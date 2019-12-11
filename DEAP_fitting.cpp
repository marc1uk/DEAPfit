/* vim:set noexpandtab tabstop=4 wrap */

#include <string>
//#include <ctime>
//#include <future>
#include <thread>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <cstring>
//#include <bitset>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TText.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TKey.h"
#include "TTimer.h"
#include "Math/MinimizerOptions.h"

// fit globals
TF1* pedestal_func=nullptr; // holds the shape of the pedestal, a TMath::Gaus
double ped_range_min=-1;    // lower range of the pedestal gaussian function
double ped_range_max=-1;    // upper range of the pedestal gaussian function
TF1* spe_func=nullptr;      // holds the shape of the SPE peak, 
double spe_range_min=-1;    // lower range of the SPE function
double spe_range_max=-1;    // upper range of the SPE function
std::vector<TF1*> npe_funcs;
std::vector<TF1Convolution*> npe_convolns;  // used to obtain the TF1s
double npe_func_min=-1;     // each npe function requires a range, but specifying many ranges seems overkill
double npe_func_max=-1;     // they'll probably all need to have a range covering the full span of charges
int max_pes=2;              // we can probably fix this?

// fitting functions
double Gamma(double* x, double* gamma_pars);
double SPE_Func(double* x, double* SPE_pars);
double FullFitFunction(double* x, double* all_pars);

int main(int argc, const char* argv[]){
    
    if(argc<2){
        std::cerr<<"usage: "<<argv[0]<<" 'file' 'source'"<<std::endl;
        std::cerr<<"file must be a full absolute path. It may be a a root macro containing a histogram, "
                 <<" or a .root file containing many histograms."<<std::endl
                 <<"The file source may be either ToolAnalysis calibration runs (0, default), "
                 <<" or the PMT Test Stand (1), specified by the 'source' argument" // TODO
                 <<std::endl;
        return 0;
    }
    
    // make the TApplication for drawing
    int myargc = 1;
    char arg1[] = "potato";
    char* myargv[]={arg1};
    TApplication *FitApp = new TApplication("FitApp",&myargc,myargv);
    
    // get the filename
    std::string filename = argv[1];
    // get it's extension to get the type
    std::string extension;
    if (filename.find(".") != std::string::npos){
        extension = filename.substr(filename.find_last_of('.'), filename.length()-filename.find_last_of('.'));
    } else {
        std::cerr<<"could not locate '.' in filename!? Can't identify filetype!"<<std::endl;
        return -1;
    }
    std::cout<<"filename is: "<<filename<<std::endl;
    std::cout<<"extension is: "<<extension<<std::endl;
    
    // get the source of file: PMT test stand or ToolAnalysis
    // this is used to determine the x-axis units, and scaling required to put in pC
    int filesource = (argc<3) ? 1 : atoi(argv[2]);
    std::string filesourcestring = (filesource) ? "ToolAnalysis" : "TestStand";
    
    // XXX DEBUG
    //filename = "/home/marc/Documents/PhD/PMTs/PMT_Test_Station/results/DEAP_refit/gain1300V_WB012.C";
    //filename = "/home/marc/Documents/PhD/PMTs/PMT_Test_Station/results/DEAP_refit/LEDRun1171S0LED6And10_PulseWindowOnly_PMTStability_Run32540.root";
    // XXX DEBUG
    
    std::map<int,TH1*> histos_to_fit;
    TCanvas* c1=nullptr;
    TH1* thehist=nullptr;
    TFile* infile=nullptr; // close when we're done
    // handle standalone histogram files
    if(extension==".C"||extension==".c"){
        // plot the histogram from file
        std::string rootcommand=".x " + filename;
        std::cout<<"Executing command "<<rootcommand<<std::endl;
        // why does this take so long to open??
        // if i move on immediately, the canvas is still loading, and then it fails to find the histogram...!
        gROOT->ProcessLine(rootcommand.c_str());
        for(int i=0; i<5; ++i){
            gSystem->ProcessEvents();
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
        gSystem->ProcessEvents();
        // now we need to pull the histogram...
        std::cout<<"searching for histogram"<<std::endl;
        c1 = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
        if(c1==nullptr){
            std::cerr<<"ROOT had no canvases after drawing the histogram"<<std::endl;
            return 0;
        }
        // extract the histogram from the canvas
        TIter nextobject(c1->GetListOfPrimitives());
        TObject* obj=nullptr;
        // loop over objects in the canvas (may include many functions owned by the histo, for example)
        for(int i=0; i<200; ++i){  // a while loop should fine, but it seems to hang sometimes
            auto&& aret = (obj=nextobject());
            //std::cout<<"ret="<<aret;
            if(not aret) break;
            //std::cout<<"primitive: "<<obj->GetName()<<" at "<<obj
            //                 <<", inherits from TH1="<<obj->InheritsFrom("TH1")<<std::endl;
            if(obj->InheritsFrom("TH1")){
                //std::cout<<"it's the histogram!"<<obj->GetName()<<std::endl;
                TH1* thehist=(TH1*)obj;
                histos_to_fit.emplace(0,thehist);
                break;
            }
        }
    // handle root files
    } else if(extension==".root"){
        // Open root file
        infile = TFile::Open(filename.c_str(),"READ");
        // pull out all the pulse charge distribution histograms
        TKey *key=nullptr;
        TIter getnextkey(infile->GetListOfKeys());
        // loop over objects in file
        std::cout<<"looping over keys"<<std::endl;
        while ((key = (TKey*)getnextkey())){
            if(key==nullptr){
                std::cerr<<"null key!"<<std::endl;
                continue;
            }
            //std::cout<<"processing key " << key->GetName()<<std::endl;
            TClass *cl = gROOT->GetClass(key->GetClassName());
            // skip unless it's a histogram
            if (!cl->InheritsFrom("TH1")) continue;
            //std::cout<<"it's a histo"<<std::endl;
            // check if it's a charge distribution
            // histogram names are "hist_charge_xxx" where xxx is the ToolAnalysis detectorkey for that PMT.
            std::string histname = key->GetName();
            int detectorkey;
            int n = 0;
            int cnt = sscanf(histname.c_str(),"hist_charge_%d %n",&detectorkey, &n);
            int success = ((n >0) && ((histname.c_str())[n]=='\0'));
            if(success){
                std::cout<<"Found charge histogram "<<histname<<std::endl;
                // read the histogram from the file into memory(?) and get a pointer
                thehist = (TH1*)key->ReadObj();
                histos_to_fit.emplace(detectorkey,thehist);
            }
        } // end loop over keys in file
        
        // also make a canvas
        c1 = new TCanvas("c1");
    } // end rootfile case
    std::cout<<"Found "<<histos_to_fit.size()<<" histograms to fit"<<std::endl;
    
    if(c1==nullptr){ std::cout<<"No canvas"<<std::endl; return 0; }
    
    // make an output file in which to save results
    TFile* outfile = new TFile("./DEAPfitter_outfile.root","RECREATE");
    outfile->cd();
    TTree* outtree = new TTree("fit_parameters","fit_parameters");
    int file_detectorkey;
    TBranch* bDetKey = outtree->Branch("DetectorKey",&file_detectorkey);
    int file_fit_success;
    TBranch* bfit_success = outtree->Branch("fit_success", &file_fit_success);
    double file_prescaling;
    double file_ped_scaling;
    double file_ped_mean;
    double file_ped_sigma;
    double file_spe_firstgamma_scaling;
    double file_spe_firstgamma_mean;
    double file_spe_firstgamma_shape;
    double file_spe_secondgamma_scaling;
    double file_spe_secondgamma_mean_scaling;
    double file_spe_secondgamma_shape_scaling;
    double file_spe_expl_scaling;
    double file_spe_expl_charge_scaling;
    double file_mean_npe;
    double file_max_pes;
    TBranch* bprescaling = outtree->Branch("prescaling", &file_prescaling);
    TBranch* bped_scaling = outtree->Branch("ped_scaling", &file_ped_scaling);
    TBranch* bped_mean = outtree->Branch("ped_mean", &file_ped_mean);
    TBranch* bped_sigma = outtree->Branch("ped_sigma", &file_ped_sigma);
    TBranch* bspe_firstgamma_scaling = outtree->Branch("spe_firstgamma_scaling", &file_spe_firstgamma_scaling);
    TBranch* bspe_firstgamma_mean = outtree->Branch("spe_firstgamma_mean", &file_spe_firstgamma_mean);
    TBranch* bspe_firstgamma_shape = outtree->Branch("spe_firstgamma_shape", &file_spe_firstgamma_shape);
    TBranch* bspe_secondgamma_scaling = outtree->Branch("spe_secondgamma_scaling", &file_spe_secondgamma_scaling);
    TBranch* bspe_secondgamma_mean_scaling = outtree->Branch("spe_secondgamma_mean_scaling", &file_spe_secondgamma_mean_scaling);
    TBranch* bspe_secondgamma_shape_scaling = outtree->Branch("spe_secondgamma_shape_scaling", &file_spe_secondgamma_shape_scaling);
    TBranch* bspe_expl_scaling = outtree->Branch("spe_expl_scaling", &file_spe_expl_scaling);
    TBranch* bspe_expl_charge_scaling = outtree->Branch("spe_expl_charge_scaling", &file_spe_expl_charge_scaling);
    TBranch* bmean_npe = outtree->Branch("mean_npe", &file_mean_npe);
    TBranch* bmax_pes = outtree->Branch("max_pes", &file_max_pes);
    // OK, file made
    
    // =========================================
    // Initialization Complete, Move to Fitting
    // =========================================
    
    // define our fit function:
    TF1* full_fit_func = new TF1("full_fit_func",FullFitFunction,0,100,14);
    full_fit_func->SetNpx(2000); // good resolution on our draws
    // it takes a *lot* of parameters...
    full_fit_func->SetParName(0,"prescaling");
    full_fit_func->SetParName(1,"ped_scaling");
    full_fit_func->SetParName(2,"ped_mean");
    full_fit_func->SetParName(3,"ped_sigma");
    full_fit_func->SetParName(4,"spe_firstgamma_scaling");
    full_fit_func->SetParName(5,"spe_firstgamma_mean");
    full_fit_func->SetParName(6,"spe_firstgamma_shape");
    full_fit_func->SetParName(7,"spe_secondgamma_scaling");
    full_fit_func->SetParName(8,"spe_secondgamma_mean_scaling");
    full_fit_func->SetParName(9,"spe_secondgamma_shape_scaling");
    full_fit_func->SetParName(10,"spe_expl_scaling");
    full_fit_func->SetParName(11,"spe_expl_charge_scaling");
    full_fit_func->SetParName(12,"mean_npe");
    full_fit_func->SetParName(13,"max_pes");
    
    // Loop over histos to fit
    int offset= (argc<4) ? 0 : atoi(argv[3]);
    int offsetcount=0;
    for(auto&& ahisto : histos_to_fit){
        // skip until requested histo
        if(offsetcount<offset){
            ++offsetcount;
            continue;
        }
        // get next histogram
        c1->Clear();
        int detkey = ahisto.first;
        thehist = ahisto.second;
        
        // check the hist and get the canvas
        if(thehist==nullptr){ std::cerr<<"No hist"<<std::endl; return 0; }
        thehist->GetListOfFunctions()->Clear(); // remove any existing fit functions
        
        // scale the hist X axis to be in picocoulombs
        if(filesourcestring=="ToolAnalysis"){
            // for these files the histograms are in nC
            double nano_to_pico = 1000;
            thehist->GetXaxis()->Set(thehist->GetNbinsX(),
                                     thehist->GetXaxis()->GetXmin()*nano_to_pico,
                                     thehist->GetXaxis()->GetXmax()*nano_to_pico);
        } else if(filesourcestring=="TestStand"){
            // for test stand results the distributions are in electron counts
            double ELECTRON_CHARGE = 1.6E-19;
            double ELECTRON_CHARGE_IN_PICOCOULOMBS = ELECTRON_CHARGE*1E12;
            thehist->GetXaxis()->Set(thehist->GetNbinsX(),
                                     thehist->GetXaxis()->GetXmin()*ELECTRON_CHARGE_IN_PICOCOULOMBS,
                                     thehist->GetXaxis()->GetXmax()*ELECTRON_CHARGE_IN_PICOCOULOMBS);
        }
        
        // We must set the ranges of all functions, and for TF1Convolution, they MUST tend to 0 at their ends!
        double histogram_minimum = thehist->GetXaxis()->GetXmin();  // this ISN'T ENOUGH
        double histogram_maximum = thehist->GetXaxis()->GetXmax();  // convolutions will be CHOPPED
        histogram_minimum = -histogram_maximum;
        //histogram_maximum *=2.;
        ped_range_min=histogram_minimum;  // lower range of the pedestal gaussian function
        ped_range_max=histogram_maximum;  // upper range of the pedestal gaussian function
        spe_range_min=histogram_minimum;  // lower range of the SPE function
        spe_range_max=histogram_maximum;  // upper range of the SPE function
        npe_func_min=histogram_minimum;   // each npe function could potentially have a different range
        npe_func_max=histogram_maximum;   // but i don't think reducing these has any benefit?
        
        // Smooth the histogram. This seems to help the fit hit the broad SPE position better
        thehist->Smooth(1.);
        
        // TODO derive this from occupancy
        double mean_pe_guess = 0.1;
        
        // FIXME dependant on the bin counts (doesn't work after we scale the histogram)
        // some very crude carry-over stuff from the triple-gaussian fit used for the PMT test stand.
        // designed to identify a few salient points in the distribution, for prior calculations
        // Try to find a second distinct peak after the pedestal
        int maxpos1=-1, maxpos2=-1, interpos=-1;
        int max1=-1, max2=-1;
        int peaktovalleymin=5; // SPE peak must be 5+ counts > the minima, to ignore noise XXX tuneme
        int intermin=std::numeric_limits<int>::max();
        // scan over X-axis
        for(int bini=1; bini<thehist->GetNbinsX(); bini++){
            // get bin content
            int bincont = thehist->GetBinContent(bini);
            // if new global maximum, set this as the pedestal location and reset any SPE/valley locations
            if(bincont>max1){ max1=bincont; maxpos1=bini; intermin=std::numeric_limits<int>::max();}
            // if we've found the pedestal peak & this is lower than our current minima, mark as valley
            else if((bincont<intermin)&&(max1>0)&&(max2==-1)){ intermin=bincont; interpos=bini; }
                     // if we've found a pedestal peak & a valley
            else if( (interpos>0) &&
                     // and this is 'sqrt(minima)' above the minima, i.e. 'large enough' to be SPE peak
                     (bincont>(intermin+2.*sqrt(intermin))) &&
                     // and this has at least 8% of the total counts, or 20 counts (for high histograms IIRC)
                     ((bincont>(thehist->GetBinContent(thehist->GetMaximumBin())*0.08)) || (bincont>20)) &&
                     // and it's higher than our current SPE peak bin
                     (bincont>max2)){
                        // record this as the new SPE bin
                        max2=bincont; maxpos2=bini;
            }
        }
        std::cout<<"Found pedestal at "<<maxpos1<<" with counts "<<max1<<std::endl;
        std::cout<<"Found valley at "<<interpos<<" with counts "<<intermin<<std::endl;
        std::cout<<"Found a SPE peak at "<<maxpos2<<" with counts "<<max2<<std::endl;
        bool spe_peak_found = ((max2>10)&&(maxpos2>0));
        if(not spe_peak_found) continue; // just don't even try
        
        // normalise the histogram to unit area? or amplitude?
        double scaling = thehist->Integral();
        //double scaling = thehist->GetBinContent(thehist->GetMaximumBin());
        thehist->Scale(1./scaling);
        double max2scaled = max2 / scaling;               // used for estimating spe amplitude
        //double max1scaled = max1/scaling;               // unused
        //double interminscaled = intermin / scaling;     // unused
        
        // clear the canvas
        c1->Clear();
        c1->cd();
        
        // draw the histogram
        thehist->Draw();
        gPad->SetLogy();
        gSystem->ProcessEvents();
        
        double pedestal_mean_guess = thehist->GetBinCenter(thehist->GetMaximumBin());
        double max_bin_count = thehist->GetBinContent(thehist->GetMaximumBin());
        double pedestal_amp_guess = max_bin_count/100.;  // TODO do better
        
        double spe_amp_guess;
        double pedestal_sigma_guess;
        // if we found a SPE bump, use it's location to estimate pedestal width and SPE amplitude
        // otherwise ... i dunno.
        if(spe_peak_found){
            pedestal_sigma_guess=thehist->GetBinCenter(interpos)/10.;  // overwrite
            spe_amp_guess = max2scaled;
        } else {
            spe_amp_guess = pedestal_amp_guess/5;
            pedestal_sigma_guess = 0.1; // based on being in pC ... TODO find a better way
        }
        std::cout<<"spe_amp_guess = "<<spe_amp_guess<<std::endl;
        
        // if we found an SPE peak use its location, otherwise take for prior just a bit above mean
        double spe_mean_guess = (spe_peak_found) ? thehist->GetBinCenter(maxpos2) : thehist->GetMean()*1.5;
        // rough order of magnitude for the SPE width
        double spe_sigma_guess = pedestal_sigma_guess*4;
        
        // Set the initial fit priors
        // ==========================
        // update the fit range
        full_fit_func->SetRange(thehist->GetXaxis()->GetXmin(),thehist->GetXaxis()->GetXmax());
        
        std::cout<<"Setting parameter limits"<<std::endl;
        // we can't set more than 11 parameters inidividually, we need to build them into an array
        std::vector<double> fit_parameters(full_fit_func->GetNpar());
        
        // general parameters
        // ------------------
        fit_parameters.at(0) = 5.;              // overall scaling of the entire function
        // TODO tune from histogram?
        // TODO this seems redundant given that we also have scaling for the individual components;
        // try removing/fixing this parameter to speed up fitting
        fit_parameters.at(12) = mean_pe_guess;  // mean number of pe
        full_fit_func->SetParLimits(12,0,0.5);  // >1 we're not in SPE regime any more!
        fit_parameters.at(13) = max_pes;
        // TODO A better prior would be the occupancy: number of flashes with a hit over total num flashes
        // we'd need to do a threshold check on each waveform before adding to the pulse charge distribution
        
        // pedestal parameters
        // -------------------
        fit_parameters.at(1) = pedestal_amp_guess;   // pedestal amplitude (gaussian scaling)
        // TODO: maybe we should normalise our histogram by integral, to better constrain our input ranges
        full_fit_func->SetParLimits(1,0,max_bin_count);
        fit_parameters.at(2) = pedestal_mean_guess;  // 0, provided baseline-subtracted, but it isn't always
        full_fit_func->SetParLimits(2,ped_range_min,ped_range_max);
        fit_parameters.at(3) = pedestal_sigma_guess; // pedestal width
        full_fit_func->SetParLimits(3,0.,ped_range_max/10);
        
        // spe parameters
        // --------------
        // the SPE consists of two 'Gaus' functions and an exponential
        
        // first gamma parameters
        std::cout<<"setting spe_firstgamma_scaling to "<<spe_amp_guess<<std::endl;
        fit_parameters.at(4) = spe_amp_guess;  // controls amplitude of main SPE function bump, indirectly
        // since we have two contributing gauss functions, and a Gaus isn't normalized
        full_fit_func->SetParLimits(4,0,max_bin_count); // TODO check upper limit isn't too restrictive
        fit_parameters.at(5) = spe_mean_guess;     // controls position of peak via a stretching from the LH edge
        full_fit_func->SetParLimits(5,0,spe_range_max);
        fit_parameters.at(6) = 10;  // an arbitrary shape factor. (1/b in the paper)
        // controls how skewed to the left (low values) or symmetric (higher values) the peak is.
        // more symmetric (higher) values also narrow the distribution.
        // minimum of 1, almost a triangle vertical at the left edge, max >100, for nearly gaussian
        full_fit_func->SetParLimits(6,1.,100);
        
        // secondary gamma parameters
        // second gamma's parameters are scaling factors to those of the first
        // secondary gamma should be at LOWER mean and slightly WIDER shape: both scaling factors should be <1.
        fit_parameters.at(7) = 0.3;                              // amplitude scaling of secondary Gamma
        full_fit_func->SetParLimits(7,0,0.75);                   // amp↑ as q↓, so even scaling 1 is too high
        fit_parameters.at(8) = 0.3;                              // mean scaling. seemed to work...
        fit_parameters.at(9) = 0.6;                              // shape scaling. seemed to work...
        full_fit_func->SetParLimits(8,0.1,1);                    // secondary gamma should be lower charge
        full_fit_func->SetParLimits(9,0.1,2);                    // secondary gammma should be wider
        
        // exponential parameters
        // SPE expl on a log-plot is a straight line dropping off from pedestal to fill the dip region
        // at most this should not exceed the dip where it fills in
        double exp_charge_scaling = 10;   // expl charge scaling: defines steepness of the expl fill-in
        fit_parameters.at(11) = exp_charge_scaling;
        full_fit_func->SetParLimits(11,3,50);
        // expl amplitude scaling defines a vertical shift of the expl fill-in
        // zero would probably be a fair prior... but let's try to get something if we have a dip? 
        if(spe_peak_found){
            double dip_charge  = thehist->GetBinCenter(interpos);
            double exp_scaling = (1./exp_charge_scaling)*exp(-dip_charge*exp_charge_scaling); // ~ dip counts
            exp_scaling /= 10;  // back it off
            fit_parameters.at(10) = exp_scaling;
        } else {
            fit_parameters.at(10) = 0;
        }
        full_fit_func->SetParLimits(10,0,max_bin_count);
        
        // notes on above calculation:
        // main scaling should be such that it is close to the number of counts at the minium of the dip region
        // i.e. at the q of the dip then (exp_scaling*charge_scaling)exp(-q_dip*charge_scaling) = dip_counts
        // so q_dip = 0.1, dip_counts = 200, with charge_scaling = 10, we have 10*exp_scaling*exp(-1)=200
        // or exp_scaling = 20/0.36 = 55. Back it off a bit to 40....
        // BUT actually, if this contributes too much, then its cutoff at the primary Gamma mean gives rise to a
        // noticable discontinuity, so back it off a lot to 10 or even 5
        // (the rest can be filled in by the secondary Gamma)
        
        // Set the fit priors
        std::cout<<"Setting parameter priors"<<std::endl;
        full_fit_func->SetParameters(fit_parameters.data());
        full_fit_func->FixParameter(13, full_fit_func->GetParameter(13));    // fix max_pes
        
        std::cout<<"Piors and their limits are:"<<std::endl;
        for(int pari=0; pari<full_fit_func->GetNpar(); ++pari){
            double lbound, ubound, value;
            full_fit_func->GetParLimits(pari,lbound,ubound);
            value = full_fit_func->GetParameter(pari);
            std::string parname = full_fit_func->GetParName(pari);
            std::cout<<pari<<": ("<<parname<<") "<<lbound<<" < "<<value<<" < "<<ubound<<std::endl;
        }
        std::cout<<std::endl;
        
        // draw, to see our prior
        std::cout<<"Drawing prior"<<std::endl;
        full_fit_func->Draw("same");
        c1->Modified();
        c1->Update();
        gSystem->ProcessEvents();
        /*  wait for user to view and fiddle
        while(gROOT->FindObject("c1")!=nullptr){
            gSystem->ProcessEvents();
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        c1 = new TCanvas("c1");
        thehist->Draw();
        */
        /*  this worked at one point in toolanalysis, but now it segfaults... :(
        TTimer* timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
        timer->TurnOn();
        timer->Reset();
        std::cout<<"Type <return> to exit:"<<std::endl;
        std::string trashme;
        std::getline (std::cin,trashme);
        timer->TurnOff();
        delete timer;
        timer = nullptr;
        */
        //std::this_thread::sleep_for(std::chrono::milliseconds(500));
        
        // FIXME doesn't add the fit function...
//        // write the histo to file, with prior fit function, for debug
//        std::cout<<"Writing prior fit to file"<<std::endl;
//        thehist->Write(TString::Format("%s_prior",thehist->GetName()));
        
        // We may be able to speed up fitting by reducing tolerances or changing out integration strategy
        // from https://root-forum.cern.ch/t/speeding-up-fitting-to-a-landau-distribution/25140/2
        //ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
        //ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);
        
        // Try to do the fit
        std::cout<<"Fitting"<<std::endl;
        auto tstart = std::chrono::high_resolution_clock::now();
        int fit_success = thehist->Fit(full_fit_func);
        auto tend = std::chrono::high_resolution_clock::now();
        double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();
        time_taken *= 1e-9; // nano seconds to seconds
        std::cout<<"Fitting took "<<time_taken<<" sec"<<std::endl;
        std::cout<<"Fit returned: "<<fit_success<<std::endl;
        std::vector<double> fitted_parameters(full_fit_func->GetNpar());
        if(fit_success==0){
            std::cout<<"Fit success on histogram "<<thehist->GetName()<<std::endl;
            std::cout<<"Fit parameters were: {";
            for(int pari=0; pari<full_fit_func->GetNpar(); ++pari){
                std::cout<<full_fit_func->GetParameter(pari);
                if((pari+1)<full_fit_func->GetNpar()) std::cout<<", ";
                fitted_parameters.push_back(full_fit_func->GetParameter(pari));
            }
            std::cout<<"}"<<std::endl;
        }
        
        // mark the canvas as modified by the fit
        std::cout<<"Updating canvas with fit"<<std::endl;
        c1->Modified();
        c1->Update();
        gSystem->ProcessEvents();
        /* Debug: allow user to fiddle
        while(gROOT->FindObject("c1")!=nullptr){
            gSystem->ProcessEvents();
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        c1 = new TCanvas("c1");
        thehist->Draw();
        */
        /*  this worked at one point in toolanalysis, but now it segfaults... :(
        timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
        timer->TurnOn();
        timer->Reset();
        std::cout<<"Type <return> to exit:"<<std::endl;
        std::getline (std::cin,trashme);
        timer->TurnOff();
        delete timer;
        timer = nullptr;
        */
        //std::this_thread::sleep_for(std::chrono::milliseconds(500));
        
        // write the fitted histo to file, with fit function
        std::cout<<"Writing to file"<<std::endl;
        thehist->Write(thehist->GetName());
        
        // copy it all into our ROOT file
        file_detectorkey = detkey;
        file_fit_success = fit_success;
        file_prescaling = fitted_parameters.at(0);
        file_ped_scaling = fitted_parameters.at(1);
        file_ped_mean = fitted_parameters.at(2);
        file_ped_sigma = fitted_parameters.at(3);
        file_spe_firstgamma_scaling = fitted_parameters.at(4);
        file_spe_firstgamma_mean = fitted_parameters.at(5);
        file_spe_firstgamma_shape = fitted_parameters.at(6);
        file_spe_secondgamma_scaling = fitted_parameters.at(7);
        file_spe_secondgamma_mean_scaling = fitted_parameters.at(8);
        file_spe_secondgamma_shape_scaling = fitted_parameters.at(9);
        file_spe_expl_scaling = fitted_parameters.at(10);
        file_spe_expl_charge_scaling = fitted_parameters.at(11);
        file_mean_npe = fitted_parameters.at(12);
        file_max_pes = fitted_parameters.at(13);
        outtree->Fill();
        outtree->Write("fit_parameters",TObject::kOverwrite);
        
        //break; // XXX DEBUG so we can check
        
    } // end loop over histograms in file
    
    // Cleanup
    if(infile){
        infile->Close();
        delete infile;
        infile = nullptr;
    }
    if(outfile){
        outfile->Write("",TObject::kOverwrite);
        outtree->ResetBranchAddresses();
        outfile->Close();
        delete outfile;
        outfile = nullptr;
    }
    
    delete FitApp;
}

// ========================================================================
// FIT FUNCTION FROM DEAP-3600
// arXiv:1705.10183v3
// ========================================================================

double Gamma(double* x, double* gamma_pars){
    if((*x)<0) return 0; // goes nuts below 0... :/
    double gamma_mean = gamma_pars[0];
    double gamma_shape = gamma_pars[1]; // this is 1/b in the paper; only it's inverse is ever used
    double gamma_product = gamma_mean/gamma_shape;
    double x_over_gamma_product = (*x)/gamma_product;
    double returnval = 1./(gamma_product*TMath::Gamma(gamma_shape));
    returnval *= pow(x_over_gamma_product,gamma_shape-1.);
    returnval *= exp(-x_over_gamma_product);
    if(returnval<0) returnval=0;
    return returnval;
}

double SPE_Func(double* x, double* SPE_pars){
    if((*x)<0) return 0; // safer, don't do anything weird...
    // the SPE function is a sum of a primary Gamma function for normal-path electrons
    // plus a second Gamma function for electrons that miss the first dynode and hit the second
    // plus an exponential that describes "scattering of the photoelectron on the dynode structure"
    // The two gamma functions are linked, such that the mean and shape of the second Gamma
    // is defined by a scaling factor multiplied by those of the first
    double firstgamma_scaling = SPE_pars[0];
    double firstgamma_mean = SPE_pars[1];
    double firstgamma_shape = SPE_pars[2];
    double secondgamma_scaling = SPE_pars[3];
    double secondgamma_mean_scaling = SPE_pars[4];
    double secondgamma_shape_scaling = SPE_pars[5];
    double expl_scaling = SPE_pars[6];
    double expl_charge_scaling = SPE_pars[7];
    
    double secondgamma_pars[]{secondgamma_mean_scaling*firstgamma_mean,secondgamma_shape_scaling*firstgamma_shape};
    
    double returnval = 0;
    returnval  = firstgamma_scaling*Gamma(x, &SPE_pars[1]);
    returnval += secondgamma_scaling*firstgamma_scaling*Gamma(x, &secondgamma_pars[0]);
    if( ((*x)<firstgamma_mean) && ((*x)>0) )
        returnval += (expl_scaling*expl_charge_scaling)*exp(-(*x)*expl_charge_scaling);
    // ^ this 'cutoff at the mean value of the primary SPE gamma' gives a discontinuity which is
    // noticable if the expl function is too large of a component... use small parameters
    // added a cutoff of the expl at 0; we need to be able to extend the SPE function out to
    // somewhere that it should tail out to 0 for the convolution. Seems to tie up with paper.
    // Note this still gives a discontinuity at the end of the function, which TF1Convolution
    // says could be a problem, but it seems to be okay.
    
    return returnval;
}

double FullFitFunction(double* x, double* all_pars){
    double prescaling_B = all_pars[0];
    double ped_scaling = all_pars[1];
    double* ped_pars = &all_pars[2]; // parameters 2 and 3 are pedestal gaussian mean and sigma
    double* spe_pars = &all_pars[4]; // parameters 4 to 11 are parameters that describe the SPE shape
    double mean_npe = all_pars[12];  // can we pull initial fit value from occupancy?
    double max_pes = all_pars[13];
    
    // build the TF1 representing the pedestal component, if we have not yet done so
    if(not pedestal_func){
        std::cout<<"Creating Pedestal TF1"<<std::endl;
        // TODO 'true' at end of TMath::Gaus adds 1/2*pi*sigma scaling to gaus, which matches paper function.
        // however for small widths, this results in a VERY large pedestal amplitude....
        // maybe we should set this to false (maximum=1) and let the firstgamma_scaling do the work...
        pedestal_func = new TF1("pedestal_func","TMath::Gaus(x, [0], [1], true)", ped_range_min, ped_range_max);
        pedestal_func->SetParName(0,"pedestal_mean");  // true adds normalization by 1/(sqrt(2pi)*sigma)
        pedestal_func->SetParName(1,"pedestal_width"); // matches DEAP paper
    }
    // update it's parameters
    pedestal_func->SetRange(ped_range_min, ped_range_max);
    pedestal_func->SetParameters(ped_pars);
    
    // the contribution from pedestal alone has its own scaling
    double pedestal_contribution = ped_scaling*pedestal_func->Eval(*x);
    
    // build the TF1 representing a SPE peak, if we have not yet done so
    if(not spe_func){
        std::cout<<"Creating SPE TF1"<<std::endl;
        spe_func = new TF1("spe_func",SPE_Func,spe_range_min,spe_range_max,8); // 8 is num parameters it takes...!
        spe_func->SetParName(0,"spe_firstgamma_scaling");
        spe_func->SetParName(1,"spe_firstgamma_mean");
        spe_func->SetParName(2,"spe_firstgamma_shape");
        spe_func->SetParName(3,"spe_secondgamma_scaling");
        spe_func->SetParName(4,"spe_secondgamma_mean_scaling");
        spe_func->SetParName(5,"spe_secondgamma_shape_scaling");
        spe_func->SetParName(6,"spe_expl_scaling");
        spe_func->SetParName(7,"spe_expl_charge_scaling");
    }
    // update it's parameters
    spe_func->SetRange(spe_range_min,spe_range_max);
    spe_func->SetParameters(spe_pars);
    
    double running_npe_contribution=0;
    // loop over the NPE peaks
    for(int i=0; i<max_pes; ++i){
        
        // if we do not have a TF1 corresponding to this NPE peak - build it
        if(npe_funcs.size()<(i+1)){
            // each Npe peak is a convolution of the pedestal, together with the SPE peak convolved with itself N times
            // we re-use the results of previous convolutions to build on for the next
            TF1* n_minus_one_func = (i==0) ? pedestal_func : npe_funcs.at(i-1);
            // convolve it with the SPE function again
            std::cout<<"Creating NPE convolution "<<i<<std::endl;
            TF1Convolution* new_conv = new TF1Convolution(n_minus_one_func, spe_func, npe_func_min, npe_func_max, true);
            //new_conv->SetNofPointsFFT(1000);  // default, should be at least this
            npe_convolns.push_back(new_conv);
            // make the function from the convolution
            std::cout<<"Creating NPE TF1 "<<i<<std::endl;
            TString funcname = TString::Format("npe_func_%d",i);
            TF1* npe_func = new TF1(funcname,*new_conv, npe_func_min, npe_func_max, new_conv->GetNpar());
            npe_funcs.push_back(npe_func);
        }
        
        // retrieve the TF1 descrbing this NPE peak
        TF1Convolution* this_npe_conv = npe_convolns.at(i);
        TF1* this_npe_func = npe_funcs.at(i);
        this_npe_conv->SetRange(npe_func_min, npe_func_max);
        this_npe_func->SetRange(npe_func_min, npe_func_max);
        
        // set it's parameters - the TF1Convolution of 2 functions (or the TF1 derived from it)
        // combines the parameters of both TF1s that went into it
        // the first couple of parameters are for the pedestal function
        this_npe_func->SetParameters(0,ped_pars[0]);
        this_npe_func->SetParameter(1,ped_pars[1]);
        // after that we keep re-convolving with the same SPE function, with the same parameters
        int n_spe_pars = spe_func->GetNpar();
        for(int j=0; j<(i+1); ++j){
            for(int pari=0; pari<n_spe_pars; ++pari){
                this_npe_func->SetParameter(2+(j*n_spe_pars)+pari,spe_pars[pari]);
                this_npe_func->SetParameter(2+(j*n_spe_pars)+pari,spe_pars[pari]);
            }
        }
        // this should now describe the shape of the distribution of the NPE peak
        double this_npe_peak_contribution = this_npe_func->Eval(*x);
        
        // each NPE peak will have an occupancy scaling
        // based on a Poisson distribution parameterised by the mean number of PE on the tube
        double poisson_scaling = TMath::Poisson(i+1,mean_npe);  // +1 because i=0 is 1PE, not 0PE
        // evaluate this npe's contribution at the specified charge
        running_npe_contribution += (poisson_scaling*this_npe_peak_contribution);
        
    }
    
    // the total number of events with a given charge is then the sum of contributions from pedestal
    // plus SPE plus 2PE plus 3PE ... 
    double total_contribution = pedestal_contribution + running_npe_contribution;
    // and there's an overall scaling factor
    return (prescaling_B*total_contribution);
}
