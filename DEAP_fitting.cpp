/* vim:set expandtab tabstop=4 wrap */

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
#include "TFitResult.h"
#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TKey.h"
#include "TTimer.h"
//#include "Math/MinimizerOptions.h"
#include "TVirtualFitter.h"

#include "DEAPFitFunction.h"
#include "sarge.h"

#ifndef STORE_PRIORS
#define STORE_PRIORS 1
#endif

#ifndef STORE_DEBUG_EXTRAS
#define STORE_DEBUG_EXTRAS 1
#endif

constexpr double ELECTRON_CHARGE = 1.6E-19;
constexpr double ELECTRON_CHARGE_IN_PICOCOULOMBS = ELECTRON_CHARGE*1E12;
constexpr int max_pes = 3;

int main(int argc, const char* argv[]){
    
    // =========================================
    // INPUT PARSING WITH SARGE
    // =========================================
    
    // Create a Sarge to parse arguments
    Sarge sarge;
    // configure sarge with this programs arguments
    sarge.setArgument("h", "help", "Get help.", false);
    sarge.setArgument("t", "type", "Did the input files come from the PMT Test Stand (0) or ToolAnalysis (1, default)", true);
    sarge.setArgument("s", "start", "Skip this many histograms before starting (default 0)", true);
    sarge.setArgument("n", "numhistos", "Process this many histograms before exiting (0, default=all from start)", true);
    sarge.setArgument("p", "precut", "Should offset and numhistos consider all histograms (0) or only those passing the PeakScan pre-cut (1, default)", true);
    sarge.setArgument("o", "output", "Name of output file (default 'DEAPfitter_outfile.root')", true);
    sarge.setArgument("j", "justscan", "Just do PeakScan and print histograms viable for fitting", false);
    sarge.setDescription("Program to try to fit PMT pulse charge distributions in order to extract gain");
    sarge.setUsage("deapfit <options> <inputfiles>");
    
    if(!sarge.parseArguments(argc, argv)){
        std::cerr << "Couldn't parse arguments..." << std::endl;
        return 1;
    }
    
    if(sarge.exists("help")){ sarge.printHelp(); return 0; }
    
    // get the type of input file
    // this is used to determine the x-axis units, and scaling required to put in pC
    std::string filesourcestring = "ToolAnalysis";
    std::string requestedtype;
    if(sarge.getFlag("type", requestedtype)){
        if     (requestedtype=="0"){ filesourcestring = "TestStand";    }
        else if(requestedtype=="1"){ filesourcestring = "ToolAnalysis"; }
        else { std::cerr<<"Unrecognised input type: "<<requestedtype
                        <<", options are 0: PMT Test Stand, 1: ToolAnalysis"<<std::endl;
               return 1;
        }
    }
    
    // get number of histograms to analyse; default all of them
    int histos_to_analyse = 0;  // all of them
    std::string requestednhistos;
    if(sarge.getFlag("numhistos", requestednhistos)){
        histos_to_analyse = stoi(requestednhistos);
    }
    std::cout<<"Will analyse ";
    if(histos_to_analyse>0) std::cout<< histos_to_analyse << " histograms ";
    else                    std::cout<<"all histograms ";
    
    // get the number of the first histogram from which to start processing
    int offset = 0;
    std::string requestedoffset;
    if(sarge.getFlag("start", requestedoffset)){
        offset = stoi(requestedoffset);
    }
    std::cout<<"starting from ";
    if(offset>0) std::cout<< "histogram " << offset;
    else         std::cout<< "the first histogram";
    
    // get whether offset and histos_to_analyse should consider only those passing PeakScan pre-cut
    int precut = 1; // true, only those passing precut
    std::string requestedprecut;
    if(sarge.getFlag("precut", requestedprecut)){
        precut = stoi(requestedprecut);
    }
    std::cout<<", considering ";
    if(precut) std::cout<<"only those that pass PeakScan check"<<std::endl;
    else       std::cout<<"all histograms, regardless of their PreScan result"<<std::endl;
    
    // get output filename
    std::string outputfilename = "DEAPfitter_outfile.root"; // default
    sarge.getFlag("output", outputfilename);
    
    // check if we're just doing a pre-scan
    int prescan_only = sarge.exists("justscan");
    if(prescan_only){
        // strip off extention
        std::string filenamewithoutextention = outputfilename;
        if(outputfilename.find(".") != std::string::npos){
            filenamewithoutextention = outputfilename.substr(0, outputfilename.find_last_of('.'));
        }
        // replace with .txt
        outputfilename = filenamewithoutextention+".txt";
        std::cout<<"Will output a list of histograms passing PeakScan to "<<outputfilename<<std::endl;
    } else {
        std::cout<<"Results will be written to "<<outputfilename<<std::endl;
    }
    
    // get input files
    std::vector<std::string> inputfiles;
    std::string inputfile;
    int filei=0;
    while(sarge.getTextArgument(filei, inputfile)){
        inputfiles.push_back(inputfile);
        filei++;
    }
    
    // =========================================
    // END OF INPUT PARSING
    // MAKE THE OUTPUT FILE
    // =========================================
    
    // declare variables for proper fitting output
    TFile* outfile=nullptr;
    TTree* outtree=nullptr;
    
    int RunNum=0;
    int SubRun=0;
    int Voltage=0;
    int file_detectorkey;
    int file_just_gaussian;
    int file_found_spe_peak;
    int file_fit_success;
    double file_fit_chi2;
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
    double file_mean_spe_charge;
    double file_gain;
    double file_spe_firstgamma_gain; // for reference...
#ifdef STORE_PRIORS
    double file_prescaling_prior;
    double file_ped_scaling_prior;
    double file_ped_mean_prior;
    double file_ped_sigma_prior;
    double file_spe_firstgamma_scaling_prior;
    double file_spe_firstgamma_mean_prior;
    double file_spe_firstgamma_shape_prior;
    double file_spe_secondgamma_scaling_prior;
    double file_spe_secondgamma_mean_scaling_prior;
    double file_spe_secondgamma_shape_scaling_prior;
    double file_spe_expl_scaling_prior;
    double file_spe_expl_charge_scaling_prior;
    double file_mean_npe_prior;
    double file_max_pes_prior;
#endif
#ifdef STORE_DEBUG_EXTRAS
    int file_maxpos1;
    int file_maxpos2;
    int file_interpos;
    int file_max1;
    int file_max2;
    int file_peaktovalleymin;
    int file_intermin;
    double file_yscaling;
#endif
    
    // or in the case of just prescan
    std::ofstream passingHistos;
    
    if(not prescan_only){
        // make an output file in which to save results
        outfile = new TFile(outputfilename.c_str(),"RECREATE");
        outfile->cd();
        outtree = new TTree("fit_parameters","fit_parameters");
        TBranch* bDetKey = outtree->Branch("DetectorKey",&file_detectorkey);
        TBranch* bRunNum = outtree->Branch("Run",&RunNum);
        TBranch* bSubRun = outtree->Branch("SubRun",&SubRun);
        TBranch* bVoltage = outtree->Branch("Voltage",&Voltage);
        TBranch* bjust_gaussian = outtree->Branch("just_gaussian", &file_just_gaussian);
        TBranch* bfound_spe_peak = outtree->Branch("found_spe_peak", &file_found_spe_peak);
        TBranch* bfit_success = outtree->Branch("fit_result", &file_fit_success);
        TBranch* bfit_chi2 = outtree->Branch("fit_chi2", &file_fit_chi2);
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
        
        // the money numbers
        TBranch* bmean_spe_charge = outtree->Branch("mean_spe_charge", &file_mean_spe_charge);
        TBranch* bgain = outtree->Branch("gain", &file_gain);
        TBranch* bspe_firstgamma_gain = outtree->Branch("spe_firstgamma_gain", &file_spe_firstgamma_gain);
        
        // also store all our priors for debug
#ifdef STORE_PRIORS
        TBranch* bprescaling_prior = outtree->Branch("prescaling_prior", &file_prescaling_prior);
        TBranch* bped_scaling_prior = outtree->Branch("ped_scaling_prior", &file_ped_scaling_prior);
        TBranch* bped_mean_prior = outtree->Branch("ped_mean_prior", &file_ped_mean_prior);
        TBranch* bped_sigma_prior = outtree->Branch("ped_sigma_prior", &file_ped_sigma_prior);
        TBranch* bspe_firstgamma_scaling_prior = outtree->Branch("spe_firstgamma_scaling_prior", &file_spe_firstgamma_scaling_prior);
        TBranch* bspe_firstgamma_mean_prior = outtree->Branch("spe_firstgamma_mean_prior", &file_spe_firstgamma_mean_prior);
        TBranch* bspe_firstgamma_shape_prior = outtree->Branch("spe_firstgamma_shape_prior", &file_spe_firstgamma_shape_prior);
        TBranch* bspe_secondgamma_scaling_prior = outtree->Branch("spe_secondgamma_scaling_prior", &file_spe_secondgamma_scaling_prior);
        TBranch* bspe_secondgamma_mean_scaling_prior = outtree->Branch("spe_secondgamma_mean_scaling_prior", &file_spe_secondgamma_mean_scaling_prior);
        TBranch* bspe_secondgamma_shape_scaling_prior = outtree->Branch("spe_secondgamma_shape_scaling_prior", &file_spe_secondgamma_shape_scaling_prior);
        TBranch* bspe_expl_scaling_prior = outtree->Branch("spe_expl_scaling_prior", &file_spe_expl_scaling_prior);
        TBranch* bspe_expl_charge_scaling_prior = outtree->Branch("spe_expl_charge_scaling_prior", &file_spe_expl_charge_scaling_prior);
        TBranch* bmean_npe_prior = outtree->Branch("mean_npe_prior", &file_mean_npe_prior);
        TBranch* bmax_pes_prior = outtree->Branch("max_pes_prior", &file_max_pes_prior);
#endif  // STORE_PRIORS
#ifdef STORE_DEBUG_EXTRAS
        // also store some of the parameters used in determining priors, in case this is useful
        TBranch* bmaxpos1 = outtree->Branch("maxpos1", &file_maxpos1);
        TBranch* bmaxpos2 = outtree->Branch("maxpos2", &file_maxpos2);
        TBranch* binterpos = outtree->Branch("interpos", &file_interpos);
        TBranch* bmax1 = outtree->Branch("max1", &file_max1);
        TBranch* bmax2 = outtree->Branch("max2", &file_max2);
        TBranch* bpeaktovalleymin = outtree->Branch("peaktovalleymin", &file_peaktovalleymin);
        TBranch* bintermin = outtree->Branch("intermin", &file_intermin);
        TBranch* byscaling = outtree->Branch("scaling", &file_yscaling);
#endif // STORE_DEBUG_EXTRAS
        // OK, file made
    } else {
        // prescan only
        passingHistos.open(outputfilename.c_str(), std::ios::out);
    }
    
    // ==================================
    // END OF OUTPUT FILE CREATION
    // MOVE TO LOOPING OVER INPUT FILES
    // ==================================
    
    // make the TApplication for drawing
    int myargc = 1;
    char arg1[] = "potato";
    char* myargv[]={arg1};
    TApplication *FitApp = new TApplication("FitApp",&myargc,myargv);
    
    // misc variables to keep track of things in each file
    std::map<int,TH1*> histos_to_fit;
    TCanvas* c1=nullptr;
    TH1* thehist=nullptr;
    TFile* infile=nullptr;    // close when we're done
    std::map<int,int> types;  // types of input files, .C or .root. Just to check we don't mix.
    int analysed_histo_count=0;
    int offsetcount=0;
    int histos_with_a_peak=0;
    
    // construct the fit function
    DEAPFitFunction deapfitter(max_pes);
    deapfitter.MakeFullFitFunction();
    
    // scan over input files and pull the histograms we're going to analyse
    for(std::string& filename : inputfiles){
        
        histos_to_fit.clear();
        
        // ==================================
        // SCAN FOR HISTOGRAMS IN THIS FILE
        // ==================================
        
        // get next file's extension to get the type
        std::string extension;
        if(filename.find(".") != std::string::npos){
            extension = filename.substr(filename.find_last_of('.'), filename.length()-filename.find_last_of('.'));
        } else {
            std::cerr<<"could not locate '.' in filename!? Can't identify filetype!"<<std::endl;
            return -1;
        }
        std::cout<<"filename is: "<<filename<<std::endl;
        std::cout<<"extension is: "<<extension<<std::endl;
        // handle standalone histogram files
        if(extension==".C"||extension==".c"){
            if(types.count(1)){
                std::cerr<<"Please don't mix input files of type .C and .root, detector keys will be meaningless."
                         <<std::endl;
                return 1;
            }
            if(types.size()==0) types.emplace(0,1);
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
                    // stop ROOT from deleting it once we draw the next histogram over it...
                    thehist->SetDirectory(0);
                    thehist->ResetBit(kCanDelete);
                    // but then no-one deletes these. That's probably fine? They'll be cleaned up on terminate...
                    // we could track them and delete them as we go, but does it offer
                    // any advantage since we pool all our histos at the start? TODO improve
                    histos_to_fit.emplace(histos_to_fit.size(),thehist); // XXX detkey is just histogram number!
                    break;
                }
            }
        // handle root files
        } else if(extension==".root"){
            if(types.count(0)){
                std::cerr<<"Please don't mix input files of type .C and .root, detector keys will be meaningless."
                         <<std::endl;
                return 1;
            }
            if(types.size()==0) types.emplace(1,1);
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
            
            // also make a canvas if we don't have one
            if(c1==nullptr) c1 = new TCanvas("c1");
        } // end rootfile case
        std::cout<<"Found "<<histos_to_fit.size()<<" histograms to fit"<<std::endl;
        
        if(c1==nullptr){ c1 = new TCanvas("c1"); }
        
        // try to pull the voltage and run number from the ToolAnalysis filename
        // expect something like: R1214S0_1000V_PMTStability_Run0.root
        int cnt = sscanf(filename.c_str(),"R%dS%d_%dV",&RunNum,&SubRun,&Voltage);
        
        // =========================================
        // END OF HISTOGRAM SCAN
        // MOVE TO FITTING HISTOGRAMS
        // =========================================
        
        // Loop over histos to fit
        for(auto&& ahisto : histos_to_fit){
            // get next histogram
            int detkey = ahisto.first;
            thehist = ahisto.second;
            if(thehist==nullptr){ std::cerr<<"No hist"<<std::endl; return 0; }
            
            if(precut==0){ // offset & histos_to_analyse consider ALL histograms, not just those viable for fit
                // in this case increment our current histo index and skip if less than our starting index
                if(offsetcount<offset){
                    ++offsetcount;
                    continue;
                }
            }
            
            // First Viability Check
            // =====================
            // matt's check for pedestal only histos:
            // integrate a gaussian, of width equal to the typical baseline sigma, from mean + 1sigma upwards
            // if <150 tag it as "empty"
            bool just_gaussian=true;
            // for the initial baseline 1-sigma, we need to do a gaussian fit.
            // To try to fit only that pedestal gaussian and not be pulled out by any
            // potential NPE peaks or extended tail, only fit from the peak out to
            // the point where bin counts drop to 10% the maximum on either side.
            // This is a little arbitrary, but hopefully will get the job done.
            // 
            // first find that fitting cutoff point
            int maxbin = thehist->GetMaximumBin();
            double maxbincount=thehist->GetBinContent(maxbin);
            double rightthresh=0;
            int rightthreshbin=0;
            for(int bini=maxbin; bini<thehist->GetNbinsX(); ++bini){
                double thisbincount = thehist->GetBinContent(bini);
                if(thisbincount<(0.1*maxbincount)){
                    //cout<<"break at "<<bini<<endl;
                    rightthresh = thehist->GetBinCenter(bini);
                    rightthreshbin = bini; // for later
                    break;
                }
            }
            double leftthresh=0;
            for(int bini=maxbin; bini<thehist->GetNbinsX(); --bini){
                double thisbincount = thehist->GetBinContent(bini);
                if(thisbincount<(0.1*maxbincount)){
                    //cout<<"break at "<<bini<<endl;
                    leftthresh = thehist->GetBinCenter(bini);
                    break;
                }
            }
            // do the fit, hopefully fitting only the pedestal.
            TF1* thegaus = new TF1("pedgaus", "gaus",leftthresh,rightthresh);
            // "gaus" is [0]*exp(-0.5*((x-[1])/[2])**2) 
            thegaus->SetNpx(2000);
            thehist->Fit(thegaus,"RQN"); // function range, quiet, do not store
            
            // Now we want to see how much of the histogram is 'within' this fit
            // and how much of it leaks out (i.e. any counts in potential NPE peaks).
            // The tricky bit here is how to suitably measure the amount 'in' and 'out'.
            // Number of counts are stats dependant. Percentage of total counts depdend on
            // the mean number of pe's, and low % doesn't necessarily mean no peaks.
            // 
            // As a shot in the dark, i'm going to count how many bins have more than X entries,
            // counting from the point where we leave the gaussian, until the first empty bin.
            // Any reasonable binning should have an SPE peak that spans at least Y bins,
            // (very coarse binning), while very fine binning should be unlikely to have X entries in
            // Y consecutive bins from just noise (also cutting off at first 0 should help).
            // I choose X = 5, Y = 10;
            int PED_COUNT_THRESH=5; // bins must have at least this many counts to be included
            int MIN_SPE_BINS=10;    // we must have this many such bins to suspect we have an SPE bump
            // do the scan
            int binsoutside=0;
            // scan over bins, might as well start from RH edge of pedestal
            for(int i=rightthreshbin; i<thehist->GetNbinsX(); ++i){
                // check if we're 'outside' the gaussian yet.
                // gonna do this by scaling the pedestal function up a bit, and checking if we're below it.
                // this captures bins in the pedestal region, but cuts off very sharply at the pedestal edge
                double binconts = thehist->GetBinContent(i);
                double fitcounts = 10.*thegaus->Eval(thehist->GetBinCenter(i));
                if(fitcounts>binconts) continue;  // still under the gaussian
                else if(binconts==0) break;       // we've hit the floor
                else if(binconts>PED_COUNT_THRESH) binsoutside++;  // another valid bin outside
                if(binsoutside>MIN_SPE_BINS){ // if we have enough bins with something in
                    just_gaussian=false;          // move to the next stage of viability checks
                    break;
                }
            }
            delete thegaus;
            
            // Second Viability Check
            // ======================
            // if we have enough bins to work with, try to find a peak & valley in them
            bool found_spe_peak=false;
            std::vector<double> precheck_pars; // access to PeakScan results if we pass the scan
            if(not just_gaussian){
                // pass to the fitter
                deapfitter.SetHisto(thehist);
                deapfitter.SmoothHisto();  // This seems to help the fit hit the broad SPE position better
                // ↑ FIXME for histos with a small number of bins TH1::Smooth ADDS valleys!
                // hopefully these are filtered out by the previous scan....
                
                // Scan for a second peak after pedestal
                found_spe_peak = deapfitter.PeakScan(&precheck_pars);
            }
            
            // update our counter of viable histograms if it passes both pre-test checks
            if(found_spe_peak){
                histos_with_a_peak++;
                // if all we're doing is counting viable histograms, we can move on
                if(prescan_only){
                    // if just counting viable histograms, print details, increase count in this file and continue
                    std::cout<<"Try to fit histo "<<thehist->GetName()<<" in file "<<filename<<std::endl;
                    passingHistos<<"Try to fit histo "<<thehist->GetName()<<" in file "<<filename<<std::endl;
                    continue;
                }
            }
            
            // having updated our count of viable histograms, check if we've hit our start point yet
            if(precut){
                // in this case offset & histos_to_analyse only count histograms after they've passed PeakScan
                if(offsetcount<offset){
                    ++offsetcount;
                    continue;
                }
            }
            
            // if we *didn't* pass both our pre-flight checks, but we *are* supposed to analyse this histo,
            // write the histo and TTree entry to file before moving on
            if(not found_spe_peak){
                // save histograms we skip, so we have a complete set and we can see what we missed
                outfile->cd();
                thehist->Write(thehist->GetName());
                
                // write an appropriate entry to the ROOT tree
                file_detectorkey = detkey;
                file_just_gaussian = just_gaussian;
                file_found_spe_peak = found_spe_peak;
                file_fit_success = false;
                file_prescaling = 0;
                file_ped_scaling = 0;
                file_ped_mean = 0;
                file_ped_sigma = 0;
                file_spe_firstgamma_scaling = 0;
                file_spe_firstgamma_mean = 0;
                file_spe_firstgamma_shape = 0;
                file_spe_secondgamma_scaling = 0;
                file_spe_secondgamma_mean_scaling = 0;
                file_spe_secondgamma_shape_scaling = 0;
                file_spe_expl_scaling = 0;
                file_spe_expl_charge_scaling = 0;
                file_mean_npe = 0;
                file_max_pes = 0;
                file_mean_spe_charge = 0;
                file_gain = 0;
                file_spe_firstgamma_gain = 0;
                
                outtree->Fill();
                outtree->Write("fit_parameters",TObject::kOverwrite);
                
                continue; // sorry, if we can't find a second bump, i can't generate suitable priors
            }
            
            //only scale the axes AFTER calling PeakScan (PeakScan should be count-independent: FIXME)
            double xscaling = 1;
            if(filesourcestring=="ToolAnalysis"){
                // for these files the histograms are in nC
                xscaling = 1000;
            } else if(filesourcestring=="TestStand"){
                // for test stand results the distributions are in electron counts
                xscaling = ELECTRON_CHARGE_IN_PICOCOULOMBS;
            }
            deapfitter.ScaleHistoXrange(xscaling);     // fix X-axis to picoCoulombs
            deapfitter.ScaleHistoYrange();             // normalise bin counts to unit integral
            deapfitter.SetRanges();                    // update internal ranges
            deapfitter.GeneratePriors();               // generate a "suitable" set of priors (needs improvement)
                                                       // Also sets the internal parameter limits, FYI
            
            /*
            // TODO derive this from occupancy
            deapfitter.FixParameter("mean_npe", 0.1);
            
            // TODO derive this from pedestal-only waveform fits
            deapfitter.FixParameter("ped_scaling", ped_scaling);
            deapfitter.FixParameter("ped_mean", ped_mean);
            deapfitter.FixParameter("ped_sigma", ped_sigma);
            */
            
#ifdef STORE_PRIORS
            // retrieve priors and add to the file for debugging
            std::vector<double> fit_parameters = deapfitter.GetParameters();
            file_prescaling_prior = fit_parameters.at(0);
            file_ped_scaling_prior = fit_parameters.at(1);
            file_ped_mean_prior = fit_parameters.at(2);
            file_ped_sigma_prior = fit_parameters.at(3);
            file_spe_firstgamma_scaling_prior = fit_parameters.at(4);
            file_spe_firstgamma_mean_prior = fit_parameters.at(5);
            file_spe_firstgamma_shape_prior = fit_parameters.at(6);
            file_spe_secondgamma_scaling_prior = fit_parameters.at(7);
            file_spe_secondgamma_mean_scaling_prior = fit_parameters.at(8);
            file_spe_secondgamma_shape_scaling_prior = fit_parameters.at(9);
            file_spe_expl_scaling_prior = fit_parameters.at(10);
            file_spe_expl_charge_scaling_prior = fit_parameters.at(11);
            file_mean_npe_prior = fit_parameters.at(12);
            file_max_pes_prior = fit_parameters.at(13);
            
            // also note even more uesless parameters used in deriving the priors, for debugging
            file_maxpos1 = precheck_pars.at(0);
            file_max1 = precheck_pars.at(1);
            file_maxpos2 = precheck_pars.at(2);
            file_max2 = precheck_pars.at(3);
            file_interpos = precheck_pars.at(4);
            file_intermin = precheck_pars.at(5);
            file_peaktovalleymin = precheck_pars.at(6);
            file_yscaling = deapfitter.GetYscaling();
#endif
            
            // draw the histogram and the prior
            std::cout<<"Drawing prior"<<std::endl;
            c1->Clear();
            c1->cd();
            thehist->Draw();
            thehist->ResetBit(kCanDelete);  // ROOT, plz... stop deleting my stuff
            thehist->SetDirectory(0);
            deapfitter.GetFullFitFunction()->Draw("same");
            gPad->SetLogy();
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
//            // write the histo to file, with prior fit function, for debug
//            std::cout<<"Writing prior fit to file"<<std::endl;
//            outfile->cd();
//            thehist->Write(TString::Format("%s_prior",thehist->GetName()));
            
            // We may be able to speed up fitting by reducing tolerances or changing out integration strategy
            // from https://root-forum.cern.ch/t/speeding-up-fitting-to-a-landau-distribution/25140/2
            //ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
            //ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100);
            //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100);
            //ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.1);
            // The above didn't seem to work, using TVirtualFitter did, but only on gpvm (ROOT 6.18),
            // not my laptop (ROOT 6.06)
            // Limiting max iterations does seem to limit maximum fitting time,
            // although the default (5000) should not have resulted in job hangs...
            // If MaxIterations is too low, the fitter may break early, thinking it's reached tolerance!
            // The returned EDM is incorrect!
            // Maybe it overruns by a certain % to check it's found the global min?
            // Bear this in mind when balaning MaxIterations vs EDM:
            // MaxIterations needs a suitable buffer above what seems necesary!
            // Maybe lowering the EDM more could also help protect against terminating at false minima...
            TVirtualFitter::SetMaxIterations(3000); // 2000 seems sufficient at a brief scan, allow more
            std::cout<<"Using at most "<<TVirtualFitter::GetMaxIterations()<<" fitting iterations"<<std::endl;
            TVirtualFitter::SetPrecision(0.1);  // 0.001 default, 20 seems sufficient...mostly. bad histograms fit badly.
            std::cout<<"Using a fit tolerance of "<<TVirtualFitter::GetPrecision()<<std::endl;
            // dunno if it helps, but it doesn't seem to hurt
            // Minuit2 is thread-safe so we could potentially fit multiple histos @ once...
            ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
            gSystem->Load("libMinuit2");
            
            // Try to do the fit
            std::cout<<"Fitting"<<std::endl;
            auto tstart = std::chrono::high_resolution_clock::now();
            TFitResultPtr fit_result = deapfitter.FitTheHisto("S");
            int fit_success = int(fit_result);
            double fit_chi2 = fit_result->Chi2();
            auto tend = std::chrono::high_resolution_clock::now();
            double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();
            time_taken *= 1e-9; // nano seconds to seconds
            std::cout<<"Fitting took "<<time_taken<<" sec"<<std::endl;
            std::cout<<"Fitting required: "<<fit_result->NCalls()<<" function calls"<<std::endl;
            std::cout<<"Fit returned: "<<fit_success<<" for histogram "<<thehist->GetName()<<std::endl;
            std::cout<<"Fit chi2 was: "<<fit_chi2<<std::endl;
            
            // get the fit parameters
            fit_parameters = deapfitter.GetParameters();
            std::cout<<"Fit parameters were: {";
            for(int pari=0; pari<fit_parameters.size(); ++pari){
                std::cout<<fit_parameters.at(pari);
                if((pari+1)<fit_parameters.size()) std::cout<<", ";
            }
            std::cout<<"}"<<std::endl;
            
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
            outfile->cd();
            thehist->Write(thehist->GetName());
            
            // copy it all into our ROOT file
            file_detectorkey = detkey;
            file_just_gaussian = just_gaussian;
            file_found_spe_peak = found_spe_peak;
            file_fit_success = fit_success;
            file_fit_chi2 = fit_chi2;
            file_prescaling = fit_parameters.at(0);
            file_ped_scaling = fit_parameters.at(1);
            file_ped_mean = fit_parameters.at(2);
            file_ped_sigma = fit_parameters.at(3);
            file_spe_firstgamma_scaling = fit_parameters.at(4);
            file_spe_firstgamma_mean = fit_parameters.at(5);
            file_spe_firstgamma_shape = fit_parameters.at(6);
            file_spe_secondgamma_scaling = fit_parameters.at(7);
            file_spe_secondgamma_mean_scaling = fit_parameters.at(8);
            file_spe_secondgamma_shape_scaling = fit_parameters.at(9);
            file_spe_expl_scaling = fit_parameters.at(10);
            file_spe_expl_charge_scaling = fit_parameters.at(11);
            file_mean_npe = fit_parameters.at(12);
            file_max_pes = fit_parameters.at(13);
            
            // Extract gain
            file_mean_spe_charge = deapfitter.GetMeanSPECharge();
            file_gain = (file_mean_spe_charge-file_ped_mean)/ELECTRON_CHARGE_IN_PICOCOULOMBS;
            file_spe_firstgamma_gain = (file_spe_firstgamma_mean-file_ped_mean)/ELECTRON_CHARGE_IN_PICOCOULOMBS;
            
            outtree->Fill();
            outtree->Write("fit_parameters",TObject::kOverwrite);
            
            // increment the counter of how many histograms we've fitted, and break if this is to be our last
            analysed_histo_count++;
            if(((analysed_histo_count+1)>histos_to_analyse)&&(histos_to_analyse>0)) break;
            
            //break; // XXX DEBUG so we can check
            
        } // end loop over histograms in file
        
        std::cout<<"We had "<<histos_with_a_peak<<" histograms that passed PeakScan in file "<<filename<<std::endl;
        
        // Cleanup
        if(infile){
            std::cout<<"Closing input file"<<std::endl;
            infile->Close();
            delete infile;
            infile = nullptr;
        }
        
    } // end loop over input files
    
    // close ROOT file if we had one
    if(outfile){
        std::cout<<"Closing output file"<<std::endl;
        //outfile->Write("",TObject::kOverwrite);
        outtree->ResetBranchAddresses();
        outfile->Close();
        delete outfile;
        outfile = nullptr;
    }
    // close prescan output file if we had one
    if(passingHistos.is_open()){
        passingHistos.close();
    }
    
    std::cout<<"Deleting canvas"<<std::endl;
    if(gROOT->FindObject("c1")!=nullptr) delete c1;
    
// gROOT maintains a list of all TF1s, and WILL DELETE THEM when it terminates.
// If a DEAPFitFunction object creates some functions, and then tries to delete them when
// it goes out of scope, it will segfault if the TApplication got there first.
// So either you rely on your user ONLY creating DEAPFitFunction objects on the heap
// so that you can call their destructor before terminating the TApplication,
// or you let DEAPFitFunction leak memory and rely on there being a TApplication
// to clean up later (i don't see why there necessarily will be one...)
//    std::cout<<"Deleting TApplication"<<std::endl;
//    if(FitApp) delete FitApp;   ...for now at least commenting this out works, though i don't know why, or if it will keep working
    std::cout<<"done, goodbye"<<std::endl;
}

