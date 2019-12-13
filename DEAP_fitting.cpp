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

#include "DEAPFitFunction.h"

#ifndef STORE_PRIORS
#define STORE_PRIORS 1
#endif

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
    // also store all our priors for debug
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
    
    // also store some of the parameters used in determining priors, in case this is useful
    int file_maxpos1;
    int file_maxpos2;
    int file_interpos;
    int file_max1;
    int file_max2;
    int file_peaktovalleymin;
    int file_intermin;
    bool file_spe_peak_found;
    double file_yscaling;
    TBranch* bmaxpos1 = outtree->Branch("maxpos1", &file_maxpos1);
    TBranch* bmaxpos2 = outtree->Branch("maxpos2", &file_maxpos2);
    TBranch* binterpos = outtree->Branch("interpos", &file_interpos);
    TBranch* bmax1 = outtree->Branch("max1", &file_max1);
    TBranch* bmax2 = outtree->Branch("max2", &file_max2);
    TBranch* bpeaktovalleymin = outtree->Branch("peaktovalleymin", &file_peaktovalleymin);
    TBranch* bintermin = outtree->Branch("intermin", &file_intermin);
    TBranch* bspe_peak_found = outtree->Branch("spe_peak_found", &file_spe_peak_found);
    TBranch* byscaling = outtree->Branch("scaling", &file_yscaling);
#endif
    // OK, file made
    
    // =========================================
    // Initialization Complete, Move to Fitting
    // =========================================
    
    // define our fit function:
    int max_pes = 3;
    DEAPFitFunction deapfitter(max_pes);
    //TTF1* full_fit_func = new TF1("full_fit_func",FullFitFunction,0,100,14);  // you can do it this way, but it's more effort
    deapfitter.MakeFullFitFunction();
    
    // Loop over histos to fit
    int analysed_histo_count=0;
    int histos_to_analyse = (argc<5) ? 0 : atoi(argv[4]);
    int offset = (argc<4) ? 0 : atoi(argv[3]);
    int offsetcount=0;
    for(auto&& ahisto : histos_to_fit){
        // skip until requested histo
        if(offsetcount<offset){
            ++offsetcount;
            continue;
        }
        analysed_histo_count++;
        if((analysed_histo_count>histos_to_analyse)&&(histos_to_analyse>0)) break;
        // get next histogram
        int detkey = ahisto.first;
        thehist = ahisto.second;
        if(thehist==nullptr){ std::cerr<<"No hist"<<std::endl; return 0; }
        
        // pass to the fitter
        deapfitter.SetHisto(thehist);
        deapfitter.SmoothHisto();  // This seems to help the fit hit the broad SPE position better
        
        // Scan for a second peak after pedestal
        std::vector<double> precheck_pars;
        bool found_spe_peak = deapfitter.PeakScan(&precheck_pars);
        if(not found_spe_peak) continue; // sorry, if we can't find a second bump, i can't generate suitable priors
        
        //only scale the axes AFTER calling PeakScan FIXME
        double xscaling = 1;
        if(filesourcestring=="ToolAnalysis"){
            // for these files the histograms are in nC
            xscaling = 1000;
        } else if(filesourcestring=="TestStand"){
            // for test stand results the distributions are in electron counts
            double ELECTRON_CHARGE = 1.6E-19;
            double ELECTRON_CHARGE_IN_PICOCOULOMBS = ELECTRON_CHARGE*1E12;
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
        file_spe_peak_found = found_spe_peak;
        file_yscaling = deapfitter.GetYscaling();
#endif
        
        // draw the histogram and the prior
        std::cout<<"Drawing prior"<<std::endl;
        c1->Clear();
        c1->cd();
        thehist->Draw();
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
        int fit_success = deapfitter.FitTheHisto();
        auto tend = std::chrono::high_resolution_clock::now();
        double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(tend-tstart).count();
        time_taken *= 1e-9; // nano seconds to seconds
        std::cout<<"Fitting took "<<time_taken<<" sec"<<std::endl;
        std::cout<<"Fit returned: "<<fit_success<<" for histogram "<<thehist->GetName()<<std::endl;
        
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
        thehist->Write(thehist->GetName());
        
        // copy it all into our ROOT file
        file_detectorkey = detkey;
        file_fit_success = fit_success;
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

