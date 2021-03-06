#include <string>
#include <algorithm>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TText.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TKey.h"

int main(int argc, const char* argv[]){
    // make the TApplication for drawing
    //int myargc = 1;
    //char arg1[] = "potato";
    //char* myargv[]={arg1};
    //TApplication *PlotApp = new TApplication("PlotApp",&myargc,myargv);
    if(argc<2){
      std::cout<<"Usage: "<<argv[0]<<" <input_base> <num_files>"<<std::endl;
      std::cout<<"<input_base> should be such that input files have names '<input_base>_*.root'"
               <<", where * runs from 0 to <num_files>-1"<<std::endl;
      return 0;
    }
    
    // XXX we we should take these by argument or something...
    std::vector<int> OffTubes{332, 338, 342, 344, 345, 346, 349, 352, 359, 394, 431, 444, 445};

    //std::string filenamebase = "LEDRun1178S0LEDs5And6And10And28_PulseWindowOnly_PMTStability_Run0";
    std::string filenamebase = "DEAPfitter_outfile";
    if(argc>1) filenamebase = argv[1];
    if(  (filenamebase.find_last_of(".") != std::string::npos) &&
         (! ((filenamebase.find_last_of(".")==1)&&(filenamebase.find_first_of('.')==0)) )){
        filenamebase = filenamebase.substr(0, filenamebase.find_last_of('.'));
    }
    std::string outputfilename = filenamebase+"_fitted.pdf";
    int num_jobs = 1; // how many files did we split the results into - how many grid jobs processed the fits?
    if(argc>2) num_jobs = atoi(argv[2]);
    // loop over jobs
    std::map<std::string,TH1*> thehistos;
    // this happens when there were no histograms in it's analysed set that had an SPE peak
    for(int jobi=0; jobi<num_jobs; ++jobi){
        // load it's output file
        std::string jobfilename = filenamebase+"_"+std::to_string(jobi)+".root";
        TFile* infile = new TFile(jobfilename.c_str(),"READ"); // we leak these files...
        if(infile==nullptr){
            std::cout<<"Could not find output file "<<filenamebase<<", skipping"<<std::endl;
        }
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
            int voltage;
            int n = 0;
            int cnt = sscanf(histname.c_str(),"hist_charge_%d_%dV %n",&detectorkey, &voltage, &n);
            int success = ((n >0) && ((histname.c_str())[n]=='\0'));
            if(success){
                std::cout<<"Found charge histogram "<<histname<<std::endl;
                std::vector<int>::iterator it = std::find(OffTubes.begin(),OffTubes.end(),detectorkey);
                if( it != OffTubes.end()) continue; // channel is turned off, skip...
                
                // read the histogram from the file into memory(?) and get a pointer
                TH1* thehist = (TH1*)key->ReadObj();
                thehist->SetName(histname.c_str());
                thehist->SetTitle(TString::Format("detkey %d at %dV",detectorkey,voltage));
                thehistos.emplace(histname,thehist);
            }
        } // end loop over keys in file
    } // end loop over results files
    std::cout<<"Found "<<thehistos.size()<<" histograms to plot"<<std::endl;
    
    std::vector<TCanvas*> canvases;
    int n_pages=0;
    int histcounter=-1;
    // loop over channels
    for(auto&& apair : thehistos){
      histcounter++;
      TH1* thehist = apair.second;
      
      // move to a new canvas each 12 histos
      if(histcounter%12==0){
        TString cname = "page_";
        cname+=n_pages;
        canvases.push_back(new TCanvas(cname,cname,1400,800));
        canvases.back()->Divide(4,3);
        n_pages++;
      }
      
      // draw the histo
      canvases.back()->cd((histcounter%12)+1);
      gPad-> SetLogy();
      thehist->SetLineColor(4);
      thehist->Draw();
      
    }
    
    std::cout<<"writing "<<canvases.size()<<" pages"<<std::endl;
    for(int pagei=0; pagei<canvases.size(); pagei++){
      // write each page as a separate file
      //std::cout << "WE HERE " << std::endl;
      //TString ofname = "charge_dists_";
      //ofname+=pagei;
      //ofname+=".pdf";
      //canvases.at(pagei)->SaveAs(ofname);
      // 
      // combine all pages into a single file
      if(pagei==0){
        canvases.at(pagei)->Print((outputfilename+"(").c_str(),TString::Format("Title:%d",pagei));
      } else if((pagei+1)==canvases.size()){
        canvases.at(pagei)->Print((outputfilename+")").c_str(),TString::Format("Title:%d",pagei));
      } else {
        canvases.at(pagei)->Print(outputfilename.c_str(),TString::Format("Title:%d",pagei));
      }
    } // end loop over writing to pdf
    std::cout<<"Finished. Goodbye"<<std::endl;
    // so why is the application still not terminating after this...?
    // these don't help...
    gROOT->CloseFiles();
    gROOT->EndOfProcessCleanups();
    gROOT->ShutDown();
    if(gApplication) gApplication->Terminate(0);
    return 0;
    // maybe we can force it with
    std::cout<<"exiting"<<std::endl;
    exit(1);
}
