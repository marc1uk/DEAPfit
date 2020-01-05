/* vim:set expandtab tabstop=4 wrap */

#include <string>
#include <thread>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iostream>
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

int main(int argc, const char* argv[]){
    if(argc<3){
        std::cout<<"Usage: "<<argv[0]<<" <input file> <output file>"<<std::endl;
        return 0;
    }
    std::string filename = argv[1];
    std::string outputfilename = argv[2];
    
    std::map<std::string, TGraph*> tgraphs_to_fit;
    TGraph* thehist=nullptr;
    
    // Open root file
    TFile* infile = TFile::Open(filename.c_str(),"READ");
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
        // skip unless it's a TGraph
        if (!cl->InheritsFrom("TGraph")) continue;
        //std::cout<<"it's a tgraph"<<std::endl;
        // graph names are "HV_vs_Gain_XXX" where xxx is the ToolAnalysis detectorkey for that PMT.
        std::string histname = key->GetName();
        int detectorkey;
        int n = 0;
        int cnt = sscanf(histname.c_str(),"HV_vs_Gain_%d %n",&detectorkey, &n);
        int success = ((n >0) && ((histname.c_str())[n]=='\0'));
        if(success){
            std::cout<<"Found tgraph "<<histname<<std::endl;
            // read the histogram from the file into memory(?) and get a pointer
            thehist = (TGraph*)key->ReadObj();  // should we call thehist->Delete() when we're done with it?
            tgraphs_to_fit.emplace(histname,thehist);
        }
    } // end loop over keys in file
    
    std::vector<TCanvas*> canvases;
    int n_pages=0;
    int summarypdf=0; // whether to do 12 plots (T) per page or 1 plot per page (0)
    int histcounter=-1;
    // loop over channels
    for(auto&& adetector : tgraphs_to_fit){
        ++histcounter;
        std::string histname = adetector.first;
        TGraph* thegraph = adetector.second;
        
        if(summarypdf==1){
            // move to a new canvas each 12 histos
            TString cname = "page_";
            cname+=n_pages;
            if(histcounter%12==0){
              canvases.push_back(new TCanvas(cname,cname,1400,800));
              canvases.back()->Divide(4,3);
              n_pages++;
            }
            canvases.back()->cd((histcounter%12)+1);
        } else {
            TString cname = "page_";
            cname+=n_pages;
            canvases.push_back(new TCanvas(cname,cname,1400,800));
            canvases.back()->cd();
            n_pages++;
        }
        
        // draw the histo
        thegraph->Draw("A*P");
        gPad->SetLogy();
    }
    
    for(int pagei=0; pagei<canvases.size(); pagei++){
      // combine all pages into a single file
      if(pagei==0){
        canvases.at(pagei)->Print((outputfilename+"(").c_str(),TString::Format("Title:%d",pagei));
      } else if((pagei+1)==canvases.size()){
        canvases.at(pagei)->Print((outputfilename+")").c_str(),TString::Format("Title:%d",pagei));
      } else{
        canvases.at(pagei)->Print(outputfilename.c_str(),TString::Format("Title:%d",pagei));
      }
    } // end loop over writing to pdf
    
    infile->Close();
    delete infile; infile=nullptr;
}
