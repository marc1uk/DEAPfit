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
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TText.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TKey.h"

int main(int argc, const char* argv[]){
	
	if(argc<2){
		std::cout<<"Usage: "<<argv[0]<<" <input file pattern> <outfile name>"<<std::endl;
		return 0;
	}
	std::string inpattern = argv[1];
	
	std::string outfilename = "DEAPFitter_HV_vs_Gain.root"; // default
	if(argc>2){
		outfilename = argv[2];
	}
	
	TChain* inchain = new TChain("fit_parameters");
	
	inchain->SetBranchStatus("*",0); // disable all branches, we only need a few
	inchain->SetBranchStatus("DetectorKey",1);          // 
	inchain->SetBranchStatus("Voltage",1);              // 
	inchain->SetBranchStatus("found_spe_peak",1);       // did we try to fit it
	inchain->SetBranchStatus("fit_result",1);           // did we think we were successful?
	inchain->SetBranchStatus("fit_chi2",1);             // were we successful?
	inchain->SetBranchStatus("gain",1);                 // what was the alleged gain?
	inchain->SetBranchStatus("spe_firstgamma_gain",1);  // debug, for comparison
	inchain->SetBranchStatus("mean_npe",1);             // for reference, what was the alleged illumination
	
	// add our list of files
	std::cout<<"Adding files from pattern "<<inpattern<<std::endl;
	int ret = inchain->Add(inpattern.c_str());
	if((ret==0)||(inchain->GetListOfFiles()->GetEntries()==0)){
		std::cerr<<"Could not find any files matching pattern "<<inpattern<<std::endl;
		return 0;
	}
	std::cout<<"Added "<<ret<<" files, loading first tree"<<std::endl;
	
	// check they have the desired tree
	ret = inchain->LoadTree(0);
	if(ret<0){
		std::cerr<<"No matching tree in file "<<inchain->GetListOfFiles()->At(0)->GetTitle()<<std::endl;
		return 0;
	}
	std::cout<<"First tree has "<<inchain->GetTree()->GetEntries()<<" entries"<<std::endl;
	
	std::cout<<"Setting branch addresses"<<std::endl;
	// connect enabled branches
	int voltage=0;
	int detectorkey;
	int found_spe_peak;
	int fit_result;
	double fit_chi2;
	double mean_npe;
	double gain;
	double spe_firstgamma_gain; // for reference...
	inchain->SetBranchAddress("DetectorKey",&detectorkey);
	inchain->SetBranchAddress("Voltage",&voltage);
	inchain->SetBranchAddress("found_spe_peak",&found_spe_peak);
	inchain->SetBranchAddress("fit_result",&fit_result);
	inchain->SetBranchAddress("fit_chi2",&fit_chi2);
	inchain->SetBranchAddress("gain",&gain);
	inchain->SetBranchAddress("spe_firstgamma_gain",&spe_firstgamma_gain);
	inchain->SetBranchAddress("mean_npe",&mean_npe);
	
	// loop over the input file, building the set of gain arrays
	// note that the chain may not be in order, and not all channels may have the same voltages
	std::map<int,std::vector<std::pair<double,double>>> gain_data;
	
	int entryi=0;
	std::cout<<"Scanning TChain"<<std::endl;
	while(true){
		int localentry = inchain->LoadTree(entryi);
		if(localentry<0) break;     // ran off end of TChain
		inchain->GetEntry(entryi);  // get next entry
		++entryi;
		
		// only plot gain values with a successful fit FIXME we should note these, somehow...
		if(not found_spe_peak) continue;
		
		// TODO should we not plot points with an invalid fit_result? Sometimes fit's good even so...
		
		// TODO use chi2 of fit to add error bars to this datapoint...
		// need to think about how to convert chi2 to an error...
		// maybe we need to store the error on the fit parameter in DEAP_fitting.cpp
		// we'd need to store these in another map...
		
		// check if this is the first time we've seen this detkey
		std::cout<<"Adding voltage "<<voltage<<" gain "<<gain<<" to detectorkey "<<detectorkey<<std::endl;
		if(gain_data.count(detectorkey)==0){
			// add a new entry to the map
			gain_data.emplace(detectorkey,std::vector<std::pair<double,double>>{std::pair<double,double>{voltage,gain}});
		} else {
			// add the new data point to the vectors
			gain_data.at(detectorkey).push_back(std::pair<double,double>{voltage,gain});
		}
	}
	std::cout<<"Scanned "<<entryi<<" entries in TChain"<<std::endl;
	std::cout<<"Found "<<gain_data.size()<<" PMTs to make plots for"<<std::endl;
	
	// close the input file
	inchain->ResetBranchAddresses();
	delete inchain; inchain=nullptr;
	
	// make an output file in which to save results
	TFile* outfile = new TFile(outfilename.c_str(),"RECREATE");
	outfile->cd();
	
	//make the TGraphs XXX TODO TGraphErrors once we have the error
	TGraph* thegraph = new TGraph();
	thegraph->GetXaxis()->SetTitle("Voltage [V]");
	thegraph->GetYaxis()->SetTitle("Gain [no units]");
	for(auto&& pmt_dataset : gain_data){
		int detectorkey = pmt_dataset.first;
		 // get N points in TGraph from size of points vector
		int npoints = pmt_dataset.second.size();
		std::cout<<"Making plot of "<<npoints<<" points"<<std::endl;
		// update size of TGraph
		thegraph->Set(npoints);
		// update points in TGraph
		for(int pointi=0; pointi<npoints; ++pointi){
			thegraph->SetPoint(pointi, pmt_dataset.second.at(pointi).first, pmt_dataset.second.at(pointi).second);
		}
		// update axes
		double y_lim_lo = *std::min_element(thegraph->GetY(),thegraph->GetY()+thegraph->GetN());
		double y_lim_hi = *std::max_element(thegraph->GetY(),thegraph->GetY()+thegraph->GetN());
		thegraph->SetMinimum(y_lim_lo*0.9);
		thegraph->SetMaximum(y_lim_hi*1.1);
		double x_lim_lo = *std::min_element(thegraph->GetX(),thegraph->GetX()+thegraph->GetN());
		double x_lim_hi = *std::max_element(thegraph->GetX(),thegraph->GetX()+thegraph->GetN());
		thegraph->GetHistogram()->GetXaxis()->SetLimits(x_lim_lo*0.9, x_lim_hi*1.1);
		
		// set title
		std::string titlestring = "HV_vs_Gain_"+std::to_string(detectorkey);
		thegraph->SetName(titlestring.c_str());
		thegraph->SetTitle(titlestring.c_str());
		thegraph->Write(titlestring.c_str(),TObject::kOverwrite);
	}
	std::cout<<"TGraphs constructed, output written to "<<outfilename<<std::endl;
	
	// close output file
	outfile->Write(); // write header of TFile, necessary even if all objects had Write called... maybe?
	outfile->Close();
	delete outfile; outfile=nullptr;
	std::cout<<"Finished. Goodbye"<<std::endl;
}
