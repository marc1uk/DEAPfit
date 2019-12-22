#include "DEAPFitFunction.h"

DEAPFitFunction deapfitter();

// Make the fit function:
deapfitter.MakeFullFitFunction();  // Do not delete the returned TF1, it is owned by the DEAPFitFunction object!
// or
TF1* full_fit_func = new TF1("full_fit_func",&deapfitter,&DEAPFitFunction::FullFitFunction,0,100,14)
// or if you want to fit just components of the full fitting function
//TF1* full_fit_func = new TF1("full_fit_func",&deapfitter,&DEAPFitFunction::Gamma,0,100,14)
//TF1* full_fit_func = new TF1("full_fit_func",&deapfitter,&DEAPFitFunction::SPE_Func,0,100,14)
deapfitter.NameParameters(full_fit_func);   // only if you're using the FullFitFunction
full_fit_func->FixParameter("max_pes", 3);

// pass the histogram to the class so we can use it to generate priors
deapfitter.SetHisto(thehist);
deapfitter.SmoothHisto();  // this seems to help improve the fit

// Scan for a second peak after pedestal
bool found_spe_peak = deapfitter.PeakScan();
if(not found_spe_peak) continue; // sorry, if we can't find a second bump, i can't generate suitable priors

//only scale the axes AFTER calling PeakScan FIXME
deapfitter.ScaleHistoXrange(nano_to_pico);  // fix X-axis to picoCoulombs
deapfitter.ScaleHistoYrange();              // normalise bin counts to unit integral
deapfitter.SetRanges();                     // update internal ranges
deapfitter.GeneratePriors();                // generate a "suitable" set of priors (needs improvement)

// set the fit range (maybe redundant?)
full_fit_func->SetRange(thehist->GetXaxis()->GetXmin(),thehist->GetXaxis()->GetXmax());

// set any priors you have more info about than i do:
// if you have mean npe from occupancy:
deapfitter.FixMeanNpe(0.1);
// or
full_fit_func->FixParameter("mean_npe", 0.1);

// if you have pedestal from pedestal-only waveform fits
deapfitter.FixPedScaling(ped_scaling);
deapfitter.FixPedMean(ped_mean);
deapfitter.FixPedSigma(ped_sigma);
// or
full_fit_func->FixParameter("ped_scaling", ped_scaling);
full_fit_func->FixParameter("ped_mean", ped_mean);
full_fit_func->FixParameter("ped_sigma", ped_sigma);

int fit_success = thehist->Fit(full_fit_func);
// or
int fit_success = deapfitter.FitTheHisto();
std::cout<<"Fit returned: "<<fit_success<<std::endl;

// get parameters
std::vector<double> fitpars(full_fit_func->GetNpar());
full_fit_func->GetParameters(fitpars.data());
// or
std::vector<double> fitpars = deapfitter.GetParameters();

// clear the canvas
c1->Clear();
c1->cd();

// draw the histogram
thehist->Draw();
gPad->SetLogy();
// draw the fit result
full_fit_func->Draw("same");
// or
deapfitter.GetFunction()->Draw("same");

gSystem->ProcessEvents();


