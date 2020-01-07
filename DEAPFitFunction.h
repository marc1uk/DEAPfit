/* vim:set noexpandtab tabstop=4 wrap */

#include <string>
#include <vector>
#include <utility>
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
#include "Math/MinimizerOptions.h"

class DEAPFitFunction {
	
	public:
	// constructor & destructor
	DEAPFitFunction(int max_pes);
	~DEAPFitFunction();
	
	// misc setters
	void SetHisto(TH1* histo, TH1* bg_histo=nullptr);  // the histogram(s) to fit
	void SetBackgroundHisto(TH1* histo);               // background-only histogram
	// ranges
	void SetPedestalRange(double min=0, double max=0);
	void SetSpeRange(double min=0, double max=0);
	void SetNpeRange(double min=0, double max=0);
	std::pair<double,double> SetRanges(double min=0, double max=0);
	int LoadRanges();                                  // internally used, retrieve min-max of histogram
	
	// pre-processing functions
	TH1* ScaleHistoXrange(double inscaling);
	TH1* ScaleHistoYrange(double inscaling=0);
	TH1* SmoothHisto(int smoothing=1);
	int FitPedestal(std::vector<double>* ped_pars=nullptr);       // fit pedestal-only histogram
	bool PeakScan(std::vector<double>* precheck_pars=nullptr);    // search for a peak for determining priors
	double EstimateMeanNPE();                                     // estimate occupancy
	bool EstimateSPE(std::vector<double>* specheck_pars=nullptr); // fermi method to estimate posn of SPE peak
	bool GeneratePriors(std::vector<double>* fit_pars=nullptr,
						std::vector<std::pair<double,double>>* par_limits=nullptr, bool force_fit=false);
	
	// parameter setters
	int SetParameter(std::string par_name, double par_value);
	int SetParameter(int par_i, double par_value);
	void SetPrescaling(const double &prescaling_in);
	void SetPedScaling(const double &ped_scaling_in);
	void SetPedMean(const double &ped_mean_in);
	void SetPedSigma(const double &ped_sigma_in);
	void SetSpeFirstGammaScaling(const double &spe_firstgamma_scaling_in);
	void SetSpeFirstGammaMean(const double &spe_firstgamma_mean_in);
	void SetSpeFirstGammaShape(const double &spe_firstgamma_shape_in);
	void SetSpeSecondGammaScaling(const double &spe_secondgamma_scaling_in);
	void SetSpeSecondGammaMeanScaling(const double &spe_secondgamma_mean_scaling_in);
	void SetSpeSecondGammaShapeScaling(const double &spe_secondgamma_shape_scaling_in);
	void SetSpeExplScaling(const double &spe_expl_scaling_in);
	void SetSpeExplChargeScaling(const double &spe_expl_charge_scaling_in);
	void SetMeanNpe(const double &mean_npe_in);
	void SetMaxNpe(const double &max_pes_in);
	// group parameter setters
	int SetParameters(std::vector<double> fit_parameters);  // returns 1 on success
	void SetParameters(double* fit_parameters);
	int SetParameterLimits(std::vector<std::pair<double,double>> ranges_in);
	int SetParameterLimits(int par_i, std::pair<double,double> range_in);
	void RefreshParameters();  // pass internal member parameters to internal functions
	
	// parameter fixers
	void FixPrescaling(const double &prescaling_in);
	void FixPedScaling(const double &ped_scaling_in);
	void FixPedMean(const double &ped_mean_in);
	void FixPedSigma(const double &ped_sigma_in);
	void FixSpeFirstGammaScaling(const double &spe_firstgamma_scaling_in);
	void FixSpeFirstGammaMean(const double &spe_firstgamma_mean_in);
	void FixSpeFirstGammaShape(const double &spe_firstgamma_shape_in);
	void FixSpeSecondGammaScaling(const double &spe_secondgamma_scaling_in);
	void FixSpeSecondGammaMeanScaling(const double &spe_secondgamma_mean_scaling_in);
	void FixSpeSecondGammaShapeScaling(const double &spe_secondgamma_shape_scaling_in);
	void FixSpeExplScaling(const double &spe_expl_scaling_in);
	void FixSpeExplChargeScaling(const double &spe_expl_charge_scaling_in);
	void FixMeanNpe(const double &mean_npe_in);
	void FixMaxNpe(const double &max_pes_in);
	void FixingWarning();
	
	// getters
	std::vector<double> GetParameters();
	double GetParameter(std::string parname);
	std::vector<std::pair<double,double>> GetParameterLimits();
	std::pair<double,double> GetParameterLimits(std::string par_name);
	std::pair<double,double> GetParameterLimits(int par_i);
	double GetXscaling();
	double GetYscaling();
	
	// mathematical functions used in the fits
	double Gamma(double* x, double* gamma_pars);
	double SPE_Func(double* x, double* SPE_pars);
	double FullFitFunction(double* x, double* all_pars);
	double operator()(double *x, double *p);      // for using with a TF1 directly, calls FullFitFunction
	
	// retrieving the internals: Do not delete the returned objects!
	TF1* GetFullFitFunction();
	TF1* GetPedFunc();
	TF1* GetSPEFunc();
	TF1* GetBackgroundFunc();
	std::vector<TF1*> GetNPEFuncs();
	std::vector<TF1Convolution*> GetNPEConvs();
	std::vector<double> GetNPEPars();
	
	// Misc additional functions: for extracting gain from the fit
	double GetMeanSPECharge();
	
	// Operators on the TF1
	TF1* NameParameters(TF1* thefunc=nullptr);       // set names of parameters
	std::string ParameterNumberToName(int par_i);    // misc
	int ParameterNameToNumber(std::string par_name); // misc
	void ConstructFunctions();                       // construct the TF1s needed by the FullFitFunction
	TF1* MakeFullFitFunction();                      // construct TF1 internally. Do not delete the returned object!
	TFitResultPtr FitTheHisto(std::string opts="");  // call TH1::Fit with internal FullFitFunction
	
	private:
	// misc variables
	TH1* thehist=nullptr;                         // the histogram being worked on
	TH1* bghist=nullptr;                          // background-only histogram for pedestal estimation
	TF1* pedestal_standalone=nullptr;             // holds pedestal fit of background-only histogram, a "gaus"
	TF1* pedestal_func=nullptr;                   // holds the shape of the pedestal, a TMath::Gaus
	TF1* spe_func=nullptr;                        // holds the shape of the SPE peak
	TF1* full_fit_func=nullptr;                   // we may do the entire fitting internally
	std::vector<double> NPE_pars;                 // buffer to hold parameters as required by NPE peak TF1s
	std::vector<TF1*> npe_funcs;                  // functions used to fit the NPE peaks
	std::vector<TF1Convolution*> npe_convolns;    // used to obtain the TF1s
	int fit_success=-1;                           // result of the fit
	double xscaling=1;                            // keep track of it, just in case the user wants it back
	double yscaling=1;                            // needed to adjust a value from PreCheck for use in a prior.
	
	// static parameters involved in the fitting process
	// use valid ranges to prevent warnings on creation, they get set later
	double histogram_maximum=100;                 // ranges of the internal histogram we're fitting
	double histogram_minimum=0;                   // 
	double ped_range_max=100;                     // upper range of the pedestal gaussian function
	double ped_range_min=0;                       // lower range of the pedestal gaussian function
	double spe_range_max=100;                     // upper range of the SPE function
	double spe_range_min=0;                       // lower range of the SPE function
	double npe_func_max=100;                      // upper range of NPE convoluted functions
	double npe_func_min=0;                        // lower range of NPE convoluted functions
	int max_pes=3;                                // we'd need to use simulated annealing to fit this
	int n_spe_pars=8;                             // 
	
	// parameters used in PreCheck scan for some second peak beyond pedestal
	int maxpos1=-1;                               // position of first peak found (pedestal)
	int maxpos2=-1;                               // position of 2nd peak found (SPE)
	int interpos=-1;                              // position of valley
	int max1=-1;                                  // counts in pedestal
	int max2=-1;                                  // counts in SPE peak
	int intermin=-1;                              // counts in valley
	int peaktovalleymin=5;                        // peak must be this much above minima to reject noise XXX tuneme
	bool spe_peak_found=false;
	
	// parameters used to keep track of what we've called, which affects how we generate priors
	bool ran_fit_pedestal=false;
	bool ran_estimate_spe=false;
	bool ran_estimate_mean_npe=false;
	
	// fit parameters
	double prescaling=-1;
	double ped_scaling=-1;
	double ped_mean=-1;
	double ped_sigma=-1;
	double spe_firstgamma_scaling=-1;
	double spe_firstgamma_mean=-1;
	double spe_firstgamma_shape=-1;
	double spe_secondgamma_scaling=-1;
	double spe_secondgamma_mean_scaling=-1;
	double spe_secondgamma_shape_scaling=-1;
	double spe_expl_scaling=-1;
	double spe_expl_charge_scaling=-1;
	double mean_npe=-1;
	
	// other derived numbers
	double mean_spe_charge=-1;                   // estimated by EstimateSPE, or from fit by GetMeanSPECharge
	double spe_charge_variance=-1;               // estimated by EstimateSPE
	
	// fit parameter ranges
	std::vector<std::pair<double,double>> fit_parameter_ranges;
	
};
