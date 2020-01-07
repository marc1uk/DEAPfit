#include "DEAPFitFunction.h"

DEAPFitFunction::DEAPFitFunction(int max_pes_in){
	// constructor
	max_pes = max_pes_in;
	ConstructFunctions();
	fit_parameter_ranges = std::vector<std::pair<double,double>>(14,std::pair<double,double>{0,0});
}

DEAPFitFunction::~DEAPFitFunction(){
	//std::cout<<"DEAPFitFunction destructor"<<std::endl;
	//std::cout<<"Deleting standalone pedestal function"<<std::endl;
	if(pedestal_standalone) delete pedestal_standalone; pedestal_standalone=nullptr;
	//std::cout<<"Deleting pedestal function"<<std::endl;
	if(pedestal_func) delete pedestal_func; pedestal_func=nullptr;
	//std::cout<<"Deleting spe function"<<std::endl;
	if(spe_func) delete spe_func; spe_func=nullptr;
	//std::cout<<"Deleting npe functions"<<std::endl;
	for(TF1* afunc : npe_funcs) delete afunc;  // do functions built from convolutions own them??? who knows
	npe_funcs.clear();
	//std::cout<<"Deleting convolutions"<<std::endl;
	for(TF1Convolution* aconv : npe_convolns) delete aconv;
	npe_convolns.clear();
	//std::cout<<"Deleting full fit function"<<std::endl;
	if(full_fit_func) delete full_fit_func; full_fit_func=nullptr;
}

void DEAPFitFunction::SetPedestalRange(double min, double max){
	// if no values passed (default) pull ranges from histogram
	if((min==0)&&(max==0)){
		int rangesok = LoadRanges();
		if(not rangesok) return;
		min = -histogram_maximum;
		max =  histogram_maximum;
	}
	ped_range_min = min;
	ped_range_max = max;
	if(pedestal_func) pedestal_func->SetRange(ped_range_min, ped_range_max);
}
void DEAPFitFunction::SetSpeRange(double min, double max){
	// if no values passed (default) pull ranges from histogram
	if((min==0)&&(max==0)){
		int rangesok = LoadRanges();
		if(not rangesok) return;
		min = -histogram_maximum;
		max =  histogram_maximum;
	}
	spe_range_min = min;
	spe_range_max = max;
	if(spe_func) spe_func->SetRange(spe_range_min,spe_range_max);
}
void DEAPFitFunction::SetNpeRange(double min, double max){
	// if no values passed (default) pull ranges from histogram
	double min_func=0;
	double min_conv=0;  // XXX FIXME the convolution range should probably be less than the input func ranges
	if((min==0)&&(max==0)){
		int rangesok = LoadRanges();
		if(not rangesok) return;
		min_conv = histogram_minimum;
		min_func = -histogram_maximum;
		min = min_func;
		max =  histogram_maximum;
	} else {
		min_func = min;
		min_conv = min;
	}
	npe_func_min = min;  // XXX we don't have enough parameters... i don't know how this works...
	npe_func_max = max;
	for(TF1Convolution* aconv : npe_convolns) aconv->SetRange(min_conv,npe_func_max);
	for(TF1* afunc : npe_funcs) afunc->SetRange(min_func,npe_func_max);
}
std::pair<double,double> DEAPFitFunction::SetRanges(double min, double max){
	// if no values passed (default) pull ranges from histogram
	if((min==0)&&(max==0)){
		int rangesok = LoadRanges();
		if(not rangesok) return std::pair<double,double>{0,0};
		min = -histogram_maximum;
		max =  histogram_maximum;
	}
	SetPedestalRange(min,max);
	SetSpeRange(min,max);
	SetNpeRange(min,max);
	if(full_fit_func) full_fit_func->SetRange(histogram_minimum,histogram_maximum);
	return std::pair<double,double>{min,max};
}
int DEAPFitFunction::LoadRanges(){
	if((thehist->GetXaxis()->GetXmax()>0)&&(thehist->GetXaxis()->GetXmax()>thehist->GetXaxis()->GetXmin())){
		// retrieve histogram range for default setting of fit function ranges
		histogram_minimum = thehist->GetXaxis()->GetXmin();  // THIS ALONE ISN'T ENOUGH
		histogram_maximum = thehist->GetXaxis()->GetXmax();  // convolutions will be CHOPPED
		// at the least we need to use -hist_max to hist_max for TF1 ranges
		return 1;
	} else {
		std::cerr<<"DEAPFitFunction::LoadRanges found invalid histogram range "
			 <<thehist->GetXaxis()->GetXmax()<<" to "<<thehist->GetXaxis()->GetXmin()<<std::endl;
		return 0;
	}
}

void DEAPFitFunction::SetHisto(TH1* histo, TH1* bg_histo){
	if(histo==nullptr){
		std::cerr<<"DEAPFitFunction::SetHisto called with nullptr!"<<std::endl;
		return;
	}
	thehist = histo;
	if(thehist->GetListOfFunctions()->GetSize()>0){
		//std::cerr<<"WARNING: DEAPFitFunction::SetHisto clears any functions "
		//	 <<"owned by input histogram!"<<std::endl;
		thehist->GetListOfFunctions()->Clear();
	}
	SetRanges();
	
	if(bg_histo) SetBackgroundHisto(bg_histo);
	else bghist = nullptr; // clear it so that we don't accidentally use it on the wrong histo
	
	// clear internal flags that record whether we've fit a background histogram
	// and whether we've used EstimateSPE. These are used in GeneratePriors.
	ran_fit_pedestal=false; ran_estimate_spe=false; ran_estimate_mean_npe=false;
}

void DEAPFitFunction::SetBackgroundHisto(TH1* histo){
	if(histo==nullptr){
		std::cerr<<"DEAPFitFunction::SetBackgroundHisto called with nullptr!"<<std::endl;
		return;
	}
	bghist = histo;
	if(bghist->GetListOfFunctions()->GetSize()>0){
		//std::cerr<<"WARNING: DEAPFitFunction::SetBackgroundHisto clears any functions "
		//	 <<"owned by input histogram!"<<std::endl;
		bghist->GetListOfFunctions()->Clear();
	}
}

TH1* DEAPFitFunction::ScaleHistoXrange(double inscaling){
	thehist->GetXaxis()->Set(thehist->GetNbinsX(),
				 thehist->GetXaxis()->GetXmin()*inscaling,
				 thehist->GetXaxis()->GetXmax()*inscaling);
	SetRanges();
	return thehist;
}
TH1* DEAPFitFunction::ScaleHistoYrange(double inscaling){
	if(inscaling==0){
		// by default scale to unit area
		//inscaling = thehist->Integral();
		// we could also consider unit peak?
		inscaling = thehist->GetBinContent(thehist->GetMaximumBin());
	}
	// keep track of it for calculating priors
	yscaling = inscaling;
	thehist->Scale(1./yscaling);
	return thehist;
}

TH1* DEAPFitFunction::SmoothHisto(int smoothing){
	//thehist->Smooth(smoothing);
	
	// ROOT's TH1::Smooth is fine when the histogram spans enough bins,
	// but for histograms with only a few populated bins it can add large peaks and valleys
	// to a previously smooth part of the histogram!
	// In extreme cases it even creates valleys that split the histogram such that
	// bins in a previously highly populated region are empty!!!
	// The problem may be related to poor handling of sharp gradients in histograms
	// dominated by a narrow gaussian
	// (this epiphany came from attempts to convert to a TGraph and use TGraphSmooth,
	//  which fit the tail well but then the errors would explode at the sharp pedestal rise)
	// 
	// SO: since we're fitting it in log form anyway, convert it to log, smooth it, then convert it back!
	TH1D smoothhist("smoothhist","smoothhist",
		thehist->GetNbinsX(),thehist->GetXaxis()->GetXmin(),thehist->GetXaxis()->GetXmax());
	for(int bini=0; bini<thehist->GetNbinsX(); bini++){
		double bincont = (thehist->GetBinContent(bini)==0) ? 0 : log(thehist->GetBinContent(bini));
		smoothhist.SetBinContent(bini,bincont);
	}
	smoothhist.Smooth(smoothing);
	
	// transfer smoothed bin contents back, accounting for our temporary log
	for(int bini=0; bini<thehist->GetNbinsX(); bini++){
		double bincont = (smoothhist.GetBinContent(bini)==0) ? 0 : exp(smoothhist.GetBinContent(bini));
		thehist->SetBinContent(bini,bincont);
	}
	
	return thehist;
}

bool DEAPFitFunction::PeakScan(std::vector<double>* precheck_pars){
	// FIXME dependant on the bin counts (doesn't work after we scale the histogram)
	// some very crude carry-over stuff from the triple-gaussian fit used for the PMT test stand.
	// designed to identify a few salient points in the distribution, for prior calculations
	// Try to find a second distinct peak after the pedestal
	
	maxpos1=-1, maxpos2=-1, interpos=-1;
	max1=-1, max2=-1;
	intermin=std::numeric_limits<int>::max();
	
	// ========================================================================
	// STEP 1: SCAN OVER HISTOGRAM AND TRY TO FIND A SECOND PEAK AFTER PEDESTAL
	// ========================================================================
	
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
			 // it's position is 'sufficiently far' from the valley:
			 // don't consider a SPE peak that's in an adjacent bin to the valley
			 ((bini-interpos)>3) &&
			 // don't consider situations where the 'valley' has 0 counts
			 // (typically we've run off the end of the distribution and hit noise
			 (intermin>0) &&
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
	std::cout<<"Found pedestal at "<<maxpos1<<" ("<<thehist->GetBinCenter(maxpos1)<<")"
		 <<" with counts "<<max1<<std::endl;
	std::cout<<"Found valley at "<<interpos<<" ("<<thehist->GetBinCenter(interpos)<<")"
		 <<" with counts "<<intermin<<std::endl;
	std::cout<<"Found a SPE peak at "<<maxpos2<<" ("<<thehist->GetBinCenter(maxpos2)<<")"
		<<" with counts "<<max2<<std::endl;
	spe_peak_found = ((max2>10)&&(maxpos2>0));
	
	// copy the intermediate pars in case user wants them
	if(precheck_pars!=nullptr){
		precheck_pars->clear();
		precheck_pars->push_back(maxpos1);
		precheck_pars->push_back(max1);
		precheck_pars->push_back(maxpos2);
		precheck_pars->push_back(max2);
		precheck_pars->push_back(interpos);
		precheck_pars->push_back(intermin);
		precheck_pars->push_back(peaktovalleymin);
	}
	
	return spe_peak_found;
}

int DEAPFitFunction::FitPedestal(std::vector<double>* ped_pars){
	if(bghist==nullptr){
		std::cerr<<"DEAPFitFunction::FidPedestal called with no background histogram set!"<<std::endl
			 <<"Call DEAPFitFunction::SetBackgroundHisto before this, or pass it "
			 <<"with DEAPFitFunction::SetHisto(TH1* signal_histo, TH1* bg_histo)"<<std::endl;
		return 0;
	}
	
	// fit a gaussian to the provided background-only histogram
	// pedestal_standalone is "gaus" == [0]*exp(-0.5*((x-[1])/[2])**2)
	int fit_success = bghist->Fit(pedestal_standalone,"QN"); // quiet, do not store
	
	if(fit_success){
		std::cerr<<"DEAPFitFunction::FidPedestal fit failed!"<<std::endl
			 <<"Pedestal parameters will not be set!"<<std::endl;
		return fit_success;
	}
	
	// retrieve the fit parameters
	// we don't take the amplitude because the counts in background and signal histos may be different
	ped_mean = pedestal_standalone->GetParameter(1);
	ped_sigma = pedestal_standalone->GetParameter(2);
	
	// pass to the internal pedestal function to keep in sync
	if(pedestal_func){
		pedestal_func->SetParameter("ped_mean",ped_mean);
		pedestal_func->SetParameter("ped_sigma",ped_sigma);
	}
	
	// if we were given a vector, fill the components with the fit parameters
	*ped_pars = std::vector<double>{pedestal_standalone->GetParameter(0),ped_mean,ped_sigma};
	
	ran_fit_pedestal = true; // note that we ran this, for use in GeneratePriors
	return fit_success;
}

double DEAPFitFunction::EstimateMeanNPE(){
	// We can use a histogram of background-only events to estimate the mean_npe (aka occupancy)
	// We do this by placing a low-charge cut on the background-only distribution,
	// such that in our total distribution it would only capture background events.
	// We apply this cut to the blank distribution and calculate what fraction of background it contains,
	// in order to determine the efficiency of the cut.
	// We then apply the cut to the total distribution, and correct for the efficiency to estimate
	// the number of 0-pe events it contains.
	// Assuming a poisson distribution for the number of pe's in an event, the number of 0-pe events
	// can be used to determine the mean number of pes via:
	// mean_npe = -ln(n_0pe_triggers/n_total_triggers)
	
	if(thehist==nullptr){
		std::cerr<<"DEAPFitFunction::EstimateMeanNPE called but we have no histogram!"
			 <<" Call SetHisto first"<<std::endl;
		return -1;
	}
	
	if(bghist==nullptr){
		std::cerr<<"DEAPFitFunction::EstimateMeanNPE requires a charge distribution"
			 <<"  of background-only events! Call SetBackgroundHisto first"<<std::endl;
		return -1;
	}
	
	// fit the pedestal if we haven't, and then use the fit to derive a cut-off charge
	if(not ran_fit_pedestal){
		int fpsucess = FitPedestal();
		if(fpsucess!=0) return -1;
	}
	double ped_cutoff = ped_mean + ped_sigma*5.; // teal uses 5, which seems like a lot...
	
	// first determine the cut efficiency with the background histogram
	// count num events below the cutoff
	double low_q_background_events=0;
	for(int bini=1; bini<bghist->GetNbinsX(); bini++){ // bin 0 is the underflow bin
		double bincenter = bghist->GetBinCenter(bini);
		if(bincenter>ped_cutoff) break;
		low_q_background_events += bghist->GetBinContent(bini);
	}
	// determine the cut efficiency
	double cut_efficiency = low_q_background_events/bghist->GetEntries();
	
	// count the number of events with charge below the cut threshold
	double n_0pe_triggers=0;
	for(int bini=1; bini<thehist->GetNbinsX(); bini++){ // bin 0 is the underflow bin
		double bincenter = thehist->GetBinCenter(bini);
		if(bincenter>ped_cutoff) break;
		n_0pe_triggers += thehist->GetBinContent(bini);
	}
	// scale it up by the efficiency
	n_0pe_triggers /= cut_efficiency;
	
	// get the occupancy
	mean_npe = -log(n_0pe_triggers/thehist->GetEntries());
	
	ran_estimate_mean_npe = true;
	return mean_npe;
}

bool DEAPFitFunction::EstimateSPE(std::vector<double>* specheck_pars){
	// ========================================================================
	// STEP 1b: USE FERMI METHOD TO ESTIMATE SPE PEAK POSITION
	// ========================================================================
	// mean_spe_charge = (mean_total_dist - mean_bg_dist) / mean_npe
	// where
	// * mean_total_dist is the mean of the charge distribution for signal+background
	// * mean_bg_dist is the mean of the charge distribution for background only
	// * mean_npe is the mean of the discrete probability distribution of the number of pe's in a flash
	
	if(not ran_estimate_mean_npe){
		double mnpeok = EstimateMeanNPE();
		if(mnpeok<0) return false;
	}
	
	// finally determine the mean spe charge
	double mean_total_dist = thehist->GetMean();
	double mean_bg_dist = bghist->GetMean();
	mean_spe_charge = (mean_total_dist - mean_bg_dist) / mean_npe;
	
	// to derive some estimate of the width of the SPE peak, we can calculate the spe variance
	// spe_charge_variance = (tot_dist_variance - bg_dist_variance)/mean_npe - mean_spe_charge^2
	double bg_dist_variance = pow(bghist->GetStdDev(),2.);
	double tot_dist_variance = pow(thehist->GetStdDev(),2.);
	spe_charge_variance = (tot_dist_variance - bg_dist_variance)/mean_npe - pow(mean_spe_charge,2.);
	
	// pass back to the user
	if(specheck_pars){
		// TODO it would make sense to pass these back as a map<name,value> so they have meaning
		*specheck_pars = std::vector<double>{mean_total_dist,
						     mean_bg_dist,
						     mean_npe,
						     mean_spe_charge,
						     bg_dist_variance,
						     tot_dist_variance,
						     spe_charge_variance};
	}
	
	ran_estimate_spe = true;  // note that we ran this, for use in GeneratePriors
	return true; // TODO error checking?
}

bool DEAPFitFunction::GeneratePriors(std::vector<double>* fit_pars, std::vector<std::pair<double,double>>* par_limits, bool force_fit){
	
	if(not spe_peak_found || force_fit){
		std::cerr<<"ERROR: Sorry, current generation of suitable priors relies on successful finding of "
			 <<" a second bump in the charge distribution. Please run PeakScan first, if you haven't."
			 <<"If you have, no suitable bump was found."<<std::endl
			 <<"To proceed, you'll need to generate your own priors "
			 <<"and pass them via DEAPFitFunction::SetParameters."
			 <<std::endl;
		return false;
	}
	
	// there are a few potential sources that may have already set prior values
	// that we may not wish to override:
	// * if we called FitPedestal, we have ped_mean and ped_sigma
	// * if we called EstimateMeanNPE, we have mean_npe
	// * if we called EstimateSPE, we have mean_spe_charge and spe_charge_variance,
	//   which can be used to determine spe_firstgamma_mean and spe_firstgamma_shape
	// 
	// we track calls to these functions with booleans that get reset on SetHisto:
	// ran_estimate_mean_npe: determines if we should leave mean_npe as is
	// ran_estimate_spe:      determines if we should use mean_spe_charge and spe_charge_variance
	//                           for calculating spe_firstgamma parameters (shape: TODO),
	// ran_fit_pedestal:      determines if we should leave ped_mean and ped_sigma as is
	
	// ========================================================================
	// STEP 2: ESTIMATE STARTING PARAMETER VALUES
	// ========================================================================
	// general parameters
	// ------------------
	prescaling = 0.1;           // overall scaling of the entire function TODO tune from histogram?
	// TODO this seems redundant given that we also have scaling for the individual components;
	// try removing/fixing this parameter to speed up fitting?
	
	// if we didn't use EstimateSPE, guess some generic low-value mean_npe...
	if(not ran_estimate_mean_npe){
		mean_npe = 0.1;
	}
	
	// highest bin should be pedestal; use it to estimate pedestal gaussian amplitude
	double max_bin_count = thehist->GetBinContent(thehist->GetMaximumBin());
	double pedestal_amp_guess = max_bin_count;
	
	// if PeakScan found an SPE peak use its location as the SPE peak position prior,
	// and to derive estimates of pedestal width and SPE amplitude
	// otherwise take for prior just a bit above mean, and guesses for the others.
	// n.b. we'll choose whether to use these values or others a little later
	double spe_mean_guess;
	double spe_amp_guess;
	double pedestal_sigma_guess;
	if(spe_peak_found){
		spe_amp_guess = max2 / yscaling;
		spe_mean_guess = thehist->GetBinCenter(maxpos2)*1.1;
		pedestal_sigma_guess = thehist->GetBinCenter(interpos)/10.;  // overwrite
	} else {
		// really crude fallback options
		spe_amp_guess = pedestal_amp_guess;
		spe_mean_guess = thehist->GetMean()*1.5;
		pedestal_sigma_guess = 0.1; // based on being in pC ... TODO find a better way
	}
	
	// pedestal parameters
	// -------------------
	ped_scaling = pedestal_amp_guess;                            // pedestal amplitude (gaussian scaling)
	// if we fit a background-only histogram, we already have priors for the mean and sigma,
	// otherwise estimate them now
	if(not ran_fit_pedestal){
		ped_mean = thehist->GetBinCenter(thehist->GetMaximumBin());  // often >0 even after baseline sub
		ped_sigma = pedestal_sigma_guess;                            // pedestal width
	}
	
	// spe parameters
	// --------------
	// the SPE consists of two 'Gaus' functions and an exponential
	
	// first gamma parameters
	if(ran_estimate_spe){
		// if we ran EstimateSPE we have a measurement of the mean spe charge,
		// which we can use to calculate the prior of spe_firstgamma_mean.
		// We need to scale it up a bit as the mean SPE charge isn't the same as the SPE peak charge,
		// but the translation between the two isn't that simple, sooo... fudge it up by 10%?
		spe_firstgamma_scaling = thehist->GetBinContent(thehist->FindBin(mean_spe_charge));
		spe_firstgamma_mean = mean_spe_charge*1.1;
		spe_firstgamma_shape = 10; // TODO derive it from spe_charge_variance
	} else {
		// use our estimate from PeakScan results
		spe_firstgamma_scaling = spe_amp_guess; // controls amplitude of main SPE function bump - but
		// indirectly since we have two contributing gauss functions whose amplitudes are linked to this
		spe_firstgamma_mean = spe_mean_guess;   // controls posn of peak by stretching from the LH edge
		spe_firstgamma_shape = 10;              // an arbitrary shape factor. (1/b in the paper)
	}
	// a little fine tuning
	spe_firstgamma_scaling /= 23.;
	spe_firstgamma_scaling *= (1./spe_firstgamma_mean);
	// the shape parameter controls how skewed left (low values) or symmetric (higher values) the peak is.
	// more symmetric (higher) values also narrows the distribution.
	// minimum of 1, almost a triangle vertical at the left edge, max >100, for nearly gaussian
	
	// secondary gamma parameters
	// second gamma's parameters are scaling factors to those of the first
	// secondary gamma should be at LOWER mean and slightly WIDER shape: both scaling factors should be <1.
	spe_secondgamma_scaling = 0.3;                // amplitude scaling of secondary Gamma
	spe_secondgamma_mean_scaling = 0.3;           // mean scaling. Don't confuse with NPE peak separation!
	spe_secondgamma_shape_scaling = 0.3;          // shape scaling. lower is wider
	
	// exponential parameters
	// SPE expl on a log-plot is a straight line dropping off from pedestal to fill the dip region
	// at most this should not exceed the dip where it fills in
	spe_expl_charge_scaling = 10;                 // defines steepness of expl fill-in
	// expl amplitude scaling defines a vertical shift of the expl fill-in
	// zero would probably be a fair prior... but let's try to get something if we have a dip?
	if(spe_peak_found){
		double dip_charge  = thehist->GetBinCenter(interpos);
		// estimate dip counts
		double exp_scaling = (1./spe_expl_charge_scaling)*exp(-dip_charge*spe_expl_charge_scaling);
		exp_scaling /= 10;  // back it off
		spe_expl_scaling = exp_scaling;
	} else {
		spe_expl_scaling = 0;
	}
	// notes on above calculation:
	// main scaling should be such that it is close to the number of counts at the minium of the dip region
	// i.e. at the q of the dip then (exp_scaling*charge_scaling)exp(-q_dip*charge_scaling) = dip_counts
	// so q_dip = 0.1, dip_counts = 200, with charge_scaling = 10, we have 10*exp_scaling*exp(-1)=200
	// or exp_scaling = 20/0.36 = 55. Back it off a bit to 40....
	// BUT actually, if this contributes too much, then its cutoff at the primary Gamma mean gives rise to a
	// noticable discontinuity, so back it off a lot to 10 or even 5
	// (the rest can be filled in by the secondary Gamma)
	
	// ========================================================================
	// STEP 3: ESTIMATE PARAMETER FITTING RANGES
	// ========================================================================
	
	// build the vector of fit parameter ranges
	// prescaling
	fit_parameter_ranges.resize(14);
	fit_parameter_ranges.at(0) = std::pair<double,double>{0,0};              // no constraints
	// pedestal amplitude
	fit_parameter_ranges.at(1) = std::pair<double,double>{0,max_bin_count*2};
	// pedestal mean
	fit_parameter_ranges.at(2) = std::pair<double,double>{0,spe_mean_guess};
	// pedestal sigma
	fit_parameter_ranges.at(3) = std::pair<double,double>{0.,spe_mean_guess}; // TODO check not too restrictive
	// spe first gamma amplitude
	fit_parameter_ranges.at(4) = std::pair<double,double>{0,max_bin_count};   // TODO check not too restrictive
	// spe first gamma mean
	fit_parameter_ranges.at(5) = std::pair<double,double>{0,spe_range_max};   // TODO check not too restrictive
	// spe first gamma shape
	fit_parameter_ranges.at(6) = std::pair<double,double>{1.,100};
	// spe second gamma amplitude scaling
	fit_parameter_ranges.at(7) = std::pair<double,double>{0,0.75}; // amp↑ as q↓, so even scaling 1 is too high
	// spe second gamma mean scaling
	fit_parameter_ranges.at(8) = std::pair<double,double>{0.1,1};  // secondary gamma should be lower charge
	// spe second gamma shape scaling
	fit_parameter_ranges.at(9) = std::pair<double,double>{0.1,2};  // secondary gammma *should* be wider
	// spe expl amplitude
	fit_parameter_ranges.at(10) = std::pair<double,double>{0,max_bin_count};
	// spe expl charge scaling
	fit_parameter_ranges.at(11) = std::pair<double,double>{3,50};
	// mean npe
	fit_parameter_ranges.at(12) = std::pair<double,double>{0,3.0};  // >1 isn't SPE regime any more!
	// max npes
	fit_parameter_ranges.at(13) = std::pair<double,double>{1,5};
	
	// assume that if we have an internal fit function, we wish to set these as priors
	std::cout<<"Setting parameter priors"<<std::endl;
	RefreshParameters();
	
	// update internal TF1 limits
	SetParameterLimits(fit_parameter_ranges);
	
	// print just to check and for info
	std::cout<<"Piors and their limits are:"<<std::endl;
	for(int pari=0; pari<full_fit_func->GetNpar(); ++pari){
		double lbound, ubound, value;
		full_fit_func->GetParLimits(pari,lbound,ubound);
		value = full_fit_func->GetParameter(pari);
		std::string parname = full_fit_func->GetParName(pari);
		std::cout<<pari<<": ("<<parname<<") "<<lbound<<" < "<<value<<" < "<<ubound<<std::endl;
	}
	std::cout<<std::endl;
	
	// update the user's vector
	if(fit_pars!=nullptr) *fit_pars = GetParameters();
	// update the user's limits
	if(par_limits!=nullptr) *par_limits = fit_parameter_ranges;
	
	return true;
}

TF1* DEAPFitFunction::MakeFullFitFunction(){
	if(full_fit_func!=nullptr){
		std::cerr<<"WARNING: You CANNOT create multiple TF1s "
			 <<"from a single DEAPFitFunction object!"<<std::endl
			 <<"Returning existing fit function."<<std::endl;
	} else {
		ConstructFunctions(); // construct the TF1s that will go into it
		full_fit_func = new TF1("full_fit_func", this, &DEAPFitFunction::FullFitFunction, histogram_minimum, histogram_maximum, 14);
		NameParameters(full_fit_func);
		full_fit_func->FixParameter(13, max_pes);    // fix max_pes
	}
	return full_fit_func;
}

TF1* DEAPFitFunction::NameParameters(TF1* thefunc){
	// if user didn't provide a function, see if we have one internally
	if(thefunc==nullptr){
		// if we have no internal function either, error
		if(full_fit_func==nullptr){
			std::cerr<<"DEAPFitFunction::NameParameters called with nullptr "
				 <<"and we have no internal fitting TF1!"<<std::endl;
			return nullptr;
		} else {
			thefunc=full_fit_func;
		}
	}
	thefunc->SetNpx(2000); // good resolution on our draws
	// it takes a *lot* of parameters...
	thefunc->SetParName(0,"prescaling");
	thefunc->SetParName(1,"ped_scaling");
	thefunc->SetParName(2,"ped_mean");
	thefunc->SetParName(3,"ped_sigma");
	thefunc->SetParName(4,"spe_firstgamma_scaling");
	thefunc->SetParName(5,"spe_firstgamma_mean");
	thefunc->SetParName(6,"spe_firstgamma_shape");
	thefunc->SetParName(7,"spe_secondgamma_scaling");
	thefunc->SetParName(8,"spe_secondgamma_mean_scaling");
	thefunc->SetParName(9,"spe_secondgamma_shape_scaling");
	thefunc->SetParName(10,"spe_expl_scaling");
	thefunc->SetParName(11,"spe_expl_charge_scaling");
	thefunc->SetParName(12,"mean_npe");
	if(thefunc->GetNpar()>13) thefunc->SetParName(13,"max_pes");
	return thefunc;
}

int DEAPFitFunction::ParameterNameToNumber(std::string par_name){
	int par_i=-1;
	if(par_name=="prescaling"){ par_i = 0; }
	if(par_name=="ped_scaling"){ par_i = 1; }
	if(par_name=="ped_mean"){ par_i = 2; }
	if(par_name=="ped_sigma"){ par_i = 3; }
	if(par_name=="spe_firstgamma_scaling"){ par_i = 4; }
	if(par_name=="spe_firstgamma_mean"){ par_i = 5; }
	if(par_name=="spe_firstgamma_shape"){ par_i = 6; }
	if(par_name=="spe_secondgamma_scaling"){ par_i = 7; }
	if(par_name=="spe_secondgamma_mean_scaling"){ par_i = 8; }
	if(par_name=="spe_secondgamma_shape_scaling"){ par_i = 9; }
	if(par_name=="spe_expl_scaling"){ par_i = 10; }
	if(par_name=="spe_expl_charge_scaling"){ par_i = 11; }
	if(par_name=="mean_npe"){ par_i = 12; }
	if(par_name=="max_pes"){ par_i = 13; }
	return par_i;
}

std::string DEAPFitFunction::ParameterNumberToName(int par_i){
	std::string par_name = "unknown";
	if(par_i==0){ par_name = "prescaling"; }
	if(par_i==1){ par_name = "ped_scaling"; }
	if(par_i==2){ par_name = "ped_mean"; }
	if(par_i==3){ par_name = "ped_sigma"; }
	if(par_i==4){ par_name = "spe_firstgamma_scaling"; }
	if(par_i==5){ par_name = "spe_firstgamma_mean"; }
	if(par_i==6){ par_name = "spe_firstgamma_shape"; }
	if(par_i==7){ par_name = "spe_secondgamma_scaling"; }
	if(par_i==8){ par_name = "spe_secondgamma_mean_scaling"; }
	if(par_i==9){ par_name = "spe_secondgamma_shape_scaling"; }
	if(par_i==10){ par_name = "spe_expl_scaling"; }
	if(par_i==11){ par_name = "spe_expl_charge_scaling"; }
	if(par_i==12){ par_name = "mean_npe"; }
	if(par_i==13){ par_name = "max_pes"; }
	return par_name;
}

// Parameter Setters
// TODO these should change the limits and expand them as necessary
void DEAPFitFunction::SetPrescaling(const double &prescaling_in){
	prescaling = prescaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("prescaling", prescaling);
}
void DEAPFitFunction::SetPedScaling(const double &ped_scaling_in){
	ped_scaling = ped_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("ped_scaling", ped_scaling);
}
void DEAPFitFunction::SetPedMean(const double &ped_mean_in){
	ped_mean = ped_mean_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("ped_mean", ped_mean);
}
void DEAPFitFunction::SetPedSigma(const double &ped_sigma_in){
	ped_sigma = ped_sigma_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("ped_sigma", ped_sigma);
}
void DEAPFitFunction::SetSpeFirstGammaScaling(const double &spe_firstgamma_scaling_in){
	spe_firstgamma_scaling = spe_firstgamma_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_firstgamma_scaling", spe_firstgamma_scaling);
}
void DEAPFitFunction::SetSpeFirstGammaMean(const double &spe_firstgamma_mean_in){
	spe_firstgamma_mean = spe_firstgamma_mean_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_firstgamma_mean", spe_firstgamma_mean);
}
void DEAPFitFunction::SetSpeFirstGammaShape(const double &spe_firstgamma_shape_in){
	spe_firstgamma_shape = spe_firstgamma_shape_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_firstgamma_shape", spe_firstgamma_shape);
}
void DEAPFitFunction::SetSpeSecondGammaScaling(const double &spe_secondgamma_scaling_in){
	spe_secondgamma_scaling = spe_secondgamma_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_secondgamma_scaling", spe_secondgamma_scaling);
}
void DEAPFitFunction::SetSpeSecondGammaMeanScaling(const double &spe_secondgamma_mean_scaling_in){
	spe_secondgamma_mean_scaling = spe_secondgamma_mean_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_secondgamma_mean_scaling", spe_secondgamma_mean_scaling);
}
void DEAPFitFunction::SetSpeSecondGammaShapeScaling(const double &spe_secondgamma_shape_scaling_in){
	spe_secondgamma_shape_scaling = spe_secondgamma_shape_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_secondgamma_shape_scaling", spe_secondgamma_shape_scaling);
}
void DEAPFitFunction::SetSpeExplScaling(const double &spe_expl_scaling_in){
	spe_expl_scaling = spe_expl_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_expl_scaling", spe_expl_scaling);
}
void DEAPFitFunction::SetSpeExplChargeScaling(const double &spe_expl_charge_scaling_in){
	spe_expl_charge_scaling = spe_expl_charge_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("spe_expl_charge_scaling", spe_expl_charge_scaling);
}
void DEAPFitFunction::SetMeanNpe(const double &mean_npe_in){
	mean_npe = mean_npe_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("mean_npe", mean_npe);
}
void DEAPFitFunction::SetMaxNpe(const double &max_pes_in){
	max_pes = max_pes_in;
	if(full_fit_func!=nullptr) full_fit_func->SetParameter("max_pes", max_pes);
	ConstructFunctions(); // in case we need to construct more TF1s
}

int DEAPFitFunction::SetParameter(std::string par_name, double par_value){
	// check we know this parameter
	int par_i = ParameterNameToNumber(par_name);
	if(par_i<0){
		std::cerr<<"DEAPFitFunction::SetParameter unknown parameter "<<par_name<<std::endl;
		return -1;
	}
	
	return SetParameter(par_i, par_value);
}

int DEAPFitFunction::SetParameter(int par_i, double par_value){
	// check parameter number in range
	if((par_i<0) || (par_i>fit_parameter_ranges.size())){
		std::cerr<<"DEAPFitFunction::SetParameter parameter number "<<par_i<<" out of range"<<std::endl;
		return -1;
	}
	
	// update the limits first if we need to
	std::pair<double,double> the_limits = fit_parameter_ranges.at(par_i);
	if(par_value<the_limits.first){
		SetParameterLimits(par_i, std::pair<double,double>{par_value,the_limits.second});
	} else if(par_value>the_limits.second){
		SetParameterLimits(par_i, std::pair<double,double>{the_limits.first,par_value});
	}
	
	// convert to name for ease
	std::string par_name = ParameterNumberToName(par_i);
	if(par_name=="unknown"){
		std::cerr<<"DEAPFitFunction::Unknown parameter number "<<par_i<<std::endl;
		return -1;
	}
	
	// set accordingly
	if(par_name=="prescaling"){ prescaling = par_value; }
	if(par_name=="ped_scaling"){ ped_scaling = par_value; }
	if(par_name=="ped_mean"){ ped_mean = par_value; }
	if(par_name=="ped_sigma"){ ped_sigma = par_value; }
	if(par_name=="spe_firstgamma_scaling"){ spe_firstgamma_scaling = par_value; }
	if(par_name=="spe_firstgamma_mean"){ spe_firstgamma_mean = par_value; }
	if(par_name=="spe_firstgamma_shape"){ spe_firstgamma_shape = par_value; }
	if(par_name=="spe_secondgamma_scaling"){ spe_secondgamma_scaling = par_value; }
	if(par_name=="spe_secondgamma_mean_scaling"){ spe_secondgamma_mean_scaling = par_value; }
	if(par_name=="spe_secondgamma_shape_scaling"){ spe_secondgamma_shape_scaling = par_value; }
	if(par_name=="spe_expl_scaling"){ spe_expl_scaling = par_value; }
	if(par_name=="spe_expl_charge_scaling"){ spe_expl_charge_scaling = par_value; }
	if(par_name=="mean_npe"){ mean_npe = par_value; }
	if(par_name=="max_pes"){ max_pes = par_value; }
	RefreshParameters();
	
	return 1;
}

// Parameter fixers: set the parameter and do not vary it during fitting
void DEAPFitFunction::FixPrescaling(const double &prescaling_in){
	prescaling = prescaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(0, prescaling);
	else FixingWarning();
}
void DEAPFitFunction::FixPedScaling(const double &ped_scaling_in){
	ped_scaling = ped_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(1, ped_scaling);
	else FixingWarning();
}
void DEAPFitFunction::FixPedMean(const double &ped_mean_in){
	ped_mean = ped_mean_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(2, ped_mean);
	else FixingWarning();
}
void DEAPFitFunction::FixPedSigma(const double &ped_sigma_in){
	ped_sigma = ped_sigma_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(3, ped_sigma);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeFirstGammaScaling(const double &spe_firstgamma_scaling_in){
	spe_firstgamma_scaling = spe_firstgamma_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(4, spe_firstgamma_scaling);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeFirstGammaMean(const double &spe_firstgamma_mean_in){
	spe_firstgamma_mean = spe_firstgamma_mean_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(5, spe_firstgamma_mean);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeFirstGammaShape(const double &spe_firstgamma_shape_in){
	spe_firstgamma_shape = spe_firstgamma_shape_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(6, spe_firstgamma_shape);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeSecondGammaScaling(const double &spe_secondgamma_scaling_in){
	spe_secondgamma_scaling = spe_secondgamma_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(7, spe_secondgamma_scaling);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeSecondGammaMeanScaling(const double &spe_secondgamma_mean_scaling_in){
	spe_secondgamma_mean_scaling = spe_secondgamma_mean_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(8, spe_secondgamma_mean_scaling);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeSecondGammaShapeScaling(const double &spe_secondgamma_shape_scaling_in){
	spe_secondgamma_shape_scaling = spe_secondgamma_shape_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(9, spe_secondgamma_shape_scaling);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeExplScaling(const double &spe_expl_scaling_in){
	spe_expl_scaling = spe_expl_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(10, spe_expl_scaling);
	else FixingWarning();
}
void DEAPFitFunction::FixSpeExplChargeScaling(const double &spe_expl_charge_scaling_in){
	spe_expl_charge_scaling = spe_expl_charge_scaling_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(11, spe_expl_charge_scaling);
	else FixingWarning();
}
void DEAPFitFunction::FixMeanNpe(const double &mean_npe_in){
	mean_npe = mean_npe_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(12, mean_npe);
	else FixingWarning();
}
void DEAPFitFunction::FixMaxNpe(const double &max_pes_in){
	max_pes = max_pes_in;
	if(full_fit_func!=nullptr) full_fit_func->FixParameter(13, max_pes);
	else FixingWarning();
}

void DEAPFitFunction::FixingWarning(){
	std::cerr<<"WARNING: DEAPFitFunction::FixParameter called, but we have no internal TF1 object"<<std::endl
		 <<"If you constructed your TF1 using this class as a functor, "
		 <<"use TF1::FixParameter instead"<<std::endl;
}

// Group Fit Parameter Setters
int DEAPFitFunction::SetParameters(std::vector<double> fit_parameters){
	if(fit_parameters.size()<12){
		std::cerr<<"WARNING: Too few parameters passed to DEAPFitFunction::SetParameters!"<<std::endl
			 <<"Expected 12 (13 including max_npes), got "<<fit_parameters.size()<<std::endl
			 <<"NO PARAMETERS WILL BE SET"<<std::endl;
		return 0;
	}
	prescaling = fit_parameters.at(0);
	ped_scaling = fit_parameters.at(1);
	ped_mean = fit_parameters.at(2);
	ped_sigma = fit_parameters.at(3);
	spe_firstgamma_scaling = fit_parameters.at(4);
	spe_firstgamma_mean = fit_parameters.at(5);
	spe_firstgamma_shape = fit_parameters.at(6);
	spe_secondgamma_scaling = fit_parameters.at(7);
	spe_secondgamma_mean_scaling = fit_parameters.at(8);
	spe_secondgamma_shape_scaling = fit_parameters.at(9);
	spe_expl_scaling = fit_parameters.at(10);
	spe_expl_charge_scaling = fit_parameters.at(11);
	mean_npe = fit_parameters.at(12);
	if(fit_parameters.size()==14) max_pes = fit_parameters.at(13); // optional fit parameter
	if(fit_parameters.size()>14){
		std::cerr<<"WARNING: Too many parameters passed to DEAPFitFunction::SetParameters!"<<std::endl
			 <<"Expected 13 (14 including max_npes), got "<<fit_parameters.size()<<std::endl
			 <<"Only parameters up to the first 13 (14) will be used"<<std::endl;
	}
	
	// update the internal TF1s
	RefreshParameters();
	return 1;
}
void DEAPFitFunction::SetParameters(double* fit_parameters){
	// sure hope this doesn't segfault
	// set members
	prescaling = fit_parameters[0];
	ped_scaling = fit_parameters[1];
	ped_mean = fit_parameters[2];
	ped_sigma = fit_parameters[3];
	spe_firstgamma_scaling = fit_parameters[4];
	spe_firstgamma_mean = fit_parameters[5];
	spe_firstgamma_shape = fit_parameters[6];
	spe_secondgamma_scaling = fit_parameters[7];
	spe_secondgamma_mean_scaling = fit_parameters[8];
	spe_secondgamma_shape_scaling = fit_parameters[9];
	spe_expl_scaling = fit_parameters[10];
	spe_expl_charge_scaling = fit_parameters[11];
	mean_npe = fit_parameters[12];
	//max_pes = fit_parameters[13];  // better to assume we're not using it, since we can't tell
	
	// update the internal TF1s
	RefreshParameters();
}

int DEAPFitFunction::SetParameterLimits(std::vector<std::pair<double,double>> ranges_in){
	if(ranges_in.size()<13){
		std::cerr<<"WARNING: Too few parameter ranges passed to "
			 <<"DEAPFitFunction::SetParameterLimits!"<<std::endl
			 <<"Expected 13 (14 including max_npes), got "<<ranges_in.size()<<std::endl
			 <<"NO PARAMETERS WILL BE SET"<<std::endl;
		return -1;
	} else if(ranges_in.size()>14){
		std::cerr<<"WARNING: Too many parameters passed to DEAPFitFunction::SetParameterLimits!"<<std::endl
			 <<"Expected 13 (14 including max_npes), got "<<ranges_in.size()<<std::endl
			 <<"Only parameters up to the first 13 (14) will be used"<<std::endl;
	}
	fit_parameter_ranges = ranges_in;
	if(ranges_in.size()==14) fit_parameter_ranges.emplace_back(std::pair<double,double>{1,5});
	
	if(full_fit_func!=nullptr){
		for(int pari=0; pari<fit_parameter_ranges.size(); pari++){
			full_fit_func->SetParLimits(pari,
						    fit_parameter_ranges.at(pari).first,
						    fit_parameter_ranges.at(pari).second);
		}
	}
	return 1;
}

// TODO: check for parameters outside new limits and coerce into range
// TODO: other setters for specific parameter ranges
int DEAPFitFunction::SetParameterLimits(int par_i, std::pair<double,double> range_in){
	if((par_i<0) || (par_i>=fit_parameter_ranges.size())){
		std::cerr<<"DEAPFitFunction::SetParameterLimits parameter number "<<par_i
			 <<"out of range 0-"<<fit_parameter_ranges.size()<<std::endl;
		return -1;
	}
	fit_parameter_ranges.at(par_i) = range_in;
	
	if(full_fit_func!=nullptr){
		for(int pari=0; pari<fit_parameter_ranges.size(); pari++){
			full_fit_func->SetParLimits(pari,
						    fit_parameter_ranges.at(pari).first,
						    fit_parameter_ranges.at(pari).second);
		}
	}
	return 1;
}

std::pair<double,double> DEAPFitFunction::GetParameterLimits(std::string par_name){
	
	// convert name into number, since limits are in a vector not a map
	int par_i = ParameterNameToNumber(par_name);
	if(par_i<0){
		std::cerr<<"DEAPFitFunction::GetParameterLimits for unknown parameter "<<par_name<<std::endl;
		return std::pair<double,double>{0,0};
	}
	
	if(par_i>=fit_parameter_ranges.size()){
		std::cerr<<"DEAPFitFunction::GetParameterLimits internal limits vector too small?"<<std::endl;
		return std::pair<double,double>{0,0};
	}
	
	return fit_parameter_ranges.at(par_i);
}

std::pair<double,double> DEAPFitFunction::GetParameterLimits(int par_i){
	if((par_i<0) || (par_i>=fit_parameter_ranges.size()) ){
		std::cerr<<"DEAPFitFunction::GetParameterLimits; parameter number "<<par_i
			 <<" is out of range 0-"<<(fit_parameter_ranges.size()-1)<<std::endl;
		return std::pair<double,double>{0,0};
	}
	return fit_parameter_ranges.at(par_i);
}

std::vector<double> DEAPFitFunction::GetParameters(){
	std::vector<double> fit_parameters(14);
	fit_parameters.at(0) = prescaling;
	fit_parameters.at(1) = ped_scaling;
	fit_parameters.at(2) = ped_mean;
	fit_parameters.at(3) = ped_sigma;
	fit_parameters.at(4) = spe_firstgamma_scaling;
	fit_parameters.at(5) = spe_firstgamma_mean;
	fit_parameters.at(6) = spe_firstgamma_shape;
	fit_parameters.at(7) = spe_secondgamma_scaling;
	fit_parameters.at(8) = spe_secondgamma_mean_scaling;
	fit_parameters.at(9) = spe_secondgamma_shape_scaling;
	fit_parameters.at(10) = spe_expl_scaling;
	fit_parameters.at(11) = spe_expl_charge_scaling;
	fit_parameters.at(12) = mean_npe;
	fit_parameters.at(13) = max_pes;
	return fit_parameters;
}

double DEAPFitFunction::GetParameter(std::string parname){
	if(parname=="prescaling") return prescaling;
	if(parname=="ped_scaling") return ped_scaling;
	if(parname=="ped_mean") return ped_mean;
	if(parname=="ped_sigma") return ped_sigma;
	if(parname=="spe_firstgamma_scaling") return spe_firstgamma_scaling;
	if(parname=="spe_firstgamma_mean") return spe_firstgamma_mean;
	if(parname=="spe_firstgamma_shape") return spe_firstgamma_shape;
	if(parname=="spe_secondgamma_scaling") return spe_secondgamma_scaling;
	if(parname=="spe_secondgamma_mean_scaling") return spe_secondgamma_mean_scaling;
	if(parname=="spe_secondgamma_shape_scaling") return spe_secondgamma_shape_scaling;
	if(parname=="spe_expl_scaling") return spe_expl_scaling;
	if(parname=="spe_expl_charge_scaling") return spe_expl_charge_scaling;
	if(parname=="mean_npe") return mean_npe;
	if(parname=="max_pes") return max_pes;
	if(parname=="mean_spe_charge") return mean_spe_charge;
	if(parname=="spe_charge_variance") return spe_charge_variance;
	else return std::numeric_limits<double>::max();
}

void DEAPFitFunction::RefreshParameters(){
	// pedestal
	if(pedestal_func){
		pedestal_func->SetParameter("ped_mean",ped_mean);
		pedestal_func->SetParameter("ped_sigma",ped_sigma);
	}
	
	// spe - these parameters aren't actually used; only the npe TF1s are used in the fitting
	// and they make copies of the spe TF1 on construction. Maintain in case the users wants to draw it?
	if(spe_func){
		spe_func->SetParameter("spe_firstgamma_scaling",spe_firstgamma_scaling);
		spe_func->SetParameter("spe_firstgamma_mean",spe_firstgamma_mean);
		spe_func->SetParameter("spe_firstgamma_shape",spe_firstgamma_shape);
		spe_func->SetParameter("spe_secondgamma_scaling",spe_secondgamma_scaling);
		spe_func->SetParameter("spe_secondgamma_mean_scaling",spe_secondgamma_mean_scaling);
		spe_func->SetParameter("spe_secondgamma_shape_scaling",spe_secondgamma_shape_scaling);
		spe_func->SetParameter("spe_expl_scaling",spe_expl_scaling);
		spe_func->SetParameter("spe_expl_charge_scaling",spe_expl_charge_scaling);
	}
	
	//npe
	NPE_pars.resize(npe_funcs.back()->GetNpar());
	// load pars into NPE_pars
	NPE_pars.at(0) = ped_mean;
	NPE_pars.at(1) = ped_sigma;
	for(int j=0; j<max_pes; ++j){
		for(int pari=0; pari<n_spe_pars; ++pari){
			NPE_pars.at(2+(j*n_spe_pars)+pari) = spe_func->GetParameter(pari);
		}
	}
	// pass into functions - these are the parameters of the internal functions to the TF1Convolutions
	for(TF1* apefunc : npe_funcs){
		apefunc->SetParameters(NPE_pars.data());
	}
	
	// pass into the full fit function
	if(full_fit_func!=nullptr){
		// load into a vector for returning or passing to the internal TF1
		// (we need to use a vector because SetParameter(X,value) only allows
		// setting up to 12 parameters at most
		std::vector<double> fit_parameters(14);
		fit_parameters.at(0) = prescaling;
		fit_parameters.at(1) = ped_scaling;
		fit_parameters.at(2) = ped_mean;
		fit_parameters.at(3) = ped_sigma;
		fit_parameters.at(4) = spe_firstgamma_scaling;
		fit_parameters.at(5) = spe_firstgamma_mean;
		fit_parameters.at(6) = spe_firstgamma_shape;
		fit_parameters.at(7) = spe_secondgamma_scaling;
		fit_parameters.at(8) = spe_secondgamma_mean_scaling;
		fit_parameters.at(9) = spe_secondgamma_shape_scaling;
		fit_parameters.at(10) = spe_expl_scaling;
		fit_parameters.at(11) = spe_expl_charge_scaling;
		fit_parameters.at(12) = mean_npe;
		fit_parameters.at(13) = max_pes;
		full_fit_func->SetParameters(fit_parameters.data());
	}
}

std::vector<std::pair<double,double>> DEAPFitFunction::GetParameterLimits(){
	return fit_parameter_ranges;
}

TF1* DEAPFitFunction::GetFullFitFunction(){
	return full_fit_func;
}

TF1* DEAPFitFunction::GetPedFunc(){
	return pedestal_func;
}

TF1* DEAPFitFunction::GetSPEFunc(){
	return spe_func;
}

std::vector<TF1*> DEAPFitFunction::GetNPEFuncs(){
	return npe_funcs;
}

std::vector<TF1Convolution*> DEAPFitFunction::GetNPEConvs(){
	return npe_convolns;
}

std::vector<double> DEAPFitFunction::GetNPEPars(){
	return NPE_pars;
}

TF1* DEAPFitFunction::GetBackgroundFunc(){
	return pedestal_standalone;
}

double DEAPFitFunction::GetXscaling(){
	return xscaling;
}

double DEAPFitFunction::GetYscaling(){
	return yscaling;
}

TFitResultPtr DEAPFitFunction::FitTheHisto(std::string opts){
	if(opts.find("S") == std::string::npos){ opts.push_back('S'); }
	TFitResultPtr fitresult = thehist->Fit(full_fit_func, opts.c_str());  // 0 means sucess.
	SetParameters(full_fit_func->GetParameters());     // update the internal members with the fit results
	return fitresult;
}

double DEAPFitFunction::GetMeanSPECharge(){
	// the mean SPE charge is not directly the spe_firstgamma_mean, due since the SPE peak
	// contains other components. We therefore need to extract the mean from the SPE function
	// as a whole
	// first make sure our SPE TF1 is up-to-date with the latest parameters
	RefreshParameters();
	// Now measure it's mean
	mean_spe_charge = spe_func->Mean(histogram_minimum,histogram_maximum);
	return mean_spe_charge;
}

// ========================================================================
// FIT FUNCTION FROM DEAP-3600
// arXiv:1705.10183v3
// ========================================================================

double DEAPFitFunction::Gamma(double* x, double* gamma_pars){
	if((*x)<0) return 0; // goes nuts below 0... :/
	double gamma_mean = gamma_pars[0];
	double gamma_shape = gamma_pars[1]; // this is 1/b in the paper; only it's inverse is ever used
	double gamma_product = gamma_mean/gamma_shape;
	double x_over_gamma_product = (*x)/gamma_product;
	double returnval = 1./(gamma_product*TMath::Gamma(gamma_shape));
	returnval *= pow(x_over_gamma_product,gamma_shape-1.);
	// maybe protecting exp against silly values will prevent it from hanging...
	// exp(30) ~ 1E13, probably don't want to go any further than that.
	if(x_over_gamma_product>30)  returnval=0;
	if(x_over_gamma_product<-30) returnval *= 1E13;
	else                         returnval *= exp(-x_over_gamma_product);
	if(returnval<0) returnval=0;
	return returnval;
}

double DEAPFitFunction::SPE_Func(double* x, double* SPE_pars){
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
	
	// calculate exponential fill-in contribution
	double expcontrib=0;
	double exparg = (*x)*expl_charge_scaling;  // note we do exp(-exparg)
	if( ((*x)>firstgamma_mean) || ((*x)<0) ) expcontrib=0;
	// ^ this 'cutoff at the mean value of the primary SPE gamma' gives a discontinuity which is
	// noticable if the expl function is too large of a component... use small parameters
	// added a cutoff of the expl at 0; we need to be able to extend the SPE function out to
	// somewhere that it should tail out to 0 for the convolution. Seems to tie up with paper.
	// Note this still gives a discontinuity at the end of the function, which TF1Convolution
	// says could be a problem, but it seems to be okay.
	// 
	// maybe protecting exp against silly values will prevent it from hanging...
	// exp(30) = 1E13, probably don't want to go any further than that.
	else if(exparg>30)  expcontrib=0;
	else if(exparg<-30) expcontrib = 1E13;
	else                expcontrib = (expl_scaling*expl_charge_scaling)*exp(-exparg);
	
	returnval += expcontrib;
	
	if(returnval<0) returnval=0;
	return returnval;
}

void DEAPFitFunction::ConstructFunctions(){
	gROOT->cd();
	
	// firstly a standalone gaussian for fitting background-only histograms
	// the reason we don't use ped_func is that ped_func doesn't include an amplitude parameter,
	// and since we convolute with it a bunch, it's better not to add unnecessary parameters
	if(not pedestal_standalone){
		pedestal_standalone = new TF1("pedestal_standalone", "gaus");
		// "gaus" is [0]*exp(-0.5*((x-[1])/[2])**2)
		pedestal_standalone->SetParName(0,"ped_amp");
		pedestal_standalone->SetParName(1,"ped_mean");
		pedestal_standalone->SetParName(2,"ped_sigma");
		pedestal_standalone->SetNpx(1000);
	}
	
	// build the TF1 representing the pedestal component, if we have not yet done so
	if(not pedestal_func){
		std::cout<<"Creating Pedestal TF1"<<std::endl;
		// TODO 'true' arg to TMath::Gaus adds 1/2*pi*sigma scaling to gaus, which matches paper function.
		// however for small widths, this results in a VERY large pedestal amplitude....
		// maybe we should set this to false (maximum=1) and let the firstgamma_scaling do the work...
		pedestal_func = new TF1("pedestal_func","TMath::Gaus(x, [0], [1], true)", ped_range_min, ped_range_max);
		pedestal_func->SetParName(0,"ped_mean");  // true adds normalization by 1/(sqrt(2pi)*sigma)
		pedestal_func->SetParName(1,"ped_sigma"); // matches DEAP paper
	}
	
	// build the TF1 representing a SPE peak, if we have not yet done so
	if(not spe_func){
		std::cout<<"Creating SPE TF1"<<std::endl;
		spe_func = new TF1("spe_func",this,&DEAPFitFunction::SPE_Func,spe_range_min,spe_range_max,8); // 8 is num parameters it takes...!
		spe_func->SetParName(0,"spe_firstgamma_scaling");
		spe_func->SetParName(1,"spe_firstgamma_mean");
		spe_func->SetParName(2,"spe_firstgamma_shape");
		spe_func->SetParName(3,"spe_secondgamma_scaling");
		spe_func->SetParName(4,"spe_secondgamma_mean_scaling");
		spe_func->SetParName(5,"spe_secondgamma_shape_scaling");
		spe_func->SetParName(6,"spe_expl_scaling");
		spe_func->SetParName(7,"spe_expl_charge_scaling");
		n_spe_pars = spe_func->GetNpar();
	}
	
	// loop over the NPE peaks
	for(int i=0; i<max_pes; ++i){
		// if we do not have a TF1 corresponding to this NPE peak - build it
		if(npe_funcs.size()<(i+1)){
			// each Npe peak is a convolution of the pedestal,
			// together with the SPE peak convolved with itself N times
			// we re-use the results of previous convolutions to build on for the next
			TF1* n_minus_one_func = (i==0) ? pedestal_func : npe_funcs.at(i-1);
			// convolve it with the SPE function again
			std::cout<<"Creating NPE convolution "<<i<<std::endl;
			TF1Convolution* new_conv = new TF1Convolution(n_minus_one_func, spe_func, npe_func_min, npe_func_max, true);
			// internally this copy-constructs two NEW TF1s based on it's arguments
			// so the passed-in functions have no further influence over it
			// (i.e. we do not need to maintain spe_func parameters beyond this method)
			//new_conv->SetNofPointsFFT(1000);  // default, should be at least this
			npe_convolns.push_back(new_conv);
			// make the function from the convolution
			std::cout<<"Creating NPE TF1 "<<i<<std::endl;
			TString funcname = TString::Format("npe_func_%d",i);
			TF1* npe_func = new TF1(funcname,*new_conv, npe_func_min, npe_func_max, new_conv->GetNpar());
			// this has no special overload for TF1Convolutions, it's just that a TF1Convolution
			// implements the operator()(double* x, double* p) method.
			// EvalPar calls will update the internal functions, and if necessary, re-calculate
			// the convolution.
			//
			// The hitch here? is that in the paper, the SPE function is "normalized to unity".
			// There is TF1::SetNormalized which will calculate the integral
			// and normalize the return values... but it doesn't account for changes in parameters!!
			// (the code even implies it *should* call TF1::Update, but it doesn't... 
			// Without hacking ROOT, i don't see how this can be achieved...
			// As a result, the SPE functions do not just vary in shape, they vary in integral,
			// such that e.g. increasing the SPE peak location decreases its amplitude
			// This is all the kind of mess that makes calculating priors a nightmare.
			npe_funcs.push_back(npe_func);
		}
	}
	
	// increase parameter buffer size if we don't have enough space
	NPE_pars.resize(npe_funcs.back()->GetNpar());
}

double DEAPFitFunction::FullFitFunction(double* x, double* all_pars){
	double prescaling_B = all_pars[0];
	double ped_scaling = all_pars[1];
	double* ped_pars = &all_pars[2]; // parameters 2 and 3 are pedestal gaussian mean and sigma
	double* spe_pars = &all_pars[4]; // parameters 4 to 11 are parameters that describe the SPE shape
	double mean_npe = all_pars[12];  // can we pull initial fit value from occupancy?
	
	// For the NPE peak TF1's, many of the parameters are redundant due to the way they're built
	// That means we need to pass it the same parameter value multiple times.
	// TODO can we remove this redundancy somehow...? Fixing parameter e.g. 3 to be the same as parameter 1?
	// build the necessary array of parameters
	// the first couple are always the pedestal
	NPE_pars.at(0) = ped_pars[0];
	NPE_pars.at(1) = ped_pars[1];
	// after that we keep re-convolving with the same SPE function, with the same parameters
	for(int j=0; j<max_pes; ++j){
		for(int pari=0; pari<n_spe_pars; ++pari){
			NPE_pars.at(2+(j*n_spe_pars)+pari) = spe_pars[pari];
		}
	}
	
	// the contribution from pedestal alone has its own scaling
	double pedestal_contribution = ped_scaling*pedestal_func->EvalPar(x, ped_pars);
	
	// because we can't properly normalize the internals of TF1Convolution,
	// compensate for it by calculating the internal scaling owing from the convolutions
	spe_func->SetParameters(spe_pars);
	double rel_tolerance = 0.0001; // default 1.e-12.
	double spe_integral = spe_func->Integral(histogram_minimum,histogram_maximum,rel_tolerance);
	/* Lowering the relative tolerance probably speeds it up, and suppresses GSL warnings:
	Error in <GSLError>: Error 18 in qags.c at 548 : cannot reach tolerance because of roundoff error
	Warning in <TF1::IntegralOneDim>: Error found in integrating function spe_func ...
	*/
	
//	// it may be quicker to use IntegralFast?
//	// not sure under what circumstances the GaussLegendreSamplingPoints need to be recalculated though...
//	Int_t np = 1000;
//	double *x=new double[np]; these would need to be moved outside so we don't leak memory
//	double *w=new double[np];
//	spe_func->CalcGaussLegendreSamplingPoints(np,x,w,rel_tolerance); // 0 is params pointer (default 0)
//	double spe_integral = spe_func->IntegralFast(np,x,w,histogram_minimum,histogram_maximum);
	
	double running_npe_contribution=0;
	// loop over the NPE peaks
	for(int i=0; i<max_pes; ++i){
		// retrieve the TF1 descrbing this NPE peak
		TF1* this_npe_func = npe_funcs.at(i);
		
		// evaluate the contribution
		double this_npe_peak_contribution = this_npe_func->EvalPar(x, NPE_pars.data());
		
		// compensate for failure to internally normalize
		double normalization_scaling = pow(spe_integral,-(i+1));
		
		// each NPE peak will have an occupancy scaling
		// based on a Poisson distribution parameterised by the mean number of PE on the tube
		double poisson_scaling = TMath::Poisson(i+1,mean_npe);  // +1 because i=0 is 1PE, not 0PE
		
		// evaluate this npe's contribution at the specified charge
		running_npe_contribution += (normalization_scaling*poisson_scaling*this_npe_peak_contribution);
	}
	
	// the total number of events with a given charge is then the sum of contributions from pedestal
	// plus SPE plus 2PE plus 3PE ...
	double total_contribution = pedestal_contribution + running_npe_contribution;
	// and there's an overall scaling factor
	return (prescaling_B*total_contribution);
}

// for using this class to create a TF1 externally
double DEAPFitFunction::operator()(double *x, double *p){
	return FullFitFunction(x,p);
}

// We may be able to speed up fitting by reducing tolerances or changing out integration strategy
// from https://root-forum.cern.ch/t/speeding-up-fitting-to-a-landau-distribution/25140/2
//ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
//ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);


//	spe_func->SetParameters(spe_pars);
//	
//	double running_npe_contribution=0;
//	// loop over the NPE peaks
//	for(int i=0; i<max_pes; ++i){
//		// retrieve the TF1 descrbing this NPE peak
//		TF1* this_npe_func = npe_funcs.at(i);
//		
//		// set it's parameters - the TF1Convolution of 2 functions (or the TF1 derived from it)
//		// combines the parameters of both TF1s that went into it
//		// the first couple of parameters are for the pedestal function
//		this_npe_func->SetParameters(0,ped_pars[0]);
//		this_npe_func->SetParameter(1,ped_pars[1]);
//		
//		// after that we keep re-convolving with the same SPE function, with the same parameters
//		int n_spe_pars = spe_func->GetNpar();
//		for(int j=0; j<i; ++j){
//			for(int pari=0; pari<n_spe_pars; ++pari){
//				this_npe_func->SetParameter(2+(j*n_spe_pars)+pari,spe_pars[pari]);
//				this_npe_func->SetParameter(2+(j*n_spe_pars)+pari,spe_pars[pari]);
//			}
//		}
//		// this should now describe the shape of the distribution of the NPE peak
//		double this_npe_peak_contribution = this_npe_func->Eval(*x);
//		
//		// each NPE peak will have an occupancy scaling
//		// based on a Poisson distribution parameterised by the mean number of PE on the tube
//		double poisson_scaling = TMath::Poisson(i+1,mean_npe);  // +1 because i=0 is 1PE, not 0PE
//		// evaluate this npe's contribution at the specified charge
//		running_npe_contribution += (poisson_scaling*this_npe_peak_contribution);
//	}

