{

// IMPORTANT: BEFORE WE CAN DEFINE FUNCTIONS IN THE ROOT INTERPRETER
// WE NEED TO ENTER 'RAW INPUT' MODE by typing '.rawInput' on the prompt
// then doing the same to exit once we're done. Note: we cannot use
// gROOT->ProcessLine to execute this...
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
  return returnval;
}

double SPE_Func(double* x, double* SPE_pars){
    //if((*x)<0) return 0; // safer, don't do anything weird...
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
    pedestal_func = new TF1("pedestal_func","TMath::Gaus(x, [0], [1], true)", ped_range_min, ped_range_max);
    pedestal_func->SetParName(0,"pedestal_mean");
    pedestal_func->SetParName(1,"pedestal_width");
  }
  // update it's parameters
  pedestal_func->SetParameters(ped_pars);
  
  // the contribution from pedestal alone has its own scaling
  double pedestal_contribution = ped_scaling*pedestal_func->Eval(*x);
  
  // build the TF1 representing a SPE peak, if we have not yet done so
  if(not spe_func){
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
      TF1Convolution* new_conv = new TF1Convolution(n_minus_one_func, spe_func, npe_func_min, npe_func_max, true);
      new_conv->SetRange(npe_func_min, npe_func_max);
      //new_conv->SetNofPointsFFT(1000);  // default, should be at least this
      npe_convolns.push_back(new_conv);
      // make the function from the convolution
      TString funcname = TString::Format("npe_func_%d",i);
      TF1* npe_func = new TF1(funcname,*new_conv, npe_func_min, npe_func_max, new_conv->GetNpar());
      npe_funcs.push_back(npe_func);
    }
    
    // retrieve the TF1 descrbing this NPE peak
    TF1* this_npe_func = npe_funcs.at(i);
    
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
    double poisson_scaling = TMath::Poisson(i,mean_npe);
    // evaluate this npe's contribution at the specified charge
    running_npe_contribution += (poisson_scaling*this_npe_peak_contribution);
    
  }
  
  // the total number of events with a given charge is then the sum of contributions from pedestal
  // plus SPE plus 2PE plus 3PE ... 
  double total_contribution = pedestal_contribution + running_npe_contribution;
  // and there's an overall scaling factor
  return (prescaling_B*total_contribution);
}

}
