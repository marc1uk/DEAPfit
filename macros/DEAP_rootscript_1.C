{
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

// We must set the ranges of all functions, and for TF1Convolution, they MUST tend to 0 at their ends!
// suitable for this particular histogram we're looking at
ped_range_min=-0.2;    // lower range of the pedestal gaussian function
ped_range_max=0.5;     // upper range of the pedestal gaussian function
spe_range_min=-0.2;    // lower range of the SPE function
spe_range_max=1.0;    // upper range of the SPE function
npe_func_min=-0.2;     // each npe function requires a range, but specifying many ranges seems overkill
npe_func_max=3.0;     // they'll probably all need to have a range covering the full span of charges

//std::string filename = "/home/marc/Documents/PhD/PMTs/PMT_Test_Station/results/DEAP_refit/gain1300V_WB012.C";
//std::string filename = "/annie/app/users/moflaher/gain1300V_WB012.C";
std::string filename = "/home/marc/Documents/PhD/PMTs/PMT_Test_Station/results/DEAP_refit/detkey_334.C";
std::string rootcommand=".x " + filename;
std::cout<<"Executing command "<<rootcommand<<std::endl;
// why does this take so long to open??
gROOT->ProcessLine(rootcommand.c_str());
//for(int i=0; i<5; ++i){
//  gSystem->ProcessEvents();
//  std::this_thread::sleep_for(std::chrono::seconds(1));
//}
gSystem->ProcessEvents();

// now we need to pull the histogram...
std::cout<<"searching for histogram"<<std::endl;
TCanvas* c1 = (TCanvas*)gROOT->GetListOfCanvases()->At(0);
TH1* thehist=nullptr;
TIter nextobject(c1->GetListOfPrimitives());
TObject* obj=nullptr;
for(int i=0; i<50; ++i){
  auto&& aret = (obj=nextobject());
  std::cout<<"ret="<<aret;
  if(not aret) break;
  std::cout<<"primitive: "<<obj->GetName()<<" at "<<obj
           <<", inherits from TH1="<<obj->InheritsFrom("TH1")<<std::endl;
  if(obj->InheritsFrom("TH1")){
    std::cout<<"it's the histogram!"<<obj->GetName()<<std::endl;
    thehist=(TH1*)obj;
    break;
  }
}

// check the hist and get the canvas
if(thehist==nullptr){ std::cerr<<"No hist"<<std::endl; return 0; }
gPad->SetLogy();

//// scale the hist X axis to be in picocoulombs - Test Stand ver
//double ELECTRON_CHARGE = 1.6E-19;
//double ELECTRON_CHARGE_IN_PICOCOULOMBS = ELECTRON_CHARGE*1E12;
//thehist->GetXaxis()->Set(thehist->GetNbinsX(),
//  thehist->GetXaxis()->GetXmin()*ELECTRON_CHARGE_IN_PICOCOULOMBS,
//  thehist->GetXaxis()->GetXmax()*ELECTRON_CHARGE_IN_PICOCOULOMBS);

// scale the x axis to be in picocoulombs - ToolAnalysis ver
double nano_to_pico = 1000;
thehist->GetXaxis()->Set(thehist->GetNbinsX(),
                     thehist->GetXaxis()->GetXmin()*nano_to_pico,
                     thehist->GetXaxis()->GetXmax()*nano_to_pico);

// We must set the ranges of all functions, and for TF1Convolution, they MUST tend to 0 at their ends!
double histogram_minimum = thehist->GetXaxis()->GetXmin();
double histogram_maximum = thehist->GetXaxis()->GetXmax();
histogram_minimum = -histogram_maximum;
histogram_maximum *=2.;
ped_range_min=histogram_minimum;  // lower range of the pedestal gaussian function
ped_range_max=histogram_maximum;  // upper range of the pedestal gaussian function
spe_range_min=histogram_minimum;  // lower range of the SPE function
spe_range_max=histogram_maximum;  // upper range of the SPE function
npe_func_min=histogram_minimum;   // each npe function could potentially have a different range
npe_func_max=histogram_maximum;   // but i don't think reducing these has any benefit?

// TODO derive this from occupancy
double mean_pe_guess = 0.1;

// maybe we could try smoothing the histogram to make it easier to find the SPE peak?
//TH1F acopyhist(*thehist);
//acopyhist.Smooth(1);

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

// normalise the histogram to unit area? or amplitude?
double scaling = thehist->Integral();
//double scaling = thehist->GetBinContent(thehist->GetMaximumBin());
thehist->Scale(1./scaling);
c1->Modified();
double max2scaled = max2 / scaling;               // used for estimating spe amplitude

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

double spe_exp_charge_scaling = 10;   // expl charge scaling: defines steepness of the expl fill-in
double dip_charge  = thehist->GetBinCenter(interpos);
double exp_scaling = (1./spe_exp_charge_scaling)*exp(-dip_charge*spe_exp_charge_scaling); // ~ dip counts
exp_scaling /= 10;  // back it off

}
