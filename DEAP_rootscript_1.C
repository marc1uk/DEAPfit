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

std::string filename = "/home/marc/Documents/PhD/PMTs/PMT_Test_Station/results/DEAP_refit/gain1300V_WB012.C";
//std::string filename = "/annie/app/users/moflaher/gain1300V_WB012.C";
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

// scale the hist X axis to be in picocoulombs
double ELECTRON_CHARGE = 1.6E-19;
double ELECTRON_CHARGE_IN_PICOCOULOMBS = ELECTRON_CHARGE*1E12;
thehist->GetXaxis()->Set(thehist->GetNbinsX(),
  thehist->GetXaxis()->GetXmin()*ELECTRON_CHARGE_IN_PICOCOULOMBS,
  thehist->GetXaxis()->GetXmax()*ELECTRON_CHARGE_IN_PICOCOULOMBS);

}
