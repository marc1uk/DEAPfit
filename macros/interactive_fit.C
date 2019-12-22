TFile* _file0 = TFile::Open("LEDRun1178S0LEDs5And6And10And28_PulseWindowOnly_PMTStability_Run0.root");
TH1F* thehist = (TH1F*)_file0->Get("hist_charge_332");
gSystem->Load("./libdeapfit.so");
DEAPFitFunction deapfitter(3);
deapfitter.SetHisto(thehist);
deapfitter.ScaleHistoXrange(1000);
double histmax = thehist->GetXaxis()->GetXmax();
double histmin = thehist->GetXaxis()->GetXmin();
deapfitter.SmoothHisto();
deapfitter.PeakScan();
deapfitter.ScaleHistoYrange();
std::vector<double> fitpars;
std::vector<std::pair<double,double>> parlimits;
deapfitter.GeneratePriors(&fitpars, &parlimits);
//TF1* full_fit_func = new TF1("full_fit_func",&deapfitter,&DEAPFitFunction::FullFitFunction,histmin,histmax,14);
deapfitter.MakeFullFitFunction();
TF1* full_fit_func = deapfitter.GetFullFitFunction();
deapfitter.NameParameters(full_fit_func);
//full_fit_func->SetParameters(fitpars.data());
deapfitter.SetParameters(fitpars);
//for(int pari=0; pari<full_fit_func->GetNpar(); pari++) std::cout<<pari<<": "<<full_fit_func->GetParName(pari)<<" = "<<fitpars.at(pari)<<std::endl;
std::cout<<"{"; for(auto&& apar : fitpars) std::cout<<apar<<", "; std::cout<<"\b\b}"<<std::endl;
for(int pari=0; pari<parlimits.size(); pari++){
  full_fit_func->SetParLimits(pari,parlimits.at(pari).first,parlimits.at(pari).second);
}
thehist->Draw();
thehist->ResetBit(kCanDelete);
thehist->SetDirectory(0);
gPad->SetLogy();
full_fit_func->Draw("same");


TF1* spefunc = deapfitter.GetSPEFunc();
spefunc->SetLineColor(kBlue);
spefunc->Draw("same")
TF1* pedfunc = deapfitter.GetPedFunc();
pedfunc->SetLineColor(kViolet);
pedfunc->Draw("same");
std::vector<TF1*> npefuncs = deapfitter.GetNPEFuncs();
npefuncs.at(0)->SetLineColor(kGreen);
npefuncs.at(0)->Draw("same");
npefuncs.at(1)->SetLineColor(kBlue);
npefuncs.at(1)->Draw("same");
