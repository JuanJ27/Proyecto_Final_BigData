using namespace RooFit; 
void CodigoGor(){
  // I m p o r t   T T r e e   i n t o   a   R o o D a t a S e t
  // -----------------------------------------------------------
  TFile f("reducetree_Bujk_AOD_3_best1.root");
  TTree* t = (TTree*) f.Get ("butree");
  Double_t B_mass, J_mass, bcpt, pion1pt, etatkl, rapidityB, pdl, pdle, p1charg, dxysig1, etamu1, etamu2, mu1pt, mu2pt, muptmin, muptmax;
  UInt_t trijk;

  t->SetBranchAddress("B_mass",&B_mass);
  t->SetBranchAddress("J_mass",&J_mass);
  t->SetBranchAddress("bcpt",&bcpt);
  t->SetBranchAddress("pion1pt",&pion1pt);
  t->SetBranchAddress("etatkl",&etatkl);
  t->SetBranchAddress("rapidityB",&rapidityB);
  t->SetBranchAddress("pdl",&pdl);
  t->SetBranchAddress("p1charg",&p1charg);
  t->SetBranchAddress("dxysig1",&dxysig1);
  t->SetBranchAddress("etamu1",  &etamu1);
  t->SetBranchAddress("etamu2",  &etamu2);
  t->SetBranchAddress("mu1pt",  &mu1pt);
  t->SetBranchAddress("mu2pt",  &mu2pt); 
  t->SetBranchAddress("muptmin", &muptmin);
  t->SetBranchAddress("muptmax", &muptmax); 
  t->SetBranchAddress("trijk",  &trijk);

  std::cout << "These are the columns B_mass, J_mass, and bcpt:" << std::endl;
  //for (auto branch : *t->GetListOfBranches()) {
    //std::cout << "Branchs: " << branch->GetName() << std::endl;
  //}

  RooRealVar Bm("B_mass","B_mass",5.0,5.6);
  
  RooDataSet *datatre = new RooDataSet("datatre", "datatre", RooArgSet(Bm));

    for (int i = 0; i < t->GetEntries(); i++) {
    t->GetEvent(i);
    //Mass windows cuts
    if(B_mass<=5.0 || B_mass>=5.6) continue;
    if(J_mass<=2.9 || J_mass>=3.3) continue;
    if(pion1pt<1.2) continue;
    if(bcpt<20.0 || bcpt>23.0) continue;
    //if(mu1pt>4.0 && mu2pt>4.0) continue;

    Bm = B_mass;
    datatre->add(RooArgSet(Bm));
    }

  datatre->Print("***\n V \n ***");
  ///
  RooRealVar mean("mean","mean of gaussians",5.2,5.35) ;
  RooRealVar sigma1("sigma1","sigma1 of gaussians",0.01,0.15) ;
  RooGaussian sig1("sig1","Signal component 1",Bm,mean,sigma1) ;
  
  // Build Chebychev 1st polynomial pdf
  RooRealVar a0("a0","a0",-10,10) ;
  RooChebychev bkg("bkg","Background",Bm,a0) ;
    
  // Associated nsig/nbkg as expected number of events with sig/bkg _in_the_range_ "signalRange"
  RooRealVar nsig("nsig","number of signal events",0.,15000) ;
  RooRealVar nbkg("nbkg","number of background events",0,15000) ;

  //model
  RooAddPdf  model("model","sig+bkg+bkg2", RooArgList(sig1,bkg), RooArgList(nsig,nbkg)) ;
  
  
  // S a m p l e   d a t a ,   f i t   m o d e l
  // -------------------------------------------

  RooFitResult* result = model.fitTo(*datatre,Extended(kTRUE),Save(kTRUE)); 

  // Plot datatre and PDF overlaid
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.42, 1, 1);
  //pad1->SetBottomMargin(0); // Sin margen inferior en el primer pad
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.43);
  pad2->SetTopMargin(0); // Sin margen superior en el segundo pad
  pad2->SetBottomMargin(0.25); // Margen inferior más grande en el segundo pad
  pad2->Draw();

  pad1->cd();

  RooPlot* xframe = Bm.frame();

  datatre->plotOn(xframe,Name("Datos"));
  model.plotOn(xframe,Name("Fit"));
  RooHist* hpullm2 = xframe->pullHist();
  model.plotOn(xframe,Components(sig1),LineColor(kRed),LineWidth(2),Name("Signal")); 
  model.plotOn(xframe, Components(bkg), LineColor(kGreen), LineStyle(kDashed),Name("Background")); 
  xframe->GetXaxis()->CenterTitle();
  xframe->GetYaxis()->SetNdivisions(505,1);
  xframe->GetXaxis()->SetNdivisions(505,1);   
  xframe->Draw();

  TLegend *legend = new TLegend(0.17, 0.67, 0.37, 0.87);
  legend->AddEntry(xframe->findObject("Signal"), "Signal");
  legend->AddEntry(xframe->findObject("Datos"), "Data");
  legend->AddEntry(xframe->findObject("Fit"), "Fit");
  legend->AddEntry(xframe->findObject("Background"), "Background");
  legend->AddEntry(xframe->findObject("Background2"), "Background2");
  legend->SetBorderSize(0);
  legend->Draw();

  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;

  Double_t G = sigma1.getVal()*1000.0;
  Double_t GE = sigma1.getError()*1000.0;

  TLegend *legpar = new TLegend(0.55,0.68,0.75,0.88);
  legpar->SetTextSize(0.035);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("Mass(B_{c}^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B_{c}^{+}} = %1.2f #pm %1.2f",nsig.getVal(),nsig.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.2f #pm %1.2f",nbkg.getVal(),nbkg.getError()),"");
  legpar->Draw();
  xframe->SetXTitle("M(J/\\psi\\pi^{+}) [GeV]");
  xframe->SetYTitle("Events/10 MeV");

  datatre->Print("v");
  result->Print("v");

  pad2->cd();
  RooPlot* framem2 = Bm.frame(Title("Pull Distribution"));
  framem2->addPlotable(hpullm2, "P");
  framem2->SetYTitle("(Data-Fit)/#sigma");
  framem2->SetXTitle("M(J/\\psi\\pi^{+}) [GeV]");
  framem2->GetXaxis()->CenterTitle();
  framem2->GetXaxis()->SetNdivisions(505,1);   
  framem2->Draw();
  TLine *line = new TLine(5.0, 0, 5.6, 0);
  line->SetLineColor(kRed);
  line->SetLineStyle(kDashed);
  line->Draw();


  canvas->SaveAs("mass&pull.pdf");

  //parameters output file
  //ofstream salida_TotalFit("output_resulTFit_1_2.txt");
  //salida_TotalFit.is_open();
  //save_resultTfit(salida_TotalFit, nsig, nbkg, a0, mean, sigma1);

}