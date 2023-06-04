using namespace RooFit; 

void analisis(){
  TFile *f = new TFile("./../../../ROOTS_Bdkstar_AOD_2018_3_best1.root");
  TTree* t = (TTree*) f->Get ("treeS");

  Double_t        massB;
  Double_t        P2;
  Double_t        Bpt;
  Double_t        Brapidity;
  Double_t        Bphi;
  Double_t        Bdl;
  Double_t        BdlE;
  Double_t        BdlIP;
  Double_t        BdlIPE;
  Double_t        masskstar;
  Double_t        massJ;
  Double_t        Jpsipt;
  Double_t        Jrapidity;
  Double_t        Jphi;
  Double_t        sigLxyJ;
  Double_t        cosalfaJ;
  Double_t        sigLxyB;
  Double_t        cosalfaB;
  Int_t           pi1charge;
  Int_t           pi2charge;
  Double_t        mu1pt;
  Double_t        mu2pt;
  Double_t        mu1phi;
  Double_t        mu2phi;
  Double_t        mu1eta;
  Double_t        mu2eta;
  Double_t        pi1pt;
  Double_t        pi1eta;
  Double_t        pi1phi;
  Double_t        pi2pt;
  Double_t        pi2eta;
  Double_t        pi2phi;
  Double_t        Bpro;
  Double_t        Jpro;
  Double_t        dxypi1;
  Double_t        dxyEpi1;
  Double_t        dxysigpi1;
  Double_t        DRmumu;
  UInt_t          nPV;
  Int_t           runn;
  Int_t           evtn;
  Int_t           Lumiblock;
  UInt_t          Triggerjk;
  UInt_t          Triggerjkjk;
  Double_t        maxPtpi;
  Double_t        minPtpi;
  Double_t        ptpipi;
  Double_t        masspipi;
  Double_t        masskk;

  t->SetBranchAddress("massB", &massB);
  t->SetBranchAddress("P2", &P2);
  t->SetBranchAddress("Bpt", &Bpt);
  t->SetBranchAddress("Brapidity", &Brapidity);
  t->SetBranchAddress("Bphi", &Bphi);
  t->SetBranchAddress("Bdl", &Bdl);
  t->SetBranchAddress("BdlE", &BdlE);
  t->SetBranchAddress("BdlIP", &BdlIP);
  t->SetBranchAddress("BdlIPE", &BdlIPE);
  t->SetBranchAddress("masskstar", &masskstar);
  t->SetBranchAddress("massJ", &massJ);
  t->SetBranchAddress("Jpsipt", &Jpsipt);
  t->SetBranchAddress("Jrapidity", &Jrapidity);
  t->SetBranchAddress("Jphi", &Jphi);
  t->SetBranchAddress("sigLxyJ", &sigLxyJ);
  t->SetBranchAddress("cosalfaJ", &cosalfaJ);
  t->SetBranchAddress("sigLxyB", &sigLxyB);
  t->SetBranchAddress("cosalfaB", &cosalfaB);
  t->SetBranchAddress("pi1charge", &pi1charge);
  t->SetBranchAddress("pi2charge", &pi2charge);
  t->SetBranchAddress("mu1pt", &mu1pt);
  t->SetBranchAddress("mu2pt", &mu2pt);
  t->SetBranchAddress("mu1phi", &mu1phi);
  t->SetBranchAddress("mu2phi", &mu2phi);
  t->SetBranchAddress("mu1eta", &mu1eta);
  t->SetBranchAddress("mu2eta", &mu2eta);
  t->SetBranchAddress("pi1pt", &pi1pt);
  t->SetBranchAddress("pi1eta", &pi1eta);
  t->SetBranchAddress("pi1phi", &pi1phi);
  t->SetBranchAddress("pi2pt", &pi2pt);
  t->SetBranchAddress("pi2eta", &pi2eta);
  t->SetBranchAddress("pi2phi", &pi2phi);
  t->SetBranchAddress("Bpro", &Bpro);
  t->SetBranchAddress("Jpro", &Jpro);
  t->SetBranchAddress("dxypi1", &dxypi1);
  t->SetBranchAddress("dxyEpi1", &dxyEpi1);
  t->SetBranchAddress("dxysigpi1", &dxysigpi1);
  t->SetBranchAddress("DRmumu", &DRmumu);
  t->SetBranchAddress("nPV", &nPV);
  t->SetBranchAddress("runn", &runn);
  t->SetBranchAddress("evtn", &evtn);
  t->SetBranchAddress("Lumiblock", &Lumiblock);
  t->SetBranchAddress("Triggerjk", &Triggerjk);
  t->SetBranchAddress("Triggerjkjk", &Triggerjkjk);
  t->SetBranchAddress("maxPtpi", &maxPtpi);
  t->SetBranchAddress("minPtpi", &minPtpi);
  t->SetBranchAddress("ptpipi", &ptpipi);
  t->SetBranchAddress("masspipi", &masspipi);
  t->SetBranchAddress("masskk", &masskk);

  RooRealVar Bm("B_mass","B_mass",5.1,5.5);

  RooDataSet *massDS = new RooDataSet("massDS", "massDS", RooArgSet(Bm));

  for (int i = 0; i < t->GetEntries(); i++) {
    t->GetEvent(i);
    //Mass windows cuts
    if(massB<5.1 || massB>5.5) continue;
    if(massJ<=3.0969-0.150 || massJ>=3.0969+0.150) continue;
    if(Bpt<20.0 || Bpt>23.0) continue;
    if(mu1pt<4.0) continue;
    if(mu2pt<4.0) continue;
    if( abs(mu1eta)>2.4)continue;
    if( abs(mu2eta)>2.4)continue;
    if( abs(Brapidity)>2.4) continue;

    Bm = massB;
    massDS->add(RooArgSet(Bm));
  }
  massDS->Print("V");
  //Fit

  //Signal
  RooRealVar mean("mean", "mean", 5.28, 5.1, 5.5);
  RooRealVar sigma("sigma", "sigma", .02, 0.0, 0.04);
  RooRealVar alpha("alpha", "alpha", 1, 0.0, 5);
  RooRealVar n("n", "n", 3, 0, 5);

  RooCBShape sig("sig", "sig", Bm, mean, sigma, alpha, n);
  //Bkg
  RooRealVar a0("a0", "a0", -1,1);
  RooRealVar a1("a1", "a1", -1,1);
  RooChebychev bkg("bkg", "bkg", Bm, RooArgList(a0, a1));

  //Num eventos
  RooRealVar nsig("nsig","nsig",1000,0.,148085);
  RooRealVar nbkg("nbkg","nbkg",2000,0.,148085);

  RooAddPdf model("model","sig+bkg", RooArgList(sig,bkg), RooArgList(nsig,nbkg));
  RooFitResult* result = model.fitTo(*massDS, Extended(kTRUE), Save(kTRUE), NumCPU(4));


  //Canvas 1
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
  c1->cd();
  TPad *pad = new TPad("p", "", 0.01,0.01,0.99,0.99);
  pad->SetLeftMargin(0.13);   
  pad->SetRightMargin(0.019);
  pad->SetBottomMargin(0.15);
  pad->SetTopMargin(0.09);
  pad->Draw();
  pad->cd();
  gPad->SetLogy();
  RooPlot* frame = Bm.frame(5.1, 5.5, 70);
  frame->SetTitle("");
  frame->SetXTitle("M(J/#psiK^{*0})(GeV)"); frame->GetXaxis()->CenterTitle();
  frame->SetYTitle("Events"); frame->GetYaxis()->CenterTitle();
  frame->SetTitleSize(0.06, "XY"); 
  frame->SetLabelSize(0.05,"XY");

  massDS->plotOn(frame, Name("data"), MarkerSize(1.0),XErrorSize(0));
  model.plotOn(frame, Components(sig), LineColor(kRed), LineWidth(2), Name("signal"));
  model.plotOn(frame, Components(bkg), LineColor(kBlack), LineWidth(2), Name("bkg"));
  model.plotOn(frame, LineColor(kBlue), LineWidth(1), Name("model"));

  TLegend *legMass = new TLegend(0.7,0.85,0.9,0.5);
  legMass->SetTextSize(0.08);
  legMass->SetFillColor(0);
  legMass->SetBorderSize(0);
  legMass->SetFillStyle(0);
  legMass->AddEntry(frame->findObject("data"), "Data", "pe");
  legMass->AddEntry(frame->findObject("model"), "Fit", "l");
  legMass->AddEntry(frame->findObject("signal"),"Signal","l");
  legMass->AddEntry(frame->findObject("bkg"),"Bkg","l");

  frame->SetMinimum(1E2);
  frame->SetMaximum(6E4);
  frame->Draw();
  legMass->Draw();

  auto mtext = new TLatex();
  mtext->SetTextSize(0.11);
  mtext->SetTextFont(42);
  mtext->DrawLatex(5.25, 9000, "B^{0}");

  c1->Draw();
  c1->SaveAs("massFit.png");

  //Canvas 2
  TCanvas *c = new TCanvas("c", "c", 1000, 700);
  c->cd();

  TPad *pad1 = new TPad("p1", "", 0.01,0.4,0.9903769, 0.99);
  pad1->SetLeftMargin(0.13);   
  pad1->SetRightMargin(0.019);
  pad1->SetBottomMargin(0.15);
  pad1->SetTopMargin(0.09);
  pad1->Draw();

  TPad *pad2 = new TPad("p2", "", 0.01,0.01,0.9903769,0.39);
  pad2->SetLeftMargin(0.13);
  pad2->SetRightMargin(0.019);
  pad2->SetBottomMargin(0.1);
  pad2->SetTopMargin(0.09);
  pad2->Draw();

  pad1->cd();
  gPad->SetLogy();
  frame->SetTitle("");
  frame->SetXTitle("M(J/#psiK^{*0})(GeV)"); frame->GetXaxis()->CenterTitle();
  frame->SetYTitle("Events"); frame->GetYaxis()->CenterTitle();
  frame->SetTitleSize(0.06, "XY"); 
  frame->SetLabelSize(0.05,"XY");

  massDS->plotOn(frame, Name("data"), MarkerSize(1.0),XErrorSize(0));
  model.plotOn(frame, Components(sig), LineColor(kRed), LineWidth(2), Name("signal"));
  model.plotOn(frame, Components(bkg), LineColor(kBlack), LineWidth(2), Name("bkg"));
  model.plotOn(frame, LineColor(kBlue), LineWidth(1), Name("model"));

  frame->SetMinimum(1E2);
  frame->SetMaximum(6E4);
  frame->Draw();
  legMass->Draw();

  mtext->SetTextSize(0.11);
  mtext->SetTextFont(42);
  mtext->DrawLatex(5.25, 9000, "B^{0}");

  RooHist* massPull = frame->pullHist();
  pad2->cd();
  RooPlot* frame2 = Bm.frame(Title("Pull Distribution")) ;
  frame2->addPlotable(massPull,"P") ;
  frame2->SetYTitle("(Data-Fit)/#sigma");
  frame2->SetTitleSize(0.05, "XY");
  frame2->GetYaxis()->CenterTitle();
  frame2->SetLabelSize(0.06,"XY");
  frame2->Draw();

  c->Draw();
  c->SaveAs("massFitPull.png");

  ofstream salida("fitData.txt");
  salida << a0.getVal() << " " << a0.getError() << " " << a1.getVal() << 
  " " << a1.getError() << " " << mean.getVal() << " " << mean.getError() <<
  " " << sigma.getVal() << " " << sigma.getError() << " " << alpha.getVal() << 
  " " << alpha.getError() << " " << n.getVal() << " " << n.getError() << 
  " " << nbkg.getVal() << " " << nbkg.getError() << 
  " " << nsig.getVal() << " " << nsig.getError();
  salida.close();
}