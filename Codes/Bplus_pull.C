#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
//#include "RooArgusBG.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"

using namespace RooFit;

void save_resultTfit(ofstream& salida, RooRealVar Ns, RooRealVar Nb, RooRealVar mean,  RooRealVar width, RooRealVar width2, RooRealVar d0, RooRealVar d1, RooRealVar c, RooRealVar mean2, RooRealVar sigma3, RooRealVar alpha, RooRealVar n)
{
    salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() << " " << mean.getVal() << " " << mean.getError() << " " <<  width.getVal() << " " << width.getError() << " " <<  width2.getVal() << " " << width2.getError() << " " <<  d0.getVal() << " " << d0.getError() << " " << d1.getVal() << " " << d1.getError() << " " << c.getVal() << " " << c.getError() << " " << mean2.getVal() << " " << mean2.getError() << " " << sigma3.getVal() << " " << sigma3.getError() << " " << alpha.getVal() << " " << alpha.getError() << " " << n.getVal() << " " << n.getError();
    cout << " el archivo se escribio bien" << endl; 
    return;
}

TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel, RooAddPdf sumgau, RooAddPdf bkg, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth) 
{
 Double_t nbin = ((supM-infM)/0.0055);

 int H = 670;
 int W = 450;
 TCanvas *c1 = new TCanvas("c1","",50,50,W,H);
 c1->cd();
 c1->SetLeftMargin(0.13);
 //c1->SetRightMargin(0.02);
 c1->SetRightMargin(0.03);
 c1->SetTopMargin(0.09);
 c1->SetBottomMargin(0.13);
 
 gPad->SetLogy(); 
 
 RooPlot* Mframe = M.frame(infM,supM,nbin);
 data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(0.6),XErrorSize(0));
 MassModel.plotOn(Mframe,DrawOption("l"),FillColor(0),LineWidth(1),Name("fittotal"));
 MassModel.plotOn(Mframe,Components(sumgau),LineColor(kRed),LineWidth(1),Name("Signal")); 
 MassModel.plotOn(Mframe,Components(bkg),LineColor(kBlack),LineWidth(1),Name("bkg")); 
 data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(0.8),XErrorSize(0),Name("Data"));
 MassModel.plotOn(Mframe);
 
 Mframe->SetTitle(" ");
 Mframe->SetYTitle("Events / 5 MeV");    
 Mframe->SetLabelSize(0.03,"XY");
 Mframe->SetTitleSize(0.04,"XY");
 Mframe->GetYaxis()->CenterTitle();   
 Mframe->GetXaxis()->CenterTitle();
 Mframe->GetYaxis()->SetNdivisions(505,1);
 Mframe->GetXaxis()->SetNdivisions(505,1);   
 Mframe->GetXaxis()->SetDecimals(1); 
 Mframe->SetTitleOffset(0.9,"X");
 Mframe->SetTitleOffset(1.1,"Y");
 Mframe->SetTitleSize(0.04,"XY");
 //Mframe->SetMinimum(10e2);
 //Mframe->SetMaximum(10e5);    
 Mframe->Draw();

 TLegend *leg = new TLegend(0.68,0.75,0.88,0.9); 
 leg->SetTextSize(0.038);
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->SetFillStyle(0);
 leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
 leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
 leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
 leg->AddEntry(Mframe->findObject("bkg"),"Comb. bkg.","l");
 leg->Draw();
 
 Double_t Mpsi = mean.getVal()*1000.0;
 Double_t MpsiE = mean.getError()*1000.0;
 Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
 Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
 

 TLegend *legMass = new TLegend(0.39,0.67,0.4,0.75);
 legMass->SetTextFont(47); 
 legMass->SetTextSize(38);  
 legMass->SetFillColor(0); 
 legMass->SetBorderSize(0);
 legMass->SetFillStyle(0); 
 legMass->SetHeader(Form("B^{+}"));
 legMass->Draw(); 
 
/*
 TLatex *tex2 = new TLatex(0.2,0.826,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.05); 
 tex2->SetLineWidth(2);
 tex2->Draw();
*/
     auto legend = new TLegend(.15,.8,.2,.9);
    //legend->AddEntry((TObject*)nullptr, " #bf{CMS}", "");
    legend->AddEntry((TObject*)nullptr, " L = 61.6 fb^{-1}", "");
    legend->AddEntry((TObject*)nullptr, "#sqrt{s} = 13 TeV", "");
    legend->SetTextSize(0.045);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

 c1->Modified();
 gPad->Update();
 gPad->RedrawAxis();
 return c1; 
 
}

TCanvas* CreateCanvasNomPull(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel,RooGaussian Sig2,RooGaussian Sig, RooAddPdf sumgau, RooExponential exp, RooGenericPdf genpdf, RooCBShape Sig3, RooAddPdf bkg, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth) 
{
  Double_t nbin = ((supM-infM)/0.0055);

  int H = 920;
  int W = 550;

  TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
  c1->cd() ;  
  c1->SetLeftMargin(0.005);
  c1->SetRightMargin(0.01);
  c1->SetTopMargin(0.09);
  c1->SetBottomMargin(0.1);

  TPad *pad1 = new TPad("pad1", "padi",0.01,0.382,0.9903769, 0.99 );
  pad1->SetLeftMargin(0.09);
  pad1->SetRightMargin(0.019);
  pad1->SetTopMargin(0.09);
  pad1->SetBottomMargin(0.0);

  TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.381);
  pad2->SetLeftMargin(0.09);
  pad2->SetRightMargin(0.019);  
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.25);
  pad2->SetFillColor(0);
  pad2->SetGridx(0);
  pad2->SetGridy(0);

  pad1->Draw();
  pad2->Draw();
  pad1->cd(); 
  gPad->SetLogy();

  RooPlot* Mframe = M.frame(infM,supM,nbin);
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
  MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  
  RooHist* hpullm2 = Mframe->pullHist() ;
  
  MassModel.plotOn(Mframe,Components(Sig),LineColor(kOrange-3),LineStyle(kDashed),Name("Sig_{1}")); 

  MassModel.plotOn(Mframe,Components(sumgau),LineColor(kRed),LineWidth(2),Name("Signal Total")); 
  MassModel.plotOn(Mframe,Components(Sig2),LineColor(kGreen),LineStyle(kDashed),Name("Sig_{2}")); 
  MassModel.plotOn(Mframe,Components(bkg),LineColor(kBlack),LineWidth(2),Name("bkg"));
  MassModel.plotOn(Mframe,Components(exp),LineColor(kBlue + 3),LineStyle(kDashed),Name("exp")); 
  MassModel.plotOn(Mframe,Components(genpdf),LineColor(kMagenta),LineStyle(kDashed),Name("ERF"));
  MassModel.plotOn(Mframe,Components(Sig3),LineColor(kGreen + 3),LineStyle(kDashed),Name("Sig_{3}"));
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel.plotOn(Mframe);
  
 Mframe->SetTitle(" ");
 Mframe->SetYTitle("Events / 5 MeV");    
 //Mframe->SetLabelSize(0.03,"XY");
 Mframe->SetLabelSize(0.0375,"XY");
 Mframe->GetYaxis()->CenterTitle();   
 Mframe->GetXaxis()->CenterTitle();
 Mframe->GetYaxis()->SetNdivisions(505,1);
 Mframe->GetXaxis()->SetNdivisions(505,1);   
 Mframe->GetXaxis()->SetDecimals(1); 
 Mframe->SetTitleOffset(0.9,"X");
 //Mframe->SetTitleOffset(1.1,"Y");
 Mframe->SetTitleOffset(0.8,"Y");
 //Mframe->SetTitleSize(0.04,"XY");
 Mframe->SetTitleSize(0.05,"XY");
 Mframe->SetMinimum(1e2);
 Mframe->SetMaximum(1e5);    
 Mframe->Draw();
  

 TLegend *leg = new TLegend(0.68,0.5,0.97,0.88); 
 leg->SetTextSize(0.038);
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->SetFillStyle(0);
 leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
 leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
 leg->AddEntry(Mframe->findObject("Signal Total")," Signal Total","l");
 leg->AddEntry(Mframe->findObject("Sig_{1}")," Sig_{1}","l");
 leg->AddEntry(Mframe->findObject("Sig_{2}")," Sig_{2}","l");
 leg->AddEntry(Mframe->findObject("bkg")," Comb. bkg.","l");
 leg->AddEntry(Mframe->findObject("exp")," Exp.","l");
 leg->AddEntry(Mframe->findObject("ERF")," B^{+} #rightarrow J/#psi #pi^{+}","l");
 leg->AddEntry(Mframe->findObject("Sig_{3}")," B #rightarrow J/#psi K^{+} X","l");
 leg->Draw();

    auto legend = new TLegend(.1,.75,.15,.89);
    //legend->AddEntry((TObject*)nullptr, " #bf{CMS}", "");
    legend->AddEntry((TObject*)nullptr, " L = 61.6 fb^{-1}", "");
    legend->AddEntry((TObject*)nullptr, "#sqrt{s} = 13 TeV", "");
    legend->SetTextSize(0.045);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();
  
  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;
  Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
  Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
  
  //TLegend *legpar = new TLegend(0.05,0.43,0.3,0.55);
  TLegend *legpar = new TLegend(0.65,0.35,0.9,0.5);
  legpar->SetTextSize(0.03);
  legpar->SetTextFont(42);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B^{+}) = %1.1f #pm %1.1f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.1f #pm %1.1f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{sig} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();

 TLegend *legMass = new TLegend(0.39,0.67,0.4,0.75);
 legMass->SetTextFont(47); 
 legMass->SetTextSize(38);  
 legMass->SetFillColor(0); 
 legMass->SetBorderSize(0);
 legMass->SetFillStyle(0); 
 legMass->SetHeader(Form("B^{+}"));
 legMass->Draw(); 
 
/*
 TLatex *tex2 = new TLatex(0.2,0.826,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.05); 
 tex2->SetLineWidth(2);
  
 tex2->Draw();
*/
  
  pad2->cd();
  
  // Create a new frame to draw the pull distribution 
  RooPlot* framem2 = M.frame(infM,supM,nbin) ;
  framem2->addPlotable(hpullm2,"P") ;
  
  framem2->SetTitle(" ");
  framem2->SetYTitle(" (Data-Fit)/#sigma");
  framem2->SetLabelSize(0.075,"XY");
  framem2->SetTitleSize(0.08,"XY");
  framem2->GetYaxis()->CenterTitle();   
  framem2->GetXaxis()->CenterTitle();
  framem2->GetYaxis()->SetNdivisions(505,1);
  framem2->GetXaxis()->SetNdivisions(505,1);   
  framem2->GetXaxis()->SetDecimals(1); 
  framem2->SetTitleOffset(0.9,"X");
  framem2->SetTitleOffset(0.5,"Y");
  framem2->SetMaximum(4.9);
  framem2->SetMinimum(-4.9);
  
  framem2->Draw();

  TLine *line1 = new TLine(infM,0.0,supM,0.0);
  line1->SetLineColor(2);
  line1->SetLineStyle(kDashed);
  line1->SetLineWidth(1);
  line1->Draw();
  
  c1->Modified();
  gPad->Update();
  gPad->RedrawAxis();
  return c1; 
}

void Bplus_pull(Double_t ptl=20.0, Double_t pth=23.0){
    // DataSet
    Double_t B_mass, J_mass, mu1pt, mu2pt, bcpt, mu1eta, mu2eta, Brapidity, dxysig1;
    UInt_t trijk;

    TChain *ch = new TChain("bstree","");
    ch->Add("../rootfiles/reducetree_Bujk_AOD_3_best1.root/butree");

    TTree *tree = (TTree*) ch;
    Long64_t nentries = tree->GetEntries();
    //Long64_t nentries = 1000;

    tree->SetBranchAddress ("B_mass", &B_mass);
    tree->SetBranchAddress ("J_mass", &J_mass);
    tree->SetBranchAddress ("mu1pt", &mu1pt);
    tree->SetBranchAddress ("mu2pt", &mu2pt);
    tree->SetBranchAddress ("etamu1", &mu1eta);
    tree->SetBranchAddress ("etamu2", &mu2eta);
    tree->SetBranchAddress ("rapidityB", &Brapidity);
    tree->SetBranchAddress ("trijk", &trijk);
    tree->SetBranchAddress ("dxysig1", &dxysig1);

    Double_t Mmin = 5.0;
    Double_t Mmax = 5.6;

    // Declare observable m
    RooRealVar m("m","M(J/#psi K^{+}) (GeV)",Mmin,Mmax);
    //RooDataSet *data = new RooDataSet("data", "data", RooArgSet(m));
    RooDataSet data("data", "data", RooArgSet(m));

    Int_t nTen = nentries/10;
    Int_t k=0;
    Int_t nbytes = 0, nb = 0;

    for(Long64_t jentry=0; jentry<nentries;jentry++)
        {
        Long64_t ientry = tree->LoadTree(jentry);
        if (ientry < 0) break;
        nb = tree->GetEntry(jentry);   nbytes += nb;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        //Mass windows cuts
        if(B_mass<=Mmin || B_mass>=Mmax) continue;
        if(J_mass<=3.0969-0.150 || J_mass>=3.0969+0.150) continue;
        if(mu1pt<4.0) continue;
        if(mu2pt<4.0) continue;
        if(bcpt<ptl || bcpt>=pth)continue; 

        if( abs(mu1eta)>2.4)continue;
        if( abs(mu2eta)>2.4)continue;
        if( abs(Brapidity)>2.4) continue;
        if(trijk>10.0 || trijk<1.0)continue;
        //if(trijk>10.0)continue;
        if(abs(dxysig1)<2.0)continue;
        m=B_mass;    
        data.add(RooArgSet(m));
        }

    // Signal model
    // Parameters for the Gaussian signal
    RooRealVar mean("mean","mean of gaussians",5.279,5.2,5.32) ;
    RooRealVar sigma1("sigma1","width of gaussians",0.010,0.001,0.015) ;
    RooRealVar sigma2("sigma2","width of gaussians",0.020,0.015,0.05) ;
    //RooRealVar alpha("alpha", "alpha", 1, 0.0, 5);
    //RooRealVar n("n", "n", 3, 0, 5);

    // Gaussians
    //RooCBShape sig1("sig1", "Signal component 1", m, mean, sigma1, alpha, n);
    RooGaussian sig1("sig1","Signal component 1",m,mean,sigma1) ;
    //RooCBShape sig2("sig2", "Signal component 2", m, mean, sigma2, alpha, n);
    RooGaussian sig2("sig2","Signal component 2",m,mean,sigma2) ;

    // Sum of both gaussians
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.);
    //RooRealVar sig0frac("sig0frac","fraction of component 0 in signal",0.8,0.,1.);
    //RooAddPdf sig_("sig_","Signal_",RooArgList(sig1,sig2),RooArgList(sig0frac));
    RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),RooArgList(sig1frac));

    /*
    //RooRealVar mean2("mean2","mean of gaussian",5.355,5.3,5.4) ;
    RooRealVar mean2("mean2","mean of gaussian",5.395,5.385,5.415);
    RooRealVar sigma3("sigma3","width of gaussian",0.010,0.001,0.015) ;
    RooRealVar alpha("alpha", "alpha", 1, 0.0, 5);
    RooRealVar n("n", "n", 3, 0, 5);
    
    RooCBShape sig3("sig3", "Signal component 3", m, mean, sigma3, alpha, n);

    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.);
    RooAddPdf sig("sig","Signal_",RooArgList(sig_,sig3),RooArgList(sig1frac));
    */
    // Background model
    // ERF
    RooRealVar d0("d0","d0",0.035,0.01,0.085);// another parameters related to "shoulder" shape on B+ mass BG                                                  
    //RooRealVar d0("d0","d0",0.05,0.01,0.105);
    RooRealVar d1("d1","d1",5.14,5.05,5.2);// describes the "shoulder" steepness
    RooGenericPdf genpdf("genpdf","genpdf","(TMath::Erf((-m + d1)/d0)+1)",RooArgSet(d0,d1,m));
    
    // Chebychev polynomial.
    //RooRealVar a0("a0","a0",-0.5,-1.,1.);
    //RooChebychev chev("chev","Chebychev polynomial degree 1",m,a0) ;

    // Exponential
    RooRealVar c("c","c",-10.0,10.0);
    RooExponential exp("exp","Exp. Background",m,c);

    //RooRealVar bkg1frac("bkg1frac","fraction of component 1 in background",0.8,0.,1.) ;
    RooRealVar bkg0frac("bkg0frac","fraction of component 0 in background",0.8,0.,1.);
    //RooAddPdf bkg("bkg","Background",RooArgList(genpdf,chev),RooArgList(bkg1frac));
    //RooAddPdf bkg_("bkg_","Background_",RooArgList(genpdf,chev),RooArgList(bkg0frac)) ;
    RooAddPdf bkg_("bkg_","Background_",RooArgList(genpdf,exp),RooArgList(bkg0frac));
    //RooAddPdf bkg("bkg","Background",RooArgList(genpdf,exp),RooArgList(bkg1frac));

    
    RooRealVar mean2("mean2","mean of gaussian",5.355,5.3,5.4) ;
    //RooRealVar mean2("mean2","mean of gaussian",5.45,5.35,5.55);
    RooRealVar sigma3("sigma3","width of gaussian",0.010,0.001,0.015) ;
    //RooRealVar alpha("alpha", "alpha", 1, 0.0, 10.0);
    RooRealVar alpha("alpha", "alpha", -1, -10.0, 0.0);
    RooRealVar n("n", "n", 3, 0, 10);
    
    RooCBShape sig3("sig3", "Signal component 3", m, mean2, sigma3, alpha, n);

    RooRealVar bkg1frac("bkg1frac","fraction of component 1 in background",0.8,0.,1.);
    RooAddPdf bkg("bkg","Background",RooArgList(bkg_,sig3),RooArgList(bkg1frac)) ;
    
    /*
    //RooRealVar mean2("mean2","mean of gaussian",5.355,5.3,5.4) ;
    RooRealVar mean2("mean2","mean of gaussian",5.45,5.35,5.55);
    RooRealVar sigma3("sigma3","width of gaussian",0.010,0.001,0.015) ;
    RooRealVar alpha("alpha", "alpha", 1, 0.0, 5);
    RooRealVar n("n", "n", 3, 0, 5);
    
    RooCBShape sig3("sig3", "Signal component 3", m, mean, sigma3, alpha, n);
    */
    // Number of events
    RooRealVar nsig("nsig","number of signal events",510000,0.,800000) ;
    RooRealVar nbkg("nbkg","number of background events",255000,0.,800000) ;
    //RooRealVar ndec("ndec","number of decayment events",80000,0.,800000) ;
  
    // Use AddPdf to extend the model. Giving as many coefficients as pdfs switches on extension.
    RooAddPdf  model("model","sig+bkg", RooArgList(sig,bkg), RooArgList(nsig,nbkg));
    RooFitResult* result = model.fitTo(data, Extended(kTRUE), Save(kTRUE));
    //data->Print("V");

    // Fit model to data
    //model.fitTo(*data, Extended(kTRUE)) ;
    /*
    //Graph
    TCanvas *canvas = new TCanvas("Mass B_u", "Mass B_u",50,50, 800, 600);
    canvas->SetLeftMargin(0.09);   
    canvas->SetRightMargin(0.019);
    canvas->SetTopMargin(0.09);
    canvas->SetBottomMargin(0.09);

    // First pad
    TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();

    // Second pad
    TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    //pad2->SetFillColor(0);
    //pad2->SetGridx(0);
    //pad2->SetGridy(0);
    pad2->Draw();

    pad1->cd();
    gPad->SetLogy();
    */
    RooPlot* frame = m.frame(Title("Mass Distribution"));
    data.plotOn(frame);//,DataError(RooAbsData::SumW2),MarkerSize(.65),XErrorSize(0));
    model.plotOn(frame,LineColor(kBlue),LineWidth(2),Name("Fit")); 
    model.plotOn(frame,Components(bkg),LineColor(kBlack),LineWidth(2),Name("bkg")); 
    model.plotOn(frame,Components(sig),LineColor(kRed),LineWidth(2),Name("Signal")); 

    // Filling pad1
    /*
    frame->SetTitle("Fit to the B_{u} mass");
    frame->SetTitleSize(0.8,"");
    frame->SetYTitle("Events");
    frame->SetLabelSize(0.06,"XY");
    frame->SetTitleSize(0.19,"XY");
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->CenterTitle();
    //frame->GetYaxis()->SetNdivisions(505,1);
    //frame->GetXaxis()->SetNdivisions(505,1);
    //frame->GetXaxis()->SetDecimals(1);
    //frame->GetXaxis()->SetTickLength(0.06); 
    frame->SetTitleSize(0.11,"X");
    frame->SetTitleSize(0.1,"Y");
    frame->SetTitleOffset(0.9,"X");
    frame->SetTitleOffset(0.42,"Y");
    frame->SetMinimum(-2.0);
    */
    /*
    frame->Draw();

    TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
    leg->SetTextSize(0.06);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(frame->findObject("Data")," Data","ep"); 
    leg->AddEntry(frame->findObject("Fit")," Fit result","l");
    leg->AddEntry(frame->findObject("Signal")," Signal","l");
    leg->AddEntry(frame->findObject("bkg"),"Comb. backg.","l");
    leg->Draw();

    Double_t Mpsi = mean.getVal()*1000.0;
    Double_t MpsiE = mean.getError()*1000.0;
    Double_t G = sqrt( sig1frac.getVal()*sigma1.getVal()*sigma1.getVal() + (1-sig1frac.getVal())*sigma2.getVal()*sigma2.getVal() )*1000.0;
    Double_t GE = (1/G)*sqrt( (sig1frac.getVal()*sig1frac.getVal())*(sigma1.getVal()*sigma1.getVal())*(sigma1.getError()*sigma1.getError()) + ((1-sig1frac.getVal())*(1-sig1frac.getVal()))*(sigma2.getVal()*sigma2.getVal())*(sigma2.getError()*sigma2.getError()) )*1000.0*1000.0;

    TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
    legpar->SetTextSize(0.06);
    legpar->SetTextFont(42);
    legpar->SetFillColor(0);
    legpar->SetBorderSize(0);
    legpar->SetFillStyle(0);
    legpar->AddEntry("",Form("M(B_{u}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
    legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
    legpar->AddEntry("",Form("N_{sig} = %1.0f #pm %1.0f",nsig.getVal(),nsig.getError()),"");
    legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",nbkg.getVal(),nbkg.getError()),"");
    legpar->Draw();

    auto legend = new TLegend(.2,.25,.4,.45);
    legend->AddEntry((TObject*)nullptr, " #bf{CMS}", "");
    legend->AddEntry((TObject*)nullptr, " L = 61.6 fb^{-1}", "");
    legend->AddEntry((TObject*)nullptr, "#sqrt{s} = 13 TeV", "");
    legend->SetTextSize(0.06);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->Draw();

    pad2->cd();

    RooHist* masspull = frame->pullHist();
    RooPlot* framem2 = m.frame(Title("Pull Distribution"));

    framem2->addPlotable(masspull,"P") ;
    framem2->SetTitle("");
    framem2->SetYTitle(" (N^{fit}-N^{true})/#sigma_{fit}");
    framem2->SetXTitle("Mass (Gev)");
    framem2->SetLabelSize(0.1,"XY");
    framem2->SetTitleSize(0.11,"X");
    framem2->SetTitleSize(0.135,"Y");  
    framem2->GetYaxis()->CenterTitle();   
    framem2->GetXaxis()->CenterTitle();
    framem2->GetYaxis()->SetNdivisions(505,1);
    framem2->GetXaxis()->SetNdivisions(505,1);
    framem2->GetXaxis()->SetTickLength(0.07);   
    framem2->SetTitleOffset(0.9,"X");
    framem2->SetTitleOffset(0.31,"Y");
    //framem2->SetMaximum(4.9);
    //framem2->SetMinimum(-4.9);
    framem2->Draw();

    canvas->Draw();
    //canvas->SetLogy(true);
    //canvas->SaveAs("fit_Bplus.png");
    canvas->SaveAs("fit.png");
    */
    /*
    ofstream salida_TotalFit("Data/output_BuFit.txt");
    salida_TotalFit.is_open();
    //save_resultTfit(salida_TotalFit, nsig, nbkg, mean, sigma1, sigma2, d0, d1, a0);
    save_resultTfit(salida_TotalFit, nsig, nbkg, mean, sigma1, sigma2, d0, d1, c, mean2, sigma3, alpha, n);
    
    TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", result, data, m, Mmax, Mmin, model, sig, bkg, nsig, nbkg, sigma1, sigma2, sig1frac, mean, ptl, pth);
    canv_nominal->Print(Form("Plots_fit/mass_BuFit_ptbins_%1.0f_%1.0f_NOpull.png",ptl,pth));
    //canv_nominal->Print(Form("plots_ptandeta/mass_BsFit_ptbins_%1.0f_%1.0f_NOpull.pdf",ptl,pth));
    TCanvas* canv_nominalpull = CreateCanvasNomPull("canv_nominalpull", result, data, m, Mmax, Mmin, model, sig2, sig1, sig, exp, genpdf, sig3, bkg, nsig, nbkg, sigma1, sigma2, sig1frac, mean, ptl, pth);  
    canv_nominalpull->Print(Form("Plots_fit/mass_BuFitpull_ptbins_%1.0f_%1.0f.png",ptl,pth));
    //canv_nominalpull->Print(Form("plots_ptandeta/mass_BsFitpull_ptbins_%1.0f_%1.0f.pdf",ptl,pth));
    */
}