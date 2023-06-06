#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooVoigtian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

#include "RooCBShape.h"

#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"
#include "RooRandom.h"

using namespace RooFit;
using namespace std;

TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, RooDataSet *data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel, RooAddPdf sumgau, RooAddPdf bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth)
{

  //Double_t nbin = ((supM-infM)/0.010) + 1;
  Double_t nbin = ((supM-infM)/0.010);
  
  int H = 600;
  int W = 800;
  TCanvas *c1 = new TCanvas("c1","",50,50,W,H);
  c1->cd();
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.09);
  c1->SetBottomMargin(0.13); 
  
  RooPlot* Mframe = M.frame(infM,supM,nbin);
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
  //MassModel.plotOn(Mframe);
  MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  MassModel.plotOn(Mframe,Components(sumgau),LineColor(kRed),LineWidth(2),Name("Signal")); 
  MassModel.plotOn(Mframe,Components(bkg1),LineColor(kGreen),LineWidth(2),Name("bkg")); 
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel.plotOn(Mframe);
  
  Mframe->SetYTitle("Events / 10 MeV"); 
  Mframe->SetLabelSize(0.04,"XY");
  Mframe->SetTitleSize(0.05,"XY");
  Mframe->GetYaxis()->CenterTitle();   
  Mframe->GetXaxis()->CenterTitle();
  Mframe->GetYaxis()->SetNdivisions(505,1);
  Mframe->GetXaxis()->SetNdivisions(505,1);   
  Mframe->GetXaxis()->SetDecimals(1); 
  Mframe->SetTitleOffset(0.9,"X");
  Mframe->SetTitleOffset(1.1,"Y");
  Mframe->SetTitleSize(0.06,"XY");
  Mframe->SetMinimum(1.0); 
  Mframe->Draw();
  
  //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88); 
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
  leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
  leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
  leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
  leg->Draw();
  
  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;
  //Double_t G = width.getVal()*1000.0;
  //Double_t GE = width.getError()*1000.0;
  Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
  Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
  
  TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);
  legpar->SetTextSize(0.035);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B_{c}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B_{c}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();
  
  TLegend *legMass = new TLegend(0.64,0.57,0.83,0.65);
  legMass->SetTextFont(43); 
  legMass->SetTextSize(20);  
  legMass->SetFillColor(0); 
  legMass->SetBorderSize(0);
  legMass->SetFillStyle(0); 
  legMass->SetHeader(Form("%1.1f #leq p_{T}(B_{c}) < %1.1f GeV ",ptl,pth));
  legMass->Draw(); 
  
  //TLatex *   tex1 = new TLatex(0.92,0.926,"61.2 fb^{-1} (13 TeV)");
  TLatex *   tex1 = new TLatex(0.92,0.926,"ToyMC");
  
  tex1->SetNDC();
  tex1->SetTextAlign(31);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.05); 
  tex1->SetLineWidth(2);
  
  TLatex *tex2 = new TLatex(0.2,0.926,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.05); 
  tex2->SetLineWidth(2);
  
  TLatex *tex3 = new TLatex(0.29,0.926,"Preliminary");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.05); 
  tex3->SetLineWidth(2);
  
  tex1->Draw();  
  tex2->Draw();
  tex3->Draw();
  
  c1->Modified();
  gPad->Update();
  gPad->RedrawAxis();
  TLine l;
  l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
  return c1; 
 
}

void Bplus_MC2(Int_t ntotal=1, Int_t seed=5, Double_t ptl=20.0, Double_t pth=23.0)
{
    /*
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);
  */

    Double_t Mmin = 5.0; 
    Double_t Mmax = 5.6; 

    RooRealVar m("m","M(J/#psi K^{+}) (GeV)",Mmin,Mmax);
  //RooDataSet data("data","data",RooArgSet(M));
 
  //RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
  //RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;
  
    //---- Mass model -------
    Double_t supM = Mmax;
    Double_t infM = Mmin;
  
    // Signal model
    // Parametros
    RooRealVar mean("mean","mean of gaussians",5.279,5.2,5.32);
    RooRealVar sigma1("sigma1","width of gaussians",0.010,0.001,0.015);
    RooRealVar sigma2("sigma2","width of gaussians",0.020,0.015,0.05);

    // Gaussianas
    RooGaussian sig1("sig1","Signal component 1",m,mean,sigma1);
    RooGaussian sig2("sig2","Signal component 2",m,mean,sigma2);

    // Suma de las gaussianas
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.);
    RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),RooArgList(sig1frac));

    // Background model
    // ERF
    RooRealVar d0("d0","d0",0.035,0.01,0.085);                 
    RooRealVar d1("d1","d1",5.14,5.05,5.2);
    RooGenericPdf genpdf("genpdf","genpdf","(TMath::Erf((-m + d1)/d0)+1)",RooArgSet(d0,d1,m));
    
    // Exponencial
    RooRealVar c("c","c",-10.0,10.0);
    RooExponential exp("exp","Exp. Background",m,c);

    // Suma de la ERF y la exponencial
    RooRealVar bkg0frac("bkg0frac","fraction of component 0 in background",0.8,0.,1.);
    RooAddPdf bkg_("bkg_","Background_",RooArgList(genpdf,exp),RooArgList(bkg0frac));
    
    // Crystal Ball function
    RooRealVar mean2("mean2","mean of gaussian",5.355,5.3,5.4) ;
    RooRealVar sigma3("sigma3","width of gaussian",0.010,0.001,0.015) ;
    RooRealVar alpha("alpha", "alpha", -1, -10.0, 0.0);
    RooRealVar n("n", "n", 3, 0, 10);
    RooCBShape sig3("sig3", "Signal component 3", m, mean2, sigma3, alpha, n);

    // Sumas de la Crystal ball con el resto del background
    RooRealVar bkg1frac("bkg1frac","fraction of component 1 in background",0.8,0.,1.);
    RooAddPdf bkg("bkg","Background",RooArgList(bkg_,sig3),RooArgList(bkg1frac)) ;
    
    // Número de eventos
    RooRealVar nsig("nsig","number of signal events",510000,0.,800000) ;
    RooRealVar nbkg("nbkg","number of background events",255000,0.,800000) ;
  
    // Creamos el modelo y lo fiteamos con la data
    RooAddPdf  model("model","sig+bkg", RooArgList(sig,bkg), RooArgList(nsig,nbkg));

    //------------ Fit procedure -------------------
    
    //RooRandom::randomGenerator()->SetSeed(3);
    RooRandom::randomGenerator()->SetSeed(seed);
    // esto que sigue solo es para que no imprima en pantalla todo lo que hace del ajuste
    
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().getStream(1).removeTopic(Integration) ;  
    RooMsgService::instance().getStream(1).removeTopic(Minimization) ;  
    RooMsgService::instance().getStream(1).removeTopic(Fitting) ;  
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
    RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
    RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
    RooMsgService::instance().getStream(1).removeTopic(Eval) ;
    RooMsgService::instance().Print() ;
  
    Double_t Nspull,Nsdif,Nbpull,Mupull,Nsig,NsigE,Nbkg,NbkgE;
    //TFile *file = new TFile("ToyMC_Bc.root","RECREATE");
    //TFile *file = new TFile(Form("ToyMC_Bs_ptbins_%1i_%1.0f_%1.0f.root",seed,ptl,pth),"RECREATE");
    TFile *file = TFile::Open(Form("ToyMC_Bc_%1i_%1.0f_%1.0f.root",seed,ptl,pth),"RECREATE");
    TTree *tree = new TTree("tree","tree");
    file->cd();

    tree->Branch("Nspull",&Nspull);
    tree->Branch("Nsdif",&Nsdif);
    tree->Branch("Nbpull",&Nbpull);    
    tree->Branch("Mupull",&Mupull);
    
    tree->Branch("Nsig",&Nsig);
    tree->Branch("NsigE",&NsigE);
    
    tree->Branch("Nbkg",&Nbkg);
    tree->Branch("NbkgE",&NbkgE);

  //*****************************************
  //    Input values from DATA fit in pt bins
  //*****************************************
    Double_t nsi, nsie, nbi, nbie;
    Double_t mui, muie, w1i, w1ie, w2i, w2ie;
    Double_t d0i, d0ie, d1i, d1ie, ci, cie;
    Double_t mu2i, mu2ie, w3i, w3ie, ali, alie, ni, nie;
    
    //Double_t Ns, Nb, mean, width1, width2, d0, d1, a0;
    ifstream entrada("../Data/output_BuFit.txt");
    //entrada.is_open();
    if ( !entrada ) 
        {
        cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
        exit( 1 );
        }
    entrada >> nsi >> nsie >> nbi >> nbie >> mui >> muie >> w1i >> w1ie >> w2i >> w2ie >> d0i >> d0ie >> d1i >> d1ie >> ci >> cie >> mu2i >> mu2ie >> w3i >> w3ie >> ali >> alie >> ni >> nie;
    entrada.close();
    
    cout << "Ya leí la data" << endl;
    
    Int_t nruns=0;
    
    for(Int_t n_=0;n_<ntotal;n_++)
    {
    //Ns.setConstant(kTRUE);
    nsig.setVal(nsi);
    nsig.setError(nsie);
    
    nbkg.setVal(nbi);
    nbkg.setError(nbie);
        
    mean.setVal(mui);
    mean.setError(muie);
    
    sigma1.setVal(w1i);
    sigma1.setError(w1ie);
    
    sigma2.setVal(w2i);
    sigma2.setError(w2ie);

    d0.setVal(d0i);
    d0.setError(d0ie);
    
    d1.setVal(d1i);
    d1.setError(d1ie);

    c.setVal(ci);
    c.setError(cie);

    mean2.setVal(mu2i);
    mean2.setError(mu2ie);

    sigma3.setVal(w3i);
    sigma3.setError(w3ie);

    alpha.setVal(ali);
    alpha.setError(alie);

    n.setVal(ni);
    n.setError(nie);

    cout << "Ya almacené los datos" << endl;

    //RooRandom::randomGenerator()->SetSeed(n);
    RooDataSet *dataToy = model.generate(RooArgSet(m), Extended(kTRUE));
    //RooDataSet *dataToy = model.generate(RooArgSet(m),Ns.getVal() +  Nb.getVal() );
    //RooDataSet dataToy = model.generate(RooArgSet(m), Ns.getVal() + Nb.getVal() );
    
    RooFitResult* fitres = model.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
    //RooFitResult* fitres = model.fitTo(data,Extended(),Minos(kTRUE),Save(kTRUE), NumCPU(4));
    //dataToy->Print("v");
    //fitres->Print("v");

    Nspull = (nsig.getVal()-nsi)/nsig.getError();
    Nsdif = (nsig.getVal()-nsi);
    Nbpull = (nbkg.getVal()-nbi)/nbkg.getError();   
    Mupull = (mean.getVal()-mui)/mean.getError();

    Nsig = nsig.getVal();
    NsigE = nsig.getError();
    
    Nbkg = nbkg.getVal();
    NbkgE = nbkg.getError();

    cout << "Ya creé Mu, Ns y Nb" << endl;

    tree->Fill();

    cout << "Ya llené el árbol" << endl;

    nruns++;
    cout << nruns << endl;
    
    TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, dataToy, m, supM, infM, model, sig, bkg, nsig, nbkg, sigma1, sigma2, sig1frac, mean, ptl, pth); 
    canv_nominal->Print(Form("plots/mass_BcFit_ptbins_ToyMC_%1.0f_%1.0f.png",ptl,pth));
    //canv_nominal->Print(Form("plots/mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.pdf",ptl,pth));
    
    cout << "Ya hice el plot" << endl;

    delete dataToy;
    delete fitres;
}

tree->Write();

cout<<"Ya escribí el arbol" << endl;

//file->Write("",TObject::kOverwrite);

cout<<" for end: "<<nruns<<endl;

}//End analysis