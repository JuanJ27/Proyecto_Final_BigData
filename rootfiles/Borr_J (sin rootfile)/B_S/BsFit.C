#include "classreduce.C"
#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
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

using namespace RooFit;
using namespace std;


void save_resultTfit(ofstream& salida, RooRealVar Ns, RooRealVar Nb, RooRealVar a0, RooRealVar fs, RooRealVar mean,  RooRealVar width, RooRealVar width2)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() <<" " <<  a0.getVal() << " " << a0.getError() << " "  
	 << fs.getVal() << " " << fs.getError() << " " << mean.getVal() << " " << mean.getError() << " "
	 <<  width.getVal() << " " << width.getError() << " " <<  width2.getVal() << " " << width2.getError()
    ;
  cout << " el archivo se escribio bien" << endl; 
  return;
}


void save_result(ofstream& salida, RooRealVar Ns, RooRealVar Nb)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError();
  cout << " el archivo se escribio bien" << endl; 
  return;
}

TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel, RooAddPdf sumgau, RooChebychev bkg, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth) 
{
 Double_t nbin = ((supM-infM)/0.0055);

 int H = 670;
 int W = 450;
 TCanvas *c1 = new TCanvas("c1","",50,50,W,H);
 c1->cd();
 c1->SetLeftMargin(0.13);
 c1->SetRightMargin(0.02);
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
 Mframe->SetMinimum(2.0);
 Mframe->SetMaximum(1390.1);    
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
 legMass->SetHeader(Form("B_{s}^{0}"));
 legMass->Draw(); 
 

 TLatex *tex2 = new TLatex(0.2,0.826,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.05); 
 tex2->SetLineWidth(2);
  
 tex2->Draw();


 c1->Modified();
 gPad->Update();
 gPad->RedrawAxis();
 return c1; 
 
}

TCanvas* CreateCanvasNomPull(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel,RooGaussian Sig2,RooGaussian Sig, RooAddPdf sumgau, RooChebychev bkg, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth) 
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
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel.plotOn(Mframe);
  
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
 Mframe->SetMinimum(2.0);
 Mframe->SetMaximum(1390.1);    
 Mframe->Draw();
  

 TLegend *leg = new TLegend(0.68,0.55,0.97,0.88); 
 leg->SetTextSize(0.038);
 leg->SetFillColor(0);
 leg->SetBorderSize(0);
 leg->SetFillStyle(0);
 leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
 leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
 leg->AddEntry(Mframe->findObject("Signal Total")," Signal Total","l");
 leg->AddEntry(Mframe->findObject("Sig_{1}")," Sig_{1}","l");
 leg->AddEntry(Mframe->findObject("Sig_{2}")," Sig_{2}","l");
 leg->AddEntry(Mframe->findObject("bkg"),"Comb. bkg.","l");
 leg->Draw();
  
  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;
  Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
  Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
  
  TLegend *legpar = new TLegend(0.05,0.43,0.3,0.65);
  legpar->SetTextSize(0.03);
  legpar->SetTextFont(42);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B_{s}^{0}) = %1.1f #pm %1.1f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.1f #pm %1.1f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B_{s}^{0}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();

 TLegend *legMass = new TLegend(0.39,0.67,0.4,0.75);
 legMass->SetTextFont(47); 
 legMass->SetTextSize(38);  
 legMass->SetFillColor(0); 
 legMass->SetBorderSize(0);
 legMass->SetFillStyle(0); 
 legMass->SetHeader(Form("B_{s}^{0}"));
 legMass->Draw(); 
 

 TLatex *tex2 = new TLatex(0.2,0.826,"CMS");
 tex2->SetNDC();
 tex2->SetTextFont(61);
 tex2->SetTextSize(0.05); 
 tex2->SetLineWidth(2);
  
 tex2->Draw();

  
  pad2->cd();
  
  // Create a new frame to draw the pull distribution 
  RooPlot* framem2 = M.frame(infM,supM,nbin) ;
  framem2->addPlotable(hpullm2,"P") ;
  
  framem2->SetYTitle(" (Data-Fit)/#sigma");
  framem2->SetLabelSize(0.03,"XY");
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


void BsFit(Double_t ptl=20.0, Double_t pth=23.0)
{
gStyle->SetOptTitle(0);
gStyle->SetOptFit(0);
gStyle->SetOptStat(0);

Double_t Mmin = 5.240; 
Double_t Mmax = 5.490; 
Double_t Mmin0 = 5.350; 
 
 
TChain *ch = new TChain("bstree","");
ch->Add("reducetree_Bsphi_AOD_3_best1.root/bstree");

TTree *tree = (TTree*) ch;
classreduce t(tree);
Long64_t nentries = t.fChain->GetEntries();
cout<<" Entries : "<<nentries<<endl;
  
//------------------------------
RooRealVar M("M"," M(J/#psi #phi) (GeV)",Mmin,Mmax);
RooRealVar M0("M0"," M0(J/#psi #phi) (GeV)",Mmin0,Mmax);
RooDataSet data("data","data",RooArgSet(M));

Int_t nTen = nentries/10;
Int_t k=0;
Int_t nbytes = 0, nb = 0;
for(Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
      if(jentry==nentries-1) cout<<endl;

      //Mass windows cuts
      if(t.B_mass<=Mmin || t.B_mass>=Mmax) continue;	
      if(t.J_mass<=3.0969-0.150 || t.J_mass>=3.0969+0.150) continue;
      if(t.phi_mass<=1.01946-0.010 || t.phi_mass>=1.01946+0.010) continue;// pdgM = 1.019461 

      if(t.bcpt<ptl || t.bcpt>=pth)continue; 

      if(abs(t.rapidityB)>2.4)continue;
      
      if((t.pdl/t.pdle)<5.0)continue;
      if(t.pion1pt<1.2)continue;
      if(t.pion2pt<1.2)continue;

      if(t.trijk>10 || t.trijk<1)continue;
      if(abs(t.dxysig1)<2.0)continue;


      M=t.B_mass;    
      data.add(RooArgSet(M));

    }

 cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
 data.Print("v");
 
 //---- Mass model -------
 Double_t supM = Mmax;
 Double_t infM = Mmin;
 
 // ****define nominal background ****
 RooRealVar a0("a0","a0",0.23,-10.0,10.0);
 RooChebychev bkg("bkg","Background",M,RooArgList(a0)); 
 
 //gaussians       
 RooRealVar mean("mean"," Mass mean",5.366,5.320,5.400,"GeV");
 RooRealVar width("width"," Mass width",0.010,0.001,0.015,"GeV");
 RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

 RooRealVar width2("width2"," Mass width2 ",0.020,0.015,0.05,"GeV");
 RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

 //********final PDF ********
 RooRealVar Ns("Ns","Ns",0.,120000);
 RooRealVar Nb("Nb","Nb",0.,120000);   
 RooRealVar fs("fs","fs",0.8,0.,1.);

 RooAddPdf sumgau("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));
 //model
 RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,bkg),RooArgList(Ns,Nb));
 
 //------------ Fit procedure -------------------
 Ns.setVal(10000.0);
 Nb.setVal(10000.0);
 
 RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
 data.Print("v"); 
 fitres->Print("v");

 Double_t Gt = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
 Double_t GtE = (1/Gt)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;

 //parameters output file
 ofstream salida_TotalFit(Form("plots_ptandeta/output_TotalFit_BsFit_nominal_%1.0f_%1.0f.txt",ptl,pth));
 salida_TotalFit.is_open();
 save_resultTfit(salida_TotalFit, Ns, Nb, a0, fs, mean, width, width2);


 //parameters output file
 ofstream salida_nominal(Form("plots_ptandeta/output_BsFit_nominal_%1.0f_%1.0f.txt",ptl,pth));
 salida_nominal.is_open();
 save_result(salida_nominal, Ns, Nb);// OJO esta es la funciona a la que se le pasa el archivo por referencia
 
 //made canvas
 TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, data, M, supM, infM, MassModel, sumgau, bkg, Ns, Nb, width, width2, fs, mean,ptl,pth);
 canv_nominal->Print(Form("plots_ptandeta/mass_BsFit_ptbins_%1.0f_%1.0f_NOpull.png",ptl,pth));
 canv_nominal->Print(Form("plots_ptandeta/mass_BsFit_ptbins_%1.0f_%1.0f_NOpull.pdf",ptl,pth));
 TCanvas* canv_nominalpull = CreateCanvasNomPull("canv_nominalpull", fitres, data, M, supM, infM, MassModel,Sig2 ,Sig, sumgau, bkg, Ns, Nb, width, width2, fs, mean,ptl,pth);  
 canv_nominalpull->Print(Form("plots_ptandeta/mass_BsFitpull_ptbins_%1.0f_%1.0f.png",ptl,pth));
 canv_nominalpull->Print(Form("plots_ptandeta/mass_BsFitpull_ptbins_%1.0f_%1.0f.pdf",ptl,pth));
 
}//End analysis
