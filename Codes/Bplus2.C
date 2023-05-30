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

void Bplus2(){
    // DataSet
    Double_t B_mass;

    TChain *ch = new TChain("bstree","");
    ch->Add("../rootfiles/reducetree_Bujk_AOD_3_best1.root/butree");

    TTree *tree = (TTree*) ch;
    Long64_t nentries = tree->GetEntries();
    //Long64_t nentries = 1000;

    tree->SetBranchAddress ("B_mass", &B_mass);

    Double_t Mmin = 5.0; 
    Double_t Mmax = 5.6; 

    RooRealVar m("m","m",Mmin,Mmax);
    RooDataSet *data = new RooDataSet("data", "data", RooArgSet(m));
    // Create a Histogram
    //int num_bins = 100;
    //auto hist = new TH1F("Histogram", "Histogram Title", num_bins, min_value, max_value);

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
        m=B_mass;    
        data->add(RooArgSet(m));
        }
        /*
//    for (int i = 0; i < tree->GetEntries(); ++i) {
    for (int i = 0; i < 1000; ++i) {

        tree->GetEntry(i); // Output: 132
        hist->Fill(B_mass);
        m = B_mass;
        data->add(RooArgSet(m));  
    }
*/
    // Declare observable m
    

    // Signal model
    // Parameters for the Gaussian signal
    RooRealVar mean("mean","mean of gaussians",5.279,5.2,5.32) ;
    RooRealVar sigma1("sigma1","width of gaussians",0.010,0.001,0.015) ;
    RooRealVar sigma2("sigma2","width of gaussians",0.020,0.015,0.05) ;

    // Gaussians
    RooGaussian sig1("sig1","Signal component 1",m,mean,sigma1) ;
    RooGaussian sig2("sig2","Signal component 2",m,mean,sigma2) ;

    // Sum of both gaussians
    RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
    RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;

    // Background model
    // ERF
    /*
    RooRealVar mean_erf("mean_erf", "mean_erf", 5.15,5.1,5.2);
    RooRealVar sigma_erf("sigma_erf", "sigma_erf", 0.035, 0.01, 0.085);
    RooGenericPdf erf("erf", "Error Function", "TMath::Erf((m - mean_erf) / sigma_erf)+1", RooArgSet(m, mean_erf, sigma_erf));
    */
    RooRealVar d0("d0","d0",0.035,0.01,0.085);// another parameters related to "shoulder" shape on B+ mass BG                                                  
    RooRealVar d1("d1","d1",5.14,5.05,5.2);// describes the "shoulder" steepness
    RooGenericPdf genpdf("genpdf","genpdf","(TMath::Erf((-m + d1)/d0)+1)",RooArgSet(d0,d1,m));
    
    // Define the parameter of the Chebychev polynomial.
    RooRealVar a0("a0","a0",-0.5,-1.,1.);
    // Create the error function.
    RooChebychev chev("chev","Chebychev polynomial degree 1",m,RooArgSet(a0)) ;

    RooRealVar bkg1frac("bkg1frac","fraction of component 1 in background",0.8,0.,1.) ;
    RooAddPdf bkg("bkg","Background",RooArgList(genpdf,chev),bkg1frac) ;
    //RooAddPdf bkg("bkg","Background",RooArgList(erf,chev),bkg1frac) ;

    // Number of events
    RooRealVar nsig("nsig","number of signal events",255000,0.,800000) ;
    RooRealVar nbkg("nbkg","number of background events",510000,0.,800000) ;
  
    // Use AddPdf to extend the model. Giving as many coefficients as pdfs switches on extension.
    RooAddPdf  model("model","sig+bkg", RooArgList(sig,bkg), RooArgList(nsig,nbkg));

    //data->Print("V");

    // Fit model to data
    model.fitTo(*data, Extended(kTRUE)) ;

    RooPlot* frame = m.frame();
    data->plotOn(frame,DataError(RooAbsData::SumW2),MarkerSize(.65),XErrorSize(0));
    model.plotOn(frame,LineColor(kRed),LineWidth(2),Name("Fit")); 
    model.plotOn(frame,Components(bkg),LineColor(kGreen),LineWidth(2),Name("bkg")); 
    model.plotOn(frame,Components(sig),LineColor(kBlue),LineWidth(2),Name("Signal")); 
    //frame->Draw();
    RooHist* hpullm2 = frame->pullHist() ;


    //Graph

    TCanvas *canvas = new TCanvas("basic", "basic",50,50, 800, 600);
    canvas->SetLeftMargin(0.09);   
    canvas->SetRightMargin(0.019);
    canvas->SetTopMargin(0.09);
    canvas->SetBottomMargin(0.09); 
    TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);  

    TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
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

    frame->SetTitle("Fit to the B_{c} mass");
    frame->SetTitleSize(0.8,"");
    frame->SetYTitle("Events / 5 MeV");
    frame->SetXTitle("M(J/ψ K^{+}) (GeV)");
    frame->SetLabelSize(0.06,"XY");
    frame->SetTitleSize(0.19,"XY");
    frame->GetYaxis()->CenterTitle();   
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->SetNdivisions(505,1);
    frame->GetXaxis()->SetNdivisions(505,1);   
    frame->GetXaxis()->SetDecimals(1);
    frame->GetXaxis()->SetTickLength(0.03); 
    frame->SetTitleSize(0.04,"X");
    frame->SetTitleSize(0.04,"Y");  
    frame->SetTitleOffset(0.9,"X");
    frame->SetTitleOffset(0.42,"Y");
    frame->SetMinimum(-2.0);   
    frame->Draw();

    TLegend *leg = new TLegend(0.75,0.63,0.95,0.88);
    leg->SetTextSize(0.045);
    leg->SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(frame->findObject("Data")," Data","ep"); 
    leg->AddEntry(frame->findObject("Fit")," Fit result","l");
    leg->AddEntry(frame->findObject("Signal")," Signal","l");
    leg->AddEntry(frame->findObject("bkg"),"Comb. backg.","l");
    leg->Draw();

    pad2->cd();
    RooPlot* framem2 = m.frame() ;


    framem2->addPlotable(hpullm2,"P") ;
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
    framem2->SetMaximum(4.9);
    framem2->SetMinimum(-4.9);

    framem2->Draw();
    canvas->Modified();
    //canvas->SetLogy(true);
    canvas->SaveAs("fit_Bplus.png");
    // Close the input file
    //infile->Close();
}