#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;
void montecarlo(){
    Float_t a0f, a0fErr, a1f, a1fErr, meanf, meanfErr, sigmaf, sigmafErr, alphaf, alphafErr, nf, nfErr, nbkgf, nbkgfErr, nsigf, nsigfErr;

    ifstream entrada("fitData.txt");
    if(!entrada) 
    {
      cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
      exit( 1 );
    }
    entrada >> a0f >> a0fErr >> a1f >> a1fErr >> meanf >> meanfErr >> sigmaf >> sigmafErr >> alphaf >> alphafErr >> nf >> nfErr >> nbkgf >> nbkgfErr >> nsigf >> nsigfErr;
    entrada.close();

    RooRealVar m("m", "m", 5.1, 5.5);

    //Signal
    RooRealVar mean("mean", "mean", 5.1, 5.5);
    RooRealVar sigma("sigma", "sigma", 0.0, 0.04);
    RooRealVar alpha("alpha", "alpha", 0.0, 5);
    RooRealVar n("n", "n", 0, 5);
    //RooGaussian sig("sig", "sig", Bm, mean, sigma);
    RooCBShape sig("sig", "sig", m, mean, sigma, alpha, n);
    //Bkg
    RooRealVar a0("a0", "a0", -1,1);
    RooRealVar a1("a1", "a1", -1,1);
    RooChebychev bkg("bkg", "bkg", m, RooArgList(a0, a1));

    //Num eventos
    RooRealVar nsig("nsig","nsig",0,110000);
    RooRealVar nbkg("nbkg","nbkg",0,50000);
    
    //Pdf total
    RooAddPdf model("model","sig+bkg", RooArgList(sig,bkg), RooArgList(nsig,nbkg));

    //DataSets
    RooRealVar massPull("massPUll", "massPull", -6, 6);
    RooRealVar nsigPull("nsigPull", "nsigPull", -6, 6);
    RooRealVar nbkgPull("nbkgPull", "nbkgPull", -6, 6);

    RooDataSet mpDset("massPUll", "massPull", RooArgSet(massPull));
    RooDataSet nsigDset("nsigPULL", "nsigPULL", RooArgSet(nsigPull));
    RooDataSet nbkgDset("nbkgPULL", "nbkgPULL", RooArgSet(nbkgPull));

    RooRandom::randomGenerator()->SetSeed(84329746);
    for(int i=0; i<100; i++){

      mean.setVal(meanf); mean.setError(meanfErr);
      sigma.setVal(sigmaf); sigma.setError(sigmafErr);
      a0.setVal(a0f); a0.setError(a0fErr);
      nsig.setVal(nsigf); nsig.setError(nsigfErr);
      a1.setVal(a1f); a1.setError(a1fErr);
      alpha.setVal(alphaf); alpha.setError(alphafErr);
      n.setVal(nf); n.setError(nfErr);
      nbkg.setVal(nbkgf); nbkg.setError(nbkgfErr);
      
      RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
      RooMsgService::instance().setSilentMode(kTRUE);
      RooMsgService::instance().setStreamStatus(1,false);
      RooMsgService::instance().getStream(1).removeTopic(Integration);  
      RooMsgService::instance().getStream(1).removeTopic(Minimization);  
      RooMsgService::instance().getStream(1).removeTopic(Fitting);  
      RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
      RooMsgService::instance().getStream(1).removeTopic(Optimization);
      RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
      RooMsgService::instance().getStream(1).removeTopic(Eval);
      RooMsgService::instance().Print();
      cout << i << endl;

      RooDataSet* toyData = model.generate(RooArgSet(m), Extended(kTRUE));
      RooFitResult* result = model.fitTo(*toyData, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), NumCPU(4));

      massPull.setVal((mean.getVal() - meanf) / mean.getError());
      mpDset.add(RooArgSet(massPull));

      nsigPull.setVal((nsig.getVal() - nsigf) / nsig.getError());
      nsigDset.add(RooArgSet(nsigPull));

      nbkgPull.setVal((nbkg.getVal() - nbkgf) / nbkg.getError());
      nbkgDset.add(RooArgSet(nbkgPull));
    }

    //Fit Pull masa
    RooRealVar masPullMean("masPullMean", "masPullMean", 0, -.5, .5);
    RooRealVar massPullSigma("massPullSigma", "massPullSigma", 1, 0.5, 1.5);
    RooGaussian massPullGauss("massPullGauss", "massPullGauss", massPull, masPullMean, massPullSigma);
    massPullGauss.fitTo(mpDset, Extended(kTRUE), Save(kTRUE));

    TCanvas* c2 = new TCanvas("c2", "c2", 900, 700);
    TPad *pad1 = new TPad("p1", "", 0.01,0.01,0.99, 0.99); pad1->Draw();

    pad1->cd();
    
    RooPlot* mframe2 = massPull.frame(-6, 6, 30);
    mpDset.plotOn(mframe2);
    massPullGauss.plotOn(mframe2, LineColor(kBlue));
    
    mframe2->SetTitle("Mass pull distribution with Toy Montecarlo");
    mframe2->Draw();

    auto mtext = new TLatex();
    mtext->SetTextSize(0.04);
    mtext->SetTextFont(42);
    mtext->DrawLatex(1, 16,Form("#mu = %1.4f #pm %1.4f", masPullMean.getVal(), masPullMean.getError()));
    mtext->DrawLatex(1, 14,Form("#sigma = %1.4f #pm %1.4f", massPullSigma.getVal(), massPullSigma.getError()));

    c2->Draw();
    c2->SaveAs("monteCarlo.png");
}