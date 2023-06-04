//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 15 19:56:13 2023 by ROOT version 6.28/00
// from TTree treeS/signal
// found on file: ROOTS_Bsphi_AOD_2018_3_best1.root
//////////////////////////////////////////////////////////

#ifndef classmake_h
#define classmake_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class classmake {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        massB;
   Double_t        Bpt;
   Double_t        Brapidity;
   Double_t        Bphi;
   Double_t        Bdl;
   Double_t        BdlE;
   Double_t        BdlIP;
   Double_t        BdlIPE;
   Double_t        massphi;
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
   Double_t        Mkstar1;
   Double_t        MPK;
   Double_t        MBd;
   Double_t        MLambda;
   Double_t        costheta;
   Double_t        cospsi;
   Double_t        cosphi;

   // List of branches
   TBranch        *b_massB;   //!
   TBranch        *b_Bpt;   //!
   TBranch        *b_Brapidity;   //!
   TBranch        *b_Bphi;   //!
   TBranch        *b_Bdl;   //!
   TBranch        *b_BdlE;   //!
   TBranch        *b_BdlIP;   //!
   TBranch        *b_BdlIPE;   //!
   TBranch        *b_massphi;   //!
   TBranch        *b_massJ;   //!
   TBranch        *b_Jpsipt;   //!
   TBranch        *b_Jrapidity;   //!
   TBranch        *b_Jphi;   //!
   TBranch        *b_sigLxyJ;   //!
   TBranch        *b_cosalfaJ;   //!
   TBranch        *b_sigLxyB;   //!
   TBranch        *b_cosalfaB;   //!
   TBranch        *b_pi1charge;   //!
   TBranch        *b_pi2charge;   //!
   TBranch        *b_mu1pt;   //!
   TBranch        *b_mu2pt;   //!
   TBranch        *b_mu1phi;   //!
   TBranch        *b_mu2phi;   //!
   TBranch        *b_mu1eta;   //!
   TBranch        *b_mu2eta;   //!
   TBranch        *b_pi1pt;   //!
   TBranch        *b_pi1eta;   //!
   TBranch        *b_pi1phi;   //!
   TBranch        *b_pi2pt;   //!
   TBranch        *b_pi2eta;   //!
   TBranch        *b_pi2phi;   //!
   TBranch        *b_Bpro;   //!
   TBranch        *b_Jpro;   //!
   TBranch        *b_dxypi1;   //!
   TBranch        *b_dxyEpi1;   //!
   TBranch        *b_dxysigpi1;   //!
   TBranch        *b_DRmumu;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_runn;   //!
   TBranch        *b_evtn;   //!
   TBranch        *b_Lumiblock;   //!
   TBranch        *b_Triggerjk;   //!
   TBranch        *b_Triggerjkjk;   //!
   TBranch        *b_maxPtpi;   //!
   TBranch        *b_minPtpi;   //!
   TBranch        *b_ptpipi;   //!
   TBranch        *b_Mkstar1;   //!
   TBranch        *b_MPK;   //!
   TBranch        *b_MBd;   //!
   TBranch        *b_MLambda;   //!
   TBranch        *b_costheta;   //!
   TBranch        *b_cospsi;   //!
   TBranch        *b_cosphi;   //!

   classmake(TTree *tree=0);
   virtual ~classmake();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef classmake_cxx
classmake::classmake(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ROOTS_Bsphi_AOD_2018_3_best1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ROOTS_Bsphi_AOD_2018_3_best1.root");
      }
      f->GetObject("treeS",tree);

   }
   Init(tree);
}

classmake::~classmake()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t classmake::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t classmake::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void classmake::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("massB", &massB, &b_massB);
   fChain->SetBranchAddress("Bpt", &Bpt, &b_Bpt);
   fChain->SetBranchAddress("Brapidity", &Brapidity, &b_Brapidity);
   fChain->SetBranchAddress("Bphi", &Bphi, &b_Bphi);
   fChain->SetBranchAddress("Bdl", &Bdl, &b_Bdl);
   fChain->SetBranchAddress("BdlE", &BdlE, &b_BdlE);
   fChain->SetBranchAddress("BdlIP", &BdlIP, &b_BdlIP);
   fChain->SetBranchAddress("BdlIPE", &BdlIPE, &b_BdlIPE);
   fChain->SetBranchAddress("massphi", &massphi, &b_massphi);
   fChain->SetBranchAddress("massJ", &massJ, &b_massJ);
   fChain->SetBranchAddress("Jpsipt", &Jpsipt, &b_Jpsipt);
   fChain->SetBranchAddress("Jrapidity", &Jrapidity, &b_Jrapidity);
   fChain->SetBranchAddress("Jphi", &Jphi, &b_Jphi);
   fChain->SetBranchAddress("sigLxyJ", &sigLxyJ, &b_sigLxyJ);
   fChain->SetBranchAddress("cosalfaJ", &cosalfaJ, &b_cosalfaJ);
   fChain->SetBranchAddress("sigLxyB", &sigLxyB, &b_sigLxyB);
   fChain->SetBranchAddress("cosalfaB", &cosalfaB, &b_cosalfaB);
   fChain->SetBranchAddress("pi1charge", &pi1charge, &b_pi1charge);
   fChain->SetBranchAddress("pi2charge", &pi2charge, &b_pi2charge);
   fChain->SetBranchAddress("mu1pt", &mu1pt, &b_mu1pt);
   fChain->SetBranchAddress("mu2pt", &mu2pt, &b_mu2pt);
   fChain->SetBranchAddress("mu1phi", &mu1phi, &b_mu1phi);
   fChain->SetBranchAddress("mu2phi", &mu2phi, &b_mu2phi);
   fChain->SetBranchAddress("mu1eta", &mu1eta, &b_mu1eta);
   fChain->SetBranchAddress("mu2eta", &mu2eta, &b_mu2eta);
   fChain->SetBranchAddress("pi1pt", &pi1pt, &b_pi1pt);
   fChain->SetBranchAddress("pi1eta", &pi1eta, &b_pi1eta);
   fChain->SetBranchAddress("pi1phi", &pi1phi, &b_pi1phi);
   fChain->SetBranchAddress("pi2pt", &pi2pt, &b_pi2pt);
   fChain->SetBranchAddress("pi2eta", &pi2eta, &b_pi2eta);
   fChain->SetBranchAddress("pi2phi", &pi2phi, &b_pi2phi);
   fChain->SetBranchAddress("Bpro", &Bpro, &b_Bpro);
   fChain->SetBranchAddress("Jpro", &Jpro, &b_Jpro);
   fChain->SetBranchAddress("dxypi1", &dxypi1, &b_dxypi1);
   fChain->SetBranchAddress("dxyEpi1", &dxyEpi1, &b_dxyEpi1);
   fChain->SetBranchAddress("dxysigpi1", &dxysigpi1, &b_dxysigpi1);
   fChain->SetBranchAddress("DRmumu", &DRmumu, &b_DRmumu);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("runn", &runn, &b_runn);
   fChain->SetBranchAddress("evtn", &evtn, &b_evtn);
   fChain->SetBranchAddress("Lumiblock", &Lumiblock, &b_Lumiblock);
   fChain->SetBranchAddress("Triggerjk", &Triggerjk, &b_Triggerjk);
   fChain->SetBranchAddress("Triggerjkjk", &Triggerjkjk, &b_Triggerjkjk);
   fChain->SetBranchAddress("maxPtpi", &maxPtpi, &b_maxPtpi);
   fChain->SetBranchAddress("minPtpi", &minPtpi, &b_minPtpi);
   fChain->SetBranchAddress("ptpipi", &ptpipi, &b_ptpipi);
   fChain->SetBranchAddress("Mkstar1", &Mkstar1, &b_Mkstar1);
   fChain->SetBranchAddress("MPK", &MPK, &b_MPK);
   fChain->SetBranchAddress("MBd", &MBd, &b_MBd);
   fChain->SetBranchAddress("MLambda", &MLambda, &b_MLambda);
   fChain->SetBranchAddress("costheta", &costheta, &b_costheta);
   fChain->SetBranchAddress("cospsi", &cospsi, &b_cospsi);
   fChain->SetBranchAddress("cosphi", &cosphi, &b_cosphi);
   Notify();
}

Bool_t classmake::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void classmake::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t classmake::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef classmake_cxx
