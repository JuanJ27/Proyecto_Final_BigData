//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 15 20:03:56 2023 by ROOT version 6.28/00
// from TTree bstree/bstree
// found on file: reducetree_Bsphi_AOD_3_best1.root
//////////////////////////////////////////////////////////

#ifndef classreduce_h
#define classreduce_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class classreduce {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        B_mass;
   Double_t        J_mass;
   Double_t        phi_mass;
   Double_t        pion1pt;
   Double_t        pion2pt;
   Double_t        ptmin;
   Double_t        ptmax;
   Double_t        pipipt;
   Double_t        etatkl;
   Double_t        etatksl;
   Double_t        bcpt;
   Double_t        rapidityB;
   Double_t        Bpro;
   Double_t        pdl;
   Double_t        pdle;
   Double_t        p1charg;
   Double_t        p2charg;
   Double_t        dxysig1;
   Double_t        etamu1;
   Double_t        etamu2;
   Double_t        mu1pt;
   Double_t        mu2pt;
   Double_t        muptmin;
   Double_t        muptmax;
   UInt_t          trijk;
   UInt_t          trijkjk;

   // List of branches
   TBranch        *b_B_mass;   //!
   TBranch        *b_J_mass;   //!
   TBranch        *b_phi_mass;   //!
   TBranch        *b_pion1pt;   //!
   TBranch        *b_pion2pt;   //!
   TBranch        *b_ptmin;   //!
   TBranch        *b_ptmax;   //!
   TBranch        *b_pipipt;   //!
   TBranch        *b_etatkl;   //!
   TBranch        *b_etatksl;   //!
   TBranch        *b_bcpt;   //!
   TBranch        *b_rapidityB;   //!
   TBranch        *b_Bpro;   //!
   TBranch        *b_pdl;   //!
   TBranch        *b_pdle;   //!
   TBranch        *b_p1charg;   //!
   TBranch        *b_p2charg;   //!
   TBranch        *b_dxysig1;   //!
   TBranch        *b_etamu1;   //!
   TBranch        *b_etamu2;   //!
   TBranch        *b_mu1pt;   //!
   TBranch        *b_mu2pt;   //!
   TBranch        *b_muptmin;   //!
   TBranch        *b_muptmax;   //!
   TBranch        *b_trijk;   //!
   TBranch        *b_trijkjk;   //!

   classreduce(TTree *tree=0);
   virtual ~classreduce();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef classreduce_cxx
classreduce::classreduce(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("reducetree_Bsphi_AOD_3_best1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("reducetree_Bsphi_AOD_3_best1.root");
      }
      f->GetObject("bstree",tree);

   }
   Init(tree);
}

classreduce::~classreduce()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t classreduce::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t classreduce::LoadTree(Long64_t entry)
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

void classreduce::Init(TTree *tree)
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

   fChain->SetBranchAddress("B_mass", &B_mass, &b_B_mass);
   fChain->SetBranchAddress("J_mass", &J_mass, &b_J_mass);
   fChain->SetBranchAddress("phi_mass", &phi_mass, &b_phi_mass);
   fChain->SetBranchAddress("pion1pt", &pion1pt, &b_pion1pt);
   fChain->SetBranchAddress("pion2pt", &pion2pt, &b_pion2pt);
   fChain->SetBranchAddress("ptmin", &ptmin, &b_ptmin);
   fChain->SetBranchAddress("ptmax", &ptmax, &b_ptmax);
   fChain->SetBranchAddress("pipipt", &pipipt, &b_pipipt);
   fChain->SetBranchAddress("etatkl", &etatkl, &b_etatkl);
   fChain->SetBranchAddress("etatksl", &etatksl, &b_etatksl);
   fChain->SetBranchAddress("bcpt", &bcpt, &b_bcpt);
   fChain->SetBranchAddress("rapidityB", &rapidityB, &b_rapidityB);
   fChain->SetBranchAddress("Bpro", &Bpro, &b_Bpro);
   fChain->SetBranchAddress("pdl", &pdl, &b_pdl);
   fChain->SetBranchAddress("pdle", &pdle, &b_pdle);
   fChain->SetBranchAddress("p1charg", &p1charg, &b_p1charg);
   fChain->SetBranchAddress("p2charg", &p2charg, &b_p2charg);
   fChain->SetBranchAddress("dxysig1", &dxysig1, &b_dxysig1);
   fChain->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   fChain->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   fChain->SetBranchAddress("mu1pt", &mu1pt, &b_mu1pt);
   fChain->SetBranchAddress("mu2pt", &mu2pt, &b_mu2pt);
   fChain->SetBranchAddress("muptmin", &muptmin, &b_muptmin);
   fChain->SetBranchAddress("muptmax", &muptmax, &b_muptmax);
   fChain->SetBranchAddress("trijk", &trijk, &b_trijk);
   fChain->SetBranchAddress("trijkjk", &trijkjk, &b_trijkjk);
   Notify();
}

Bool_t classreduce::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void classreduce::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t classreduce::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef classreduce_cxx
