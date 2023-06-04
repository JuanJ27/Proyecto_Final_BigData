#include "classmake.C"

#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include  <vector>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodBDT.h"
#include <TVector.h>

using namespace std;

void reduce_rootuples_Bsphi(int era=3, int vq=1)
{

TChain *ch = new TChain("treeS",""); 
ch->Add(Form("ROOTS_Bsphi_AOD_2018_%1i_best%1i.root/treeS",era,vq)); 

TTree *tree = (TTree*) ch;
classmake t(tree);
Long64_t nentries = t.fChain->GetEntries();
cout<<" Entries : "<<nentries<<endl;
cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
cout<<" Reading data ..."<<nentries<<endl;

//------------------------------

TString outfileName( Form("reducetree_Bsphi_AOD_%1i_best%1i.root",era,vq)  ); 
    
TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
 
//*******************************
//*     VARIABLES TO SAVE       *
//*******************************

Double_t B_mass, J_mass, phi_mass;
Double_t pion1pt, pion2pt, bcpt, ptmin,ptmax, pipipt,  etatkl, etatksl ;
 Double_t rapidityB, pdl, pdle, cphaB, Bpro;
Double_t p1charg, p2charg;
Double_t dxysig1;
Double_t etamu1, etamu2, mu1pt, mu2pt, muptmin, muptmax;
UInt_t   trijk, trijkjk;
 
TTree *bstree =  new TTree("bstree","bstree");
bstree->Branch("B_mass", &B_mass,  "B_mass/D");
bstree->Branch("J_mass", &J_mass,  "J_mass/D");
bstree->Branch("phi_mass", &phi_mass,  "phi_mass/D");

bstree->Branch("pion1pt", &pion1pt,  "pion1pt/D");
bstree->Branch("pion2pt", &pion2pt,  "pion2pt/D");
bstree->Branch("ptmin", &ptmin,  "ptmin/D");
bstree->Branch("ptmax", &ptmax,  "ptmax/D");
bstree->Branch("pipipt", &pipipt,  "pipipt/D");
bstree->Branch("etatkl",  &etatkl,    "etatkl/D"); 
bstree->Branch("etatksl",  &etatksl,    "etatksl/D");
bstree->Branch("bcpt", &bcpt,  "bcpt/D");
bstree->Branch("rapidityB", &rapidityB,  "rapidityB/D");
//bstree->Branch("cphaB", &cphaB,  "cphaB/D");
bstree->Branch("Bpro", &Bpro,  "Bpro/D");
bstree->Branch("pdl", &pdl,  "pdl/D");
bstree->Branch("pdle", &pdle,  "pdle/D");  
bstree->Branch("p1charg", &p1charg,  "p1charg/D"); 
bstree->Branch("p2charg", &p2charg,  "p2charg/D");
bstree->Branch("dxysig1", &dxysig1,  "dxysig1/D");
bstree->Branch("etamu1",  &etamu1,    "etamu1/D");
bstree->Branch("etamu2",  &etamu2,    "etamu2/D");
bstree->Branch("mu1pt",  &mu1pt,    "mu1pt/D");
bstree->Branch("mu2pt",  &mu2pt,    "mu2pt/D"); 
bstree->Branch("muptmin", &muptmin,  "muptmin/D");
bstree->Branch("muptmax", &muptmax,  "muptmax/D"); 
bstree->Branch("trijk",  &trijk,    "trijk/i");
bstree->Branch("trijkjk",  &trijkjk,    "trijkjk/i");

Int_t nTen = nentries/10;
Int_t k=0;
Int_t nbytes = 0, nb = 0;
for(Long64_t jentry=0; jentry<nentries;jentry++)
    {
      //esto es lo que mmuestra el porcentaje
      Long64_t ientry = t.LoadTree(jentry);
      if (ientry < 0) break;
      nb = t.fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
      if(jentry==nentries-1) cout<<endl;
      //aqui termina lo del porcentaje
      //Mass windows cuts
      if(t.massB<=5.0 || t.massB>=5.8) continue;
      if(t.massJ<=2.9 || t.massJ>=3.3) continue;
      //if(t.massphi<=1.01946-0.010 || t.massphi>=1.01946+0.010) continue;

      if(t.DRmumu>1.2) continue;
      if(t.mu1pt<4.0) continue;
      if(t.mu2pt<4.0) continue;
      if(t.Jpro<0.1)continue;   
      //if(t.Jpsipt<7.9) continue;
      
      //if(t.Brapidity<-2.2 || t.Brapidity>2.2) continue;
      if(t.Bpt<10.0) continue;
      if(t.Bpro<0.1)continue;  
      //if(t.BdlIP<0.01) continue;
      //if((t.BdlIP/t.BdlIPE)<5.0) continue;

      if(t.pi1pt<1.2) continue;
      if(t.pi2pt<1.2) continue;      
      
      if(t.sigLxyJ<3.0)continue;
      if(t.cosalfaJ<0.9)continue;
      //if(t.cosalfaB<0.999)continue;

      if( abs(t.mu1eta)>2.4)continue;
      if( abs(t.mu2eta)>2.4)continue;
      if( abs(t.pi1eta)>2.4)continue;
      if( abs(t.pi2eta)>2.4)continue;
      if( abs(t.Brapidity)>2.4) continue;

      if(t.pi1charge==t.pi2charge)continue;         

      //*******************
      //*   fill tree    *
      //*******************

      if( t.mu1pt > t.mu2pt ){
	muptmax = t.mu1pt;
	muptmin = t.mu2pt;
	etamu1  = t.mu1eta;
	etamu2  = t.mu2eta;
      }
      else{
	muptmax = t.mu2pt;
	muptmin = t.mu1pt;
	etamu1  = t.mu2eta;
	etamu2  = t.mu1eta;
      }

      if( t.pi1pt > t.pi2pt ){
	etatkl  = t.pi1eta;
	etatksl  = t.pi2eta;
      }
      else{
	etatkl  = t.pi2eta;
	etatksl  = t.pi1eta;
      }
      
      B_mass = t.massB;
      J_mass  = t.massJ;
      phi_mass = t.massphi;   
      pion1pt = t.pi1pt;
      pion2pt = t.pi2pt;
      ptmin = t.minPtpi;
      ptmax = t.maxPtpi;
      pipipt= t.ptpipi;
      bcpt    = t.Bpt;
      rapidityB = t.Brapidity;
      Bpro = t.Bpro;
      //cphaB = t.cosalfaB;
      pdl = t.BdlIP;
      pdle = t.BdlIPE;
      dxysig1 = t.dxysigpi1;
      mu1pt = t.mu1pt;
      mu2pt = t.mu2pt;
      //etamu1  = t.mu1eta;
      //etamu2  = t.mu2eta;
      trijk  = t.Triggerjk;
      trijkjk  = t.Triggerjkjk;      
      p1charg = t.pi1charge;
      p2charg = t.pi2charge;

      bstree->Fill();

    }

 outputFile->Write("",TObject::kOverwrite);

}//End analysis
