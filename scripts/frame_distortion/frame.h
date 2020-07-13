//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 14 21:49:31 2019 by ROOT version 6.18/04
// from TChain velocity/Event/
//////////////////////////////////////////////////////////

#ifndef frame_h
#define frame_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class frame {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Double_t        evttime;
   Int_t           run;
   Int_t           Nactivefembs[6];
   Int_t           subrun;
   vector<double>  *T0_values;
   vector<float>   *xprojectedlen;
   vector<float>   *trackthetaxz;
   vector<float>   *trackthetayz;
   vector<float>   *trkstartx;
   vector<float>   *trkstarty;
   vector<float>   *trkstartz;
   vector<float>   *trkendx;
   vector<float>   *trkendy;
   vector<float>   *trkendz;
   vector<float>   *trklen;
   vector<int>     *TrkID;
   vector<int>     *tot_trks;
   vector<vector<float> > *hit_peakT0;
   vector<vector<int> > *hit_tpc0;
   vector<vector<int> > *hit_wire0;
   vector<vector<int> > *hit_channel0;
   vector<vector<float> > *trkhitx0;
   vector<vector<float> > *trkhity0;
   vector<vector<float> > *trkhitz0;
   vector<vector<float> > *trkdq_int0;
   vector<vector<float> > *trkdq_amp0;
   vector<vector<float> > *hit_peakT1;
   vector<vector<int> > *hit_tpc1;
   vector<vector<int> > *hit_wire1;
   vector<vector<int> > *hit_channel1;
   vector<vector<float> > *trkhitx1;
   vector<vector<float> > *trkhity1;
   vector<vector<float> > *trkhitz1;
   vector<vector<float> > *trkdq_int1;
   vector<vector<float> > *trkdq_amp1;
   vector<vector<float> > *hit_peakT2;
   vector<vector<int> > *hit_tpc2;
   vector<vector<int> > *hit_wire2;
   vector<vector<int> > *hit_channel2;
   vector<vector<float> > *trkhitx2;
   vector<vector<float> > *trkhity2;
   vector<vector<float> > *trkhitz2;
   vector<vector<float> > *trkhitz_wire2;
   vector<vector<float> > *trkdq_int2;
   vector<vector<float> > *trkdq_amp2;
   vector<vector<float> > *trkstartcosxyz;
   vector<vector<float> > *trkendcosxyz;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_run;   //!
   TBranch        *b_Nactivefembs;   //!
   TBranch        *b_surbrun;   //!
   TBranch        *b_T0_values;   //!
   TBranch        *b_xprojectedlen;   //!
   TBranch        *b_trackthetaxz;   //!
   TBranch        *b_trackthetayz;   //!
   TBranch        *b_trkstartx;   //!
   TBranch        *b_trkstarty;   //!
   TBranch        *b_trkstartz;   //!
   TBranch        *b_trkendx;   //!
   TBranch        *b_trkendy;   //!
   TBranch        *b_trkendz;   //!
   TBranch        *b_trklen;   //!
   TBranch        *b_TrkID;   //!
   TBranch        *b_tot_trks;   //!
   TBranch        *b_hit_peakT0;   //!
   TBranch        *b_hit_tpc0;   //!
   TBranch        *b_hit_wire0;   //!
   TBranch        *b_hit_channel0;   //!
   TBranch        *b_trkhitx0;   //!
   TBranch        *b_trkhity0;   //!
   TBranch        *b_trkhitz0;   //!
   TBranch        *b_trkdq_int0;   //!
   TBranch        *b_trkdq_amp0;   //!
   TBranch        *b_hit_peakT1;   //!
   TBranch        *b_hit_tpc1;   //!
   TBranch        *b_hit_wire1;   //!
   TBranch        *b_hit_channel1;   //!
   TBranch        *b_trkhitx1;   //!
   TBranch        *b_trkhity1;   //!
   TBranch        *b_trkhitz1;   //!
   TBranch        *b_trkdq_int1;   //!
   TBranch        *b_trkdq_amp1;   //!
   TBranch        *b_hit_peakT2;   //!
   TBranch        *b_hit_tpc2;   //!
   TBranch        *b_hit_wire2;   //!
   TBranch        *b_hit_channel2;   //!
   TBranch        *b_trkhitx2;   //!
   TBranch        *b_trkhity2;   //!
   TBranch        *b_trkhitz2;   //!
   TBranch        *b_trkhitz_wire2;   //!
   TBranch        *b_trkdq_int2;   //!
   TBranch        *b_trkdq_amp2;   //!
   TBranch        *b_trkstartcosxyz;   //!
   TBranch        *b_trkendcosxyz;   //!

   frame(TTree *tree=0);
   virtual ~frame();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef frame_cxx
frame::frame(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("velocity/Event",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("velocity/Event","");
      chain->Add("/dune/data2/users/tianlel/reco/v08_32_01/combined.root/velocity/Event");
      chain->Add("/dune/data/users/tianlel/reco/v08_32_01/combined2.root/velocity/Event");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

frame::~frame()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t frame::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t frame::LoadTree(Long64_t entry)
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

void frame::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   T0_values = 0;
   xprojectedlen = 0;
   trackthetaxz = 0;
   trackthetayz = 0;
   trkstartx = 0;
   trkstarty = 0;
   trkstartz = 0;
   trkendx = 0;
   trkendy = 0;
   trkendz = 0;
   trklen = 0;
   TrkID = 0;
   tot_trks = 0;
   hit_peakT0 = 0;
   hit_tpc0 = 0;
   hit_wire0 = 0;
   hit_channel0 = 0;
   trkhitx0 = 0;
   trkhity0 = 0;
   trkhitz0 = 0;
   trkdq_int0 = 0;
   trkdq_amp0 = 0;
   hit_peakT1 = 0;
   hit_tpc1 = 0;
   hit_wire1 = 0;
   hit_channel1 = 0;
   trkhitx1 = 0;
   trkhity1 = 0;
   trkhitz1 = 0;
   trkdq_int1 = 0;
   trkdq_amp1 = 0;
   hit_peakT2 = 0;
   hit_tpc2 = 0;
   hit_wire2 = 0;
   hit_channel2 = 0;
   trkhitx2 = 0;
   trkhity2 = 0;
   trkhitz2 = 0;
   trkhitz_wire2 = 0;
   trkdq_int2 = 0;
   trkdq_amp2 = 0;
   trkstartcosxyz = 0;
   trkendcosxyz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("Nactivefembs", Nactivefembs, &b_Nactivefembs);
   fChain->SetBranchAddress("subrun", &subrun, &b_surbrun);
   fChain->SetBranchAddress("T0_values", &T0_values, &b_T0_values);
   fChain->SetBranchAddress("xprojectedlen", &xprojectedlen, &b_xprojectedlen);
   fChain->SetBranchAddress("trackthetaxz", &trackthetaxz, &b_trackthetaxz);
   fChain->SetBranchAddress("trackthetayz", &trackthetayz, &b_trackthetayz);
   fChain->SetBranchAddress("trkstartx", &trkstartx, &b_trkstartx);
   fChain->SetBranchAddress("trkstarty", &trkstarty, &b_trkstarty);
   fChain->SetBranchAddress("trkstartz", &trkstartz, &b_trkstartz);
   fChain->SetBranchAddress("trkendx", &trkendx, &b_trkendx);
   fChain->SetBranchAddress("trkendy", &trkendy, &b_trkendy);
   fChain->SetBranchAddress("trkendz", &trkendz, &b_trkendz);
   fChain->SetBranchAddress("trklen", &trklen, &b_trklen);
   fChain->SetBranchAddress("TrkID", &TrkID, &b_TrkID);
   fChain->SetBranchAddress("tot_trks", &tot_trks, &b_tot_trks);
   fChain->SetBranchAddress("hit_peakT0", &hit_peakT0, &b_hit_peakT0);
   fChain->SetBranchAddress("hit_tpc0", &hit_tpc0, &b_hit_tpc0);
   fChain->SetBranchAddress("hit_wire0", &hit_wire0, &b_hit_wire0);
   fChain->SetBranchAddress("hit_channel0", &hit_channel0, &b_hit_channel0);
   fChain->SetBranchAddress("trkhitx0", &trkhitx0, &b_trkhitx0);
   fChain->SetBranchAddress("trkhity0", &trkhity0, &b_trkhity0);
   fChain->SetBranchAddress("trkhitz0", &trkhitz0, &b_trkhitz0);
   fChain->SetBranchAddress("trkdq_int0", &trkdq_int0, &b_trkdq_int0);
   fChain->SetBranchAddress("trkdq_amp0", &trkdq_amp0, &b_trkdq_amp0);
   fChain->SetBranchAddress("hit_peakT1", &hit_peakT1, &b_hit_peakT1);
   fChain->SetBranchAddress("hit_tpc1", &hit_tpc1, &b_hit_tpc1);
   fChain->SetBranchAddress("hit_wire1", &hit_wire1, &b_hit_wire1);
   fChain->SetBranchAddress("hit_channel1", &hit_channel1, &b_hit_channel1);
   fChain->SetBranchAddress("trkhitx1", &trkhitx1, &b_trkhitx1);
   fChain->SetBranchAddress("trkhity1", &trkhity1, &b_trkhity1);
   fChain->SetBranchAddress("trkhitz1", &trkhitz1, &b_trkhitz1);
   fChain->SetBranchAddress("trkdq_int1", &trkdq_int1, &b_trkdq_int1);
   fChain->SetBranchAddress("trkdq_amp1", &trkdq_amp1, &b_trkdq_amp1);
   fChain->SetBranchAddress("hit_peakT2", &hit_peakT2, &b_hit_peakT2);
   fChain->SetBranchAddress("hit_tpc2", &hit_tpc2, &b_hit_tpc2);
   fChain->SetBranchAddress("hit_wire2", &hit_wire2, &b_hit_wire2);
   fChain->SetBranchAddress("hit_channel2", &hit_channel2, &b_hit_channel2);
   fChain->SetBranchAddress("trkhitx2", &trkhitx2, &b_trkhitx2);
   fChain->SetBranchAddress("trkhity2", &trkhity2, &b_trkhity2);
   fChain->SetBranchAddress("trkhitz2", &trkhitz2, &b_trkhitz2);
   fChain->SetBranchAddress("trkhitz_wire2", &trkhitz_wire2, &b_trkhitz_wire2);
   fChain->SetBranchAddress("trkdq_int2", &trkdq_int2, &b_trkdq_int2);
   fChain->SetBranchAddress("trkdq_amp2", &trkdq_amp2, &b_trkdq_amp2);
   fChain->SetBranchAddress("trkstartcosxyz", &trkstartcosxyz, &b_trkstartcosxyz);
   fChain->SetBranchAddress("trkendcosxyz", &trkendcosxyz, &b_trkendcosxyz);
   Notify();
}

Bool_t frame::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void frame::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t frame::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef frame_cxx
