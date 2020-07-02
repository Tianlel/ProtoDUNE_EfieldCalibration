#define frame_cxx
#include "frame.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;


/********** constant variables ***********/

/* output file */
string output_PATH = "/dune/app/users/tianlel/protoDUNE/E_field/ProtoDUNE_EfieldCalibration/output_ROOTtree/frame_distortion/";
string output_file_name = "test_plots.root";

/* optional variables */
Long64_t select_nentries = 0; // 0 -> use all nentries

/* frame parameters */
doublt frame_Zpositions[7] = {0,120,230,345,465,575,695};

/********** constant variables end ***********/

void frame::Loop()
{
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    if (select_nentries != 0) nentries = select_nentries;   

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if (jentry % 10000 == 0)
            std::cout << jentry << "/" << nentries << std::endl;

        if(!trkhitz_wire2->size()) continue;


    
    } // end of nentries loop
}
