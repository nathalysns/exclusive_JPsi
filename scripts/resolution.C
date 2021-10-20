#include "plot.h"
#include "../common/ECCEStyle.C"
#include "cuts_list.h"

void resolution(
	TString inFile            	= "./../rootfiles/jpsi_final.root"
){

gStyle->SetOptStat(0);
// load tree
TChain *const tt_event = new TChain("T");

    if (inFile.EndsWith(".root")) {                     // are we loading a single root tree?
        std::cout << "loading a single root file" << std::endl;
        tt_event->AddFile(inFile);
    }
    else {                                              // or are we loading a bunch?
        std::cout << "loading a list of files" << std::endl;
        std::ifstream files(inFile);
        std::string filePath;

        while (std::getline(files, filePath)) {
            tt_event->AddFile(filePath.c_str());
        }
        files.close();
    }
 
SetECCEStyle();

TCanvas *c1 = new TCanvas("c1","",0,0,1000,500);
c1->Divide(2,1);
c1->cd(1);
TString Q2resoform = "(Q2-Q2MC)";
TH1F* Q2reso = new TH1F("Q2reso","",bins,-0.5,0.5);
SetStyleHistoTH1ForGraphs(Q2reso,"","(Q^{2}_{reco} - Q^{2}_{gen})",  "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
tt_event->Draw(Form("%s>>Q2reso",Q2resoform.Data()), totalcuts,"goff");
Q2reso->Draw();
c1->cd(2);
TH2F* Q2reso2 = new TH2F("Q2reso2","", bins, -8, 1, bins, 0.5, 1.5);
SetStyleHistoTH2ForGraphs(Q2reso2, "#eta electron _{gen}", "Q^{2}_{reco}/Q^{2}_{gen}", 0.04, 0.05, 0.04, 0.05, 1.1, 1.1);
tt_event->Draw("(Q2MC/Q2):e_eta_Hepmc>>Q2reso2", totalcuts,"goff");
Q2reso2->Draw("colz");
c1->SaveAs(Form("%sQ2resolution.pdf",outdir.Data()));

TCanvas *c2 = new TCanvas("c2","",0,0,1000,500);
c2->Divide(2,1);
c2->cd(1);
TString tresoform = "(t-tMC)";
TH1F* treso = new TH1F("treso","",bins,-0.5,0.5);
SetStyleHistoTH1ForGraphs(treso,"", "(t_{reco} - t_{gen})", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
tt_event->Draw(Form("%s>>treso",tresoform.Data()), totalcuts,"goff");
treso->Draw();
c2->cd(2);
TH2F* treso2 = new TH2F("treso2","", bins,4, 10, bins, 0.5, 1.5);
SetStyleHistoTH2ForGraphs(treso2, "#eta proton _{gen}", "t_{reco}/t_{gen}", 0.04, 0.05, 0.04, 0.05, 1.1, 1.1);
tt_event->Draw("(t/tMC):p_eta_Hepmc>>treso2", totalcuts,"goff");
treso2->Draw("colz");
c2->SaveAs(Form("%stresolution.pdf",outdir.Data()));


TCanvas *c3 = new TCanvas("c3","",0,0,1000,500);
c3->Divide(2,1);
c3->cd(1);
TString xbjresoform = "(xbj-xbjMC)";
TH1F* xbjreso = new TH1F("xbjreso","",bins,-0.5,0.5);
SetStyleHistoTH1ForGraphs(xbjreso,"", "(x_{reco} - x_{gen})", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
tt_event->Draw(Form("%s>>xbjreso",xbjresoform.Data()), totalcuts,"goff");
xbjreso->Draw();
c3->cd(2);
TH2F* xbjreso2 = new TH2F("xbjreso2","", bins,-8, 1, bins, 0.5, 1.5);
SetStyleHistoTH2ForGraphs(xbjreso2, "#eta proton _{gen}", "x_{reco}/x_{gen}", 0.04, 0.05, 0.04, 0.05, 1.1, 1.1);
tt_event->Draw("(xbj/xbjMC):e_eta_Hepmc>>xbjreso2", totalcuts,"goff");
xbjreso2->Draw("colz");
c3->SaveAs(Form("%sxbjresolution.pdf",outdir.Data()));



}