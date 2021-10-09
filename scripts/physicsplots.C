#include "plot.h"
#include "../common/ECCEStyle.C"
#include "cuts_list.h"

void physicsplots(
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

TCanvas *c2 = new TCanvas("c2","",0,0,500,700);
auto *p2 = new TPad("p2","p3",0.,0.,1.,0.3); p2->Draw();
p2->SetTopMargin(0.001);
p2->SetBottomMargin(0.3);
p2->SetGrid();
auto *p1 = new TPad("p1","p1",0.,0.35,1.,1.);  p1->Draw();
p1->SetBottomMargin(0.001);
p1->cd();
p1->SetGrid();
gPad->SetLogy();
double protetamin = 4.1;
double protetamax = 12;
double electronmin = -5;
double electronmax = 5;
TH1F* Q2hist = new TH1F("Q2hist","",bins,0,50);
TH1F* Q2histHEPMC = new TH1F("Q2histHEPMC","",bins,0,50);
SetStyleHistoTH1ForGraphs(Q2hist,"", "Q^{2}", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(Q2histHEPMC,"", "Q^{2}", "counts", 0.7*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.7*textSizeSinglePad,1.1,1.1);
tt_event->Draw("Q2>>Q2hist",totalcuts,"goff");
tt_event->Draw("Q2MC>>Q2histHEPMC","Q2MC>1","goff");
auto legend = new TLegend(0.5,0.6,0.9,0.9);
legend->SetHeader("ECCE","C");
legend->AddEntry((TObject*)0, "e+P 18x275", "");
legend->AddEntry((TObject*)0, "J/#Psi #rightarrow ee", "");
legend->AddEntry((TObject*)0, "Q^{2} > 1 GeV^{2} ");
legend->AddEntry(Q2hist, "Generated", "f");
legend->AddEntry(Q2histHEPMC, "Reconstructed", "ep");
TH1F* Q2original = (TH1F*) Q2hist->Clone();
Q2histHEPMC->SetFillColor(kCyan+2);
Q2histHEPMC->Draw();
Q2hist->Draw("pe same");
legend->Draw();
p2->cd();
Q2original->Divide(Q2histHEPMC);
SetStyleHistoTH1ForGraphs(Q2original,"", "Q^{2} [GeV^{2}] ", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
Q2original->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
Q2original->Draw("p");
Q2original->GetYaxis()->SetRangeUser(0,1);
c2->Update();
c2->SaveAs(Form("%sQ2.pdf",outdir.Data()));

TCanvas *c3 = new TCanvas("c3","",0,0,500,700);
auto *p3 = new TPad("p3","",0.,0.,1.,0.3); p3->Draw();
p3->SetTopMargin(0.001);
p3->SetBottomMargin(0.3);
p3->SetGrid();
auto *p4 = new TPad("p1","p1",0.,0.35,1.,1.);  p4->Draw();
p4->SetBottomMargin(0.001);
p4->cd();
p4->SetGrid();
gPad->SetLogy();
//gPad->SetLogx();


TH1F* Xbjhist = new TH1F("Xbjhist","",bins,0,0.5);
TH1F* XbjhistHEPMC = new TH1F("XbjhistHEPMC","",bins,0,0.5);
SetStyleHistoTH1ForGraphs(Xbjhist,"", "x_{bj}", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(XbjhistHEPMC,"", "x_{bj}", "counts", 0.7*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.7*textSizeSinglePad,1.1,1.1);
tt_event->Draw("xbj>>Xbjhist",totalcuts,"goff");
tt_event->Draw("xbjMC>>XbjhistHEPMC","Q2MC>1","goff");
TH1F* Xbjoriginal = (TH1F*) Xbjhist->Clone();
XbjhistHEPMC->SetFillColor(kCyan+2);
XbjhistHEPMC->Draw();
Xbjhist->Draw("pe same");
legend->Draw();
p3->cd();
Xbjoriginal->Divide(XbjhistHEPMC);
SetStyleHistoTH1ForGraphs(Xbjoriginal,"", "x_{bj}", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
Xbjoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
Xbjoriginal->Draw("p");
Xbjoriginal->GetYaxis()->SetRangeUser(0,1);
c3->Update();
c3->SaveAs(Form("%sxbj.pdf",outdir.Data()));

}