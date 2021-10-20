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
SetStyleHistoTH1ForGraphs(Q2hist,"", "Q^{2}", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(Q2histHEPMC,"", "Q^{2}", "counts/fb^{-1}", 0.7*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.7*textSizeSinglePad,1.1,1.1);
tt_event->Draw("Q2>>Q2hist",Form("%f",1/lumi)*(totalcuts),"goff");
tt_event->Draw("Q2MC>>Q2histHEPMC",Form("%f",1/lumi),"goff");
auto legend = new TLegend(0.5,0.6,0.9,0.9);
legend->SetHeader("ECCE","C");
legend->AddEntry((TObject*)0, "e+P 18x275", "");
legend->AddEntry((TObject*)0, "J/#Psi #rightarrow ee", "");
legend->AddEntry((TObject*)0, "Q^{2} > 1 GeV^{2} ");
legend->AddEntry(Q2hist, "Generated", "f");
legend->AddEntry(Q2histHEPMC, "Reconstructed", "ep");
TH1F* Q2original = (TH1F*) Q2hist->Clone();
Q2histHEPMC->SetFillColor(kCyan+2);
Q2histHEPMC->Draw("hist");
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
auto *p4 = new TPad("p4","p4",0.,0.35,1.,1.);  p4->Draw();
p4->SetBottomMargin(0.001);
p4->cd();
p4->SetGrid();
gPad->SetLogy();
TH1F* Xbjhist = new TH1F("Xbjhist","",bins,0,0.1);
TH1F* XbjhistHEPMC = new TH1F("XbjhistHEPMC","",bins,0,0.1);
SetStyleHistoTH1ForGraphs(Xbjhist,"", "x_{bj}", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(XbjhistHEPMC,"", "x_{bj}", "counts/fb^{-1}", 0.7*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.7*textSizeSinglePad,1.1,1.1);
tt_event->Draw("xbj>>Xbjhist",Form("%f",1/lumi)*(totalcuts),"goff");
tt_event->Draw("xbjMC>>XbjhistHEPMC",Form("%f",1/lumi)*(Q2cut),"goff");
TH1F* Xbjoriginal = (TH1F*) Xbjhist->Clone();
XbjhistHEPMC->SetFillColor(kCyan+2);
XbjhistHEPMC->Draw("hist");
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

TCanvas *c4 = new TCanvas("c4","",0,0,500,700);
auto *p5 = new TPad("p5","",0.,0.,1.,0.3); p5->Draw();
p5->SetTopMargin(0.001);
p5->SetBottomMargin(0.3);
p5->SetGrid();
auto *p6 = new TPad("p6","",0.,0.35,1.,1.);  p6->Draw();
p6->SetBottomMargin(0.001);
p6->cd();
p6->SetGrid();
gPad->SetLogy();
TH1F* thist = new TH1F("thist","",bins80,0,2);
TH1F* thistHEPMC = new TH1F("thistHEPMC","",bins80,0,2);
SetStyleHistoTH1ForGraphs(thist,"", "-t [GeV^{2}]", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.5);
SetStyleHistoTH1ForGraphs(thistHEPMC,"", "-t [GeV^{2}]", "counts/fb^{-1}", 0.7*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.7*textSizeSinglePad,1.1,1.5);
tt_event->Draw("-t>>thist",Form("%f",1/lumi)*(totalcuts),"goff");
tt_event->Draw("-tHepmc>>thistHEPMC",Form("%f",1/lumi)*(Q2cut),"goff");
TH1F* toriginal = (TH1F*) thist->Clone();
thistHEPMC->SetFillColor(kCyan+2);
thistHEPMC->Draw("hist");
thist->Draw("pe same");
legend->Draw();
p5->cd();
toriginal->Divide(thistHEPMC);
SetStyleHistoTH1ForGraphs(toriginal,"", "-t [GeV^{2}]", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
toriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
toriginal->Draw("p");
toriginal->GetYaxis()->SetRangeUser(0,1);
c4->Update();
c4->SaveAs(Form("%sxbj.pdf",outdir.Data()));

}