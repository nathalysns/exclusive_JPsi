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

tt_event->Draw("Q2>>Q2hist",Form("%f",1/(2*lumi))*(totalcuts),"goff");
tt_event->Draw("Q2MC>>Q2histHEPMC",Form("%f",1/lumi)*(Q2cut),"goff");

auto legend = new TLegend(0.5,0.6,0.9,0.9);
legend->SetHeader("ECCE","C");
legend->AddEntry((TObject*)0, Form("e+P %s",config.Data()), "");
legend->AddEntry((TObject*)0, "J/#Psi #rightarrow ee", "");
legend->AddEntry((TObject*)0, "Q^{2} > 1 GeV^{2} ","");
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

tt_event->Draw("xbj>>Xbjhist",Form("%f",1/(2*lumi))*(totalcuts),"goff");
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

tt_event->Draw("-t>>thist",Form("%f",1/(2.*lumi))*(totalcuts),"goff");
tt_event->Draw("-tMC>>thistHEPMC",Form("%f",1/lumi),"goff");

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

//===============================================================================
TCanvas *c6 = new TCanvas("c6","",0,0,500,700);
auto *p9 = new TPad("p9","",0.,0.,1.,0.3); p9->Draw();
p9->SetTopMargin(0.001);
p9->SetBottomMargin(0.3);
p9->SetGrid();
auto *p10 = new TPad("p10","",0.,0.35,1.,1.);  p10->Draw();
p10->SetBottomMargin(0.001);
p10->cd();
p10->SetGrid();
TH1F* ptproton = new TH1F("ptproton","",bins,0,2);
TH1F* ptprotoncutt = new TH1F("ptprotoncutt","",bins,0,2);

tt_event->Draw("p_pT_Hepmc>>ptproton", Form("%f",1/lumi),"goff");
tt_event->Draw("p_pT>>ptprotoncutt",  Form("%f",1/(2*lumi)), "goff");

SetStyleHistoTH1ForGraphs(ptproton,"", "#eta ", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(ptprotoncutt,"", "#eta ", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
ptproton->SetFillColor(kCyan+2);
ptproton->Draw("hist");
ptprotoncutt->Draw("ep same");
legend->Draw();
TH1F* ptprotoncutoriginal = (TH1F*) ptprotoncutt->Clone();

p9->cd();
ptprotoncutoriginal->Divide(ptproton);
SetStyleHistoTH1ForGraphs(ptprotoncutoriginal,"", "Proton P_{t} [GeV]", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
ptprotoncutoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
ptprotoncutoriginal->GetYaxis()->SetRangeUser(0,1.2);
ptprotoncutoriginal->Draw("p");
c6->SaveAs(Form("%sptproton.pdf",outdir.Data()));

//===============================================================================
TCanvas *c7 = new TCanvas("c7","",0,0,500,700);
auto *p11 = new TPad("p11","",0.,0.,1.,0.3); p11->Draw();
p11->SetTopMargin(0.001);
p11->SetBottomMargin(0.3);
p11->SetGrid();

auto *p12 = new TPad("p12","",0.,0.35,1.,1.);  p12->Draw();
p12->SetBottomMargin(0.001);
p12->cd();
p12->SetGrid();
gPad->SetLogy();
TH1F* xlproton = new TH1F("xlproton","",2*bins,0.8,1.01);
TH1F* xlprotoncutt = new TH1F("xlprotoncutt","",2*bins,0.8,1.01);

tt_event->Draw(Form("p_E_Hepmc/%f>>xlproton",energyT), Form("%f",1/lumi),"goff");
tt_event->Draw(Form("p_E/%f>>xlprotoncutt",energyT),  Form("%f",1/(2*lumi))*(totalcuts), "goff");

SetStyleHistoTH1ForGraphs(xlproton,"", "#eta ", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(xlprotoncutt,"", "#eta ", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
xlproton->SetFillColor(kCyan+2);
xlproton->Draw("hist");
xlprotoncutt->Draw("ep same");
legend->Draw();
TH1F* xlprotoncutoriginal = (TH1F*) xlprotoncutt->Clone();

auto legend2 = new TLegend(0.5,0.6,0.9,0.9);
legend2->SetHeader("ECCE","C");
legend2->AddEntry((TObject*)0, Form("e+P %s",config.Data()), "");
legend2->AddEntry((TObject*)0, "J/#Psi #rightarrow ee", "");
legend2->AddEntry((TObject*)0, "Q^{2} > 1 GeV^{2} ","");
legend2->AddEntry(xlproton, "Generated", "f");
legend2->AddEntry(xlprotoncutt, "Reconstructed", "ep");

p11->cd();
xlprotoncutoriginal->Divide(xlproton);
SetStyleHistoTH1ForGraphs(xlprotoncutoriginal,"", "Proton x_{L}=E'_{p}/E_{p}", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
xlprotoncutoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
xlprotoncutoriginal->GetYaxis()->SetRangeUser(0,1.2);
xlprotoncutoriginal->Draw("p");
c7->SaveAs(Form("%sxlproton.pdf",outdir.Data()));

//===============================================================================
TCanvas *c8 = new TCanvas("c8","",0,0,500,700);
auto *p13 = new TPad("p13","",0.,0.,1.,0.3); p13->Draw();
p13->SetTopMargin(0.001);
p13->SetBottomMargin(0.3);
p13->SetGrid();

auto *p14 = new TPad("p14","",0.,0.35,1.,1.);  p14->Draw();
p14->SetBottomMargin(0.001);
p14->cd();
p14->SetGrid();
gPad->SetLogy();
TH1F* xvproton = new TH1F("xvproton","",2*bins,0.,0.05);
TH1F* xvprotoncutt = new TH1F("xvprotoncutt","",2*bins,0.,0.05);

tt_event->Draw("xv>>xvproton", Form("%f",1/lumi),"goff");
tt_event->Draw("xvMC>>xvprotoncutt",  Form("%f",1/(2*lumi))*(totalcuts), "goff");

SetStyleHistoTH1ForGraphs(xvproton,"", "x_{v} ", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(xvprotoncutt,"", "x_{v}", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
xvproton->SetFillColor(kCyan+2);
xvproton->Draw("hist");
xvprotoncutt->Draw("ep same");
legend->Draw();
TH1F* xvprotoncutoriginal = (TH1F*) xvprotoncutt->Clone();

legend->Draw();

p13->cd();
xvprotoncutoriginal->Divide(xvproton);
SetStyleHistoTH1ForGraphs(xvprotoncutoriginal,"", "x_{v}", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
xvprotoncutoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
xvprotoncutoriginal->GetXaxis()->SetNdivisions(506, 4, 0, kTRUE);
xvprotoncutoriginal->GetYaxis()->SetRangeUser(0,1.2);
xvprotoncutoriginal->Draw("p");
c8->SaveAs(Form("%sxv.pdf",outdir.Data()));

//===============================================================================
TCanvas *c9 = new TCanvas("c9","",0,0,500,700);
auto *p15 = new TPad("p15","",0.,0.,1.,0.3); p15->Draw();
p15->SetTopMargin(0.001);
p15->SetBottomMargin(0.3);
p15->SetGrid();

auto *p16 = new TPad("p16","",0.,0.35,1.,1.);  p16->Draw();
p16->SetBottomMargin(0.001);
p16->cd();
p16->SetGrid();
//gPad->SetLogy();
TH1F* phihist = new TH1F("phihist","",bins,0.,6.3);
TH1F* phicuthist = new TH1F("phicuthist","",bins,0.,6.3);

tt_event->Draw("phiJpsiMC>>phihist", Form("%f",1/lumi),"goff");
tt_event->Draw("phiJpsi>>phicuthist",  Form("%f",1/(2*lumi))*(totalcuts), "goff");

SetStyleHistoTH1ForGraphs(phicuthist,"", "#phi _{reco}", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(phihist,"", "#phi _{gen}", "counts/fb^{-1}", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);

phihist->SetFillColor(kCyan+2);
phihist->GetYaxis()->SetRangeUser(0,1.5*phihist->GetMaximum());
phihist->Draw("hist");
phicuthist->Draw("ep same");

legend->Draw();

TH1F* phicutoriginal = (TH1F*) phicuthist->Clone();

legend->Draw();

p15->cd();
phicutoriginal->Divide(phihist);
SetStyleHistoTH1ForGraphs(phicutoriginal,"", "#phi ", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
phicutoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
phicutoriginal->GetXaxis()->SetNdivisions(506, 4, 0, kTRUE);
phicutoriginal->GetYaxis()->SetRangeUser(0,1.2);
phicutoriginal->Draw("p");
c8->SaveAs(Form("%sphi.pdf",outdir.Data()));

}