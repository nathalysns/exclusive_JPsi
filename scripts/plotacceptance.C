#include "plot.h"
#include "../common/ECCEStyle.C"
#include "cuts_list.h"

void plotacceptance(
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

double protetamin = 4.1;
double protetamax = 12;

double electronmin = -5;
double electronmax = 5;


TH1F* etaproton = new TH1F("etaproton","",bins,protetamin,protetamax);
TH1F* etaprotoncut = new TH1F("etaprotoncut","",bins,protetamin,protetamax);
SetStyleHistoTH1ForGraphs(etaproton,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(etaprotoncut,"", "#eta ", "counts", 0.7*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.7*textSizeSinglePad,1.1,1.1);
tt_event->Draw("p_eta_Hepmc>>etaproton", RPpid + RPhits,"goff");
tt_event->Draw("p_eta_Hepmc>>etaprotoncut",allcuts,"goff");
auto legend = new TLegend(0.5,0.6,0.9,0.9);
legend->SetHeader("ECCE","C");
legend->AddEntry((TObject*)0, "e+P 18x275", "");
legend->AddEntry((TObject*)0, "J/#Psi #rightarrow ee", "");
legend->AddEntry((TObject*)0, "Q^{2} > 1 GeV^{2} ");
legend->AddEntry(etaproton, "Generated", "f");
legend->AddEntry(etaprotoncut, "Reconstructed", "ep");
TH1F* etaprotonoriginal = (TH1F*) etaprotoncut->Clone();
etaproton->SetFillColor(kCyan+2);
etaproton->Draw();
etaprotoncut->Draw("pe same");
legend->Draw();
p2->cd();
etaprotonoriginal->Divide(etaproton);
SetStyleHistoTH1ForGraphs(etaprotonoriginal,"", "Protons #eta ", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
etaprotonoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
etaprotonoriginal->Draw("p");
c2->Update();
c2->SaveAs(Form("%sProtons_accepted.pdf",outdir.Data()));
TH1F* etaelectron = new TH1F("etaelectron","",bins,electronmin,electronmax);
TH1F* etaelectroncut = new TH1F("etaelectroncut","",bins,electronmin,electronmax);


TCanvas *c3 = new TCanvas("c3","",0,0,500,700);
auto *p3 = new TPad("p3","",0.,0.,1.,0.3); p3->Draw();
p3->SetTopMargin(0.001);
p3->SetBottomMargin(0.3);
p3->SetGrid();
auto *p4 = new TPad("p4","",0.,0.35,1.,1.);  p4->Draw();
p4->SetBottomMargin(0.001);
p4->cd();
p4->SetGrid();
tt_event->Draw("e_eta_Hepmc>>etaelectron", elec_eta_ccut ,"goff");
tt_event->Draw("e_eta>>etaelectroncut", nTrackcut + elec_eta_ccut + elec_etarecp_ccut,"goff");
SetStyleHistoTH1ForGraphs(etaelectron,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(etaelectroncut,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
etaelectron->SetFillColor(kCyan+2);
etaelectron->Draw("");
etaelectroncut->Draw("ep  same");
TH1F* etaelectronoriginal = (TH1F*) etaelectroncut->Clone();
legend->Draw();
p3->cd();
etaelectronoriginal->Divide(etaelectron);
SetStyleHistoTH1ForGraphs(etaelectronoriginal,"", "Electrons #eta ", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
etaelectronoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
etaelectronoriginal->GetYaxis()->SetRangeUser(0,1.2);
etaelectronoriginal->Draw("p");
c3->SaveAs(Form("%sElectrons_accepted.pdf",outdir.Data()));

//===============================================================================
TCanvas *c4 = new TCanvas("c4","",0,0,500,700);
auto *p5 = new TPad("p5","",0.,0.,1.,0.3); p5->Draw();
p5->SetTopMargin(0.001);
p5->SetBottomMargin(0.3);
p5->SetGrid();
auto *p6 = new TPad("p4","",0.,0.35,1.,1.);  p6->Draw();
p6->SetBottomMargin(0.001);
p6->cd();
p6->SetGrid();
TH1F* etaelectropn = new TH1F("etaelectropn","",bins,electronmin,electronmax+3);
TH1F* etaelectronpncut = new TH1F("etaelectronpncut","",bins,electronmin,electronmax+3);
tt_event->Draw("ep_eta_Hepmc>>etaelectropn", "" ,"goff");
tt_event->Draw("ep_eta>>etaelectronpncut", nTrackcut + posi_jpsi_etareco_ccut,"goff");
SetStyleHistoTH1ForGraphs(etaelectropn,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(etaelectronpncut,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
etaelectropn->SetFillColor(kCyan+2);
etaelectropn->Draw("");
etaelectronpncut->Draw("ep  same");
TH1F* etaelectronporiginal = (TH1F*) etaelectronpncut->Clone();
legend->Draw();
p5->cd();
etaelectronporiginal->Divide(etaelectropn);
SetStyleHistoTH1ForGraphs(etaelectronporiginal,"", "e^{+} (from J/#psi) #eta ", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
etaelectronporiginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
etaelectronporiginal->GetYaxis()->SetRangeUser(0,1.2);
etaelectronporiginal->Draw("p");
c4->SaveAs(Form("%spositron_jpsi_accepted.pdf",outdir.Data()));

//===============================================================================
TCanvas *c5 = new TCanvas("c5","",0,0,500,700);
auto *p7 = new TPad("p7","",0.,0.,1.,0.3); p7->Draw();
p7->SetTopMargin(0.001);
p7->SetBottomMargin(0.3);
p7->SetGrid();
auto *p8 = new TPad("p4","",0.,0.35,1.,1.);  p8->Draw();
p8->SetBottomMargin(0.001);
p8->cd();
p8->SetGrid();
TH1F* etaelectromn = new TH1F("etaelectromn","",bins,electronmin,electronmax+4);
TH1F* etaelectromncut = new TH1F("etaelectromncut","",bins,electronmin,electronmax+4);
tt_event->Draw("em_eta_Hepmc>>etaelectromn", "" ,"goff");
tt_event->Draw("em_eta>>etaelectromncut", nTrackcut + elec_jpsi_etareco_ccut, "goff");

SetStyleHistoTH1ForGraphs(etaelectromn,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(etaelectromncut,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
etaelectromn->SetFillColor(kCyan+2);
etaelectromn->Draw("");
etaelectromncut->Draw("ep same");
legend->Draw();
TH1F* etaelectronmoriginal = (TH1F*) etaelectromncut->Clone();

p7->cd();
etaelectronmoriginal->Divide(etaelectromn);
SetStyleHistoTH1ForGraphs(etaelectronmoriginal,"", "e^{-} (from J/#psi) #eta ", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
etaelectronmoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
etaelectronmoriginal->GetYaxis()->SetRangeUser(0,1.2);
etaelectronmoriginal->Draw("p");
c5->SaveAs(Form("%selectron_jpsi_accepted.pdf",outdir.Data()));

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
tt_event->Draw("p_pT_Hepmc>>ptproton", RPpid + RPhits,"goff");
tt_event->Draw("p_pT_Hepmc>>ptprotoncutt",  allcuts, "goff");

SetStyleHistoTH1ForGraphs(ptproton,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(ptprotoncutt,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
ptproton->SetFillColor(kCyan+2);
ptproton->Draw("");
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
gPad->SetLogy();
auto *p12 = new TPad("p11","",0.,0.35,1.,1.);  p12->Draw();
p12->SetBottomMargin(0.001);
p12->cd();
p12->SetGrid();
TH1F* xlproton = new TH1F("xlproton","",bins,0,2);
TH1F* xlprotoncutt = new TH1F("xlprotoncutt","",bins,0,2);
tt_event->Draw("p_E_Hepmc/275.>>xlproton", RPpid + RPhits,"goff");
tt_event->Draw("p_E_Hepmc/275.>>xlprotoncutt",  allcuts, "goff");

SetStyleHistoTH1ForGraphs(xlproton,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(xlprotoncutt,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
xlproton->SetFillColor(kCyan+2);
xlproton->Draw("");
xlprotoncutt->Draw("ep same");
legend->Draw();
TH1F* xlprotoncutoriginal = (TH1F*) xlprotoncutt->Clone();

p11->cd();
xlprotoncutoriginal->Divide(xlproton);
SetStyleHistoTH1ForGraphs(xlprotoncutoriginal,"", "Proton x_{L}=E'_{p}/E_{p}", "ratio", 2*textSizeSinglePad,2*textSizeSinglePad, 1.7*textSizeSinglePad,2*textSizeSinglePad,1.1,0.5);
xlprotoncutoriginal->GetYaxis()->SetNdivisions(506, 4, 0, kTRUE);
xlprotoncutoriginal->GetYaxis()->SetRangeUser(0,1.2);
xlprotoncutoriginal->Draw("p");
c7->SaveAs(Form("%sxlproton.pdf",outdir.Data()));

}