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



TH1F* etaproton = new TH1F("etaproton","",bins,protetamin,protetamax);
TH1F* etaprotoncut = new TH1F("etaprotoncut","",bins,protetamin,protetamax);
SetStyleHistoTH1ForGraphs(etaproton,"", "#eta ", "counts", 0.7*textSizeSinglePad,textSizeSinglePad,0.7*textSizeSinglePad,textSizeSinglePad,1.1,1.4);
SetStyleHistoTH1ForGraphs(etaprotoncut,"", "#eta ", "counts", 0.7*textSizeSinglePad,0.7*textSizeSinglePad, 0.7*textSizeSinglePad,0.7*textSizeSinglePad,1.1,1.1);

tt_event->Draw("p_eta_MC>>etaproton", RPpid + RPhits,"goff");
tt_event->Draw("p_eta_MC>>etaprotoncut",allcuts,"goff");


auto legend = new TLegend(0.5,0.6,0.9,0.9);
legend->SetHeader("ECCE","C");
legend->AddEntry((TObject*)0, "e+P 18x275", "");
legend->AddEntry((TObject*)0, "J/#Psi #rightarrow ee", "");
legend->AddEntry((TObject*)0, "Q^{2} > 1 GeV^2 ");
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
c2->SaveAs(Form("%s/Protons_accepted.pdf",outdir.Data()));

}