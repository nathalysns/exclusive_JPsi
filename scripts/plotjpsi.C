#include "plot.h"

void plotjpsi(
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
 
if(!tt_event){ std::cout << "tree not found... returning!"<< std::endl; return;}
TCut cut = "e_E > 4. && Q2 > 1. && e_eta < 0";
double emin = 0; double emax = 300;
//==================== HEpmc
Int_t bins1 =100;
TH2F *kinepHepmc    = new TH2F("kinepHepmc", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinemHepmc    = new TH2F("kinemHepmc", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinjpsiHepmc  = new TH2F("kinjpsiHepmc", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinelHepmc    = new TH2F("kinelHepmc", "", bins1, 0., 20, bins1, -20, 20);
TH2F *kinpHepmc     = new TH2F("kinpHepmc", "", bins1, emin, emax, bins1, -10, 20);

tt_event->Draw("ep_eta_Hepmc:ep_p_Hepmc>>kinepHepmc",cut,"goff");
tt_event->Draw("em_eta_Hepmc:em_p_Hepmc>>kinemHepmc",cut,"goff");
tt_event->Draw("Jpsi_eta_Hepmc:Jpsi_p_Hepmc>>kinjpsiHepmc",cut,"goff");
tt_event->Draw("e_eta_Hepmc:e_p_Hepmc>>kinelHepmc",cut,"goff");
tt_event->Draw("p_eta_Hepmc:p_p_Hepmc>>kinpHepmc",cut,"goff");

//========================================================================//
TCanvas *c1 = new TCanvas("c1","",0,0,1100,1000);
DrawGammaCanvasSettings( c1, 0.095, 0.01, 0.01, 0.105);
SetStyleHistoTH1ForGraphs(kinepHepmc, "e^{-} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinemHepmc, "e^{+} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinjpsiHepmc, "J/#psi", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinelHepmc, "Electron", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinpHepmc, "Proton", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
c1->Divide(3,2);
c1->cd(1);
kinepHepmc->Draw("colz");
c1->cd(2);
kinemHepmc->Draw("colz");
c1->cd(3);
kinjpsiHepmc->Draw("colz");
c1->cd(4);
kinelHepmc->Draw("colz");
c1->cd(5);
kinpHepmc->Draw("colz");

c1->SaveAs("plots/Hepmc_eta_p.pdf");

//==================== MC
TH2F *kinepMC  = new TH2F("kinepMC", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinemMC  = new TH2F("kinemMC", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinjpsiMC  = new TH2F("kinjpsiMC", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinelMC  = new TH2F("kinelMC", "", bins1, 0., 20, bins1, -20, 20);
TH2F *kinpMC  = new TH2F("kinpMC", "", bins1, emin, emax, bins1, -10, 20);

tt_event->Draw("ep_eta_MC:ep_p_MC>>kinepMC",cut,"goff");
tt_event->Draw("em_eta_MC:em_p_MC>>kinemMC",cut,"goff");
tt_event->Draw("Jpsi_eta_MC:Jpsi_p_MC>>kinjpsiMC",cut,"goff");
tt_event->Draw("e_eta_MC:e_p_MC>>kinelMC",cut,"goff");
tt_event->Draw("p_eta_MC:p_p_MC>>kinpMC",cut,"goff");

//========================================================================//
TCanvas *c2 = new TCanvas("c2","",0,0,1100,1000);
DrawGammaCanvasSettings( c2, 0.095, 0.01, 0.01, 0.105);
SetStyleHistoTH1ForGraphs(kinepMC, "e^{-} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinemMC, "e^{+} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinjpsiMC, "J/#psi", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinelMC, "Electron", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinpMC, "Proton", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
c2->Divide(3,2);
c2->cd(1);
kinepMC->Draw("colz");
c2->cd(2);
kinemMC->Draw("colz");
c2->cd(3);
kinjpsiMC->Draw("colz");
c2->cd(4);
kinelMC->Draw("colz");
c2->cd(5);
kinpMC->Draw("colz");

c2->SaveAs("plots/MC_eta_p.pdf");

//==================== Track

TH2F *kinep  = new TH2F("kinep", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinem  = new TH2F("kinem", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinjpsi  = new TH2F("kinjpsi", "", bins1, 0., 20, bins1, -10, 10);
TH2F *kinel  = new TH2F("kinel", "", bins1, 0., 20, bins1, -20, 20);
TH2F *kinp  = new TH2F("kinp", "", bins1, emin, emax, bins1, -10, 20);


tt_event->Draw("ep_eta:ep_p>>kinep",cut,"goff");
tt_event->Draw("em_eta:em_p>>kinem",cut,"goff");
tt_event->Draw("Jpsi_eta:Jpsi_p>>kinjpsi",cut,"goff");
tt_event->Draw("e_eta:e_p>>kinel",cut,"goff");
tt_event->Draw("p_eta:p_p>>kinp",cut,"goff");

//========================================================================//
TCanvas *c3 = new TCanvas("c3","",0,0,1100,1000);
DrawGammaCanvasSettings( c3, 0.095, 0.01, 0.01, 0.105);
SetStyleHistoTH1ForGraphs(kinep, "e^{-} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinem, "e^{+} from (J/#psi)", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinjpsi, "J/#psi", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinel, "Electron", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(kinp, "Proton", "p [GeV]", "#eta", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
c3->Divide(3,2);
c3->cd(1);
kinep->Draw("colz");
c3->cd(2);
kinem->Draw("colz");
c3->cd(3);
kinjpsi->Draw("colz");
c3->cd(4);
kinel->Draw("colz");
c3->cd(5);
kinp->Draw("colz");

c3->SaveAs("plots/Track_eta_p.pdf");


//==================== PiMass
TH1F *massJPSIMC  = new TH1F("massJPSIMC", "", bins1, 1.9, 5.1);
TH1F *massJPSI    = new TH1F("massJPSI", "", bins1, 1.9, 5.1);

tt_event->Draw("Jpsi_mass_MC>>massJPSIMC",cut,"goff");
tt_event->Draw("Jpsi_M>>massJPSI",cut,"goff");
//========================================================================//
TCanvas *c4 = new TCanvas("c4","",0,0,1000,500);
DrawGammaCanvasSettings( c4, 0.095, 0.01, 0.01, 0.105);
c4->Divide(2,1);
c4->cd(1);
SetStyleHistoTH1ForGraphs(massJPSIMC, "M_(J/#psi)", "M_(J/#psi) [GeV]", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(massJPSI, "M_(J/#psi)", "M_(J/#psi) [GeV]", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
massJPSIMC->SetLineColor(1);
massJPSI->SetLineColor(2);
massJPSIMC->Draw();
massJPSI->Draw("same");
c4->cd(2);
massJPSI->Draw();
c4->SaveAs("plots/JpsiMass.pdf");

//==================== Q2-xbj
TH1F *Q2MC      = new TH1F("Q2MC", "", bins1, 0, 50);
TH1F *xbjMC     = new TH1F("xbjMC", "", bins1, 0, 0.5);
TH1F *Q2        = new TH1F("Q2", "", bins1, 0, 50);
TH1F *xbj       = new TH1F("xbj", "", bins1, 0, 0.5);

tt_event->Draw("Q2MC>>Q2MC",cut,"goff");
tt_event->Draw("xbjMC>>xbjMC",cut,"goff");
tt_event->Draw("Q2>>Q2",cut,"goff");
tt_event->Draw("xbj>>xbj",cut,"goff");
//========================================================================//
TCanvas *c5 = new TCanvas("c5","",0,0,1000,500);
DrawGammaCanvasSettings( c5, 0.095, 0.01, 0.01, 0.105);

c5->Divide(2,1);
c5->cd(1);
c5->cd(1)->SetLogy();
SetStyleHistoTH1ForGraphs(xbjMC, "", "x_{bj}", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(xbj, "", "x_{bj}", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(Q2MC, "", "Q^{2} [GeV^{2}]", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(Q2, "", "Q^{2} [GeV^{2}]", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);

xbjMC->SetLineColor(1);
xbj->SetLineColor(2);
Q2MC->SetLineColor(1);
Q2->SetLineColor(2);

Q2MC->Draw();
Q2->Draw("same");
c5->cd(2);
c5->cd(2)->SetLogy();
xbjMC->Draw();
xbj->Draw("same");
c5->SaveAs("plots/kin.pdf");

//==================== eta
bins1 = 80;
TH1F *ep_eta         = new TH1F("ep_eta", "", bins1, -8, 8);
TH1F *em_eta          = new TH1F("em_eta", "", bins1, -8, 8);
TH1F *e_eta           = new TH1F("e_eta", "", bins1,  -8, 8);
TH1F *jpsi_eta        = new TH1F("jpsi_eta", "", bins1,  -8, 8);
TH1F *ep_etaMC        = new TH1F("ep_etaMC", "", bins1, -8, 8);
TH1F *em_etaMC        = new TH1F("em_etaMC", "", bins1, -8, 8);
TH1F *e_etaMC         = new TH1F("e_etaMC", "", bins1,  -8, 8);
TH1F *jpsi_etaMC      = new TH1F("jpsi_etaMC", "", bins1,  -8, 8);

tt_event->Draw("ep_eta>>ep_eta",cut,"goff");
tt_event->Draw("em_eta>>em_eta",cut,"goff");
tt_event->Draw("e_eta>>e_eta",cut,"goff");
tt_event->Draw("Jpsi_eta>>jpsi_eta",cut,"goff");

tt_event->Draw("ep_eta_MC>>ep_etaMC",cut,"goff");
tt_event->Draw("em_eta_MC>>em_etaMC",cut,"goff");
tt_event->Draw("e_eta_MC>>e_etaMC",cut,"goff");
tt_event->Draw("Jpsi_eta_MC>>jpsi_etaMC",cut,"goff");
//========================================================================//
TCanvas *c6 = new TCanvas("c6","",0,0,1000,1000);
DrawGammaCanvasSettings( c6, 0.095, 0.01, 0.01, 0.105);

SetStyleHistoTH1ForGraphs(ep_eta, "", "#eta(e^{+} from (J/#psi))", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(em_eta, "", "#eta(e^{-} from (J/#psi))", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(e_eta, "", "#eta(e^{-})", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(jpsi_eta, "", "#eta(J/#psi)", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(ep_etaMC, "", "#eta(e^{+} from (J/#psi))", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(em_etaMC, "", "#eta(e^{-} from (J/#psi))", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(e_etaMC, "", "#eta(e^{-})", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);
SetStyleHistoTH1ForGraphs(jpsi_etaMC, "", "#eta(J/#psi)", "", 0.85*textSizeSinglePad,textSizeSinglePad, 0.85*textSizeSinglePad,textSizeSinglePad, 0.9,0.81);

jpsi_etaMC->SetLineColor(1);
jpsi_eta->SetLineColor(2);
e_etaMC->SetLineColor(1);
e_eta->SetLineColor(2);
em_etaMC->SetLineColor(1);
em_eta->SetLineColor(2);
ep_etaMC->SetLineColor(1);
ep_eta->SetLineColor(2);

c6->Divide(2,2);
c6->cd(1);
ep_etaMC->Draw();
ep_eta->Draw("same");
c6->cd(2);
em_etaMC->Draw();
em_eta->Draw("same");
c6->cd(3);
jpsi_etaMC->Draw();
jpsi_eta->Draw("same");
c6->cd(4);
e_etaMC->Draw();
e_eta->Draw("same");

c6->SaveAs("plots/eta.pdf");


}