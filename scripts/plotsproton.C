#include "plot.h"
#include "cuts_list.h"


void plotsproton(

	TString inFile            	= "./../rootfiles/jpsi_final.root"

){


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


TCut allcuts = RPpid + RPhits + !RPcut1 + RPcut2;// + RPcutxmax + RPcuty + RPcutymax;


TCut tcut1 ="RPpid[0]==2212";
TCut tcut2 ="RPpid[1]==2212";
TH3F* RPhist1 = new TH3F("RPhist1","",100,2500,3000,100,-100,-50,100,-10,10);
TH3F* RPhist2 = new TH3F("RPhist2","",100,2500,3000,100,-100,-50,100,-10,10);
tt_event->Draw("RPy[0]:RPx[0]:RPz[0]>>RPhist1", allcuts, "goff");
tt_event->Draw("RPy[1]:RPx[1]:RPz[1]>>RPhist2", allcuts, "goff");
SetStyleHistoTH3ForGraphs(RPhist1,"Z [cm]", "X [cm]", "Y [cm]", 0.03, 0.05, 0.03, 0.05, 0.03, 0.05);
TCanvas *c1 = new TCanvas("c1","",0,0,800,800);
RPhist1->SetFillColor(9);
RPhist1->Draw("");
RPhist2->SetFillColor(2);
RPhist2->Draw("SAME");

TPolyLine3D *l = new TPolyLine3D(4);
l->SetPoint(0,2600,xini-xcutouter,-ycutouter);
l->SetPoint(1,2600,xini-xcutouter,ycutouter);
l->SetPoint(2,2600,xini+xcutouter,ycutouter);
l->SetPoint(3,2600,xini+xcutouter,-ycutouter);
l->SetPoint(4,2600,xini-xcutouter,-ycutouter);
l->SetLineWidth(2);
l->Draw("same");
TPolyLine3D *l2 = new TPolyLine3D(4);
l2->SetPoint(0,2600,xini-xcut,-ycut);
l2->SetPoint(1,2600,xini-xcut,ycut);
l2->SetPoint(2,2600,xini+xcut,ycut);
l2->SetPoint(3,2600,xini+xcut,-ycut);
l2->SetPoint(4,2600,xini-xcut,-ycut);
l2->SetLineColor(2);
l2->SetLineWidth(2);
l2->Draw("same");
double xini2 = -92.2032;
double zrp2 = 2800;
TPolyLine3D *l3 = new TPolyLine3D(4);
l3->SetPoint(0,zrp2,xini2-xcutouter,-ycutouter);
l3->SetPoint(1,zrp2,xini2-xcutouter,ycutouter);
l3->SetPoint(2,zrp2,xini2+xcutouter,ycutouter);
l3->SetPoint(3,zrp2,xini2+xcutouter,-ycutouter);
l3->SetPoint(4,zrp2,xini2-xcutouter,-ycutouter);
l3->SetLineWidth(2);
l3->Draw("same");
TPolyLine3D *l4 = new TPolyLine3D(4);
l4->SetPoint(0,zrp2,xini2-xcut,-ycut);
l4->SetPoint(1,zrp2,xini2-xcut,ycut);
l4->SetPoint(2,zrp2,xini2+xcut,ycut);
l4->SetPoint(3,zrp2,xini2+xcut,-ycut);
l4->SetPoint(4,zrp2,xini2-xcut,-ycut);
l4->SetLineColor(2);
l4->SetLineWidth(2);
l4->Draw("same");
c1->SaveAs(Form("%s/RP_nocuts.pdf",outdir.Data()));

TCanvas *c2 = new TCanvas("c2","",0,0,800,800);

ECCE
Sartre e+Pb 18x108.4
 ee→ ψJ/ 2 < 10 GeV21 < Q
HEPMC level
G4 input
Reconstructed

TH1F* etaproton = new TH1F("etaproton","",bins,4,5);
TH1F* etaprotoncut = new TH1F("etaprotoncut","",bins,4,5);

tt_event->Draw("p_eta_MC>>etaproton", RPpid + RPhits,"goff");
tt_event->Draw("p_eta_MC>>etaprotoncut",allcuts,"goff");

etaproton->SetLineColor(kBlack);
etaproton->Draw();
etaprotoncut->Draw("same");




}