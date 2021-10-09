

double xini = -83.224;
double ytight = 2;
double xtight = 5;
double ycut = ytight;
double xcut = xtight;
double ycutouter = 5;
double xcutouter = 12.5;

const int bins = 100;

TCut RPindex = "RPind == 1";
TCut RPpid = "RPpid == 2212";
TCut RPhits = "RPhits < 3";
TCut nTrackcut = "nTracks == 3"; 


TCut RPcut1    =  Form("((RPx[0] > %f && RPx[0] < %f) && (RPy[0] > %f  && RPy[0] < %f))", xini-xcut, xini+xcut, -ycut, ycut);
TCut RPcut2    =  Form("((RPx[0] > %f && RPx[0] < %f) && (RPy[0] > %f  && RPy[0] < %f))", xini-xcutouter, xini+xcutouter, -ycutouter, ycutouter);


TCut allcuts = RPpid + RPhits + !RPcut1 + RPcut2;
TString outdir = "../plots/";

TCut elec_eta_ccut = "e_eta_Hepmc < 0";
TCut elec_etarecp_ccut = "e_eta < 0";

TCut posi_jpsi_etareco_ccut = "ep_eta > -4 && ep_eta <4";
TCut elec_jpsi_etareco_ccut = "em_eta > -4 && em_eta <4";

TCut Q2cut = "Q2>1";

TCut totalcuts = allcuts  + elec_eta_ccut  +  elec_etarecp_ccut + posi_jpsi_etareco_ccut + elec_jpsi_etareco_ccut + Q2cut;