

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


TCut RPcut1    =  Form("((RPx[0] > %f && RPx[0] < %f) && (RPy[0] > %f  && RPy[0] < %f))", xini-xcut, xini+xcut, -ycut, ycut);
TCut RPcut2    =  Form("((RPx[0] > %f && RPx[0] < %f) && (RPy[0] > %f  && RPy[0] < %f))", xini-xcutouter, xini+xcutouter, -ycutouter, ycutouter);

TCut RPcutxmax =  Form("RPx > %f && RPx < %f", xini-xcutouter, xini+xcutouter);
TCut RPcuty    =  Form("RPy < %f || RPy > %f", -ycut, ycut);
TCut RPcutymax =  Form("RPy > %f && RPy < %f", -ycutouter, ycutouter);

TCut allcuts = RPpid + RPhits + !RPcut1 + RPcut2;
TString outdir = "../plots/";