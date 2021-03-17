/*
This program is written to calculate MOLLER envelopes between Col4 and hybrid wall 
Chandan Ghosh March 2021

*/
#include "remolltypes.hh"
void set_plot_style();
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid, std::vector<int> &col4trid, std::vector<int> &PbWalltrid);
void RotateXY(double &x,double &y);
double getAngle(double x, double y);
const double pi = acos(-1);
const double septant = (2*pi/7.);
const double septantStart = 3*septant;
const double septantStop = septantStart+septant;
int envelope_region2(TString source, TString out)
{
TChain T("T");
T.Add(Form("%s",source.Data()));
/*const Int_t nfile=10;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out/envelope_ee_newgeo/remollout_ee%d.root",ij);
//string1<<Form("remollout_beam%d.root",ij);
add_file[ij-1]=string1.str();
//cout<<string1.str()<<endl;
T.Add(add_file[ij-1]);
}*/
//T.Add("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/ExistingColl/remollout_beam*.root");
Int_t nEvents = T.GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
Double_t weight = 1.;///nEvents;
//weight = 1.;
//TFile f("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield_ep/DSConfig0Shield_ep_disk_28.root","RECREATE");
TFile f(Form("%s",out.Data()),"RECREATE");
set_plot_style();

std::map<TString,TH2D*> h_d_xy;//d==disk with out any cut radial redial distributions
std::map<TString,TH1D*> h_d_r;


TString part,part1,part2,part3,part4;
const Int_t nEnergy=1;
//0:Total; 1:hit.p<1; 2: hit.p>=1 && hit.p<10;  3: hit.p>=10 && hit.p<100; 4: hit.p>=100;
const Int_t nParticle=1;
const Int_t nDet=2;

string sParticle[nParticle] = {"electron"};//,"positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{65,"Col4Ent"},{66,"PbWall"}};
string sDet[nDet] = {"Col4Ent","PbWall"};


for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  part=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
  if(k==0){
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("xy dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),230,-250.,-20.,75,-5,70);
  h_d_r[part] = new TH1D(part+"_d_r",Form("R dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),50,0,250);
 }
  else if(k==1){
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("xy dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),400,-390.,10.,155,-5,150);
  h_d_r[part] = new TH1D(part+"_d_r",Form("R dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),80,0,400);
 }
 }
}


Double_t fRate=0;
remollEvent_t *fEvent=0;
std::vector<remollGenericDetectorHit_t> *fHit=0;
std::vector<remollEventParticle_t> *fPart=0;

T.SetBranchAddress("hit", &fHit);
T.SetBranchAddress("rate", &fRate);

for(size_t j=0;j<nEvents;j++){
//for(size_t j=0;j<1000;j++){
  T.GetEntry(j);

  std::vector<int> col2trid;
  std::vector<int>::iterator col2it;
  std::vector<int> col4trid;
  std::vector<int>::iterator col4it;
  std::vector<int> PbWalltrid;
  std::vector<int>::iterator PbWallit;
  if(fRate==0) fRate=1;
  //isValid(fHit,&test1,&test2,&test3,&test4);
  isValid1(fHit,col2trid,col4trid,PbWalltrid);
  /*for(it=det28trid.begin();it!=det28trid.end();++it)
        std::cout<<"Inside main event number "<<j<<" track id "<<*it<<std::endl;*/
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 

    Int_t det = hit.det;
    Int_t pid = abs(hit.pid);
    Int_t trid = hit.trid;
    if(hit.e<=10 && hit.vz>-3875)continue;
    //if((det==28 && pid ==11 && hit.r>640)||(det==59 && pid==11&&hit.r<=101 &&hit.x<=0)||(det==65 && pid==11 && hit.r<200 && hit.x<=0))std::cout<<"Inside Continue loop::Event "<<j<<" hit "<<i<<" det "<<det<<" r "<<hit.r<<" det28trid size "<<det28trid.size()<<" col2 size "<<col2trid.size()<<" col4exit size "<<col4trid.size()<<std::endl;
    col2it = find(col2trid.begin(),col2trid.end(),trid);
    col4it = find(col4trid.begin(),col4trid.end(),trid);
    PbWallit = find(PbWalltrid.begin(),PbWalltrid.end(),trid);
    if(!(col2it!=col2trid.end() && PbWallit!=PbWalltrid.end() && col4it!=col4trid.end())) {/*std::cout<<"Inside Continue loop::Event "<<j<<" hit "<<i<<" det "<<det<<" det28trid size "<<det28trid.size()<<" col2 size "<<col2trid.size()<<" col4exit size "<<col4trid.size()<<std::endl;*/continue;}
    //std::cout<<"After Continue :: Event "<<j<<" hit "<<i<<" trid "<<trid<<" iterator find "<<std::endl;  
    if((det>=65 && det<=66)){
//	std::cout<<"Hi Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<" phi "<< (180./acos(-1))*atan(hit.y/hit.x)<<std::endl;
		double hitx = hit.x; 
		double hity = hit.y;
		double hite = hit.e;
		//hite=1.;
		double angle = getAngle(hitx,hity);
		//std::cout<<"Before rotation x "<<hitx<<" y  "<<hity<<" angle "<<angle<<std::endl; 
		//if (angle<=septantStart || angle>=septantStop)
		if(hity<0)hity=-hity;
		RotateXY(hitx,hity);
		//std::cout<<"After rotation x "<<hitx<<" y  "<<hity<<" angle "<<(180./acos(-1))*atan2(hity,hitx)<<std::endl; 
        	part=Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());
		h_d_xy[part]->Fill(hitx,hity,(fRate)*weight*hite);
		h_d_r[part]->Fill(hit.r,(fRate)*weight*hite);
//		std::cout<<"Event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
		//h_t_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);
    }
  }
}


for(Int_t k=0;k<nDet;k++){
  for(Int_t i=0;i<nParticle;i++){
      part=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
      h_d_xy[part]->Write("", TObject::kOverwrite);
      h_d_r[part]->Write("", TObject::kOverwrite);
  }
}


return 0;
}
void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid, std::vector<int> &col4trid, std::vector<int> &PbWalltrid)
{
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid; 
    int hite = hit.e; 
    //if((det==59 || det==65) && x<=0) std::cout<<"det "<<det<<" x "<<x<<" phi "<<phi<<std::endl;
    if(det==59 && hit.r>=35 && hit.r<=101 && abs(pid)==11 && hite>10 && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col2trid.push_back(hit.trid);
    }
    if(det==65 && hit.r>=50 && hit.r<=197 && abs(pid)==11 && hite>10 && hit.vz<=-3875) {
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col4trid.push_back(hit.trid);
    }
    if(det==66 && abs(pid)==11 && hite>10 && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	PbWalltrid.push_back(hit.trid);
    }
 }
 // std::cout<<"end of isvalid"<<std::endl;
}

void RotateXY(double &x, double &y){
double angle = getAngle(x,y);
const double s1=sin(septant);
const double c1=cos(septant);
const double s2=sin(2*septant);
const double c2=cos(2*septant);
const double s3=sin(3*septant);
const double c3=cos(3*septant);
 double tx,ty;
if (angle>=0 && angle<septant){
	tx = x*c3-y*s3;
	ty = x*s3+y*c3;
	x=tx;
	y=(ty<0) ? -ty : ty;
	}
if (angle>=septant && angle<2*septant){
	tx = x*c2-y*s2;
	ty = x*s2+y*c2;
	x=tx;
	y=(ty<0) ? -ty : ty;
	}
if (angle>=2*septant && angle<3*septant){
	tx = x*c1-y*s1;
	ty = x*s1+y*c1;
	x=tx;
	y=(ty<0) ? -ty : ty;
	}
 return;
}

double getAngle(double x, double y){
double angle=atan2(y,x);
return (angle<0) ? (2*pi+angle) : angle;
}
