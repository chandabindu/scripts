/*
This program is written to see what outer radius of collimator 2 is really needed.
Chandan Ghosh April 2020

*/
#include "remolltypes.hh"
void set_plot_style();
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid, std::vector<int> &col2trid1, std::vector<int> &col2trid2, std::vector<int> &col2trid3);
void RotateXY(double &x,double &y);
double getAngle(double x, double y);
const double pi = acos(-1);
const double septant = (2*pi/7.);
const double midangle = 360./14.;
const double septantStart = 3*septant;
const double septantStop = septantStart+septant;
int Coll2_design(TString source, TString out)
{
TChain T("T");
T.Add(Form("%s",source.Data()));
/*const Int_t nfile=10;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out/Collar_moller_NoShld/remollout_moller%d.root",ij);
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

std::map<TString,TH2D*> h_d_xy;
std::map<TString,TH2D*> h_d_rz;
std::map<TString,TH1D*> h_d_r;
std::map<TString,TH1D*> h_d_p;
std::map<TString,TH1D*> h_d_th;

std::map<TString,TH2D*> h_d1_xy;
std::map<TString,TH2D*> h_d1_rz;
std::map<TString,TH2D*> h_d1_rph;
std::map<TString,TH1D*> h_d1_r;
std::map<TString,TH1D*> h_d1_p;
std::map<TString,TH1D*> h_d_ph;


TString part,part1,part2,part3,part4;
const Int_t nEnergy=1;
//0:Total; 1:hit.p<1; 2: hit.p>=1 && hit.p<10;  3: hit.p>=10 && hit.p<100; 4: hit.p>=100;
const Int_t nParticle=1;
const Int_t nDet=2;

string sParticle[nParticle] = {"electron"};//,"positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{76,"Col2Ent"},{77,"Col2Exit"}};
string sDet[nDet] = {"Col2Ent","Col2Exit"};


for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  part=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
  h_d_p[part] = new TH1D(part+"_d_p",Form("Momemtum dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),1200,0,12000);
  h_d_rz[part] = new TH2D(part+"_d_rz",Form("Vertex rz dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),3800,-6000,32000,200,-200,200);
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("xy dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),210,-220.,-10.,75,-5,70);
  h_d_r[part] = new TH1D(part+"_d_r",Form("R dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),44,20,240);
  h_d_th[part] = new TH1D(part+"_d_th",Form("Scatt. Angle at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),100,0,2.5);

  h_d1_p[part] = new TH1D(part+"_d1_p",Form("Momemtum dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),1200,0,12000);
  h_d1_rz[part] = new TH2D(part+"_d1_rz",Form("Vertex rz dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),3800,-6000,32000,200,-200,200);
  h_d1_rph[part] = new TH2D(part+"_d1_rph",Form("r vs Ph at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),130,0,130,100,-0.5,0.5);
  h_d1_xy[part] = new TH2D(part+"_d1_xy",Form("xy dist. at %s for %s",sDet[k].c_str(),sParticle[i].c_str()),300,-150.,150.,300,-150,150);
  h_d1_r[part] = new TH1D(part+"_d1_r",Form("R dist. at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),44,20,240);
  h_d_ph[part] = new TH1D(part+"_d_ph",Form("Phi. Angle at  %s for %s",sDet[k].c_str(),sParticle[i].c_str()),100,-0.5,0.5);
}
}


Double_t fRate=0;
remollEvent_t *fEvent=0;
std::vector<remollGenericDetectorHit_t> *fHit=0;
std::vector<remollEventParticle_t> *fPart=0;

T.SetBranchAddress("hit", &fHit);
T.SetBranchAddress("rate", &fRate);

for(size_t j=0;j<nEvents;j++){
  T.GetEntry(j);

  std::vector<int> col2trid;
  std::vector<int>::iterator col2it;
  std::vector<int> col2trid1;
  std::vector<int>::iterator col2it1;

  std::vector<int> col2trid2;//col2 exit
  std::vector<int>::iterator col2it2;
  std::vector<int> col2trid3;//coll2 exit
  std::vector<int>::iterator col2it3;

  if(fRate==0) fRate=1;
  isValid1(fHit,col2trid,col2trid1,col2trid2,col2trid3);
  /*for(it=det28trid.begin();it!=det28trid.end();++it)
        std::cout<<"Inside main event number "<<j<<" track id "<<*it<<std::endl;*/
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 

    Int_t det = hit.det;
    Int_t pid = abs(hit.pid);
    Int_t trid = hit.trid;
    if(!((det==76 || det==77) && hit.vz<=-3875)) continue;
    col2it = find(col2trid.begin(),col2trid.end(),trid);
    col2it1 = find(col2trid1.begin(),col2trid1.end(),trid);
    col2it2 = find(col2trid2.begin(),col2trid2.end(),trid);
    col2it3 = find(col2trid3.begin(),col2trid3.end(),trid);

    //if(!(col2it!=col2trid.end())) {/*std::cout<<"Inside Continue loop::Event "<<j<<" hit "<<i<<" det "<<det<<" det28trid size "<<det28trid.size()<<" col2 size "<<col2trid.size()<<" col4exit size "<<col4trid.size()<<std::endl;*/continue;}
    //std::cout<<"After Continue :: Event "<<j<<" hit "<<i<<" trid "<<trid<<" iterator find "<<std::endl;  
    if((det==76 || det==77) && pid==11){
//	std::cout<<"Hi Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<" phi "<< (180./acos(-1))*atan(hit.y/hit.x)<<std::endl;
		double hitx = hit.x; 
		double hity = hit.y;
		double hite = hit.e;
		double px = hit.px; 
		double py = hit.py;
		double pz = hit.pz;
		//int trid = hit.trid;
		hite=1.;
        	part=Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());
		if(col2it1!=col2trid1.end() && col2it3!=col2trid3.end()){
		if((trid==*col2it1) && (trid==*col2it3)){
		h_d1_xy[part]->Fill(hitx,hity,(fRate)*weight);
		h_d1_rz[part]->Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy),(fRate)*weight);
		h_d1_r[part]->Fill(hit.r,(fRate)*weight);
		h_d1_p[part]->Fill(hit.p,(fRate)*weight);
		double Phi1=-1000;
		//if(px<0) px=-px;
		//if(py<0) py=-py;
		double phi = atan2(py,px)*(180./pi);
		/*if(phi<0)phi=phi*(-1);
		if(abs(phi)>=0 && abs(phi)<2*midangle) Phi1=abs(phi)-midangle;
		if(abs(phi)>=2*midangle && abs(phi)<4*midangle) Phi1=abs(phi)-3*midangle;
		if(abs(phi)>=4*midangle && abs(phi)<6*midangle) Phi1=abs(phi)-5*midangle;
		if(abs(phi)>=6*midangle && abs(phi)<=8*midangle) Phi1=abs(phi)-7*midangle;*/
		Phi1=phi;
		double Phi = atan(hit.r*sin(Phi1*acos(-1)/180.)/hit.z)*(180/pi);
		Phi = (180./pi)*(-hity*px+hitx*py)/(hit.r*pz);
		h_d1_rph[part]->Fill(hit.r,Phi,(fRate)*weight);
		h_d_ph[part]->Fill(Phi,(fRate)*weight);}
		}

		if(hity<0)hity=-hity;
		RotateXY(hitx,hity);
		if(col2it!=col2trid.end() && col2it2!=col2trid2.end()){
		if((trid==*col2it) && (trid==*col2it2)){
		h_d_xy[part]->Fill(hitx,hity,(fRate)*weight);
		h_d_rz[part]->Fill(hit.vz,sqrt(hit.vx*hit.vx+hit.vy*hit.vy),(fRate)*weight);
		h_d_r[part]->Fill(hit.r,(fRate)*weight);
		h_d_p[part]->Fill(hit.p,(fRate)*weight);
		double theta = (180./pi)*atan2(sqrt(px*px+py*py),pz);
		h_d_th[part]->Fill(theta,(fRate)*weight);}
		}	
		}
	}

}


for(Int_t k=0;k<nDet;k++){
  for(Int_t i=0;i<nParticle;i++){
      part=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
      h_d_xy[part]->Write("", TObject::kOverwrite);
      h_d_p[part]->Write("", TObject::kOverwrite);
      h_d_rz[part]->Write("", TObject::kOverwrite);
      h_d_r[part]->Write("", TObject::kOverwrite);
      h_d_th[part]->Write("", TObject::kOverwrite);

      h_d1_xy[part]->Write("", TObject::kOverwrite);
      h_d1_p[part]->Write("", TObject::kOverwrite);
      h_d1_rz[part]->Write("", TObject::kOverwrite);
      h_d1_rph[part]->Write("", TObject::kOverwrite);
      h_d1_r[part]->Write("", TObject::kOverwrite);
      h_d_ph[part]->Write("", TObject::kOverwrite);
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

void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &col2trid, std::vector<int> &col2trid1, std::vector<int> &col2trid2, std::vector<int> &col2trid3)
{
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid; 
    int hite = hit.e;
    double hitx = hit.x;
    double hity = hit.y;
    if(hity<0)hity=-hity;
    RotateXY(hitx,hity);
    double angle = getAngle(hitx,hity);
    //double angle = atan2(hity,hitx); 
    //if(det==68 && angle<2.888) std::cout<<"det "<<det<<" x "<<hitx<<" y "<<hity<<" phi "<<angle<<std::endl;
    /*if(det==68 && angle<2.888 && abs(pid)==11){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	CutPlanetrid.push_back(hit.trid);
    }
    if(det==28 && hit.r>=640 &&hit.r<=1195 && pid==11){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	det28trid.push_back(hit.trid);
    }*/
    if(det==76 && hit.r>=98 && hit.r<=104 && pid==11 && hit.vz<=-3875){//97.7 -33.7-35
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col2trid.push_back(hit.trid);
    }
    if(det==76 && angle<2.9325 && angle>2.90 && hit.r>=35 && hit.r<=104 && pid==11 && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col2trid1.push_back(hit.trid);
    }
    if(det==77 && hit.r>=98 && hit.r<=120 && pid==11 && hit.vz<=-3875){//97.7 -33.7-35
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col2trid2.push_back(hit.trid);
    }
    if(det==77 && angle<2.96 && angle>2.8 && hit.r>=35 && hit.r<=104 && pid==11  && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col2trid3.push_back(hit.trid);
    }
    /*if(det==65 && hit.r>=180 && hit.r<=197 && pid==11){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	col4trid.push_back(hit.trid);
    }*/
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
