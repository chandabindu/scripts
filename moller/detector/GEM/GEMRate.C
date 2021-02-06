/*
This program is written to calculate the theta and phi angles at det 28 of high energy electrons for engineers to design the detector support structure. The angles are calculated based on the three momentum components.
Chandan Ghosh Jan 2021

*/
#include "remolltypes.hh"
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &det28trid, std::vector<int> &det68trid);
const double pi = acos(-1);
const double rad2deg = 180./pi;
const double midangle=360./14.;
const int nRing = 6;
const int nSector = 3;
int whichRing(Double_t,Double_t, Int_t, Double_t &RMax);
int whichSector(Double_t,Double_t);
int GEMRate(TString source, TString out)
{
TChain T("T");
T.Add(Form("%s",source.Data()));
/*const Int_t nfile=3;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
//string1<<Form("/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out/SymmetricCoil3_ep/remollout_ep%d.root",ij);
string1<<Form("remollout_beam%d.root",ij);
//string1<<Form("remollout_beam%d.root",ij);
add_file[ij-1]=string1.str();
//cout<<string1.str()<<endl;
T.Add(add_file[ij-1]);
}*/
//T.Add("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/ExistingColl/remollout_beam*.root");
Int_t nEvents = T.GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
Double_t weight = 1./nEvents;
weight = 1.;
//TFile f("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield_ep/DSConfig0Shield_ep_disk_28.root","RECREATE");
TFile f(Form("%s",out.Data()),"RECREATE");

std::map<TString,TH2D*> h_xy;
std::map<TString,TH1D*> h_r;
std::map<TString,TH2D*> h_xy1;
std::map<TString,TH1D*> h_phi;
TString part,part1,part2,part3,part4;
const Int_t nParticle=1;
const Int_t nDet=2;

string sParticle[nParticle] = {"electron"};//,"positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{28,"maindet"},{91,"GEMPlane"}};
string sDet[nDet] = {"maindet","GEMPlane"};
string SRing[nRing] = {"Ring1","Ring2","Ring3","Ring4","Ring5","Ring6"};
string SSector[nSector] = {"closed","transition","open"};

for(Int_t ij=0;ij<nDet;ij++)
{
   part=Form("XY_%s",sDet[ij].c_str());
   part1=Form("r_%s",sDet[ij].c_str());
   //part1=Form("EP_phi_%s_%s",SRing[ij].c_str(),SSector[pk].c_str());
   h_xy[part] = new TH2D(part,Form("xy distribution at %s",sDet[ij].c_str()),2400,-1200,1200,2400,-1200,1200);
   h_r[part1] = new TH1D(part1,Form("Radial distribution at %s",sDet[ij].c_str()),700,500,1200);
   //h_xy1[part1] = new TH2D(part,Form("Sector summed xy distribution at %s",sDet[ij].c_str()),2400,-1200,1200,2400,-1200,1200);
   //h_phi[part1] = new TH1D(part1,Form("Phi distribution at main detector for %s %s",SRing[ij].c_str(),SSector[pk].c_str()),100,-5,5);
}

Double_t fRate=0;
remollEvent_t *fEvent=0;
std::vector<remollGenericDetectorHit_t> *fHit=0;
std::vector<remollEventParticle_t> *fPart=0;

T.SetBranchAddress("ev", &fEvent);
T.SetBranchAddress("hit", &fHit);
T.SetBranchAddress("rate", &fRate);
T.SetBranchAddress("part", &fPart);

for(size_t j=0;j<nEvents;j++){
//for(size_t j=0;j<1000;j++){
  T.GetEntry(j);
  std::vector<int> det28trid;
  std::vector<int>::iterator det28it;
  std::vector<int> det91trid;
  std::vector<int>::iterator det91it;
  isValid1(fHit,det28trid,det91trid);
  /*for(it=det28trid.begin();it!=det28trid.end();++it)
        std::cout<<"Inside main event number "<<j<<" track id "<<*it<<std::endl;*/

  if(fRate==0) fRate=1;
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 

    Int_t det = hit.det;
    Int_t pid = hit.pid;
    double e = hit.e;
    Double_t RMax;
    double Phi1;
    //std::cout<<"Before first cont. Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<std::endl;
    if(!((det==28 || det==91) && (abs(pid)==11 || abs(pid)==211 || pid==22)))continue;
    //std::cout<<"After first cont. Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" det28 size "<<det28trid.size()<<" size "<<det91trid.size()<<std::endl;
    Int_t trid = hit.trid;
    det28it = find(det28trid.begin(),det28trid.end(),trid);
    det91it = find(det91trid.begin(),det91trid.end(),trid);
    if(!(det28it!=det28trid.end() && det91it!=det91trid.end()))continue;

    part=Form("XY_%s",snDet[det].c_str());
    part1=Form("r_%s",snDet[det].c_str());
    //part1=Form("XY_sum%s",snDet[det].c_str());
    double x = hit.x;
    double y = hit.y;
    //std::cout<<"After second cont. Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" x "<<x<<" y "<<y<<std::endl;
    h_xy[part]->Fill(hit.x,hit.y,fRate*weight);
    h_r[part1]->Fill(hit.r,fRate*weight);
    /*double phi = atan2(y,x)*rad2deg;
    if(phi<0)phi=phi+2*pi;
    //double phi = atan(py/px)*rad2deg;
    if(abs(phi)>=0 && abs(phi)<2*midangle) Phi1=abs(phi)-midangle;
    if(abs(phi)>=2*midangle && abs(phi)<4*midangle) Phi1=abs(phi)-3*midangle;
    if(abs(phi)>=4*midangle && abs(phi)<6*midangle) Phi1=abs(phi)-5*midangle;
    if(abs(phi)>=6*midangle && abs(phi)<8*midangle) Phi1=abs(phi)-7*midangle;
    if(abs(phi)>=9*midangle && abs(phi)<10*midangle) Phi1=abs(phi)-9*midangle;
    if(abs(phi)>=10*midangle && abs(phi)<12*midangle) Phi1=abs(phi)-11*midangle;
    if(abs(phi)>=12*midangle && abs(phi)<14*midangle) Phi1=abs(phi)-13*midangle;*/
    //std::cout<<"Event "<<j<<" hit "<<i<<" x "<<x<<" y "<<y<<" r "<<hit.r<<" px "<<px<<" py "<<py<<" pz "<<pz<<" sector "<<sector<<" ring "<<ring<<" theta "<<theta<<" Phi1 "<<Phi1<<" RMax "<<RMax<<" Phi "<<Phi<<std::endl;
    }
}

for(Int_t ij=0;ij<nDet;ij++){
    part=Form("XY_%s",sDet[ij].c_str());
    part1=Form("r_%s",sDet[ij].c_str());
    h_xy[part]->Write("",TObject::kOverwrite);    
    h_r[part1]->Write("",TObject::kOverwrite);    
}

return 0;
}

int whichRing(Double_t x, Double_t y, Int_t sector, Double_t &RMax ){
	std::vector<std::vector<double> > rMin, rMax; 
	RMax=-200;
        rMin={
    	  { 640.0,  640.0,  640.0},
      	  { 680.0,  680.0,  680.0},
          { 750.0,  750.0,  750.0},
          { 855.0,  847.5,  825.0},
          { 935.0,  920.0,  875.0},
          {1075.0, 1080.0, 1090.0},
        };
       rMax={
          { 680.0,  680.0,  680.0},
          { 750.0,  750.0,  750.0},
          { 855.0,  847.5,  825.0},
          { 935.0,  920.0,  875.0},
          {1075.0, 1060.0, 1055.0},
          {1190.0, 1190.0, 1190.0},
       };

	double r = sqrt(x*x+y*y);
	int ring=-10;
	for(int ij=0;ij<nRing;ij++)
		if(r >= rMin[ij][sector] && r <=rMax[ij][sector]){
			RMax=rMax[ij][sector];
			return ij;}

	return ring;
}

int whichSector(Double_t x,Double_t y){
	int sector=-10;
	double phi = atan2(y,x);
	if(phi<0) phi+=2*pi;
	const double secPhi = fmod(phi, 2.*pi/7.);
        //0,1,2 == closed, transition, open
	if(secPhi < pi/28.)
		sector = 0;
	else if (secPhi < 3.*pi/28.)
		sector = 1; 
	else if (secPhi < 5.*pi/28.)
		sector = 2; 
	else if (secPhi < 7.*pi/28.)
		sector = 1; 
	else if (secPhi < 8.*pi/28.)
		sector = 0; 

	return sector;
}

void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &det28trid, std::vector<int> &det91trid)
{
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid;
    if(det==28 && hit.r>=640 && hit.r<=1190 && hit.e>100){// && (abs(pid)==11 || abs(pid)==211 || pid==22) && hit.e>500){
    //std::cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	det28trid.push_back(hit.trid);
    }
    if(det==91 && hit.e>100 ){// && abs(pid)==211 && pid==22) && hit.e>500){
    //std::cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	det91trid.push_back(hit.trid);
    }
 }
  //std::cout<<"end of isvalid"<<std::endl;
}
