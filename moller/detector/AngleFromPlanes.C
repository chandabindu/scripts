/*
This program is written to calculate the theta and phi angles at det 28 of high energy electrons for engineers to design the detector support structure. The angles are calculated based on the three momentum components.
Chandan Ghosh Jan 2021

*/
#include "remolltypes.hh"
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &det28trid, std::vector<int> &det68trid);
const double pi = acos(-1);
const int nRing = 6;
const int nSector = 3;
const double rad2deg = 180./pi;
int whichRing(Double_t,Double_t, Int_t, Double_t &RMax);
int whichSector(Double_t,Double_t);
int AngleFromPlanes(TString source, TString out)
{
TChain T("T");
T.Add(Form("%s",source.Data()));
/*const Int_t nfile=100;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/moller12gev/chandan/sim_out/SymmetricCoil3_ep/remollout_ep%d.root",ij);
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

std::map<TString,TH1D*> h_theta;
std::map<TString,TH1D*> h_phi;
TString part,part1;

string SRing[nRing] = {"Ring1","Ring2","Ring3","Ring4","Ring5","Ring6"};
string SSector[nSector] = {"closed","transition","open"};

for(Int_t ij=0;ij<nRing;ij++)
{
  for(Int_t pk=0;pk<nSector;pk++){
   part=Form("EP_theta_%s_%s",SRing[ij].c_str(),SSector[pk].c_str());
   part1=Form("EP_phi_%s_%s",SRing[ij].c_str(),SSector[pk].c_str());
   h_theta[part] = new TH1D(part,Form("Theta distribution at main detector for %s %s",SRing[ij].c_str(),SSector[pk].c_str()),100,0,5);
   h_phi[part1] = new TH1D(part1,Form("Phi distribution at main detector for %s %s",SRing[ij].c_str(),SSector[pk].c_str()),90,-15,15);
   }
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
  T.GetEntry(j);
  double x1{0},y1{0},z1{0},x2{0},y2{0},z2{0};
  double phi1{-100},phi2{-100};
  std::vector<int> det28trid;
  std::vector<int>::iterator det28it;
  std::vector<int> det68trid;
  std::vector<int>::iterator det68it;
  isValid1(fHit,det28trid,det68trid);
  /*for(it=det28trid.begin();it!=det28trid.end();++it)
        std::cout<<"Inside main event number "<<j<<" track id "<<*it<<std::endl;*/
  bool flag1=false;
  bool flag2=false;
  int ring{-20},sector{-20};
  if(fRate==0) fRate=1;
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    Bool_t ene_cut;

    Int_t det = hit.det;
    Int_t pid = hit.pid;
    double e = hit.e;
    Double_t RMax;
    if(!(det==28 || det==68 || pid==11))continue;
    Int_t trid = hit.trid;
    det28it = find(det28trid.begin(),det28trid.end(),trid);
    det68it = find(det68trid.begin(),det68trid.end(),trid);
    if(!(det28it!=det28trid.end() && det68it!=det68trid.end()))continue;
    if(det==28){
    x2 = hit.x;
    y2 = hit.y;
    z2 = hit.z;
    phi2=atan2(y2,x2);
    double px = hit.px;
    double py = hit.py;
    double pz = hit.pz;
    sector = whichSector(x2,y2);
    ring = whichRing(x2,y2,sector,RMax);
    std::cout<<"Det "<<det<<"  Event "<<j<<" hit "<<i<<" x2 "<<x2<<" y2 "<<y2<<" z2 "<<z2<<" r "<<hit.r<<" sector "<<sector<<" ring "<<ring<<" phi2 "<<phi2<<" RMax "<<RMax<<std::endl;
    if(ring==-10 || sector==-10){flag1=false; continue;}
    flag1=true;}
    if(det==68){
    x1 = hit.x;
    y1 = hit.y;
    z1 = hit.z;
    phi1=atan2(y1,x1);
    std::cout<<"Det "<<det<<"  Event "<<j<<" hit "<<i<<" x1 "<<x1<<" y1 "<<y1<<" z1 "<<z1<<" r "<<hit.r<<" sector "<<sector<<" ring "<<ring<<" phi1 "<<phi1<<std::endl;
    flag2=true;}
    if(flag1&&flag2){
    double theta = atan2((sqrt(x2*x2+y2*y2)-sqrt(x1*x1+y1*y1)),(z2-z1));
    double Phi = phi2-phi1;
    //double theta = atan2(sqrt(px*px+py*py),pz);
    //double Phi = atan2(py,px);
    //std::cout<<"Ring "<<ring<<" sector "<<sector<<std::endl;
    std::cout<<"Inside AND condition  Event "<<j<<" hit "<<i<<" r "<<hit.r<<" sector "<<sector<<" ring "<<ring<<" theta "<<theta<<" Phi "<<Phi<<" phi2 "<<phi2<<" phi1 "<<phi1<<std::endl;
    part=Form("EP_theta_%s_%s",SRing[ring].c_str(),SSector[sector].c_str());
    part1=Form("EP_phi_%s_%s",SRing[ring].c_str(),SSector[sector].c_str());
    h_theta[part]->Fill(theta*rad2deg,fRate);
    h_phi[part1]->Fill(Phi*rad2deg,fRate);
    flag1=false;flag2=false;
    }
    }
}

for(Int_t ij=0;ij<nRing;ij++){
  for(Int_t pk=0;pk<nSector;pk++){
    part=Form("EP_theta_%s_%s",SRing[ij].c_str(),SSector[pk].c_str());
    part1=Form("EP_phi_%s_%s",SRing[ij].c_str(),SSector[pk].c_str());
    h_theta[part]->Write("",TObject::kOverwrite);    
    h_phi[part1]->Write("",TObject::kOverwrite);    
    }
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

void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &det28trid, std::vector<int> &det68trid)
{
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid;
    if(det==28 && hit.r>=640 && hit.r<=1190 && pid==11 && hit.e>500){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	det28trid.push_back(hit.trid);
    }
    if(det==68 && pid==11 && hit.e>500){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	det68trid.push_back(hit.trid);
    }
 }
 // std::cout<<"end of isvalid"<<std::endl;
}
