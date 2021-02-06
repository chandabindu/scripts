/*
This program is written to calculate the theta and phi angles at det 28 of high energy electrons for engineers to design the detector support structure. The angles are calculated based on the three momentum components.
Chandan Ghosh Jan 2021

*/
#include "remolltypes.hh"
const double pi = acos(-1);
const double rad2deg = 180./pi;
const double midangle=360./14.;
const int nRing = 6;
const int nSector = 3;
int whichRing(Double_t,Double_t, Int_t, Double_t &RMax);
int whichSector(Double_t,Double_t);
int Angle(TString source, TString out)
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
TString part,part1,part2,part3,part4;
const Int_t nParticle=1;
const Int_t nDet=1;

string sParticle[nParticle] = {"electron"};//,"positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{28,"maindet"}};
string sDet[nDet] = {"maindet"};
string SRing[nRing] = {"Ring1","Ring2","Ring3","Ring4","Ring5","Ring6"};
string SSector[nSector] = {"closed","transition","open"};

for(Int_t ij=0;ij<nRing;ij++)
{
  for(Int_t pk=0;pk<nSector;pk++){
   part=Form("EP_theta_%s_%s",SRing[ij].c_str(),SSector[pk].c_str());
   part1=Form("EP_phi_%s_%s",SRing[ij].c_str(),SSector[pk].c_str());
   h_theta[part] = new TH1D(part,Form("Theta distribution at main detector for %s %s",SRing[ij].c_str(),SSector[pk].c_str()),100,0,5);
   h_phi[part1] = new TH1D(part1,Form("Phi distribution at main detector for %s %s",SRing[ij].c_str(),SSector[pk].c_str()),100,-5,5);
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
//for(size_t j=0;j<1000;j++){
  T.GetEntry(j);

  if(fRate==0) fRate=1;
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    Bool_t ene_cut;

    Int_t det = hit.det;
    Int_t pid = hit.pid;
    double e = hit.e;
    Double_t RMax;
    double Phi1;
    if(det==28 && e>500 && hit.r>=640 && hit.r<1190 && pid==11){ 
    double x = hit.x;
    double y = hit.y;
    double px = hit.px;
    double py = hit.py;
    double pz = hit.pz;
    int sector = whichSector(x,y);
    int ring = whichRing(x,y,sector,RMax);
    if(ring==-10 || sector==-10) continue;
    double theta = atan2(sqrt(px*px+py*py),pz);
    double phi = atan2(py,px)*rad2deg;
    if(phi<0)phi=phi*(-1);
    //double phi = atan(py/px)*rad2deg;
    if(abs(phi)>=0 && abs(phi)<2*midangle) Phi1=abs(phi)-midangle;
    if(abs(phi)>=2*midangle && abs(phi)<4*midangle) Phi1=abs(phi)-3*midangle;
    if(abs(phi)>=4*midangle && abs(phi)<6*midangle) Phi1=abs(phi)-5*midangle;
    if(abs(phi)>=6*midangle && abs(phi)<=7*midangle) Phi1=abs(phi)-7*midangle;
    /*if(Phi>=0 && Phi<2*midangle) Phi=Phi-midangle;
    if(Phi>=2*midangle && Phi<4*midangle) Phi=Phi-3*midangle;
    if(Phi>=4*midangle && Phi<6*midangle) Phi=Phi-5*midangle;
    if(Phi>=6*midangle && Phi<7*midangle) Phi=Phi-7*midangle;
    if(Phi<0 && Phi>=-2*midangle) Phi=Phi+midangle;
    if(Phi<-2*midangle && Phi>=-4*midangle) Phi=Phi+3*midangle;
    if(Phi<-4*midangle && Phi>=-6*midangle) Phi=Phi+5*midangle;
    if(Phi<-6*midangle && Phi>=-7*midangle) Phi=Phi+7*midangle;*/
    //std::cout<<"Ring "<<ring<<" sector "<<sector<<std::endl;
    //Angle this phi section makes at the end of downstream magnet =10m distance   
    double Phi = atan(RMax*sin(Phi1*acos(-1)/180.)/10000.)*rad2deg;
    //std::cout<<"Event "<<j<<" hit "<<i<<" x "<<x<<" y "<<y<<" r "<<hit.r<<" px "<<px<<" py "<<py<<" pz "<<pz<<" sector "<<sector<<" ring "<<ring<<" theta "<<theta<<" Phi1 "<<Phi1<<" RMax "<<RMax<<" Phi "<<Phi<<std::endl;
    part=Form("EP_theta_%s_%s",SRing[ring].c_str(),SSector[sector].c_str());
    part1=Form("EP_phi_%s_%s",SRing[ring].c_str(),SSector[sector].c_str());
    h_theta[part]->Fill(theta*rad2deg,fRate);
    h_phi[part1]->Fill(Phi,fRate);
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
