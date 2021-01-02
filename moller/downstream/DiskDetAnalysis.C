/*
This program is written to calculate power deposition on the downstream coil and epoxy for MOLLER experiment with cuts at different detectors at different z-locations
Chandan Ghosh Dec 2020

*/
#include "remolltypes.hh"
void set_plot_style();
float CalRad(Double_t);
bool isValid(std::vector<remollGenericDetectorHit_t> *fHit=0);
const double pi = acos(-1);
int DiskDetAnalysis()
{
TChain T("T");
const Int_t nfile=1100;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield0_ep/remollout_ep%d.root",ij);
//string1<<Form("remollout_beam%d.root",ij);
add_file[ij-1]=string1.str();
//cout<<string1.str()<<endl;
T.Add(add_file[ij-1]);
}
//T.Add("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/ExistingColl/remollout_beam*.root");
Int_t nEvents = T.GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
Double_t weight = 1./nEvents;
weight = 1.;
//TFile f("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield_ep/DSConfig0Shield_ep_disk_28.root","RECREATE");
TFile f("DSDiskDet_Config0Shield0_ep.root","RECREATE");
set_plot_style();

std::map<TString,TH2D*> h_d_xy;//d==disk
std::map<TString,TH1D*> h_d_e;
std::map<TString,TH2D*> h_t_phiz;//t==tube
std::map<TString,TH2D*> h_t_rz;
std::map<TString,TH1D*> h_t_e;
TString part,part1,part2;
const Int_t nEnergy=5;
//0:Total; 1:hit.p<1; 2: hit.p>=1 && hit.p<10;  3: hit.p>=10 && hit.p<100; 4: hit.p>=100;
const Int_t nParticle=5;
const Int_t nDet=17;
const Int_t nTDet=1;

string sParticle[nParticle] = {"All","electron","positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"},{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{63,"afterCol4"},{64,"beforeCoil1"},{65,"afterCoil1"},{66,"afterCoil2"},{72,"beforeLintel"},{67,"afterCoil3"},{73,"At13M"},{77,"At13p1M"},{76,"At13p4M"},{78,"At13p5M"},{74,"At14M"},{79,"At14p1M"},{75,"At14p5M"},{71,"At14p6M"},{68,"beforeCollar"},{69,"driftRegion"},{28,"maindet"}};
string sDet[nDet] = {"afterCol4","beforeCoil1","afterCoil1","afterCoil2","beforeLintel","afterCoil3","At13M","At13p1M","At13p4M","At13p5M","At14M","At14p1M","At14p5M","At14p6M","beforeCollar","driftRegion","maindet"};
string sTDet[nTDet] = {"cylinder"};
std::map<int,string> TDet  {{57,"cylinder"}};
string sEnergy[nEnergy] = {"Total", "p<1", "p>=1 && p<10", "p>=10 && p<100","p>=100"};


for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  if(i==0){//this is to get total energy deposition from all particles
  part=Form("h_%s",sDet[k].c_str());
  //std::cout<<"Test inside defination loop:i==0  "<<part<<std::endl;
  if(k<=3)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s on %s",sParticle[i].c_str(),sDet[k].c_str()),200,-50, 50, 200,-50,50);
  else if(k>=4 && k<=13)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s on %s",sParticle[i].c_str(),sDet[k].c_str()),600,-150, 150, 600,-150,150);
  else if(k==14)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s on %s",sParticle[i].c_str(),sDet[k].c_str()),1400,-700, 700, 1400,-700,700);
  else if(k==15)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s on %s",sParticle[i].c_str(),sDet[k].c_str()),2000,-1000, 1000, 2000,-1000,1000);
  if(k==16)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s on %s",sParticle[i].c_str(),sDet[k].c_str()),2400,-1200, 1200, 2400,-1200,1200);
  }
  for(int j=0;j<nEnergy && i>0;j++){
  part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
  if(k<=3)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),200,-50, 50, 200,-50,50);
  else if(k>=4 && k<=13)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),600,-150, 150, 600,-150,150);
  else if(k==14)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),1400,-700, 700, 1400,-700,700);
  else if(k==15)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),2000,-1000, 1000, 2000,-1000,1000);
  else if(k==16)
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),2400,-1200, 1200, 2400,-1200,1200);
  }
  part1=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
  //std::cout<<"Test inside defination loop:1D  "<<part1<<std::endl;
  h_d_e[part1] = new TH1D(part1+"_d_e",Form("%s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 110,0,11000);
 }
}

for(Int_t k=0;k<nTDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  if(i==0){//this is to get total energy deposition from all particles
  part=Form("h_%s",sTDet[k].c_str());
  //std::cout<<"Test inside defination loop:i==0  "<<part<<std::endl;
  h_t_rz[part] = new TH2D(part+"_t_rz",Form("%s on %s",sParticle[i].c_str(),sTDet[k].c_str()),625,3300, 15800, 90,0,100);
 /* if(k>=7){
   part=Form("h_%s",sDet[k].c_str());
   h_ub_rz[part]= new TH2D(part+"_ub_rz",Form("Under belly epoxy:%s on %s",sParticle[i].c_str(),sDet[k].c_str()),375,4500, 12000, 5,25,50);
   }*/
  }
  for(int j=0;j<nEnergy && i>0;j++){
  part=Form("h_%s_%s_E%d",sTDet[k].c_str(),sParticle[i].c_str(),j);
  //std::cout<<"Test inside defination loop:i!=0 "<<part<<std::endl;
//  std::cout<<"Inside the histogram defination "<<part<<std::endl;
  h_t_rz[part] = new TH2D(part+"_t_rz",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sTDet[k].c_str()),625,3300, 15800, 90,0,100);
  /*if(k>=7)h_ub_rz[part] = new TH2D(part+"_ub_rz",Form("Under belly epoxy: %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,800, 3200, 5,25,50);*/
  if(k==0){//-180 to +180deg phi angle covers all the coils. So define the histograms once is enough!!
 	part2=Form("h_%s_E%d",sParticle[i].c_str(),j);
        h_t_phiz[part2] = new TH2D(part2+"_t_phiz",Form("On Cylindrical %s for %s",sEnergy[j].c_str(),sParticle[i].c_str()),625,3300, 15800, 380,-190,190);
  }
  }
  part1=Form("h_%s_%s",sTDet[k].c_str(),sParticle[i].c_str());
  //std::cout<<"Test inside defination loop:1D  "<<part1<<std::endl;
  h_t_e[part1] = new TH1D(part1+"_t_e",Form("%s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sTDet[k].c_str()), 110,0,11000);
  //if(k>=7)h_ub_e[part1] = new TH1D(part1+"_ub_e",Form("Under belly epoxy: %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 110,0,11000);
 // part=Form("pr_E%d",k);
 // h_u_xy[part] = new TH2D(part+"_u_xy",Form("%s upstream einc",part.Data()),180,-450, 450, 180, -450,450);
 // h_ue_xy[part] = new TH2D(part+"_ue_xy",Form("%s upstream epoxy einc",part.Data()),180,-450, 450, 180, -450,450);
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
//for(size_t j=0;j<100;j++){
  T.GetEntry(j);
  if(fRate==0) fRate=1;
  bool test = isValid(fHit);
  if(!test) continue;
  //std::cout<<"event number "<<j<<" test "<<test<<std::endl;
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    Bool_t ene_cut;

    std::map<TString, Bool_t> hit_pid;
    Int_t det = hit.det;
    Int_t pid = hit.pid;
   // std::cout<<"Before Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<std::endl;
    //if(((det>=63 && det<=69) || det==28 || (det>=57 && det<=62)) && (abs(pid)==11 || abs(pid)==2112 || pid==22)){
    if(((det>=57 && det<=79) || det==28 ) && (abs(pid)==11) && (hit.trid==1 && hit.mtrid==0) && hit.pz>0){
//std::cout<<"Hi Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<std::endl;
        Int_t hitz = hit.z;
	if((det>=63 && det<=79) || det==28 ){
	//std::cout<<"Event "<<j<<"  filled disk dets Total"<<std::endl;
        part=Form("h_%s",snDet[det].c_str());
	h_d_xy[part]->Fill(hit.y,hit.x,(fRate)*weight);}
	if(det>=57 && det<=62){
	//std::cout<<"Event "<<j<<"  filled tube dets Total"<<std::endl;
        part=Form("h_%s",TDet[57].c_str());
	h_t_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);}
	/*float hitr = CalRad(hitz);
	bool hit_belly = ((det==4008 ||det==4009 ||det==4010 ||det==4011 ||det==4012 ||det==4013 ||det==4014)&&(hit.r<hitr)); 
	//if(hit_belly){part=Form("h_%s",snDet[det].c_str());  h_ub_rz[part]->Fill(hit.z,hit.r,hit.edep*(fRate)*weight);}
	if(hit_belly) h_ub_rz[part]->Fill(hit.z,hit.r,hit.edep*(fRate)*weight);*/
        
//	std::cout<<"hit det  "<<det<<"  hitz  "<<hitz<<"  hit.r  "<<hit.r<<" hitr  "<<hitr<<"  hit_belly  "<<hit_belly<<std::endl;
        for(Int_t k=0;k<nEnergy;k++){
          if(k==0){ene_cut=hit.p<12000;}
          else if(k==1){ene_cut=hit.p<1;}
          else if (k==2) {ene_cut=hit.p>=1 && hit.p<10;}
          else if (k==3) {ene_cut=hit.p>=10 && hit.p<100;}
          else {ene_cut=hit.p>=100;}
	     if(ene_cut){
		if((det>=63 && det<=79)||det==28 ){
        	part=Form("h_%s_%s_E%d",snDet[det].c_str(),snParticle[pid].c_str(),k);
//		std::cout<<"Event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
		//std::cout<<"Event "<<j<<"  Hit  "<<i<<" Energy "<<k<<" part inside hit loop  "<<part<<std::endl;
		h_d_xy[part]->Fill(hit.y,hit.x,(fRate)*weight);}
		if(det>=57 && det<=62){
        	part=Form("h_%s_%s_E%d",TDet[57].c_str(),snParticle[pid].c_str(),k);
//		std::cout<<"Event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
		h_t_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);
		//if(hit_belly){//std::cout<<"Inside hit_belly loop for k= "<<k<<std::endl;
		//h_ub_rz[part]->Fill(hit.z,hit.r,hit.edep*(fRate)*weight);
 	        part2=Form("h_%s_E%d",snParticle[pid].c_str(),k);
		double phi =(180./pi)*atan2(hit.y,hit.x); 
		h_t_phiz[part2]->Fill(hit.z,phi,(fRate)*weight);
		//}
		}
	     }
   	}
	if((det>=63 && det<=79)||det==28){
        part1=Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());
	h_d_e[part1]->Fill(hit.edep,(fRate)*weight);}
	if(det>=57 && det<=62){
        part1=Form("h_%s_%s",TDet[57].c_str(),snParticle[pid].c_str());
	h_t_e[part1]->Fill(hit.edep,(fRate)*weight);
	/*if(hit_belly){//std::cout<<"Inside hit_belly loop"<<std::endl;
	h_ub_e[part1]->Fill(hit.edep,(fRate)*weight);}*/
	}
    }
  }
}


for(Int_t k=0;k<nDet;k++){
  for(Int_t i=0;i<nParticle;i++){
    if(i==0){
    part=Form("h_%s",sDet[k].c_str());
    h_d_xy[part]->Write("", TObject::kOverwrite);
    //std::cout<<"Test inside write loop:2D 0 ::"<<nEnergy<<std::endl;
    //if(k>=7)h_ub_rz[part]->Write("", TObject::kOverwrite);
    }
    
    for(Int_t j=0;j<nEnergy && i>0;j++){
      part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
      h_d_xy[part]->Write("", TObject::kOverwrite);
      //std::cout<<"Test inside write loop:2D"<<std::endl;
      /*if(k>=7)
      h_ub_rz[part]->Write("", TObject::kOverwrite);
      if(k==0){
	      part2=Form("h_%s_E%d",sParticle[i].c_str(),j);
	      h_ub_phiz[part2]->Write("", TObject::kOverwrite);
      }*/
    }
    if(i>0){
    part1=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
    h_d_e[part1]->Write("",TObject::kOverwrite);
    //std::cout<<"Test inside write loop:1D"<<std::endl;
    /*if(k>=7)
    h_ub_e[part1]->Write("",TObject::kOverwrite);*/
    }
  }
}

for(Int_t k=0;k<1;k++){
  for(Int_t i=0;i<nParticle;i++){
    if(i==0){
    part=Form("h_%s",sTDet[k].c_str());
    h_t_rz[part]->Write("", TObject::kOverwrite);
    //std::cout<<"Test inside write loop:2D 0 ::"<<nEnergy<<std::endl;
    //if(k>=7)h_ub_rz[part]->Write("", TObject::kOverwrite);
    }
    
    for(Int_t j=0;j<nEnergy && i>0;j++){
      part=Form("h_%s_%s_E%d",sTDet[k].c_str(),sParticle[i].c_str(),j);
      h_t_rz[part]->Write("", TObject::kOverwrite);
      //std::cout<<"Test inside write loop:2D"<<std::endl;
      /*if(k>=7)
      h_ub_rz[part]->Write("", TObject::kOverwrite);*/
      if(k==0){
	      part2=Form("h_%s_E%d",sParticle[i].c_str(),j);
	      h_t_phiz[part2]->Write("", TObject::kOverwrite);
      }
    }
    if(i>0){
    part1=Form("h_%s_%s",sTDet[k].c_str(),sParticle[i].c_str());
    h_t_e[part1]->Write("",TObject::kOverwrite);
      //std::cout<<"Test inside write loop:1D"<<std::endl;
    /*if(k>=7)
    h_ub_e[part1]->Write("",TObject::kOverwrite);*/
    }
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

bool isValid(std::vector<remollGenericDetectorHit_t> *fHit=0)
{
  bool test1=false;	
  bool test2=false;	
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid;
    if(abs(pid)!=11 || hit.trid!=1 || hit.mtrid!=0 || hit.pz<=0) continue;
    if(det==63 && hit.r<31) test1=true; 
    if(det==28 && hit.r>=640) test2=true;// && hit.r<=1075) test2=true; 
    if(test1&&test2){return true;}
   }

  return false;
}




float CalRad(Double_t hitz)
{
	float slope = 4.615/1792.73;
        float z_offset = 2000; 

	float hitr = (float)(32.466+1.1+(hitz-z_offset)*slope);
	return hitr;
}
