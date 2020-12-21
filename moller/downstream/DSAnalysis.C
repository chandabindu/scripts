/*
This program is written to calculate power deposition on the upstream coil and epoxy for MOLLER experiment
Chandan Ghosh Oct 2020

*/
#include "remolltypes.hh"
void set_plot_style();
float CalRad(Double_t);
const double pi = acos(-1);
int DSAnalysis()
{
TChain T("T");
const Int_t nfile=1100;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield_ep/remollout_ep%d.root",ij);
//string1<<Form("remollout_beam%d.root",ij);
add_file[ij-1]=string1.str();
//cout<<string1.str()<<endl;
T.Add(add_file[ij-1]);
}
//T.Add("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/ExistingColl/remollout_beam*.root");
Int_t nEvents = T.GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
Double_t weight = 1./nEvents;
TFile f("DSConfig0Shield_ep.root","RECREATE");
set_plot_style();

std::map<TString,TH2D*> h_u_rz;
std::map<TString,TH1D*> h_u_e;
//std::map<TString,TH2D*> h_ub_rz;
std::map<TString,TH2D*> h_t_phiz;
std::map<TString,TH2D*> h_t_rz;
std::map<TString,TH1D*> h_t_e;
//std::map<TString,TH2D*> h_u_xy;
//std::map<TString,TH2D*> h_ue_xy;
//std::map<TString,TH1D*> h_ub_e;
TString part,part1,part2;
const Int_t nEnergy=5;
//0:Total; 1:hit.p<1; 2: hit.p>=1 && hit.p<10;  3: hit.p>=10 && hit.p<100; 4: hit.p>=100;
const Int_t nParticle=5;
const Int_t nDet=14;
const Int_t nTDet=6;

string sParticle[nParticle] = {"All","electron","positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"},{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{3001,"coil1"},{3002,"coil2"},{3003,"coil3"},{3004,"coil4"},{3005,"coil5"},{3006,"coil6"},{3007,"coil7"},{3008,"epoxy1"},{3009,"epoxy2"},{3010,"epoxy3"},{3011,"epoxy4"},{3012,"epoxy5"},{3013,"epoxy6"},{3014,"epoxy7"}};
string sDet[nDet] = {"coil1","coil2","coil3","coil4","coil5","coil6","coil7","epoxy1","epoxy2","epoxy3","epoxy4","epoxy5","epoxy6","epoxy7"};
string sTDet[nTDet] = {"cylinder","tube2","tube3","tube4","tube5","tube6"};
std::map<int,string> TDet  {{57,"cylinder"},{58,"tube2"},{59,"tube3"},{60,"tube4"},{61,"tube5"},{62,"tube6"}};
string sEnergy[nEnergy] = {"Total", "p<1", "p>=1 && p<10", "p>=10 && p<100","p>=100"};


for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  if(i==0){//this is to get total energy deposition from all particles
  part=Form("h_%s",sDet[k].c_str());
  //std::cout<<"Test inside defination loop:i==0  "<<part<<std::endl;
  h_u_rz[part] = new TH2D(part+"_u_rz",Form("%s on %s",sParticle[i].c_str(),sDet[k].c_str()),375,4500, 12000, 90,0,450);
 /* if(k>=7){
   part=Form("h_%s",sDet[k].c_str());
   h_ub_rz[part]= new TH2D(part+"_ub_rz",Form("Under belly epoxy:%s on %s",sParticle[i].c_str(),sDet[k].c_str()),375,4500, 12000, 5,25,50);
   }*/
  }
  for(int j=0;j<nEnergy && i>0;j++){
  part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
  //std::cout<<"Test inside defination loop:i!=0 "<<part<<std::endl;
//  std::cout<<"Inside the histogram defination "<<part<<std::endl;
  h_u_rz[part] = new TH2D(part+"_u_rz",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),375,4500, 12000, 90,0,450);
  /*if(k>=7)h_ub_rz[part] = new TH2D(part+"_ub_rz",Form("Under belly epoxy: %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,800, 3200, 5,25,50);
  if(k==0){//-180 to +180deg phi angle covers all the coils. So define the histograms once is enough!!
 	part2=Form("h_%s_E%d",sParticle[i].c_str(),j);
        h_ub_phiz[part2] = new TH2D(part2+"_ub_phiz",Form("Under belly epoxy: %s for %s",sEnergy[j].c_str(),sParticle[i].c_str()),120,800, 3200, 380,-190,190);
  }*/
  }
  part1=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
  //std::cout<<"Test inside defination loop:1D  "<<part1<<std::endl;
  h_u_e[part1] = new TH1D(part1+"_u_e",Form("%s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 110,0,11000);
  //if(k>=7)h_ub_e[part1] = new TH1D(part1+"_ub_e",Form("Under belly epoxy: %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 110,0,11000);
 // part=Form("pr_E%d",k);
 // h_u_xy[part] = new TH2D(part+"_u_xy",Form("%s upstream einc",part.Data()),180,-450, 450, 180, -450,450);
 // h_ue_xy[part] = new TH2D(part+"_ue_xy",Form("%s upstream epoxy einc",part.Data()),180,-450, 450, 180, -450,450);
 }
}

for(Int_t k=0;k<1;k++){
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
  fRate=1;
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    Bool_t ene_cut;

    std::map<TString, Bool_t> hit_pid;
    Int_t det = hit.det;
    Int_t pid = hit.pid;
    if((det==3001 ||det==3002 ||det==3003 ||det==3004 ||det==3005 ||det==3006 ||det==3007 ||
        det==3008 ||det==3009 ||det==3010 ||det==3011 ||det==3012 ||det==3013 ||det==3014 || 
	det==57 || det==58 || det==59 || det==60 || det==61 ||det==62)&&( 
        abs(pid)==11 || abs(pid)==2112 || pid==22)){
        Int_t hitz = hit.z;
	if(det>=3001 && det<=3014){
        part=Form("h_%s",snDet[det].c_str());
	h_u_rz[part]->Fill(hit.z,hit.r,hit.edep*(fRate)*weight);}
	if(det>=57 && det<=62){
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
		if(det>=3001 && det<=3014){
        	part=Form("h_%s_%s_E%d",snDet[det].c_str(),snParticle[pid].c_str(),k);
		//std::cout<<"hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<std::endl;
		//std::cout<<"Event "<<j<<"  Hit  "<<i<<" Energy "<<k<<" part inside hit loop  "<<part<<std::endl;
		h_u_rz[part]->Fill(hit.z,hit.r,hit.edep*(fRate)*weight);}
		if(det>=57 && det<=62){
        	part=Form("h_%s_%s_E%d",TDet[57].c_str(),snParticle[pid].c_str(),k);
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
	if(det>=3001 && det<=3014){
        part1=Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());
	h_u_e[part1]->Fill(hit.edep,(fRate)*weight);}
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
    h_u_rz[part]->Write("", TObject::kOverwrite);
    //std::cout<<"Test inside write loop:2D 0 ::"<<nEnergy<<std::endl;
    //if(k>=7)h_ub_rz[part]->Write("", TObject::kOverwrite);
    }
    
    for(Int_t j=0;j<nEnergy && i>0;j++){
      part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
      h_u_rz[part]->Write("", TObject::kOverwrite);
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
    h_u_e[part1]->Write("",TObject::kOverwrite);
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


float CalRad(Double_t hitz)
{
	float slope = 4.615/1792.73;
        float z_offset = 2000; 

	float hitr = (float)(32.466+1.1+(hitz-z_offset)*slope);
	return hitr;
}
