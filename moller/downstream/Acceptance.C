/*
This program is written to calculate power deposition on the downstream coil and epoxy for MOLLER experiment with cuts at different detectors at different z-locations
Chandan Ghosh Dec 2020

*/
#include "remolltypes.hh"
void set_plot_style();
const double pi = acos(-1);
int Acceptance()
{
TChain T("T");
const Int_t nfile=1100;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield5_ep/remollout_ep%d.root",ij);
//string1<<Form("remollout_beam%d.root",ij);
add_file[ij-1]=string1.str();
//cout<<string1.str()<<endl;
T.Add(add_file[ij-1]);
}
//T.Add("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/ExistingColl/remollout_beam*.root");
Int_t nEvents = T.GetEntries();
std::cout<<"Number of entries  "<<nEvents<<std::endl;
Double_t weight = 1./nEvents;
weight = 1;
//TFile f("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/DSConfig0Shield_ep/DSConfig0Shield_ep_disk_28.root","RECREATE");
TFile f("Acceptance_Config0Shield5_ep.root","RECREATE");
set_plot_style();

std::map<TString,TH2D*> h_d_xy;//d==disk
std::map<TString,TH2D*> h_d_RTheta;//d==disk
std::map<TString,TH1D*> h_d_e;
std::map<TString,TH1D*> h_d_r;
std::map<TString,TH1D*> h_d_thlab;
TString part,part1,part2;
const Int_t nEnergy=5;
//0:Total; 1:hit.p<1; 2: hit.p>=1 && hit.p<10;  3: hit.p>=10 && hit.p<100; 4: hit.p>=100;
const Int_t nParticle=4;
const Int_t nDet=2;

string sParticle[nParticle] = {"electron","positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"},{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{68,"beforeCollar"},{28,"maindet"}};
string sDet[nDet] = {"beforeCollar","maindet"};
string sEnergy[nEnergy] = {"Total", "p<1", "p>=1 && p<10", "p>=10 && p<100","p>=100"};


for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  for(int j=0;j<nEnergy ;j++){
  part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
  //std::cout<<"Test inside defination loop:2D  "<<part<<std::endl;
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),3000,-1500, 1500, 3000,-1500,1500);
  h_d_RTheta[part] = new TH2D(part+"_d_RTheta",Form("%s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,0, 0.05, 500,0,1500);
  h_d_r[part] = new TH1D(part+"_d_r",Form("Radial:: %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),1500,0,1500);
  h_d_thlab[part] = new TH1D(part+"_d_thlab",Form("theta lab:: %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,0,0.05);
  }
  part1=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
  //std::cout<<"Test inside defination loop:1D  "<<part1<<std::endl;
  h_d_e[part1] = new TH1D(part1+"_d_e",Form("%s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 1200,0,12000);
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
//for(size_t j=0;j<10;j++){
  T.GetEntry(j);
  if(fRate==0) fRate=1;
  //std::cout<<"event number "<<j<<" test "<<test<<std::endl;
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    Bool_t ene_cut;
    
    //if(!(hit.trid==1 && hit.mtrid==0)) continue;
    std::map<TString, Bool_t> hit_pid;
    Int_t det = hit.det;
    Int_t pid = hit.pid;
   // std::cout<<"Before Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<std::endl;
    if((det==28 ||det==68) && (abs(pid)==11 || pid==22 || pid==2112 )) {
//std::cout<<"Hi Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<std::endl;
        Int_t hitz = hit.z;
//	std::cout<<"hit det  "<<det<<"  hitz  "<<hitz<<"  hit.r  "<<hit.r<<" hitr  "<<hitr<<"  hit_belly  "<<hit_belly<<std::endl;
        for(Int_t k=0;k<nEnergy;k++){
          if(k==0){ene_cut=hit.p<12000;}
          else if(k==1){ene_cut=hit.p<1;}
          else if (k==2) {ene_cut=hit.p>=1 && hit.p<10;}
          else if (k==3) {ene_cut=hit.p>=10 && hit.p<100;}
          else {ene_cut=hit.p>=100;}
	     if(ene_cut){
        	part=Form("h_%s_%s_E%d",snDet[det].c_str(),snParticle[pid].c_str(),k);
//		std::cout<<"Event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<" part.th "<<par.th<<std::endl;
//		std::cout<<"Event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
//		std::cout<<"Event "<<j<<"  Hit  "<<i<<" Energy "<<k<<" part inside hit loop  "<<part<<std::endl;
		h_d_xy[part]->Fill(hit.y,hit.x,(fRate)*weight);
		h_d_r[part]->Fill(hit.r,(fRate)*weight);
		if(hit.mtrid==0 && hit.trid==1){
		remollEventParticle_t par=fPart->at(0); 
		if(det==28){
		h_d_thlab[part]->Fill(par.th,fRate*weight*(hit.r>=640));
		h_d_RTheta[part]->Fill(par.th,hit.r,fRate*weight*(hit.r>=640));
		}
		if(det==68){
		h_d_thlab[part]->Fill(par.th,fRate*weight*(hit.r>=70));
		h_d_RTheta[part]->Fill(par.th,hit.r,fRate*weight*(hit.r>=70));
		}
		}

	     }
	}	
        part1=Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());
	h_d_e[part1]->Fill(hit.e,(fRate)*weight);
	}
    }
  }



for(Int_t k=0;k<nDet;k++){
  for(Int_t i=0;i<nParticle;i++){
    for(Int_t j=0;j<nEnergy;j++){
      part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
      h_d_xy[part]->Write("", TObject::kOverwrite);
      h_d_r[part]->Write("", TObject::kOverwrite);
      h_d_thlab[part]->Write("", TObject::kOverwrite);
      h_d_RTheta[part]->Write("", TObject::kOverwrite);
      //std::cout<<"Test inside write loop:2D"<<std::endl;
      /*if(k>=7)
      h_ub_rz[part]->Write("", TObject::kOverwrite);
      if(k==0){
	      part2=Form("h_%s_E%d",sParticle[i].c_str(),j);
	      h_ub_phiz[part2]->Write("", TObject::kOverwrite);
      }*/
    }
    part1=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
    h_d_e[part1]->Write("",TObject::kOverwrite);
    //std::cout<<"Test inside write loop:1D"<<std::endl;
    /*if(k>=7)
    h_ub_e[part1]->Write("",TObject::kOverwrite);*/
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
