/*
This program is written to calculate power deposition on the downstream coil and epoxy for MOLLER experiment with cuts at different detectors at different z-locations
Chandan Ghosh Dec 2020

*/
#include "remolltypes.hh"
void set_plot_style();
void isValid(std::vector<remollGenericDetectorHit_t> *fHit, bool *test1, bool *test2, bool *test3, bool *test4);
void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &det28trid, std::vector<int> &Ring5trid, std::vector<int> &Boretrid, std::vector<int> &Acceptancetrid);
const double pi = acos(-1);
int ShiftedCoil_BoreVsAcceptance_2dets(TString source, TString out)
{
TChain T("T");
T.Add(Form("%s",source.Data()));
/*const Int_t nfile=1100;
TString add_file[nfile]={""};
for(int ij=1;ij<=nfile;ij++){
ostringstream string1;
string1<<Form("/lustre19/expphy/volatile/halla/parity/chandan/sim_out/OffsetCoil2_ep/remollout_ep%d.root",ij);
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
set_plot_style();

std::map<TString,TH2D*> h_d_xy;//d==disk with out any cut radial redial distributions
std::map<TString,TH2D*> h_d_re;
std::map<TString,TH1D*> h_d_e;
std::map<TString,TH1D*> h_d_r;
std::map<TString,TH2D*> h_t_phiz;//t==tube
std::map<TString,TH2D*> h_t_rz;
std::map<TString,TH1D*> h_t_e;

std::map<TString,TH2D*> h_d1_xy;//d==disk; Maindet acceptance all radii
std::map<TString,TH2D*> h_d1_re;
std::map<TString,TH1D*> h_d1_e;
std::map<TString,TH1D*> h_d1_r;
std::map<TString,TH2D*> h_t1_phiz;//t==tube
std::map<TString,TH2D*> h_t1_rz;
std::map<TString,TH1D*> h_t1_e;

std::map<TString,TH2D*> h_d2_xy;//d==disk; Maindet acceptance  -through bore
std::map<TString,TH2D*> h_d2_re;
std::map<TString,TH1D*> h_d2_e;
std::map<TString,TH1D*> h_d2_r;
std::map<TString,TH2D*> h_t2_phiz;//t==tube
std::map<TString,TH2D*> h_t2_rz;
std::map<TString,TH1D*> h_t2_e;

std::map<TString,TH2D*> h_d3_xy;//d==disk; Ring 5 - all radii
std::map<TString,TH2D*> h_d3_re;
std::map<TString,TH1D*> h_d3_e;
std::map<TString,TH1D*> h_d3_r;
std::map<TString,TH2D*> h_t3_phiz;//t==tube
std::map<TString,TH2D*> h_t3_rz;
std::map<TString,TH1D*> h_t3_e;

std::map<TString,TH2D*> h_d4_xy;//d==disk; Ring 5 - trhough bore
std::map<TString,TH2D*> h_d4_re;
std::map<TString,TH1D*> h_d4_e;
std::map<TString,TH1D*> h_d4_r;
std::map<TString,TH2D*> h_t4_phiz;//t==tube
std::map<TString,TH2D*> h_t4_rz;
std::map<TString,TH1D*> h_t4_e;

TString part,part1,part2,part3,part4;
const Int_t nEnergy=1;
//0:Total; 1:hit.p<1; 2: hit.p>=1 && hit.p<10;  3: hit.p>=10 && hit.p<100; 4: hit.p>=100;
const Int_t nParticle=1;
const Int_t nDet=2;
const Int_t nTDet=1;

string sParticle[nParticle] = {"electron"};//,"positron","photon","neutron"};
std::map<int,string> snParticle  {{11,"electron"}};//,{-11,"positron"},{22,"photon"},{2112,"neutron"}};
std::map<int,string> snDet  {{63,"afterCol4"},{28,"maindet"}};
string sDet[nDet] = {"afterCol4","maindet"};
string sTDet[nTDet] = {"cylinder"};
std::map<int,string> TDet  {{57,"cylinder"}};
string sEnergy[nEnergy] = {"Total"};//, "p<1", "p>=1 && p<10", "p>=10 && p<100","p>=100"};


for(Int_t k=0;k<nDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  for(int j=0;j<nEnergy ;j++){
  part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
  if(k<=0){
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("No cut %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,-500, 500, 100,-500,500);
  h_d1_xy[part] = new TH2D(part+"_d1_xy",Form("Maindet all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,-500, 500, 100,-500,500);
  h_d2_xy[part] = new TH2D(part+"_d2_xy",Form("Maindet from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,-100, 100, 100,-100,100);
  h_d3_xy[part] = new TH2D(part+"_d3_xy",Form("Ring 5 all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,-500, 500, 100,-500,500);
  h_d4_xy[part] = new TH2D(part+"_d4_xy",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,-100, 100, 100,-100,100);
  h_d_re[part] = new TH2D(part+"_d_re",Form("No cut %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 50,0,500);
  h_d1_re[part] = new TH2D(part+"_d1_re",Form("Maindet all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 50,0,500);
  h_d2_re[part] = new TH2D(part+"_d2_re",Form("Maindet from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 50,0,500);
  h_d3_re[part] = new TH2D(part+"_d3_re",Form("Ring 5 all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 100,0,500);
  h_d4_re[part] = new TH2D(part+"_d4_re",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 50,0,500);
  h_d_r[part] = new TH1D(part+"_d_r",Form("No cut %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,0,500);
  h_d1_r[part] = new TH1D(part+"_d1_r",Form("Maindet all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,0,500);
  h_d2_r[part] = new TH1D(part+"_d2_r",Form("Maindet from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,0,100);
  h_d3_r[part] = new TH1D(part+"_d3_r",Form("Ring 5 all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,0,500);
  h_d4_r[part] = new TH1D(part+"_d4_r",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),100,0,100);}
  else if(k==1){
  h_d_xy[part] = new TH2D(part+"_d_xy",Form("No cut %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),1500,-1500, 1500, 1500,-1500,1500);
  h_d1_xy[part] = new TH2D(part+"_d1_xy",Form("Maindet all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),1500,-1500, 1500, 1500,-1500,1500);
  h_d2_xy[part] = new TH2D(part+"_d2_xy",Form("Maindet from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),1500,-1500, 1500, 1500,-1500,1500);
  h_d3_xy[part] = new TH2D(part+"_d3_xy",Form("Ring 5 all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),1500,-1500, 1500, 1500,-1500,1500);
  h_d4_xy[part] = new TH2D(part+"_d4_xy",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),1500,-1500, 1500, 1500,-1500,1500);
  h_d_re[part] = new TH2D(part+"_d_re",Form("No cut %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 70,600,1300);
  h_d1_re[part] = new TH2D(part+"_d1_re",Form("Maindet all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 70,600,1300);
  h_d2_re[part] = new TH2D(part+"_d2_re",Form("Maindet frim col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 70,600,1300);
  h_d3_re[part] = new TH2D(part+"_d3_re",Form("Ring5 all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 70,600,1300);
  h_d4_re[part] = new TH2D(part+"_d4_re",Form("Ring5 fromcol4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()),120,0, 12000, 170,600,1300);
  h_d_r[part] = new TH1D(part+"_d_r",Form("No cut %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 800,500,1300);
  h_d1_r[part] = new TH1D(part+"_d1_r",Form("Maindet all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 800,500,1300);
  h_d2_r[part] = new TH1D(part+"_d2_r",Form("Maindet from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 800,500,1300);
  h_d3_r[part] = new TH1D(part+"_d3_r",Form("Ring 5 all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 800,500,1300);
  h_d4_r[part] = new TH1D(part+"_d4_r",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 800,500,1300);}
  }
  part1=Form("h_%s_%s",sDet[k].c_str(),sParticle[i].c_str());
  //std::cout<<"Test inside defination loop:1D  "<<part1<<std::endl;
  h_d_e[part1] = new TH1D(part1+"_d_e",Form("No cut %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 1200,0,12000);
  h_d1_e[part1] = new TH1D(part1+"_d1_e",Form("Maindet all %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 1200,0,12000);
  h_d2_e[part1] = new TH1D(part1+"_d2_e",Form("Maindet from col4 bore %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 1200,0,12000);
  h_d3_e[part1] = new TH1D(part1+"_d3_e",Form("Ring 5 all %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 1200,0,12000);
  h_d4_e[part1] = new TH1D(part1+"_d4_e",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sDet[k].c_str()), 1200,0,12000);
 }
}

for(Int_t k=0;k<nTDet;k++){
 for(Int_t i=0;i<nParticle;i++){
  for(int j=0;j<nEnergy ;j++){
  part=Form("h_%s_%s_E%d",sTDet[k].c_str(),sParticle[i].c_str(),j);
  //std::cout<<"Test inside defination loop:i!=0 "<<part<<std::endl;
  h_t_rz[part] = new TH2D(part+"_t_rz",Form("No cut %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sTDet[k].c_str()),625,3300, 15800, 100,0,100);
  h_t1_rz[part] = new TH2D(part+"_t1_rz",Form("Maindet all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sTDet[k].c_str()),625,3300, 15800, 100,0,100);
  h_t2_rz[part] = new TH2D(part+"_t2_rz",Form("Maindet from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sTDet[k].c_str()),625,3300, 15800, 100,0,100);
  h_t3_rz[part] = new TH2D(part+"_t3_rz",Form("Ring 5 all %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sTDet[k].c_str()),625,3300, 15800, 90,0,100);
  h_t4_rz[part] = new TH2D(part+"_t4_rz",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[j].c_str(),sParticle[i].c_str(),sTDet[k].c_str()),625,3300, 15800, 100,0,100);
  if(k==0){//-180 to +180deg phi angle covers all the coils. So define the histograms once is enough!!
 	part2=Form("h_%s_E%d",sParticle[i].c_str(),j);
        h_t_phiz[part2] = new TH2D(part2+"_t_phiz",Form("No cut On Cylindrical %s for %s",sEnergy[j].c_str(),sParticle[i].c_str()),625,3300, 15800, 380,-190,190);
        h_t1_phiz[part2] = new TH2D(part2+"_t1_phiz",Form("Maindet all On Cylindrical %s for %s",sEnergy[j].c_str(),sParticle[i].c_str()),625,3300, 15800, 380,-190,190);
        h_t2_phiz[part2] = new TH2D(part2+"_t2_phiz",Form("Maindet from col4 bore On Cylindrical %s for %s",sEnergy[j].c_str(),sParticle[i].c_str()),625,3300, 15800, 380,-190,190);
        h_t3_phiz[part2] = new TH2D(part2+"_t3_phiz",Form("Ring 5 all On Cylindrical %s for %s",sEnergy[j].c_str(),sParticle[i].c_str()),625,3300, 15800, 380,-190,190);
        h_t4_phiz[part2] = new TH2D(part2+"_t4_phiz",Form("Ring 5 from col4 bore On Cylindrical %s for %s",sEnergy[j].c_str(),sParticle[i].c_str()),625,3300, 15800, 380,-190,190);
  }
  }
  part1=Form("h_%s_%s",sTDet[k].c_str(),sParticle[i].c_str());
  //std::cout<<"Test inside defination loop:1D  "<<part1<<std::endl;
  h_t_e[part1] = new TH1D(part1+"_t_e",Form("No cut %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sTDet[k].c_str()), 1200,0,12000);
  h_t1_e[part1] = new TH1D(part1+"_t1_e",Form("Maindet all %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sTDet[k].c_str()), 1200,0,12000);
  h_t2_e[part1] = new TH1D(part1+"_t2_e",Form("Maindet from col4 bore %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sTDet[k].c_str()), 1200,0,12000);
  h_t3_e[part1] = new TH1D(part1+"_t3_e",Form("Ring 5 all %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sTDet[k].c_str()), 1200,0,12000);
  h_t4_e[part1] = new TH1D(part1+"_t4_e",Form("Ring 5 from col4 bore %s %s for %s",sEnergy[0].c_str(),sParticle[i].c_str(),sTDet[k].c_str()), 1200,0,12000);
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

  std::vector<int> det28trid;
  std::vector<int>::iterator it;
  std::vector<int> Ring5trid;
  std::vector<int>::iterator Ring5it;
  std::vector<int> Boretrid;
  std::vector<int>::iterator Boreit;
  std::vector<int> Acceptancetrid;
  std::vector<int>::iterator Acceptanceit;
  if(fRate==0) fRate=1;
  //isValid(fHit,&test1,&test2,&test3,&test4);
  isValid1(fHit,det28trid,Ring5trid,Boretrid,Acceptancetrid);
  /*for(it=det28trid.begin();it!=det28trid.end();++it)
        std::cout<<"Inside main event number "<<j<<" track id "<<*it<<std::endl;*/
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    Bool_t ene_cut;

    std::map<TString, Bool_t> hit_pid;
    Int_t det = hit.det;
    Int_t pid = hit.pid;
    Int_t trid = hit.trid;
    it = find(det28trid.begin(),det28trid.end(),trid);
    Ring5it = find(Ring5trid.begin(),Ring5trid.end(),trid);
    Boreit = find(Boretrid.begin(),Boretrid.end(),trid);
    Acceptanceit = find(Acceptancetrid.begin(),Acceptancetrid.end(),trid);
    //Ring5it=it;
    if(!(it!=det28trid.end() && Acceptanceit!=Acceptancetrid.end())) {/*std::cout<<"Inside Continue loop::Event "<<j<<" hit "<<i<<" det "<<det<<" trid "<<trid<<std::endl;*/continue;}
  //  std::cout<<"After Continue :: Event "<<j<<" hit "<<i<<" trid "<<trid<<" iterator find "<<*it<<std::endl;  
    //std::cout<<"Before Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<std::endl;
    //if(((det>=63 && det<=69) || det==28 || (det>=57 && det<=62)) && (abs(pid)==11 || abs(pid)==2112 || pid==22))
    //if(((det>=57 && det<=79) || det==28 ) && (abs(pid)==11 || pid==22 || abs(pid)==2112) && (hit.trid==1 && hit.mtrid==0) && hit.pz>0){
    if(((det>=57 && det<=63) || det==28 ) && pid==11 && hit.vz<=-3875){
//	std::cout<<"Hi Event "<<j<<" hit "<<i<<" det "<<det<<" pid "<<pid<<" r "<<hit.r<<" trid "<<hit.trid<<" mtrid  "<<hit.mtrid<<std::endl;
        
//	std::cout<<"hit det  "<<det<<"  hitz  "<<hit.z<<"  hit.r  "<<hit.r<<" hitr  "<<hitr<<"  hit_belly  "<<hit_belly<<std::endl;
        for(Int_t k=0;k<nEnergy;k++){
          ene_cut=hit.p<12000 && hit.p>=0;
          /*if(k==0){ene_cut=hit.p<12000 && hit.p>=0;}
          else if(k==1){ene_cut=hit.p<1;}
          else if (k==2) {ene_cut=hit.p>=1 && hit.p<10;}
          else if (k==3) {ene_cut=hit.p>=10 && hit.p<100;}
          else {ene_cut=hit.p>=100;}*/
	     if(ene_cut){
		if(det==63 ||det==28 ){
        	part=Form("h_%s_%s_E%d",snDet[det].c_str(),snParticle[pid].c_str(),k);
//		std::cout<<"Event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
		//std::cout<<"Event "<<j<<"  Hit  "<<i<<" Energy "<<k<<" part inside hit loop  "<<part<<std::endl;
		h_d_xy[part]->Fill(hit.x,hit.y,(fRate)*weight);
		h_d_r[part]->Fill(hit.r,(fRate)*weight);
		h_d_re[part]->Fill(hit.e,hit.r,(fRate)*weight);
		if(Acceptanceit!=Acceptancetrid.end() && it!=det28trid.end())
		if((trid==*Acceptanceit) && (trid==*it)){
		h_d1_xy[part]->Fill(hit.x,hit.y,(fRate)*weight);
		h_d1_r[part]->Fill(hit.r,(fRate)*weight);
		h_d1_re[part]->Fill(hit.e,hit.r,(fRate)*weight);}
//		std::cout<<"Event1  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
		if(Boreit!=Boretrid.end() && it!=det28trid.end())
		if((trid==*Boreit) && (trid==*it)){
		h_d2_xy[part]->Fill(hit.x,hit.y,(fRate)*weight);
		h_d2_r[part]->Fill(hit.r,(fRate)*weight);
		h_d2_re[part]->Fill(hit.e,hit.r,(fRate)*weight);}
		if(Acceptanceit!=Acceptancetrid.end() && Ring5it!=Ring5trid.end())
		if((trid==*Acceptanceit) && (trid==*Ring5it)){
		h_d3_xy[part]->Fill(hit.x,hit.y,(fRate)*weight);
		h_d3_r[part]->Fill(hit.r,(fRate)*weight);
		h_d3_re[part]->Fill(hit.e,hit.r,(fRate)*weight);}
		if(Boreit!=Boretrid.end() && Ring5it!=Ring5trid.end())
		if((trid==*Boreit) && (trid==*Ring5it)){
		h_d4_xy[part]->Fill(hit.x,hit.y,(fRate)*weight);
		h_d4_r[part]->Fill(hit.r,(fRate)*weight);
		h_d4_re[part]->Fill(hit.e,hit.r,(fRate)*weight);}
		}
		if(det>=57 && det<=62){
        	part=Form("h_%s_%s_E%d",TDet[57].c_str(),snParticle[pid].c_str(),k);
 	        part2=Form("h_%s_E%d",snParticle[pid].c_str(),k);
		double phi =(180./pi)*atan2(hit.y,hit.x); 
//		std::cout<<"Event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
		//h_t_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);
		h_t_phiz[part2]->Fill(hit.z,phi,(fRate)*weight);
		if(Acceptanceit!=Acceptancetrid.end() && it!=det28trid.end())
		if((trid==*Acceptanceit) && (trid==*it)){
		h_t1_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);
		h_t1_phiz[part2]->Fill(hit.z,phi,(fRate)*weight);}
		if(Boreit!=Boretrid.end() && it!=det28trid.end())
		if((trid==*Boreit) && (trid==*it)){
		h_t2_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);
		h_t2_phiz[part2]->Fill(hit.z,phi,(fRate)*weight);}
		if(Acceptanceit!=Acceptancetrid.end() && Ring5it!=Ring5trid.end())
		if((trid==*Acceptanceit) && (trid==*Ring5it)){
		h_t3_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);
		h_t3_phiz[part2]->Fill(hit.z,phi,(fRate)*weight);}
		if(Boreit!=Boretrid.end() && Ring5it!=Ring5trid.end())
		if((trid==*Boreit) && (trid==*Ring5it)){
		h_t4_rz[part]->Fill(hit.z,hit.r,(fRate)*weight);
		h_t4_phiz[part2]->Fill(hit.z,phi,(fRate)*weight);}
		//std::cout<<"Inside hit_belly loop for k= "<<k<<std::endl;
		//h_ub_rz[part]->Fill(hit.z,hit.r,hit.edep*(fRate)*weight);
		//}
		}
//		std::cout<<"End of event  "<<j<<"  hit.det  "<<det<<"  pid   "<<pid<<"  p  "<<hit.p<<"  r  "<<hit.r<<std::endl;
	     }
   	}
	if(det==63 || det==28){
        part1=Form("h_%s_%s",snDet[det].c_str(),snParticle[pid].c_str());
	h_d_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Acceptanceit!=Acceptancetrid.end() && it!=det28trid.end())
	if((trid==*Acceptanceit) && (trid==*it))
	h_d1_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Boreit!=Boretrid.end() && it!=det28trid.end())
	if((trid==*Boreit) && (trid==*it))
	h_d2_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Acceptanceit!=Acceptancetrid.end() && Ring5it!=Ring5trid.end())
	if((trid==*Acceptanceit) && (trid==*Ring5it))
	h_d3_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Boreit!=Boretrid.end() && Ring5it!=Ring5trid.end())
	if((trid==*Boreit) && (trid==*Ring5it))
	h_d4_e[part1]->Fill(hit.e,(fRate)*weight);}
	if(det>=57 && det<=62){
        part1=Form("h_%s_%s",TDet[57].c_str(),snParticle[pid].c_str());
	h_t_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Acceptanceit!=Acceptancetrid.end() && it!=det28trid.end())
	if((trid==*Acceptanceit) && (trid==*it))
	h_t1_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Boreit!=Boretrid.end() && it!=det28trid.end())
	if((trid==*Boreit) && (trid==*it))
	h_t2_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Acceptanceit!=Acceptancetrid.end() && Ring5it!=Ring5trid.end())
	if((trid==*Acceptanceit) && (trid==*Ring5it))
	h_t3_e[part1]->Fill(hit.e,(fRate)*weight);
	if(Boreit!=Boretrid.end() && Ring5it!=Ring5trid.end())
	if((trid==*Boreit) && (trid==*Ring5it))
	h_t4_e[part1]->Fill(hit.e,(fRate)*weight);
	}
    }
  }
}


for(Int_t k=0;k<nDet;k++){
  for(Int_t i=0;i<nParticle;i++){
    for(Int_t j=0;j<nEnergy;j++){
      part=Form("h_%s_%s_E%d",sDet[k].c_str(),sParticle[i].c_str(),j);
      h_d_xy[part]->Write("", TObject::kOverwrite);
      h_d1_xy[part]->Write("", TObject::kOverwrite);
      h_d2_xy[part]->Write("", TObject::kOverwrite);
      h_d3_xy[part]->Write("", TObject::kOverwrite);
      h_d4_xy[part]->Write("", TObject::kOverwrite);
      h_d_re[part]->Write("", TObject::kOverwrite);
      h_d1_re[part]->Write("", TObject::kOverwrite);
      h_d2_re[part]->Write("", TObject::kOverwrite);
      h_d3_re[part]->Write("", TObject::kOverwrite);
      h_d4_re[part]->Write("", TObject::kOverwrite);
      h_d_r[part]->Write("", TObject::kOverwrite);
      h_d1_r[part]->Write("", TObject::kOverwrite);
      h_d2_r[part]->Write("", TObject::kOverwrite);
      h_d3_r[part]->Write("", TObject::kOverwrite);
      h_d4_r[part]->Write("", TObject::kOverwrite);
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
    h_d1_e[part1]->Write("",TObject::kOverwrite);
    h_d2_e[part1]->Write("",TObject::kOverwrite);
    h_d3_e[part1]->Write("",TObject::kOverwrite);
    h_d4_e[part1]->Write("",TObject::kOverwrite);
    //std::cout<<"Test inside write loop:1D"<<std::endl;
    /*if(k>=7)
    h_ub_e[part1]->Write("",TObject::kOverwrite);*/
   
  }
}

for(Int_t k=0;k<1;k++){
  for(Int_t i=0;i<nParticle;i++){
    for(Int_t j=0;j<nEnergy;j++){
      part=Form("h_%s_%s_E%d",sTDet[k].c_str(),sParticle[i].c_str(),j);
      h_t_rz[part]->Write("", TObject::kOverwrite);
      h_t1_rz[part]->Write("", TObject::kOverwrite);
      h_t2_rz[part]->Write("", TObject::kOverwrite);
      h_t3_rz[part]->Write("", TObject::kOverwrite);
      h_t4_rz[part]->Write("", TObject::kOverwrite);
      //std::cout<<"Test inside write loop:2D"<<std::endl;
      /*if(k>=7)
      h_ub_rz[part]->Write("", TObject::kOverwrite);*/
      if(k==0){
	      part2=Form("h_%s_E%d",sParticle[i].c_str(),j);
	      h_t_phiz[part2]->Write("", TObject::kOverwrite);
	      h_t1_phiz[part2]->Write("", TObject::kOverwrite);
	      h_t2_phiz[part2]->Write("", TObject::kOverwrite);
	      h_t3_phiz[part2]->Write("", TObject::kOverwrite);
	      h_t4_phiz[part2]->Write("", TObject::kOverwrite);
      }
    }
    
    part1=Form("h_%s_%s",sTDet[k].c_str(),sParticle[i].c_str());
    h_t_e[part1]->Write("",TObject::kOverwrite);
    h_t1_e[part1]->Write("",TObject::kOverwrite);
    h_t2_e[part1]->Write("",TObject::kOverwrite);
    h_t3_e[part1]->Write("",TObject::kOverwrite);
    h_t4_e[part1]->Write("",TObject::kOverwrite);
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

void isValid(std::vector<remollGenericDetectorHit_t> *fHit, bool *test1, bool *test2, bool *test3, bool *test4)
{
  bool test11=false;	
  bool test21=false;	
  bool test12=false;	
  bool test22=false;	
  bool test13=false;	
  bool test23=false;	
  bool test14=false;	
  bool test24=false;	
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid;
   //std::cout<<"Inside isValid before Det "<<det<<" pid "<< pid <<" radius "<<hit.r<<" mtrid "<<hit.mtrid<<"  trid "<<hit.trid<<" pz "<<hit.pz<<" test1 "<<*test1<<" test2 "<<*test2<<" test3 "<<*test3<<" test4 "<<*test4<<std::endl;  
    //if((det!=28||det!=63) && (abs(pid)!=11 || pid!=22 || hit.trid!=1 || hit.mtrid!=0 || hit.pz<=0)) continue;
    //if(!((det==28 || det==63)&&(abs(pid)==11 || abs(pid)==22|| pid==2112||hit.mtrid==0||hit.trid==1||hit.pz>0))){/*std::cout<<"strange det "<<det<<" condiftion "<<(bool)(!(det==28||det==63))<<std::endl;*/ continue;}
    if(!((det==28 || det==63) && abs(pid)==11 && hit.mtrid==0 && hit.trid==1 && hit.pz>0)){/*std::cout<<"strange det "<<det<<" condiftion "<<(bool)(!(det==28||det==63))<<std::endl;*/ continue;}
    //if(abs(pid)!=11 || pid!=22 || hit.trid!=1 || hit.mtrid!=0 || hit.pz<=0) continue;
    if(det==63 && hit.r<200) test11=true; 
    if(det==28 && hit.r>=640 && hit.r<=1200) test21=true;// && hit.r<=1075) test2=true; 

    if(det==63 && hit.r<31) test12=true; 
    if(det==28 && hit.r>=640 && hit.r<=1200) test22=true;// && hit.r<=1075) test2=true; 

    if(det==63 && hit.r<200) test13=true; 
    if(det==28 && hit.r>=935 && hit.r<=1075) test23=true; 

    if(det==63 && hit.r<31) test14=true; 
    if(det==28 && hit.r>=935 && hit.r<=1075) test24=true;
   //std::cout<<"Inside isValid after Det "<<det<<" pid "<< pid <<" radius "<<hit.r<<" mtrid "<<hit.mtrid<<"  trid "<<hit.trid<<" pz "<<hit.pz<<" test11 "<<test11<<" test21 "<<test21<<" test12 "<<test12<<" test22 "<<test22<<std::endl;  
   }
  *test1 = test11&&test21;
  *test2 = test12&&test22;
  *test3 = test13&&test23;
  *test4 = test14&&test24;
  //std::cout<<"test1 "<<*test1<<" test2 "<<*test2<<" test3 "<<*test3<<" test4 "<<*test4<<std::endl;

  return;
}

void isValid1(std::vector<remollGenericDetectorHit_t> *fHit, std::vector<int> &det28trid, std::vector<int> &Ring5trid, std::vector<int> &Boretrid, std::vector<int> &Acceptancetrid)
{
  for(size_t i=0;i<fHit->size();i++){
    remollGenericDetectorHit_t hit=fHit->at(i); 
    int det = hit.det; 
    int pid = hit.pid;
    if(det==28 && hit.r>=640 &&hit.r<=1200 && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	det28trid.push_back(hit.trid);
    }
    if(det==28 && hit.r>=935 &&hit.r<=1075 && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	Ring5trid.push_back(hit.trid);
    }
    if(det==63 && hit.r<=31 && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	Boretrid.push_back(hit.trid);
    }
    if(det==63 && hit.r<=200 && hit.vz<=-3875){
	//std:;cout<<"Inside isvalid1 det "<<det<<" r "<<hit.r<<"   track id "<<hit.trid<<std::endl;
	Acceptancetrid.push_back(hit.trid);
    }
 }
 // std::cout<<"end of isvalid"<<std::endl;
}
