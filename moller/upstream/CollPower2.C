/*
This program for calculating the power on the collimator 1,2,4,5 for MOLLER 
experiment New Design
Chandan Ghosh
*/
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <new>
#include <cstdlib>
#include <math.h>

#include <TRandom.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TSystem.h>

#include <TH2F.h>
#include <TTree.h>
#include <TF1.h>
#include <TProfile.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TDatime.h>
#include <TStopwatch.h>
#include <stdexcept>
#include <time.h>
#include <cstdio>
#include <map>
#include <cassert>

#include <TMath.h>
#include <TStyle.h>
#include <TPaveStats.h>

#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TFrame.h>
#include <TObjArray.h>
#include <TVector2.h>
#include "remolltypes.hh"

using namespace std;

void set_plot_style();

void CollPower2()
{
	  gROOT->SetStyle("Plain");
	  gStyle->SetOptStat(0); 
	  //gSystem->Load("libremoll.so");
	  //gStyle->SetOptStat("eMR");
	  set_plot_style();

	
	  const Int_t nofile =1;
	  const Int_t ndet =10;
	  Int_t collimator[ndet]={2001,2002,2008,2004,2005,2009,2010,2011,2012,2007};
	  TString sDet[ndet]={"Col1","Col2_CW","Col2_Cu","Col4","Col5","Col1_H2O_US","Col1_H2O_DS","Col1_H2O_CW","Col1_Jacket","lintel"};
	  Int_t volume;
	  Double_t hitx,hity,hitz;
	  Double_t dumyenergy=0.;
	  Int_t particle;
	  Double_t esum[ndet]={0.};

	  //Histograms
	  TH2F* h_2d[ndet];
	  TH2F* h_2d_rz[ndet];
	  h_2d[0] = new TH2F("h_col1_CW","Coll 1 CW  power",400,-4,4,400,-4,4);//
	  h_2d[1] = new TH2F("h_col2_CW","Coll 2 CW Power",700,-35,35,700,-35,35);
	  h_2d[2] = new TH2F("h_col2_Cu","Coll 2 Cu power",700,-35,35,700,-35,35);
	  h_2d[3] = new TH2F("h_col4","Coll 4 power",700,-35,35,700,-35,35);
	  h_2d[4] = new TH2F("h_col5","Coll 5 power",400,-20,20,400,-20,20);
	  h_2d[5] = new TH2F("h_col1_H2O_US","Col1 H2O US Power",400,-4,4,400,-4,4);
	  h_2d[6] = new TH2F("h_col1_H2O_DS","Col1 H2O DS Power",400,-4,4,400,-4,4);
	  h_2d[7] = new TH2F("h_col1_H2O_CW","Col1 H2O+CW Power",400,-4,4,400,-4,4);
	  h_2d[8] = new TH2F("h_col1_jacket","Col1 Jacket Power",400,-4,4,400,-4,4);
	  h_2d[9] = new TH2F("h_lintel","lintel Power",140,-70,70,140,-70,70);
	  TString sXY[ndet] = {"25W/mm^{2}","W/mm^{2}","W/mm^{2}","W/mm^{2}","W/mm^{2}","25W/mm^{2}","25W/mm^{2}","25W/mm^{2}","25W/mm^{2}","W/cm^{2}"};

	  h_2d_rz[0] = new TH2F("Col1_CW_rz","Coll_1_rz",600,0.35,0.95,300,1.0,4.0);//z in m; r in cm //0.1mm^2
	  h_2d_rz[1] = new TH2F("Col2_CW_rz","Coll2_CW_rz",250,0.700,0.95,300,2.0,32.0);//z in m; r in cm //1 mm^2
	  h_2d_rz[2] = new TH2F("Col2_Cu_rz","Coll2_Cu_rz",250,0.700,0.95,300,2.0,32.0);//z in m; r in cm //mm^2
	  h_2d_rz[3] = new TH2F("Col_4_rz","Coll_4_rz",200,3.20,3.4,300,2,32);//z in m; r in cm //mm^2
	  h_2d_rz[4] = new TH2F("Col_5_rz","Coll_5_rz",110,7.78,7.89,140,6,20);//z in m; r in cm //mm^2
	  h_2d_rz[5] = new TH2F("Col1_H2O_US_rz","Col1_H20_US_rz",600,0.35,0.95,300,1.0,4.0);//z in m; r in cm //0.1mm^2
	  h_2d_rz[6] = new TH2F("Col1_H2O_DS_rz","Col1_H20_DS_rz",600,0.35,0.95,300,1.0,4.0);//z in m; r in cm //0.1mm^2
	  h_2d_rz[7] = new TH2F("Col1_H2O_CW_rz","Col1_H20_CW_rz",600,0.35,0.95,300,1.0,4.0);//z in m; r in cm //0.1mm^2
	  h_2d_rz[8] = new TH2F("Col1_jacket_rz","Col1_jacket_rz",600,0.35,0.95,300,1.0,4.0);//z in m; r in cm //0.1mm^2
	  h_2d_rz[9] = new TH2F("lintel_rz","lintel_rz",200,7.735,7.935,300,40.0,70.0);//z in m; r in cm //mm^2
	  TString sRZ[ndet] = {"10W/mm^{2}","W/mm^{2}","W/mm^{2}","W/mm^{2}","W/mm^{2}","10W/mm^{2}","10W/mm^{2}","10W/mm^{2}","10W/mm^{2}","W/mm^{2}"};
	  //TH1F* KE_lintel = new TH1F("lintel_ke","lintel_ke",1000,0,10000);
	  //collimator edep
	  TChain * Tmol =new TChain("T");
	  TString added_file_array[nofile]={""};
	  for(int v=1;v<=nofile;v++)
	  {
		ostringstream temp_string1;
                ostringstream temp_string2;
                temp_string1<<v;
                TString vS;
		int dumyint = 2;
                vS = temp_string1.str();
                temp_string2<<Form("remollout%d",dumyint)<<".root";
                added_file_array[v-1]=temp_string2.str();
                cout<<temp_string2.str()<<endl;
                //cout<<"file "<<v<<endl;
                //Tmol->Add(added_file_array[v-1]);
                //cout<<"file next "<<v<<endl;
         }
	  //Tmol->Add("/home/lead/G4WORK/remoll/remollout_v10.root");
	  Tmol->Add("/w/halla-scifs17exp/parity/disk1/chandan/gitdir/remoll/g10_6_100k_merged.root");
	  std::vector< remollGenericDetectorHit_t > *fGenDetHit =0;
	  Tmol->SetBranchAddress("hit",&fGenDetHit);
	  
	  Int_t nentries = (Int_t)Tmol->GetEntries();
	  Double_t MeV2Watt = 70./(1.*nentries);//for 70 uA
	  std::cout<<"Number of entries in TChain "<<nentries<<std::endl;
	  //for(int ij=0;ij<10000;ij++)//Event loop
	  for(int ij=0;ij<nentries;ij++)//Event loop
	  {
		  Tmol->GetEntry(ij);
		  for(size_t pk=0;pk<fGenDetHit->size();pk++)
		  {
			  volume = fGenDetHit->at(pk).det;
			  dumyenergy = fGenDetHit->at(pk).edep;
			  hitx =fGenDetHit->at(pk).x/10;
			  hity =fGenDetHit->at(pk).y/10;
			  hitz =fGenDetHit->at(pk).z/1000;
			  //std::cout<<"Det "<<volume<<" hitx "<<hitx<<" hity  "<<hity<<"  hitz  "<<hitz<<"  edep  "<<dumyenergy<<std::endl; 
			  if(volume==collimator[0])//collimator 1
			  {
				  esum[0] +=dumyenergy;
			  	  h_2d_rz[0]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[0]->Fill(hitx,hity,MeV2Watt*dumyenergy);
				  
			  }
			  if(volume==collimator[1])//collimator 2 CW 
			  {
				  esum[1] +=dumyenergy;
			  	  h_2d_rz[1]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[1]->Fill(hitx,hity,MeV2Watt*dumyenergy);
				  
			  }
			  if(volume==collimator[2])//collimator 2 Cu
			  {
				  esum[2] +=dumyenergy;
			  	  h_2d_rz[2]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[2]->Fill(hitx,hity,MeV2Watt*dumyenergy);
			  }
			  if(volume==collimator[3])//collimator 4
			  {
				  esum[3] +=dumyenergy;
			  	  h_2d_rz[3]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[3]->Fill(hitx,hity,MeV2Watt*dumyenergy);
			  }
			  if(volume==collimator[4])//collimator 5
			  {
				  esum[4] +=dumyenergy;
			  	  h_2d_rz[4]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[4]->Fill(hitx,hity,MeV2Watt*dumyenergy);
				  
			  }
			  if(volume==collimator[5])//collimator 1 US water part
                          {
                                  esum[5] +=dumyenergy;
                                  h_2d_rz[5]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
                                  h_2d[5]->Fill(hitx,hity,MeV2Watt*dumyenergy);
                          }
			  if(volume==collimator[6])//collimator 1 DS water part
                          {
                                  esum[6] +=dumyenergy;
                                  h_2d_rz[6]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
                                  h_2d[6]->Fill(hitx,hity,MeV2Watt*dumyenergy);
                          }
			  if(volume==collimator[7])//collimator 1 water+CW part
			  {
				  esum[7] +=dumyenergy;
			  	  h_2d_rz[7]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[7]->Fill(hitx,hity,MeV2Watt*dumyenergy);
			  }
			  if(volume==collimator[8])//collimator 1 jacket
			  {
				  esum[8] +=dumyenergy;
			  	  h_2d_rz[8]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[8]->Fill(hitx,hity,MeV2Watt*dumyenergy);
			  }
			  if(volume==collimator[9])//collimator 1 jacket
			  {
				  esum[9] +=dumyenergy;
			  	  h_2d_rz[9]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[9]->Fill(hitx,hity,MeV2Watt*dumyenergy);
			  }
		  }
	  }
	  for(int pq=0;pq<ndet;pq++)
	  {
		  std::cout<<Form("Power deposition inside %s  :",sDet[pq].Data())<<esum[pq]*70./nentries<<"W/70 uA"<<std::endl;
	  }
	  std::cout<<"Total energy deposited inside collimator 1  "<<(esum[0]+esum[5]+esum[6]+esum[7]+esum[8])*70./nentries<<"W/70 uA"<<std::endl;//sum of collimator 1 and water and jacket
	  std::cout<<"Total energy deposited inside collimator 2  "<<(esum[1]+esum[2])*70./nentries<<"W/70 uA"<<std::endl;//sum of collimator 1 and water and jacket

	  TCanvas *plot[ndet];
	  for(int ij=0;ij<ndet;ij++)
	  {
		  plot[ij] = new TCanvas(Form("C_%s",sDet[ij].Data()),Form("C_%s",sDet[ij].Data()),600,600);

		  plot[ij]->SetRightMargin(0.22);
		  plot[ij]->SetLeftMargin(0.125);
		  h_2d[ij]->GetXaxis()->SetTitle("X (cm)");
		  h_2d[ij]->GetXaxis()->SetTitleSize(0.05);
		  h_2d[ij]->GetXaxis()->SetLabelSize(0.05);
		  h_2d[ij]->GetYaxis()->SetTitle("Y (cm)");
		  h_2d[ij]->GetYaxis()->SetTitleSize(0.05);
		  h_2d[ij]->GetYaxis()->SetLabelSize(0.05);
		  //h_2d[ij]->GetZaxis()->SetTitle("Power (arb unit)");
		  h_2d[ij]->GetZaxis()->SetTitle(Form("Power (%s/70#muA)",sXY[ij].Data()));
		  h_2d[ij]->GetZaxis()->SetLabelSize(0.05);
		  h_2d[ij]->GetZaxis()->SetTitleOffset(2.0);

		  h_2d[ij]->Draw("COLZ");
		  
		  plot[ij]->SaveAs(Form("%s_xy.png",sDet[ij].Data()));
	  }

	  TCanvas *plot1[ndet];
	  for(int ij=0;ij<ndet;ij++)
	  {
			plot1[ij]= new TCanvas(Form("C_%s_rz",sDet[ij].Data()),Form("C_%s_rz",sDet[ij].Data()),600,600);
		 plot1[ij]->SetRightMargin(0.22);
		 plot1[ij]->SetLeftMargin(0.125);
		 h_2d_rz[ij]->GetXaxis()->SetTitle("Z (m)");
		 h_2d_rz[ij]->GetXaxis()->SetTitleSize(0.05);
		 h_2d_rz[ij]->GetXaxis()->SetNdivisions(505);
		 h_2d_rz[ij]->GetXaxis()->SetLabelSize(0.05);
		 h_2d_rz[ij]->GetYaxis()->SetTitle("r (cm)");
		 h_2d_rz[ij]->GetYaxis()->SetTitleSize(0.05);
		 h_2d_rz[ij]->GetYaxis()->SetLabelSize(0.05);
		 //h_2d_rz[ij]->GetZaxis()->SetTitle("Power (W/(mm^2) 70#muA)");
		 h_2d_rz[ij]->GetZaxis()->SetTitle(Form("Power (%s/70#muA)",sRZ[ij].Data()));
		 h_2d_rz[ij]->GetZaxis()->SetLabelSize(0.05);
		 h_2d_rz[ij]->GetZaxis()->SetNdivisions(505);
		 h_2d_rz[ij]->GetZaxis()->SetTitleOffset(2.0);
	 	 h_2d_rz[ij]->Draw("COLZ");
		 plot1[ij]->SaveAs(Form("%s_rz.png",sDet[ij].Data()));
	  }


}
void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    // See class TColor documentation and SetPalette() command
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

