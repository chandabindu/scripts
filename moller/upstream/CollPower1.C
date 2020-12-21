/*
This program for calculating the power on the collimator 1,2,4,5 for MOLLER 
experiment old design
Chandan Ghosh
part of this is taken from Rakitha's code coll_scan.cc
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

void CollPower1()
{
	  gROOT->SetStyle("Plain");
	  gStyle->SetOptStat(0); 
	  //gSystem->Load("libremoll.so");
	  //gStyle->SetOptStat("eMR");
	  set_plot_style();

	
	  const Int_t nofile =1;
	  Int_t collimator[6]={2001,2006,2002,2004,2005,2007};
	  Int_t volume;
	  Double_t hitx,hity,hitz;
	  Double_t dumyenergy=0.;
	  Int_t particle;
	  Double_t esum[6]={0.};

	  //Histograms
	  TH2F* h_2d[6];
	  TH2F* h_2d_rz[6];
	  h_2d[0] = new TH2F("h_coll_1","Coll 1 power",400,-4,4,400,-4,4);
	  h_2d[1] = new TH2F("h_coll_1_fins","Coll 1 fins",650,-6.5,6.5,600,-6,6);
	  h_2d[2] = new TH2F("h_coll_2","Coll 2 power",700,-35,35,700,-35,35);
	  h_2d[3] = new TH2F("h_coll_4","Coll 4 power",700,-35,35,700,-35,35);
	  h_2d[4] = new TH2F("h_coll_5","Coll 5 power",400,-20,20,400,-20,20);
	  h_2d[5] = new TH2F("h_coll_6","Lintel Power",360,-90,90,360,-90,90);

	  h_2d_rz[0] = new TH2F("Coll_1_rz","Coll_1_rz",600,0.1,0.7,300,1.0,4.0);//z in m; r in cm //0.1mm^2
	  h_2d_rz[1] = new TH2F("Coll_1_fins_rz","Coll_1_fins_rz",170,5.150,5.320,400,2.5,6.5);//z in m; r in cm //0.1 mm^2
	  h_2d_rz[2] = new TH2F("Coll_2_rz","Coll_2_rz",250,0.70,0.950,300,2,32);//z in m; r in cm //mm^2
	  h_2d_rz[3] = new TH2F("Coll_4_rz","Coll_4_rz",200,3.20,3.400,300,2,32);//z in m; r in cm //mm^2
	  h_2d_rz[4] = new TH2F("Coll_5_rz","Coll_5_rz",110,7.780,7.890,140,6,20);//z in m; r in cm //mm^2
	  h_2d_rz[5] = new TH2F("Coll_6_rz","lintel_rz",140,12.780,12.920,600,30,90);//z in m; r in cm //mm^2
	  TH1F* KE_lintel = new TH1F("lintel_ke","lintel_ke",1000,0,10000);
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
	  //Tmol->Add("/home/lead/G4WORK/remoll/remollout1.root");
	  Tmol->Add("/w/halla-scifs17exp/parity/disk1/chandan/gitdir/remoll/g10_6_100k.root");
	  std::vector< remollGenericDetectorHit_t > *fGenDetHit =0;
	  Tmol->SetBranchAddress("hit",&fGenDetHit);
	  
	  Int_t nentries = (Int_t)Tmol->GetEntries();
	  Double_t MeV2Watt = 65./(1.*nentries);//for 65 uA
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
			  
			  if(volume==collimator[0])//collimator 1
			  {
				  esum[0] +=dumyenergy;
			  	  h_2d_rz[0]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[0]->Fill(hitx,hity,MeV2Watt*dumyenergy);
				  
			  }
			  if(volume==collimator[1])//collimator 1 fins
			  {
				  esum[1] +=dumyenergy;
			  	  h_2d_rz[1]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
			  	  h_2d[1]->Fill(hitx,hity,MeV2Watt*dumyenergy);
				  
			  }
			  if(volume==collimator[2])//collimator 2
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
			  if(volume==collimator[5])//collimator 5
                          {
                                  esum[5] +=dumyenergy;
                                  h_2d_rz[5]->Fill(hitz,sqrt(hitx*hitx+hity*hity),MeV2Watt*dumyenergy);
                                  h_2d[5]->Fill(hitx,hity,MeV2Watt*dumyenergy);
     				  KE_lintel->Fill(dumyenergy);
                          }


		  }
	  }
	  for(int pq=0;pq<6;pq++)
	  {
		  std::cout<<"Collimator  "<<collimator[pq]<<"  Deposited energy "<<esum[pq]*65./nentries<<"W/65 uA"<<std::endl;
	  }
	  std::cout<<"Total energy deposited inside collimator "<<esum[0]*65./nentries+esum[1]*65./nentries<<"W/65 uA"<<std::endl;//sum of collimator 1 and fins

	  TCanvas *plot[6];
	  for(int ij=0;ij<6;ij++)
	  {
		  if(ij==0 || ij==3 || ij==4)
		  plot[ij] = new TCanvas(Form("Collimator_%d",ij+1),Form("Collimator_%d",ij+1),600,600);
		  if(ij==1)
		  plot[ij] = new TCanvas(Form("Collimator_%d_fins",ij),Form("Collimator_%d_fins",ij),600,600);
		  if(ij==2)
		  plot[ij] = new TCanvas(Form("Collimator_%d",ij),Form("Collimator_%d",ij),600,600);
		   if(ij==5)
                  plot[ij] = new TCanvas(Form("lintel_%d",ij),Form("lintel_%d",ij),600,600);

		  plot[ij]->SetRightMargin(0.22);
		  plot[ij]->SetLeftMargin(0.125);
		  h_2d[ij]->GetXaxis()->SetTitle("X (cm)");
		  h_2d[ij]->GetXaxis()->SetTitleSize(0.05);
		  h_2d[ij]->GetXaxis()->SetLabelSize(0.05);
		  h_2d[ij]->GetYaxis()->SetTitle("Y (cm)");
		  h_2d[ij]->GetYaxis()->SetTitleSize(0.05);
		  h_2d[ij]->GetYaxis()->SetLabelSize(0.05);
		  h_2d[ij]->GetZaxis()->SetTitle("Power (arb unit)");
		  h_2d[ij]->GetZaxis()->SetLabelSize(0.05);
		  h_2d[ij]->GetZaxis()->SetTitleOffset(2.0);

		  h_2d[ij]->Draw("COLZ");
		  
		  if(ij==0 || ij==3 || ij==4)
		  plot[ij]->SaveAs(Form("Collimator_%d_xy.png",ij+1));
		  if(ij==1)
		  plot[ij]->SaveAs(Form("Collimator_%d_fins_xy.png",ij));
		  if(ij==2)
		  plot[ij]->SaveAs(Form("Collimator_%d_xy.png",ij));
		  if(ij==5)
		  plot[ij]->SaveAs(Form("lintel_%d_xy.png",ij));
	  }

	  TCanvas *plot1[6];
	  for(int ij=0;ij<6;ij++)
	  {
		 if(ij==0 ||ij==3|ij==4)
			plot1[ij]= new TCanvas(Form("Coll_%d_rz",ij+1),Form("Collimator %d rz power distribution",ij+1),600,600);
		 if(ij==1)
			plot1[ij]= new TCanvas(Form("Coll_%d_fins_rz",ij),Form("Collimator %d fins rz power distribution",ij),600,600);
		 if(ij==2)
			plot1[ij]= new TCanvas(Form("Coll_%d_rz",ij),Form("Collimator %d rz power distribution",ij),600,600);
		 if(ij==5)
			plot1[ij]= new TCanvas(Form("lintel_%d_rz",ij),Form("lintel %d rz power distribution",ij),600,600);
		 plot1[ij]->SetRightMargin(0.22);
		 plot1[ij]->SetLeftMargin(0.125);
		 h_2d_rz[ij]->GetXaxis()->SetTitle("Z (m)");
		 h_2d_rz[ij]->GetXaxis()->SetTitleSize(0.05);
		 h_2d_rz[ij]->GetXaxis()->SetNdivisions(505);
		 h_2d_rz[ij]->GetXaxis()->SetLabelSize(0.05);
		 h_2d_rz[ij]->GetYaxis()->SetTitle("r (cm)");
		 h_2d_rz[ij]->GetYaxis()->SetTitleSize(0.05);
		 h_2d_rz[ij]->GetYaxis()->SetLabelSize(0.05);
		 h_2d_rz[ij]->GetZaxis()->SetTitle("Power (W/(mm^2) 65#muA)");
		 if(ij<2)
		 h_2d_rz[ij]->GetZaxis()->SetTitle("Power (W/(0.1mm^2) 65#muA)");
		 h_2d_rz[ij]->GetZaxis()->SetLabelSize(0.05);
		 h_2d_rz[ij]->GetZaxis()->SetNdivisions(505);
		 h_2d_rz[ij]->GetZaxis()->SetTitleOffset(2.0);
	 	 h_2d_rz[ij]->Draw("COLZ");
		 if(ij==0 || ij==3 || ij==4)
		 plot1[ij]->SaveAs(Form("Collimator_%d_rz.png",ij+1));
		 if(ij==1)
		 plot1[ij]->SaveAs(Form("Collimator_%d_fins_rz.png",ij));
		 if(ij==2)
		 plot1[ij]->SaveAs(Form("Collimator_%d_rz.png",ij));
		 if(ij==5)
		 plot1[ij]->SaveAs(Form("lintel_%d_rz.png",ij));
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

