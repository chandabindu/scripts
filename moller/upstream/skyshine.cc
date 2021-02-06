#include <vector>
#include <string>
#include <sstream>
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
#include <TH2D.h>
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
#include <TLatex.h>
#include "remolltypes.hh"
#define MAXHIT 100000
//#pragma cling load("libremoll.so");
void set_plot_style();
const int nenergies = 3;
const int nparticles = 3;
Float_t energy_ranges[nparticles][nenergies+1] = {{0,10,100,12100},{0,10,100,12100},{0,10,30,12100}};
Float_t hallx,hallz;
Int_t GetEnergyRange(Int_t pid, Float_t e);
Int_t GetDetector(Float_t y, Float_t py, Float_t theta, Float_t phi);
void skyshine()
{       gROOT->SetStyle("Plain");
	gSystem->Load("libremoll.so");
    	gStyle->SetOptStat(0); 
        //gStyle->SetOptStat("eMR");
        gStyle->SetNumberContours(255);
        
	const Int_t nofile = 1000;
	TString Sparticles[nparticles] = {"Electron/Positron","Photon","Neutron"};
	TString Senergies[nparticles][nenergies] = {{"KE<10","10<KE<100","100<KE"}, {"KE<10","10<KE<100","100<KE"}, {"KE<10","10<KE<30","30<KE"}};

	Float_t TotEn[nparticles][nenergies] = {{0}};
	Int_t TotCount[nparticles][nenergies] = {{0}};
        Int_t CESumed[nparticles] = {0};//C: count, E:Energy,
        Int_t CDetSumed[nparticles] = {0};

	Float_t EMEnergy= 0.;

	Int_t HE_N = 0;//number of high energy neutron
	Float_t HEN_Energy= 0.;//energy of high energy neutron

        Float_t EESumed[nparticles]={0};
	Float_t EDEP = 0;

	
	//Histograms
	TH1F *h_kinE[nparticles][nenergies];
	TH1F *h_TotkinE[nparticles][nenergies];//this is for all the region for a particular detector
	TH2D *h_xz[nparticles][nenergies];
	TH2D *h_vertex[nparticles][nenergies];
	//TH2D *h_base[det][nregions+1][nparticles][nenergies];
	//TH2D *h_vertex = new TH2D("h_vertex","vertex",66,-33,33,66,-33,33);

	Int_t theenergy = -10;
	Int_t theparticle = -10;
	Int_t particle;
	


        Float_t energy_bin[nenergies] = {100,90,1200};

        TCanvas *plot;
        TCanvas *plot1;
        TCanvas *plot2;
        TCanvas *plot3;
        TCanvas *plot5;
        TCanvas *plot6;

	for(int ip=0;ip<nparticles;ip++)
	{
		for(int iq=0;iq<nenergies;iq++)
		{
			h_kinE[ip][iq] = new TH1F(Form("KinE_%d_%d",ip+1,iq+1),Form("%s from in %s range: kinE (MeV)",Sparticles[ip].Data(),Senergies[ip][iq].Data()),energy_bin[iq],energy_ranges[ip][iq],energy_ranges[ip][iq+1]);
			h_xz[ip][iq] = new TH2D(Form("h_xz_%d_%d",ip+1,iq+1),Form("XZ: for %s for %s",Sparticles[ip].Data(),Senergies[ip][iq].Data()),130,-26,39,132,-33,33);
			h_vertex[ip][iq] = new TH2D(Form("h_vertex_%d_%d",ip+1,iq+1),Form("XZ vertices: for %s for %s",Sparticles[ip].Data(),Senergies[ip][iq].Data()),130,-26,39,132,-33,33);
		}
	}
	
/*	
	for(int ip=0;ip<nparticles;ip++)
	{
		h_TotkinE[ip] = new TH1F(Form("TotKinE_%d",ip+1),Form("Tot kinE for %s )",Sparticles[ip].Data()),energy_bin[iq],energy_ranges[ip][iq],energy_ranges[ip][iq+1]);
	}
*/		

	TChain *T = new TChain("T");
    	TString added_file_array[nofile]={""};
//    	TString outrootfile;
//    	TString outtxtfile;
    	for(int v =1; v<=nofile;v++)
	{
		ostringstream temp_string1;
		ostringstream temp_string2;
		temp_string1<<v;
		TString vS;
		vS = temp_string1.str();
		temp_string2<<"remollout"<<Form("%d",v)<<"_101.root";
		added_file_array[v-1]=temp_string2.str();
		cout<<temp_string2.str()<<endl;
		//cout<<"file "<<v<<endl;
	        T->Add(added_file_array[v-1]);
		//cout<<"file next "<<v<<endl;
	 }
        
     	TString outrootfile = "moller_targetwithconcrete.root";
    	TString outtxtfile = "moller_targetwithconcrete.txt";
    	ofstream fh_out;
    	fh_out.open(outtxtfile,std::ofstream::out);
    	TFile *rootfile = new TFile(outrootfile,"RECREATE");
    	rootfile->cd();

	//T->Add("remollout1.root");
   	std::vector< remollGenericDetectorHit_t > *fGenDetHit =0;
    	Float_t volume,px,py,pz,p,mass;
    	Float_t energy, pdgId, kinE,edep;
    	Float_t vx, vy, vz, vr;
   	Float_t vz_min, vz_max, vy_max;
	Float_t hally;
	Float_t hitx,hity,hitz,hitr;
	Float_t distance,dummy,dummy1;
	Float_t theta,phi;
	Float_t rad2deg = 180./TMath::Pi();
     
        T->SetBranchAddress("hit",&fGenDetHit);
    	Int_t entries = T->GetEntries();
    	cout<<"Entries  "<<entries<<endl;
        set_plot_style();
        for(int ij =0;ij<entries;ij++)
        //for(int ij =0;ij<10000;ij++)
    	{
	    if(ij%500000==0)cout<<ij<<"  Events analysed"<<endl;
	    T->GetEntry(ij);
	    for(size_t pk = 0; pk<fGenDetHit->size(); pk++)
	    {		    
		    volume = fGenDetHit->at(pk).det;
                    
		    if(volume==101){
		    px = fGenDetHit->at(pk).px;
                    py = fGenDetHit->at(pk).py;
                    pz = fGenDetHit->at(pk).pz;
                    energy = fGenDetHit->at(pk).e;
		    edep = fGenDetHit->at(pk).edep;
                    pdgId = fGenDetHit->at(pk).pid;
                    mass = fGenDetHit->at(pk).m;
                    vz = fGenDetHit->at(pk).vz;
                    vx = fGenDetHit->at(pk).vx;
                    vy = fGenDetHit->at(pk).vy;
                    hitx = fGenDetHit->at(pk).x;
                    hity = fGenDetHit->at(pk).y;
                    hitz = fGenDetHit->at(pk).z;
                    kinE = energy - mass;
		    particle = (Int_t)TMath::Abs(pdgId);// pidmap[(Int_t)TMath::Abs(fGenDetHit->at(pk).pid)];
		    //h_vertex->Fill(vz/1000.,vx/1000.);
		   
		    hitr = sqrt(hitx*hitx+hity*hity+hitz*hitz);
		    //cout<<"hitr "<<hitr<<endl;
		    theta = acos((hitz)/hitr);//*rad2deg;
		    phi = TMath::ATan2(hity,hitx);//*rad2deg;
		    p=sqrt(px*px+py*py+pz*pz);
		    
		    if(particle==11) theparticle = 0;
		    else if(particle==22) theparticle = 1;
		    else if(particle==2112) theparticle=2;
		    else continue;
		    
		    theenergy = GetEnergyRange(theparticle,kinE);

		    if(theparticle>=0 && theenergy>=0)
		    {
			    h_kinE[theparticle][theenergy]->Fill(kinE);
			   // h_TotkinE[theparticle]->Fill(kinE);
			    h_xz[theparticle][theenergy]->Fill(hitz/1000.,hitx/1000.);
			    h_vertex[theparticle][theenergy]->Fill(vz/1000.,vx/1000.);
		    }
		    
		    if((theparticle==0 || theparticle==1) && kinE>=100)
		    {
			    EMEnergy += kinE;
		    }
		    if(theparticle==2 && kinE>=30)
		    {
			HE_N += 1;
			HEN_Energy +=kinE;
		    } 
		    

		    TotEn[theparticle][theenergy] += kinE;
		}
}
	}
	
	
	plot=new TCanvas("C_KinE","C_kineE",800,800);
	plot->Divide(nparticles,nenergies);
	for(int ip=0;ip<nparticles;ip++)
	{
		for(int iq=0;iq<nenergies;iq++)
		{
			TotCount[ip][iq] = h_kinE[ip][iq]->Integral();
			plot->cd(nenergies*ip+iq+1);
			plot->cd(nenergies*ip+iq+1)->SetLogy();
			h_kinE[ip][iq]->GetXaxis()->SetTitle("Energy (MeV)");
			h_kinE[ip][iq]->GetYaxis()->SetTitle("Counts");
			h_kinE[ip][iq]->GetXaxis()->SetTitleSize(0.05);
			h_kinE[ip][iq]->GetYaxis()->SetTitleSize(0.05);
			h_kinE[ip][iq]->GetXaxis()->SetLabelSize(0.05);
			h_kinE[ip][iq]->GetYaxis()->SetLabelSize(0.05);

			h_kinE[ip][iq]->Draw();
			//cout<<"Region "<<ik<<" particle  "<<ip<<" energy range "<<iq<<" Integral "<<TotCount[ij][ik][ip][iq]<<endl;
		}

	}


	plot1=new TCanvas("C_XZ","C_XZ_distribution",800,800);
		
	plot1->Divide(nparticles,nenergies);
	for(int ip=0;ip<nparticles;ip++)
	{
		for(int iq=0;iq<nenergies;iq++)
		{
			plot1->cd(nenergies*ip+iq+1);
			plot1->cd(nenergies*ip+iq+1)->SetLogz();
			h_xz[ip][iq]->GetXaxis()->SetTitle("Z (m)");
                        h_xz[ip][iq]->GetYaxis()->SetTitle("X (m)");
                        h_xz[ip][iq]->GetXaxis()->SetTitleSize(0.05);
                        h_xz[ip][iq]->GetYaxis()->SetTitleSize(0.05);
                        h_xz[ip][iq]->GetXaxis()->SetLabelSize(0.05);
                        h_xz[ip][iq]->GetYaxis()->SetLabelSize(0.05);
			h_xz[ip][iq]->Draw("COLZ");
		}
	}
	/*	
	plot3=new TCanvas("Tot KinE distribution","Tot kinE distribution",800,800);

        plot3->Divide(nparticles);
        for(int ip=0;ip<nparticles;ip++)
        {
		plot3->cd(nenergies*ip+1);
                plot3->cd(nenergies*ip+1)->SetLogy();
                h_TotkinE[ip]->GetXaxis()->SetTitle("Energy (MeV)");
                h_TotkinE[ip]->GetYaxis()->SetTitle("Counts");
                h_TotkinE[ip]->GetXaxis()->SetTitleSize(0.05);
                h_TotkinE[ip]->GetYaxis()->SetTitleSize(0.05);
                h_TotkinE[ip]->Draw();
         }
         */
                
        fh_out<<"For 100 million particle on target"<<endl;
        fh_out<<"Particle:: 0 - electron, 1- photon, 2 -neutron"<<endl;
	for(int ip=0;ip<nparticles;ip++)
	{
		for(int iq =0;iq<nenergies;iq++)
		{
			CESumed[ip] += TotCount[ip][iq];
			EESumed[ip] += TotEn[ip][iq];
					
		}
		fh_out<<"Summed over energy ranges for particle "<<ip<<" counts "<<CESumed[ip]<<"  energy  "<<EESumed[ip]<<" MeV"<<endl; 
	}
	
	fh_out<<"Number of high energy (>30 MeV) neutron going towards :roof "<<HE_N<<"  Energy  "<<HEN_Energy<<" MeV"<<endl;
	
	/*
	plot3->Write();
	for(int ip=0;ip<nparticles;ip++)
	{
		h_TotkinE[ip]->Write();
	}
	*/

	plot6=new TCanvas("C_XZ_vertices","vertices_ XZ",800,800);
		
	plot6->Divide(nparticles,nenergies);
	for(int ip=0;ip<nparticles;ip++)
	{
		for(int iq=0;iq<nenergies;iq++)
		{
			plot6->cd(nenergies*ip+iq+1);
			plot6->cd(nenergies*ip+iq+1)->SetLogz();
			h_vertex[ip][iq]->GetXaxis()->SetTitle("Z (m)");
                        h_vertex[ip][iq]->GetYaxis()->SetTitle("X (m)");
                        h_vertex[ip][iq]->GetXaxis()->SetTitleSize(0.05);
                        h_vertex[ip][iq]->GetXaxis()->SetLabelSize(0.05);
                        h_vertex[ip][iq]->GetYaxis()->SetTitleSize(0.05);
                        h_vertex[ip][iq]->GetYaxis()->SetLabelSize(0.05);
			h_vertex[ip][iq]->Draw("COLZ");
		}
	}
/*
	TCanvas *plot6 = new TCanvas("Vertex","Vertex",800,800);
	plot6->SetLogz();
        h_vertex->GetXaxis()->SetTitle("Z (m)");
        h_vertex->GetYaxis()->SetTitle("X (m)");
	h_vertex->GetXaxis()->SetTitleSize(0.04);
        h_vertex->GetYaxis()->SetTitleSize(0.04);
        h_vertex->GetYaxis()->SetTitleOffset(1.20);
	h_vertex->Draw("COLZ");
	

	h_vertex->Write();
	plot6->Write();
*/
	for(int ip=0;ip<nparticles;ip++)
	{
		for(int iq=0;iq<nenergies;iq++)
		{
			h_kinE[ip][iq]->Write();
			h_xz[ip][iq]->Write();
			h_vertex[ip][iq]->Write();
		}
	}
	
	plot->Write();
	plot1->Write();
	plot6->Write();
	rootfile->Write();
	rootfile->Close();
	fh_out.close();

				
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

Int_t GetEnergyRange(Int_t pid, Float_t e)
{
	Int_t theregion =-10;
	for(int ij=0;ij<nenergies;ij++)
	{
		if(e>=energy_ranges[pid][ij] &&  e<energy_ranges[pid][ij+1])
		{
		theregion = ij;
		break;
		//return ij;
		}
		else
			theregion=  -10;
	}

	return theregion;
}

