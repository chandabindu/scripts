#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
void set_plot_style();
const int nparticle = 3;
void beamAna1()
{
	gROOT->SetStyle("Plain");
	//gSystem->Load("libremoll.so");
	gStyle->SetOptStat(1000011);
	TString fin;
	TString fout;
	TString Sparticle[nparticle]={"Electron/Positron","Photon","Neutron"};
	Float_t current = 70.;//uA
	Float_t factor = 6.25e12;
	Double_t rate=0;
	ostringstream temp_string1;
	ostringstream temp_string2;
	ostringstream temp_string3;
	ostringstream temp_string4;
	temp_string1<<"beam_bare";
	//temp_string1<<"beam_default";
	temp_string3<<"remoll_"<<temp_string1.str()<<"_out.root";
	fout = temp_string3.str();
	TFile *fh_out = new TFile(fout,"RECREATE");
	TChain *T = new TChain("T");
	for(int nfile=0;nfile<250;nfile++)
	{
		temp_string2<<"/lustre/expphy/volatile/halla/parity/chandan/sim_out/default2/remollout_beam"<<Form("%d",nfile+1)<<".root";
		fin = temp_string2.str();
		T->Add(fin);
		std::cout<<temp_string2.str()<<std::endl;
		temp_string2.str("");
		temp_string2.clear();
		//std::cout<<temp_string2.str()<<std::endl;
	}
	std::vector< remollGenericDetectorHit_t > *fGenDetHit = 0;
	T->SetBranchAddress("hit", &fGenDetHit);
	//T->SetBranchAddress("rate", &rate);
	Long64_t entries = T->GetEntries();
	std::cout<<"Number of entries  "<<entries<<std::endl;
	
	TH1F* hit_r[nparticle];
	TH2F* h_zphi94[nparticle];
	TH2F* h_zphi95[nparticle];
	TH2F* hVxVz94[nparticle];
	TH2F* hVxVz95[nparticle];
	TH1F *r = new TH1F("r","radial hit",100,0,2100);
	TH1F *rphoton = new TH1F("rphoton","radial hit by photon",100,0,2100);
	TH1F *rRate = new TH1F("rRate","rate weighted radial hit",100,0,2100);
	//TH1F *z = new TH1F("z","z",500,4000,20000);

	TH2F *zphi94 = new TH2F("zphi94","zphi94",200,4000,20000,380,-190,190);
	TH2F *zphi95 = new TH2F("zphi95","zphi95",200,19000,22000,380,-190,190);
	
	for(int ip =0;ip<nparticle;ip++)
	{
		hit_r[ip] = new TH1F(Form("hitr_%d",ip+1),Form("hit r for det 28 %s above 1 MeV",Sparticle[ip].Data()),500,0,2000);	
		h_zphi94[ip] = new TH2F(Form("h_zphi94_%d",ip+1),Form("Phi vz Z  for det 94 for %s above 1 MeV",Sparticle[ip].Data()),200,4000,20000,380,-190,190);	
		h_zphi95[ip] = new TH2F(Form("h_zphi95_%d",ip+1),Form("Phi vz Z  for det 95 for %s above 1 MeV",Sparticle[ip].Data()),200,19000,22000,380,-190,190);	
		hVxVz94[ip] = new TH2F(Form("h_VxVz94_%d",ip+1),Form("Vx vz Vz for det 94 for %s above 1 MeV",Sparticle[ip].Data()),400,-6000,20000,500,-500,500);	
		hVxVz95[ip] = new TH2F(Form("h_VxVz95_%d",ip+1),Form("Vx vz Vz for det 95 for %s above 1 MeV",Sparticle[ip].Data()),400,-6000,22000,500,-500,500);	
	}

	Float_t energy(-1.e-12),vx(-1.e-12),vz(-1.e-12),vy(-1.e-12),detector(1.e-12),hitr,hitz(-1.e-12);
	Int_t pid(0);
        set_plot_style();	
	for(Long64_t ij=0; ij<entries;ij++)
	//for(Long64_t ij=0; ij<1000;ij++)
	{
		if(ij%100000 == 1)
		std::cout<<"Analyzed "<<ij<<" events"<<std::endl;
		T->GetEntry(ij);
		for(size_t pk = 0; pk<fGenDetHit->size();pk++)
		{
			pid = (Int_t)TMath::Abs(fGenDetHit->at(pk).pid);
			detector = fGenDetHit->at(pk).det;
			energy = fGenDetHit->at(pk).e;
			vx = fGenDetHit->at(pk).vx;
			vz = fGenDetHit->at(pk).vz;
			hitr = fGenDetHit->at(pk).r;
			hitz = fGenDetHit->at(pk).z;
			double phi = (180./TMath::Pi())*TMath::ATan2(fGenDetHit->at(pk).y,fGenDetHit->at(pk).x);

			if(detector==28 && hitr>100 )
			{
				r->Fill(hitr);
				rRate->Fill(hitr,factor/entries);
				if(energy>1 && pid==11)
				{
					hit_r[0]->Fill(hitr);
				}
				if(energy>1 && pid==22)
				{
					hit_r[1]->Fill(hitr);
				}
				if(energy>1 && pid==2112)
				{
					hit_r[2]->Fill(hitr);
				}
				//std::cout<<"hit.x "<<fGenDetHit->at(pk).x<<"  hit.y "<<fGenDetHit->at(pk).y<<" radius "<<sqrt(fGenDetHit->at(pk).x*fGenDetHit->at(pk).x+fGenDetHit->at(pk).y*fGenDetHit->at(pk).y)<<" hit.r  "<<fGenDetHit->at(pk).r<<std::endl;
			}
			
			if(detector==28 && pid==22)
			{
				rphoton->Fill(hitr);
			}
			if(detector==94 && hitr>100)
			{
				zphi94->Fill(hitz,phi);
				if(energy>1 && pid==11)
				{
					h_zphi94[0]->Fill(hitz,phi);
					hVxVz94[0]->Fill(vz,vx);	
				}
				if(energy>1 && pid==22)
				{
					h_zphi94[1]->Fill(hitz,phi);
					hVxVz94[1]->Fill(vz,vx);	
				}
				if(energy>1 && pid==2112)
				{
					h_zphi94[2]->Fill(hitz,phi);
					hVxVz94[2]->Fill(vz,vx);	
				}
				//z->Fill(fGenDetHit->at(pk).z);

			}
			if(detector==95 && hitr>100)
			{
				zphi95->Fill(hitz,phi);
				if(energy>1 && pid==11)
				{
					h_zphi95[0]->Fill(hitz,phi);
					hVxVz95[0]->Fill(vz,vx);	
				}
				if(energy>1 && pid==22)
				{
					h_zphi95[1]->Fill(hitz,phi);
					hVxVz95[1]->Fill(vz,vx);	
				}
				if(energy>1 && pid==2112)
				{
					h_zphi95[2]->Fill(hitz,phi);
					hVxVz95[2]->Fill(vz,vx);	
				}

			}
		}
	}
	TCanvas *radius = new TCanvas("radius","radius",800,600);
	radius->SetLogy();
	r->GetYaxis()->SetTitle("Counts");
	r->GetXaxis()->SetTitle("r[mm]");
	r->Draw();
	temp_string4<<"remoll_"<<temp_string1.str()<<"_radial_28.gif";
	TString fgif = temp_string4.str();
	radius->SaveAs(fgif);
	//radius->SaveAs("ep_radial_default_92.gif");

	temp_string4.str("");
	temp_string4.clear();
	TCanvas *r_photon = new TCanvas("r_photon","radius for photon",800,600);
	r_photon->SetLogy();
	rphoton->GetYaxis()->SetTitle("Counts");
	rphoton->GetXaxis()->SetTitle("r[mm]");
	rphoton->Draw();
	temp_string4<<"remoll_"<<temp_string1.str()<<"_radial_photon_28.gif";
	fgif = temp_string4.str();
	r_photon->SaveAs(fgif);

	temp_string4.str("");
	temp_string4.clear();
	TCanvas *radRate = new TCanvas("radRate","rate weighted radius",800,600);
	radRate->SetLogy();
	r->GetYaxis()->SetTitle("Counts");
	rRate->GetYaxis()->SetTitle("Hz/#muA");
	rRate->GetXaxis()->SetTitle("r[mm]");
	rRate->Draw();
	temp_string4<<"remoll_"<<temp_string1.str()<<"_rRate_28.gif";
	fgif = temp_string4.str();
	radRate->SaveAs(fgif);
	//radRate->SaveAs("ep_rRate_default_28.gif");


	TCanvas *ZPhi_94 = new TCanvas("ZPhi_94","ZPhi_94",800,600);
	ZPhi_94->SetLogy(0);
	zphi94->GetYaxis()->SetTitle("#phi");
	zphi94->GetXaxis()->SetTitle("z[mm]");
	zphi94->Draw("colz");
	temp_string4.str("");
	temp_string4.clear();
	temp_string4<<"remoll_"<<temp_string1.str()<<"_ZPhi_94.gif";
	fgif = temp_string4.str();
	ZPhi_94->SaveAs(fgif);
	//ZPhi_92->SaveAs("ep_rRate_default_92.gif");

	TCanvas *ZPhi_95 = new TCanvas("ZPhi_95","ZPhi_95",800,600);
	ZPhi_95->SetLogy(0);
	zphi95->GetYaxis()->SetTitle("#phi");
	zphi95->GetXaxis()->SetTitle("z[mm]");
	zphi95->Draw("colz");
	temp_string4.str("");
	temp_string4.clear();
	temp_string4<<"remoll_"<<temp_string1.str()<<"_ZPhi_95.gif";
	fgif = temp_string4.str();
	ZPhi_95->SaveAs(fgif);
	//ZPhi_95->SaveAs("ep_rRate_default_95.gif");

	temp_string4.str("");
	temp_string4.clear();
	TCanvas *c_r = new TCanvas("c_r","radius distribution",800,600);
	c_r->SetLogy();
	r->SetLineColor(kMagenta);
	r->Draw();
	hit_r[0]->SetLineColor(kBlue);
	hit_r[0]->Draw("SAME");
	hit_r[1]->SetLineColor(kRed);
	hit_r[1]->Draw("SAME");
	hit_r[2]->SetLineColor(kGreen);
	hit_r[2]->Draw("SAME");
	r->GetYaxis()->SetTitle("Counts");
	r->GetXaxis()->SetTitle("r[mm]");
	temp_string4<<"remoll_"<<temp_string1.str()<<"_hitr_28.gif";
	fgif = temp_string4.str();
	c_r->SaveAs(fgif);

	temp_string4.str("");
	temp_string4.clear();
	TCanvas *c_zphi94 = new TCanvas("c_zphi94","Phi vs Z dist on det 94",800,600);
	c_zphi94->Divide(1,4);
	c_zphi94->cd(1);
	zphi94->Draw("COLZ");
	c_zphi94->cd(2);
	h_zphi94[0]->Draw("colz");
	c_zphi94->cd(3);
	h_zphi94[1]->Draw("colz");
	c_zphi94->cd(4);
	h_zphi94[2]->Draw("colz");
	zphi94->GetYaxis()->SetTitle("#phi");
	zphi94->GetXaxis()->SetTitle("z[mm]");
	temp_string4<<"remoll_"<<temp_string1.str()<<"_zphi_94.gif";
	fgif = temp_string4.str();
	c_zphi94->SaveAs(fgif);

	temp_string4.str("");
	temp_string4.clear();
	TCanvas *c_zphi95 = new TCanvas("c_zphi95","Phi vs Z dist on det 95",800,600);
	c_zphi95->Divide(1,4);
	c_zphi95->cd(1);
	zphi95->Draw("COLZ");
	c_zphi95->cd(2);
	h_zphi95[0]->Draw("colz");
	c_zphi95->cd(3);
	h_zphi95[1]->Draw("colz");
	c_zphi95->cd(4);
	h_zphi95[2]->Draw("colz");
	zphi95->GetYaxis()->SetTitle("#phi");
	zphi95->GetXaxis()->SetTitle("z[mm]");
	temp_string4<<"remoll_"<<temp_string1.str()<<"_zphi_95.gif";
	fgif = temp_string4.str();
	c_zphi95->SaveAs(fgif);

	temp_string4.str("");
	temp_string4.clear();
	TCanvas *c_vxvz94 = new TCanvas("c_vxvz94","Vertex dist for det 94",800,600);
	c_vxvz94->Divide(1,3);
	c_vxvz94->cd(1);
	hVxVz94[0]->Draw("COLZ");
	c_vxvz94->cd(2);
	hVxVz94[1]->Draw("colz");
	c_vxvz94->cd(3);
	hVxVz94[2]->Draw("colz");
	hVxVz94[0]->GetYaxis()->SetTitle("hit.vx[mm]");
	hVxVz94[0]->GetXaxis()->SetTitle("hit.z[mm]");
	temp_string4<<"remoll_"<<temp_string1.str()<<"_VxVz_94.gif";
	fgif = temp_string4.str();
	c_vxvz94->SaveAs(fgif);

	temp_string4.str("");
	temp_string4.clear();
	TCanvas *c_vxvz95 = new TCanvas("c_vxvz95","Vertex dist for det 95",800,600);
	c_vxvz95->Divide(1,3);
	c_vxvz95->cd(1);
	hVxVz95[0]->Draw("COLZ");
	c_vxvz95->cd(2);
	hVxVz95[1]->Draw("colz");
	c_vxvz95->cd(3);
	hVxVz95[2]->Draw("colz");
	hVxVz95[0]->GetYaxis()->SetTitle("#phi");
	hVxVz95[0]->GetXaxis()->SetTitle("z[mm]");
	temp_string4<<"remoll_"<<temp_string1.str()<<"_VxVz_95.gif";
	fgif = temp_string4.str();
	c_vxvz95->SaveAs(fgif);

	r->Write();
	rphoton->Write();
	rRate->Write();
	zphi94->Write();
	zphi95->Write();
	hit_r[0]->Write();
	hit_r[1]->Write();
	hit_r[2]->Write();
	h_zphi94[0]->Write();
	h_zphi94[1]->Write();
	h_zphi94[2]->Write();
	h_zphi95[0]->Write();
	h_zphi95[1]->Write();
	h_zphi95[2]->Write();
	hVxVz94[0]->Write();
	hVxVz94[1]->Write();
	hVxVz94[2]->Write();
	hVxVz95[0]->Write();
	hVxVz95[1]->Write();
	hVxVz95[2]->Write();
	radius->Write();
	r_photon->Write();
	radRate->Write();
	ZPhi_94->Write();
	ZPhi_95->Write();
	c_r->Write();
	c_zphi94->Write();
	c_zphi95->Write();
	c_vxvz94->Write();
	c_vxvz95->Write();
	fh_out->Close();

	//z->Draw();
// return 0;
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
