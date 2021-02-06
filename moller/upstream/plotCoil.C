void setcolor(){

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

}

void plotCoil()
{
gROOT->Reset();
   //TChanin *T = new TChain("T");
   //T->Add("outfile2.root");
   TFile *f = new TFile("outfile10.root");
   gStyle->SetOptStat(0);
   gStyle->SetStatH(0.3);
   gStyle->SetStatW(0.3);
   gStyle->SetTitleH(0.09);
   gStyle->SetTitleW(0.3);
   gStyle->SetLabelSize(0.04,"x");
   gStyle->SetLabelSize(0.04,"y");
   gROOT->ForceStyle();


   const int nEnergy = 5;
   const int nParticle = 4;
   const int nDet = 14;
   TH2F *h = new TH2F("h","h",480,800,3200,600,0,300);
   TH2F *h1 = new TH2F("h1","h1",480,800,3200,5,25,50);
   TH2F *h2 = new TH2F("h2","h2",480,800,3200,380,-190,190);
   TH1F *h4 = new TH1F("h4","h4",110,0,11000);
   //h->GetZaxis()->SetRangeUser(0,0.01);
   //h->SetMinimum(0);
   //h->SetMinimum(0.1);
  string out = "outfile.pdf";

  TCanvas c1;

  c1.Print(Form("%s[",out.c_str()),"pdf");
  //c1.Divide(1,7);
  string sParticle[nParticle] = {"electron","positron","photon","neutron"};
  string Ecut[nEnergy]={"Total","p<1","p>=1 && p<10","p>=10 && p<100","p>=100"};
  string sDet[nDet] = {"coil1","coil2","coil3","coil4","coil5","coil6","coil7","epoxy1","epoxy2","epoxy3","epoxy4","epoxy5","epoxy6","epoxy7"};

  //string septant[4]={"p<1","p>=1 && p<10","p>=10 && p<100","p>=100"};

  for(int ij=0;ij<nDet;ij++){
	for(int ip=0;ip<nParticle;ip++){
		for(int ik=0;ik<nEnergy;ik++){
			c1.cd(1);
			setcolor();
			h=(TH2F*)f->Get(Form("h_%s_%s_E%d_u_rz",sDet[ij].c_str(),sParticle[ip].c_str(),ik));
			h->GetZaxis()->SetTitleOffset(0.5);
			h->SetTitle(Form("%s %s on %s",Ecut[ik].c_str(),sParticle[ip].c_str(),sDet[ij].c_str()));
			h->SetTitleSize(0.05);
			h->Draw("colz");
			c1.Print(out.c_str(),"pdf");
			}
	}
   }

  for(int ij=7;ij<nDet;ij++){
	for(int ip=0;ip<nParticle;ip++){
		for(int ik=0;ik<nEnergy;ik++){
			c1.cd(1);
			setcolor();
			h1=(TH2F*)f->Get(Form("h_%s_%s_E%d_ub_rz",sDet[ij].c_str(),sParticle[ip].c_str(),ik));
			h1->GetZaxis()->SetTitleOffset(0.5);
			h1->SetTitle(Form("Under Belly:%s %s on %s",Ecut[ik].c_str(),sParticle[ip].c_str(),sDet[ij].c_str()));
			h1->SetTitleSize(0.05);
			h1->Draw("colz");
			c1.Print(out.c_str(),"pdf");
			}
	}
   }

  for(int ip=0;ip<nParticle;ip++){
	for(int ik=0;ik<nEnergy;ik++){
		c1.cd(1);
		setcolor();
		h2=(TH2F*)f->Get(Form("h_%s_E%d_ub_phiz",sParticle[ip].c_str(),ik));
		h2->GetZaxis()->SetTitleOffset(0.5);
		h2->SetTitle(Form("Under Belly:%s on %s",Ecut[ik].c_str(),sParticle[ip].c_str()));
		h2->SetTitleSize(0.05);
		h2->Draw("colz");
		c1.Print(out.c_str(),"pdf");
		}
  }
  for(int ij=0;ij<nDet;ij++){
	for(int ip=0;ip<nParticle;ip++){
		c1.cd(1);
		setcolor();
		h4=(TH1F*)f->Get(Form("h_%s_%s_u_e",sDet[ij].c_str(),sParticle[ip].c_str()));
		h4->SetTitle(Form("%s on %s",sParticle[ip].c_str(),sDet[ij].c_str()));
		h4->SetTitleSize(0.05);
		h4->Draw("colz");
		c1.Print(out.c_str(),"pdf");
		}
	}
   
  for(int ij=7;ij<nDet;ij++){
	for(int ip=0;ip<nParticle;ip++){
		c1.cd(1);
		setcolor();
		h4=(TH1F*)f->Get(Form("h_%s_%s_ub_e",sDet[ij].c_str(),sParticle[ip].c_str()));
		h4->SetTitle(Form("%s on %s",sParticle[ip].c_str(),sDet[ij].c_str()));
		h4->SetTitleSize(0.05);
		h4->Draw("colz");
		c1.Print(out.c_str(),"pdf");
		}
	}

  c1.Print(Form("%s]",out.c_str()),"pdf");



}
