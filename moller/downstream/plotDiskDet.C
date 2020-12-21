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

void plotDiskDet()
{
gROOT->Reset();
   //TChanin *T = new TChain("T");
   //T->Add("outfile2.root");
   TFile *f = new TFile("DSConfig0_ep_disk_28_ring5.root");
   gStyle->SetOptStat(0);
   gStyle->SetStatH(0.3);
   gStyle->SetStatW(0.3);
   gStyle->SetTitleH(0.09);
   gStyle->SetTitleW(0.3);
   gStyle->SetLabelSize(0.04,"x");
   gStyle->SetLabelSize(0.04,"y");
   gROOT->ForceStyle();


   const int nEnergy = 5;
   const int nParticle = 5;
   const int nDet = 8;
   const int nTDet = 1;
   TH2F *h = new TH2F("h","h",200,-1900,-1900,200,-1900,1900);
   TH2F *h1 = new TH2F("h1","h1",120,800,3200,5,25,50);
   TH2F *h2 = new TH2F("h2","h2",120,800,3200,380,-190,190);
   TH1F *h4 = new TH1F("h4","h4",110,0,11000);
   //h->GetZaxis()->SetRangeUser(0,0.01);
   //h->SetMinimum(0);
   //h->SetMinimum(0.1);
  string out = "DSConfig0_ep_disk_28_ring5.pdf";
  ofstream fh_out;
  fh_out.open("DSConfig0_ep_disk_28_ring5.txt");
  

  TCanvas c1;
  c1.SetRightMargin(0.15);

  c1.Print(Form("%s[",out.c_str()),"pdf");
  //c1.Divide(1,7);
  string sParticle[nParticle] = {"All","electron","positron","photon","neutron"};
  string Ecut[nEnergy]={"Total","p<1","p>=1 && p<10","p>=10 && p<100","p>=100"};
  string sDet[nDet] = {"afterCol4","beforeCoil1","afterCoil1","afterCoil2","afterCoil3","beforeCollar","driftRegion","maindet"};
  //string sTDet[nTDet] = {"tube1","tube2","tube3","tube4","tube5","tube6"};
  string sTDet[nTDet] = {"cylinder"};

  //string septant[4]={"p<1","p>=1 && p<10","p>=10 && p<100","p>=100"};

  for(int ij=0;ij<nDet;ij++){
	for(int ip=0;ip<nParticle;ip++){
          if(ip==0){
         c1.cd(1);
       	 //setcolor();
         h=(TH2F*)f->Get(Form("h_%s_d_xy",sDet[ij].c_str()));
         //h->GetZaxis()->SetTitleOffset(0.5);
         h->SetTitle(Form("Total Energy on %s",sDet[ij].c_str()));
         h->SetTitleSize(0.05);
	 fh_out<<"Det "<<Form("%s",sDet[ij].c_str())<<" Max energy "<<h->GetMaximum()<<"Z for Max. "<<h->GetXaxis()->GetBinCenter(h->GetMaximumBin())<<std::endl;
         h->Draw("colz");
         c1.Print(out.c_str(),"pdf");

     	}
             
		for(int ik=0;ik<nEnergy && ip>0;ik++){
			c1.cd(1);
			//setcolor();
			h=(TH2F*)f->Get(Form("h_%s_%s_E%d_d_xy",sDet[ij].c_str(),sParticle[ip].c_str(),ik));
			h->GetZaxis()->SetTitleOffset(0.5);
			h->SetTitle(Form("%s %s on %s",Ecut[ik].c_str(),sParticle[ip].c_str(),sDet[ij].c_str()));
			h->SetTitleSize(0.05);
			h->Draw("colz");
			c1.Print(out.c_str(),"pdf");
			}
	}
   }

  for(int ij=0;ij<nTDet;ij++){
	for(int ip=0;ip<nParticle;ip++){
          if(ip==0){
         c1.cd(1);
       	 //setcolor();
         h=(TH2F*)f->Get(Form("h_%s_t_rz",sTDet[ij].c_str()));
         h->GetZaxis()->SetTitleOffset(0.5);
         h->SetTitle(Form("Total Energy on %s",sTDet[ij].c_str()));
         h->SetTitleSize(0.05);
	 fh_out<<"Det "<<Form("%s",sTDet[ij].c_str())<<" Max energy "<<h->GetMaximum()<<"Z for Max. "<<h->GetXaxis()->GetBinCenter(h->GetMaximumBin())<<std::endl;
         h->Draw("colz");
         c1.Print(out.c_str(),"pdf");

     	}
             
		for(int ik=0;ik<nEnergy && ip>0;ik++){
			c1.cd(1);
			//setcolor();
			h=(TH2F*)f->Get(Form("h_%s_%s_E%d_t_rz",sTDet[ij].c_str(),sParticle[ip].c_str(),ik));
			h->GetZaxis()->SetTitleOffset(0.5);
			h->SetTitle(Form("%s %s on %s",Ecut[ik].c_str(),sParticle[ip].c_str(),sTDet[ij].c_str()));
			h->SetTitleSize(0.05);
			h->Draw("colz");
			c1.Print(out.c_str(),"pdf");
			}
	}
   }
   fh_out.close();

  for(int ip=0;ip<nParticle;ip++){
	for(int ik=0;ik<nEnergy && ip>0;ik++){
		c1.cd(1);
		//setcolor();
		h2=(TH2F*)f->Get(Form("h_%s_E%d_t_phiz",sParticle[ip].c_str(),ik));
		h2->GetZaxis()->SetTitleOffset(0.5);
		h2->SetTitle(Form("On cylindricl det::%s on %s",Ecut[ik].c_str(),sParticle[ip].c_str()));
		h2->SetTitleSize(0.05);
		h2->Draw("colz");
		c1.Print(out.c_str(),"pdf");
		}
  }

  for(int ij=0;ij<nDet;ij++){
	for(int ip=1;ip<nParticle;ip++){
		c1.cd(1);
	///	setcolor();
		h4=(TH1F*)f->Get(Form("h_%s_%s_d_e",sDet[ij].c_str(),sParticle[ip].c_str()));
		h4->SetTitle(Form("%s on %s",sParticle[ip].c_str(),sDet[ij].c_str()));
		h4->SetTitleSize(0.05);
		h4->Draw("colz");
		c1.Print(out.c_str(),"pdf");
		}
	}

  for(int ij=0;ij<nTDet;ij++){
	for(int ip=1;ip<nParticle;ip++){
		c1.cd(1);
	///	setcolor();
		h4=(TH1F*)f->Get(Form("h_%s_%s_t_e",sTDet[ij].c_str(),sParticle[ip].c_str()));
		h4->SetTitle(Form("%s on %s",sParticle[ip].c_str(),sTDet[ij].c_str()));
		h4->SetTitleSize(0.05);
		h4->Draw("colz");
		c1.Print(out.c_str(),"pdf");
		}
	}

  c1.Print(Form("%s]",out.c_str()),"pdf");



}
