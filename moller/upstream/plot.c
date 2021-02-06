int plot()
{
	gStyle->SetOptStat(0);
        TH1F *h_1 = new TH1F("h_1","Photon:with 1mm thich beam pipe", 500,0,2000);
        TH1F *h_2 = new TH1F("h_2","Photon:with 2mm thich beam pipe", 500,0,2000);
        TH1F *h_3 = new TH1F("h_3","Photon:with 3mm thich beam pipe", 500,0,2000);
        TH1F *h_4 = new TH1F("h_4","Photon:with 4mm thich beam pipe", 500,0,2000);
        TH1F *h_5 = new TH1F("h_5","Photon:with 5mm thich beam pipe", 500,0,2000);
        TH1F *h_nopipe = new TH1F("h_nopipe","Photon:with no beam pipe", 500,0,2000);
        TH1F *ha_1 = new TH1F("ha_1","all:with 1mm thich beam pipe", 500,0,2000);
        TH1F *ha_2 = new TH1F("ha_2","all:with 2mm thich beam pipe", 500,0,2000);
        TH1F *ha_3 = new TH1F("ha_3","all:with 3mm thich beam pipe", 500,0,2000);
        TH1F *ha_4 = new TH1F("ha_4","all:with 4mm thich beam pipe", 500,0,2000);
        TH1F *ha_5 = new TH1F("ha_5","all:with 5mm thich beam pipe", 500,0,2000);
        TH1F *ha_nopipe = new TH1F("ha_nopipe","all:with no beam pipe", 500,0,2000);
        TFile f1("/lustre/expphy/volatile/halla/parity/chandan/sim_out/default11/remoll_beam_bare_out.root");//thickeness 1mm
        TFile f2("/lustre/expphy/volatile/halla/parity/chandan/sim_out/default2/remoll_beam_bare_out.root");//thickeness 2mm
        TFile f3("/lustre/expphy/volatile/halla/parity/chandan/sim_out/default1/remoll_beam_bare_out.root");//thickeness 3mm
        TFile f4("/lustre/expphy/volatile/halla/parity/chandan/sim_out/default4/remoll_beam_bare_out.root");//thickeness 4mm
        TFile f5("/lustre/expphy/volatile/halla/parity/chandan/sim_out/default5/remoll_beam_bare_out.root");//thickeness 5mm
        TFile f6("/lustre/expphy/volatile/halla/parity/chandan/sim_out/default3/remoll_beam_bare_out.root");//thickeness nopipe


	//Int_t photon_int[6]={0};
	//Int_t all_int[6]={0};
	Float_t factor = (6.25e12)*(1.e-9);
	Float_t nbeam = 1.25e7;
	Int_t xmin=600,xmax=1250;
        h_1=(TH1F*)f1.Get("rphoton");
        h_2=(TH1F*)f2.Get("rphoton");
        h_3=(TH1F*)f3.Get("rphoton");
        h_4=(TH1F*)f4.Get("rphoton");
        h_5=(TH1F*)f5.Get("rphoton");
        h_nopipe=(TH1F*)f6.Get("rphoton");
        ha_1=(TH1F*)f1.Get("r");
        ha_2=(TH1F*)f2.Get("r");
        ha_3=(TH1F*)f3.Get("r");
        ha_4=(TH1F*)f4.Get("r");
        ha_5=(TH1F*)f5.Get("r");
        ha_nopipe=(TH1F*)f6.Get("r");

	Int_t bmin = h_nopipe->GetXaxis()->FindBin(xmin);
	Int_t bmax = h_nopipe->GetXaxis()->FindBin(xmax);
        std::cout<<"bin min "<<bmin <<"  bin max  "<<bmax<<std::endl; 	
	std::cout<<"Photon counts between 600 -1250 mm (GHz/microA):for no pipe - "<<(factor/nbeam)*h_nopipe->Integral(bmin,bmax)<<", With 1mm pipe - "<<(factor/nbeam)*h_1->Integral(bmin,bmax)<<", With 2mm pipe - "<<(factor/nbeam)*h_2->Integral(bmin,bmax)<<", With 3mm pipe - "<<(factor/nbeam)*h_3->Integral(bmin,bmax)<<", With 4mm pipe - "<<(factor/nbeam)*h_4->Integral(bmin,bmax)<<", With 5mm pipe - "<<(factor/nbeam)*h_5->Integral(bmin,bmax)<<std::endl;
	std::cout<<"All particle counts between 600 -1250 mm (GHz/microA):for no pipe - "<<(factor/nbeam)*ha_nopipe->Integral(bmin,bmax)<<", With 1mm pipe - "<< (factor/nbeam)*ha_1->Integral(bmin,bmax)<<", With 2mm pipe - "<<(factor/nbeam)*ha_2->Integral(bmin,bmax)<<", With 3mm pipe - "<<(factor/nbeam)*ha_3->Integral(bmin,bmax)<<", With 4mm pipe - "<<(factor/nbeam)*ha_4->Integral(bmin,bmax)<<", With 5mm pipe - "<<(factor/nbeam)*ha_5->Integral(bmin,bmax)<<std::endl;
	//std::cout<<"Without entries "<<without_w->GetEntries()<<"  with "<<with_w->GetEntries()<<std::endl;
        h_1->SetLineColor(kRed);
        h_2->SetLineColor(kBlue);
        h_3->SetLineColor(kGreen);
        h_4->SetLineColor(kBlack);
        h_5->SetLineColor(kMagenta);
        h_nopipe->SetLineColor(kGreen+4);
        ha_1->SetLineColor(kRed);
        ha_2->SetLineColor(kBlue);
        ha_3->SetLineColor(kGreen);
        ha_4->SetLineColor(kBlack);
        ha_5->SetLineColor(kMagenta);
        ha_nopipe->SetLineColor(kGreen+4);
        TCanvas *plot = new TCanvas("Plot","plot",800,600);
	gStyle->SetOptStat(0);
	TLatex latex;
	latex.SetTextSize(0.05);
	latex.SetTextAlign(13);

	//latex.DrawLatex(0.2,0.2,"Red - with old geometry");
/*
	TLatex latex1;
	latex1.SetTextSize(0.05);
	latex1.SetTextAlign(13);
	latex1.DrawLatex(0.5,without_w->GetMaximum()/1.5,"Blue - with new geometry");
*/	
	plot->SetLogy();
        h_nopipe->Draw();
        h_nopipe->GetYaxis()->SetRangeUser(2.0e3,1.0e5);
        h_1->Draw("same");
        h_2->Draw("same");
        h_3->Draw("same");
        h_4->Draw("same");
        h_5->Draw("same");
	latex.SetTextColor(kGreen+4);
	latex.DrawLatex(400,90000,Form("without beampipe - %0.1f GHz/#muA",(factor/nbeam)*h_nopipe->Integral(bmin,bmax)));
	latex.SetTextColor(kRed);
	latex.DrawLatex(400,65000,Form("with 1mm thickness - %0.1f GHz/#muA",(factor/nbeam)*h_1->Integral(bmin,bmax)));
	latex.SetTextColor(kBlue);
	latex.DrawLatex(400,48000,Form("with 2mm thickness - %0.1f GHz/#muA",(factor/nbeam)*h_2->Integral(bmin,bmax)));
	latex.SetTextColor(kGreen);
	latex.DrawLatex(400,35000,Form("with 3mm thickness - %0.1f GHz/#muA",(factor/nbeam)*h_3->Integral(bmin,bmax)));
	latex.SetTextColor(kBlack);
	latex.DrawLatex(400,25000,Form("with 4mm thickness - %0.1f GHz/#muA",(factor/nbeam)*h_4->Integral(bmin,bmax)));
	latex.SetTextColor(kMagenta);
	latex.DrawLatex(400,18000,Form("with 5mm thickness - %0.1f GHz/#muA",(factor/nbeam)*h_5->Integral(bmin,bmax)));
	plot->cd();
	plot->SaveAs("upstream_beampipe_ThicknessStudy_photon.png");

        TCanvas *plot1 = new TCanvas("Plot1","Radial distributions for all particles",800,600);
	gStyle->SetOptStat(0);
	plot1->SetLogy();
        ha_nopipe->Draw();
        ha_nopipe->GetYaxis()->SetRangeUser(5.0e3,1.0e5);
        ha_1->Draw("same");
        ha_2->Draw("same");
        ha_3->Draw("same");
        ha_4->Draw("same");
        ha_5->Draw("same");
	latex.SetTextColor(kGreen+4);
	latex.DrawLatex(400,90000,Form("without beampipe - %0.1f GHz/#muA",(factor/nbeam)*ha_nopipe->Integral(bmin,bmax)));
	latex.SetTextColor(kRed);
	latex.DrawLatex(400,76000,Form("with 1mm thickness - %0.1f GHz/#muA",(factor/nbeam)*ha_1->Integral(bmin,bmax)));
	latex.SetTextColor(kBlue);
	latex.DrawLatex(400,64000,Form("with 2mm thickness - %0.1f GHz/#muA",(factor/nbeam)*ha_2->Integral(bmin,bmax)));
	latex.SetTextColor(kGreen);
	latex.DrawLatex(400,53000,Form("with 3mm thickness - %0.1f GHz/#muA",(factor/nbeam)*ha_3->Integral(bmin,bmax)));
	latex.SetTextColor(kBlack);
	latex.DrawLatex(400,45000,Form("with 4mm thickness - %0.1f GHz/#muA",(factor/nbeam)*ha_4->Integral(bmin,bmax)));
	latex.SetTextColor(kMagenta);
	latex.DrawLatex(400,38000,Form("with 5mm thickness - %0.1f GHz/#muA",(factor/nbeam)*ha_5->Integral(bmin,bmax)));
	plot1->cd();
	plot1->SaveAs("upstream_beampipe_ThicknessStudy_all.png");
return 0;
}


