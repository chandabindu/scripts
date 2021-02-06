int plot()
{
	gStyle->SetOptStat(0);
        TH2D *h_ee = new TH2D();//"h_ee","MOLLLER::X vs Y", 2400,-1200,1200,2400,-1200,1200);
        TH2D *h_ep = new TH2D();//"h_ep","ep::X vs Y", 2400,-1200,1200,2400,-1200,1200);
        TH2D *h_beam = new TH2D();//"h_beam","beam::X vs Y", 2400,-1200,1200,2400,-1200,1200);
        TH1D *ha_ee = new TH1D("ha_ee","EE::Radial rate distributins", 700,500,1200);
        TH1D *ha_ep = new TH1D("ha_ep","EP::Radial rate distributins", 700,500,1200);
        TH1D *ha_beam = new TH1D("ha_beam","Beam::Radial rate distributins", 700,500,1200);
        TFile f1("./added_ee_withR.root");
        TFile f2("./added_ep_withR.root");
        TFile f3("./added_beam_withR.root");
        
	float factor1 = 1.e-9/70./208.;
	float factor2 = 1.e-9/70./210.;
	float factor3= 1.e-9/70./1018.;

        h_ee=(TH2D*)f1.Get("XY_GEMPlane");
        h_ep=(TH2D*)f2.Get("XY_GEMPlane");
        h_beam=(TH2D*)f3.Get("XY_GEMPlane");
        ha_ee=(TH1D*)f1.Get("r_GEMPlane");
        ha_ep=(TH1D*)f2.Get("r_GEMPlane");
        ha_beam=(TH1D*)f3.Get("r_GEMPlane");


	h_ee->Scale(factor1);//ha_1->Scale(factor1);htheta_1->Scale(factor1);
	h_ep->Scale(factor2);//ha_2->Scale(factor2);htheta_2->Scale(factor2);
	h_beam->Scale(factor3);//ha_3->Scale(factor3);htheta_3->Scale(factor3);
	ha_ee->Scale(factor1);
	ha_ep->Scale(factor2);
	ha_beam->Scale(factor3);

	TH2D *h4 = new TH2D("summed","ep+ee",2400,-1200,1200,2400,1200,1200);
	h4->Add(h_ee);
	h4->Add(h_ep);

	TH1D *h5 = new TH1D("summed1D","1D summed ee+ep",700,500,1200);
	h5->Add(ha_ee);
	h5->Add(ha_ep);
	ha_ee->Rebin(5.);
	ha_ep->Rebin(5.);
	ha_beam->Rebin(5.);
	h5->Rebin(5.);

	double SUM = 0;
	for(int binx=1;binx<=2400;binx++){
		for(int biny=1;biny<=2400;biny++){
			double r = sqrt(h4->GetXaxis()->GetBinCenter(binx)*h4->GetXaxis()->GetBinCenter(binx)
					+h4->GetYaxis()->GetBinCenter(biny)*h4->GetYaxis()->GetBinCenter(biny));
			if(r>=580 && r<=1182.5){
				double phi = atan(h4->GetYaxis()->GetBinCenter(biny)/h4->GetXaxis()->GetBinCenter(binx))*180./acos(-1);
				if(abs(phi)<=12.25 && (h4->GetXaxis()->GetBinCenter(binx)<0)){
					SUM=SUM+h4->GetBinContent(binx,biny);
				}
			}
		}
	}
	std::cout<<"Integral "<<SUM<<std::endl;
	//std::cout<<"Without entries "<<without_w->GetEntries()<<"  with "<<with_w->GetEntries()<<std::endl;
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
        //h_nopipe->GetYaxis()->SetRangeUser(2.0e3,1.0e5);
	plot->Divide(2,2);
	plot->cd(1);
	gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
	gPad->SetLogz();
        h_ee->Draw("colz");
	plot->cd(2);
	gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
	gPad->SetLogz();
        h_ep->Draw("colz");
	plot->cd(3);
	gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
	gPad->SetLogz();
        h_beam->Draw("colz");
	plot->cd(4);
	gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
	gPad->SetLogz();
        h4->Draw("colz");
	h_ee->GetXaxis()->SetTitle("x (mm)");
	h_ee->GetYaxis()->SetTitle("y (mm)");
	h_ep->GetXaxis()->SetTitle("x (mm)");
	h_ep->GetYaxis()->SetTitle("y (mm)");
	h_beam->GetXaxis()->SetTitle("x (mm)");
	h_beam->GetYaxis()->SetTitle("y (mm)");
	h4->GetXaxis()->SetTitle("x (mm)");
	h4->GetYaxis()->SetTitle("y (mm)");
	/*latex.SetTextColor(kRed);
	latex.DrawLatex(1000,30,Form("Beam"));
	latex.SetTextColor(kBlue);
	latex.DrawLatex(1000,10,Form("EP"));
	latex.SetTextColor(kBlack);
	latex.DrawLatex(1000,3,Form("MOLLER"));
	latex.SetTextColor(kGreen);
	latex.DrawLatex(1000,1,Form("MOLLER+ep"));*/
	//plot->cd();
	plot->SaveAs("Hits2D_1.png");

        TCanvas *plot2 = new TCanvas("Plot2","plot2",800,600);
	gStyle->SetOptStat(0);
	latex.SetTextSize(0.05);
	latex.SetTextAlign(13);

	//latex.DrawLatex(0.2,0.2,"Red - with old geometry");
/*
	TLatex latex1;
	latex1.SetTextSize(0.05);
	latex1.SetTextAlign(13);
	latex1.DrawLatex(0.5,without_w->GetMaximum()/1.5,"Blue - with new geometry");
*/	
        //h_nopipe->GetYaxis()->SetRangeUser(2.0e3,1.0e5);
	plot2->Divide(2,2);
	plot2->cd(1);
	ha_ee->GetXaxis()->SetTitle("R (mm)");
	ha_ee->GetYaxis()->SetTitle("Rate (GHz/(5mm.#muA.Sector))");
	gPad->SetLogy();
	ha_ee->SetLineColor(kRed);
        ha_ee->Draw("HIST");
	latex.SetTextColor(kRed);
	latex.DrawLatex(800,1.e-6,Form("MOLLER"));
	plot2->cd(2);
	gPad->SetLogy();
	ha_ep->SetLineColor(kBlue);
	ha_ep->GetXaxis()->SetTitle("R (mm)");
	ha_ep->GetYaxis()->SetTitle("Rate (GHz/(5mm.#muA.Sector))");
        ha_ep->Draw("HIST");
	latex.SetTextColor(kBlue);
	latex.DrawLatex(800,1.e-5,Form("EP"));
	plot2->cd(3);
	gPad->SetLogy();
	ha_beam->GetXaxis()->SetTitle("R (mm)");
	ha_beam->GetYaxis()->SetTitle("Rate (GHz/(5mm.#muA.Sector))");
	ha_beam->SetLineColor(kBlack);
        ha_beam->Draw("HIST");
	latex.SetTextColor(kBlack);
	latex.DrawLatex(800,1.e-4,Form("Beam"));
	plot2->cd(4);
	gPad->SetLogy();
	h5->GetXaxis()->SetTitle("R (mm)");
	h5->GetYaxis()->SetTitle("Rate (GHz/(5mm.#muA.Sector))");
	h5->SetLineColor(kMagenta);
	ha_beam->Draw("HIST");
	ha_ee->Draw("HIST && SAME");
	ha_ep->Draw("HIST && SAME");
        h5->Draw("HIST && SAME");
	h_ee->GetXaxis()->SetTitle("x (mm)");
	h_ee->GetYaxis()->SetTitle("y (mm)");
	latex.SetTextColor(kMagenta);
	latex.DrawLatex(800,1.e-4,Form("MOLLER+EP"));
	plot->cd();
	plot2->SaveAs("Hits1D.png");

/*        TCanvas *plot1 = new TCanvas("Plot1","Radial distributions for all particles",800,600);
	gStyle->SetOptStat(0);
	plot1->SetLogy();
        //ha_nopipe->GetYaxis()->SetRangeUser(5.0e3,1.0e5);
        ha_1->Draw("HIST");
        ha_2->Draw("same HIST");
        ha_3->Draw("same HIST");
	ha_1->GetXaxis()->SetTitle("Radius at main det. plane (mm)");
	ha_1->GetYaxis()->SetTitle("Rate (GHz/(1mm.#muA.Sector))");
	latex.SetTextColor(kRed);
	latex.DrawLatex(1200,100,Form("Beam"));
	latex.SetTextColor(kBlue);
	latex.DrawLatex(1200,40,Form("EP"));
	latex.SetTextColor(kBlack);
	latex.DrawLatex(1200,15,Form("MOLLER"));
	plot1->cd();
	plot1->SaveAs("plot/Unshielded_Photon_Rate_acceptance.png");

        TCanvas *plot2 = new TCanvas("Plot2","Vertex theta distributons",800,600);
	gStyle->SetOptStat(0);
	plot2->SetLogy();
        //ha_nopipe->GetYaxis()->SetRangeUser(5.0e3,1.0e5);
        htheta_1->Draw("HIST");
        htheta_2->Draw("same HIST");
        htheta_3->Draw("same HIST");
	htheta_1->GetXaxis()->SetTitle("#theta_{vertex} (rad) ");
	htheta_1->GetYaxis()->SetTitle("Rate (GHz/(1mm.#muA.Sector))");
	latex.SetTextColor(kRed);
	latex.DrawLatex(0.005,0.01,Form("Beam"));
	latex.SetTextColor(kBlue);
	latex.DrawLatex(0.005,0.006,Form("EP"));
	latex.SetTextColor(kBlack);
	latex.DrawLatex(0.005,0.0035,Form("MOLLER"));
	plot2->cd();
	plot2->SaveAs("plot/Unshielded_Theta_vertex_acceptance.png");*/
return 0;
}


