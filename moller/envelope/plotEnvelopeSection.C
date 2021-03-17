void plotEnvelopeSection()
{
	gStyle->SetOptStat(0);
	string Dstr = "Collar1";
	float xmin1=-120;
	float xmax1=-770;
	float ymin1=-5;
	float ymax1=260;

	//For inner photon pad4
	float xmin2=5;
	float xmax2=-45;
	float ymin2=-10;
	float ymax2=20;

	TFile fee("envelope_ee_ew_region4.root");
	TFile fep("envelope_ep_ew_region4.root");
	TH2D *histee = new TH2D();
	TH2D *histep = new TH2D();
	TGraph *gcee,*gcep,*gceeS,*gcepS,*gcop,*gcip;
	histee=(TH2D*)fee.Get(Form("h_%s_electron_d_xy",Dstr.c_str()));
	histep=(TH2D*)fep.Get(Form("h_%s_electron_d_xy",Dstr.c_str()));
	histee->RebinX(2);
	histee->RebinY(2);
	histep->RebinX(2);
	histep->RebinY(2);
	TString pngfile =Form("./plot/CSection_region4_1_%s.png",Dstr.c_str());
	TString ftxtee = Form("./txtfile/CSignal_ee_rw_newgeo_%s.txt",Dstr.c_str());
	TString ftxtep = Form("./txtfile/CSignal_ep_rw_newgeo_%s.txt",Dstr.c_str());
	TString ftxteeS = Form("./txtfile/CEE_EW_Region4_%s.txt",Dstr.c_str());
	TString ftxtepS = Form("./txtfile/CEP_EW_Region4_%s.txt",Dstr.c_str());
	TString ftxtop = Form("./txtfile/CPhoton_EW_Region4_%s.txt",Dstr.c_str());
	ifstream txtfileee,txtfileep,txtfileop,txtfileip,txtfileeeS,txtfileepS;
	txtfileee.open(ftxtee);
	txtfileep.open(ftxtep);
	txtfileeeS.open(ftxteeS);//S for Segmented
	txtfileepS.open(ftxtepS);
	txtfileop.open(ftxtop);
	int Nee,Nep,Nop,Nip,NeeS,NepS;
	txtfileee>>Nee;
	txtfileep>>Nep;
	txtfileeeS>>NeeS;
	txtfileepS>>NepS;
	txtfileop>>Nop;
	txtfileip>>Nip;
	//txtfileip>>Nip;
	std::cout<<"Nee "<<Nee<<" Nep "<<Nep<<" Nop  "<<Nop<<" Nip "<<Nip<<" NeeS "<<NeeS<<" NepS "<<NepS<<std::endl;
	const int NoEee = (const int) Nee;
	const int NoEep = (const int) Nep;
	const int NoEeeS = (const int) NeeS;
	const int NoEepS = (const int) NepS;
	const int NoEop = (const int) Nop;
	const int NoEip = (const int) Nip;
	//float x[NoE]{},y[NoE]{},z[NoE]{};
	float xee[NoEee],yee[NoEee],zee[NoEee],xep[NoEep],yep[NoEep],zep[NoEep];
	float xeeS[NoEeeS],yeeS[NoEeeS],zeeS[NoEeeS],xepS[NoEepS],yepS[NoEepS],zepS[NoEepS];
	float xop[NoEop],yop[NoEop],zop[NoEop],xip[NoEip],yip[NoEip],zip[NoEip];
	int ij=0;
	while(txtfileeeS>>xeeS[ij]>>yeeS[ij]>>zeeS[ij]){
			//std::cout<<"ij "<<ij<<" x "<<x[ij]<<" y "<<y[ij]<<" z "<<z[ij]<<std::endl;
			ij++;
			}

	txtfileeeS.close();
	ij=0;
	while(txtfileepS>>xepS[ij]>>yepS[ij]>>zepS[ij]){
			//std::cout<<"ij "<<ij<<" x "<<x[ij]<<" y "<<y[ij]<<" z "<<z[ij]<<std::endl;
			ij++;
			}

	txtfileepS.close();
	ij=0;
	while(txtfileee>>xee[ij]>>yee[ij]>>zee[ij]){
			//std::cout<<"ij "<<ij<<" x "<<x[ij]<<" y "<<y[ij]<<" z "<<z[ij]<<std::endl;
			ij++;
			}

	txtfileee.close();
	ij=0;
	while(txtfileep>>xep[ij]>>yep[ij]>>zep[ij]){
			//std::cout<<"ij "<<ij<<" x "<<x[ij]<<" y "<<y[ij]<<" z "<<z[ij]<<std::endl;
			ij++;
			}

	txtfileep.close();
	ij=0;
	while(txtfileop>>xop[ij]>>yop[ij]>>zop[ij]){
			//std::cout<<"ij "<<ij<<" x "<<x[ij]<<" y "<<y[ij]<<" z "<<z[ij]<<std::endl;
			ij++;
			}

	txtfileop.close();
	TLatex tt1 = TLatex();
       	tt1.SetNDC();
	tt1.SetTextSize(0.1);
	TCanvas *plot= new TCanvas("envelope","envelope",1000,600);
	plot->Divide(1,2);
	plot->cd(1);
	//gPad->SetTitle(Form("Moller at %s",Dstr.c_str()));
	//gPad->SetTitleSize(0.05);
	gPad->SetGrid(1,1);
	gPad->SetRightMargin(0.13);
	gPad->SetLeftMargin(0.055);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.08);
	histee->GetXaxis()->SetTitle("x (mm)");
	histee->GetXaxis()->SetTitleSize(0.05);
	histee->GetXaxis()->SetLabelSize(0.05);
	histee->GetXaxis()->SetTitleOffset(0.8);
	histee->GetYaxis()->SetTitle("y (mm)");
	histee->GetYaxis()->SetTitleSize(0.05);
	histee->GetYaxis()->SetLabelSize(0.05);
	histee->GetYaxis()->SetTitleOffset(0.5);
	histee->GetZaxis()->SetLabelSize(0.05);
	histee->GetZaxis()->CenterTitle(1);
	histee->GetZaxis()->SetTitle("Rate*e (arb. unit)");
	histee->GetZaxis()->SetTitleSize(0.05);
	histee->GetZaxis()->SetTitleOffset(.6);
	histee->SetTitle(Form("Moller at %s",Dstr.c_str()));
	histee->Draw("COLZ");
	//hist->SetTitleOffset(-1.0);
	histee->GetXaxis()->SetRangeUser(xmax1,xmin1);
	histee->GetYaxis()->SetRangeUser(ymin1,ymax1);
	gcee=new TGraph(Nee,xee,yee);
	gcee->SetLineColor(kRed);
	gcee->SetLineWidth(2);
	gcee->Draw("SAME L");
	tt1.SetTextColor(kRed);
       	tt1.DrawLatex(0.7, 0.7, "Signal");
	plot->Update();
	gceeS=new TGraph(NeeS,xeeS,yeeS);
	gceeS->SetLineColor(kMagenta);
	gceeS->SetLineWidth(2);
	gceeS->Draw("SAME L");
	tt1.SetTextColor(kMagenta);
       	tt1.DrawLatex(0.7, 0.8, "Segmented");
	plot->Update();
        gcop=new TGraph(Nop,xop,yop);
        gcop->SetLineColor(kBlack);
        gcop->SetLineWidth(2);
        gcop->Draw("SAME L");
	tt1.SetTextColor(kBlack);
       	tt1.DrawLatex(0.7, 0.6, "Photon");
	/*TLine *line1 = new TLine(-53.3,10.3,-88.60,17.68);
	line1->SetLineColor(kGreen);
	line1->SetLineWidth(2);
	line1->Draw();
	TLine *line2 = new TLine(-88.60,17.68,-123.25,28.48);
	line2->SetLineColor(kGreen);
	line2->SetLineWidth(2);
	line2->Draw();
	TLine *line3 = new TLine(-123.25,28.48,-187.22,45.04);
	line3->SetLineColor(kGreen);
	line3->SetLineWidth(2);
	line3->Draw();
	tt1.SetTextColor(kGreen);
       	tt1.DrawLatex(0.7, 0.5, "Col4 boundary");*/
        plot->Update();

	plot->cd(2);
	//gPad->SetTitle(Form("EP at %s",Dstr.c_str()));
	gPad->SetGrid(1,1);
	gPad->SetRightMargin(0.13);
	gPad->SetLeftMargin(0.055);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.08);
	histep->GetXaxis()->SetTitle("x (mm)");
	histep->GetXaxis()->SetTitleSize(0.05);
	histep->GetXaxis()->SetLabelSize(0.05);
	histep->GetXaxis()->SetTitleOffset(0.8);
	histep->GetYaxis()->SetTitle("y (mm)");
	histep->GetYaxis()->SetTitleSize(0.05);
	histep->GetYaxis()->SetLabelSize(0.05);
	histep->GetYaxis()->SetTitleOffset(0.5);
	histep->GetZaxis()->SetLabelSize(0.05);
	histep->GetZaxis()->CenterTitle(1);
	histep->GetZaxis()->SetTitle("Rate*e (arb. unit)");
	histep->GetZaxis()->SetTitleSize(0.05);
	histep->GetZaxis()->SetTitleOffset(.6);
	histep->Draw("COLZ");
	histep->SetTitle(Form("EP at %s",Dstr.c_str()));
	//hist->SetTitleOffset(-1.0);
	histep->GetXaxis()->SetRangeUser(xmax1,xmin1);
	histep->GetYaxis()->SetRangeUser(ymin1,ymax1);
	gcep=new TGraph(Nep,xep,yep);
	gcep->SetLineColor(kRed);
	gcep->SetLineWidth(2);
	gcep->Draw("SAME L");
	plot->Update();
	gcepS=new TGraph(NepS,xepS,yepS);
	gcepS->SetLineColor(kMagenta);
	gcepS->SetLineWidth(2);
	gcepS->Draw("SAME L");
	plot->Update();
        gcop->Draw("SAME L");
/*	line1->Draw();
	line2->Draw();
	line3->Draw();*/
        plot->Update();
	/*plot->cd(3);
	gPad->SetRightMargin(0.13);
	gPad->SetLeftMargin(0.075);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.08);
	gPad->SetGrid(1,1);
	histop->GetXaxis()->SetTitle("x (mm)");
	histop->GetXaxis()->SetTitleSize(0.05);
	histop->GetXaxis()->SetLabelSize(0.05);
	histop->GetXaxis()->SetTitleOffset(0.8);
	histop->GetYaxis()->SetTitle("y (mm)");
	histop->GetYaxis()->SetTitleSize(0.05);
	histop->GetYaxis()->SetLabelSize(0.05);
	histop->GetYaxis()->SetTitleOffset(0.8);
	histop->GetZaxis()->SetLabelSize(0.05);
	histop->GetZaxis()->CenterTitle(1);
	histop->GetZaxis()->SetTitle("Rate (arb. unit)");
	histop->GetZaxis()->SetTitleSize(0.05);
	histop->GetZaxis()->SetTitleOffset(.8);
	histop->Draw("COLZ");
	histop->SetTitle(Form("Outer Photon at %s",Dstr.c_str()));
	//hist->SetTitleOffset(-1.0);
	histop->GetXaxis()->SetRangeUser(xmax1,xmin1);
	histop->GetYaxis()->SetRangeUser(ymin1,ymax1);
	gcop=new TGraph(Nop,xop,yop);
	gcop->SetLineColor(kRed);
	gcop->SetLineWidth(2);
	gcop->Draw("SAME L");
	plot->Update();
	plot->cd(4);
	gPad->SetRightMargin(0.13);
	gPad->SetLeftMargin(0.075);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.08);
	gPad->SetGrid(1,1);
	histip->GetXaxis()->SetTitle("x (mm)");
	histip->GetXaxis()->SetTitleSize(0.05);
	histip->GetXaxis()->SetLabelSize(0.05);
	histip->GetXaxis()->SetTitleOffset(0.8);
	histip->GetYaxis()->SetTitle("y (mm)");
	histip->GetYaxis()->SetTitleSize(0.05);
	histip->GetYaxis()->SetLabelSize(0.05);
	histip->GetYaxis()->SetTitleOffset(0.8);
	histip->GetZaxis()->SetLabelSize(0.05);
	histip->GetZaxis()->CenterTitle(1);
	histip->GetZaxis()->SetTitle("Rate (arb. unit)");
	histip->GetZaxis()->SetTitleSize(0.05);
	histip->GetZaxis()->SetTitleOffset(.8);
	histip->Draw("COLZ");
	histip->SetTitle(Form("Inner Photon at %s",Dstr.c_str()));
	//hist->SetTitleOffset(-1.0);
	histip->GetXaxis()->SetRangeUser(xmax2,xmin2);
	histip->GetYaxis()->SetRangeUser(ymin2,ymax2);*/
	/*gcip=new TGraph(Nip,xip,yip);
	gcip->SetLineColor(kRed);
	gcip->SetLineWidth(2);
	gcip->Draw("SAME L");*/
	plot->Update();
	plot->SaveAs(pngfile);
	plot->SaveAs("testCombined.png");

}
