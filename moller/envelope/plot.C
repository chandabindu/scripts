void plot()
{
	gStyle->SetOptStat(0);
	TFile f1("envelope_OP_R4.root");
	TH2D *hist = new TH2D();
	TGraph *gc;
	string Dstr = "Collar1";
	hist=(TH2D*)f1.Get(Form("h_%s_electron_d_xy",Dstr.c_str()));
	hist->RebinX(2);
	hist->RebinY(2);
	TString pngfile =Form("./plot/CPhoton_EW_Region4_%s.png",Dstr.c_str());
	TString ftxt = Form("./txtfile/CPhoton_EW_Region4_%s.txt",Dstr.c_str());
	ifstream txtfile;
	txtfile.open(ftxt);
	int N;
	txtfile>>N;
	std::cout<<"N "<<N<<std::endl;
	const int NoE = (const int) N;
	//float x[NoE]{},y[NoE]{},z[NoE]{};
	float x[NoE],y[NoE],z[NoE];
	int ij=0;
	while(txtfile>>x[ij]>>y[ij]>>z[ij]){
			//std::cout<<"ij "<<ij<<" x "<<x[ij]<<" y "<<y[ij]<<" z "<<z[ij]<<std::endl;
			ij++;
			}

	txtfile.close();
	TCanvas *plot= new TCanvas("envelope","envelope",800,600);
	gPad->SetRightMargin(0.125);
	gPad->SetLeftMargin(0.125);
	gPad->SetGrid(1,1);
	hist->Draw("COLZ");
	hist->SetTitle(Form("Photon envelope at %s",Dstr.c_str()));
	//hist->SetTitle(Form("Moller envelope at Col4Exit"));
	//hist->SetTitleOffset(-1.0);
	hist->GetXaxis()->SetTitle("x (mm)");
	hist->GetYaxis()->SetTitle("y (mm)");
	hist->GetXaxis()->SetRangeUser(-370,-120);
	hist->GetYaxis()->SetRangeUser(-10,90);
	gc=new TGraph(N,x,y);
	gc->SetLineColor(kRed);
	gc->SetLineWidth(2);
	gc->Draw("SAME L");
	plot->Update();
	plot->SaveAs(pngfile);
	plot->SaveAs("test.png");

}
