double YIntegrate(int, TH2D*);
int CalculateBoundary()
{
	gStyle->SetOptStat(0);
	TLatex latex;
	latex.SetTextSize(0.05);
	latex.SetTextAlign(13);
        TFile f1("./envelope_ee.root");
        //TFile f2("./envelope_ep.root");
        TH2D *hist = new TH2D();//"ha_ee","EE::Radial rate distributins", 700,500,1200);
        //TH2D *hist1 = new TH2D();//"ha_ee","EE::Radial rate distributins", 700,500,1200);
	TGraph *gc;
        
	string out = "Boundary_ee.pdf";
	TCanvas c1;
	c1.Print(Form("%s[",out.c_str()),"pdf");
	const int nDet=32;
	//string sDet[nDet] = {"tgtPbWall","tgtPbCollar","Col2Exit","USCoil1","USCoil2"};
	string sDet[nDet] = {"tgtPbWall","tgtPbCollar","Col2Exit","USCoil1","USCoil2","USCoil3","USCoil4","USCoil5","Col4Exit","DSCoil1","DSCoil2","DSCoil3","DSCoil4","DSCoil5","DSCoil6","DSCoil7","DSCoil8","DSCoil9","DSCoil10","DSCoil11","DSCoil12","DSCoil13","DSCoil14","DSCoil15","DSCoil16","Collar1","Drift","Collar2","GEM4","maindet","piondet","HallExit"};
	for(int ndet=0;ndet<nDet;ndet++){
	//c1.cd(1);
	//c1.Clear();
        hist=(TH2D*)f1.Get(Form("h_%s_electron_d_xy",sDet[ndet].c_str()));
        //hist1=(TH2D*)f2.Get(Form("h_%s_electron_d_xy",sDet[ndet].c_str()));
	//hist->Add(hist1);
        //hist=(TH2D*)f1.Get("h_Collar1_electron_d_xy");
	std::cout<<"Before rebin: X bins "<<hist->GetNbinsX()<<" Y bins "<<hist->GetNbinsY()<<std::endl;
	double precision=0.000001*hist->Integral();//1% precision
	int NBins = hist->GetNbinsX();
	double threshold = precision/NBins;
	std::cout<<"Integral "<<hist->Integral()<<" precision "<<precision<<" NBins "<<NBins<<" threshold "<<threshold<<std::endl;
	if(ndet<8){
	hist->RebinX(1);
	hist->RebinY(1);}
	else{
	hist->RebinX(5);
	hist->RebinY(5);}
	std::cout<<"After rebin: X bins "<<hist->GetNbinsX()<<" Y bins "<<hist->GetNbinsY()<<std::endl;
	//hist->Rebin(1);
	std::vector<double> xpos;
	std::vector<double> ypos;
	std::vector<double> zvalue;
	double dummyx{10.},dummyy{-10};

	for(int ij=1;ij<=NBins;ij++){
		double YSum = YIntegrate(ij,hist); 
		if(YSum<=precision) continue;
		double cumulative=0;
		for(int ip=1;ip<=hist->GetNbinsY();ip++){
			cumulative+=hist->GetBinContent(ij,ip);
			//std::cout<<" ij x "<<ij<<" ip y "<<ip<<" binc "<<hist->GetBinContent(ij,ip)<<" cumulative "<<cumulative<<" YSum "<<YSum<<std::endl;
			if(cumulative>=0.99*YSum){
				dummyx=((TAxis*)hist->GetXaxis())->GetBinCenter(ij);
				dummyy=((TAxis*)hist->GetYaxis())->GetBinCenter(ip);
			//	std::cout<<" binx "<<ij<<" biny "<<ip<<" x "<<dummyx<<" y "<<dummyy<<std::endl;
				xpos.push_back(dummyx);
				ypos.push_back(dummyy);
				break;
			}
		}
	}
	const int Points = xpos.size();
	double xpos1[Points+2];
	double ypos1[Points+2];
	//std::cout<<"Size of x "<<xpos.size()<<" Size of y "<<ypos.size()<<std::endl;
	for(int iz=0;iz<Points;iz++){
		//if(iz==0){xpos1[iz]=xpos[0]-1;ypos1[iz]=0;}
		//else if(iz==(Points-1)){xpos1[iz]=xpos[Points-2-1]+1;ypos1[iz]=0;}
		//else if (iz>=1 && iz<Points-2){
		xpos1[iz+1]=xpos[iz];
		ypos1[iz+1]=ypos[iz];//}
		std::cout<<"Index: x bin "<<iz<<" xpos "<<xpos[iz]<<" ypos "<<ypos[iz]<<std::endl;}
	xpos1[0]=xpos1[1]-1;ypos1[0]=0;
	xpos1[Points+1]=xpos1[Points]+1;ypos1[Points+1]=0;
	std::cout<<"Index: first xpos1 "<<xpos1[0]<<" y0 "<<ypos1[0]<< " last xpos1 "<<xpos1[Points+1]<<"  ylast "<<ypos1[Points+1]<<std::endl;
	gc = new TGraph(Points+2,xpos1,ypos1);
	TH2D *h2 = (TH2D*) hist->Clone();
	
	double factor1 = hist->Integral();
	double factor2 = hist->GetBinContent(hist->GetMaximumBin());
	std::cout<<"Before Scaling Integral "<<factor1<<" Max Bin "<<hist->GetMaximumBin()<<" Maxiumum binc "<<factor2<<std::endl;
	h2->Scale(1./factor2);
	factor1 = h2->Integral();
	factor2 = h2->GetBinContent(hist->GetMaximumBin());
	std::cout<<"After Scaling Integral "<<factor1<<" Max Bin "<<h2->GetMaximumBin()<<" Maxiumum binc "<<factor2<<std::endl;
	c1.Clear();
	c1.cd(1);
	c1.Divide(2,1);
	c1.cd(1);
	gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
        hist->Draw("COLZ");
	hist->GetXaxis()->SetTitle("x (mm)");
	hist->GetYaxis()->SetTitle("y (mm)");
	c1.Update();
	c1.cd(2);
	//std::cout<<"Without entries "<<without_w->GetEntries()<<"  with "<<with_w->GetEntries()<<std::endl;
	
	gPad->SetRightMargin(0.125);
        gPad->SetLeftMargin(0.125);
	//double contours[8]={0.01,0.25,0.5,0.75,0.9,0.95,0.99,1.0};
	double contours[2]={0.001,0.2};
	//h2->SetContour(2, contours);

        //h_nopipe->GetYaxis()->SetRangeUser(2.0e3,1.0e5);
//	gPad->SetLogy();
	//hist->SetLineColor(kRed);
	//hist->SetLineWidth(2);
        h2->Draw("COLZ");
        //h2->Draw("CONT Z LIST");

        //hist->GetXaxis()->SetRangeUser(1,5);
	h2->GetXaxis()->SetTitle("x (mm)");
	h2->GetYaxis()->SetTitle("y (mm)");
	gc->SetLineColor(kRed);
	gc->Draw("SAME L");
	c1.Update();
	//c1.cd(1);
	

   c1.Print(out.c_str(),"pdf");
   //c2->SaveAs("Contour.png");
	}

        c1.Print(Form("%s]",out.c_str()),"pdf");
return 0;
}


double YIntegrate(int binx, TH2D* hist)
{
	double sum=0;
	for(int ij=1;ij<=hist->GetNbinsY();ij++){
		//std::cout<<"binx "<<binx<<" biny "<<ij<<" binc "<<hist->GetBinContent(binx,ij)<<std::endl;
		sum=sum+hist->GetBinContent(binx,ij);
	}
	return sum;
}
