#include <list>
std::list<double> listTemp;
std::list<double> listTemp1;
double YIntegrate(int, TH2D*);//, std::vector<double>*);
double XIntegrate(int , TH2D*, int &);// , std::vector<double> *);
double GetMovingAvg(double binc);
double GetMovingAvg1(double binc);
int CalculateBoundaryWithMovingAverage()
{
	gStyle->SetOptStat(0);
	TLatex latex;
	latex.SetTextSize(0.05);
	latex.SetTextAlign(13);
        TFile f1("./envelope_ee_ew_newMap.root");
        TFile f2("./envelope_ep.root");
        TH2D *hist = new TH2D();//"ha_ee","EE::Radial rate distributins", 700,500,1200);
        TH2D *hist1 = new TH2D();//"ha_ee","EE::Radial rate distributins", 700,500,1200);
	TGraph *gc;
	TGraph *gc1;
	TGraph *gc2;
        
	string out = "MovingAverage_ee_ew_newMap.pdf";
	TCanvas c1;
	c1.Print(Form("%s[",out.c_str()),"pdf");
	const int nDet=32;
	//string sDet[nDet] = {"tgtPbWall","tgtPbCollar","Col2Ent","USCoil1","USCoil2"};
	string sDet[nDet] = {"tgtPbWall","tgtPbCollar","Col2Ent","USCoil1","USCoil2","USCoil3","USCoil4","USCoil5","Col4Ent","PbWall","DSCoil2","DSCoil3","DSCoil4","DSCoil5","DSCoil6","DSCoil7","DSCoil8","DSCoil9","DSCoil10","Lintel","DSCoil12","DSCoil13","DSCoil14","DSCoil15","DSCoil16","Collar1","Drift","Collar2","GEM4","maindet","piondet","HallExit"};
	float zloc1[]={-2900,-1449,750,1000,1500,2000,2500,3000,3225,4430.5,5150,5400,5650,6000,6250,6500,6750,7000,7300,7785,8000,9000,10000,11000,11600,12250,15550,18850,20694,22000,25350,26200};
	float zloc[]={5150};
	//for(int ndet=0;ndet<1;ndet++){
	for(int ndet=0;ndet<1;ndet++){
	//c1.cd(1);
	//c1.Clear();
        //hist=(TH2D*)f1.Get(Form("h_%s_electron_d_xy",sDet[ndet].c_str()));
        //hist1=(TH2D*)f2.Get(Form("h_%s_electron_d_xy",sDet[ndet].c_str()));
	//hist->Scale(1./2399.);
	//hist->Scale(1./2409.);
	string Dstr = "DSCoil2";
        hist=(TH2D*)f1.Get(Form("h_%s_electron_d_xy",Dstr.c_str()));
        hist1=(TH2D*)f2.Get(Form("h_%s_electron_d_xy",Dstr.c_str()));
	TString ftxt = Form("./txtfile/Signal_ee_ew_newgeo_%s.txt",Dstr.c_str());
	ofstream txtfile;
	txtfile.open(ftxt);
	//hist->Add(hist1);
        //hist=(TH2D*)f1.Get("h_maindet_electron_d_xy");
	//std::cout<<"Before rebin: X bins "<<hist->GetNbinsX()<<" Y bins "<<hist->GetNbinsY()<<std::endl;
	
	double precision=0.000005*hist->Integral();
	int NBins = hist->GetNbinsX();
	double threshold = precision/NBins;
	//std::cout<<"Integral "<<hist->Integral()<<" precision "<<precision<<" NBins "<<NBins<<" threshold "<<threshold<<std::endl;
	/*if(ndet<=4){
	hist->RebinX(1);
	hist->RebinY(1);}
	else if(ndet==8){
	hist->RebinX(2);
	hist->RebinY(2);}
	else{
	hist->RebinX(5);
	hist->RebinY(5);}*/
	hist->RebinX(2);
	hist->RebinY(2);
	std::vector<double> xpos;
	std::vector<double> ypos;
	double dummyx{10.},dummyy{-10};
	double minBinC = hist->GetBinContent(hist->GetMinimumBin());
	if(minBinC==0) minBinC = 5.e11;
	std::cout<<"After rebin: X bins "<<hist->GetNbinsX()<<" Y bins "<<hist->GetNbinsY()<<" Min Binc "<<minBinC<<std::endl;

	for(int ij=1;ij<=hist->GetNbinsX();ij++){
		std::vector<double> zvalue;
		double YSum = YIntegrate(ij,hist);//,&zvalue); 
		if(YSum<=precision) continue;
		double cumulative=0;
		for(int ip=1;ip<=hist->GetNbinsY();ip++){
			if((dummyy=((TAxis*)hist->GetYaxis())->GetBinCenter(ip))<0)continue;
			double dumybinc=hist->GetBinContent(ij,ip);
			double flag1 = GetMovingAvg(dumybinc);
			cumulative+=dumybinc;
	//		std::cout<<"x bin "<<ij<<" y bin "<<ip<<" binc "<<dumybinc<<" Moving avg. "<<flag1<<" cumulative "<<cumulative<<" YSum "<<YSum<<" dummy y "<<dummyy<<std::endl;
			//if(cumulative>=1.*YSum){
			if(flag1<minBinC){
				dummyx=((TAxis*)hist->GetXaxis())->GetBinCenter(ij);
				dummyy=((TAxis*)hist->GetYaxis())->GetBinCenter(ip-10);
				//std::cout<<"Y uppor scan binx "<<ij<<" biny "<<ip<<" x "<<dummyx<<" y "<<dummyy<<std::endl;
				xpos.push_back(dummyx);
				ypos.push_back(dummyy);
				break;
			}
		}
		listTemp.clear();
	}
	listTemp.clear();
	std::vector<double> xscan1,xscan2;
	std::vector<double> yscan1,yscan2;
	double xBinMin =5.e6;
	int MaxBin;
	for(int ij=1;ij<=hist->GetNbinsY();ij++){
		double XSum = XIntegrate(ij,hist,MaxBin);//&zvalue); 
		//std::cout<<"Xscan before continue x bin "<<ij<<" xsum "<<XSum<<" precision "<<precision<<std::endl;
		if(XSum<=precision) continue;
		//std::cout<<"Xscan Y bin "<<ij<<" xsum "<<XSum<<" precision "<<precision<<std::endl;
		double cumulative=0;
		xBinMin = 1.e11;
		bool flag_lowerbin = false;
		bool flag_upperbin = false;
		for(int ip=MaxBin;ip<hist->GetNbinsX();ip++){
			//std::cout<<"Xscan Y bin "<<ij<<" X bin "<<ip <<" xsum "<<XSum<<" precision "<<precision<<std::endl;
			double dumybinc=hist->GetBinContent(ip,ij);
			double flag1 = GetMovingAvg(dumybinc);
			cumulative+=dumybinc;
	//		std::cout<<"x bin "<<ij<<" y bin "<<ip<<" binc "<<dumybinc<<" Moving avg. "<<flag1<<" cumulative "<<cumulative<<" YSum "<<YSum<<" dummy y "<<dummyy<<std::endl;
			//if(cumulative>=1.*YSum){
			if(flag1<xBinMin && (!flag_upperbin)){
				dummyx=((TAxis*)hist->GetXaxis())->GetBinCenter(ip-10);
				dummyy=((TAxis*)hist->GetYaxis())->GetBinCenter(ij);
				//std::cout<<"X lower binx "<<ip<<" biny "<<ij<<" x "<<dummyx<<" y "<<dummyy<<" Max bin "<<MaxBin<<std::endl;
				xscan1.push_back(dummyx);
				yscan1.push_back(dummyy);
				flag_upperbin=true;
				//break;
			}
		}

		xBinMin = 5.e10;
		for(int ip=MaxBin;ip>0;ip--){
			//std::cout<<"Xscan Y bin "<<ij<<" X bin "<<ip <<" xsum "<<XSum<<" precision "<<precision<<std::endl;
			double dumybinc=hist->GetBinContent(ip,ij);
			double flag1 = GetMovingAvg1(dumybinc);
			cumulative+=dumybinc;
			if(flag1<xBinMin && (!flag_lowerbin)){
				dummyx=((TAxis*)hist->GetXaxis())->GetBinCenter(ip+10);
				dummyy=((TAxis*)hist->GetYaxis())->GetBinCenter(ij);
				//std::cout<<"X upper binx "<<ip<<" biny "<<ij<<" x "<<dummyx<<" y "<<dummyy<<std::endl;
				xscan2.push_back(dummyx);
				yscan2.push_back(dummyy);
				flag_lowerbin=true;
				//break;
			}
		}
		listTemp.clear();
		listTemp1.clear();
	}
	/*const int Points = xpos.size();
	double xpos1[Points+2];
	double ypos1[Points+2];
	//std::cout<<"Size of x "<<xpos.size()<<" Size of y "<<ypos.size()<<std::endl;
	for(int iz=0;iz<Points;iz++){
		//if(iz==0){xpos1[iz]=xpos[0]-1;ypos1[iz]=0;}
		//else if(iz==(Points-1)){xpos1[iz]=xpos[Points-2-1]+1;ypos1[iz]=0;}
		//else if (iz>=1 && iz<Points-2){
		xpos1[iz+1]=xpos[iz];
		ypos1[iz+1]=ypos[iz];//}
		//std::cout<<"Index: x bin "<<iz<<" xpos "<<xpos[iz]<<" ypos "<<ypos[iz]<<std::endl;}
	}
	xpos1[0]=xpos1[1]-1;ypos1[0]=0;
	xpos1[Points+1]=xpos1[Points]+1;ypos1[Points+1]=0;
	//std::cout<<"Index: first xpos1 "<<xpos1[0]<<" y0 "<<ypos1[0]<< " last xpos1 "<<xpos1[Points+1]<<"  ylast "<<ypos1[Points+1]<<std::endl;
	gc = new TGraph(Points+2,xpos1,ypos1);*/
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
		//std::cout<<"Index: x bin "<<iz<<" xpos "<<xpos[iz]<<" ypos "<<ypos[iz]<<std::endl;}
	}
	//std::cout<<"Index: first xpos1 "<<xpos1[0]<<" y0 "<<ypos1[0]<< " last xpos1 "<<xpos1[Points+1]<<"  ylast "<<ypos1[Points+1]<<std::endl;
	xpos1[0]=xpos[0];
	ypos1[0]=0;
	xpos1[Points+1]=xpos[Points-1];
	ypos1[Points+1]=0;
	gc = new TGraph(Points+2,xpos1,ypos1);
	const int NXscan1 = xscan1.size();
	double XScan1[NXscan1+1];
	double YScan1[NXscan1+1];
	//std::cout<<"Size of x "<<xpos.size()<<" Size of y "<<ypos.size()<<std::endl;
	for(int iz=0;iz<NXscan1;iz++){
		//if(iz==0){xpos1[iz]=xpos[0]-1;ypos1[iz]=0;}
		//else if(iz==(Points-1)){xpos1[iz]=xpos[Points-2-1]+1;ypos1[iz]=0;}
		//else if (iz>=1 && iz<Points-2){
		XScan1[iz+1]=xscan1[iz];
		YScan1[iz+1]=yscan1[iz];//}
		//std::cout<<"Index: x bin "<<iz<<" xpos "<<xpos[iz]<<" ypos "<<ypos[iz]<<std::endl;}
	}
	//std::cout<<"Index: first xpos1 "<<xpos1[0]<<" y0 "<<ypos1[0]<< " last xpos1 "<<xpos1[Points+1]<<"  ylast "<<ypos1[Points+1]<<std::endl;
	XScan1[0]=xscan1[0];
	YScan1[0]=0;
	gc1 = new TGraph(NXscan1+1,XScan1,YScan1);
	const int NXscan2 = xscan2.size();
	double XScan2[NXscan2+1];
	double YScan2[NXscan2+1];
	//std::cout<<"Size of x "<<xpos.size()<<" Size of y "<<ypos.size()<<std::endl;
	for(int iz=0;iz<NXscan2;iz++){
		//if(iz==0){xpos1[iz]=xpos[0]-1;ypos1[iz]=0;}
		//else if(iz==(Points-1)){xpos1[iz]=xpos[Points-2-1]+1;ypos1[iz]=0;}
		//else if (iz>=1 && iz<Points-2){
		XScan2[iz+1]=xscan2[iz];
		YScan2[iz+1]=yscan2[iz];//}
		//std::cout<<"Index: x bin "<<iz<<" xpos "<<xpos[iz]<<" ypos "<<ypos[iz]<<std::endl;}
	}
	XScan2[0]=xscan2[0];
	YScan2[0]=0;
	//std::cout<<"Index: first xpos1 "<<xpos1[0]<<" y0 "<<ypos1[0]<< " last xpos1 "<<xpos1[Points+1]<<"  ylast "<<ypos1[Points+1]<<std::endl;
	//ypos1[0]=0;
	//ypos1[Points-1]=0;
	gc2 = new TGraph(NXscan2+1,XScan2,YScan2);
	txtfile<<NXscan1+NXscan2+2<<std::endl;
	for(int ijk=0;ijk<(NXscan1+1);ijk++)
	{
		txtfile<<XScan1[ijk]<<"\t"<<YScan1[ijk]<<"\t"<<zloc[0]<<std::endl;
	}
	for(int ijk=NXscan2;ijk>=0;ijk--)
	{
		txtfile<<XScan2[ijk]<<"\t"<<YScan2[ijk]<<"\t"<<zloc[0]<<std::endl;
	}
	txtfile.close();
	TH2D *h2 = (TH2D*) hist->Clone();
	
	double factor1 = hist->Integral();
	double factor2 = hist->GetBinContent(hist->GetMaximumBin());
	//std::cout<<"Before Scaling Integral "<<factor1<<" Max Bin "<<hist->GetMaximumBin()<<" Maxiumum binc "<<factor2<<std::endl;
	h2->Scale(1./factor2);
	factor1 = h2->Integral();
	factor2 = h2->GetBinContent(hist->GetMaximumBin());
	//std::cout<<"After Scaling Integral "<<factor1<<" Max Bin "<<h2->GetMaximumBin()<<" Maxiumum binc "<<factor2<<std::endl;
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
	//gc->Draw("SAME L");
	gc1->SetLineWidth(2);
	gc1->SetLineColor(kMagenta);
	gc1->Draw("SAME L");
	gc2->SetLineColor(kOrange);
	gc2->SetLineWidth(2);
	gc2->Draw("SAME L");
	c1.Update();
	//c1.cd(1);
	

   c1.Print(out.c_str(),"pdf");
   //c2->SaveAs("Contour.png");
	}

        c1.Print(Form("%s]",out.c_str()),"pdf");
return 0;
}


double YIntegrate(int binx, TH2D* hist)//, std::vector<double> *zvalue)
{
	double sum=0;
	for(int ij=1;ij<=hist->GetNbinsY();ij++){
		//std::cout<<"binx "<<binx<<" biny "<<ij<<" binc "<<hist->GetBinContent(binx,ij)<<std::endl;
		double dummy = hist->GetBinContent(binx,ij);
		sum=sum+dummy;
		//zvalue->push_back(dummy);
	}
	return sum;
}

double XIntegrate(int biny, TH2D* hist, int &MaxBin)//, std::vector<double> *zvalue)
{
	double sum=0;
	double centroid=0;
	for(int ij=1;ij<=hist->GetNbinsX();ij++){
		//std::cout<<"binx "<<binx<<" biny "<<ij<<" binc "<<hist->GetBinContent(binx,ij)<<std::endl;
		double dummy = hist->GetBinContent(ij,biny);
		centroid += ij*hist->GetBinContent(ij,biny);
		sum=sum+dummy;
	}
	double ratio =centroid/sum;
	MaxBin =  (int) ratio;
	//std::cout<<"Inside X integral centroid "<<centroid<<" sum "<<sum<<" MaxBin "<<MaxBin<<" ratio "<<ratio<<std::endl;
	return sum;
}

double GetMovingAvg(double binc){
	listTemp.push_back(binc);
	int avgSize = 10;
	if(listTemp.size()>avgSize) listTemp.pop_front();
        double sum = 0;
	for(std::list<double>::iterator p = listTemp.begin(); p!=listTemp.end();++p){
	sum += (double)*p;
	}
	return 	sum/(1.*avgSize);
}
double GetMovingAvg1(double binc){
	listTemp1.push_back(binc);
	int avgSize = 10;
	if(listTemp1.size()>avgSize) listTemp1.pop_front();
        double sum = 0;
	for(std::list<double>::iterator p = listTemp1.begin(); p!=listTemp1.end();++p){
	sum += (double)*p;
	}
	return 	sum/(1.*avgSize);
}
