int Contour()
{
	gStyle->SetOptStat(0);
	TLatex latex;
	latex.SetTextSize(0.05);
	latex.SetTextAlign(13);
        TFile f1("./envelope_OP_R4.root");
        //TFile f1("./10MeV_envelope_ep_region1.root");
        //TFile f2("./envelope_ep.root");
        TH2D *hist = new TH2D();//"ha_ee","EE::Radial rate distributins", 700,500,1200);
        TH2D *hist1 = new TH2D();//"ha_ee","EE::Radial rate distributins", 700,500,1200);
        
	string out = "Contour_ep_rw_newgeo.pdf";
	TCanvas c1;
	c1.Print(Form("%s[",out.c_str()),"pdf");
	const int nDet=1;
	string sDet[nDet] = {"Collar1"};
	float zloc = 12250;
	float zloc1[]={-2900,-1449,750,1000,1500,2000,2500,3000,3225,4430.5,5150,5400,5650,6000,6250,6500,6750,7000,7300,7785,8000,9000,10000,11000,11600,12250,15550,18850,20694,22000,25350,26200};

	//string sDet[nDet] = {"tgtPbWall","tgtPbCollar","Col2Ent","USCoil1","USCoil2","USCoil3","USCoil4","USCoil5","Col4Ent","PbWall","DSCoil2","DSCoil3","DSCoil4","DSCoil5","DSCoil6","DSCoil7","DSCoil8","DSCoil9","DSCoil10","Lintel","DSCoil12","DSCoil13","DSCoil14","DSCoil15","DSCoil16","Collar1","Drift","Collar2","GEM4","maindet","piondet","HallExit"};
	for(int ndet=0;ndet<nDet;ndet++){
	//c1.cd(1);
	//c1.Clear();
        hist=(TH2D*)f1.Get(Form("h_%s_electron_d_xy",sDet[ndet].c_str()));
        //hist1=(TH2D*)f2.Get(Form("h_%s_electron_d_xy",sDet[ndet].c_str()));
	TString ftxt = Form("./txtfile/CPhoton_EW_Region4_%s.txt",sDet[ndet].c_str());
        ofstream txtfile;
        txtfile.open(ftxt);

	//hist->Add(hist1);
        //hist=(TH2D*)f1.Get("h_maindet_electron_d_xy");
	std::cout<<"Before rebin: X bins "<<hist->GetNbinsX()<<" Y bins "<<hist->GetNbinsY()<<std::endl;
	if(ndet<8){
	hist->RebinX(2);
	hist->RebinY(2);}
	else{
	hist->RebinX(5);
	hist->RebinY(5);}
	std::cout<<"After rebin: X bins "<<hist->GetNbinsX()<<" Y bins "<<hist->GetNbinsY()<<std::endl;
	TH2D *h2 = (TH2D*) hist->Clone();
	/*int MaxX,MaxY,MaxZ;
	h2->GetMaximumBin(MaxX,MaxY,MaxZ);
	std::cout<<"MaxX "<<MaxX<<" MaxY "<<MaxY<<std::endl;*/
	
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
	h2->SetContour(2, contours);

        //h_nopipe->GetYaxis()->SetRangeUser(2.0e3,1.0e5);
//	gPad->SetLogy();
	//hist->SetLineColor(kRed);
	//hist->SetLineWidth(2);
        //hist->Draw("COLZ");
        h2->Draw("CONT Z LIST");

        //hist->GetXaxis()->SetRangeUser(1,5);
	h2->GetXaxis()->SetTitle("x (mm)");
	h2->GetYaxis()->SetTitle("y (mm)");
	c1.Update();
	//c1.cd(1);
	

	//Get Contours
	
   TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
   TList* contLevel = NULL;
   TGraph* curv     = NULL;
   TGraph* gc       = NULL;

   Int_t nGraphs    = 0;
   Int_t TotalConts = 0;

   if (conts == NULL){
      printf("*** No Contours Were Extracted!\n");
      TotalConts = 0;
      return 0;
   } else {
      TotalConts = conts->GetSize();
   }

   printf("TotalConts = %d\n", TotalConts);

   for(int i = 0; i < TotalConts; i++){
      contLevel = (TList*)conts->At(i);
      printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
      nGraphs += contLevel->GetSize();
   }

   nGraphs = 0;
    //c1.cd(1);
    //c1.Clear();
   //TCanvas* c2 = new TCanvas("c2","Contour List",610,0,600,600);
   //c2->SetTopMargin(0.15);
   //c2->cd(1);
   /*TH2F *hr = new TH2F("hr",
   "#splitline{Negative contours are returned first (highest to lowest). Positive contours are returned from}{lowest to highest. On this plot Negative contours are drawn in red and positive contours in blue.}",
   200, -220, -20, 75, -5, 70);*/

   //hr->Draw();
   Double_t xval0, yval0, zval0;
   TLatex l;
   l.SetTextSize(0.03);
   char val[20];

   for(int i = 0; i < TotalConts; i++){
      contLevel = (TList*)conts->At(i);
      zval0 = contours[i];
      printf("Z-Level Passed in as:  Z = %f\n", zval0);

      // Get first graph from list on curves on this level
      curv = (TGraph*)contLevel->First();
      for(int j = 0; j < contLevel->GetSize(); j++){
      //for(int j = 0; j < 1; j++){
         curv->GetPoint(0, xval0, yval0);
         curv->SetLineColor(kRed);
         nGraphs ++;
         printf("\tGraph: %d  -- %d Elements\n", nGraphs,curv->GetN());

         // Draw clones of the graphs to avoid deletions in case the 1st
         // pad is redrawn.
         gc = (TGraph*)curv->Clone();
	 Double_t xij[4000],yij[4000];
	 if(j==2){
	 txtfile<<curv->GetN()<<std::endl;
	 for(int ijk=0;ijk<curv->GetN();ijk++){
		 int dummmys=curv->GetPoint(ijk,xij[ijk],yij[ijk]);
		 int yint = (int) (yij[ijk]*10);
		 //int yint = (int) (yij[ijk]);
		 if((yint%5)==0){
		 //if((yint%2)==1){
		 std::cout<<"Graph "<<nGraphs<<" entry "<<ijk<<" X "<<xij[ijk]<<" Y "<<yij[ijk]<<"  yint "<<yint<<std::endl;
	 txtfile<<xij[ijk]<<"\t"<<yij[ijk]<<"\t"<<zloc<<std::endl;}}
         txtfile.close();}
         gc->Draw("C");

         //sprintf(val,"%g",zval0);
         //l.DrawLatex(xval0,yval0,val);
         curv = (TGraph*)contLevel->After(curv); // Get Next graph
      }
   }
   c1.Update();
   printf("\n\n\tExtracted %d Contours and %d Graphs \n", TotalConts, nGraphs );
   gStyle->SetTitleW(0.);
   gStyle->SetTitleH(0.);
   c1.Print(out.c_str(),"pdf");
   //c2->SaveAs("Contour.png");
	}

        c1.Print(Form("%s]",out.c_str()),"pdf");
return 0;
}


