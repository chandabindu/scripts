/* This program is written to calculate the bin center (of a non-zero bin)
 *  of the radial distributions of the inner photon envelope
 *  Chandan March 2021
 */
void GetMaxRadius()
{
        gStyle->SetOptStat(0);
        TFile f1("envelope_photon_inner.root");
        TH1D *hist = new TH1D();
        TGraph *gc;
        string Dstr = "DSCoil9";
        hist=(TH1D*)f1.Get(Form("h_%s_electron_d_r",Dstr.c_str()));
	Int_t NBin = hist->GetNbinsX();
	double bincenter=-10;
	double binc=-10;
	std::cout<<"# of bins "<<NBin<<std::endl;
	for(int ij=1;ij<NBin;ij++){
		if(hist->GetBinContent(ij)==0){
			bincenter = hist->GetBinCenter(ij-1);
			break;
		}
	}
	std::cout<<"Max. Radius at "<<Dstr.c_str()<<" is "<<bincenter+0.5<<std::endl;
}

