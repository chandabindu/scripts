//This program is written to get the particle rates at different detectors planes
//using data from a 2D histogram in a root file
//
//Chandan Ghosh dec 26 2020
void RingRates(){

TFile *fin = new TFile("/w/halla-scifs17exp/parity/disk1/chandan/gitdir/rad_analysis/script/DSConfig0_ep_disk.root");
TH2F *h1 = new TH2F("h1","h1",2400,-1200,1200,2400,-1200,1200);

const int nRing = 6;
const int nSector = 3;
//h_open = 

int counter=0;
h1=(TH2F*)fin->Get("h_maindet_d_xy");
cout<<"Number of bins X:: "<<h1->GetNbinsX()<<" Y  "<<h1->GetNbinsY()<<" Entries "<<h1->GetEntries()<<"  Integral "<<h1->Integral()<<endl;
for(int ij=1;ij<=h1->GetNbinsX();ij++)
{
for(int pk=1;pk<=h1->GetNbinsY();pk++){
if(h1->GetBinContent(ij,pk)==0) continue;
counter++;
cout<<"Counter  :: "<<counter<<"  Bin x :"<<h1->GetXaxis()->GetBinCenter(ij)<<" Bin Y "<<h1->GetYaxis()->GetBinCenter(pk)<<"  linear bin number  : "<<h1->GetBin(ij,pk)<<"  Bin content "<<h1->GetBinContent(h1->GetBin(ij,pk))<<" From Hist Bin Content "<<h1->GetBinContent(ij,pk)<<endl;
}
}
}
