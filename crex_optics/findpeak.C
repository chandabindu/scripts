void findpeak(int run, float adccut)
{
int ij=0;
TChain *T = new TChain("T");
if(run<10000){
T->Add(Form("./rootfiles/prexLHRS_%d_-1.root",run));
TH1F *hp = new TH1F("hp","momentum distributions",200,2.17,2.19);
TCut cut=Form("P.upQadcL>%f",adccut);
T->Project(hp->GetName(),"L.gold.p",cut);
hp->Draw("Q");
int binmax=hp->GetMaximumBin();
double peak = hp->GetXaxis()->GetBinCenter(binmax);
std::cout<<"run "<<run<<"  binmax "<<binmax<<"  peak "<<peak<<std::endl;
}
else{
T->Add(Form("./rootfiles/prexRHRS_%d_-1.root",run));
TH1F *hp = new TH1F("hp","momentum distributions",500,2.17,2.19);
TCut cut=Form("P.upQadcR>%f",adccut);
T->Project(hp->GetName(),"R.gold.p",cut);
int binmax=hp->GetMaximumBin();
double peak = hp->GetXaxis()->GetBinCenter(binmax);
std::cout<<"run "<<run<<"  binmax "<<binmax<<"  peak "<<peak<<std::endl;
}
}
