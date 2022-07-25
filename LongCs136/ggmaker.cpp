#include <iostream>
#include <string>
#include <ctime>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

using namespace std;
using namespace ROOT;

int main() {
    string file = "data/Run157.root";
    unique_ptr<TFile> myFile(TFile::Open(file.c_str()));
    TTree *t = (TTree *) myFile->Get("ids");
    auto entries = t -> GetEntries();
    double_t Clov_En[16];
    Int_t mult;
    Int_t Clov_Time[16];
    t->SetBranchAddress("Time_Clov", &Clov_Time);
    t->SetBranchAddress("Energy_Clov", &Clov_En);
    t->SetBranchAddress("Multiplicity", &mult);


    double currentEn;
    double secondEn;
    string histname = "histogram";
    TH2F *ggHist = new TH2F(histname.c_str(),histname.c_str(),1800,0,1799,1800,0,1799);
    //auto timeHist = new TH1F(histname.c_str(),histname.c_str(),5000,0,4999);
    cout << entries << endl;
    for(int i = 0; i < entries; i++){
        if(i%1000000 == 0){cout << 1.*i/entries << endl;}
        t ->GetEntry(i);
        if(mult < 2) continue;
        for(int j = 0; j < 16; j++){
            currentEn = Clov_En[j];
            if(currentEn > 0.1){
                for(int k = j+1; k < 16; k++){
                    secondEn = Clov_En[k];
                    if(secondEn > 0.1){
                        //don't want coincidences between adjacent clovers; could just be compton events
                        if(j/4 == k/4) continue;
                        if(abs(Clov_Time[j]-Clov_Time[k]) > 50) continue;
                        if(secondEn > currentEn){
                            ggHist -> Fill(secondEn,currentEn);
                        }
                        if(secondEn < currentEn){
                            ggHist -> Fill(currentEn,secondEn);
                        }
                    }
                }
            }
        }
    }
    string filename = "ggCs136.root";
    TFile *histFile = new TFile(filename.c_str(),"RECREATE");
    histFile -> cd();
    ggHist -> Write();
    histFile -> Close();
}