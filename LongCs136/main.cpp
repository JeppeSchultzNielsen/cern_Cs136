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
    TH2F *energy_timeHist = new TH2F(histname.c_str(),histname.c_str(),1600,0,1599,500,0,499);
    //auto timeHist = new TH1F(histname.c_str(),histname.c_str(),5000,0,4999);
    for(int i = 0; i < 100000000; i++){
        t ->GetEntry(i);
        if(mult < 2) continue;
        for(int j = 0; j < 16; j++){
            currentEn = Clov_En[j];
            if(currentEn > 0.1){
                for(int k = j+1; k < 16; k++){
                    secondEn = Clov_En[k];
                    if(secondEn > 0.1){
                        if( (secondEn > 810 && secondEn < 825)){
                            energy_timeHist -> Fill(currentEn,abs(Clov_Time[j]-Clov_Time[k]));
                        }
                        if((currentEn > 810 && currentEn < 825)){
                            energy_timeHist -> Fill(secondEn,abs(Clov_Time[j]-Clov_Time[k]));
                        }
                    }
                }
            }
        }
    }
    string filename = "time.root";
    TFile *histFile = new TFile(filename.c_str(),"RECREATE");
    histFile -> cd();
    energy_timeHist -> Write();
    histFile -> Close();
}