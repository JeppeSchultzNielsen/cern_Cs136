#include <iostream>
#include <string>
#include <ausa/json/IO.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/StringUtil.h>
#include <ausa/util/FileUtil.h>
#include <ausa/setup/DoubleSidedSiliconDetector.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/eloss/Default.h>
#include <ausa/constants/Mass.h>
#include <ausa/output/OutputConvenience.h>
#include <ctime>
#include <string>
#include <TFile.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <TFitResult.h>
#include <TF1.h>
#include <fstream>
#include <iterator>
#include <sstream>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLine.h>
#include "TROOT.h"
#include "TSystem.h"

using namespace::std;

class DecayFinder{
public:
    DecayFinder(){};

    void createHistogram(string filename, int dt, int lowerBound, int upperBound){
        unique_ptr<TFile> myFile(TFile::Open(filename.c_str()));
        TTree *t = (TTree *) myFile->Get("ids");
        auto entries = t -> GetEntries();
        double_t Clov_En[16];
        ULong64_t timestamp;
        t->SetBranchAddress("Energy_Clov", &Clov_En);
        t->SetBranchAddress("Timestamp", &timestamp);
        t->GetEntry(1);
        ULong64_t t0 = timestamp;
        ULong64_t currentTime = timestamp;

        string saveto = "decay/dt:"+ to_string(dt)+"lower"+ to_string(lowerBound)+"upper"+ to_string(upperBound)+".txt";
        ofstream mytxt (saveto);

        int sum = 0;
        int factor = 1;

        for(int i = 1; i < 10000000; i++){
            t->GetEntry(i);
            if(i%1000000 == 0){cout << 1.*i/entries << endl;}
            if(timestamp >= t0 + factor*dt){
                t->GetEntry(i);
                //write time and sum to file
                mytxt << timestamp-t0 << "\t" <<  sum << endl;
                //reset sum
                sum = 0;
                factor++;
            }

            for(int j = 0; j < 16; j++){
                if(Clov_En[j] > lowerBound){
                    if(Clov_En[j] < upperBound){
                        sum++;
                    }
                }
            }
        }
        mytxt.close();
    }
};

int main() {
    auto df = new DecayFinder();
    df ->createHistogram("data/Run118.root",500,515,521);
}