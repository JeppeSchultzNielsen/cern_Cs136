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

        for(int i = 1; i < entries; i++){
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

    static double gaussLin(double *x, double *par) {
        return par[0] + par[1]*x[0] + par[2]*TMath::Gaus(x[0],par[3],par[4], true);
    }

    static vector<double> findGuesses(double centroidGuess, double fittingWidth, TH1F *histogram){
        auto centroidBin = histogram -> FindBin(centroidGuess);
        auto lowerBin = histogram -> FindBin(centroidGuess - fittingWidth);
        auto upperBin = histogram -> FindBin(centroidGuess + fittingWidth);
        double ymax = 0;
        double xmax = 0;

        int sum = 0;
        int iter = 0;

        for(int k = lowerBin; k < upperBin; k++){
            double y = histogram->GetBinContent(k);
            if(y > ymax){
                ymax = y;
                xmax = histogram->GetBinCenter(k);
            }
            sum += y;
            iter++;
        }

        double avgpoint;
        bool unbroken = true;

        for(int k = histogram->GetXaxis()->FindBin(xmax); unbroken; k++) {
            if (histogram->GetBinContent(k) < sum / iter) {
                avgpoint = (histogram->GetBinCenter(k) - xmax);
                unbroken = false;
            }
            if (k > histogram->GetXaxis()->FindBin(xmax) + 800) {
                avgpoint = fittingWidth;
                unbroken = false;
            }
        }
        return {ymax-histogram->GetBinContent(lowerBin), xmax, histogram->GetBinContent(lowerBin), avgpoint};
    }

    void fit818(string filename, int startTime, int endTime, string saveTo){
        //create a histogram from the startTime to the endTime.
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

        auto lateHist = new TH1F("lateHist","lateHist",12000,0,3000);

        for(int i = 1; i < entries; i++){
            t->GetEntry(i);
            if(i%1000000 == 0){cout << 1.*i/entries << endl;}
            t->GetEntry(i);
            if(timestamp > t0 + endTime){
                break;
            }

            if(timestamp < t0 + startTime){
                continue;
            }

            for(int j = 0; j < 16; j++){
                if(Clov_En[j] > 0.1){
                    lateHist ->Fill(Clov_En[j]);
                }
            }
        }

        double fittingWidth = 3;

        auto guesses = findGuesses(818,fittingWidth,lateHist);
        TF1 *func = new TF1("fit", gaussLin,818-fittingWidth,818+fittingWidth+2,5);
        func ->SetParameters(0,0,3*guesses[0],guesses[1],0.5);
        TFitResultPtr fp = lateHist->Fit("fit","+ && Q && S && E","",818-fittingWidth,818+fittingWidth+2);

        TFile *histFile = TFile::Open(saveTo.c_str(), "RECREATE");
        lateHist -> Write();
        cout << fp ->Parameter(2) << endl;
        cout << fp ->Error(2) << endl;
        histFile -> Close();
        myFile -> Close();
    }
};

int main() {
    auto df = new DecayFinder();
    df -> fit818("data/Run118.root",250000,1000000,"decay/818decay.root");
    //df ->createHistogram("data/Run118.root",500,515,521);
}