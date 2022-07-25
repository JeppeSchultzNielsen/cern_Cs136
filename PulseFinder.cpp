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

class PulseFinder{
public:
    TH1F *histogram;

    PulseFinder(){};

    void createHistogram(string filename, double peakEnergy, double width){
        unique_ptr<TFile> myFile(TFile::Open(filename.c_str()));
        //TFile *myFile = TFile::Open(filename.c_str());
        //Hent træet
        TTree *t = (TTree *) myFile->Get("ids");
        auto entries = t -> GetEntries();
        double_t Clov_En[16];
        ULong64_t timestamp;
        t->SetBranchAddress("Energy_Clov", &Clov_En);
        t->SetBranchAddress("Timestamp", &timestamp);
        t->GetEntry(1);
        ULong64_t t_i = timestamp;
        t->GetEntry(entries-1);
        ULong64_t t_f = timestamp;

        string name = "pulsefinding";
        histogram = new TH1F(name.c_str(),name.c_str(),(t_f-t_i)/5,t_i,t_f);
        for(int i = 0; i < entries; i++){
            if(i%1000000 == 0){cout << 1.*i/entries << endl;}
            t->GetEntry(i);
            for(int j = 0; j < 16; j++)
                if(abs(peakEnergy-Clov_En[j]) < width){
                    histogram -> Fill(timestamp);
                }
        }
        string filename2 = "pulseHistogram.root";
        TFile *file = new TFile(filename2.c_str(), "RECREATE");
        file->cd();

        histogram -> Write();
        file -> Close();

        myFile -> Close();
    }

    void loadHistogram(){
        //TFile *myFile = TFile::Open("pulseHistogram.root");
        unique_ptr<TFile> myFile(TFile::Open("pulseHistogram.root"));
        histogram = (TH1F *) myFile->Get("pulsefinding");
    }

    vector<int> findPulses(){
        //iterer over all bins. Når der er en stor stigning, er det nok en puls.
        vector<int> timeStamps = {};
        bool recentPulse = false;
        int ticksSinceLastPulse = 0;
        for(int i = 0; i < histogram -> GetNbinsX(); i++){
            ticksSinceLastPulse++;
            if(ticksSinceLastPulse > 1000){ recentPulse = false;}
            if(recentPulse){continue;}
            if(ticksSinceLastPulse < 1000){continue;}

            double prev10avg = 0;
            int k = 0;
            if(i > 10){
                for(int j = 1; j < 31; j++){
                    prev10avg += histogram -> GetBinContent(i-j);
                    k++;
                }
            }
            prev10avg = prev10avg/k;

            k = 0;
            double next10avg = 0;
            for(int j = 0; j < 10; j++){
                next10avg += histogram -> GetBinContent(i+j);
                k++;
            }
            next10avg = next10avg/k;

            if(next10avg > prev10avg + 10){
                recentPulse = true;
                ticksSinceLastPulse = 0;
                timeStamps.push_back(histogram ->GetBinCenter(i));
            }
        }
        TCanvas *c1= new TCanvas;
        histogram -> Draw();
        for(int i = 0; i < timeStamps.size(); i++){
            TLine *l=new TLine(timeStamps[i],0,timeStamps[i],60);
            l->SetLineColor(kBlack);
            l->Draw();
        }
        string filename2 = "canvas.root";
        TFile *file = new TFile(filename2.c_str(), "RECREATE");
        c1 -> Write();
        file -> Close();
        return timeStamps;
    }
};

class SortedPulseFileCreator{
public:
    string filename;
    vector<int> timeStamps;
    SortedPulseFileCreator(string filename){
        this->filename = filename;
        auto puls = new PulseFinder();
        puls -> loadHistogram();
        timeStamps = puls -> findPulses();
    }

    void createSortedPulseFile(string sortedname){
        unique_ptr<TFile> myFile(TFile::Open(filename.c_str()));
        //Hent træet
        TTree *t = (TTree *) myFile->Get("ids");
        auto entries = t -> GetEntries();
        double_t Clov_En[16];
        ULong64_t timestamp;
        double_t energyBeta[4];
        t->SetBranchAddress("Energy_Clov", &Clov_En);
        t->SetBranchAddress("Timestamp", &timestamp);
        t->SetBranchAddress("Energy_VETO",&energyBeta);

        vector<double_t> energySorted;
        vector<double_t> energyBetaSorted;
        ULong64_t timestampSorted;

        TFile *file = TFile::Open(sortedname.c_str(), "RECREATE");
        TTree *sorted = new TTree("Sorted", "Sorted");
        sorted->Branch("energy", &energySorted);
        sorted->Branch("time", &timestampSorted);
        sorted->Branch("beta", &energyBetaSorted);

        int currentTimeStampIndex = 0;
        for(int i = 0; i < entries; i++){
            if(i%1000000 == 0){cout << 1.*i/entries << endl;}
            t->GetEntry(i);
            energySorted = {};
            energyBetaSorted = {};
            //dont want to start in middle of collection: start at first pulse
            if(timestamp < timeStamps[0]) continue;
            if(timestamp > timeStamps[currentTimeStampIndex+1]){currentTimeStampIndex++;}
            for(int j = 0; j < 16; j++){
                if(Clov_En[j] > 0.1){
                    energySorted.push_back(Clov_En[j]);
                }
            }
            for(int j = 0; j < 4; j++){
                if(energyBeta[j] > 0.1){
                    energyBetaSorted.push_back(energyBeta[j]);
                }
            }
            timestampSorted = timestamp - timeStamps[currentTimeStampIndex];
            sorted->Fill();
        }
        sorted -> Write();
        myFile -> Close();
        file -> Close();
    }
};

class DelayedHistogramCreator{
public:
    int t_1;
    int t_2;
    int dt;

    DelayedHistogramCreator(int t1, int t2, int dt) : t_1(t1), t_2(t2), dt(dt) {
    }

    static double gauss2nd(double *x, double *par){
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*TMath::Gaus(x[0],par[4],par[5], true);
    }

    static vector<double> findGuesses(double centroidGuess, double fittingWidth, TH1F *histogram){
        auto centroidBin = histogram -> FindBin(centroidGuess);
        auto lowerBin = histogram -> FindBin(centroidGuess - fittingWidth);
        auto upperBin = histogram -> FindBin(centroidGuess + fittingWidth);
        double ymax = 0;
        double xmax = 0;
        for(int k = lowerBin; k < upperBin; k++){
            double y = histogram->GetBinContent(k);
            if(y > ymax){
                ymax = y;
                xmax = histogram->GetBinCenter(k);
            }
        }
        return {ymax-histogram->GetBinContent(lowerBin), xmax};
    }

    void generateHistograms(string filename, bool betaGate,double fittingWidth){
        //TFile *myFile = TFile::Open(filename.c_str());
        unique_ptr<TFile> myFile(TFile::Open(filename.c_str()));
        //Hent træet

        TTree *t = (TTree *) myFile->Get("Sorted");

        auto entries = t -> GetEntries();
        vector<double_t> *energy = 0;
        vector<double_t> *beta = 0;
        ULong64_t timestamp = 0;

        t->SetBranchAddress("energy", &energy);
        t->SetBranchAddress("beta", &beta);
        t->SetBranchAddress("time", &timestamp);
        string histfilename = "dHistograms/t1:" + to_string(t_1) + "t2:" + to_string(t_2) + "dt:" + to_string(dt) + "betaGate:" + to_string(betaGate) + ".root";
        TFile *histFile = TFile::Open(histfilename.c_str(), "RECREATE");

        auto earlyHist = new TH1F("earlyHist","earlyHist",12000,0,3000);
        auto laterHist = new TH1F("laterHist","laterHist",12000,0,3000);
        auto dHist = new TH1F("dHist","dHist",12000,0,3000);
        for(int i = 1; i < entries; i++){
            t ->GetEntry(i);
            if(betaGate){
                if(beta -> size() == 0){continue;}
            }
            if(timestamp > t_1 && timestamp < t_1 + dt){
                for(int j = 0; j < energy->size(); j++){
                    earlyHist ->Fill(energy -> at(j));
                }
            }

            if(timestamp > t_2 && timestamp < t_2 + dt){
                for(int j = 0; j < energy->size(); j++){
                    laterHist->Fill(energy->at(j));
                }
            }
        }

        string saveto = "dHistograms/t1:" + to_string(t_1) + "t2:" + to_string(t_2) + "dt:" + to_string(dt) + "betaGate:" + to_string(betaGate) + ".txt";
        ofstream mytxt (saveto);

        vector<double> expectedPeaks = {105,340,518,818,1048,1235};
        //histograms are now created, can fit the peaks in each of them.
        for(int i = 0; i < expectedPeaks.size(); i++){
            //first fit on the early histogram
            auto max = findGuesses(expectedPeaks[i],fittingWidth,earlyHist);
            //fit Gauss with 2nd degree background
            TF1 *func = new TF1("fit", gauss2nd,expectedPeaks[i]-fittingWidth,expectedPeaks[i]+fittingWidth,6);
            func ->SetParameters(0,0,0,max[0],max[1],1);
            TFitResultPtr fp = earlyHist->Fit("fit","+ && Q && S","",expectedPeaks[i]-fittingWidth,expectedPeaks[i]+fittingWidth);
            auto counts1 = fp -> Parameter(3);

            //fit also the the laterhist
            max = findGuesses(expectedPeaks[i],fittingWidth,laterHist);
            func = new TF1("fit", gauss2nd,expectedPeaks[i]-fittingWidth,expectedPeaks[i]+fittingWidth,6);
            func ->SetParameters(0,0,0,max[0],max[1],1);
            fp = laterHist->Fit("fit","+ && Q && S","",expectedPeaks[i]-fittingWidth,expectedPeaks[i]+fittingWidth);

            //save results to txt
            mytxt << expectedPeaks[i] << "\t" << counts1 << "\t" << fp -> Parameter(3) << "\n";
        }

        mytxt.close();
        //dHist -> Add(earlyHist,laterHist, 1, -1);
        histFile -> cd();
        earlyHist -> Write();
        laterHist -> Write();
        //dHist ->Write();
        myFile -> Close();
        histFile -> Close();
    }


};

int main() {
    /*auto puls = new PulseFinder();
    puls -> loadHistogram();
    puls -> findPulses();
    auto sortedCreator = new SortedPulseFileCreator("data/Run103.root");
    sortedCreator ->createSortedPulseFile("sorted.root");*/

    auto histCreator = new DelayedHistogramCreator(3000,19000,15000);
    histCreator ->generateHistograms("sorted.root",false,10);
    histCreator = new DelayedHistogramCreator(3000,34000,5000);
    histCreator ->generateHistograms("sorted.root",false,10);
    return 0;
}
