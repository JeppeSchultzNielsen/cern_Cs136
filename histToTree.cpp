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
#include "TH2.h"

using namespace::std;

class HistToTree{
public:
    HistToTree(){};

    void create1DHist(string filename, string histname, int lower, int higher){
        //open the file
        TFile *myFile = TFile::Open(filename.c_str());
        //retrieve the histogram
        string clover_toload = histname;
        int point1;
        int point2;
        int sum;
        TH2F *histogram = (TH2F*)myFile->Get(histname.c_str());
        string hist1dname = "histUpper" + to_string(higher) + "lower" + to_string(lower) + ".root";
        TH1F *histogram1d = new TH1F(hist1dname.c_str(),hist1dname.c_str(),1600,1,1600);
        for(int i = lower; i < higher+1; i++){
            for(int j = 0; j < 1600; j++){
                point1 = histogram ->GetBinContent(i,j);
                point2 = histogram ->GetBinContent(j,i);
                sum = point1 + point2;
                cout << sum << endl;
                cout << j << endl;
                histogram1d ->SetBinContent(j,sum);
            }
        }
        string hist1dfilename = "1dhists/"+hist1dname;
        TFile *file = new TFile(hist1dfilename.c_str(), "RECREATE");
        file->cd();
        histogram1d -> Write();
        file -> Close();
        myFile -> Close();
    }
};

int main() {
    auto hTT = new HistToTree();
    hTT ->create1DHist("data/fit_gg.root","h1",1320,1324);
}