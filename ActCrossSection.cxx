#include "ActCrossSection.h"

#include <string>
#include <fstream>




void ActPhysics::CrossSection::ReadData(const std::string& file)
{
    std::ifstream streamer {file};
    if(!streamer)
        throw std::runtime_error("CrossSection::ReadData(): could not open input file named " + file);
    double angle {};
    double xs {};
    while(streamer >> angle >> xs)
    {
        fAngleData.push_back(angle * TMath::DegToRad());
        fXSData.push_back(xs);
        fTotalXS += xs * (angle * TMath::DegToRad()) * TMath::Sin(angle * TMath::DegToRad());
    }
    streamer.close();
    //Once is read, we initialize our object
    Init();
}

void ActPhysics::CrossSection::ReadData(const TString& file) 
{
    std::ifstream streamer {file};
    if(!streamer)
        throw std::runtime_error("CrossSection::ReadData(): could not open input file named " + file);
    double angle {};
    double xs {};
    while(streamer >> angle >> xs)
    {
        fAngleData.push_back(angle * TMath::DegToRad());
        fXSData.push_back(xs);
        fTotalXS += xs * (angle * TMath::DegToRad()) * TMath::Sin(angle * TMath::DegToRad());
    }
    streamer.close();
    //Once is read, we initialize our object
    Init();
}

void ActPhysics::CrossSection::Init()
{
    // Initialize the object after reading file, it returns the CDF
    // We need the sum of the XS to normalize
    for(const auto& xs : fXSData)
    {
        fSumXS += xs;
    }
    // Now the CDF
    for(int i = 0; i < fXSData.size(); i++ )
    {
        double termCDF {};
        for(int j = 0; j <= i; j++)
        {
            termCDF += fXSData[j];
        }
        termCDF /= fSumXS;
        fCDFData.push_back(termCDF);
    }
    // We need the Spline
    fCDF = new TSpline3 {"fCDF", &(fCDFData[0]), &(fAngleData[0]), (int)fCDFData.size(), "b2,e2", 0, 0};
    fCDF->SetTitle("CDF;r;#theta_{CM} [#circ]");
}

void ActPhysics::CrossSection::Draw() const
{
    auto* c {new TCanvas};
    fCDF->SetLineWidth(2);
    fCDF->SetLineColor(kRed);
    fCDF->Draw();
}

double ActPhysics::CrossSection::Sample(const double angle)
{
    return fCDF->Eval(angle);
}

void ActPhysics::CrossSection::Theo()
{
    auto* c1 {new TCanvas};
    fTheoXS = new TSpline3 {"theoXS", &(fAngleData[0]), &(fXSData[0]), (int)fAngleData.size(), "b2,e2", 0, 0};
    fTheoXS->SetLineColor(kRed);
    fTheoXS->Draw();
}
