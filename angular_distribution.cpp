#include "TFile.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"


void angular_distribution()
{
    auto* fin {new TFile {"./graphs.root"}};
    auto* g0 {fin->Get<TGraphErrors>("gGS")};
    auto* g1 {fin->Get<TGraphErrors>("g1")};
    auto* g2 {fin->Get<TGraphErrors>("g2")};


    double gasDensity {2.428e-4}; // g/cm3
    double gasMolarDensity {0.9 * 4.0282 + 0.1 * 58.12}; // g/mol 
    
    double Nt {(gasDensity/gasMolarDensity) * 6.022e23 * 25.6}; // particles/cm3 * ACTAR length
    double Nb {(3e3) * 6 * 24 * 3600}; // 3e5 pps 6 days


    auto* fin0 {new TFile {"./Outputs/7.5MeV/transfer_TRIUMF_Eex_0.000_nPS_0_pPS_0.root"}};
    auto* eff0 {fin0->Get<TEfficiency>("eff")};
    auto* fin1 {new TFile {"./Outputs/7.5MeV/transfer_TRIUMF_Eex_0.130_nPS_0_pPS_0.root"}};
    auto* eff1 {fin1->Get<TEfficiency>("eff")};
    auto* fin2 {new TFile {"./Outputs/7.5MeV/transfer_TRIUMF_Eex_0.435_nPS_0_pPS_0.root"}};
    auto* eff2 {fin2->Get<TEfficiency>("eff")};


    std::vector<double> theta_low {10., 30., 50., 70., 90., 110.};
    std::vector<double> theta_up {30., 50., 70., 90., 110., 130.};

    std::vector<double> solidAngle;
    std::vector<double> centralInterval;

    for(int i=0; i<theta_low.size(); i++){
        solidAngle.push_back(TMath::TwoPi() * (TMath::Cos(theta_low[i] * TMath::DegToRad()) - TMath::Cos(theta_up[i] * TMath::DegToRad())));
        centralInterval.push_back(theta_up[i]-(theta_up[i] - theta_low[i])/2);
    }

    TGraph *gOmega = new TGraph(centralInterval.size(), centralInterval.data(), solidAngle.data());

    for(int p = 0; p < g0->GetN(); p++)
    {
        auto bin {eff0->FindFixBin(centralInterval[p])};
        auto epsilon {eff0->GetEfficiency(bin)};
        auto Omega {gOmega->GetPointY(p)}; 

        double num {g0->GetPointY(p)};
        double denom {Nt * Nb * epsilon * Omega};
 

        auto xs {num / denom};
        xs *= 1e27;

 

        g0->SetPoint(p, centralInterval[p], xs);
    }
    for(int p = 0; p < g1->GetN(); p++)
    {
        auto bin {eff1->FindFixBin(centralInterval[p])};
        auto epsilon {eff1->GetEfficiency(bin)};
        auto Omega {gOmega->GetPointY(p)}; 

        double num {g1->GetPointY(p)};
        double denom {Nt * Nb * epsilon * Omega};
 

        auto xs {num / denom};
        xs *= 1e27;

 

        g1->SetPoint(p, centralInterval[p], xs);
    }
    for(int p = 0; p < g0->GetN(); p++)
    {
        auto bin {eff2->FindFixBin(centralInterval[p])};
        auto epsilon {eff2->GetEfficiency(bin)};
        auto Omega {gOmega->GetPointY(p)}; 

        double num {g2->GetPointY(p)};
        double denom {Nt * Nb * epsilon * Omega};
 

        auto xs {num / denom};
        xs *= 1e27;

        g2->SetPoint(p, centralInterval[p], xs);
    }



    //Plotting efficiency

    auto* cEff = new TCanvas("cEff", "Efficiency", 1200, 800);
    cEff->Divide(3, 1);

    cEff->cd(1);
    eff0->Draw();

    cEff->cd(2);
    eff1->Draw();
    
    cEff->cd(3);
    eff2->Draw();

    //Plotting solid angle
    auto* cSolidAngle = new TCanvas("cOmega", "Solid Angle", 1200, 800);
    gOmega->SetMarkerStyle(23);
    gOmega->Draw("AP");

    //Plotting xs
    auto* cXS = new TCanvas("cXS", "Efficiency", 1200, 800);
    cXS->Divide(3, 1);

    cXS->cd(1);
    gPad->SetLogy();
    g0->SetMarkerStyle(23);
    g0->Draw("AP");

    cXS->cd(2);
    gPad->SetLogy();
    g0->SetMarkerStyle(23);
    g1->Draw("AP");
    
    cXS->cd(3);
    gPad->SetLogy();
    g0->SetMarkerStyle(23);
    g2->Draw("AP");

}