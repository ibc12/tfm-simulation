#include "TFile.h"
#include "TEfficiency.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"
#include "TAxis.h"


void angular_distribution()
{
    auto* fin {new TFile {"./graphs.root"}};
    auto* g0 {fin->Get<TGraphErrors>("gGS")};
    auto* g1 {fin->Get<TGraphErrors>("g1")};
    auto* g2 {fin->Get<TGraphErrors>("g2")};

    auto* g0_xs = new TGraphErrors();
    auto* g1_xs = new TGraphErrors();
    auto* g2_xs = new TGraphErrors();

    double gasDensity {2.428e-4}; // g/cm3
    double gasMolarDensity {0.9 * 4.0282 + 0.1 * 58.12}; // g/mol 
    
    //double Nt {(gasDensity/gasMolarDensity) * 6.022e23 * 25.6}; // particles/cm3 * ACTAR length
    double Nt {1.082e21};
    double Nb {(3e3) * 6 * 24 * 3600}; // 3e5 pps 6 days


    auto* fin0 {new TFile {"./Outputs/7.5MeV/transfer_TRIUMF_Eex_0.000_nPS_0_pPS_0.root"}};
    auto* eff0 {fin0->Get<TEfficiency>("eff")};
    auto* fin1 {new TFile {"./Outputs/7.5MeV/transfer_TRIUMF_Eex_0.130_nPS_0_pPS_0.root"}};
    auto* eff1 {fin1->Get<TEfficiency>("eff")};
    auto* fin2 {new TFile {"./Outputs/7.5MeV/transfer_TRIUMF_Eex_0.435_nPS_0_pPS_0.root"}};
    auto* eff2 {fin2->Get<TEfficiency>("eff")};


    std::vector<double> theta_low {20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120.};
    std::vector<double> theta_up  {30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130.};



    //std::vector<double> theta_low {10., 30., 50., 70., 90., 110.};
    //std::vector<double> theta_up {30., 50., 70., 90., 110., 130.};

    std::vector<double> solidAngle;
    std::vector<double> centralInterval;

    for(int i=0; i<theta_low.size(); i++){
        solidAngle.push_back(TMath::Pi() * (TMath::Cos(theta_low[i] * TMath::DegToRad()) - TMath::Cos(theta_up[i] * TMath::DegToRad())));
        centralInterval.push_back(theta_up[i]-(theta_up[i] - theta_low[i])/2);
    }

    TGraph *gOmega = new TGraph(centralInterval.size(), centralInterval.data(), solidAngle.data());


    for(int p = 0; p < g0->GetN(); p++)
    {
        auto bmax {eff0->FindFixBin(theta_up[p])};
        auto bmin {eff0->FindFixBin(theta_low[p])};
        std::vector<double> vals;
        for(int i = bmin; i < bmax; i++)
        {
            vals.push_back(eff0->GetEfficiency(i));
        }
        auto epsilon {TMath::Mean(vals.begin(), vals.end())};
        auto Omega {gOmega->GetPointY(p)}; 


        //auto epsilon {eff0->GetEfficiency(bin)};
        

        double num {g0->GetPointY(p)};
        double denom {Nt * Nb * epsilon * Omega};
 
        auto bin {eff0->FindFixBin(centralInterval[p])};
        auto ueff {( eff0->GetEfficiencyErrorLow(bin) + eff0->GetEfficiencyErrorUp(bin) )/ 2};
        double coeffEps {-num / (Nt * Nb * Omega) / TMath::Power(epsilon, 2)};


        auto xs {num / denom};
        xs *= 1e27;

 
        double coeffNum {TMath::Sqrt(num)/denom};
        double coeffNb {-num / (Nt * Omega * epsilon) / TMath::Power(Nb, 2)};
        double uNb {TMath::Sqrt(Nb)};

        g0_xs->SetPoint(p, centralInterval[p], xs);
        g0_xs->SetPointError(p, 0, TMath::Sqrt(coeffNum * coeffNum + coeffNb * coeffNb * uNb * uNb + coeffEps * coeffEps * ueff * ueff)*1e27);
    }
    for(int p = 0; p < g1->GetN(); p++)
    {
        auto bmax {eff1->FindFixBin(theta_up[p])};
        auto bmin {eff1->FindFixBin(theta_low[p])};
        std::vector<double> vals;
        for(int i = bmin; i < bmax; i++)
        {
            vals.push_back(eff1->GetEfficiency(i));
        }
        auto epsilon {TMath::Mean(vals.begin(), vals.end())};

        //auto epsilon {eff1->GetEfficiency(bin)};
        auto Omega {gOmega->GetPointY(p)}; 

        double num {g1->GetPointY(p)};
        double denom {Nt * Nb * epsilon * Omega};
 

        auto bin {eff1->FindFixBin(centralInterval[p])};
        auto ueff {( eff1->GetEfficiencyErrorLow(bin) + eff1->GetEfficiencyErrorUp(bin) )/ 2};
        double coeffEps {-num / (Nt * Nb * Omega) / TMath::Power(epsilon, 2)};

        auto xs {num / denom};
        xs *= 1e27;


        double coeffNum {TMath::Sqrt(num)/denom};
        double coeffNb {-num / (Nt * Omega * epsilon) / TMath::Power(Nb, 2)};
        double uNb {TMath::Sqrt(Nb)};

        g1_xs->SetPoint(p, centralInterval[p], xs);
        g1_xs->SetPointError(p, 0, TMath::Sqrt(coeffNum * coeffNum + coeffNb * coeffNb * uNb * uNb + coeffEps * coeffEps * ueff * ueff)*1e27);
    }
    for(int p = 0; p < g0->GetN(); p++)
    {
        auto bmax {eff2->FindFixBin(theta_up[p])};
        auto bmin {eff2->FindFixBin(theta_low[p])};
        std::vector<double> vals;
        for(int i = bmin; i < bmax; i++)
        {
            vals.push_back(eff2->GetEfficiency(i));
        }
        auto epsilon {TMath::Mean(vals.begin(), vals.end())};

        //auto bin {eff2->FindFixBin(centralInterval[p])};
        //auto epsilon {eff2->GetEfficiency(bin)};
        auto Omega {gOmega->GetPointY(p)}; 

        double num {g2->GetPointY(p)};
        double denom {Nt * Nb * epsilon * Omega};
        
        auto bin {eff2->FindFixBin(centralInterval[p])};
        auto ueff {( eff2->GetEfficiencyErrorLow(bin) + eff2->GetEfficiencyErrorUp(bin) )/ 2};
        double coeffEps {-num / (Nt * Nb * Omega) / TMath::Power(epsilon, 2)};

        auto xs {num / denom};
        xs *= 1e27;


        double coeffNum {TMath::Sqrt(num)/denom};
        double coeffNb {-num / (Nt * Omega * epsilon) / TMath::Power(Nb, 2)};
        double uNb {TMath::Sqrt(Nb)};

        g2_xs->SetPoint(p, centralInterval[p], xs);
        g2_xs->SetPointError(p, 0, TMath::Sqrt(coeffNum * coeffNum + coeffNb * coeffNb * uNb * uNb + coeffEps * coeffEps * ueff * ueff)*1e27);
    }


    auto* gtheo0 {new TGraphErrors("./Inputs/TheoXS/7.5MeV/angs12nospin.dat", "%lg %lg")};
    auto* gtheo1 {new TGraphErrors("./Inputs/TheoXS/7.5MeV/angp12nospin.dat", "%lg %lg")};
    auto* gtheo2 {new TGraphErrors("./Inputs/TheoXS/7.5MeV/angp32nospin.dat", "%lg %lg")};


    std::vector<int> colors {6, 8, 46};

    g0_xs->SetMarkerColor(colors[0]);
    //g0_xs->SetLineColor(colors[0]);
    gtheo0->SetLineColor(colors[0]);
    gtheo0->SetLineStyle(2);
    gtheo0->SetLineWidth(4);

    g1_xs->SetMarkerColor(colors[1]);
    //g1_xs->SetLineColor(colors[1]);
    gtheo1->SetLineColor(colors[1]);
    gtheo1->SetLineStyle(2);
    gtheo1->SetLineWidth(4);

    g2_xs->SetMarkerColor(colors[2]);
    //g2_xs->SetLineColor(colors[2]);
    gtheo2->SetLineColor(colors[2]);
    gtheo2->SetLineStyle(2);
    gtheo2->SetLineWidth(4);


    //Plotting efficiency

    auto* cEff = new TCanvas("cEff", "Efficiency", 1200, 800);
    cEff->Divide(3, 1);

    cEff->cd(1);
    eff0->Draw("AP");
    eff0->SetTitle("E_{x} = 0 MeV;#theta_{CM} [deg];Eficiencia");

    cEff->cd(2);
    eff1->Draw("AP");
    eff1->SetTitle("E_{x} = 0,130 MeV;#theta_{CM} [deg];Eficiencia");
    
    cEff->cd(3);
    eff2->Draw("AP");
    eff2->SetTitle("E_{x} = 0,435 MeV;#theta_{CM} [deg];Eficiencia");

    //Plotting solid angle
    auto* cSolidAngle = new TCanvas("cOmega", "Solid Angle", 800, 800);
    gOmega->SetMarkerStyle(23);
    gOmega->Draw("AP");

    //Plotting xs
    auto* cXS = new TCanvas("cXS", "XS", 1200, 800);
    cXS->Divide(3, 1);

    cXS->cd(1);
    gPad->SetLogy();
    g0_xs->SetMarkerStyle(24);
    //g0_xs->Scale(1/TMath::TwoPi());
    g0_xs->Draw("AP");
    gtheo0->Draw("same");
    g0_xs->SetTitle("E_{x} = 0 MeV;#theta_{CM} [deg];d#sigma / d#Omega (mb/sr)");

    cXS->cd(2);
    gPad->SetLogy();
    g1_xs->SetMarkerStyle(24);
    //g1_xs->Scale(1/TMath::TwoPi());
    g1_xs->Draw("AP");
    gtheo1->Draw("same");
    g1_xs->SetTitle("E_{x} = 0,130 MeV;#theta_{CM} [deg];d#sigma / d#Omega (mb/sr)");
    
    cXS->cd(3);
    gPad->SetLogy();
    g2_xs->SetMarkerStyle(24);
    //g2_xs->Scale(1/TMath::TwoPi());
    g2_xs->Draw("AP");
    g2_xs->SetTitle("E_{x} = 0,435 MeV;#theta_{CM} [deg];d#sigma / d#Omega (mb/sr)");
    gtheo2->Draw("same");
    

}