#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "Rtypes.h"

#include "TAttLine.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include "TF1.h"
#include "TVirtualFitter.h"

#include "ActCrossSection.cxx"

#include <string>
#include <vector>
void merger()
{
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);


    ROOT::EnableImplicitMT();

    double T1 {7.5};
    std::vector<double> Exs {0., 0.130, 0.435};

 // Construct the output folder path based on T1 with fixed precision
    std::ostringstream oss;
    oss << "./Outputs/" << std::fixed << std::setprecision(1) << T1 << "MeV/";
    std::string outputPath = oss.str();

    std::vector<std::string> files {
        outputPath + "transfer_TRIUMF_Eex_0.000_nPS_0_pPS_0.root",
        outputPath + "transfer_TRIUMF_Eex_0.130_nPS_0_pPS_0.root",
        outputPath + "transfer_TRIUMF_Eex_0.435_nPS_0_pPS_0.root"
    };

    // Read dfs
    std::vector<ROOT::RDF::RNode> dfs;
    for(const auto& file : files)
    {
        dfs.push_back(ROOT::RDataFrame {"SimulationTTree", file});
    }
    // Compute scaling factors
    double gasDensity {2.428e-4}; // g/cm3
    double gasMolarDensity {0.9 * 4.0282 + 0.1 * 58.12}; // g/mol 
    double Nt {(gasDensity/gasMolarDensity) * 6.022e23 * 25.6}; // particles/cm3 * ACTAR length
    double Np {(3e3) * 6 * 24 * 3600}; // 3e5 pps 6 days
    double Nit {1.e6};

    // Set histogram
    int nbins {200};
    double xmin1 {-1};
    double xmax1 {2};
    auto* hEx {new TH1D {
        "hEx", TString::Format("Ex for all peaks;E_{x} [MeV];Counts / %.0f keV", (xmax1 - xmin1) / nbins * 1000), nbins,
        xmin1, xmax1}};
    hEx->Sumw2();
    double ymin {0};
    double ymax {40};
    double xmin2 {0};
    double xmax2 {180};
    auto* hKin {new TH2D {
        "hKin", "Kinematics for all E_{x};Theta_{Lab} [degree];E_{Lab} [MeV]", nbins,
        xmin2, xmax2, nbins, ymin, ymax}};
    hKin->Sumw2();
    std::vector<TH1D*> hs1;
    std::vector<TH2D*> hs2;

    int contador {0};
    for(auto& df : dfs)
    {
        double Ex = Exs[contador];
        auto* xs {new ActPhysics::CrossSection()};
        if(Ex == 0)
        {
            TString data_to_read {TString::Format("./Inputs/TheoXS/%.1fMeV/angs12nospin.dat", T1)};
            xs->ReadData(data_to_read);
        }
        else if(Ex == 0.130)
        {
            TString data_to_read {TString::Format("./Inputs/TheoXS/%.1fMeV/angp12nospin.dat", T1)};
            xs->ReadData(data_to_read);
        }
        else if(Ex == 0.435)
        {
            TString data_to_read {TString::Format("./Inputs/TheoXS/%.1fMeV/angp32nospin.dat", T1)};
            xs->ReadData(data_to_read);
        }

        double totalXS {xs->GetTotalXScm()};
        // compute scaling factor
        double scaling {(Nt * Np * totalXS) / Nit};
        // Get temporary histogram
        auto h1 {df.Histo1D(
            {"h1", "Ex in for loop", hEx->GetNbinsX(), hEx->GetXaxis()->GetXmin(), hEx->GetXaxis()->GetXmax()}, "Eex")};
        h1->Scale(scaling);
        hEx->Add(h1.GetPtr());
        hs1.push_back((TH1D*)h1->Clone());

        auto h2 {df.Histo2D(
            {"h2", "Kinematics in for loop", hKin->GetNbinsX(), hKin->GetXaxis()->GetXmin(), hKin->GetXaxis()->GetXmax(), 
                hKin->GetNbinsY(), hKin->GetYaxis()->GetXmin(), hKin->GetYaxis()->GetXmax()}, "theta3Lab", "EVertex")};
        h2->Scale(scaling);
        hKin->Add(h2.GetPtr());
        hs2.push_back((TH2D*)h2->Clone());

        contador += 1;

        delete xs;
    }

    auto* f {new TF1{"f", "[0] * TMath::Voigt(x - [1], [2], [3]) + [4] * TMath::Voigt(x - [5], [2], [6])  + [7] * TMath::Voigt(x - [8], [2], [9]) ", -2, 2}};
    f->SetParameters(150, 0, 0.07, 0.1, 250, 0.13, 0.1, 140, 0.4, 0.1);
    //f->FixParameter(1, 0);
    // f->FixParameter(5, 0.130);
    // f->FixParameter(9, 0.435);
    f->SetParLimits(0, 10, 350);
    f->SetParLimits(3, 0.05, 0.5);
    f->SetParLimits(4, 10, 350);
    f->SetParLimits(6, 0.05, 0.5);
    f->SetParLimits(7, 10, 350);
    f->SetParLimits(9, 0.05, 0.5);

    hEx->Fit(f, "0M+I10000");
    
    auto* fGS {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
    double* paramsGS = f->GetParameters();
    fGS->SetParameters(paramsGS);
    auto* f1st {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
    double params1st[4];
    params1st[0] = f->GetParameter(4);
    params1st[1] = f->GetParameter(5);
    params1st[2] = f->GetParameter(2);
    params1st[3] = f->GetParameter(6);
    f1st->SetParameters(params1st);
    auto* f2nd {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
    double params2nd[4];
    params2nd[0] = f->GetParameter(7);
    params2nd[1] = f->GetParameter(8);
    params2nd[2] = f->GetParameter(2);
    params2nd[3] = f->GetParameter(9);
    f2nd->SetParameters(params2nd);

    // plot
    std::vector<int> colors {6, 8, 46};
    auto* c0 {new TCanvas {"c0", "Merger canvas 1D"}};
    gStyle->SetOptStat(0);
    hEx->SetLineWidth(2);
    hEx->Draw("histe");
    f->Draw("same");
    fGS->SetLineColor(colors[0]);
    fGS->Draw("same");
    f1st->SetLineColor(colors[1]);
    f1st->Draw("same");
    f2nd->SetLineColor(colors[2]);
    f2nd->Draw("same");
    //hEx->SaveAs("hit_merger_Exs.root");
    std::vector<std::string> labels {"0", "0.130", "0.435"};
    auto* leg1 {new TLegend {0.2, 0.2}};
    leg1->SetHeader("E_{x} [MeV]");
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    for(int i = 0; i < hs1.size(); i++)
    {
        hs1[i]->SetLineColor(colors[i]);
        hs1[i]->SetLineStyle(kDashed);
        hs1[i]->SetLineWidth(2);
        leg1->AddEntry(hs1[i], labels[i].c_str());
        hs1[i]->Draw("hist same");
    }
    leg1->Draw();

    // add text
    auto* latex {new TLatex{0.5, 0.5, "#font[42]{#sigma #approx 70 keV}"}};
    latex->Draw();

    auto* c1 {new TCanvas {"c1", "Merger canvas 2D"}};
    hKin->Draw("colz");
    auto* leg2 {new TLegend {0.2, 0.2}};
    leg2->SetHeader("E_{x} [MeV]");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    for(int i = 0; i < hs2.size(); i++)
    {
        leg2->AddEntry(hs2[i], labels[i].c_str());
    }
}
