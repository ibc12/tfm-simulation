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
void merger_fit_picos()
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
        "hEx", TString::Format("Ex para todos os picos;E_{x} [MeV];Contas / %.0f keV", (xmax1 - xmin1) / nbins * 1000), nbins,
        xmin1, xmax1}};
    hEx->Sumw2();
    std::vector<TH1D*> hs1;
    std::vector<TF1*> fits;
    std::vector<double> sigmas;

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

        contador += 1;

        TF1* fit = new TF1("fit", "gaus", xmin1, xmax1);
        h1->Fit(fit, "Q");
        fits.push_back(fit);

        // Get sigma of the fit and store it
        double sigma = fit->GetParameter(2);
        sigmas.push_back(sigma);

        // Print sigma of the fit
        std::cout << "Sigma for Ex = " << Ex << " is " << fit->GetParameter(2) << std::endl;

        delete xs;
    }


    // Calculate and print the mean of the sigmas
    double sigmaSum = 0;
    for (const auto& sigma : sigmas)
    {
        sigmaSum += sigma;
    }
    double sigmaMean = sigmaSum / sigmas.size();
    std::cout << "Mean of the sigmas is " << sigmaMean << std::endl;


    // plot
    std::vector<int> colors {6, 8, 46};
    auto* c0 {new TCanvas {"c0", "Merger canvas 1D"}};
    gStyle->SetOptStat(0);
    hEx->SetLineWidth(2);
    hEx->GetXaxis()->SetRangeUser(-0.8, 1);
    hEx->Draw("hist");

    for(size_t i = 0; i < fits.size(); ++i)
    {
        fits[i]->SetLineColor(colors[i]); // Different color for each fit
        fits[i]->Draw("same");
    }

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
    //auto* latex {new TLatex{0.5, 0.5, "#font[42]{#sigma #approx 70 keV}"}};
    //latex->Draw();

}
