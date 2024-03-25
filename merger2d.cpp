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

#include "ActCrossSection.cxx"

#include <string>
#include <vector>
void merger2d()
{
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
    double xs {4 * TMath::Pi() *
               1e-27}; // not computed here, should be obtained from XSSampler as the absolute xs -> ~ 1 mbarn

    // Set histogram
    int nbins {200};
    double xmin {0};
    double xmax {180};
    double ymin {0};
    double ymax {40};
    auto* hKin {new TH2F {
        "hKin", "Kinematics for all E_{x};Theta_{Lab} [MeV];E_{Lab}", nbins,
        xmin, xmax, nbins, ymin, ymax}};
    hKin->Sumw2();
    std::vector<TH2F*> hs;
    for(auto& df : dfs)
    {
        int contador {0};
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
        auto h {df.Histo2D(
            {"h", "Kinematics in for loop", hKin->GetNbinsX(), hKin->GetXaxis()->GetXmin(), hKin->GetXaxis()->GetXmax(), 
                hKin->GetNbinsY(), hKin->GetYaxis()->GetXmin(), hKin->GetYaxis()->GetXmax()}, "theta3Lab", "EVertex")};
        h->Scale(scaling);
        hKin->Add(h.GetPtr());
        hs.push_back((TH2F*)h->Clone());

        contador += 1;
    }

    // plot
    auto* c0 {new TCanvas {"c0", "Merger canvas"}};
    //gStyle->SetOptStat(0);
    //hKin->SetLineWidth(2);
    hKin->Draw("colz");
    std::vector<int> colors {6, 8, 46};
    std::vector<std::string> labels {"0", "0.130", "0.435"};
    auto* leg {new TLegend {0.2, 0.2}};
    leg->SetHeader("E_{x} [MeV]");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    for(int i = 0; i < hs.size(); i++)
    {
        //hs[i]->SetLineColor(colors[i]);
        //hs[i]->SetLineStyle(kDashed);
        //hs[i]->SetLineWidth(2);
        leg->AddEntry(hs[i], labels[i].c_str());
        //hs[i]->Draw("hist same");
    }
    //leg->Draw();

    // add text
    //auto* latex {new TLatex{0.5, 0.5, "#font[42]{#sigma #approx 70 keV}"}};
    //latex->Draw();

}