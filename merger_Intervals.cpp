#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "Rtypes.h"

#include "TAttLine.h"
#include "TGraphErrors.h"
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
void merger_Intervals()
{

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);

    ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");

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
    //double Nt {(gasDensity/gasMolarDensity) * 6.022e23 * 25.6}; // particles/cm3 * ACTAR length
    double Nt {1.082e21};
    double Np {(3e3) * 6 * 24 * 3600}; // 3e5 pps 6 days
    double Nit {1.e6};

    // Set histogram
    int nbins {40};
    double xmin1 {-0.6};
    double xmax1 {1};

    double ymin {0};
    double ymax {40};
    double xmin2 {0};
    double xmax2 {180};
    //Intervals for merge
    double theta_up {30.};
    double theta_low {20.};

    //Create TGraph for Angular Distribution
    auto* g_GS = new TGraphErrors();
    auto* g_1st = new TGraphErrors();
    auto* g_2nd = new TGraphErrors();



    auto* canvas {new TCanvas("canvas", TString::Format("Ex for diferent theta intervals at E_{beam} = %.1fMeV", T1))};
    canvas->DivideSquare(11);

    std::vector<TH1D*> hExs;

    for(int i = 0; i < 11; i++){
        std::vector<TH1D*> hs1;

        auto* hEx {new TH1D {
        TString::Format("hEx %.0f-%.0f", theta_low, theta_up), TString::Format("Exs  #theta_{CM} #in (%.0f , %.0f) [deg];E_{x} [MeV];Counts / %.0f keV", 
        theta_low, theta_up, (xmax1 - xmin1) / nbins * 1000), nbins,xmin1, xmax1}};
        hEx->Sumw2();

        int contador {0};
        for(auto& df : dfs)
        {
            
            double Ex = Exs[contador];
            auto* xs {new ActPhysics::CrossSection()};
            double intervalXS {};
            if(Ex == 0)
            {
                TString data_to_read {TString::Format("./Inputs/TheoXS/%.1fMeV/angs12nospin.dat", T1)};
                xs->ReadData(data_to_read);
                intervalXS = xs->GetTotalXScm();
            }
            else if(Ex == 0.130)
            {
                TString data_to_read {TString::Format("./Inputs/TheoXS/%.1fMeV/angp12nospin.dat", T1)};
                xs->ReadData(data_to_read);
                intervalXS = xs->GetTotalXScm();
            }
            else if(Ex == 0.435)
            {
                TString data_to_read {TString::Format("./Inputs/TheoXS/%.1fMeV/angp32nospin.dat", T1)};
                xs->ReadData(data_to_read);
                intervalXS = xs->GetTotalXScm();
            }

            
            // compute scaling factor
            double scaling {(Nt * Np * intervalXS) / Nit};
            // Get temporary histogram
            auto str {TString::Format("%f <= theta3CM && theta3CM < %f", theta_low, theta_up)};
            // %f indicates that you are formating a float
            // Create a node: a copy of the RDF with custom cuts, new columns, etc
            auto node {df.Filter(str.Data())}; // this filters all the columns,
            // storing in the variable named node only the entries that fullfil the condition
            // You have to use the .Data() method of the TString str bc the Filter method does not allow
            // a TString argument. .Data() converts to a old-C char* type (cousas técnicas, con saber que a próxima vez
            // telo que facer así xa chega :)
            // And now get your histogram as usual (ofc using the node variable)

            auto hInner {node.Histo1D(
                {"h", str, hEx->GetNbinsX(), hEx->GetXaxis()->GetXmin(), hEx->GetXaxis()->GetXmax()}, "Eex")};
            // and push back (best idea is to store in the vector the clone of the temp histogram hInner)
            hInner->Scale(scaling);
            hEx->Add(hInner.GetPtr());
            hs1.push_back((TH1D*)hInner->Clone());

            contador += 1;

            if(contador == 2){
                hExs.push_back(hEx);
            }
        }

        auto* f {new TF1{"f", "[0] * TMath::Voigt(x - [1], [2], [3]) + [4] * TMath::Voigt(x - [5], [6], [7])  + [8] * TMath::Voigt(x - [9], [10], [11]) ", -2, 2}};
        Double_t params[12] = {150, 0, 0.1018, 0.1, 250, 0.13, 0.08895, 0.02, 140, 0.4, 0.09646, 0.08};
        f->SetParameters(params);
        //Fix sigmas
        f->FixParameter(2, 0.106039);
        f->FixParameter(6, 0.0932712);
        f->FixParameter(10, 0.0996123);
        //Fix Ex
        f->FixParameter(1, -3.43436e-03);
        f->FixParameter(5, 1.27182e-01);
        f->FixParameter(9, 4.34972e-01);
        //Fix Gammas
        f->FixParameter(3, 9.77561e-02);
        f->FixParameter(7, 2.36665e-02);
        f->FixParameter(11, 8.71715e-02);
        //Par Lims to amplitudes
        f->SetParLimits(0, 0.1, 350);
        f->SetParLimits(4, 0.1, 350);
        f->SetParLimits(8, 0.1, 350);

        hEx->Fit(f, "0M+I10000");
        
        auto* fGS {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
        double* paramsGS = f->GetParameters();
        fGS->SetParameters(paramsGS);
        auto* f1st {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
        double params1st[4];
        params1st[0] = f->GetParameter(4);
        params1st[1] = f->GetParameter(5);
        params1st[2] = f->GetParameter(6);
        params1st[3] = f->GetParameter(7);
        f1st->SetParameters(params1st);
        auto* f2nd {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
        double params2nd[4];
        params2nd[0] = f->GetParameter(8);
        params2nd[1] = f->GetParameter(9);
        params2nd[2] = f->GetParameter(10);
        params2nd[3] = f->GetParameter(11);
        f2nd->SetParameters(params2nd);

        double centerInterval {theta_up-(theta_up-theta_low)/2};
        g_GS->SetPoint(g_GS->GetN(), centerInterval, fGS->Integral(-2,2) / hEx->GetBinWidth(0));
        g_1st->SetPoint(g_1st->GetN(), centerInterval, f1st->Integral(-2,2) / hEx->GetBinWidth(0));
        g_2nd->SetPoint(g_2nd->GetN(), centerInterval, f2nd->Integral(-2,2) / hEx->GetBinWidth(0));

        g_GS->SetPointError(g_GS->GetN() - 1, 0, TMath::Sqrt(fGS->Integral(-2,2) / hEx->GetBinWidth(0)));
        g_1st->SetPointError(g_1st->GetN() - 1, 0, TMath::Sqrt(f1st->Integral(-2,2) / hEx->GetBinWidth(0)));
        g_2nd->SetPointError(g_2nd->GetN() - 1, 0, TMath::Sqrt(f2nd->Integral(-2,2) / hEx->GetBinWidth(0)));



        std::vector<int> colors {6, 8, 46};
        // plot
        canvas->cd(i+1);
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
        std::vector<std::string> labels {"0", "0.130", "0.435"};
        auto* leg1 {new TLegend {0.2, 0.2}};
        leg1->SetHeader("E_{x} [MeV]");
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        for(int j = 0; j < hs1.size(); j++)
        {
            hs1[j]->SetLineColor(colors[j]);
            hs1[j]->SetLineStyle(kDashed);
            hs1[j]->SetLineWidth(2);
            leg1->AddEntry(hs1[j], labels[j].c_str());
            hs1[j]->Draw("hist same");
        }
        leg1->Draw();

        // add text
        //auto* latex {new TLatex{0.5, 0.5, "#font[42]{#sigma #approx 70 keV}"}};
        //latex->Draw();

        theta_up += 10.;
        theta_low += 10.;
    }

    auto* canvas2 = new TCanvas("canvas2", "Integrales con errores", 1200, 800);
    canvas2->Divide(3, 1);

    canvas2->cd(1);
    gPad->SetLogy();
    g_GS->SetTitle("Integral GS;#theta_{CM} [deg];Integral");
    g_GS->SetMarkerStyle(21);
    g_GS->Draw("AP");

    canvas2->cd(2);
    gPad->SetLogy();
    g_1st->SetTitle("Integral 1st;#theta_{CM} [deg];Integral");
    g_1st->SetMarkerStyle(22);
    g_1st->Draw("AP");

    canvas2->cd(3);
    gPad->SetLogy();
    g_2nd->SetTitle("Integral 2nd;#theta_{CM} [deg];Integral");
    g_2nd->SetMarkerStyle(23);
    g_2nd->Draw("AP");



    auto* fout {new TFile {"./graphs.root", "recreate"}};
    g_GS->Write("gGS");
    g_1st->Write("g1");
    g_2nd->Write("g2");

    fout->Close();
}
