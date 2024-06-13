#include "TF1.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "ROOT/RResultPtr.hxx"


void fit_Ex(){
    auto file = "./Outputs/7.5MeV/transfer_TRIUMF_Eex_0.000_nPS_0_pPS_0.root";

    auto df = ROOT::RDataFrame {"SimulationTTree", file};

    double xmin {-1};
    double xmax {1};
    int nbins = 100;

    auto hEx {df.Histo1D({"h1", "Ex for gs", nbins, xmin, xmax}, "Eex")};   

    auto* f {new TF1{"f", "[0] * TMath::Voigt(x - [1], [2], [3]) ", -2, 2}};
    f->SetParameters(1300, 0, 0.2, 0.1);
    f->SetParLimits(0, 0, 2000);
    f->SetParLimits(3, 0.05, 0.2);

    hEx->Fit(f, "0M+");


    auto* c {new TCanvas {"c", "c"}};
    hEx->DrawClone("histe");
    f->Draw("same");

}