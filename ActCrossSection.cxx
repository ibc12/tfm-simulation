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
        fAngleDataGraph.push_back(angle);
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
        fAngleDataGraph.push_back(angle);
        fXSData.push_back(xs);
        fTotalXS += xs * (angle * TMath::DegToRad()) * TMath::Sin(angle * TMath::DegToRad());
    }
    streamer.close();
    //Once is read, we initialize our object
    Init();
}

double ActPhysics::CrossSection::xsIntervalcm(const TString& file, double minAngle, double maxAngle) 
{
    std::ifstream streamer {file};
    if(!streamer)
        throw std::runtime_error("CrossSection::ReadData(): could not open input file named " + file);

    double angle {};
    double xs {};
    double xsIntervalValue {};    
    while(streamer >> angle >> xs)
    {
        // Only process the data if the angle is within the specified range
        if (angle >= minAngle && angle <= maxAngle)
        {
            xsIntervalValue += xs * (angle * TMath::DegToRad()) * TMath::Sin(angle * TMath::DegToRad());
        }
    }
    streamer.close();
    // Once data is read, initialize our object
    return xsIntervalValue * 1e-27;
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
    fCDFGraph = new TSpline3 {"fCDFGraph", &(fCDFData[0]), &(fAngleDataGraph[0]), (int)fCDFData.size(), "b2,e2", 0, 0};
    fCDFGraph->SetTitle("CDF;r;#theta_{CM} [deg]");
}

void ActPhysics::CrossSection::Draw() const
{
    auto* c {new TCanvas};

        // Crear un TGraph para dibujar los ejes
    TGraph* graph = new TGraph(static_cast<int>(fAngleDataGraph.size()), &(fAngleDataGraph[0]), &(fXSData[0]));
    graph->SetTitle(""); // Sin título para el gráfico
    graph->Draw("AP"); // Dibuja el gráfico con puntos y ejes

    // Configurar los ejes del TGraph
    graph->GetXaxis()->SetTitle("Ángulo [deg]"); // Título del eje X
    graph->GetYaxis()->SetTitle("d#sigma / d#Omega [mb sr^{-1}]"); // Título del eje Y

    fCDFGraph->SetLineWidth(2);
    fCDFGraph->SetLineColor(kRed);
    fCDFGraph->SetLineWidth(3);
    fCDFGraph->Draw();
}

double ActPhysics::CrossSection::Sample(const double angle)
{
    return fCDF->Eval(angle);
}

void ActPhysics::CrossSection::Theo()
{
    auto* c1 = new TCanvas;

    // Crear un TGraph para dibujar los ejes
    TGraph* graph = new TGraph(static_cast<int>(fAngleDataGraph.size()), &(fAngleDataGraph[0]), &(fXSData[0]));
    graph->SetTitle("TheoXS"); // Sin título para el gráfico
    graph->Draw("AP"); // Dibuja el gráfico con puntos y ejes

    // Configurar los ejes del TGraph
    graph->GetXaxis()->SetTitle("Ángulo [deg]"); // Título del eje X
    graph->GetYaxis()->SetTitle("d#sigma / d#Omega [mb sr^{-1}]"); // Título del eje Y

    // Crear y dibujar el TSpline3
    fTheoXS = new TSpline3("theoXS", &(fAngleDataGraph[0]), &(fXSData[0]), static_cast<int>(fAngleDataGraph.size()), "b2,e2", 0, 0);
    fTheoXS->SetLineColor(kRed);
    fTheoXS->SetLineWidth(3);
    fTheoXS->Draw("same"); // Dibuja el spline en el mismo lienzo

    c1->Modified(); // Marcar el lienzo como modificado
    c1->Update(); // Volver a actualizar el lienzo para reflejar los cambios

    }