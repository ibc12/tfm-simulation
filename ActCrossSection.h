#ifndef ActCrossSection_h
#define ActCrossSection_h

#include "TSpline.h"
#include "TCanvas.h"


#include <string>
#include <vector>

namespace ActPhysics
{
    class CrossSection
    {
    private:
        std::vector<double> fAngleData {};
        std::vector<double> fAngleDataGraph {};
        std::vector<double> fXSData {};
        double fSumXS {};
        double fTotalXS {};
        std::vector<double> fCDFData {};
        TSpline3* fCDF {};
        TSpline3* fCDFGraph{};
        TSpline3* fTheoXS {};

    public:
        void ReadData(const std::string& file);
        void ReadData(const TString& file);
        double xsIntervalcm(const TString& file, double minAngle, double maxAngle);
        void Draw() const;
        double Sample(const double angle);
        void Theo();
        double GetTotalXSmbarn() { return fTotalXS; }
        double GetTotalXScm() const { return fTotalXS * 1e-27; }

    private:
        void Init();
    };
}
#endif
