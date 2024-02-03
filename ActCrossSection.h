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
        std::vector<double> fXSData {};
        double fSumXS {};
        std::vector<double> fCDFData {};
        TSpline3* fCDF {};
        TSpline3* fTheoXS {};

    public:
        void ReadData(const std::string& file);
        void ReadData(const TString& file);
        void Draw() const;
        double Sample(const double angle);
        void Theo();

    private:
        void Init();
    };
}
#endif
