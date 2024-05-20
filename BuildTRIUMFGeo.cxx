#include "ActGeometry.h"

#include <map>
#include <string>
#include <utility>

void BuildTRIUMFGeo(bool draw = true)
{
    // Define parameters
    // Remember that we work with HALF LENGTHS
    // Drift cage for ACTAR
    double driftX {25.6 / 2}; // cm
    double driftY {25.6 / 2};
    double driftZ {25.6 / 2};
    ActSim::DriftChamber actar(driftX, driftY, driftZ);
    // unit silicon size front type
    double silicon1X {5.0E-2 / 2}; // cm
    double silicon1Y {8. / 2};
    double silicon1Z {5.0 / 2};
    // unit silicon size side type
    double silicon2X {1.5E-1 / 2}; // cm
    double silicon2Y {5. / 2};
    double silicon2Z {5.0 / 2};    
    ActSim::SilUnit silUnit1(0, silicon1X, silicon1Y, silicon1Z);
    ActSim::SilUnit silUnit1_1500(0, silicon2X, silicon1Y, silicon1Z);
    ActSim::SilUnit silUnit2(0, silicon2X, silicon2Y, silicon2Z);
    // set placements for front L0
    std::map<int, std::pair<double, double>> l0Placements {{0, {+2 * silicon1Y, -2 * silicon1Z}},
                                                           {1, {0, -2 * silicon1Z}},
                                                           {2, {-2 * silicon1Y, -2 * silicon1Z}},
                                                           {3, {+1 * silicon1Y + 2.2, 0}},
                                                           {4, {-1 * silicon1Y - 2.2, 0}},
                                                           {5, {+2 * silicon1Y, +2 * silicon1Z}},
                                                           {6, {0, +2 * silicon1Z}},
                                                           {7, {-2 * silicon1Y, +2 * silicon1Z}},
                                                           {8, {+2 * silicon1Y, +4 * silicon1Z}},
                                                           {9, {0, +4 * silicon1Z}},
                                                           {10, {-2 * silicon1Y, +4 * silicon1Z}}};
    ActSim::SilAssembly l0Assembly(0, silUnit1, true, false);
    // offset from flange of ACTAR
    double l0offset {10.4}; // cm
    l0Assembly.SetOffsets(l0offset);
    l0Assembly.SetAssemblyPlacements(l0Placements);
    // L1 assembly! same as L0 but only changing placements for sils 3 and 4
    auto l1Placements {l0Placements};
    l1Placements.at(3) = {-2 * silicon1Y, 0};
    l1Placements.at(4) = {+2 * silicon1Y, 0};
    ActSim::SilAssembly l1Assembly(1, silUnit1, true, false);
    l1Assembly.SetAssemblyPlacements(l1Placements);
    l1Assembly.SetOffsets(l0offset + 2.9);

    // BACKWARDS assembly
    ActSim::SilAssembly backAssembly {2, silUnit1_1500, true, false};
    auto backOffset {-0.5 - (2 * driftX)}; // cm respect to actar rear (beam output)
    backAssembly.SetOffsets(backOffset);
    backAssembly.SetAssemblyPlacements(l0Placements);

    // SIDE assembly left
    std::map<int, std::pair<double, double>> sidePlacements {
        {0, {-3 * silicon2Y, +2 * silicon2Z}}, {1, {-1 * silicon2Y, +2 * silicon2Z}},
        {2, {+1 * silicon2Y, +2 * silicon2Z}}, {3, {+3 * silicon2Y, +2 * silicon2Z}},
        {4, {-3 * silicon2Y, +0 * silicon2Z}}, {5, {-1 * silicon2Y, +0 * silicon2Z}},
        {6, {+1 * silicon2Y, +0 * silicon2Z}}, {7, {+3 * silicon2Y, +0 * silicon2Z}},
        {8, {-3 * silicon2Y, -2 * silicon2Z}}, {9, {-1 * silicon2Y, -2 * silicon2Z}},
        {10, {+1 * silicon2Y, -2 * silicon2Z}}, {11, {+3 * silicon2Y, -2 * silicon2Z}},
    };
    ActSim::SilAssembly sideAssemblyLeft {3, silUnit2, false, true};
    auto sideOffsetRight {10.}; // cm from actar left side
    sideAssemblyLeft.SetOffsets(-1, sideOffsetRight);
    sideAssemblyLeft.SetAssemblyPlacements(sidePlacements);

    // SIDE assembly right

    ActSim::SilAssembly sideAssemblyRight {4, silUnit2, false, true};
    auto sideOffsetLeft {-10. - (2 * driftY)}; // cm from actar right side
    sideAssemblyRight.SetOffsets(-1, sideOffsetLeft);
    sideAssemblyRight.SetAssemblyPlacements(sidePlacements);

    // BUILD GEOMETRY
    ActSim::Geometry geo {};
    geo.SetDrift(actar);
    geo.AddAssemblyData(l0Assembly);
    geo.AddAssemblyData(l1Assembly);
    geo.AddAssemblyData(backAssembly);
    geo.AddAssemblyData(sideAssemblyLeft);
    geo.AddAssemblyData(sideAssemblyRight);    
    geo.Construct();
    geo.Print();

    // SAVE GEO
    std::string path {"./Geometries/"};
    geo.WriteGeometry(path, "triumf");

    // and draw it if necessary
    if(draw)
        geo.Draw();
}
