#include <TROOT.h>

void LoadTrk()
{
gROOT->Reset();
gROOT->ProcessLine(".I ~/work/Kalman-tracking/examples/classes");
gROOT->ProcessLine(".L examples/classes/SolGeom.cc+");
gROOT->ProcessLine(".L examples/classes/TrkUtil.cc+");
gROOT->ProcessLine(".L examples/classes/SolTrack.cc+");
gROOT->ProcessLine(".L examples/classes/KalmanCk.cc+");
gROOT->ProcessLine(".L examples/KalmanCheck.c+");
gROOT->ProcessLine(".L examples/KalmanComp.c+");
gROOT->ProcessLine(".L examples/PlotGeo.c+");
gROOT->ProcessLine(".L examples/PlotGeoM.c+");
gROOT->ProcessLine(".L examples/EfficiencyRes.c+");
gROOT->ProcessLine(".L examples/TrkTime.c+");
}
