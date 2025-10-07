#include <TROOT.h>

void LoadTrk()
{
gROOT->Reset();
//gROOT->ProcessLine(".L $PYTHIA8/lib/libpythia8.so");
//gROOT->ProcessLine(".L libDelphes.so");
gROOT->ProcessLine(".I ~/git/fbedesch/delphes2/delphes/examples/classes");
gROOT->ProcessLine(".L examples/classes/SolGeom.cc+");
gROOT->ProcessLine(".L examples/classes/TrkUtil.cc+");
gROOT->ProcessLine(".L examples/classes/SolTrack.cc+");
gROOT->ProcessLine(".L examples/classes/KalmanCk.cc+");
gROOT->ProcessLine(".L examples/KalmanCheck.c+");
gROOT->ProcessLine(".L examples/KalmanComp.c+");
gROOT->ProcessLine(".L examples/PlotGeo.c+");
}
