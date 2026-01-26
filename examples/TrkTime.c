
#include <iostream>
#include <TVector3.h>
#include <TString.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include "examples/classes/SolTrack.h"
//
//
void TrkTime(TString GEOM, Int_t Ncycle = 1000)
{
//
// Initialize geometry
	std::cout << "Text geometry file: " << GEOM << std::endl;
	char* GeoName = (char*) GEOM.Data();
	SolGeom* G = new SolGeom(GeoName);
	//
	// Initialize track
	Double_t mass = 139.57039e-3;	// Pion
	TVector3 x(0.0,0.0,0.0);
	TVector3 p(0.7,0.7,0.5);
	Bool_t Res = kTRUE;
	Bool_t MS  = kTRUE;
	SolTrack StdTrk(x,p,G);
	SolTrack psKlmTrk(x,p,G);
	SolTrack trKlmTrk(x,p,G);
//
//	Initialize times
//
	TStopwatch *StdTime = new TStopwatch();	// Standard method
	StdTime->Stop();
	TStopwatch *psKlmTime = new TStopwatch();	// Pseudo Kalman
	psKlmTime->Stop();
	TStopwatch *trKlmTime = new TStopwatch();	// True Kalman
	trKlmTime->Stop();
//
// Timing loop
//
	for(Int_t n=0; n<Ncycle; n++){
		// Standard
		StdTime->Start(kFALSE);
		StdTrk.CovCalc(Res, MS);
		StdTime->Stop();
		// Pseudo Kalman
		psKlmTime->Start(kFALSE);;
		psKlmTrk.KalmanCov(Res, MS,mass);
		psKlmTime->Stop();
		// True Kalman
		trKlmTime->Start(kFALSE);;
		trKlmTrk.KalmanCovT(Res, MS,mass);
		trKlmTime->Stop();
	}
//
// Print results:
//
	cout<<"=============================================="<<endl;
	cout<<"Timing of Standard Methods over "<<Ncycle<<endl;
	StdTime->Print();
//
	cout<<"=============================================="<<endl;
	cout<<"Timing of pseudo Kalman over "<<Ncycle<<endl;
	psKlmTime->Print();
//
	cout<<"=============================================="<<endl;
	cout<<"Timing of pseudo Kalman over "<<Ncycle<<endl;
	trKlmTime->Print();
}

