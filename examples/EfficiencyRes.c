#include <iostream>
#include <TString.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <vector>
#include "SolGeom.h"
#include "SolTrack.h"
//
//
void EfficiencyRes(TString GEOM, Bool_t Res=kTRUE, Bool_t MS=kTRUE)
{
//
// Initialize geometry
	std::cout << "Text geometry file: " << GEOM <<", B = 2T"<< std::endl;
	char* GeoName = (char*) GEOM.Data();
	SolGeom* G = new SolGeom(GeoName);	// Import geometry
	Double_t B = 2;
	G->SetBfield(B);	// Set Magnetic field
	//
	// Setup graphs
	//
	const Int_t nPoints = 1000; // #points/position
	const Int_t nPt = 16;		// # pt locations
	Double_t Pta[nPt] = {0.1, 0.2, 0.4, 0.6, 0.8,
						 1.0, 2.0, 4.0, 6.0, 8.0,
						 10., 20., 40., 60., 80.,100.};
	//
	TGraph *gPt_D; 		std::vector<Double_t>D;
	TGraph *gPt_Phi;	std::vector<Double_t>Phi;
	TGraph *gPt_Pt;		std::vector<Double_t>Pt;
	TGraph *gPt_Z0;		std::vector<Double_t>Z0;
	TGraph *gPt_Theta;	std::vector<Double_t>Theta;
	//
	// Kalman results
	Int_t MinMeas = 6;
	std::vector<Double_t>ptg;
	Int_t np_len = 0;
	TVector3 x(0.0,0.0,0.0);
	Double_t sm = 0.02;
	//
	// Loop on events and fill arrays
	//
	for(Int_t ip=0; ip<nPt; ip++){
		TVector3 p(Pta[ip],0.0,0.0);
		Double_t Cot = p.Z()/p.Pt();
		SolTrack Ktrack(x, p, G);		// Initialize track
		for(Int_t n=0; n<nPoints; n++){
			if(Ktrack.nMeas() >= MinMeas){	// Protect from insufficient hits
				Ktrack.KalmanCov(Res, MS);	// Calculate covariance
				// Fill resolution arrays
				D.push_back(Ktrack.s_D()*1.e6);	// Convert to microns
				Phi.push_back(Ktrack.s_phi0());
				Pt.push_back(Ktrack.s_pt());		// sigma(pt)/pt
				Z0.push_back(Ktrack.s_z0()*1.e6);	// Convert to microns
				Double_t CotSig = Ktrack.s_ct();
				Theta.push_back(CotSig/(1.+Cot*Cot));
				Double_t R = sm*gRandom->Rndm();
				ptg.push_back(Pta[ip]+R);
				np_len++;
			}
		}
	}
	//
	// Fill/configure graphs
	//
	Double_t mSize = 0.5;
	gPt_D 	= new TGraph(np_len,ptg.data(),D.data());
	gPt_D->SetTitle("Transverse impact parameter (#mum)");
	gPt_D->SetMarkerStyle(kFullCircle);
	gPt_D->SetMarkerSize(mSize);
	gPt_D->SetMarkerColor(kRed);
	gPt_D->SetMaximum(50.);
	gPt_D->SetMinimum(1.);
	gPt_Phi = new TGraph(np_len,ptg.data(),Phi.data());
	gPt_Phi->SetTitle("#varphi_{0}");
	gPt_Phi->SetMarkerStyle(kFullCircle);
	gPt_Phi->SetMarkerSize(mSize);
	gPt_Phi->SetMarkerColor(kRed);
	gPt_Phi->SetMaximum(2.0e-2);
	gPt_Phi->SetMinimum(1.0e-5);
	gPt_Pt 	= new TGraph(np_len,ptg.data(),Pt.data());
	gPt_Pt->SetTitle("Transverse momentum (GeV)");
	gPt_Pt->SetMarkerStyle(kFullCircle);
	gPt_Pt->SetMarkerSize(mSize);
	gPt_Pt->SetMarkerColor(kRed);
	gPt_Pt->SetMaximum(1.0e-2);
	gPt_Z0 	= new TGraph(np_len,ptg.data(),Z0.data());
	gPt_Z0->SetTitle("Z_{0} (#mum)");
	gPt_Z0->SetMarkerStyle(kFullCircle);
	gPt_Z0->SetMarkerSize(mSize);
	gPt_Z0->SetMarkerColor(kRed);
	gPt_Z0->SetMaximum(50.);
	gPt_Z0->SetMinimum(1.);
	gPt_Theta = new TGraph(np_len,ptg.data(),Theta.data());
	gPt_Theta->SetTitle("#theta");
	gPt_Theta->SetMarkerStyle(kFullCircle);
	gPt_Theta->SetMarkerSize(mSize);
	gPt_Theta->SetMarkerColor(kRed);
	gPt_Theta->SetMaximum(1.0e-2);
	gPt_Theta->SetMinimum(1.0e-5);
	//
	// Plot
	//
	TCanvas *cnv = new TCanvas("cnv","Track resolutions with efficiency",
							50,50, 900,600);
	cnv->Divide(3,2);
	cnv->cd(1); gPad->SetLogy(1);
	gPt_D->Draw("AP");
	cnv->cd(2); gPad->SetLogy(1);
	gPt_Phi->Draw("AP");
	cnv->cd(3); gPad->SetLogy(1);
	gPt_Pt->Draw("AP");
	cnv->cd(4); gPad->SetLogy(1);
	gPt_Z0->Draw("AP");
	cnv->cd(5); gPad->SetLogy(1);
	gPt_Theta->Draw("AP");
}

