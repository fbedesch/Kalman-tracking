
#include <iostream>
#include <TString.h>
#include <TLegend.h>
#include "examples/classes/KalmanCk.h"
//
//
void KalmanComp(TString GEOM1, TString GEOM2, Bool_t Res=kTRUE, Bool_t MS=kTRUE)
{
	//
	// PLots for B field = 2 T
	//
//
// Initialize geometry 1
	std::cout << "Text geometry file: " << GEOM1 <<", B = 2T"<< std::endl;
	char* GeoName1 = (char*) GEOM1.Data();
	SolGeom* G1 = new SolGeom(GeoName1);
	Double_t B = 2;
	G1->SetBfield(B);
	// Start plotting class
	KalmanCk * K1 = new KalmanCk(G1);
	K1->DrawPtScan(10);	// Draw sample tracks
	K1->SetMode(Res, MS);	// Set tracking mode
	// Fill graphs
	K1->Fill();
	// Display
	K1->Print();
//
// Initialize geometry 2
	std::cout << "Text geometry file: " << GEOM2 <<", B = 2T"<< std::endl;
	char* GeoName2 = (char*) GEOM2.Data();
	SolGeom* G2 = new SolGeom(GeoName2);
	G2->SetBfield(B);
	// Start plotting class
	KalmanCk * K2 = new KalmanCk(G2);
	//K1->DrawPtScan(10);	// Draw sample tracks
	K2->SetMode(Res, MS);	// Set tracking mode
	// Fill graphs
	K2->Fill();
	// Display
	//K2->Print();


	//
	// PLots for B field = 3 T
	//
//
// Initialize geometry 1
	std::cout << "Text geometry file: " << GEOM1 <<", B = 3T"<< std::endl;
	//char* GeoName1 = (char*) GEOM1.Data();
	SolGeom* G1b3 = new SolGeom(GeoName1);
	B = 3;
	G1b3->SetBfield(B);
	// Start plotting class
	KalmanCk * K1b3 = new KalmanCk(G1b3);
	//K1b3->DrawPtScan(10);	// Draw sample tracks
	K1b3->SetMode(Res, MS);	// Set tracking mode
	// Fill graphs
	K1b3->Fill();
	// Display
	K1b3->Print();
//
// Initialize geometry 2
	std::cout << "Text geometry file: " << GEOM2 <<", B = 3T"<< std::endl;
	//char* GeoName2 = (char*) GEOM2.Data();
	SolGeom* G2b3 = new SolGeom(GeoName2);
	G2b3->SetBfield(B);
	// Start plotting class
	KalmanCk * K2b3 = new KalmanCk(G2b3);
	//K1->DrawPtScan(10);	// Draw sample tracks
	K2b3->SetMode(Res, MS);	// Set tracking mode
	// Fill graphs
	K2b3->Fill();
	// Display
	//K2->Print();
	//
	// Publication plots: D Resolutions
	//
	TCanvas *ck = new TCanvas("ck","Kalman method - Pt resolutions",250,250, 900, 500);
	ck->Divide(1,1);
	ck->cd(1); gPad->SetLogy(1);
	//
	// B=2 T plots
	//
	TLegend* lg = new TLegend(0.2, 0.9, 0.6, 0.70); 
	TString lPt1 = "Baseline DCH (B = 2T)";
	TString lPt2 = "Heavier  DCH (B = 2T)";
	//
	// Draw GEOM1 pt resolutions
	//
	Int_t i90 = K1->fNangFa-1;
	Double_t ymax = 1.e-2;
	Int_t iColor = (Int_t) kBlue;
	TString title = "#sigma(p_{t})/p_{t}";
	TString title_pt = "p_{t} GeV/c";
	K1->GrSetup(K1->gk_Pt_Pt[i90], ymax, iColor,title, title_pt);
	lg->AddEntry(K1->gk_Pt_Pt[i90], lPt1, "L");
	K1->gk_Pt_Pt[i90]->SetMinimum(1.e-4);
	K1->gk_Pt_Pt[i90]->Draw("APL");
	//
	// Draw GEOM2 resolutions
	//
	Int_t j90 = K2->fNangFa-1;
	//Double_t ymax = .5e-2;
	iColor = (Int_t) kRed;
	//TString title_pt = "p_{t} GeV/c";
	K2->GrSetup(K2->gk_Pt_Pt[j90], ymax, iColor,title, title_pt);
	lg->AddEntry(K2->gk_Pt_Pt[j90], lPt2, "L");
	K2->gk_Pt_Pt[i90]->SetMinimum(1.e-4);
	K2->gk_Pt_Pt[j90]->Draw("PLSAME");
	//
	//lg->Draw("SAME");
	// Print 100 GeV resolution
	Double_t PtPr = K1->fPt_fixA[K1->fNpt_Fa-2];
	Double_t res1 = K1->gk_Pt_Pt[i90]->Eval(PtPr);
	Double_t res2 = K2->gk_Pt_Pt[j90]->Eval(PtPr);
	cout<<"Pt= "<<PtPr<<", s(pt)/pt = "<<res1
	<<" (STD Geo)"<<res2<<" (Heavy Geo)"<<endl;
	cout<<"% degradation = "<<100*(res2-res1)/res1<<"%"<<endl;


	//
	// B=3 T plots
	//
	//TLegend* lg = new TLegend(0.2, 0.9, 0.6, 0.70);
	TString lPt1b3 = "Baseline DCH (B = 3T)";
	TString lPt2b3 = "Heavier  DCH (B = 3T)";
	//
	// Draw GEOM1 pt resolutions
	//
	Int_t i903 = K1b3->fNangFa-1;
	//Double_t ymax = 1.e-2;
	iColor = (Int_t) kBlack;
	TString title3 = "#sigma(p_{t})/p_{t}";
	TString title3_pt = "p_{t} GeV/c";
	K1b3->GrSetup(K1b3->gk_Pt_Pt[i903], ymax, iColor,title3, title3_pt);
	lg->AddEntry(K1b3->gk_Pt_Pt[i903], lPt1b3, "L");
	K1b3->gk_Pt_Pt[i903]->SetMinimum(1.e-4);
	K1b3->gk_Pt_Pt[i903]->Draw("PLSAME");
	//
	// Draw GEOM2 resolutions
	//
	Int_t j903 = K2b3->fNangFa-1;
	//Double_t ymax = .5e-2;
	iColor = (Int_t) kMagenta;
	//TString title_pt = "p_{t} GeV/c";
	K2b3->GrSetup(K2b3->gk_Pt_Pt[j903], ymax, iColor,title3, title3_pt);
	lg->AddEntry(K2b3->gk_Pt_Pt[j903], lPt2b3, "L");
	K2b3->gk_Pt_Pt[i903]->SetMinimum(1.e-4);
	K2b3->gk_Pt_Pt[j903]->Draw("PLSAME");
	//
	lg->Draw("SAME");
	// Print 100 GeV resolution
	Double_t PtPrb3 = K1b3->fPt_fixA[K1b3->fNpt_Fa-2];
	Double_t res1b3 = K1b3->gk_Pt_Pt[i903]->Eval(PtPr);
	Double_t res2b3 = K2b3->gk_Pt_Pt[j903]->Eval(PtPr);
	cout<<"Pt= "<<PtPr<<", s(pt)/pt = "<<res1b3
	<<" (STD Geo)"<<res2b3<<" (Heavy Geo)"<<endl;
	cout<<"% degradation = "<<100*(res2-res1)/res1<<"%"<<endl;
}

