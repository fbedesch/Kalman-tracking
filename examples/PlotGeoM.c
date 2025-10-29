#include <iostream>
#include <TMath.h>
#include <TH1.h>
#include <TAxis.h>
#include <THStack.h>
#include <TString.h>
#include <TLegend.h>
#include "SolGeom.h"
#include <vector>
//
//
void PlotGeoM(TString GEOM)
{
	//
	// Initialize geometry
	//
	std::cout << "Text geometry file: " << GEOM << std::endl;
	char* GeoName = (char*) GEOM.Data();
	SolGeom* G = new SolGeom(GeoName);
	//***************************************
	// Plot geometry         ****************
	//***************************************
	G->Draw();
	//
	//***************************************
	// Plot material         ****************
	//***************************************
	Int_t nDet = G->Ndet();				// Max # detector types
	// Define groups for output
	std::vector<std::vector<TString *>> Groups;
	std::vector<TString>gName;
	// Beam pipe
	gName.push_back("Pipe");
	std::vector<TString *> vPipe;
	vPipe.push_back(new TString("PIPE"));
	Groups.push_back(vPipe);
	// Vertex detector
	gName.push_back("Vertex");
	std::vector<TString *> vVertex;
	vVertex.push_back(new TString("VTX"));
	Groups.push_back(vVertex);
	// Middle detectors
	gName.push_back("Middle");
	std::vector<TString *> vMiddle;
	vMiddle.push_back(new TString("INN_CYL"));
	vMiddle.push_back(new TString("MID1"));
	vMiddle.push_back(new TString("MID2"));
	Groups.push_back(vMiddle);
	// Outer detectors
	gName.push_back("Outer");
	std::vector<TString *> vOuter;
	vOuter.push_back(new TString("OUT1"));
	vOuter.push_back(new TString("OUT2"));
	vOuter.push_back(new TString("OB_SUP"));
	vOuter.push_back(new TString("OB_FLA"));
	Groups.push_back(vOuter);
	//
	// Printout groups
	Int_t ng = (Int_t) Groups.size();
	for(Int_t i=0; i<ng; i++){
		std::cout<<"Group "<<i<<" name is "<<gName[i]<<
		" and includes the following labels:"<<std::endl;
		//
		std::vector<TString *>List = Groups[i];
		Int_t n = (Int_t)List.size();
		for(Int_t j=0; j<n; j++){
			TString Label = *List[j];
			std::cout<<"Label "<<j<<" = "<<Label<<", n = "<<n<<std::endl;
		}
	}
	//
	// Define histograms
	Double_t hSMax= 0.3;
	THStack *hSth = new THStack("hSth", "Radiation length fraction vs #theta");
	hSth->SetMaximum(hSMax);
	TH1D **hX0th = new TH1D*[ng];
	THStack *hScs = new THStack("hScs", "Radiation length fraction vs cos(#theta)");
	hScs->SetMaximum(hSMax);
	TH1D **hX0cs = new TH1D*[ng];
	Int_t stat=0;
	const Int_t Ncol = 5;
	Int_t FillColor[Ncol] = {kYellow, kRed, kAzure, kGreen, kGray};
	for(Int_t i=0; i<ng; i++){
		char hTitleTh[50];
		char hTitleCs[50];
		stat = sprintf(hTitleTh,"Radiation length fraction vs #theta");
		stat = sprintf(hTitleCs,"Radiation length fraction vs cos(#theta)");
		char hNameTh[20];
		char hNameCs[20];
		stat = sprintf(hNameTh,"hTh%d",i);
		stat = sprintf(hNameCs,"hCs%d",i);
		Int_t nBins = 1000;
		hX0th[i] = new TH1D(hNameTh,hTitleTh,nBins,0.,90.);
		if(i<ng)hX0th[i]->SetFillColor(FillColor[i]);
		else hX0th[i]->SetFillColor(FillColor[Ncol-1]);
		hX0cs[i] = new TH1D(hNameCs,hTitleCs,nBins,0.,1.0);
		if(i<ng)hX0cs[i]->SetFillColor(FillColor[i]);
		else hX0cs[i]->SetFillColor(FillColor[Ncol-1]);
	}
	//std::cout<<"Done histogram definition"<<std::endl;
	//
	// Fill histograms
	for(Int_t i=0; i<ng; i++){
		//std::cout<<"Processing group "<<i<<std::endl;
		Int_t Nbt = hX0th[i]->GetNbinsX();
		Int_t Nbc = hX0cs[i]->GetNbinsX();
		std::vector<TString *>List = Groups[i];
		Int_t n = (Int_t)List.size();
		//std::cout<<"Nbt = "<<Nbt<<", Nbc = "<<Nbc<<", n = "<<n<<std::endl;
		// FIll theta
		for(Int_t j=1; j<=Nbt; j++){
			Double_t Th = hX0th[i]->GetBinCenter(j)*TMath::Pi()/180.;
			//std::cout<<"#theta = "<<Th<<std::endl;
			Double_t *Mat = new Double_t[nDet];
			Mat = G->FracX0(Th);
			Double_t x0 = 0;
			for(Int_t k1=0; k1<nDet; k1++){
				for(Int_t k2=0; k2<n; k2++){
					if(G->dType(k1) == *List[k2]){
						x0 += Mat[k1];
						//std::cout<<"dType("<<k1<<") = "<<G->dType(k1)<<
						//", *List["<<k2<<"] = "<<*List[k2]<<std::endl;

					}
				}
			}
			hX0th[i]->SetBinContent(j,x0);
			delete [] Mat;
		}
		//std::cout<<"Done histogram theta fill"<<std::endl;
		// FIll cos(theta)
		for(Int_t j=1; j<Nbc; j++){
			Double_t cs = hX0cs[i]->GetBinCenter(j);
			//std::cout<<"cos(#theta) = "<<cs<<std::endl;
			Double_t Th = TMath::ACos(cs);
			Double_t *Mat = new Double_t[nDet];
			Mat = G->FracX0(Th);
			Double_t x0 = 0;
			for(Int_t k1=0; k1<=nDet; k1++){
				for(Int_t k2=0; k2<n; k2++){
					if(G->dType(k1) == *List[k2])x0 += Mat[k1];
				}
			}
			hX0cs[i]->SetBinContent(j,x0);
			delete [] Mat;
		}
		//std::cout<<"Done histogram cos(#theta) fill"<<std::endl;
	}
	//
	// Stack histograms
	//
	for(Int_t i=0; i<ng; i++){
		hSth->Add(hX0th[i]);
		hScs->Add(hX0cs[i]);
	}
	//
	// Now plot
	TCanvas *cnvX0 = new TCanvas("cnvX0", "Radiation length summary",150,150, 900, 500);
	cnvX0->Divide(2,1);
	//
	TLegend* lgTh = new TLegend(0.2, 0.9, 0.6, 0.70);
	for(Int_t i=0; i<ng;  i++)lgTh->AddEntry(hX0th[i],gName[i],"f");
	TLegend* lgCs = new TLegend(0.2, 0.9, 0.6, 0.70);
	for(Int_t i=0; i<ng;  i++)lgCs->AddEntry(hX0cs[i],gName[i],"f");
	//
	cnvX0->cd(1);
	hSth->Draw();
	lgTh->Draw();
	//
	cnvX0->cd(2);
	hScs->Draw();
	lgCs->Draw();
}

