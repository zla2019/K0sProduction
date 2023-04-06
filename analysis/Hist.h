#ifndef HIST_H
#define HIST_H
#include <TTree.h>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include "MyTree.h"

class Hist
{
public:
	Hist() {};
	~Hist() {};
	void init();
	void Fill(MyTree::Particle& p);
	void FillLeft(MyTree::Particle& p);
	void FillRight(MyTree::Particle& p);
	void FillAll(MyTree::Particle& p, int cent9);

	void Fill(float qinv, int cent9, int isSideBand1, int isSideBand2);
	void FillMix(float qinv, int cent9, int isSideBand1, int isSideBand2, float sideBandWeight[]);

	void FillCut(MyTree::Particle& p);

	void Write(TFile* of);

	//hist list
	TH1F* hSameKPDG;
	TH1F* hVz;
	TH2F* hVr;
	TH1F* hCent9;
	TH1F* hRefMult;

	TH1F* hCosTheta;
	TH1F* hDecayLength;
	TH1F* hDgDCA;
	TH1F* hDCA;
	TH2F* hDedx;
	TH2F* hMass2;
	TH2F* hAllKRapPt;
	TH2F* hSameKRapPt;
	TH1F* hSameKMass;
	TH1F* hSameKPhi;
	TH2F* hSameKPipRapPt;
	TH2F* hSameKPimRapPt;
	TH1F* hDaughterPipDCA;
	TH1F* hDaughterPimDCA;
	TH1F* hPipNSigma;
	TH1F* hPimNSigma;
	TH1F* hPipNSigma2;
	TH1F* hPimNSigma2;
	TH3F* hSameKPtRapMass[9];

	//cut plots
	TH1F* hCutChi2Topo;
	TH1F* hCutChi2NDF;
	TH1F* hCutChi2PrimA;
	TH1F* hCutChi2PrimB;
	TH1F* hCutNHitsA;
	TH1F* hCutNHitsB;
	TH1F* hCutDCAA;
	TH1F* hCutDCAB;
	TH1F* hCutDecayLength;
	TH2F* hEtaPt;
};
#endif
