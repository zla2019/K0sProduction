#include <iostream>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include "MyTree.h"
#include "Config.h"
#include <unordered_map>
#include <TCanvas.h>
#include <TLatex.h>
#include <fstream>
#include <ctime>
#include "Hist.h"

const int nPBins = 23;
//PID calibration
const float p_low[nPBins]  = {0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8};
const float p_high[nPBins] = {0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 10.0};
float pidCalib_pion_3p9[nPBins]  = {1.33, 1.50, 1.54, 1.55, 1.55, 1.56, 1.57, 1.58, 1.60, 1.62, 1.63, 1.64, 1.66, 1.67, 1.69, 1.70, 1.72, 1.73, 1.75, 1.78, 1.80, 1.82, 1.84};
float pidCalib_pion_3p5[nPBins]  = {1.30, 1.49, 1.53, 1.54, 1.54, 1.55, 1.56, 1.57, 1.59, 1.61, 1.62, 1.63, 1.65, 1.66, 1.68, 1.69, 1.71, 1.73, 1.75, 1.77, 1.79, 1.82, 1.83};
float pidCalib_pion_3p2[nPBins]  = {1.39, 1.57, 1.61, 1.62, 1.64, 1.65, 1.66, 1.68, 1.69, 1.71, 1.73, 1.74, 1.75, 1.77, 1.79, 1.80, 1.82, 1.84, 1.86, 1.88, 1.91, 1.94, 1.94};
float pidCalib_pion_3p0[nPBins]  = {0.354, 0.354, 0.388, 0.396, 0.404, 0.414, 0.425, 0.437, 0.450, 0.461, 0.472, 0.482, 0.493, 0.499, 0.509, 0.516, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526, 0.526};
float* pidCalib_pion = nullptr;

bool passCut(double a, double lower, double upper);
int getMomBin(float mom, const float* MomLower, const float* MomUpper, const float NMom);
bool loadBadRun(std::unordered_map<int, bool>& badRunList, std::ifstream& ifBadRunList);
bool passCut(double a, Config& config, std::string cutName);
bool passCut(double a, double b, Config& config, std::string cutName);
inline bool passAllCuts(MyTree::Particle& p, Config& config);

int main(int argc, char **argv)
{
	time_t first, second;
	first = time(NULL);
	//initial{{{
	std::string ifListName, ifCutListName, ofName;
	if(argc < 4) {
		return 1;
	}
	ifListName = argv[1];
	ifCutListName = argv[2];
	ofName = argv[3];
	std::cout << "LOG: Reading file list: " << ifListName << std::endl;
	std::cout << "LOG: Reading Cut File: " << ifCutListName << std::endl;
	std::ifstream ifFilelist(ifListName);
	std::ifstream ifCutList(ifCutListName);
	Config config(&ifCutList);
	std::unordered_map<int, bool> badRunList;
	const float Kmass = 0.497611;
	//mass window info
	//const float Mean[3][4] = { { 0.4979, 0.4981, 0.4979, 0.4949 }, { 0.4978, 0.4981, 0.4979, 0.4979 }, { 0.4979, 0.4981, 0.4980, 0.4980 } };
	//const float Sigma[3][4] = { { 0.0027, 0.0035, 0.0033, 0.0033 }, { 0.0027, 0.0034, 0.0036, 0.0036 }, { 0.0027, 0.0035, 0.0039, 0.0039 } };
	const float Mean[3][4] = { { 0.4981, 0.4981, 0.4981, 0.4981 }, { 0.4981, 0.4981, 0.4981, 0.4981 }, { 0.4981, 0.4981, 0.4981, 0.4981 } };
	const float Sigma[3][4] = { { 0.0035, 0.0035, 0.0035, 0.0035 }, { 0.0035, 0.0035, 0.0035, 0.0035 }, { 0.0035, 0.0035, 0.0035, 0.0035 } };
	const float NMassSigma = std::stof(config.mSetList["NSigmaMass"]);
	const bool MuteWarning = true;	//used for debug

	float beamRapidity;
	if(config.mSetList["Energy"] == "3.0") beamRapidity = 1.045;
	else if(config.mSetList["Energy"] == "3.2") beamRapidity = 1.135;
	else if(config.mSetList["Energy"] == "3.5") beamRapidity = 1.24;
	else if(config.mSetList["Energy"] == "3.9") beamRapidity = 1.36;

	if(config.mSetList["Energy"] == "3.0") pidCalib_pion = pidCalib_pion_3p0;
	else if(config.mSetList["Energy"] == "3.2") pidCalib_pion = pidCalib_pion_3p2;
	else if(config.mSetList["Energy"] == "3.5") pidCalib_pion = pidCalib_pion_3p5;
	else if(config.mSetList["Energy"] == "3.9") pidCalib_pion = pidCalib_pion_3p9;

	if(config.mSetList["Energy"] == "3.2") {
		std::ifstream ifBadRun("/home/zla/CF/K0CF/analysis/3p2badRun.list");
		loadBadRun(badRunList, ifBadRun);
	}
	//}}}

	//Save cut info{{{
	TCanvas* caCut = new TCanvas("cuts", "cuts", 1920, 1080);
	caCut->cd();
	TLatex* ltx = new TLatex(0.0, 0.9, "Cuts:");
	//}}}

	//prepare plots{{{
	Hist hist;
	hist.init();
	TFile* ifPurity;
	TH2F* hPurity[3];
	if(config.mSwitchList["OpenPairPurity"]) {
		ifPurity = TFile::Open(config.mSetList["PurityPath"].c_str());
		if(!ifPurity->IsOpen()) {
			std::cout << "ERROR: Purity file is not open" << std::endl;
		}
		for(int icent = 0; icent < 3; ++icent) {
			hPurity[icent] = (TH2F*)ifPurity->Get(Form("hKPurity_cent%i", icent));
		}
	}
	//}}}

	//processing{{{
	std::string ifName;
	while(getline(ifFilelist, ifName)) {
		TFile* ifTree = TFile::Open(ifName.c_str());
		if(!ifTree) {
			std::cout << ifName.c_str() << " can not open" << std::endl;
			continue;
		}
		MyTree *myTree = new MyTree((TTree*)ifTree->Get("tree"));
		if(!myTree->mTree) {
			std::cout << "ERROR: Open tree fail, skipping to next file" << std::endl;
			continue;
		}
		Long64_t nEvent = myTree->getNEvent();

		std::cout << "LOG: Reading file: " << ifName << std::endl;
		std::cout << "LOG: nEvent: " << nEvent << std::endl;
		for(int ievt = 0; ievt < nEvent; ++ievt) {
			if(ievt % 10000 == 0) {
				std::cout << "LOG: ievent: " << ievt << std::endl;
			}
			myTree->getEntry(ievt);
			if(badRunList[myTree->mBufferRunId] == true) {
				if(!MuteWarning) {
					std::cout << "WARNING: Bad run " << myTree->mBufferRunId << std::endl;
				}
				continue;
			}
			float vx = myTree->mBufferVx;
			float vy = myTree->mBufferVy;
			float vz = myTree->mBufferVz;
			float vr = sqrt(vx*vx + vy*vy);
			float cent9 = myTree->mBufferCent9;
			int centForPurity;
			if(cent9 >= 7) {
				centForPurity = 0;
			} else if(cent9 >= 2) {
				centForPurity = 1;
			} else {
				centForPurity = 2;
			}
			unsigned int nK = myTree->mBufferNTrack;
			hist.hRefMult->Fill(myTree->mBufferRefMult);
			if(cent9 < 0 || cent9 > 8) {
				continue;
			}
			hist.hVz->Fill(vz);
			hist.hVr->Fill(vx, vy);
			hist.hCent9->Fill((int)cent9);
			//pre select
			std::vector<int> idxK;
			for(int icurK = 0; icurK < nK; ++icurK) {
				MyTree::Particle curKTmp = myTree->getParticle(icurK, beamRapidity);
				if(!passAllCuts(curKTmp, config)) continue;
				idxK.push_back(icurK);
			}
			int nPassCutK = idxK.size();
			for(int iK1 = 0; iK1 < nPassCutK; ++iK1) {
				int icurK = idxK[iK1];
				MyTree::Particle curK = myTree->getParticle(icurK, beamRapidity);
				TVector3 p(curK.px, curK.py, curK.pz);
				TVector3 pos(curK.bx, curK.by, curK.bz);
				TVector3 vectPmPos = pos - p;
				float cosTheta = vectPmPos.CosTheta();
				int momBinA = getMomBin(curK.pA, p_low, p_high, nPBins);
				int momBinB = getMomBin(curK.pB, p_low, p_high, nPBins);

				float nsigmaA = curK.nSigmaA;
				float nsigmaB = curK.nSigmaB;
				if(config.mSwitchList["OpenNSigmaShift"]) {
					nsigmaA -= pidCalib_pion[momBinA];
					nsigmaB -= pidCalib_pion[momBinB];
				}

				hist.FillAll(curK, (int)cent9);
			}
		}
		delete myTree;
		ifTree->Close();
		delete ifTree;
		myTree = nullptr;
		ifTree = nullptr;
	}
	//}}}

	//prepare output{{{
	TFile* ofPlots = new TFile(ofName.c_str(), "RECREATE");
	hist.Write(ofPlots);
	ofPlots->Close();
	delete ofPlots;
	//}}}
	second = time(NULL);
	std::cout << "processing time: " << difftime(second, first) << std::endl;
	return 0;
}

//my tools{{{
bool loadBadRun(std::unordered_map<int, bool>& badRunList, std::ifstream& ifBadRunList)
{
	if(!ifBadRunList.is_open()) {
		std::cout << "ERROR: Bad run list file not open" << std::endl;
	}
	while(!ifBadRunList.eof()) {
		std::string badRunStr = "";
		ifBadRunList >> badRunStr;
		if(badRunStr == "") {
			break;
		}
		int badRun = std::stoi(badRunStr);
		badRunList[badRun] = true;
	}
	std::cout << "LOG: " << badRunList.size() << " Bad Runs Loaded" << std::endl;
	return true;
}

bool passCut(double a, double lower, double upper)
{
	if(a < lower) {
		return false;
	} else if(a > upper) {
		return false;
	}
	return true;
}

int getMomBin(float mom, const float* MomLower, const float* MomUpper, const float NMom)
{
	int momBin = -1;
	for(int imom = 0; imom < NMom; ++imom) {
		if(mom > MomLower[imom] && mom < MomUpper[imom]) {
			momBin = imom;
			break;
		}
	}
	return momBin;
}

bool passCut(double a, Config& config, std::string cutName)
{
	if(a < config.mCutList[cutName].first) {
		return false;
	} else if(a > config.mCutList[cutName].second) {
		return false;
	}
	return true;
}
bool passCut(double a, double b, Config& config, std::string cutName)
{
	bool hasTHS = false;
	float ths = 0;
	bool pattern = 0;
	if(config.mThsList.find(cutName) == config.mThsList.end()) {
		hasTHS = false;
	} else {
		hasTHS = true;
		ths = config.mThsList[cutName].first;
		pattern = config.mThsList[cutName].second;
	}
	if(hasTHS) {
		bool passTHS = false;
		if(!pattern) {
			passTHS = passCut(b, -1e7, ths);
			return !(passTHS && !passCut(a, config, cutName));
		} else {
			passTHS = passCut(b, ths, 1e7);
			return !(passTHS && !passCut(a, config, cutName));
		}
	} else {
		return passCut(a, config, cutName);
	}

	return false;
}

inline bool passAllCuts(MyTree::Particle& p, Config& config)
{
	TVector3 mom(p.px, p.py, p.pz);
	TVector3 pos(p.bx, p.by, p.bz);
	TVector3 vectPmPos = pos - mom;
	float cosTheta = vectPmPos.CosTheta();
	int momBinA = getMomBin(p.pA, p_low, p_high, nPBins);
	int momBinB = getMomBin(p.pB, p_low, p_high, nPBins);

	float nsigmaA = p.nSigmaA;
	float nsigmaB = p.nSigmaB;
	if(config.mSwitchList["OpenNSigmaShift"]) {
		nsigmaA -= pidCalib_pion[momBinA];
		nsigmaB -= pidCalib_pion[momBinB];
	}

	if(!passCut(p.chi2Topo, config, "Chi2Topo")) return false;
	if(!passCut(p.chi2NDF, config, "Chi2NDF")) return false;
	if(!passCut(p.chi2PrimPip, config, "Chi2PrimPip")) return false;
	if(!passCut(p.chi2PrimPim, config, "Chi2PrimPim")) return false;
	if(!passCut(p.rap, config, "Rap")) return false;
	if(!passCut(p.pt, config, "Pt")) return false;
	if(!passCut(p.nHitsA, config, "NHitsA")) return false;
	if(!passCut(p.nHitsB, config, "NHitsB")) return false;
	if(!passCut(p.dcaA, config, "DCAA")) return false;
	if(!passCut(p.dcaB, config, "DCAB")) return false;
	if(!passCut(p.m2A, p.pA, config, "Mass2Pip")) return false;
	if(!passCut(p.m2B, p.pB, config, "Mass2Pim")) return false;
	if(!passCut(p.dca, config, "KS0DCA")) return false;
	if(!passCut(p.decayLength, config, "DecayLength")) return false;
	if(!passCut(nsigmaA, config, "NSigmaPi")) return false;
	if(!passCut(nsigmaB, config, "NSigmaPi")) return false;
	//if(p.rap > 0 && p.pt < 0.3) return false;
	//if(p.etaA < -2.0 || p.etaA > 0) return false;
	//if(p.etaB < -2.0 || p.etaB > 0) return false;
	return true;
}
//}}}
