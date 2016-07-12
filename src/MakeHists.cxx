#include "MakeHists.hpp"
#include "Functions.hpp"
#include "TMath.h"

using namespace std;

bool MakeHists::initialize(ConfigParser *config, DSHandler *ds) {
  mConfig = config;
  hs = new HistStore();
  calculator = new VariableCalculator();
  InitYields(ds);
  mTRFvariables.push_back("Mbb_MindR");
  mTRFvariables.push_back("dRbb_MaxPt");
  mTRFvariables.push_back("dRbb_MaxM");

  mVarToCalc.push_back("pT_jet1");
  mVarToCalc.push_back("pT_jet2");
  mVarToCalc.push_back("eta_jet1");
  mVarToCalc.push_back("eta_jet2");
  mVarToCalc.push_back("phi_jet1");
  mVarToCalc.push_back("phi_jet2");
  mVarToCalc.push_back("pT_bJet1");
  mVarToCalc.push_back("eta_bJet1");
  mVarToCalc.push_back("phi_bJet1");
  mVarToCalc.push_back("pT_lep1");
  mVarToCalc.push_back("pT_lep2");
  mVarToCalc.push_back("eta_lep1");
  mVarToCalc.push_back("eta_lep2");
  mVarToCalc.push_back("phi_lep1");
  mVarToCalc.push_back("phi_lep2");

  mSysName.clear();
  if (config->HasCommonSetting("TreeSysName"))
    mSysName = config->GetCommonSetting("TreeSysName");
  mSysName.push_back("nominal");
  return true;
}

bool MakeHists::run(TTree *event, map<string, float> weights,
                    TTreeFormulaContainer *formulas,
                    std::map<string, bool> bControl) {
  mEvent = event;
  mSysWeights.clear();
  std::vector<string> mRegions;
  mRegions.clear();
  string tmpTreeName = mEvent->GetName();
  if (tmpTreeName == "nominal" || tmpTreeName == "nominal_Loose")
    isNominal = true;
  else
    isNominal = false;
  string mSample;
  bool isTRF = bControl.at("isTRF");
  bool doTTbarMerge = bControl.at("doTTbarMerge");
  bool doSysmetic = bControl.at("doSysmetic");

  calculator->CalculateVariables(mEvent);
  int mcChannel =
      *(Tools::Instance().GetTreeValue<int>(mEvent, "mcChannelNumber"));
  mSample = Tools::Instance().GetSampleType(mcChannel);
  if (doSysmetic && isNominal) {
    GetSysWeights(mEvent, mSysWeights);
    for (auto w : weights) {
      if (w.first == "norm" || w.first == "weight_NNLO" ||
          w.first == "ttbb_Nominal_weight")
        continue;
      string tmpName = w.first.substr(0, w.first.find("_weight"));
      mSysWeights[tmpName] = w.second;
    }
  }
  mSysWeights["nominal"] = 1.0;

  // Get Weights!
  float tmpWeights = Tools::Instance().GetWeight(mcChannel);
  float weight_mc = 1.0;
  float weight_pileup = 1.0;
  float weight_jvt = 1.0;
  float weight_leptonSF = 1.0;
  float weight_bTagSF_77 = 1.0;

  // No weight for DATA!
  if (mcChannel != 0) {
    if (mSample == "ttbar") {
      mSample = doHeavyFlavor(mEvent);
    }

    weight_mc = *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_mc"));
    weight_pileup =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_pileup"));
    weight_jvt = *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_jvt"));
    weight_leptonSF =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_leptonSF"));
    weight_bTagSF_77 =
        *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_bTagSF_77"));

    if (mSample == "ttlight" || mSample == "ttcc" || mSample == "ttbb")
      tmpWeights = tmpWeights * weight_mc * weight_bTagSF_77 * weight_leptonSF *
                   weight_jvt * weight_pileup * weights["ttbb_Nominal_weight"] *
                   weights["weight_NNLO"];
    else
      tmpWeights = tmpWeights * weight_mc * weight_bTagSF_77 * weight_leptonSF *
                   weight_jvt * weight_pileup;
  }

  for (auto w : mSysName) {
    if (!isNominal && w != "nominal")
      continue;

    float weight = mSysWeights.at(w);
    float mWeights = tmpWeights;

    if (mcChannel == 0)
      isTRF = false;

    // Fakes!
    // bool isFake = false;
    if (mSample == "DATA" && w != "nominal")
      continue;
    if ((mSample != "ttbb" && mSample != "ttcc" && mSample != "ttlight") &&
        (w.find("ttbb") != string::npos || w.find("NNLO") != string::npos))
      continue;

    // Treat sys weight for MC!
    if (mcChannel != 0) {
      if (mSample == "ttbb" && (w == "NNLO_topPtUp" || w == "NNLO_ttbarPtUP"))
        continue;

      if (w.find("leptonSF") != string::npos)
        mWeights = mWeights * weight / weight_leptonSF;
      else if (w.find("bTagSF_77") != string::npos)
        mWeights = mWeights * weight / weight_bTagSF_77;
      else if (w.find("pileup") != string::npos)
        mWeights = mWeights * weight / weight_pileup;
      else if (w.find("jvt") != string::npos)
        mWeights = mWeights * weight / weight_jvt;

      if (mSample == "ttlight" || mSample == "ttcc") {
        if (w.find("ttbb") != string::npos)
          mWeights = mWeights * weight / weights["ttbb_Nominal_weight"];
        else if (w.find("NNLO") != string::npos)
          mWeights = mWeights * weight / weights["weight_NNLO"];
      } else if (mSample == "ttbb") {
        if (w.find("ttbb") != string::npos)
          mWeights = mWeights * weight / weights["ttbb_Nominal_weight"];
      }
    }

    // Selectioins!
    for (int ientry = 0; ientry < formulas->GetEntries(); ++ientry) {
      TTreeFormula *formula = formulas->GetFormula(ientry);
      if (formula->EvalInstance()) {
        mRegions.push_back(formula->GetName());
      }
    }
    if (mRegions.empty())
      return false;

    std::map<string, float> mTRFweights;
    if (isTRF) {
      mTRFweights["2bex"] =
          *(Tools::Instance().GetTreeValue<float>(mEvent, "trf_weight_77_2ex"));
      mTRFweights["3bex"] =
          *(Tools::Instance().GetTreeValue<float>(mEvent, "trf_weight_77_3ex"));
      mTRFweights["4bin"] =
          *(Tools::Instance().GetTreeValue<float>(mEvent, "trf_weight_77_4in"));
      std::vector<float> tmpTRFweights =
          *(Tools::Instance().GetTreeValue<vector<float>>(mEvent,
                                                          "trf_weight_77_in"));
      mTRFweights["2bin"] = tmpTRFweights[2];
      mTRFweights["3bin"] = tmpTRFweights[3];
    }
    // Fake leptons removal
    if (mcChannel != 0) {
      if (doTTbarMerge) {
        if (!doTTbarCombination(mEvent))
          return false;
      }
      if (isFake(mEvent))
        mSample = "Fakes";
    }
    // Calculate Variables
    // calculator->CalculateVariables(mEvent);
    for (auto region : mRegions) {
      if (isTRF) {
        if (region.find("2bex") != string::npos) {
          mWeights = mWeights * mTRFweights["2bex"];
        }
        if (region.find("3bex") != string::npos) {
          mWeights = mWeights * mTRFweights["3bex"];
        }
        if (region.find("4bin") != string::npos) {
          mWeights = mWeights * mTRFweights["4bin"];
        }
      }
      if (isNominal) {
        mRawYields.at(region).at(mSample) += 1;
        mWeightedYields.at(region).at(mSample) += (weights["norm"] * mWeights);
        mWeightedYieldsError.at(region).at(mSample) +=
            (weights["norm"] * mWeights) * (weights["norm"] * mWeights);
      }
      std::vector<string> vars = mConfig->GetRegionVars(region);
      for (auto var : vars) {
        float value;
        if (find(mVarToCalc.begin(), mVarToCalc.end(), var) !=
            mVarToCalc.end()) {
          value = calculator->GetVarValue(var);
        } else {
          string valueType = Tools::Instance().GetValueType(mEvent, var);
          if (!isTRF) {
            if (valueType == "int" || valueType == "Int_t") {
              value = *(Tools::Instance().GetTreeValue<int>(mEvent, var));
            } else {
              value = *(Tools::Instance().GetTreeValue<float>(mEvent, var));
            }
          } else {
            if (valueType == "int" || valueType == "Int_t") {
              value = *(Tools::Instance().GetTreeValue<int>(mEvent, var));
            } else {
              if (find(mTRFvariables.begin(), mTRFvariables.end(), var) !=
                  mTRFvariables.end()) {
                if (region.find("2bex") != string::npos) {
                  string tmpVar = var + "_77_2ex";
                  value =
                      *(Tools::Instance().GetTreeValue<float>(mEvent, tmpVar));
                } else if (region.find("3bex") != string::npos) {
                  string tmpVar = var + "_77_3ex";
                  value =
                      *(Tools::Instance().GetTreeValue<float>(mEvent, tmpVar));
                } else if (region.find("4bin") != string::npos) {
                  string tmpVar = var + "_77_4in";
                  value =
                      *(Tools::Instance().GetTreeValue<float>(mEvent, tmpVar));
                } else {
                  value = *(Tools::Instance().GetTreeValue<float>(mEvent, var));
                }
              } else {
                value = *(Tools::Instance().GetTreeValue<float>(mEvent, var));
              }
            }
          }
        }
        string SysName;
        if (isNominal)
          SysName = w;
        else
          SysName = tmpTreeName;
        string hname = Tools::Instance().GenName(var, region, mSample, SysName);
        if (hs->HasHist(hname)) {
          hs->GetHist(hname)->Fill(value, mWeights * weights["norm"]);
        } else {
          std::vector<float> bins = mConfig->GetVarBinning(region, var);
          hs->AddHist(hname, bins[0], bins[1], bins[2]);
          hs->GetHist(hname)->Fill(value, mWeights * weights["norm"]);
        }
      }
    }
  }
  return true;
}

bool MakeHists::finalize(TFile *fFile) {
  FillYields();
  hs->SaveAllHists(fFile);
  printf("MakeHists:: finalize:: MakeHists has finished running!\n");
  return true;
}

void MakeHists::InitYields(DSHandler *ds) {
  std::vector<string> regions = mConfig->GetRegions();
  std::vector<string> samples = ds->GetAllTypes();
  int nSamples = samples.size();
  int nRegions = regions.size();
  if (!(find(samples.begin(), samples.end(), "ttbar") == samples.end())) {
    nSamples += 2;
  }
  hs->AddHist2D("hist_raw_yields", nSamples, 0, nSamples, nRegions, 0,
                nRegions);
  hs->AddHist2D("hist_weighted_yields", nSamples, 0, nSamples, nRegions, 0,
                nRegions);
  for (auto region : regions) {
    for (auto sample : samples) {
      if (sample == "ttbar") {
        mRawYields[region]["ttlight"] = 0;
        mWeightedYields[region]["ttlight"] = 0;
        mWeightedYieldsError[region]["ttlight"] = 0;
        mRawYields[region]["ttcc"] = 0;
        mWeightedYields[region]["ttcc"] = 0;
        mWeightedYieldsError[region]["ttcc"] = 0;
        mRawYields[region]["ttbb"] = 0;
        mWeightedYields[region]["ttbb"] = 0;
        mWeightedYieldsError[region]["ttbb"] = 0;
      } else {
        mRawYields[region][sample] = 0;
        mWeightedYields[region][sample] = 0;
        mWeightedYieldsError[region][sample] = 0;
      }
    }
  }
  int nx = 1, ny = 1;
  for (auto region : regions) {
    hs->GetHist2D("hist_raw_yields")
        ->GetYaxis()
        ->SetBinLabel(ny, region.c_str());
    hs->GetHist2D("hist_weighted_yields")
        ->GetYaxis()
        ->SetBinLabel(ny++, region.c_str());
  }
  for (auto sample : samples) {
    if (sample != "ttbar") {
      hs->GetHist2D("hist_raw_yields")
          ->GetXaxis()
          ->SetBinLabel(nx, sample.c_str());
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, sample.c_str());
    } else {
      hs->GetHist2D("hist_raw_yields")->GetXaxis()->SetBinLabel(nx, "ttlight");
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, "ttlight");
      hs->GetHist2D("hist_raw_yields")->GetXaxis()->SetBinLabel(nx, "ttcc");
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, "ttcc");
      hs->GetHist2D("hist_raw_yields")->GetXaxis()->SetBinLabel(nx, "ttbb");
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(nx++, "ttbb");
    }
  }
}

void MakeHists::FillYields() {
  int nx = hs->GetHist2D("hist_raw_yields")->GetNbinsX();
  int ny = hs->GetHist2D("hist_raw_yields")->GetNbinsY();
  for (int ix = 1; ix < nx + 1; ix++) {
    string xname =
        hs->GetHist2D("hist_raw_yields")->GetXaxis()->GetBinLabel(ix);
    for (int iy = 1; iy < ny + 1; iy++) {
      string yname =
          hs->GetHist2D("hist_raw_yields")->GetYaxis()->GetBinLabel(iy);
      string tmpName = yname + "_" + xname;
      printf("HistsGen:: FillYields:: Filling Yields %s\n", tmpName.c_str());
      float raw, weighted, weighted_error;
      raw = mRawYields.at(yname).at(xname);
      weighted = mWeightedYields.at(yname).at(xname);
      weighted_error = mWeightedYieldsError.at(yname).at(xname);
      hs->GetHist2D("hist_raw_yields")->SetBinContent(ix, iy, raw);
      hs->GetHist2D("hist_weighted_yields")->SetBinContent(ix, iy, weighted);
      hs->GetHist2D("hist_weighted_yields")
          ->SetBinError(ix, iy, TMath::Sqrt(weighted_error));
    }
  }
}
