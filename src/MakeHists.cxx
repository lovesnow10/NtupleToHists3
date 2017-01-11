#include "Functions.hpp"
#include "MakeHists.hpp"
#include "TMath.h"

using namespace std;

bool MakeHists::initialize(ConfigParser *config, DSHandler *ds) {
  mConfig = config;
  hs = new HistStore();
  calculator = new VariableCalculator();
  // InitYields(ds);
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
  bool doFakes = bControl.at("doFakes");
  bool doFakesOnly = bControl.at("doFakesOnly");

  calculator->CalculateVariables(mEvent);
  int mcChannel =
      *(Tools::Instance().GetTreeValue<int>(mEvent, "mcChannelNumber"));
  mSample = Tools::Instance().GetSampleType(mcChannel);
  if (doSysmetic && isNominal) {
    GetSysWeights(mEvent, mSysWeights);
    for (auto w : weights) {
      if (w.first == "norm" || w.first == "weight_NNLO")
        continue;
      string tmpName = w.first.substr(7);
      mSysWeights[tmpName] = w.second;
    }
  }
  mSysWeights["nominal"] = 1.0;

  // Get Weights!
  // float tmpWeights = Tools::Instance().GetWeight(mcChannel);
  float weight_mc = 1.0;
  float weight_pileup = 1.0;
  float weight_jvt = 1.0;
  float weight_leptonSF = 1.0;
  float weight_bTagSF_70 = 1.0;
  float weight_ttbb_Nominal = 1.0;

  // No weight for DATA!
  if (mcChannel != 0) {
    if (mSample == "ttbar") {
      mSample = doHeavyFlavor(mEvent);
    }
    // Fake leptons removal
    if (doTTbarMerge) {
      if (!doTTbarCombination(mEvent))
        return false;
    }
    if (isFake(mEvent))
      mSample = "Fakes";
  }
  if (mSample == "Fakes")
  {
    if (!doFakes) return false;
  }
  else{
    if (doFakes && doFakesOnly) return false;
  }

  for (auto w : mSysName) {
    if (!isNominal && w != "nominal")
      continue;

    float weight = mSysWeights.at(w);
    float mWeights = Tools::Instance().GetWeight(mcChannel);

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

      weight_mc = *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_mc"));
      weight_pileup =
          *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_pileup"));
      weight_jvt =
          *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_jvt"));
      weight_leptonSF =
          *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_leptonSF"));
      weight_bTagSF_70 =
          *(Tools::Instance().GetTreeValue<float>(mEvent, "weight_bTagSF_70"));
      weight_ttbb_Nominal = *(
          Tools::Instance().GetTreeValue<float>(mEvent, "weight_ttbb_Nominal"));

      if (mSample != "ttbb") {
        if (w.find("leptonSF_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_jvt * weight_bTagSF_70 * weight_ttbb_Nominal *
                     weights["weight_NNLO"];
        else if (w.find("bTagSF_70_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_jvt * weight_leptonSF * weight_ttbb_Nominal *
                     weights["weight_NNLO"];
        else if (w.find("pileup_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_jvt *
                     weight_leptonSF * weight_bTagSF_70 * weight_ttbb_Nominal *
                     weights["weight_NNLO"];
        else if (w.find("jvt_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_leptonSF * weight_bTagSF_70 * weight_ttbb_Nominal *
                     weights["weight_NNLO"];
        else if (w.find("ttbb_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_leptonSF * weight_bTagSF_70 *
                     weights["weight_NNLO"];
        else if (w.find("NNLO_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_leptonSF * weight_bTagSF_70 * weight_ttbb_Nominal;
        else
          mWeights = mWeights * weight_mc * weight_pileup * weight_jvt *
                     weight_leptonSF * weight_bTagSF_70 * weight_ttbb_Nominal *
                     weights["weight_NNLO"];
      } else {
        if (w.find("leptonSF_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_jvt * weight_bTagSF_70 * weight_ttbb_Nominal;
        else if (w.find("bTagSF_70_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_jvt * weight_leptonSF * weight_ttbb_Nominal;
        else if (w.find("pileup_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_jvt *
                     weight_leptonSF * weight_bTagSF_70 * weight_ttbb_Nominal;
        else if (w.find("jvt_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_leptonSF * weight_bTagSF_70 * weight_ttbb_Nominal;
        else if (w.find("ttbb_") != string::npos)
          mWeights = mWeights * weight * weight_mc * weight_pileup *
                     weight_leptonSF * weight_bTagSF_70;
        else
          mWeights = mWeights * weight_mc * weight_pileup * weight_jvt *
                     weight_leptonSF * weight_bTagSF_70 * weight_ttbb_Nominal;
      }
    }

    // Selectioins!
    mRegions.clear();
    for (int ientry = 0; ientry < formulas->GetEntries(); ++ientry) {
      TTreeFormula *formula = formulas->GetFormula(ientry);
      if (formula->EvalInstance()) {
        mRegions.push_back(formula->GetName());
      }
    }
    if (mRegions.empty())
      return false;
    // continue;

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
        if (Tools::Instance().CheckYieldsMap(mRawYields, region, mSample)) {
          mRawYields.at(region).at(mSample) += 1;
          mWeightedYields.at(region).at(mSample) +=
              (weights["norm"] * mWeights);
          mWeightedYieldsError.at(region).at(mSample) +=
              (weights["norm"] * mWeights) * (weights["norm"] * mWeights);
        } else {
          mRawYields[region][mSample] = 1;
          mWeightedYields[region][mSample] = weights["norm"] * mWeights;
          mWeightedYieldsError[region][mSample] =
              weights["norm"] * mWeights * weights["norm"] * mWeights;
        }
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

void MakeHists::FillYields() {
  int nRegions, nSamples;
  nRegions = mRawYields.size();
  nSamples = (mRawYields.begin()->second).size();

  hs->AddHist2D("hist_raw_yields", nSamples, 0, nSamples, nRegions, 0,
                nRegions);
  hs->AddHist2D("hist_weighted_yields", nSamples, 0, nSamples, nRegions, 0,
                nRegions);
  int iRegion = 1;
  for (auto region : mRawYields) {
    hs->GetHist2D("hist_raw_yields")
        ->GetYaxis()
        ->SetBinLabel(iRegion, region.first.c_str());
    hs->GetHist2D("hist_weighted_yields")
        ->GetYaxis()
        ->SetBinLabel(iRegion, region.first.c_str());
    int iSample = 1;
    for (auto sample : region.second) {
      hs->GetHist2D("hist_raw_yields")
          ->GetXaxis()
          ->SetBinLabel(iSample, sample.first.c_str());
      hs->GetHist2D("hist_weighted_yields")
          ->GetXaxis()
          ->SetBinLabel(iSample, sample.first.c_str());
      string tmpName = region.first + "_" + sample.first;
      printf("HistsGen:: FillYields:: Filling Yields %s\n", tmpName.c_str());
      float raw = sample.second;
      float weighted = mWeightedYields.at(region.first).at(sample.first);
      float weighted_error =
          mWeightedYieldsError.at(region.first).at(sample.first);
      hs->GetHist2D("hist_raw_yields")->SetBinContent(iSample, iRegion, raw);
      hs->GetHist2D("hist_weighted_yields")
          ->SetBinContent(iSample, iRegion, weighted);
      hs->GetHist2D("hist_weighted_yields")
          ->SetBinError(iSample, iRegion, TMath::Sqrt(weighted_error));
      iSample++;
    }
    iRegion++;
  }
}
