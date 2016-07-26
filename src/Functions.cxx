#include "Functions.hpp"

bool isFake(TTree *event) {
  bool fake = false;
  int mcChannel =
      *(Tools::Instance().GetTreeValue<int>(event, "mcChannelNumber"));
  if (mcChannel == 0)
    return false;
  string mSample = Tools::Instance().GetSampleType(mcChannel);
  if (mSample == "Fakes")
    return true;

  int nEl = *(Tools::Instance().GetTreeValue<int>(event, "nElectrons"));
  int nMu = *(Tools::Instance().GetTreeValue<int>(event, "nMuons"));

  std::vector<int> el_true_type = *(
      Tools::Instance().GetTreeValue<std::vector<int>>(event, "el_true_type"));
  std::vector<int> el_true_origin =
      *(Tools::Instance().GetTreeValue<std::vector<int>>(event,
                                                         "el_true_origin"));
  std::vector<int> mu_true_type = *(
      Tools::Instance().GetTreeValue<std::vector<int>>(event, "mu_true_type"));
  std::vector<int> mu_true_origin =
      *(Tools::Instance().GetTreeValue<std::vector<int>>(event,
                                                         "mu_true_origin"));
  std::vector<int> el_true_originbkg =
      *(Tools::Instance().GetTreeValue<std::vector<int>>(event,
                                                         "el_true_originbkg"));
  int et1 = -1;
  int et2 = -1;
  int eo1 = -1;
  int eo2 = -1;
  int eb1 = -1;
  int eb2 = -1;

  if (nEl == 0 && nMu == 2) {
    et1 = mu_true_type.at(0);
    et2 = mu_true_type.at(1);
    eo1 = mu_true_origin.at(0);
    eo2 = mu_true_origin.at(1);
    if (!(et1 == 6 && et2 == 6 &&
          (eo1 == 10 || eo1 == 12 || eo1 == 13 || eo1 == 14) &&
          (eo2 == 10 || eo2 == 12 || eo2 == 13 || eo2 == 14)))
      fake = true;
  } else if (nEl == 2 && nMu == 0) {
    et1 = el_true_type.at(0);
    et2 = el_true_type.at(1);
    eo1 = el_true_origin.at(0);
    eo2 = el_true_origin.at(1);
    eb1 = el_true_originbkg.at(0);
    eb2 = el_true_originbkg.at(1);

    bool pass =
        ((et1 == 2 && (eo1 == 10 || eo1 == 12 || eo1 == 13 || eo1 == 14)) ||
         (et1 == 4 && eo1 == 5 &&
          (eb1 == 10 || eb1 == 12 || eb1 == 13 || eb1 == 14)));
    pass = pass &&
           ((et2 == 2 && (eo2 == 10 || eo2 == 12 || eo2 == 13 || eo2 == 14)) ||
            (et2 == 4 && eo2 == 5 &&
             (eb2 == 10 || eb2 == 12 || eb2 == 13 || eb2 == 14)));
    if (!pass)
      fake = true;
  } else if (nEl == 1 && nMu == 1) {
    et1 = el_true_type.at(0);
    et2 = mu_true_type.at(0);
    eo1 = el_true_origin.at(0);
    eo2 = mu_true_origin.at(0);
    eb1 = el_true_originbkg.at(0);

    bool pass = et2 == 6 && (eo2 == 10 || eo2 == 12 || eo2 == 13 || eo2 == 14);
    pass = pass &&
           ((et1 == 2 && (eo1 == 10 || eo1 == 12 || eo1 == 13 || eo1 == 14)) ||
            (et1 == 4 && eo1 == 5 &&
             (eb1 == 10 || eb1 == 12 || eb1 == 13 || eb1 == 14)));
    if (!pass)
      fake = true;
  }

  return fake;
}

bool doTTbarCombination(TTree *event) {
  int TopHeavyFlavorFilterFlag =
      *(Tools::Instance().GetTreeValue<int>(event, "TopHeavyFlavorFilterFlag"));
  bool truth_top_dilep_filter =
      *(Tools::Instance().GetTreeValue<bool>(event, "truth_top_dilep_filter"));
  int mcChannel =
      *(Tools::Instance().GetTreeValue<int>(event, "mcChannelNumber"));

  if (mcChannel == 410000 &&
      (TopHeavyFlavorFilterFlag == 5 || truth_top_dilep_filter == true))
    return false;
  if (mcChannel == 410120 && truth_top_dilep_filter == true)
    return false;
  if (mcChannel == 410009 && TopHeavyFlavorFilterFlag == 5)
    return false;
  return true;
}

string doHeavyFlavor(TTree *event) {
  int HF_Classification =
      *(Tools::Instance().GetTreeValue<int>(event, "HF_Classification"));
  if (HF_Classification == 0)
    return "ttlight";
  else if (abs(HF_Classification) > 0 && abs(HF_Classification) < 100)
    return "ttcc";
  else
    return "ttbb";
}

void GetSysWeights(TTree *event, std::map<string, float> &mSysWeights) {
  mSysWeights.clear();
  mSysWeights["ttbb_CSS_KIN"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_CSS_KIN"));
  mSysWeights["ttbb_NNPDF"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_NNPDF"));
  mSysWeights["ttbb_MSTW"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_MSTW"));
  mSysWeights["ttbb_Q_CMMPS"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_Q_CMMPS"));
  mSysWeights["ttbb_glosoft"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_glosoft"));
  mSysWeights["ttbb_defaultX05"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_defaultX05"));
  mSysWeights["ttbb_defaultX2"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_defaultX2"));
  mSysWeights["ttbb_MPI_UP"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_MPIup"));
  mSysWeights["ttbb_MPI_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_MPIdown"));
  mSysWeights["ttbb_MPIfactor"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_MPIfactor"));
  mSysWeights["ttbb_aMcAtNloHpp"] = *(
      Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_aMcAtNloHpp"));
  mSysWeights["ttbb_aMcAtNloPy8"] = *(
      Tools::Instance().GetTreeValue<float>(event, "weight_ttbb_aMcAtNloPy8"));
  mSysWeights["pileup_UP"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_pileup_UP"));
  mSysWeights["pileup_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_pileup_DOWN"));
  mSysWeights["jvt_UP"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_jvt_UP"));
  mSysWeights["jvt_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(event, "weight_jvt_DOWN"));
  mSysWeights["leptonSF_EL_SF_Trigger_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_EL_SF_Trigger_UP"));
  mSysWeights["leptonSF_EL_SF_Trigger_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_EL_SF_Trigger_DOWN"));
  mSysWeights["leptonSF_EL_SF_Reco_UP"] =
      *(Tools::Instance().GetTreeValue<float>(event,
                                              "weight_leptonSF_EL_SF_Reco_UP"));
  mSysWeights["leptonSF_EL_SF_Reco_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_EL_SF_Reco_DOWN"));
  mSysWeights["leptonSF_EL_SF_ID_UP"] = *(Tools::Instance().GetTreeValue<float>(
      event, "weight_leptonSF_EL_SF_ID_UP"));
  mSysWeights["leptonSF_EL_SF_ID_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(event,
                                              "weight_leptonSF_EL_SF_ID_DOWN"));
  mSysWeights["leptonSF_EL_SF_Isol_UP"] =
      *(Tools::Instance().GetTreeValue<float>(event,
                                              "weight_leptonSF_EL_SF_Isol_UP"));
  mSysWeights["leptonSF_EL_SF_Isol_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_EL_SF_Isol_DOWN"));
  mSysWeights["leptonSF_MU_SF_Trigger_STAT_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Trigger_STAT_UP"));
  mSysWeights["leptonSF_MU_SF_Trigger_STAT_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Trigger_STAT_DOWN"));
  mSysWeights["leptonSF_MU_SF_Trigger_SYST_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Trigger_SYST_UP"));
  mSysWeights["leptonSF_MU_SF_Trigger_SYST_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Trigger_SYST_DOWN"));
  mSysWeights["leptonSF_MU_SF_ID_STAT_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_STAT_UP"));
  mSysWeights["leptonSF_MU_SF_ID_STAT_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_STAT_DOWN"));
  mSysWeights["leptonSF_MU_SF_ID_SYST_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_SYST_UP"));
  mSysWeights["leptonSF_MU_SF_ID_SYST_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_SYST_DOWN"));
  mSysWeights["leptonSF_MU_SF_ID_STAT_LOWPT_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_STAT_LOWPT_UP"));
  mSysWeights["leptonSF_MU_SF_ID_STAT_LOWPT_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_STAT_LOWPT_DOWN"));
  mSysWeights["leptonSF_MU_SF_ID_SYST_LOWPT_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_SYST_LOWPT_UP"));
  mSysWeights["leptonSF_MU_SF_ID_SYST_LOWPT_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_ID_SYST_LOWPT_DOWN"));
  mSysWeights["leptonSF_MU_SF_Isol_STAT_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Isol_STAT_UP"));
  mSysWeights["leptonSF_MU_SF_Isol_STAT_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Isol_STAT_DOWN"));
  mSysWeights["leptonSF_MU_SF_Isol_SYST_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Isol_SYST_UP"));
  mSysWeights["leptonSF_MU_SF_Isol_SYST_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_Isol_SYST_DOWN"));
  mSysWeights["leptonSF_MU_SF_TTVA_STAT_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_TTVA_STAT_UP"));
  mSysWeights["leptonSF_MU_SF_TTVA_STAT_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_TTVA_STAT_DOWN"));
  mSysWeights["leptonSF_MU_SF_TTVA_SYST_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_TTVA_SYST_UP"));
  mSysWeights["leptonSF_MU_SF_TTVA_SYST_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_leptonSF_MU_SF_TTVA_SYST_DOWN"));
  std::vector<float> bTagSF_70_eigenvars_B_UP =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(
          event, "weight_bTagSF_70_eigenvars_B_up"));
  std::vector<float> bTagSF_70_eigenvars_B_DOWN =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(
          event, "weight_bTagSF_70_eigenvars_B_down"));
  std::vector<float> bTagSF_70_eigenvars_C_UP =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(
          event, "weight_bTagSF_70_eigenvars_C_up"));
  std::vector<float> bTagSF_70_eigenvars_C_DOWN =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(
          event, "weight_bTagSF_70_eigenvars_C_down"));
  std::vector<float> bTagSF_70_eigenvars_Light_UP =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(
          event, "weight_bTagSF_70_eigenvars_Light_up"));
  std::vector<float> bTagSF_70_eigenvars_Light_DOWN =
      *(Tools::Instance().GetTreeValue<std::vector<float>>(
          event, "weight_bTagSF_70_eigenvars_Light_down"));
  mSysWeights["bTagSF_70_extrapolation_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_bTagSF_70_extrapolation_up"));
  mSysWeights["bTagSF_70_extrapolation_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_bTagSF_70_extrapolation_down"));
  mSysWeights["bTagSF_70_extrapolation_from_charm_UP"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_bTagSF_70_extrapolation_from_charm_up"));
  mSysWeights["bTagSF_70_extrapolation_from_charm_DOWN"] =
      *(Tools::Instance().GetTreeValue<float>(
          event, "weight_bTagSF_70_extrapolation_from_charm_down"));
  mSysWeights["bTagSF_70_eigenvars_B1_UP"] = bTagSF_70_eigenvars_B_UP.at(0);
  mSysWeights["bTagSF_70_eigenvars_B2_UP"] = bTagSF_70_eigenvars_B_UP.at(1);
  mSysWeights["bTagSF_70_eigenvars_B3_UP"] = bTagSF_70_eigenvars_B_UP.at(2);
  mSysWeights["bTagSF_70_eigenvars_B4_UP"] = bTagSF_70_eigenvars_B_UP.at(3);
  mSysWeights["bTagSF_70_eigenvars_B5_UP"] = bTagSF_70_eigenvars_B_UP.at(4);
  mSysWeights["bTagSF_70_eigenvars_B1_DOWN"] = bTagSF_70_eigenvars_B_DOWN.at(0);
  mSysWeights["bTagSF_70_eigenvars_B2_DOWN"] = bTagSF_70_eigenvars_B_DOWN.at(1);
  mSysWeights["bTagSF_70_eigenvars_B3_DOWN"] = bTagSF_70_eigenvars_B_DOWN.at(2);
  mSysWeights["bTagSF_70_eigenvars_B4_DOWN"] = bTagSF_70_eigenvars_B_DOWN.at(3);
  mSysWeights["bTagSF_70_eigenvars_B5_DOWN"] = bTagSF_70_eigenvars_B_DOWN.at(4);
  mSysWeights["bTagSF_70_eigenvars_C1_UP"] = bTagSF_70_eigenvars_C_UP.at(0);
  mSysWeights["bTagSF_70_eigenvars_C2_UP"] = bTagSF_70_eigenvars_C_UP.at(1);
  mSysWeights["bTagSF_70_eigenvars_C3_UP"] = bTagSF_70_eigenvars_C_UP.at(2);
  mSysWeights["bTagSF_70_eigenvars_C4_UP"] = bTagSF_70_eigenvars_C_UP.at(3);
  mSysWeights["bTagSF_70_eigenvars_C1_DOWN"] = bTagSF_70_eigenvars_C_DOWN.at(0);
  mSysWeights["bTagSF_70_eigenvars_C2_DOWN"] = bTagSF_70_eigenvars_C_DOWN.at(1);
  mSysWeights["bTagSF_70_eigenvars_C3_DOWN"] = bTagSF_70_eigenvars_C_DOWN.at(2);
  mSysWeights["bTagSF_70_eigenvars_C4_DOWN"] = bTagSF_70_eigenvars_C_DOWN.at(3);
  mSysWeights["bTagSF_70_eigenvars_Light1_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(0);
  mSysWeights["bTagSF_70_eigenvars_Light2_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(1);
  mSysWeights["bTagSF_70_eigenvars_Light3_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(2);
  mSysWeights["bTagSF_70_eigenvars_Light4_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(3);
  mSysWeights["bTagSF_70_eigenvars_Light5_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(4);
  mSysWeights["bTagSF_70_eigenvars_Light6_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(5);
  mSysWeights["bTagSF_70_eigenvars_Light7_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(6);
  mSysWeights["bTagSF_70_eigenvars_Light8_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(7);
  mSysWeights["bTagSF_70_eigenvars_Light9_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(8);
  mSysWeights["bTagSF_70_eigenvars_Light10_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(9);
  mSysWeights["bTagSF_70_eigenvars_Light11_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(10);
  mSysWeights["bTagSF_70_eigenvars_Light12_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(11);
  mSysWeights["bTagSF_70_eigenvars_Light13_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(12);
  mSysWeights["bTagSF_70_eigenvars_Light14_UP"] =
      bTagSF_70_eigenvars_Light_UP.at(13);
  mSysWeights["bTagSF_70_eigenvars_Light1_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(0);
  mSysWeights["bTagSF_70_eigenvars_Light2_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(1);
  mSysWeights["bTagSF_70_eigenvars_Light3_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(2);
  mSysWeights["bTagSF_70_eigenvars_Light4_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(3);
  mSysWeights["bTagSF_70_eigenvars_Light5_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(4);
  mSysWeights["bTagSF_70_eigenvars_Light6_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(5);
  mSysWeights["bTagSF_70_eigenvars_Light7_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(6);
  mSysWeights["bTagSF_70_eigenvars_Light8_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(7);
  mSysWeights["bTagSF_70_eigenvars_Light9_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(8);
  mSysWeights["bTagSF_70_eigenvars_Light10_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(9);
  mSysWeights["bTagSF_70_eigenvars_Light11_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(10);
  mSysWeights["bTagSF_70_eigenvars_Light12_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(11);
  mSysWeights["bTagSF_70_eigenvars_Light13_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(12);
  mSysWeights["bTagSF_70_eigenvars_Light14_DOWN"] =
      bTagSF_70_eigenvars_Light_DOWN.at(13);
}
