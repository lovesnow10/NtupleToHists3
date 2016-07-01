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
  return false;
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
