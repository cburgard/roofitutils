#include "RooFitUtils/AbsMeasurement.h"
#include "RooFitUtils/ExtendedMinimizer.h"

#include "TArrayI.h"
#include "TH1.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TString.h"

#include "RooBinning.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"

// ____________________________________________________________________________|__________

RooFitUtils::AbsMeasurement::AbsMeasurement(const std::string &MeasurementName,
                                            RooWorkspace *ws,
                                            const std::string &ModelConfigName,
                                            const std::string &DataName)
    : AbsMeasurement(MeasurementName, ws ? ws->GetName() : "", ModelConfigName,
                     DataName) {
  // Constructor
  fWorkSpace = ws;
}

// ____________________________________________________________________________|__________

RooFitUtils::AbsMeasurement::AbsMeasurement(const std::string &MeasurementName,
                                            const std::string &FileName,
                                            const std::string &WorkspaceName,
                                            const std::string &ModelConfigName,
                                            const std::string &DataName)
    : AbsMeasurement(MeasurementName, WorkspaceName, ModelConfigName,
                     DataName) {
  // Constructor
  fFileName = FileName;
}

// ____________________________________________________________________________|__________

RooFitUtils::AbsMeasurement::AbsMeasurement(const std::string &MeasurementName,
                                            const std::string &WorkspaceName,
                                            const std::string &ModelConfigName,
                                            const std::string &DataName)
    : TNamed(MeasurementName.c_str(), MeasurementName.c_str()),
      fWorkspaceName(WorkspaceName), fModelConfigName(ModelConfigName),
      fDataName(DataName) {
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(),
                                                    fMinimizerAlgo.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);

  coutP(InputArguments) << "AbsMeasurement::AbsMeasurement(" << fName
                        << ") created" << std::endl;
}

// ____________________________________________________________________________|__________

RooFitUtils::AbsMeasurement::~AbsMeasurement() {
  // Destructor
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::enablePruning(const std::string &poi,
                                                const std::string &filter,
                                                const std::string &weight,
                                                const std::string &threshold,
                                                int additionalDigit) {
  coutP(InputArguments) << "AbsMeasurement::enablePruning(" << fName
                        << ") pruning enabled" << std::endl;
  std::cout << "\n*************************************************************"
               "***************"
            << std::endl;
  std::cout << "* WARNING: Parameters will be pruned automatically based on "
               "the change of  *"
            << std::endl;
  std::cout << "*          Hesse uncertainty on the specified parameters of "
               "interest       *"
            << std::endl;
  std::cout << "*          respecting the specified thresholds.                "
               "            *"
            << std::endl;
  std::cout << "*                                                              "
               "            *"
            << std::endl;
  std::cout << "*          It is recommended to verify the list of pruned "
               "parameters!      *"
            << std::endl;
  std::cout << "***************************************************************"
               "*************\n"
            << std::endl;

  fIsPrunable = kTRUE;
  fPruningPoi = poi;
  fPruningFilter = filter;
  fPruningWeight = weight;
  fPruningThreshold = threshold;
  fPruningAdditionalDigit = additionalDigit;
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::SetPrunedNuisanceParameters(
    std::string parameters) {
  TString allParameters = parameters;
  allParameters.ReplaceAll(" ", "");
  TObjArray *allParametersArray = allParameters.Tokenize(",");
  unsigned int numParameters = allParametersArray->GetEntries();
  for (unsigned int itrPar = 0; itrPar < numParameters; ++itrPar) {
    TString thisPar =
        ((TObjString *)allParametersArray->At(itrPar))->GetString();
    fPrunedNuisanceParameters.push_back(thisPar.Data());
  }
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::PruneNuisanceParameters() {
  coutP(ObjectHandling)
      << "AbsMeasurement::PruneNuisanceParameters(" << fName
      << ") performing fit to compute Hesse matrix for pruning." << std::endl;

  RooArgSet *allParameters = fModelConfig->GetPdf()->getParameters(*fData);
  fWorkSpace->saveSnapshot("tmpBeforePruning", *allParameters);

  for (RooLinkedListIter it =
           fModelConfig->GetParametersOfInterest()->iterator();
       RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());) {
    v->setConstant(1);
  }

  TString allPoi = fPruningPoi;
  allPoi.ReplaceAll(" ", "");
  TObjArray *allPoiArray = allPoi.Tokenize(",");
  unsigned int numPoi = allPoiArray->GetEntries();
  for (unsigned int itrPoi = 0; itrPoi < numPoi; ++itrPoi) {
    TString thisPoi = ((TObjString *)allPoiArray->At(itrPoi))->GetString();
    fWorkSpace->var(thisPoi.Data())->setConstant(0);
  }

  ExtendedMinimizer minimizer("minimizer", fPdf, fData);
  minimizer.minimize(
      RooFit::Minimizer(fMinimizerType.c_str(), fMinimizerAlgo.c_str()),
      RooFit::Strategy(fDefaultStrategy),
      RooFit::Constrain(*fModelConfig->GetNuisanceParameters()),
      RooFit::GlobalObservables(*fModelConfig->GetGlobalObservables()),
      RooFit::NumCPU(fNumCPU, 3), RooFit::Offset(1), RooFit::Optimize(2),
      RooFit::Save(), RooFit::Hesse());
  RooFitResult *fitresult = minimizer.GetFitResult();
  TMatrixDSym hesse = minimizer.GetHesseMatrix();

  fPrunedNuisanceParameters = PruneNuisanceParameters(
      hesse, fitresult, fPruningPoi, fPruningFilter, fPruningWeight,
      fPruningThreshold, fPruningAdditionalDigit, fPrunedNuisanceParameters);

  fWorkSpace->loadSnapshot("tmpBeforePruning");

  delete fitresult;
}

// ____________________________________________________________________________|__________

std::list<std::string> RooFitUtils::AbsMeasurement::PruneNuisanceParameters(
    const TMatrixDSym chesse, RooFitResult *fitresult, const std::string &poi,
    const std::string &filter, const std::string &weight,
    const std::string &threshold, int additionalDigit,
    std::list<std::string> prePrunedParameters) {
  // Ranking/pruning function. Compute symmetric uncertainty on POIs using the
  // Hessian matrix. Remove iteratively a single NP and compare uncertainty from
  // reduced Hessian matrix to initial one. Multiple POIs are combined by
  // computing
  // the weighted average and propagating the partial uncertainties and taking
  // the
  // OR of the individual rankings.

  // Decompose POI names and associated weights and thresholds
  TString allPoi = poi;
  TString allWeight = weight;
  TString allThreshold = threshold;
  TString allFilter = filter;
  allPoi.ReplaceAll(" ", "");
  allWeight.ReplaceAll(" ", "");
  allThreshold.ReplaceAll(" ", "");
  allFilter.ReplaceAll(" ", "");
  TObjArray *allPoiArray = allPoi.Tokenize(",");
  TObjArray *allWeightArray = allWeight.Tokenize(",");
  TObjArray *allThresholdArray = allThreshold.Tokenize(",");
  TObjArray *allFilterArray = allFilter.Tokenize(",");
  unsigned int numPoi = allPoiArray->GetEntries();
  unsigned int numWeight = allWeightArray->GetEntries();
  unsigned int numThreshold = allThresholdArray->GetEntries();
  unsigned int numFilter = allFilterArray->GetEntries();

  // Check if number of POIs and weights match
  if (numPoi != numWeight && numWeight != 1) {
    coutE(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Number of POIs does not match number of "
                             "specified weights. No additional pruning."
                          << std::endl;
    return prePrunedParameters;
  }

  // Check if number of POIs and thresholds match
  if (numPoi != numThreshold && numThreshold != 1) {
    coutE(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Number of POIs does not match number of "
                             "specified thresholds. No additional pruning."
                          << std::endl;
    return prePrunedParameters;
  }

  // Store information in easier usable objects
  std::vector<TString> PoiList;
  std::map<std::string, double> WeightList;
  std::map<std::string, TString> ThresholdList;
  std::map<std::string, double> ThresholdCutList;
  std::vector<TString> FilterList;
  double totalWeights = 0.0;

  for (unsigned int itrPoi = 0; itrPoi < numPoi; ++itrPoi) {
    TString thisPoi = ((TObjString *)allPoiArray->At(itrPoi))->GetString();
    TString thisWeight =
        (numWeight == 1)
            ? ((TObjString *)allWeightArray->At(0))->GetString()
            : ((TObjString *)allWeightArray->At(itrPoi))->GetString();
    TString thisThreshold =
        (numThreshold == 1)
            ? ((TObjString *)allThresholdArray->At(0))->GetString()
            : ((TObjString *)allThresholdArray->At(itrPoi))->GetString();
    double tmpWeight = atof(thisWeight.Data());

    TObjArray *thisThresholdArray = thisThreshold.Tokenize(":");
    TString thisSaveThreshold =
        ((TObjString *)thisThresholdArray->At(0))->GetString();
    double thisThresholdCut =
        (thisThresholdArray->GetEntries() == 1)
            ? 0.0
            : atof(((TObjString *)thisThresholdArray->At(1))
                       ->GetString()
                       .Data());

    PoiList.push_back(thisPoi);
    WeightList[thisPoi.Data()] = tmpWeight;
    ThresholdList[thisPoi.Data()] = thisSaveThreshold;
    ThresholdCutList[thisPoi.Data()] = thisThresholdCut;

    totalWeights += tmpWeight;

    coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Parsed POI " << thisPoi.Data()
                          << ", assigning weight " << tmpWeight
                          << " and threshold " << thisThreshold.Data()
                          << std::endl;
  }

  if (numFilter > 0) {
    coutI(ObjectHandling)
        << "AbsMeasurement::PruneNuisanceParameters(" << fName
        << ") Nuisance parameters will be pruned in the following order:";
    for (unsigned int itrSet = 0; itrSet < numFilter; ++itrSet) {
      TString thisSet = ((TObjString *)allFilterArray->At(itrSet))->GetString();
      FilterList.push_back(thisSet);
      std::cout << " " << thisSet.Data() << ";";
    }
    std::cout << std::endl;
  } else {
    coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Include all parameters in ranking" << std::endl;
    FilterList.push_back(".*");
  }

  // Clean up
  delete allPoiArray;
  delete allWeightArray;
  delete allThresholdArray;
  delete allFilterArray;

  // Normalize weights
  for (auto itrWeight = WeightList.begin(); itrWeight != WeightList.end();
       ++itrWeight) {
    WeightList[itrWeight->first] = itrWeight->second / totalWeights;
  }

  // Print the RooFitResult used for pruning
  fitresult->Print();

  // Get initial information, such as
  //   - covariance matrix
  //   - floating parameters
  //   - position of POIs in the list of floating parameters
  //   - best fit and initial Hesse error
  // TMatrixDSym cov = fitresult->covarianceMatrix();
  TMatrixDSym hes(chesse);
  TMatrixDSym cov(hes.Invert());
  RooArgList pars = fitresult->floatParsFinal();

  std::map<std::string, int> index;
  std::map<std::string, double> initErr;
  std::map<std::string, double> bestFit;
  std::map<std::string, double> initErrSig;

  // map<std::string, double> ThresholdList;
  for (std::vector<TString>::const_iterator itr = PoiList.begin(),
                                            end = PoiList.end();
       itr != end; ++itr) {
    int thisIndex = pars.index(itr->Data());
    double thisInitErr = sqrt(cov[thisIndex][thisIndex]);
    double thisBestFit = ((RooRealVar *)pars.at(thisIndex))->getVal();
    std::pair<double, double> tmpVals =
        PDGrounding(thisBestFit, thisInitErr, additionalDigit);

    index[itr->Data()] = thisIndex;
    initErr[itr->Data()] = thisInitErr;
    initErrSig[itr->Data()] = tmpVals.second;
    bestFit[itr->Data()] = thisBestFit;

    coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Harvest parameter (" << thisIndex
                          << "): " << *itr << " with Hesse error "
                          << tmpVals.first << " (" << thisBestFit << ") +/- "
                          << tmpVals.second << " (" << thisInitErr << ")"
                          << std::endl;
  }

  // Compute the weighted average of POIs
  double totalBestFit = 0.0;
  for (std::vector<TString>::const_iterator itr = PoiList.begin(),
                                            end = PoiList.end();
       itr != end; ++itr) {
    totalBestFit += WeightList[itr->Data()] * bestFit[itr->Data()];
  }

  // Compute uncertainty on the weighted average of the POIs using error
  // propagation.
  // This takes into account the two-point correlations between the POIs.
  double initTotalErrorSquared = 0.0;

  for (std::vector<TString>::const_iterator itr = PoiList.begin(),
                                            end = PoiList.end();
       itr != end; ++itr) {
    for (std::vector<TString>::const_iterator jtr = PoiList.begin(),
                                              end = PoiList.end();
         jtr != end; ++jtr) {
      initTotalErrorSquared += WeightList[itr->Data()] *
                               WeightList[jtr->Data()] *
                               cov[index[itr->Data()]][index[jtr->Data()]];
    }
  }
  double initTotalError = sqrt(initTotalErrorSquared);

  coutI(ObjectHandling)
      << "AbsMeasurement::PruneNuisanceParameters(" << fName
      << ") Initial averaged POI with propagated uncertainties: "
      << totalBestFit << " +/- " << initTotalError << std::endl;
  ;

  // Ranking starts here.
  // Iteratively remove parameters and recompute Hesse uncertainty on the POI
  // and on
  // the averaged parameter used for the final ranking. Maybe filter list of
  // parameters
  // to prune base on regexp.
  std::list<std::string> order = prePrunedParameters;
  unsigned int allParams = pars.getSize();
  for (auto itr = order.begin(); itr != order.end(); ++itr) {
    std::string name(*itr);
    int par2rem_index = pars.index(name.c_str());
    pars.remove(*pars.at(par2rem_index));
  }
  std::set<std::pair<double, std::string>> uncertsPmO;

  for (std::vector<TString>::const_iterator itrFilter = FilterList.begin(),
                                            end = FilterList.end();
       itrFilter != end; ++itrFilter) {
    TString thisFilter = *itrFilter;
    std::string stringFilter = thisFilter.Data();

    TRegexp reg(thisFilter);
    Ssiz_t dummy(0);

    while ((unsigned int)order.size() < allParams - (numPoi + 1)) {
      std::set<std::pair<double, std::string>> uncerts;
      std::map<std::string, std::set<std::pair<double, std::string>>>
          part_uncerts;

      for (RooLinkedListIter it = pars.iterator();
           RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());) {
        std::string name = v->GetName();

        // Skip removing POIs
        if (allPoi.Contains(name)) {
          continue;
        }

        // Apply filtering on a subset of parameters
        if (stringFilter != "" &&
            reg.Index(TString(v->GetName()), &dummy, 0) == -1) {
          continue;
        }

        // Attach current parameter to a temporary list of parameters that
        // should be removed from
        // Hesse matrix. This includes the previously ranked, and already
        // removed parameters.
        std::list<std::string> names = order;
        names.push_back(name);

        // Grab covariance matrix and floating parameters for doing
        // calculations. These objects will
        // be modified for every parameter, and thus need to be reverted to the
        // original ones in
        // in every iteration.
        // TMatrixDSym tmp_cov = fitresult->covarianceMatrix();
        TMatrixDSym tmp_hesse(chesse);
        RooArgList tmp_pars = fitresult->floatParsFinal();

        // Reduce temporary Hesse matrix by a set of parameters and keep track
        // of the parameters left
        RemoveParameter(tmp_hesse, tmp_pars, names);
        TMatrixDSym tmp_cov = tmp_hesse.Invert();

        // Harvest and store reduced components
        std::map<std::string, int> index_red;
        // map<string, double > errred;
        for (std::vector<TString>::const_iterator itr = PoiList.begin(),
                                                  end = PoiList.end();
             itr != end; ++itr) {
          int thisIndex = tmp_pars.index(itr->Data());
          double thisErr = sqrt(tmp_cov[thisIndex][thisIndex]);

          index_red[itr->Data()] = thisIndex;
          // errred[itr->Data()] = thisErr;

          // Store individual uncertainties for reference
          part_uncerts[itr->Data()].insert(std::make_pair(thisErr, name));
        }

        // Compute uncertainty on the weighted average of the POIs using error
        // propagation.
        // This takes into account the two-point correlations between the POIs.
        double redTotalErrorSquared = 0.0;
        for (std::vector<TString>::const_iterator i = PoiList.begin(),
                                                  end = PoiList.end();
             i != end; ++i) {
          for (std::vector<TString>::const_iterator j = PoiList.begin(),
                                                    end = PoiList.end();
               j != end; ++j) {
            redTotalErrorSquared +=
                WeightList[i->Data()] * WeightList[j->Data()] *
                tmp_cov[index_red[i->Data()]][index_red[j->Data()]];
          }
        }
        double redTotalError = sqrt(redTotalErrorSquared);

        // Store total uncertainty for the ranking for remaining parameters
        uncerts.insert(std::make_pair(redTotalError, name));
      }
      if (uncerts.empty())
        break;

      // Print ranking of remaining parameters based on uncertainty on combined
      // inclusive POI
      std::cout << std::endl;
      coutI(ObjectHandling)
          << "AbsMeasurement::PruneNuisanceParameters(" << fName
          << ") Ranking of remaining parameters based on combined inclusive POI"
          << std::endl;
      PrintRanking(uncerts, initTotalError);

      // Find the lowest ranked parameter which changes none of the
      // uncertainties of the individual
      // parameters above the given threshold
      bool foundPar2Rem = false;
      std::string par2rem = uncerts.rbegin()->second;
      double par2rem_err = uncerts.rbegin()->first;
      unsigned int tmpSum = pars.getSize() + 1 - numPoi;

      // Loop over all remaining parameters in the global ranking
      for (std::set<std::pair<double, std::string>>::reverse_iterator rankitr =
               uncerts.rbegin();
           rankitr != uncerts.rend(); ++rankitr) {
        par2rem = rankitr->second;
        par2rem_err = rankitr->first;
        bool passThreshold = true;
        int tmpIndex = 0;

        coutI(ObjectHandling)
            << "AbsMeasurement::PruneNuisanceParameters(" << fName
            << ") --> Testing " << par2rem << std::endl;
        ;

        // Loop over the individual POIs to test effect of removing proposed
        // parameter on each of them
        for (auto iterator = part_uncerts.begin();
             iterator != part_uncerts.end(); ++iterator) {
          std::string thisPoi = iterator->first;

          coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters("
                                << fName << ")  | " << thisPoi;

          std::set<std::pair<double, std::string>> thisRanking =
              iterator->second;
          unsigned int thisPosition = 1;
          double thisReduction = 1.0;
          double thisInitError = initErr[thisPoi];
          TString thisThreshold = ThresholdList[thisPoi];

          // Find the proposed parameter, the position in the individual ranking
          // and the estimated uncertainty reduction
          for (std::set<std::pair<double, std::string>>::reverse_iterator
                   subitr = thisRanking.rbegin();
               subitr != thisRanking.rend(); ++subitr) {
            std::string varName(subitr->second);
            double thisUncert = subitr->first;

            // Check the individual uncertainty reduction if the parameter is
            // found and exit search
            if (varName == par2rem) {
              thisReduction = thisUncert / thisInitError - 1;

              // Check if the threshold is passed, either the relative reduction
              // or check
              // whether quoted uncertainty would be changed in case of 'auto'
              if (thisThreshold != "auto") {
                double tmpThreshold = atof(thisThreshold.Data());
                if (fabs(thisReduction) > tmpThreshold) {
                  passThreshold = false;
                }
              } else {
                std::pair<double, double> tmpVals =
                    PDGrounding(bestFit[thisPoi], thisUncert, additionalDigit);

                std::ostringstream streamInit, streamRed;
                streamInit << initErrSig[thisPoi];
                streamRed << tmpVals.second;
                TString strInit(streamInit.str());
                TString strRed(streamRed.str());

                strInit.ReplaceAll(".", "");
                strRed.ReplaceAll(".", "");

                while (strInit.Length() < strRed.Length()) {
                  strInit.Append("0");
                }

                while (strRed.Length() < strInit.Length()) {
                  strRed.Append("0");
                }

                double lastInit = atof(strInit.Data());
                double lastRed = atof(strRed.Data());

                std::cout << ": absolute change " << initErrSig[thisPoi]
                          << " --> " << tmpVals.second;

                // Allow the specified non-significant digit to change by a
                // specified value
                if (tmpVals.second < initErrSig[thisPoi]) {
                  if (additionalDigit > 0) {
                    if (!AlmostEqualUlpsAndAbs(lastInit, lastRed,
                                               ThresholdCutList[thisPoi], 4)) {
                      std::cout << " reduced too much: " << lastInit << " - "
                                << lastRed << " > "
                                << ThresholdCutList[thisPoi];
                      passThreshold = false;
                    } else {
                      std::cout << " tolerated: " << lastInit << " - "
                                << lastRed << " < "
                                << ThresholdCutList[thisPoi];
                    }
                  } else {
                    passThreshold = false;
                  }
                }
              }
              break;
            }
            thisPosition++;
          }

          tmpIndex++;

          std::cout << ", position " << thisPosition << " of " << tmpSum << " ("
                    << Form("%09f%%", thisReduction * 100) << ")" << std::endl;
        }
        // std::cout << std::endl;

        // All individual changes of the uncertainties on the POIs to combine
        // are below
        // the specified thresholds. Accept the found parameter and leave loop
        if (!passThreshold) {
          coutI(ObjectHandling)
              << "AbsMeasurement::PruneNuisanceParameters(" << fName
              << ")  └-> Vetoed parameter because uncertainty on one POI "
                 "reduced too much!"
              << std::endl;
        } else {
          coutI(ObjectHandling)
              << "AbsMeasurement::PruneNuisanceParameters(" << fName
              << ")  └-> Removing parameter!" << std::endl;
          foundPar2Rem = true;
          break;
        }
      }

      // If no parameter is found, break the pruning procedure
      if (!foundPar2Rem) {
        coutI(ObjectHandling)
            << "AbsMeasurement::PruneNuisanceParameters(" << fName
            << ")  └-> No more parameters found to prune!" << std::endl;
        break;
      }

      // Add the found parameter to the global ranking and remove it from
      // the list of parameters present in the covariance matrix
      order.push_back(par2rem);
      int par2rem_index = pars.index(par2rem.c_str());
      pars.remove(*pars.at(par2rem_index));
      uncertsPmO.insert(std::make_pair(par2rem_err, par2rem));
    }
  }

  // Print final ranking of non-filtered parameters based on the uncertainty on
  // combined inclusive POI
  coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                        << ") Ranking of all non-filtered parameters based on "
                           "combined inclusive POI"
                        << std::endl;
  PrintRanking(uncertsPmO, initTotalError);

  // Return the set of parameters that can be pruned safely
  return order;
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::RemoveParameter(
    TMatrixDSym &hes, RooArgList &pars, std::list<std::string> names) {
  // Remove a set of parameters from a matrix and thus reduce it
  // Find rows and columns to keep and remove
  std::set<int> removeRows;
  std::vector<int> keepRows;

  for (auto it = names.begin(); it != names.end(); ++it) {
    std::string name = *it;
    int index = pars.index(name.c_str());
    removeRows.insert(index);
  }

  RooArgList redpars(pars);
  for (RooLinkedListIter it = redpars.iterator();
       RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());) {
    std::string name = v->GetName();
    int index = redpars.index(name.c_str());
    if (removeRows.find(index) == removeRows.end()) {
      keepRows.push_back(index);
    } else {
      int remindex = pars.index(name.c_str());
      pars.remove(*pars.at(remindex));
    }
  }

  // Remove specified rows and columns from Hesse matrix
  int *a = &keepRows[0];
  TArrayI keepRow(keepRows.size(), a);

  for (Int_t i = 0; i < keepRow.GetSize(); i++) {
    TMatrixDColumn(hes, i) = TMatrixDColumn(hes, keepRow[i]);
    TMatrixDRow(hes, i) = TMatrixDRow(hes, keepRow[i]);
  }

  hes.ResizeTo(keepRow.GetSize(), keepRow.GetSize());
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::PrintRanking(
    std::set<std::pair<double, std::string>> uncerts, double initTotalError) {
  // Print ranking, which is stored in a (ordered) set
  for (std::set<std::pair<double, std::string>>::reverse_iterator itr =
           uncerts.rbegin();
       itr != uncerts.rend(); ++itr) {
    std::string varName(itr->second);
    double uncert(itr->first);

    double delta = uncert / initTotalError - 1;

    coutI(ObjectHandling) << "AbsMeasurement::PrintRanking(" << fName << ") "
                          << Form("%09f%% (%05f -> %05f): %s", delta * 100,
                                  initTotalError, uncert, varName.c_str())
                          << std::endl;
  }
  std::cout << std::endl;
}

// ____________________________________________________________________________|__________

std::pair<double, double>
RooFitUtils::AbsMeasurement::PDGrounding(double value, double error,
                                         int additionalDigit) {
  // Given a value and an error, round and format them according to the PDG
  // rules
  // for significant digits

  int threeDigits = GetThreeDigits(error);
  int nSignificantDigits = GetNSigDigits(threeDigits);
  nSignificantDigits += additionalDigit;

  // extraRound is meant for the special case of threeDigits > 950
  int extraRound;
  if (threeDigits >= 950)
    extraRound = 1;
  else
    extraRound = 0;

  // Convert to mantissa + exponent representation
  int expVal, expErr;
  frexp10(value, &expVal);
  frexp10(error, &expErr);

  // Format the value and error
  double formVal = FormatValue(
      value, expVal, expVal - expErr + nSignificantDigits, extraRound);
  double formErr = FormatValue(error, expErr, nSignificantDigits, extraRound);

  return std::make_pair(formVal, formErr);
}

// ____________________________________________________________________________|__________

int RooFitUtils::AbsMeasurement::GetThreeDigits(double error) {
  // Extract the three most significant digits and return them as an integer
  std::ostringstream stream;
  stream << Form("%.2e", error);
  TString str(stream.str());
  str = ((TObjString *)(str.Tokenize("e"))->At(0))->GetString();
  str.ReplaceAll(".", "");
  str.ReplaceAll("+", "");
  str.ReplaceAll("-", "");

  int threeDigits = atoi(str.Data());

  return threeDigits;
}

// ____________________________________________________________________________|__________

int RooFitUtils::AbsMeasurement::GetNSigDigits(int threeDigits) {
  // Find the number of significant digits
  assert(threeDigits < 1000);
  int nSignificantDigits;
  if (threeDigits < 101)
    nSignificantDigits = 2;
  else if (threeDigits < 356)
    nSignificantDigits = 2;
  else if (threeDigits < 950)
    nSignificantDigits = 1;
  else
    nSignificantDigits = 2;

  return nSignificantDigits;
}

// ____________________________________________________________________________|__________

double RooFitUtils::AbsMeasurement::frexp10(double x, int *exp) {
  // Convert a number to mantissa + exponent representation in base 10
  double mantissa = .0 > x ? -x : x;
  *exp = 0;

  if (mantissa >= 10.) {
    *exp = 1;
    for (; !((mantissa /= 10.) < 10.); ++(*exp))
      ;
  } else if (!(mantissa > 1.)) {
    *exp = -1;
    for (; !((mantissa *= 10.) > 1.); --(*exp))
      ;
  }

  return mantissa;
}

// ____________________________________________________________________________|__________

double RooFitUtils::AbsMeasurement::FormatValue(double value, int exponent,
                                                int nDigits, int extraRound) {
  // Format a value correctly and remove not needed digits
  int roundAt = nDigits - 1 - exponent - extraRound;
  int nDec;
  if (exponent < nDigits)
    nDec = roundAt;
  else
    nDec = 0;

  std::ostringstream stream;
  stream << "%." << nDec << "f";

  double tmp = pow(10, roundAt);
  double formVal = atof(Form(stream.str().c_str(), round(value * tmp) / tmp));

  return formVal;
}

namespace {
union MyFloat_t {
  MyFloat_t(float num = 0.0f) : f(num) {}

  // Portable extraction of components.
  bool Negative() const { return (i >> 31) != 0; }
  int32_t RawMantissa() const { return i & ((1 << 23) - 1); }
  int32_t RawExponent() const { return (i >> 23) & 0xFF; }

  int32_t i;
  float f;
};
}

bool RooFitUtils::AbsMeasurement::AlmostEqualUlpsAndAbs(float A, float B,
                                                        float maxDiff,
                                                        int maxUlpsDiff) {
  // Check if the numbers are really close -- needed  when comparing numbers
  // near zero.
  float absDiff = fabs(A - B);
  if (absDiff <= maxDiff)
    return true;

  MyFloat_t uA(A);
  MyFloat_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative())
    return false;

  // Find the difference in ULPs.
  int ulpsDiff = abs(uA.i - uB.i);
  if (ulpsDiff <= maxUlpsDiff)
    return true;

  return false;
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::writeToFile() {
  // Save the measurement
  coutP(ObjectHandling) << "AbsMeasurement::writeToFile(" << fName
                        << ") saving measurement as " << fFileName << std::endl;
  fWorkSpace->writeToFile(fFileName.c_str());
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::writeToFile(const char *fileName) {
  // Save the measurement
  fFileName = std::string(fileName);
  this->writeToFile();
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::SetDatasetBinning(
    Int_t setNbins, const char *generateBinnedTag, const char *binnedCategories,
    const char *unbinnedCategories, const char *weightVarName) {
  // Parse configuration for binning data
  fSetBinning = kTRUE;
  fSetNbins = setNbins;
  fGenerateBinnedTag = generateBinnedTag;
  fBinnedCategories = binnedCategories;
  fUnbinnedCategories = unbinnedCategories;
  fWeightVarName = weightVarName;
}

// ____________________________________________________________________________|__________

// Generate binned data
// Method adapted from Tim Adye <T.J.Adye@rl.ac.uk>
// Originally implemented in
// https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/HiggsPhys/CombinationTools/RooStatTools/trunk/StandardHypoTestInv/WorkspaceCalculator.cxx

ClassImp(TOwnedList)

    TOwnedList::TOwnedList()
    : TList() {
  SetOwner();
}
TOwnedList::~TOwnedList() { Clear(); }
void TOwnedList::Clear(Option_t *option) {
  if (!option || strcmp(option, "nodelete") != 0)
    for (TIter it(this); TObject *obj = it();) {
      // cout << "Delete
      // "<<obj->ClassName()<<"::"<<obj->GetName()<<(obj->IsOnHeap()?"":" (not
      // on heap)")<<std::endl;
      delete obj;
    }
  TList::Clear("nodelete");
}

// ____________________________________________________________________________|__________

RooDataSet *RooFitUtils::AbsMeasurement::SetDatasetBinning(
    const RooAbsPdf *pdfIn, const RooAbsData *data, Int_t setNbins,
    const char *generateBinnedTag, const char *binnedCategories,
    const char *unbinnedCategories, const char *weightVarName) {
  WildcardList binnedCategoriesWild(binnedCategories ? binnedCategories : "");
  WildcardList unbinnedCategoriesWild(unbinnedCategories ? unbinnedCategories
                                                         : "");
  TOwnedList localObjs;
  const RooSimultaneous *pdf = dynamic_cast<const RooSimultaneous *>(pdfIn);
  if (!pdf)
    return 0;
  const RooAbsCategoryLValue &cat = pdf->indexCat();
  TList *dataList = data->split(cat, true);
  if (!dataList)
    return 0;
  localObjs.Add(dataList);
  std::map<std::string, RooDataSet *> dataMap;
  RooRealVar weightVar(weightVarName, "", 1.0);
  int dsBinned = 0;
  for (TIter nextds = dataList;
       RooAbsData *datai = dynamic_cast<RooAbsData *>(nextds());) {
    localObjs.Add(datai); // make sure we delete them
    // Make a copy just so we can match the weight variable names
    TString copyName = Form("%s_unbinned", datai->GetName());
    RooDataSet *copyData =
        new RooDataSet(copyName, copyName, RooArgSet(*datai->get(), &weightVar),
                       weightVar.GetName());
    copyData->append(*(RooDataSet *)datai); // also works for RooDataHist
    dataMap[datai->GetName()] = copyData;
    localObjs.Add(copyData);
    int isBinned = -1;
    const char *dataType = "";
    if (RooAbsPdf *pdfi = pdf->getPdf(datai->GetName())) {
      if (RooArgSet *obs = pdfi->getObservables(*datai)) {
        isBinned = pdfi->isBinnedDistribution(*obs) ? 1 : 0;
        dataType = (isBinned ? " binned" : " unbinned");
        delete obs;
      }
    }
    coutI(ObjectHandling) << "AbsMeasurement::SetDatasetBinning(" << fName
                          << ") "
                          << Form("Category %s%s dataset has %d/%g entries",
                                  datai->GetName(), dataType,
                                  datai->numEntries(), datai->sumEntries())
                          << std::endl;
    if (datai->numEntries() <= 0)
      continue;
    bool rebin = false, genBinOnly = false;
    if (isBinned == 1 && unbinnedCategoriesWild.Match(datai->GetName())) {
      if (datai->numEntries() <= datai->sumEntries())
        continue;
      rebin = true;
    } else if (isBinned == 0 && binnedCategoriesWild.Match(datai->GetName())) {
      if (setNbins <= 0)
        continue;
      genBinOnly = (datai->numEntries() <= setNbins);
      if (genBinOnly && datai->sumEntries() <= setNbins)
        continue;
    } else {
      continue;
    }
    RooAbsPdf *pdfi = pdf->getPdf(datai->GetName());
    if (!pdfi)
      continue;
    RooArgSet *obs = pdfi->getObservables(*datai);
    if (!rebin && generateBinnedTag && *generateBinnedTag) {
      pdfi->setAttribute(generateBinnedTag);
    }
    coutP(ObjectHandling) << "AbsMeasurement::SetDatasetBinning(" << fName
                          << ") "
                          << Form("%s binning on%s dataset %s with %d/%g "
                                  "entries, PDF %s, variables",
                                  (rebin ? "Prune"
                                         : genBinOnly ? "Generate" : "Set"),
                                  dataType, datai->GetName(),
                                  datai->numEntries(), datai->sumEntries(),
                                  pdfi->GetName())
                          << std::endl;

    for (RooLinkedListIter it = obs->iterator(); TObject *o = it.Next();) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(o);
      if (!v)
        continue;
      coutI(ObjectHandling) << "AbsMeasurement::SetDatasetBinning(" << fName
                            << ")   " << v->GetName() << std::endl;
      if (!rebin)
        v->setBinning(RooBinning(setNbins, v->getMin(), v->getMax()));
    }
    std::cout << std::endl;
    if (genBinOnly) {
      delete obs;
      continue;
    }
    if (obs->getSize() != 1 && !rebin) {
      std::cerr << "PDF " << pdfi->GetName() << " has " << obs->getSize()
                << " observables - can't create binned dataset" << std::endl;
      delete obs;
      continue;
    }
    RooRealVar *obsVar = dynamic_cast<RooRealVar *>(obs->first());

    TString newName = Form("%s_binned", datai->GetName());
    RooDataSet *dataiNew = new RooDataSet(
        newName, newName, RooArgSet(*obs, weightVar), weightVar.GetName());
    localObjs.Add(dataiNew);
    if (rebin) {
      const RooAbsBinning &binning = obsVar->getBinning();
      coutI(ObjectHandling)
          << "AbsMeasurement::SetDatasetBinning(" << fName << ") "
          << Form("Observable %s bins (%d,%g,%g) contain %d/%g entries:",
                  obsVar->GetName(), binning.numBins(), binning.lowBound(),
                  binning.highBound(), datai->numEntries(), datai->sumEntries())
          << std::endl;
      for (int j = 0, n = datai->numEntries(); j < n; j++) {
        const RooArgSet *row = datai->get(j);
        Double_t w = datai->weight();
        if (w == 0)
          continue;
        coutI(ObjectHandling)
            << "AbsMeasurement::SetDatasetBinning(" << fName << ") "
            << " " << dynamic_cast<RooRealVar *>(row->find(*obsVar))->getVal()
            << ":" << w << std::endl;
        dataiNew->add(*row, w);
      }
      std::cout << std::endl;
      if (dataiNew->numEntries() == 0) {
        dataiNew->add(*datai->get(0),
                      0.0); // keep at least one entry (is this necessary?)
      }
      if (dataiNew->numEntries() >= datai->numEntries()) {
        delete obs;
        continue;
      }
      coutI(ObjectHandling)
          << "AbsMeasurement::SetDatasetBinning(" << fName << ") "
          << Form("Rebin dataset from %d/%g to %d/%g in %s",
                  datai->numEntries(), datai->sumEntries(),
                  dataiNew->numEntries(), dataiNew->sumEntries(),
                  dataiNew->GetName())
          << std::endl;
    } else {
      TString histName = Form("%s_hist", datai->GetName());
      TH1 *hist = datai->createHistogram(histName, *obsVar);
      coutI(ObjectHandling)
          << "AbsMeasurement::SetDatasetBinning(" << fName << ") "
          << Form("Created histogram %s(%d,%g,%g)", hist->GetName(),
                  hist->GetNbinsX(), hist->GetXaxis()->GetXmin(),
                  hist->GetXaxis()->GetXmax())
          << std::endl;
      for (int j = 1, n = hist->GetNbinsX(); j <= n; ++j) {
        Double_t x = hist->GetXaxis()->GetBinCenter(j),
                 y = hist->GetBinContent(j);
        if (y == 0.0)
          continue;
        obsVar->setVal(x);
        dataiNew->add(*obs, y);
      }
      if (dataiNew->numEntries() ==
          0) { // keep at least one entry (is this necessary?)
        obsVar->setVal(hist->GetXaxis()->GetBinCenter(1));
        dataiNew->add(*obs, 0.0);
      }
      coutI(ObjectHandling)
          << "AbsMeasurement::SetDatasetBinning(" << fName << ") "
          << Form("in dataset %s with %d/%g entries", dataiNew->GetName(),
                  dataiNew->numEntries(), dataiNew->sumEntries())
          << std::endl;
      delete hist;
    }
    delete obs;
    dataMap[datai->GetName()] = dataiNew; // replace datai
    dsBinned++;
  }

  if (dsBinned <= 0)
    return 0;

  RooCategory *catVar = dynamic_cast<RooCategory *>(data->get()->find(cat));
  if (!catVar) {
    std::cerr << "Dataset " << data->GetName()
              << " does not have an index variable '" << cat.GetName()
              << "' - can't create binned dataset" << std::endl;
    return 0;
  }

  RooArgSet newObs, *allObs = pdf->getObservables(*data);
  for (RooFIter it = data->get()->fwdIterator();
       const RooAbsArg *dsvar = it.next();) {
    if (RooAbsArg *v = allObs->find(*dsvar))
      newObs.add(*v);
  }
  delete allObs;
  TString name = Form("%s_binnedgg", data->GetName());
  RooDataSet *newData = new RooDataSet(
      name, name, RooArgSet(newObs, weightVar), RooFit::Index(*catVar),
      RooFit::Import(dataMap), RooFit::WeightVar(weightVar));

  coutP(ObjectHandling) << "AbsMeasurement::SetDatasetBinning(" << fName << ") "
                        << Form("Replace dataset %s (%d/%g entries) with "
                                "dataset %s (%d/%g entries)",
                                data->GetName(), data->numEntries(),
                                data->sumEntries(), newData->GetName(),
                                newData->numEntries(), newData->sumEntries())
                        << std::endl;

  return newData;
}
