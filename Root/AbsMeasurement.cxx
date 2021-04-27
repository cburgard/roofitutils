#include "RooFitUtils/AbsMeasurement.h"
#include "RooFitUtils/ExtendedMinimizer.h"
#include "TKey.h"
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

//  coutP(InputArguments) << "AbsMeasurement::AbsMeasurement(" << fName
//                        << ") created" << std::endl;
}

// ____________________________________________________________________________|__________

RooFitUtils::AbsMeasurement::~AbsMeasurement() {
  // Destructor
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::enablePruning(const std::string &poi,
                                                const std::string &filter,
                                                const std::string &percentage) {

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
  fPruningPercentage = percentage;
}

// ____________________________________________________________________________|__________

void RooFitUtils::AbsMeasurement::MakeConstSnapshot(
     std::string infilename, RooFitResult* fitresult, 
     std::string parameters, std::string outwsname,
     std::string outsnapshot){
   TFile* file = new TFile(infilename.c_str(),"READ");
   TIter next(file->GetListOfKeys()); TKey *key;
   RooWorkspace* ws = NULL;
    while ((key = (TKey*)next())) {
      TString keyname = key->GetName();
      if (!keyname.Contains("ProcessID") && !keyname.Contains("nllscan") && !keyname.Contains("fitResult"))
      ws = (RooWorkspace*) file->Get(keyname);
     }
   
   TString allParameters = parameters;
   allParameters.ReplaceAll(" ", "");
   TObjArray *allParametersArray = allParameters.Tokenize(",");
   unsigned int numParameters = allParametersArray->GetEntries();
   for (unsigned int itrPar = 0; itrPar < numParameters; ++itrPar) {
     TString varname = ((TObjString *)allParametersArray->At(itrPar))->GetString();
     auto fltpars = fitresult->floatParsFinal();
     RooRealVar* par = (RooRealVar*) fltpars.find(varname);
     ws->var(par->GetName())->setVal(par->getVal());
     ws->var(par->GetName())->setConstant(kTRUE);
   }
   ws->saveSnapshot(outsnapshot.c_str(), parameters.c_str());
   ws->writeToFile(outwsname.c_str(),kTRUE);
}

// ____________________________________________________________________________|__________

std::set<std::pair<double,std::string>> RooFitUtils::AbsMeasurement::OrderNuisanceParameters(
  const TMatrixDSym chesse, RooFitResult *fitresult, const std::string &poi,
  const std::string &filter, unsigned int nlo, unsigned nhi) {
  // Ranking/pruning function 
  // Two step Procedure -
  // 1. Rank the nuisance parameters according to a criteria 
  //    current criteria is,
  //    rank(NP) = max. (fractional change in variance due NP) over all POIs
  // 

  // Decompose POI names and associated weights and thresholds
  TString allPoi = poi;
  TString allFilter = filter;
  allPoi.ReplaceAll(" ", "");
  allFilter.ReplaceAll(" ", "");
  TObjArray *allPoiArray = allPoi.Tokenize(",");
  TObjArray *allFilterArray = allFilter.Tokenize(",");
  unsigned int numPoi = allPoiArray->GetEntries();
  unsigned int numFilter = allFilterArray->GetEntries();

  // Store information in easier usable objects
  std::vector<TString> PoiList;
  std::vector<TString> FilterList;

  // Loop over the POIs to check for pruning and collect assigned percentages
  for (unsigned int itrPoi = 0; itrPoi < numPoi; ++itrPoi) {
    TString thisPoi = ((TObjString *)allPoiArray->At(itrPoi))->GetString();
    PoiList.push_back(thisPoi);
    coutP(ObjectHandling) << "AbsMeasurement::OrderNuisanceParameters(" << fName
                          << ") Parsed POI " << thisPoi.Data()
                          << std::endl;
  }

  // Check if there are filters for checking the nuisance parameters to prune
  // in a particular order
  if (numFilter > 0) {
    coutP(ObjectHandling) << "AbsMeasurement::OrderNuisanceParameters(" << fName
			  << ") Nuisance parameters will be pruned in the following order:";

    for (unsigned int itrSet = 0; itrSet < numFilter; ++itrSet) {
      TString thisSet = ((TObjString *)allFilterArray->At(itrSet))->GetString();
      FilterList.push_back(thisSet);
      std::cout << " " << thisSet.Data() << ";";
    }

  }
  // If none are provide choose all the nuisance paramter in the model to check 
  // for pruning
  else {
    coutP(ObjectHandling) << "AbsMeasurement::OrderNuisanceParameters(" << fName
                          << ") Include all parameters in ranking" << std::endl;
    FilterList.push_back(".*");
  }
  std::cout << std::endl;
  // Clean up
  delete allPoiArray;
  delete allFilterArray;

  // Print the RooFitResult used for pruning
  // fitresult->Print();

  // Get initial information, such as
  //   - covariance matrix
  //   - floating parameters
  //   - position of POIs in the list of floating parameters
  //   - best fit and initial Hesse error
  TMatrixDSym hes(chesse);
  TMatrixDSym cov(hes.Invert());
  RooArgList pars = fitresult->floatParsFinal();

  // maps to keep track of the position, initial error, best fit value
  // and the sign of the init error
  std::map<std::string, int> index;
  std::map<std::string, double> initErr;
  std::map<std::string, double> bestFit;
  std::map<std::string, double> initErrSig;

  // loop over the pois to prune and extract the necessary info.
  for (std::vector<TString>::const_iterator itr = PoiList.begin(),
                                            end = PoiList.end();
       itr != end; ++itr) {
    int thisIndex = pars.index(itr->Data());
    double thisInitErr = sqrt(cov[thisIndex][thisIndex]);
    double thisBestFit = ((RooRealVar *)pars.at(thisIndex))->getVal();
    std::pair<double, double> tmpVals =
        std::make_pair(thisBestFit, thisInitErr);

    index[itr->Data()] = thisIndex;
    initErr[itr->Data()] = thisInitErr;
    initErrSig[itr->Data()] = tmpVals.second;
    bestFit[itr->Data()] = thisBestFit;

    coutP(ObjectHandling) << "AbsMeasurement::OrderNuisanceParameters(" << fName
                          << ") Harvest parameter (" << thisIndex
                          << "): " << *itr << " "
                          << tmpVals.first << " +/- "
                          << tmpVals.second << std::endl;
  }

  // Ranking step starts here.
  // Iteratively remove parameters and recompute Hesse uncertainty
  // on the POI and obtain the max. of their fractional changes
  // final ranking.
  // Maybe filter list of parameters to prune base on regexp.


  // ranks is a set of the NP and their estimated rank rank
  std::set<std::pair<double, std::string>> ranks;
  unsigned int tmpSum = pars.getSize() + 1 - numPoi;
  // Start looping over the different NP sets if any order is mentioned
  // if everything is looked at in one go
  for (std::vector<TString>::const_iterator itrFilter = FilterList.begin(),
                                            end = FilterList.end();
       itrFilter != end; ++itrFilter) {
    TString thisFilter = *itrFilter;
    std::string stringFilter = thisFilter.Data();

    TRegexp reg(thisFilter);
    Ssiz_t dummy(0);

    // loop over all the parameters that are floating in the final fit
    // and find their rank 
    unsigned int itrpos = 0;
    nhi = (nhi > 0) ? nhi : tmpSum;
    coutP(ObjectHandling) << "AbsMeasurement::OrderNuisanceParameters(" << fName
			  << ") Estimating ranks for " << nhi - nlo + 1 
                          << " nuisance parameters." << std::endl;

    for (RooLinkedListIter it = pars.iterator();
         RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());) {
      std::string name = v->GetName();
      ++itrpos;
      if (itrpos  <= nlo || itrpos > nhi + 1) continue;
      // Skip removing POIs
      if (allPoi.Contains(name) || (stringFilter != "" &&
          reg.Index(TString(v->GetName()), &dummy, 0) == -1)) continue;

      std::list<std::string> names;
     
      // make and invert reduced hesse matrix used for ranking 
      names.push_back(name);
      TMatrixDSym tmp_hesse(chesse);
      RooArgList tmp_pars = fitresult->floatParsFinal();
      RemoveParameter(tmp_hesse, tmp_pars, names);
      TMatrixDSym tmp_cov = tmp_hesse.Invert();

      // calculating the rank of the nuisance parameter
      // -- find the absolute fractional change in 
      //    variance relative to the intial error 
      //    for all POIs and then choose
      //    the maximum from this group

      std::vector<double> poisfracErr;
      std::map<std::string, int> index_red;
      for (std::vector<TString>::const_iterator itr = PoiList.begin(),
           end = PoiList.end(); itr != end; ++itr) {
 
        int thisIndex = tmp_pars.index(itr->Data());
	double thisErr = sqrt(tmp_cov[thisIndex][thisIndex]);
 	index_red[itr->Data()] = thisIndex;

	double diff = initErr[itr->Data()] - thisErr;
	double fracErr = diff/(initErr[itr->Data()]);
	poisfracErr.push_back(std::abs(fracErr));
      }

      // find the maximum amongst the group of absolute fractional variance changes
      double max_val = *std::max_element(std::begin(poisfracErr), std::end(poisfracErr));
      ranks.insert(std::make_pair(max_val, name));
    }
  }
    coutP(ObjectHandling) << "AbsMeasurement::OrderNuisanceParameters(" << fName
			  << ") Ranked " << ranks.size() 
                          << " nuisance parameters." << std::endl;
    return ranks;
}

// ____________________________________________________________________________|__________

std::set<std::string> RooFitUtils::AbsMeasurement::PruneNuisanceParameters(
  std::set<std::pair<double,std::string>> ranks,
  const TMatrixDSym chesse, RooFitResult *fitresult,
  const std::string &poi,const std::string &percentage, const std::string &filter) {

// 2. Once ranked find the largest ranked nuisance parameter 
//    by performing a binary search
//    -- if necessary the remain NPs can be also checked using checkAll option

  // Decompose POI names and associated weights and thresholds
  TString allPoi = poi;
  TString allPercentage = percentage;
  TString allFilter = filter;
  allPoi.ReplaceAll(" ", "");
  allFilter.ReplaceAll(" ", "");
  allPercentage.ReplaceAll(" ", "");
  TObjArray *allPoiArray = allPoi.Tokenize(",");
  TObjArray *allFilterArray = allFilter.Tokenize(",");
  TObjArray *allPercentageArray = allPercentage.Tokenize(",");
  unsigned int numPoi = allPoiArray->GetEntries();
  unsigned int numFilter = allFilterArray->GetEntries();
  unsigned int numPercentage = allPercentageArray->GetEntries();

  // Check if number of POIs and percentages match
  if (numPoi != numPercentage && numPercentage != 1) {
    coutE(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Number of POIs does not match number of "
                             "specified threshold percentages. No additional pruning."
                          << std::endl;
    return std::set<std::string>();
  }
  // Store information in easier usable objects
  std::vector<TString> PoiList;
  std::map<std::string, double> PercentageList;
  std::vector<TString> FilterList;

  // Loop over the POIs to check for pruning and collect assigned percentages
  for (unsigned int itrPoi = 0; itrPoi < numPoi; ++itrPoi) {
    TString thisPoi = ((TObjString *)allPoiArray->At(itrPoi))->GetString();
    TString thisPercentage =
      	    (numPercentage == 1)
     	    ? ((TObjString *)allPercentageArray->At(0))->GetString()
            : ((TObjString *)allPercentageArray->At(itrPoi))->GetString();

    double tmpPercentage = atof(thisPercentage.Data());

    PoiList.push_back(thisPoi);
    PercentageList[thisPoi.Data()] = 0.01*tmpPercentage;
    coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Parsed POI " << thisPoi.Data()
                          << ", assigning percentage " << tmpPercentage << " %"
                          << std::endl;
  }

  // Check if there are filters for checking the nuisance parameters to prune
  // in a particular order
  if (numFilter > 0) {
    coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
              << ") Nuisance parameters will be pruned in the following order:";

    for (unsigned int itrSet = 0; itrSet < numFilter; ++itrSet) {
      TString thisSet = ((TObjString *)allFilterArray->At(itrSet))->GetString();
      FilterList.push_back(thisSet);
      std::cout << " " << thisSet.Data() << ";";
    }

  }
  // If none are provide choose all the nuisance paramter in the model to check 
  // for pruning
  else {
   coutI(ObjectHandling)  << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Include all parameters in ranking" << std::endl;
    FilterList.push_back(".*");
  }
  std::cout << std::endl;
  // Clean up
  delete allPoiArray;
  delete allPercentageArray;
  delete allFilterArray;

  // Print the RooFitResult used for pruning
  // fitresult->Print();

  // Get initial information, such as
  //   - covariance matrix
  //   - floating parameters
  //   - position of POIs in the list of floating parameters
  //   - best fit and initial Hesse error
  TMatrixDSym hes(chesse);
  TMatrixDSym cov(hes.Invert());
  RooArgList pars = fitresult->floatParsFinal();

  // maps to keep track of the position, initial error, best fit value
  // and the sign of the init error
  std::map<std::string, int> index;
  std::map<std::string, double> initErr;
  std::map<std::string, double> bestFit;
  std::map<std::string, double> initErrSig;

  // loop over the pois to prune and extract the necessary info.
  for (std::vector<TString>::const_iterator itr = PoiList.begin(),
                                            end = PoiList.end();
       itr != end; ++itr) {
    int thisIndex = pars.index(itr->Data());
    double thisInitErr = sqrt(cov[thisIndex][thisIndex]);
    double thisBestFit = ((RooRealVar *)pars.at(thisIndex))->getVal();
    std::pair<double, double> tmpVals =
        std::make_pair(thisBestFit, thisInitErr);

    index[itr->Data()] = thisIndex;
    initErr[itr->Data()] = thisInitErr;
    initErrSig[itr->Data()] = tmpVals.second;
    bestFit[itr->Data()] = thisBestFit;

    coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                          << ") Harvest parameter (" << thisIndex
                          << "): " << *itr << " "
                          << tmpVals.first << " +/- "
                          << tmpVals.second << std::endl;
  }

  // Ranking step starts here.
  // Iteratively remove parameters and recompute Hesse uncertainty
  // on the POI and obtain the max. of their fractional changes
  // final ranking.
  // Maybe filter list of parameters to prune base on regexp.

  std::set<std::string> finalNPnames;
  // Start looping over the different NP sets if any order is mentioned
  // if everything is looked at in one go
  for (std::vector<TString>::const_iterator itrFilter = FilterList.begin(),
                                            end = FilterList.end();
       itrFilter != end; ++itrFilter) {
    TString thisFilter = *itrFilter;
    std::string stringFilter = thisFilter.Data();

    TRegexp reg(thisFilter);

    // ranks is a set of the NP and their estimated rank rank
    // part_uncerts is the variance change due to the NP for all POIs
    if (ranks.empty()) break;

    // Print ranking of parameters
    coutI(ObjectHandling)
              << "AbsMeasurement::PruneNuisanceParameters(" << fName
              << ") Ranking of parameters"
              << std::endl;
    PrintRanking(ranks);

    bool foundPar2Rem = false;
    std::string par2rem = ranks.begin()->second;

    unsigned int hipos = ranks.size();
    unsigned int lopos = 0;
    // Loop over parameters in the global ranking
    for (auto rankitr:ranks){
      std::set<std::string> pruneNPnames;
      par2rem = rankitr.second;
      bool passThreshold = true;
      int tmpIndex = 0;

      coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                            << ") --> Testing " << par2rem << std::endl;
 
      std::map<std::string, std::set<std::pair<double, std::string>>> pois_uncerts;

      TMatrixDSym tmp_hesse(chesse);
      RooArgList tmp_pars = fitresult->floatParsFinal();

      auto x = ranks.begin();
      int pos = (hipos + lopos)/2;
      for(int j = 0; j < pos; ++j){
        pruneNPnames.insert(x->second);
        ++x;
      }
      std::list<std::string> nuis_names;
      for (RooLinkedListIter it = tmp_pars.iterator();
        RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());){
        std::string varname = v->GetName();
        bool found = (std::find(pruneNPnames.begin(), pruneNPnames.end(), varname) != pruneNPnames.end());
        if (found) nuis_names.push_back(varname);
        else{if (allPoi.Contains(varname)) continue;
        if (varname == rankitr.second) continue;}
      } 

      RemoveParameter(tmp_hesse, tmp_pars, nuis_names);
      TMatrixDSym tmp_cov = tmp_hesse.Invert();
      std::map<std::string, int> index_red;
      for (std::vector<TString>::const_iterator itr = PoiList.begin(), end = PoiList.end(); itr != end; ++itr) {
        int thisIndex = tmp_pars.index(itr->Data());
        double thisErr = sqrt(tmp_cov[thisIndex][thisIndex]);
        index_red[itr->Data()] = thisIndex;
        pois_uncerts[itr->Data()].insert(std::make_pair(thisErr, par2rem));
      }

      // Loop over the individual POIs to test effect of removing proposed
      // parameter on each of them
      for (auto iterator = pois_uncerts.begin();
             iterator != pois_uncerts.end(); ++iterator) {
        std::string thisPoi = iterator->first;
       
        coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters("
       	               << fName << ")  | " << thisPoi;
       
        std::set<std::pair<double, std::string>> thisRanking = iterator->second;
        unsigned int thisPosition = 1;
        double thisReduction = 1.0;
        double thisInitError = initErr[thisPoi];
       
        // Find the proposed parameter, the position in the individual ranking
        // and the estimated uncertainty reduction
        for (std::set<std::pair<double, std::string>>::reverse_iterator
              subitr = thisRanking.rbegin(); subitr != thisRanking.rend(); ++subitr) {
           std::string varName(subitr->second);
          double thisUncert = subitr->first;
        
          // Check the individual uncertainty reduction if the parameter is
          // found and exit search
          if (varName == par2rem) {
            thisReduction = thisUncert / thisInitError - 1;
        
            // Check if the threshold is passed, either the relative reduction
            // or check
            // whether quoted uncertainty would be changed in case of 'auto'
            std::pair<double, double> tmpVals = std::make_pair(bestFit[thisPoi], thisUncert);
        
            std::cout << ": absolute change " << initErrSig[thisPoi]
                      << " --> " << tmpVals.second;
        
            // Allow the specified non-significant digit to change by a
            // specified value
            if (tmpVals.second < initErrSig[thisPoi]) {
              if (!AlmostEqualUlpsAndAbs(initErrSig[thisPoi], tmpVals.second,
        			         PercentageList[thisPoi]*initErrSig[thisPoi], 4)) {
                std::cout << " reduced too much: " << initErrSig[thisPoi] << " - "
        	 	  << tmpVals.second << " > "
        	 	  << PercentageList[thisPoi]*initErrSig[thisPoi];
        	passThreshold = false;
              } else {
                std::cout << " tolerated: " << initErrSig[thisPoi] << " - "
        	  	  << tmpVals.second << " < "
                          << PercentageList[thisPoi]*initErrSig[thisPoi];
              }
            }
            break;
          }
          thisPosition++;
        }
        tmpIndex++;
        std::cout << ", ("
	      << Form("%09f%%", thisReduction * 100) << ")" << std::endl;
      }

      // All individual changes of the uncertainties on the POIs to check
      // are below
      // the specified thresholds. Accept the found parameter and leave loop
      if (!passThreshold) {
          coutI(ObjectHandling)
          << "AbsMeasurement::PruneNuisanceParameters(" << fName
          << ")  └-> Vetoed parameter because uncertainty on one POI "
             "reduced too much!"
          << std::endl;
      hipos = pos - 1;
      } else {
          coutI(ObjectHandling)
          << "AbsMeasurement::PruneNuisanceParameters(" << fName
          << ")  └-> Removing parameter "<< par2rem << std::endl;
          lopos = pos + 1;
          foundPar2Rem = true;
          pruneNPnames.insert(par2rem);
      }

      // If no parameter is found, break the pruning procedure
      if (!foundPar2Rem || lopos > hipos) {
        coutI(ObjectHandling)
        << "AbsMeasurement::PruneNuisanceParameters(" << fName
        << ")  └-> No more parameters found to prune!" << std::endl;
      
        coutI(ObjectHandling)
        << "AbsMeasurement::PruneNuisanceParameters(" << fName
        << ") Nuisance Parameters identified for the following criteria, " << std::endl;
        for (std::vector<TString>::const_iterator itr = PoiList.begin(),
                                                  end = PoiList.end();
             itr != end; ++itr) {
          int thisIndex = pars.index(itr->Data());
          double thisInitErr = sqrt(cov[thisIndex][thisIndex]);
          double thisBestFit = ((RooRealVar *)pars.at(thisIndex))->getVal();
          std::pair<double, double> tmpVals =
              std::make_pair(thisBestFit, thisInitErr);
      
          index[itr->Data()] = thisIndex;
          initErr[itr->Data()] = thisInitErr;
          initErrSig[itr->Data()] = tmpVals.second;
          bestFit[itr->Data()] = thisBestFit;
      
          coutI(ObjectHandling) << "AbsMeasurement::PruneNuisanceParameters(" << fName
                                << ") (" << thisIndex
                                << " ): " << *itr << " " << tmpVals.first << " +/- "
                                << tmpVals.second << " "
			        << " for a fractional variation of "
                                << 100*PercentageList[itr->Data()] <<" %"
                                << std::endl;
        }
     
        for (auto v : pruneNPnames){
           finalNPnames.insert(v);
        }
        break;
      }
      // Add the found parameter to the global ranking and remove it from
      // the list of parameters present in the covariance matrix
      int par2rem_index = pars.index(par2rem.c_str());
      pars.remove(*pars.at(par2rem_index));
    }
  }
  // Return the set of parameters that can be pruned safely
  return finalNPnames;
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
    std::set<std::pair<double, std::string>> uncerts) {
  // Print ranking, which is stored in a (ordered) set
  for (std::set<std::pair<double, std::string>>::reverse_iterator itr =
	   uncerts.rbegin();
       itr != uncerts.rend(); ++itr) {
    std::string varName(itr->second);
    double uncert(itr->first);

    coutI(ObjectHandling) << "AbsMeasurement::PrintRanking(" << fName << ") "
			  << Form(" (%03f %%): %s", 100*uncert, varName.c_str())
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
