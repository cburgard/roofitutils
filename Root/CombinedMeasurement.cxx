#include "RooFitUtils/CombinedMeasurement.h"
#include "RooFitUtils/ExtendedMinimizer.h"

#include "RooGaussian.h"
#include "RooLognormal.h"
#include "RooPoisson.h"
#include "RooGamma.h"
#include "RooBifurGauss.h"

// ____________________________________________________________________________|__________
// Constructor
RooFitUtils::CombinedMeasurement::CombinedMeasurement( const std::string& CombinedMeasurementName, const std::string& WorkspaceName, const std::string& ModelConfigName, const std::string& DataName )
  :
  AbsMeasurement( CombinedMeasurementName, WorkspaceName, ModelConfigName, DataName )
{
  coutP(InputArguments) << "CombinedMeasurement::CombinedMeasurement(" << fName <<") created" << endl;

  fPdfName = "combPdf";
  fCategoryName = "master_measurement";

  fWorkSpace = new RooWorkspace(fWorkspaceName.c_str());
  fWorkSpace->autoImportClassCode(kTRUE);
  fModelConfig = new ModelConfig(fModelConfigName.c_str(), fWorkSpace);
  fParametersOfInterest = new RooArgSet("ParametersOfInterest");
  fNuisanceParameters = new RooArgSet("NuisanceParameters");
  fObservables = new RooArgSet("Observables");
  fGlobalObservables = new RooArgSet("GlobalObservables");
  fAsimovData = NULL;
}

// ____________________________________________________________________________|__________
// Constructor to load existing (regularised measurement)
RooFitUtils::CombinedMeasurement::CombinedMeasurement( const std::string& CombinedMeasurementName, const std::string& FileName, const std::string& WorkspaceName, const std::string& ModelConfigName, const std::string& DataName )
  :
  AbsMeasurement( CombinedMeasurementName, FileName, WorkspaceName, ModelConfigName, DataName )
{
  fFile = TFile::Open(fFileName.c_str());
  fWorkSpace = (RooWorkspace*) fFile->Get(fWorkspaceName.c_str());
  fModelConfig = (ModelConfig*)fWorkSpace->obj(fModelConfigName.c_str());
  fData = (RooAbsData*) fWorkSpace->obj(fDataName.c_str());
  fPdf = (RooSimultaneous*)fModelConfig->GetPdf();
  fParametersOfInterest = (RooArgSet*) fModelConfig->GetParametersOfInterest();
  fNuisanceParameters = (RooArgSet*) fModelConfig->GetNuisanceParameters();
  fObservables = (RooArgSet*) fModelConfig->GetObservables();
  fGlobalObservables = (RooArgSet*) fModelConfig->GetGlobalObservables();
  fAsimovData = NULL;

  fParametersOfInterestString = "";
  TIterator* poiItr = fParametersOfInterest->createIterator();
  RooRealVar* nextPOI;
  while ((nextPOI = (RooRealVar*)poiItr->Next())) {
    fParametersOfInterestString += "," + string(nextPOI->GetName());
  }

  fPdfName = fPdf->GetName();
  fCategoryName = ((RooSimultaneous*)fPdf)->indexCat().GetName();
}

// ____________________________________________________________________________|__________
// Destructor
RooFitUtils::CombinedMeasurement::~CombinedMeasurement()
{
  coutI(InputArguments) << "CombinedMeasurement::~CombinedMeasurement(" << fName << ") cleaning up" << endl;

  fNuisanceParameters->Delete();
  delete fNuisanceParameters;

  fObservables->Delete();
  delete fObservables;

  fGlobalObservables->Delete();
  delete fGlobalObservables;

  fParametersOfInterest->Delete();
  delete fParametersOfInterest;
}

// ____________________________________________________________________________|__________
// Initialisation
void RooFitUtils::CombinedMeasurement::initialise()
{
  fIsInitialised = kTRUE;
}

// ____________________________________________________________________________|__________
// Collect all measurement after regularising their individual channels
void RooFitUtils::CombinedMeasurement::CollectMeasurements()
{
  // initialise the measurements and get nuisance parameters from individual measurements
  // to possibly add to correlation scheme
  std::map< std::string, RooArgSet* > tmpAllNuisanceParameters;
  for (auto measItr = fMeasurements.begin(); measItr != fMeasurements.end(); ++measItr) {
    Measurement* meas = measItr->second;
    meas->initialise();
    RooArgSet* thisNuisanceParameters = meas->GetNuisanceParameters();
    tmpAllNuisanceParameters[meas->GetName()] = thisNuisanceParameters;
  }

  DetermineAutoCorrelations(tmpAllNuisanceParameters);

  // actual collection of measurements and channels which includes regularisation
  for (std::map< std::string, Measurement* >::iterator measItr = fMeasurements.begin(); measItr != fMeasurements.end(); ++measItr) {
    Measurement* meas = measItr->second;
    meas->SetRenamingMap(fCorrelationScheme->GetRenamingMap(meas->GetName()));
    meas->SetCorrelationFactors(fCorrelationScheme->GetCorrelationFactors(meas->GetName()));
    meas->CollectChannels();
    fCorrelationScheme->SetRenamingMap(meas->GetName(), meas->GetRenamingMap());
  }
}

// ____________________________________________________________________________|__________
// Add missing parameters to the correlation scheme. Maybe correlate parameters
// automatically, if requested from the CorrelationSchme. Use this feature with
// caution, as meaning of parameters and thus correlation might be different
// than the name suggests!
// Current implementation works only for Gaussian and Poisson constrained
// nuisance parameters.
// TODO implementation for other constraint types.
void RooFitUtils::CombinedMeasurement::DetermineAutoCorrelations( std::map< std::string, RooArgSet* >& tmpAllNuisanceParameters )
{
  bool enableAutoCorrelation = fCorrelationScheme->GetAutoCorrelation();
  std::map< std::string, std::set< std::string > > tmpAllNuisanceParameterMap;

  // Loop over all sets, containing the nuisance parameters of a single measurement
  for (std::map< std::string, RooArgSet* >::iterator it = tmpAllNuisanceParameters.begin(); it != tmpAllNuisanceParameters.end(); ++it) {
    std::string thisMeasurementName = it->first;
    RooArgSet* thisNuisanceParameters = it->second;

    // Get the existing RenamingMap for this Measurement
    RenamingMap thisRenamingMap = fCorrelationScheme->GetRenamingMap(thisMeasurementName);
    std::map< std::string, std::string > tmpRenamingMap = thisRenamingMap.GetRenamingMap();

    // Loop over all nuisance parameters in the set
    TIterator* nuiItr = thisNuisanceParameters->createIterator();
    RooRealVar* nextNui;
    while ((nextNui = (RooRealVar*)nuiItr->Next())) {
      std::string nextNuiName = nextNui->GetName();

      // Skip the parameter if it was already added explicitly
      if (tmpRenamingMap.find(nextNuiName) != tmpRenamingMap.end()) {
        continue;
      }

      // Skip all non-Gaussian and non-Poisson constrained parameters
      // TODO loosen requirement once proper handling of other constraint types is implemented in RenamingMap, etc.
      RenamingMap::ConstraintType thisConstraintType = fMeasurements[thisMeasurementName]->DetermineConstraintType(nextNui);
      if (thisConstraintType != RenamingMap::Gaussian && thisConstraintType != RenamingMap::Poisson) {
        continue;
      }

      // Add them to the temporary map used for adding them to the CorrelationScheme later
      tmpAllNuisanceParameterMap[nextNuiName].insert(thisMeasurementName);
    }
    delete nuiItr;
  }

  // Loop over all nuisance parameters not yet added to the CorrelationScheme and identify
  // the ones appearing potentially in more than one channel.
  // NEED TO BE CAREFUL, AS ONLY MATCHING BY NAME IMPLEMENTED HERE
  // TODO implement smarter matching algorithm, taking into account ranges, global observables,
  // constraint types, ...
  for (std::map< std::string, std::set< std::string > >::iterator it = tmpAllNuisanceParameterMap.begin(); it != tmpAllNuisanceParameterMap.end(); ++it) {
    std::string thisParameterName = it->first;
    std::set< std::string > thisMeasurements = it->second;

    if (thisMeasurements.size() == 1) {
      // In case a parameter shows up only in a single measurement (and was not yet added to the CorrelationScheme),
      // add it, but make sure that it can be identified later on as automatically added and the channel it originates from
      std::string thisMeasurementName = *thisMeasurements.begin();
      fCorrelationScheme->RenameParameter(thisMeasurementName.c_str(), thisParameterName.c_str(), ("auto_" + thisParameterName + "_" + thisMeasurementName).c_str());
    } else {
      // If a parameter is part of more than one measurement, either correlate it or de-correlate it, depending on
      // specified setting of CorrelationScheme. In case they should be correlated, a correlation string is composed
      // from the measurement name and parameter name. Parameters can be identified as automatically correlated by a prefix
      if (!enableAutoCorrelation) {
        for (std::set< std::string >::iterator measItr = thisMeasurements.begin(); measItr != thisMeasurements.end(); ++measItr) {
          std::string thisMeasurementName = *measItr;
          fCorrelationScheme->RenameParameter(thisMeasurementName.c_str(), thisParameterName.c_str(), ("auto_" + thisParameterName + "_" + thisMeasurementName).c_str());
        }
      } else {
        std::string correlationString;
        for (std::set< std::string >::iterator measItr = thisMeasurements.begin(); measItr != thisMeasurements.end(); ++measItr) {
          correlationString += *measItr + "::" + thisParameterName + ",";
        }
        correlationString.erase(correlationString.size()-1);
        fCorrelationScheme->CorrelateParameter(correlationString.c_str(), thisParameterName.c_str());
      }
    }
  }
}

// ____________________________________________________________________________|__________
// Build a new combined workspace from the regularised channels and measurements
void RooFitUtils::CombinedMeasurement::CombineMeasurements()
{
  std::map< std::string, RooDataSet* > datasetMap;
  std::map< std::string, RooAbsPdf* > pdfMap;
  coutP(ObjectHandling) << "CombinedMeasurement::CombineMeasurements(" << fName << ") Making new master category " << fCategoryName << endl;
  RooCategory* category = new RooCategory(fCategoryName.c_str(), fCategoryName.c_str());

  int nrMeasurements = fMeasurements.size();
  if (nrMeasurements == 0) {
    coutF(InputArguments) << "CombinedMeasurement::CombineMeasurements(" << fName << ") no measurements for combination specified" << endl;
    exit(-1);
  } else {
    coutP(InputArguments) << "CombinedMeasurement::CombineMeasurements(" << fName << ") will combine " << nrMeasurements << " measurements" << endl;
  }

  int numTotPdf = 0;

  for (std::map< std::string, Measurement* >::iterator measItr = fMeasurements.begin(); measItr != fMeasurements.end(); ++measItr) {
    Measurement* thisMeasurement = measItr->second;
    std::vector< Channel* > thisChannels = thisMeasurement->GetChannels();

    RooArgSet* thisObservables = thisMeasurement->GetObservables();
    RooArgSet* thisGlobalObservables = thisMeasurement->GetGlobalObservables();
    RooArgSet* thisNuisanceParameters = thisMeasurement->GetNuisanceParameters();

    for (RooLinkedListIter it = thisObservables->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      fObservables->add(*v);
    }

    for (RooLinkedListIter it = thisGlobalObservables->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      TString thisName = v->GetName();
      if (thisName.BeginsWith("PRUNED_")) {
        continue;
      }
      fGlobalObservables->add(*v);
    }
    fGlobalObservables->setAttribAll("GLOBAL_OBSERVABLE");

    for (RooLinkedListIter it = thisNuisanceParameters->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      TString thisName = v->GetName();
      if (thisName.BeginsWith("PRUNED_")) {
        continue;
      }
      fNuisanceParameters->add(*v);
    }

    int nrChannels = thisChannels.size();

    for (int itrChan = 0; itrChan < nrChannels; ++itrChan) {
      Channel* thisChannel = thisChannels[itrChan];
      RooAbsPdf* pdf_tmp = (RooAbsPdf*)thisChannel->GetPdf();
      RooAbsData* data_tmp = thisChannel->GetData();

      coutP(ObjectHandling) << "CombinedMeasurement::CombineMeasurements(" << fName << ") defining category " << thisChannel->GetName() << endl;
      coutP(ObjectHandling) << "CombinedMeasurement::CombineMeasurements(" << fName << ") with PDF " << pdf_tmp->GetName() << " and dataset " << data_tmp->GetName() << endl;
      category->defineType(thisChannel->GetName(),numTotPdf);
      category->setLabel(thisChannel->GetName(), true);

      pdfMap[thisChannel->GetName()] = pdf_tmp;
      datasetMap[thisChannel->GetName()] = (RooDataSet*)data_tmp;

      numTotPdf++;
    }
  }
  fObservables->add(*category);

  coutP(ObjectHandling) << "CombinedMeasurement::CombineMeasurements(" << fName << ") Making new simultaneous pdf " << fPdfName << endl;
  fPdf = new RooSimultaneous( fPdfName.c_str(), fPdfName.c_str(), pdfMap, *category );
  fPdf->setStringAttribute( "DefaultGlobalObservablesTag","GLOBAL_OBSERVABLE" );
  RooRealVar weightVar( "weightVar", "", 1., -1e10, 1e10 );

  // making the combined dataset
  coutP(ObjectHandling) << "CombinedMeasurement::CombineMeasurements(" << fName << ") Making combined dataset " << fDataName << endl;
  RooArgSet obs_cat_weight;
  obs_cat_weight.add(*fObservables);
  obs_cat_weight.add(*category);
  obs_cat_weight.add(weightVar);
  fData = new RooDataSet(fDataName.c_str(), fDataName.c_str(), obs_cat_weight, Index(*category), Import(datasetMap), WeightVar("weightVar"));

  // putting everything together
  MakeCleanWorkspace();

  for (std::map< std::string, RooDataSet* >::iterator it = datasetMap.begin(); it != datasetMap.end(); ++it) {
    delete it->second;
  }

  for (std::map< std::string, RooAbsPdf* >::iterator it = pdfMap.begin(); it != pdfMap.end(); ++it) {
    delete it->second;
  }
}

// ____________________________________________________________________________|__________
// Re-parametrise a (combined) measurement
void RooFitUtils::CombinedMeasurement::ParametriseMeasurements()
{
  coutP(ObjectHandling) << "CombinedMeasurement::ParametriseMeasurements() trying to re-parametrise the model according to " << fParametrisationSequence->GetName() << endl;

  // Get the list of schemes for this model from the parametrisation sequence
  std::list< ParametrisationScheme > thisSequence = fParametrisationSequence->GetSequence();

  // Loop through the schemes
  std::string originalName = fPdf->GetName();

  for (std::list< ParametrisationScheme >::const_iterator itrSequence = thisSequence.begin(), end = thisSequence.end(); itrSequence != end; ++itrSequence) {
    // Get the list of expressions from the parametrisation scheme
    ParametrisationScheme thisScheme = *itrSequence;
    coutP(ObjectHandling) << "CombinedMeasurement::ParametriseMeasurements() on scheme " << thisScheme.GetName() << endl;
    std::list< std::string > thisExpressions = thisScheme.GetExpressions();

    // Loop through the expressions
    stringstream editStr;
    editStr << "EDIT::" << originalName <<"(" << originalName;

    for (std::list< std::string >::const_iterator itrScheme = thisExpressions.begin(), end = thisExpressions.end(); itrScheme != end; ++itrScheme) {
      // Get the expression to execute
      TString thisExpression = *itrScheme;

      if (thisExpression.Contains('=')) {
        // convert to proper factory replacement
        editStr << "," << thisExpression;
      } else {
        coutP(ObjectHandling) << "CombinedMeasurement::ParametriseMeasurements() trying to call " << thisExpression << endl;
        fWorkSpace->factory(thisExpression);
      }
    }

    editStr << ")";
    coutP(ObjectHandling) << "CombinedMeasurement::ParametriseMeasurements() editing command so far " << editStr.str() << endl;

    // replacement
    fWorkSpace->factory(editStr.str().c_str());
    RooAbsPdf* tmpPdf = fWorkSpace->pdf(originalName.c_str());

    if (!tmpPdf) {
      coutE(ObjectHandling) << "CombinedMeasurement::ParametriseMeasurements() something went wrong when editing the Pdf for channel " << fName << ". Skipping this replacement." << endl;
      // continue;
    } else {
      coutI(ObjectHandling) << "CombinedMeasurement::ParametriseMeasurements() successful replacement, continuing." << endl;
      fPdf = tmpPdf;
    }

    // add new nuisance parameters and global observables
    std::list< std::string > thisNewNuisanceParameter = thisScheme.GetNewNuisanceParameters();
    for (std::list< std::string >::const_iterator itrNP = thisNewNuisanceParameter.begin(), end = thisNewNuisanceParameter.end(); itrNP != end; ++itrNP) {
      // Get the expression to execute
      TString thisExpression = *itrNP;
      fNuisanceParameters->add(*fWorkSpace->var(thisExpression));
    }

    std::list< std::string > thisNewGlobalObservable = thisScheme.GetNewGlobalObservables();
    for (std::list< std::string >::const_iterator itrGlob = thisNewGlobalObservable.begin(), end = thisNewGlobalObservable.end(); itrGlob != end; ++itrGlob) {
      // Get the expression to execute
      TString thisExpression = *itrGlob;
      fGlobalObservables->add(*fWorkSpace->var(thisExpression));
    }
    fGlobalObservables->setAttribAll("GLOBAL_OBSERVABLE");
  }

  MakeCleanWorkspace();
}

// ____________________________________________________________________________|__________
// Set parameters of interest for the combined model
void RooFitUtils::CombinedMeasurement::DefineParametersOfInterest( std::string ParametersOfInterest, ModelConfig* tmpModelConfig )
{
  TString allParametersOfInterest = ParametersOfInterest;
  TObjArray* allParametersOfInterestArray = allParametersOfInterest.Tokenize(",");
  unsigned int numPOIs = allParametersOfInterestArray->GetEntries();

  RooWorkspace* tmpWorkspace = tmpModelConfig->GetWorkspace();
  RooArgSet* oldParametersOfInterest = fParametersOfInterest;
  fParametersOfInterest->removeAll();

  for (unsigned int itrPOIs = 0; itrPOIs < numPOIs; ++itrPOIs) {
    TString thisParametersOfInterest = ((TObjString*)allParametersOfInterestArray->At(itrPOIs))->GetString();

    if (tmpWorkspace->var(thisParametersOfInterest)) {
      fParametersOfInterest->add(*tmpWorkspace->var(thisParametersOfInterest));
      coutP(ObjectHandling) << "CombinedMeasurement::DefineParametersOfInterest(" << fName << ") adding POI: " << thisParametersOfInterest << endl;

      if (fNuisanceParameters->find(*tmpWorkspace->var(thisParametersOfInterest))) {
        fNuisanceParameters->remove(*fNuisanceParameters->find(*tmpWorkspace->var(thisParametersOfInterest)));
        coutP(ObjectHandling) << "CombinedMeasurement::DefineParametersOfInterest(" << fName << ") removing new POI " << thisParametersOfInterest << " from nuisance parameters" << endl;
      }
    } else {
      coutE(ObjectHandling) << "CombinedMeasurement::DefineParametersOfInterest(" << fName << ") failed adding POI: " << thisParametersOfInterest << endl;
    }
  }

  // Move old POIs to nuisance parameter in case they still exist and are not a POI any longer
  TIterator* PoiItr = oldParametersOfInterest->createIterator();
  RooRealVar* nextParameterOfInterest;
  while ((nextParameterOfInterest = (RooRealVar*)PoiItr->Next())) {
    if (tmpWorkspace->var(nextParameterOfInterest->GetName()) && !fParametersOfInterest->find(*nextParameterOfInterest)) {
      coutP(ObjectHandling) << "CombinedMeasurement::DefineParametersOfInterest(" << fName << ") transferring old POI " << nextParameterOfInterest->GetName() << " to nuisance parameters" << endl;
      fNuisanceParameters->add(*nextParameterOfInterest);
    }
  }
  delete PoiItr;

  allParametersOfInterestArray->Delete();
  delete allParametersOfInterestArray;

}

// ____________________________________________________________________________|__________
// Clean the workspace
void RooFitUtils::CombinedMeasurement::MakeCleanWorkspace()
{
  coutP(ObjectHandling) << "CombinedMeasurement::MakeCleanWorkspace(" << fName << ") Building clean workspace" << endl;

  fParametersOfInterest->sort();
  fNuisanceParameters->sort();
  fGlobalObservables->sort();
  fObservables->sort();

  RooWorkspace* tmpWorkspace = new RooWorkspace(fWorkSpace->GetName());
  tmpWorkspace->autoImportClassCode(kTRUE);
  ModelConfig* tmpModelConfig = new ModelConfig(fModelConfig->GetName(), tmpWorkspace);

  std::string tmpPdfName = fPdf->GetName();
  tmpWorkspace->import(*fPdf, RecycleConflictNodes());
  delete fPdf;
  fPdf = tmpWorkspace->pdf(tmpPdfName.c_str());

  tmpWorkspace->import(*fData);
  if (fAsimovData) tmpWorkspace->import(*fAsimovData);

  tmpModelConfig->SetPdf(*fPdf);
  DefineParametersOfInterest(fParametersOfInterestString, tmpModelConfig);
  tmpModelConfig->SetParametersOfInterest(*fParametersOfInterest);
  tmpModelConfig->SetNuisanceParameters(*fNuisanceParameters);
  tmpModelConfig->SetObservables(*fObservables);
  tmpModelConfig->SetGlobalObservables(*fGlobalObservables);

  tmpWorkspace->import(*tmpModelConfig);

  fModelConfig = tmpModelConfig;
  fWorkSpace = tmpWorkspace;

  RooArgSet* prunedNuisanceParameters = (RooArgSet*)((fWorkSpace->allVars()).selectByName("PRUNED_NUIS_*"));
  RooArgSet* prunedGlobalObservables = (RooArgSet*)((fWorkSpace->allVars()).selectByName("PRUNED_GLOB_*"));
  SetAllConstant(*prunedNuisanceParameters);
  SetAllConstant(*prunedGlobalObservables);
	fWorkSpace->defineSet("ModelConfig_PrunedNuisParams", *prunedNuisanceParameters);
	fWorkSpace->defineSet("ModelConfig_PrunedGlobalObservables", *prunedGlobalObservables);

//  fObservables->Print();
//  fGlobalObservables->Print();
//  fNuisanceParameters->Print();
//  prunedNuisanceParameters->Print();
//  prunedGlobalObservables->Print();
}

// ____________________________________________________________________________|__________
// Interface to add Asimov data
void RooFitUtils::CombinedMeasurement::MakeAsimovData( bool Conditional, CombinedMeasurement::SnapshotName profileGenerateAt )
{
  MakeAsimovData(Conditional, profileGenerateAt, profileGenerateAt);
}

// ____________________________________________________________________________|__________
// Interface to add Asimov data
void RooFitUtils::CombinedMeasurement::MakeAsimovData( bool Conditional, CombinedMeasurement::SnapshotName profileAt, CombinedMeasurement::SnapshotName generateAt )
{
  coutP(ObjectHandling) << "CombinedMeasurement::MakeAsimovData(" << fName << ") adding " << (!Conditional?"unconditional":"conditional") << " Asimov data, profile at " << profileAt << ", generate at " << generateAt  << endl;

  // Generate the needed snapshots
  MakeSnapshots(profileAt, Conditional);
  if (generateAt == ucmles && generateAt != profileAt) {
    MakeSnapshots(generateAt, Conditional);
  }

  // Set the value of the POI for Asimov data generation
  RooRealVar* mu = (RooRealVar*)fModelConfig->GetParametersOfInterest()->first();
  double muVal = mu->getVal();

  if (generateAt == background) {
    muVal = 0;
  } else if (generateAt == nominal) {
    muVal = 1;
  } else if (generateAt == ucmles) {
    fWorkSpace->loadSnapshot("ucmles");
    muVal = mu->getVal();
    fWorkSpace->loadSnapshot("nominalNuis");
    fWorkSpace->loadSnapshot("nominalGlobs");
  } else {
    coutE(ObjectHandling) << "CombinedMeasurement::MakeAsimovData(" << fName << ") Unknown value for generation requested." << endl;
    return;
  }

  // Define the snapshot name to load at which the Asimov data will be generated
  std::string muStr;
  if (profileAt == background) {
    muStr = "_0";
  } else if (profileAt == nominal) {
    muStr = "_1";
  } else if (profileAt == ucmles) {
    muStr = "_muhat";
  } else {
    coutE(ObjectHandling) << "CombinedMeasurement::MakeAsimovData(" << fName << ") Unknown value for profiling requested." << endl;
    return;
  }

  // Collect all values for Asimov data generation
  if (!Conditional) {
    fWorkSpace->loadSnapshot("nominalGlobs");
    fWorkSpace->loadSnapshot("nominalNuis");
  }
  fWorkSpace->loadSnapshot("nominalGlobs");
  fWorkSpace->loadSnapshot(("conditionalNuis" + muStr).c_str());

  mu->setVal(muVal);

  RooArgSet* genPoiValues = (RooArgSet*)fParametersOfInterest->snapshot();
  RooArgSet*  allParams = fModelConfig->GetPdf()->getParameters(fData);
  RooStats::RemoveConstantParameters(allParams);
  *allParams = *genPoiValues;


  RooAbsData* thisAsimovData = (RooAbsData*)RooStats::AsymptoticCalculator::MakeAsimovData(*fModelConfig, *allParams, *fGlobalObservables);
  thisAsimovData->SetName(("asimovData" + (profileAt == generateAt ? muStr : "")).c_str());
  thisAsimovData->SetTitle(("asimovData" + (profileAt == generateAt ? muStr : "")).c_str());
  fWorkSpace->import(*thisAsimovData);
  fAsimovData = thisAsimovData;

  delete allParams;

  fWorkSpace->loadSnapshot("nominalNuis");
  fWorkSpace->loadSnapshot("nominalGlobs");
}

// ____________________________________________________________________________|__________
// Unfold constraints, needed for snapshot making
void RooFitUtils::CombinedMeasurement::UnfoldConstraints( RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter )
{
  if (counter > 50) {
    coutF(ObjectHandling) << "CombinedMeasurement::UnfoldConstraints(" << fName << ") failed to unfold constraints! Details given below." << endl;
    coutF(ObjectHandling) << "CombinedMeasurement::UnfoldConstraints(" << fName << ") Initial: " << endl;
    initial.Print("v");
    cout << endl;
    coutF(ObjectHandling) << "CombinedMeasurement::UnfoldConstraints(" << fName << ") Final: " << endl;
    final.Print("v");
    exit(-1);
  }

  TIterator* itr = initial.createIterator();
  RooAbsPdf* thisPdf;
  while ((thisPdf = (RooAbsPdf*)itr->Next())) {
    RooArgSet nuis_tmp = nuis;
    RooArgSet constraint_set(*thisPdf->getAllConstraints(obs, nuis_tmp, false));
    if (thisPdf->IsA() != RooGaussian::Class() && thisPdf->IsA() != RooLognormal::Class() && thisPdf->IsA() != RooGamma::Class() && thisPdf->IsA() != RooPoisson::Class() && thisPdf->IsA() != RooBifurGauss::Class()) {
      counter++;
      UnfoldConstraints(constraint_set, final, obs, nuis, counter);
    } else {
      final.add(*thisPdf);
    }
  }
  delete itr;
}

// ____________________________________________________________________________|__________
// Add snapshots to the workspace
void RooFitUtils::CombinedMeasurement::MakeSnapshots( CombinedMeasurement::SnapshotName Snapshot, bool Conditional )
{
  coutP(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") making snapshots." << endl;

  // Set the value of the POI for (conditional) profiling
  RooRealVar* mu = (RooRealVar*)fModelConfig->GetParametersOfInterest()->first();
  double muVal = mu->getVal();
  if (Snapshot == background) {
    muVal = 0;
  } else if (Snapshot == nominal) {
    muVal = 1;
  } else if (Snapshot == ucmles) {
    muVal = 1; // will be set later
  } else {
    coutE(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Unknown snapshot requested." << endl;
    return;
  }

  mu->setVal(muVal);

  // Define the snapshot name for saving
  std::string muStr;
  if (Snapshot == background) {
    muStr = "_0";
  } else if (Snapshot == nominal) {
    muStr = "_1";
  } else {
    muStr = "_muhat";
  }

  RooArgSet mc_obs = *fModelConfig->GetObservables();
  RooArgSet mc_globs = *fModelConfig->GetGlobalObservables();
  RooArgSet mc_nuis = *fModelConfig->GetNuisanceParameters();

  // Pair the nuisance parameters and global observables
  RooArgSet mc_nuis_tmp = mc_nuis;
  RooArgList nui_list("ordered_nuis");
  RooArgList glob_list("ordered_globs");
  RooArgSet constraint_set_tmp(*fPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
  RooArgSet constraint_set;

  int counter_tmp = 0;
  UnfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);

  TIterator* cIter = constraint_set.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)cIter->Next())) {
    RooAbsPdf* pdf = (RooAbsPdf*)arg;
    if (!pdf) continue;

    TIterator* nIter = mc_nuis.createIterator();
    RooRealVar* thisNui = NULL;
    RooAbsArg* nui_arg;
    while ((nui_arg = (RooAbsArg*)nIter->Next())) {
      if (pdf->dependsOn(*nui_arg)) {
        thisNui = (RooRealVar*)nui_arg;
        break;
      }
    }
    delete nIter;

    // Need this incase the observable isn't fundamental.
    // In this case, see which variable is dependent on the nuisance parameter and use that.
    RooArgSet* components = pdf->getComponents();
    components->remove(*pdf);
    if (components->getSize()) {
      TIterator* itr1 = components->createIterator();
      RooAbsArg* arg1;
      while ((arg1 = (RooAbsArg*)itr1->Next())) {
        TIterator* itr2 = components->createIterator();
        RooAbsArg* arg2;
        while ((arg2 = (RooAbsArg*)itr2->Next())) {
          if (arg1 == arg2) continue;
          if (arg2->dependsOn(*arg1)) {
            components->remove(*arg1);
          }
        }
        delete itr2;
      }
      delete itr1;
    }
    if (components->getSize() > 1) {
      coutE(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Failed to isolate proper nuisance parameter." << endl;
      return;
    }
    else if (components->getSize() == 1) {
      thisNui = (RooRealVar*)components->first();
    }

    TIterator* gIter = mc_globs.createIterator();
    RooRealVar* thisGlob = NULL;
    RooAbsArg* glob_arg;
    while ((glob_arg = (RooAbsArg*)gIter->Next())) {
      if (pdf->dependsOn(*glob_arg)) {
        thisGlob = (RooRealVar*)glob_arg;
        break;
      }
    }
    delete gIter;

    if (!thisNui || !thisGlob) {
      coutW(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Failed to find nuisance parameter or global observable for constraint " << pdf->GetName() << endl;
      continue;
    }

    coutP(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Pairing nuisance parameter: " << thisNui->GetName() << " with global observable: " << thisGlob->GetName() << " from constraint " << pdf->GetName() << endl;

    nui_list.add(*thisNui);
    glob_list.add(*thisGlob);
  }
  delete cIter;

  // Save the snapshots of nominal parameters
  coutP(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Saving nominal snapshots." << endl;
  fWorkSpace->saveSnapshot("nominalGlobs", *fModelConfig->GetGlobalObservables());
  fWorkSpace->saveSnapshot("nominalNuis", *fModelConfig->GetNuisanceParameters());

  RooArgSet nuiSet_tmp(nui_list);

  mu->setVal(muVal);
  if (Snapshot == background || Snapshot == nominal) {
    mu->setConstant(1);
  } else {
    mu->setConstant(0);
  }

  if (Conditional) {
    coutP(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") performing conditional fit." << endl;

    ExtendedMinimizer minimizer("minimizer", fPdf, fData);
    minimizer.minimize(Minimizer(fMinimizerType.c_str(), fMinimizerAlgo.c_str()), Strategy(fDefaultStrategy),
                       Constrain(*fModelConfig->GetNuisanceParameters()), GlobalObservables(*fModelConfig->GetGlobalObservables()),
                       NumCPU(fNumCPU, 3), Offset(1), Optimize(2));
  }

  mu->setConstant(0);

  // loop over the nui/glob list, grab the corresponding variable from the tmp ws, and set the glob to the value of the nui
  int nrNuis = nui_list.getSize();
  if (nrNuis != glob_list.getSize()) {
    coutE(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") number of nuisance paraeters does not match number of global observables" << endl;
    return;
  }

  for (int i=0;i<nrNuis;i++) {
    RooRealVar* nui = (RooRealVar*)nui_list.at(i);
    RooRealVar* glob = (RooRealVar*)glob_list.at(i);

    coutI(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Nuisance parameter: " << nui << ", global observable: " << glob << endl;
    coutI(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Setting global observable " << glob->GetName() << " (previous value: " << glob->getVal() << ") to conditional value: " << nui->getVal() << endl;

    glob->setVal(nui->getVal());
  }

  // Save the snapshots of conditional parameters
  coutP(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Saving conditional snapshots." << endl;
  fWorkSpace->saveSnapshot(("conditionalGlobs" + muStr).c_str(), *fModelConfig->GetGlobalObservables());
  fWorkSpace->saveSnapshot(("conditionalNuis" + muStr).c_str(), *fModelConfig->GetNuisanceParameters());

  if (Snapshot == ucmles) {
    RooArgSet nuisAndPOI(*fModelConfig->GetNuisanceParameters(), *mu);
    fWorkSpace->saveSnapshot("ucmles", nuisAndPOI);
  }

  // Bring us back to nominal for exporting
  coutP(ObjectHandling) << "CombinedMeasurement::MakeSnapshots(" << fName << ") Return to nominal parameter values." << endl;
  fWorkSpace->loadSnapshot("nominalNuis");
  fWorkSpace->loadSnapshot("nominalGlobs");
}

// ____________________________________________________________________________|__________
// Add multiple snapshots to the workspace
void RooFitUtils::CombinedMeasurement::MakeSnapshots( std::set<CombinedMeasurement::SnapshotName> Snapshots, bool Conditional )
{
  for (std::set<CombinedMeasurement::SnapshotName>::iterator itr = Snapshots.begin(); itr != Snapshots.end(); ++itr) {
    MakeSnapshots(*itr, Conditional);
  }
}

// ____________________________________________________________________________|__________
// Print information of the combined measurement
void RooFitUtils::CombinedMeasurement::Print()
{
  coutI(ObjectHandling) << "CombinedMeasurement::Print(" << fName << ") printing information about the comined measurement\n" << endl;

  cout << "== COMBINED SUMMARY ==== " << endl;
  cout << "Combined measurement   : " << fName << endl;
  cout << "\nIndividual measurements: " << fMeasurements.size() << endl;
  cout << "\nParameters of interest : " << fParametersOfInterest->getSize() << endl;
  PrintCollection(fParametersOfInterest);
  cout << "\nObservables            : " << fObservables->getSize() << endl;
  PrintCollection(fObservables);
  cout << "\nNuisance parameters    : " << fNuisanceParameters->getSize() << endl;
  PrintCollection(fNuisanceParameters);
  cout << "\nGlobal observables     : " << fGlobalObservables->getSize() << endl;
  PrintCollection(fGlobalObservables);
  cout << "\n" << endl;

  std::set<std::string> inputMeasurements;
  for (std::map< std::string, Measurement* >::iterator measItr = fMeasurements.begin(); measItr != fMeasurements.end(); ++measItr) {
    RooFitUtils::Measurement* meas = measItr->second;
    inputMeasurements.insert(meas->GetName());
  }

  fCorrelationScheme->Print(inputMeasurements);
}

// ____________________________________________________________________________|__________
// Print collection, map name to title (to be used as description)
void RooFitUtils::CombinedMeasurement::PrintCollection(RooAbsCollection* collection)
{
  int nrColumns = 1;
  string* header = new string[nrColumns];
  header[0] = "Description";

  int nrEntries = collection->getSize();
  if (nrEntries == 0) {
    // coutW(ObjectHandling) << "CombinedMeasurement::PrintCollection(" << this->GetName() << ") collection " << collection->GetName() << " is empty" << endl;
    cout << "CombinedMeasurement::PrintCollection() collection " << collection->GetName() << " is empty" << endl;
    delete[] header;
    return;
  }
  string* firstCol = new string[nrEntries+1];
  firstCol[0] = "Parameter";
  std::string** matrix = new string*[nrEntries];
  for (int in = 0; in < nrEntries; in++) {
    matrix[in] = new std::string[nrColumns];
    for (int i = 0; i < nrColumns; i++) {
      matrix[in][i] = "";
    }
  }

  int irow = 1;
  for (RooLinkedListIter it = collection->iterator(); RooAbsArg* arg = dynamic_cast<RooAbsArg*>(it.Next());) {
    firstCol[irow] = arg->GetName();
    int icol = 0;
    string description = arg->GetTitle();
    matrix[irow-1][icol] = description;
    irow++;
  }

  RooFitUtils::CorrelationScheme::PrintTable(firstCol, matrix, NULL, header, nrEntries, nrColumns, 0, cout, "    ");

  delete[] header;
  delete[] firstCol;

  for (int in = 0; in < nrEntries; in++) {
    delete[] matrix[in];
  }

  delete[] matrix;
}
