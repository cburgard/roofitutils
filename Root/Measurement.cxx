#include "RooFitUtils/Measurement.h"
#include "RooFitUtils/CombinedMeasurement.h"

#include "TRegexp.h"

#include "RooDataHist.h"
#include <stdexcept>

// ____________________________________________________________________________|__________
// Constructor
Measurement::Measurement( const std::string& MeasurementName, const std::string& FileName, const std::string& WorkspaceName, const std::string& ModelConfigName, const std::string& DataName, const std::string& SnapshotName, bool BinnedLikelihood )
  :
  AbsMeasurement( MeasurementName, FileName, WorkspaceName, ModelConfigName, DataName ),
  fSnapshotName ( SnapshotName ),
  fBinnedLikelihood( BinnedLikelihood )
{
  coutP(InputArguments) << "Measurement::Measurement(" << fName <<") created" << endl;
}

// ____________________________________________________________________________|__________
// Constructor
Measurement::Measurement( const std::string& MeasurementName, RooWorkspace* ws, const std::string& ModelConfigName, const std::string& DataName, const std::string& SnapshotName, bool BinnedLikelihood )
  :
  AbsMeasurement( MeasurementName, ws, ModelConfigName, DataName ),
  fSnapshotName ( SnapshotName ),
  fBinnedLikelihood( BinnedLikelihood )
{
  coutP(InputArguments) << "Measurement::Measurement(" << fName <<") created" << endl;
}

// ____________________________________________________________________________|__________
// Destructor
Measurement::~Measurement()
{
  coutI(InputArguments) << "Measurement::~Measurement(" << fName << ") cleaning up" << endl;
  // TODO

  if (fIsInitialised) {
    fNuisanceParameters->Delete();
    delete fNuisanceParameters;

    fObservables->Delete();
    delete fObservables;

    fGlobalObservables->Delete();
    delete fGlobalObservables;

    fFile->Close();
  }
}

// ____________________________________________________________________________|__________
// Import a model for a specific channel
void Measurement::initialise()
{
  if(!fWorkSpace){
    // Initialisation
    fFile = TFile::Open(fFileName.c_str());
    if (!fFile) {
      coutF(InputArguments) << "Measurement::initialise(" << fName << ") could not open file " << fFileName << endl;
      exit(-1);
    } else {
      coutP(InputArguments) << "Measurement::initialise(" << fName << ") opened file " << fFileName << endl;
    }
    
    fWorkSpace = (RooWorkspace*) fFile->Get(fWorkspaceName.c_str());
  }

  if (!fWorkSpace) {
    coutF(InputArguments) << "Measurement::initialise(" << fName << ") could not find workspace " << fWorkspaceName << endl;
    exit(-1);
  } else {
    coutP(InputArguments) << "Measurement::initialise(" << fName << ") using workspace " << fWorkspaceName << endl;
  }

  if (fBinnedLikelihood) {
    // Activate binned likelihood calculation for binned models
    RooFIter iter = fWorkSpace->components().fwdIterator() ;
    RooAbsArg* arg ;
    while((arg = iter.next())) {
      if (arg->IsA() == RooRealSumPdf::Class()) {
        arg->setAttribute("BinnedLikelihood");
        coutP(InputArguments) << "Measurement::initialise(" << fName << ") activating BinnedLikelihood for " << arg->GetName() << endl;
      }
    }
  }

  fModelConfig = (ModelConfig*)fWorkSpace->obj(fModelConfigName.c_str());
  if (!fModelConfig) {
    coutE(InputArguments) << "Measurement::initialise(" << fName << ") could not find ModelConfig " << fModelConfigName << endl;
    exit(-1);
  } else {
    coutP(InputArguments) << "Measurement::initialise(" << fName << ") using ModelConfig " << fModelConfigName << endl;
  }

  fData = (RooAbsData*) fWorkSpace->obj(fDataName.c_str());
  if (!fData) {
    coutF(InputArguments) << "Measurement::initialise(" << fName << ") could not find data " << fDataName << endl;
    throw std::runtime_error(TString::Format("unable to find data with name '%s'",fDataName.c_str()).Data());
  } else {
    coutP(InputArguments) << "Measurement::initialise(" << fName << ") using data " << fDataName << endl;
  }

  // fix POIs at their nominal values, set nuisance parameters to their nominal values
  fWorkSpace->loadSnapshot(fSnapshotName.c_str());
  const RooArgSet* pois = fModelConfig->GetParametersOfInterest();
  if(!pois){
    coutE(InputArguments) << "Measurement::initialise(" << fName << ") could not find list of POIs " << endl;
  } else {
    for (RooLinkedListIter it = pois->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      // v->setVal(1.0);
      v->setConstant(1);
    }
  }

  // Grab all needed information: categories, (global) observables after editing the workspace
  // find top-level RooSimultaneous, which is the physics model
  // fPdf = (RooAbsPdf*)fModelConfig->GetPdf();
  fPdf = fModelConfig->GetPdf();
  RooSimultaneous* thisSim = NULL;
  if(!fPdf){
    coutF(InputArguments) << "Measurement::initialise(" << fName << ") Did not find a Pdf in ModelConfig." << endl;
    throw std::runtime_error("unable to find Pdf in ModelConfig");
  } else if (fPdf->IsA() == RooSimultaneous::Class()) {
   thisSim = (RooSimultaneous*)fPdf;
  } else if (fPdf->IsA() == RooProdPdf::Class()) {
    coutI(InputArguments) << "Measurement::initialise(" << fName << ") Found a RooProdPdf, expected a RooSimultaneous." << endl;
    coutI(InputArguments) << "Measurement::initialise(" << fName << ") Trying to find top-level RooSimultaneous." << endl;

    RooArgList pdfList = ((RooProdPdf*)fPdf)->pdfList();
    TIterator* pdfItr = pdfList.createIterator();
    RooAbsArg* nextArg;
    while ((nextArg = (RooAbsArg*)pdfItr->Next())) {
      if (nextArg->IsA() == RooSimultaneous::Class()) {
        coutI(InputArguments) << "Measurement::initialise(" << fName << ") Found RooSimultaneous " << nextArg->GetName() << endl;
        thisSim = (RooSimultaneous*)nextArg;
        thisSim->SetName(fPdf->GetName());
        thisSim->SetTitle(fPdf->GetTitle());
        break;
      }
    }
    delete pdfItr;
  } else {
    coutF(InputArguments) << "Measurement::initialise(" << fName << ") Did not find a RooSimultaneous." << endl;
    throw std::runtime_error("unable to find RooSimultaneous Pdf");
  }

  // RooCategory* cat = (RooCategory*)&((RooSimultaneous*)fPdf)->indexCat();

  // thisSim->Print();
  // fData->Print();

  // Maybe bin the data
  if (fSetBinning) {
    fData = SetDatasetBinning(thisSim, fData, fSetNbins, fGenerateBinnedTag, fBinnedCategories, fUnbinnedCategories, fWeightVarName);
  }

  RooCategory* cat = (RooCategory*)&thisSim->indexCat();
  dataList = fData->split(*cat, true);
  fObservables = (RooArgSet*)fModelConfig->GetObservables();
  fGlobalObservables = (RooArgSet*)fModelConfig->GetGlobalObservables();
  fNuisanceParameters = (RooArgSet*)fModelConfig->GetNuisanceParameters();

  // Fix for non-binned data
  if (fData->IsA() == RooDataHist::Class()) {
    coutW(InputArguments) << "Measurement::initialise(" << fName << ") data " << fDataName << " is a RooDataHist. Transform it to a RooDataSet." << endl;
    hist2dataset(thisSim, cat);
  }

  // Find all constraint terms for later usage in regularisation
  // Copies, to keep original sets in place after getAllconstraints call
  RooArgSet tmpAllNuisanceParameters2 = *fNuisanceParameters;
  RooArgSet tmpAllObservables2 = *fObservables;
  RooArgSet* AllConstraints = fPdf->getAllConstraints(tmpAllObservables2, tmpAllNuisanceParameters2, kFALSE);

  // Take care of the case where we have a product of constraint terms
  TIterator* ConstraintItr = AllConstraints->createIterator();
  RooAbsArg* nextConstraint;
  fConstraints = new RooArgSet(AllConstraints->GetName());
  while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next())) {
    if (nextConstraint->IsA() == RooProdPdf::Class()) {
      RooArgSet thisComponents;
      FindUniqueProdComponents((RooProdPdf*)nextConstraint, thisComponents);
      fConstraints->add(thisComponents);
    } else {
      fConstraints->add(*nextConstraint);
    }
  }
  delete ConstraintItr;

  // Make the fPdf the physics pdf of the model
  fPdf = (RooSimultaneous*)thisSim;

  // Prune nuisance parameters if requested
  if (fIsPrunable) {
    PruneNuisanceParameters();
  }

  // All initialised
  fIsInitialised = kTRUE;
}

// ____________________________________________________________________________|__________
// Collect all measurement after regularising their individual channels
void Measurement::CollectChannels()
{
  if (!fIsInitialised) {
    initialise();
  }

  RooSimultaneous* thisSim = (RooSimultaneous*)fPdf;
  int nrEntries = dataList->GetEntries();

  RooArgSet* NuisanceParametersPlusPOIs = (RooArgSet*)fModelConfig->GetNuisanceParameters();
  const RooArgSet* pois = fModelConfig->GetParametersOfInterest();
  if(pois){
    NuisanceParametersPlusPOIs->add(*pois);
  }

  RooArgSet* tmpNuisanceParameters = new RooArgSet("NuisanceParameters");
  RooArgSet* tmpObservables = new RooArgSet("Observables");
  RooArgSet* tmpGlobalObservables = new RooArgSet("GlobalObservables");

	// Add missing attributes to the renaming map
	fRenamingMap.AddAttributes(fModelConfig);

  // Tag parameters that should get pruned
  if (fIsPrunable){
    TagPrunableParameters();
  }

  // TODO: the union implementation only separates standalone expressions -> switch c++11 regex when using cling
  TObjArray* thisFilterArray = (TString(fChannelFilter).ReplaceAll(" ", "")).Tokenize("|");
  unsigned int numRules = thisFilterArray->GetEntries();

  for (int itrChan = 0; itrChan < nrEntries; ++itrChan) {
    RooAbsData* ds = (RooAbsData*)dataList->At(itrChan);
    std::string label(ds->GetName());
    std::string thisChannelName = label + "_" + fName.Data();

    // Filter channels based on a regular expression
    bool useChannel = false;

    for (unsigned int itrRule = 0; itrRule < numRules; ++itrRule) {
      TString thisRule = ((TObjString*)thisFilterArray->At(itrRule))->GetString();
      TRegexp reg(thisRule);
      Ssiz_t dummy(0);
      if (reg.Index(TString(thisChannelName), &dummy, 0) != -1) {
        useChannel = true;
        break;
      }
    }

    if (!useChannel) {
      coutW(ObjectHandling) << "Measurement::CollectChannels(" << fName << ") not using channel" << thisChannelName << endl;
      fExcludedChannels.push_back(thisChannelName);
      continue;
    }

    // Weight the dataset in case it is not already
    if (!ds->isWeighted()) {
      coutW(ObjectHandling) << "Measurement::CollectChannels(" << fName << ") dataset " << ds->GetName() << "is not weighted. Adding weight 1." << endl;
      RooRealVar* weightVar = new RooRealVar( "weightVar", "", 1., -1e10, 1e10 );
      ((RooDataSet*)ds)->addColumn(*weightVar);
      RooArgSet* obs_cat_weight = new RooArgSet();
      obs_cat_weight->add(*fObservables);
      obs_cat_weight->add(*weightVar);
      RooDataSet* tmpDs = new RooDataSet(ds->GetName(), ds->GetTitle(), *obs_cat_weight, Import(*(RooDataSet*)ds), WeightVar("weightVar"));
      ds = tmpDs;
    }

    // Get the pdf
    RooAbsPdf* thisPdf = (RooAbsPdf*)thisSim->getPdf(label.c_str());

    // Attach all (necessary) constraints to the current channel to make it stand-alone
    RooArgSet thistmpNuisanceParameters;
    TIterator* NuisItr = fNuisanceParameters->createIterator();
    RooRealVar* nextNuisanceParameter;
    while ((nextNuisanceParameter = (RooRealVar*)NuisItr->Next())) {
      if (thisPdf->dependsOn(*nextNuisanceParameter)) {
				thistmpNuisanceParameters.add(*nextNuisanceParameter);
			}
    }

    RooArgSet* tmpComp = new RooArgSet();
    tmpComp->add(*thisPdf);
    TIterator* ConstraintItr = fConstraints->createIterator();
    RooAbsArg* nextConstraint;
    while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next())) {
      delete NuisItr;
      NuisItr = thistmpNuisanceParameters.createIterator();
      while ((nextNuisanceParameter = (RooRealVar*)NuisItr->Next())) {
        if (nextConstraint->dependsOn(*nextNuisanceParameter)) {
          tmpComp->add(*nextConstraint);
          coutI(ObjectHandling) << "Measurement::CollectChannels(" << fName << ") using " << nextConstraint->GetName() << " for channel " <<  thisChannelName << endl;
          break;
        }
      }
    }

    // Build a new RooProdPdf from channel and all constraints, unifying the structure of all pdfs:
    // physics model x subsidiary measurements
    RooProdPdf* prodPdf = new RooProdPdf(thisPdf->GetName(), thisPdf->GetTitle(), *tmpComp);

    // Find the unique, fundamental components
    RooArgSet thisComponents;
    FindUniqueProdComponents(prodPdf, thisComponents);
    //    thisComponents.Print();

    // Build a new, clean RooProdPdf used for this channel
    RooProdPdf* thisProdPdf = new RooProdPdf(thisPdf->GetName(),thisPdf->GetName(), thisComponents);

    string newName = string(thisPdf->GetName()) + "_complete";
    thisProdPdf->SetName(newName.c_str());
    thisProdPdf->SetTitle(newName.c_str());

    // Make the channel and regularise it.
    Channel* thisChannel = new Channel(thisChannelName, (RooAbsPdf*)thisProdPdf, ds, string(fName));
    thisChannel->SetRenamingMap(fRenamingMap);
    thisChannel->SetCorrelationFactors(fCorrelationFactors);
    thisChannel->SetGlobalObservables(fGlobalObservables);
    thisChannel->SetNuisanceParameters(NuisanceParametersPlusPOIs);
    thisChannel->SetObservables(fObservables);
    thisChannel->RegulariseChannel();

    // Retrieve observables, etc.
    RooArgSet* thisObservables = thisChannel->GetObservables();
    RooArgSet* thisGlobalObservables = thisChannel->GetGlobalObservables();
    RooArgSet* thisNuisanceParameters = thisChannel->GetNuisanceParameters();

    for (RooLinkedListIter it = thisObservables->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      tmpObservables->add(*v);
    }

    for (RooLinkedListIter it = thisGlobalObservables->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      tmpGlobalObservables->add(*v);
    }
    fGlobalObservables->setAttribAll("GLOBAL_OBSERVABLE");

    for (RooLinkedListIter it = thisNuisanceParameters->iterator(); RooRealVar* v = dynamic_cast<RooRealVar*>(it.Next());) {
      tmpNuisanceParameters->add(*v);
    }

    fChannels.push_back(thisChannel);

    delete tmpComp;
    delete prodPdf;
    // delete thisProdPdf;
    delete NuisItr;
    delete ConstraintItr;
  }

  SetObservables(tmpObservables);
  SetGlobalObservables(tmpGlobalObservables);
  SetNuisanceParameters(tmpNuisanceParameters);

  delete thisSim;
}

// ____________________________________________________________________________|__________
// Find all unique components of a RooProdPdf
void Measurement::FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components )
{
  static int counter = 0;
  counter++;

  if (counter > 50) {
    coutE(ObjectHandling) << "Measurement::FindUniqueProdComponents(" << fName << ") detected infinite loop. Please check." << endl;
    exit(1);
  }

  RooArgList pdfList = Pdf->pdfList();
  if (pdfList.getSize() == 1) {
    coutI(ObjectHandling) << "Measurement::FindUniqueProdComponents(" << fName << ") " << pdfList.at(0)->GetName() << " is fundamental." << endl;
    Components.add(pdfList);
  } else {
    TIterator* pdfItr = pdfList.createIterator();
    RooAbsArg* nextArg;
    while ((nextArg = (RooAbsArg*)pdfItr->Next())) {
      RooProdPdf* Pdf = (RooProdPdf*)nextArg;
      if (string(Pdf->ClassName()) != "RooProdPdf") {
        coutI(ObjectHandling) << "Measurement::FindUniqueProdComponents(" << fName << ") " << Pdf->GetName() << " is no RooProdPdf. Adding it." << endl;
        Components.add(*Pdf);
        continue;
      }
      FindUniqueProdComponents(Pdf, Components);
    }
    delete pdfItr;
  }
  counter = 0;
}

// ____________________________________________________________________________|__________
//
RenamingMap::ConstraintType Measurement::DetermineConstraintType( RooRealVar* Parameter )
{
  std::string tmpConstraintType = "unconstrained";
  TIterator* ConstraintItr = fConstraints->createIterator();
  RooAbsArg* nextConstraint;
  bool foundConstraint = kFALSE;
  while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next()) && !foundConstraint) {
    if (nextConstraint->dependsOn(*Parameter)) {
      foundConstraint = kTRUE;

      // find type of constraint
      TString thisConstraintType = nextConstraint->ClassName();
      tmpConstraintType = (thisConstraintType.ReplaceAll("Roo", "")).Data();
    }
  }
  delete ConstraintItr;

  if (tmpConstraintType == "Gaussian") return RenamingMap::Gaussian;
  if (tmpConstraintType == "Poisson") return RenamingMap::Poisson;
  else return RenamingMap::unknown;
}

// ____________________________________________________________________________|__________
// Print information about the measurement
void Measurement::Print()
{
  cout << "== MEASUREMENT SUMMARY ==== " << endl;
  cout << "Measurement           : " << fName << endl;
  cout << "\nChannels              : " << fChannels.size() << endl;
  cout << "\nExcluded channels     : " << fExcludedChannels.size() << endl;
  if(fModelConfig && fModelConfig->GetParametersOfInterest()){
    cout << "\nParameters of interest: " << fModelConfig->GetParametersOfInterest()->getSize() << endl;
    CombinedMeasurement::PrintCollection((RooArgSet*)fModelConfig->GetParametersOfInterest());
  } else {
    cout << "\nParameters of interest: NONE" << endl;
  }
  if(fObservables){
    cout << "\nObservables           : " << fObservables->getSize() << endl;
    CombinedMeasurement::PrintCollection(fObservables);
  } else {
    cout << "\nObservables           : NONE" << endl;
  }
  if(fNuisanceParameters){
    cout << "\nNuisance parameters   : " << fNuisanceParameters->getSize() << endl;
    CombinedMeasurement::PrintCollection(fNuisanceParameters);
  } else {
    cout << "\nNuisance parameters   : NONE" << endl;
  }
  if(fGlobalObservables){
    cout << "\nGlobal observables    : " << fGlobalObservables->getSize() << endl;
    CombinedMeasurement::PrintCollection(fGlobalObservables);
  } else {
    cout << "\nGlobal observables    : NONE" << endl;
  }
  cout << "\nRenaming map          : " << fName << endl;
  fRenamingMap.Print();
  cout << "\n" << endl;
}

// ____________________________________________________________________________|__________
void Measurement::TagPrunableParameters()
{
  std::map< std::string, std::string > thisRenamingMap = fRenamingMap.GetRenamingMap();

  for (std::list< std::string >::const_iterator pruneItr = fPrunedNuisanceParameters.begin(), end = fPrunedNuisanceParameters.end(); pruneItr != end; ++pruneItr) {
    std::string thisParameter = *pruneItr;
    std::string thisNewObservableName = thisRenamingMap[thisParameter];

    std::string thisNewConstraintName        = fRenamingMap.GetAttribute(thisNewObservableName, RenamingMap::Constraint, RenamingMap::combined);
    std::string thisNewObservableRange       = fRenamingMap.GetAttribute(thisNewObservableName, RenamingMap::ObservableRange, RenamingMap::combined);
    std::string thisNewGlobalObservableName  = fRenamingMap.GetAttribute(thisNewObservableName, RenamingMap::GlobalObservable, RenamingMap::combined);
    std::string thisNewGlobalObservableRange = fRenamingMap.GetAttribute(thisNewObservableName, RenamingMap::GlobalObservableRange, RenamingMap::combined);
    std::string thisNewConstraintType        = fRenamingMap.GetAttribute(thisNewObservableName, RenamingMap::Type, RenamingMap::combined);
    std::string thisNewSigmaName             = fRenamingMap.GetAttribute(thisNewObservableName, RenamingMap::Sigma, RenamingMap::combined);
    std::string thisNewSigmaRange            = fRenamingMap.GetAttribute(thisNewObservableName, RenamingMap::SigmaRange, RenamingMap::combined);

    thisNewObservableName = "PRUNED_NUIS_" + thisNewObservableName;
    if (thisNewConstraintName != "") thisNewConstraintName = "PRUNED_" + thisNewConstraintName;
    if (thisNewGlobalObservableName != "") thisNewGlobalObservableName = "PRUNED_GLOB_" + thisNewGlobalObservableName;
    if (thisNewSigmaName != "") thisNewSigmaName = "PRUNED_" + thisNewSigmaName;

    thisRenamingMap[thisParameter] = thisNewObservableName;
    fRenamingMap.SetRenamingMap(thisRenamingMap);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::Observable, thisNewObservableName, RenamingMap::combined);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::Constraint, thisNewConstraintName, RenamingMap::combined);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::ObservableRange, thisNewObservableRange, RenamingMap::combined);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::GlobalObservable, thisNewGlobalObservableName, RenamingMap::combined);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::GlobalObservableRange, thisNewGlobalObservableRange, RenamingMap::combined);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::Type, thisNewConstraintType, RenamingMap::combined);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::Sigma, thisNewSigmaName, RenamingMap::combined);
    fRenamingMap.SetAttribute(thisNewObservableName, RenamingMap::SigmaRange, thisNewSigmaRange, RenamingMap::combined);
  }
}

// ____________________________________________________________________________|__________
// Convert a RooDataHist to a RooDataSet incl. weights.
// Courtesy of Hongtao Yang <Hongtao.Yang@cern.ch>.
void Measurement::hist2dataset( RooSimultaneous* thisSim, RooCategory* cat )
{
    RooRealVar* x[500], *w[500];
    RooDataSet* data[500];
    std::map<std::string, RooDataSet*> datasetMap;
    RooArgSet *Observables = new RooArgSet();
    int numChannels = dataList->GetEntries();

    for (int ich = 0; ich < numChannels; ich++) {
      cat->setBin(ich);
      TString channelname=cat->getLabel();
      RooAbsPdf* pdfi = thisSim->getPdf(channelname.Data());
      RooAbsData* datai = (RooAbsData*)(dataList->At(ich));
      RooRealVar *obsi = (RooRealVar*)pdfi->getObservables(datai)->first();
      x[ich] = fWorkSpace->var(obsi->GetName());
      w[ich] = new RooRealVar(Form("wt_%d", ich), Form("wt_%d", ich), 1);

      RooArgSet* args = new RooArgSet();
      args->add(RooArgSet(*x[ich], *w[ich]));
      data[ich] = new RooDataSet("combData", "combData", *args, WeightVar(*w[ich]));

      RooArgSet* obs_tmp = (RooArgSet*)datai->get();
      RooRealVar* xdata_tmp = (RooRealVar*)obs_tmp->find(obsi->GetName());

      for (int ievt = 0; ievt < datai->numEntries(); ievt++) {
        datai->get(ievt) ;
        double weight = datai->weight();
        x[ich]->setVal(xdata_tmp->getVal());
        w[ich]->setVal(weight);
        data[ich]->add(RooArgSet(*x[ich], *w[ich]), weight);
      }
      Observables->add(*x[ich]);
      datasetMap[channelname.Data()] = data[ich];
    }

    RooRealVar wt("wt", "wt", 1); // , 0, 10000);
    RooArgSet *args = new RooArgSet();
    args->add(*Observables);
    args->add(wt);
    RooDataSet* combData = new RooDataSet("combData", "combData", *args, Index(*cat), Import(datasetMap), WeightVar(wt));
    fWorkSpace->import(*combData);

    fData=fWorkSpace->data(combData->GetName());
    dataList = fData->split(*cat, true);
    delete combData;
}