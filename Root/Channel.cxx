#include "RooFitUtils/Channel.h"

// ____________________________________________________________________________|__________
// Constructor
Channel::Channel( std::string ChannelName, RooAbsPdf* Pdf, RooAbsData* Data, std::string ParentName )
  :
  TNamed( ChannelName.c_str(), ChannelName.c_str() ),
  fPdf( Pdf ),
  fData( Data ),
  fParentMeasurement( ParentName )
{
  coutP(InputArguments) << "Channel::Channel(" << fName <<") new Channel from " << fParentMeasurement << " with PDF " << fPdf->GetName() << " and data " << fData->GetName() << "." << endl;

  // Create a workspace to use the workspace factory for the regularisation
  std::string tmpWSName = fParentMeasurement + "_" + fName.Data();
  coutP(ObjectHandling) << "Channel::Channel(" << fName << ") making workspace " << tmpWSName << " for using the workspace factory." << endl;
  fWorkSpace = new RooWorkspace(tmpWSName.c_str());
}

// ____________________________________________________________________________|__________
// Destructor
Channel::~Channel()
{
  coutI(InputArguments) << "Channel::~Channel(" << fName << ") cleaning up" << endl;
  // TODO
}

// ____________________________________________________________________________|__________
// Edit the pdf and parameters according to renaming scheme
void Channel::RegulariseChannel()
{
  std::map< std::string, std::string > thisRenamingMap = fRenamingMap.GetRenamingMap();
  std::map< std::string, std::pair< TString, TMatrixDSym > > thisCorrelationFactors = fCorrelationFactors;

  // Import Pdf into the workspace, append unique suffix (i.e. measurement name to which
  // this channel belongs) and retrieve the Pdf from the workspace again.
  coutI(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") importing Pdf " << fPdf->GetName() << " in workspace " << fWorkSpace->GetName() << endl;
  fWorkSpace->import(*fPdf, RenameAllNodes(fParentMeasurement.c_str()), RenameAllVariables(fParentMeasurement.c_str()), Silence());

  // Populate workspace with the new nuisance parameters
  // and fill tmp maps with just the nuisance parameter and global observables names for later usage
  std::map< std::string, std::string > thisRenamedGlobs;
  std::map< std::string, std::string > thisRenamedNuis;
  // std::map< std::string, std::string > thisRenamedConstraints;

  for (std::map< std::string, std::string >::iterator it = thisRenamingMap.begin(); it != thisRenamingMap.end(); ++it) {
    std::string oldObsName = it->first;
    std::string newObsName = it->second;

    std::string renamedVar = oldObsName + "_" + fParentMeasurement;

    // fill local maps only if a true replacement planned and channel depends on old variable
    if (fWorkSpace->obj(renamedVar.c_str())) {
      thisRenamedNuis[oldObsName] = newObsName;

      std::string oldGlobsName = fRenamingMap.GetAttribute(oldObsName, RenamingMap::GlobalObservable, RenamingMap::individual);
      std::string newGlobsName = fRenamingMap.GetAttribute(newObsName, RenamingMap::GlobalObservable, RenamingMap::combined);

      thisRenamedGlobs[oldGlobsName] = newGlobsName;

      std::string factoryString = fRenamingMap.FactoryExpression(newObsName, RenamingMap::combined);

      fWorkSpace->factory(factoryString.c_str());
    }
  }

  // Add partial correlations to the pdf
  std::string tmpAllProducts = string(fPdf->GetName()) + "_" + fParentMeasurement;
  std::string newPdfName = string(fPdf->GetName()) + "_" + fParentMeasurement + "Corr";
  for (std::map< std::string, std::pair< TString, TMatrixDSym > >::iterator it = thisCorrelationFactors.begin(); it != thisCorrelationFactors.end(); ++it) {
    TString thisCommonName = it->first;
    cout << thisCommonName << endl;
    std::pair< TString, TMatrixDSym > thisFactor = it->second;
    fWorkSpace->import(thisFactor.second, thisCommonName + "_corr");
    fWorkSpace->factory(thisFactor.first);
    cout << thisFactor.first << endl;
    tmpAllProducts += "," + thisCommonName + "Corr";
  }

  // Retrieve the pdf
  if (!thisCorrelationFactors.empty()) {
    cout << "PROD::"+newPdfName+"({"+tmpAllProducts+"})" << endl;
    fWorkSpace->factory(("PROD::"+newPdfName+"({"+tmpAllProducts+"})").c_str());
  } else {
    newPdfName = string(fPdf->GetName()) + "_" + fParentMeasurement;
  }
  coutI(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") getting Pdf " << newPdfName << " from workspace" << fWorkSpace->GetName() << endl;
  fPdf = (RooAbsPdf*)fWorkSpace->pdf(newPdfName.c_str());

  // Regularise the Pdf
  const RooArgSet* obsSet = fData->get();

  coutI(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") started writing local editing command for PDF " << fPdf->GetName() << endl;
  stringstream editStr;
  editStr << "EDIT::" << fPdf->GetName() << "(" << fPdf->GetName();
  int nrEdits = 0;

  // RooArgSet* tmpGlobalObservables = new RooArgSet(fGlobalObservables->GetName());
  // RooArgSet* tmpObservables = new RooArgSet(fObservables->GetName());
  RooArgSet* tmpGlobalObservables = new RooArgSet("GlobalObservables");
  RooArgSet* tmpObservables = new RooArgSet("Observables");

  for (map<string,string>::iterator itr = thisRenamingMap.begin(); itr != thisRenamingMap.end(); ++itr) {
    std::string oldObsName = itr->first;
    std::string newObsName = itr->second;

    std::string oldGlobsName = fRenamingMap.GetAttribute(oldObsName, RenamingMap::GlobalObservable, RenamingMap::individual);
    std::string newGlobsName = fRenamingMap.GetAttribute(newObsName, RenamingMap::GlobalObservable, RenamingMap::combined);

    std::string oldConstraintName = fRenamingMap.GetAttribute(oldObsName, RenamingMap::Constraint, RenamingMap::individual);
    std::string newConstraintName = fRenamingMap.GetAttribute(newObsName, RenamingMap::Constraint, RenamingMap::combined);

    std::string newConstraintType = fRenamingMap.GetAttribute(newObsName, RenamingMap::Type, RenamingMap::combined);

    // Verify that the new object exists ...
    if (!fWorkSpace->var(newObsName.c_str()) && !fWorkSpace->obj(newObsName.c_str())) {
      coutE(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") variable " << newObsName  << " was not imported into workspace " << fWorkSpace->GetName() << ". Skipping it." << endl;
      continue;
    }

    // ... and that the renamed object exists, too.
    std::string renamedVar = oldObsName + "_" + fParentMeasurement;
    if (!fWorkSpace->obj(renamedVar.c_str())) {
      coutE(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") object " << renamedVar << " was not imported into workspace " << fWorkSpace->GetName() << ". Skipping it." << endl;
      continue;
    }

    // Don't edit if the Pdf doesn't use the variable.
    if (!fPdf->dependsOn(*(RooAbsArg*)fWorkSpace->obj(renamedVar.c_str())))  {
      coutW(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") skipping variable " << renamedVar << endl;
      continue;
    }

    // Write the edit command for the current parameter
    if (newConstraintType == "Gaussian" || newConstraintType == "Poisson") {
      // global observable
      editStr << "," << oldGlobsName + "_" + fParentMeasurement << "=" << newGlobsName;
      // nuisance parameter
      editStr << "," << oldObsName + "_" + fParentMeasurement << "=" << newObsName;
      // constraint Pdf
      editStr << "," << oldConstraintName + "_" + fParentMeasurement  << "=" << newConstraintName;
      // Specify new global observable
      RooRealVar* tmpglob = fWorkSpace->var(newGlobsName.c_str());
      tmpglob->setAttribute("GLOBAL_OBSERVABLE");
      tmpGlobalObservables->add(*tmpglob);
    } else if (newConstraintType == "unconstrained") {
      // (unconstrained) nuisance parameter
      editStr << "," << oldObsName + "_" + fParentMeasurement << "=" << newObsName;
    } else {
      coutE(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") did not recognize constraint type " << newConstraintType << ". Skipping." << endl;
      continue;
    }
    nrEdits++;
  }

  // Loop over the nuisance parameters for this pdf and add them for the combined model
  // RooArgSet* tmpNuisanceParameters = new RooArgSet(fNuisanceParameters->GetName());
  RooArgSet* tmpNuisanceParameters = new RooArgSet("NuisanceParameters");
  TIterator* nuiItr = fNuisanceParameters->createIterator();
  RooRealVar* nextNui;
  while ((nextNui = (RooRealVar*)nuiItr->Next())) {
    stringstream argName;

    if (thisRenamedNuis[string(nextNui->GetName())] == "")  {
      argName << nextNui->GetName() << "_" + fParentMeasurement;
    } else {
      argName << thisRenamedNuis[string(nextNui->GetName())];
    }

    coutI(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") trying to get nuisance parameter " << argName.str() << endl;
    RooRealVar* nextNui_store = nextNui;
    nextNui = (RooRealVar*)fWorkSpace->obj(argName.str().c_str());

    if (!nextNui) {
      if (fPdf->dependsOn(*nextNui_store)) {
        coutF(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") couldn't find nuisance parameter " << argName.str() << " in workspace " << fWorkSpace->GetName() << " for channel " << fName << " from measurement " << fParentMeasurement << endl;
        exit(-1);
      } else {
        coutW(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") nuisance parameter " << argName.str() <<" in list of parameters, but Pdf " << fPdf->GetName() << "from channel " << fName << " has no dependence. Ignoring it." << endl;
        continue;
      }
    }

    if (nextNui->isFundamental()) {
      coutI(ObjectHandling) << "Channel::RegulariseChannel(" << fName <<  ") adding nuisance parameter: " << nextNui->GetName() << " from channel " << fName << " from measurement " << fParentMeasurement << endl;
      tmpNuisanceParameters->add(*nextNui);
    }
  }
  delete nuiItr;
  SetNuisanceParameters(tmpNuisanceParameters);

  // loop over the observables in the dataset for this channel and add to the editing command
  TIterator* obsSetItr = obsSet->createIterator();
  RooRealVar* nextObs;
  while ((nextObs = (RooRealVar*)obsSetItr->Next())) {
    RooRealVar* fromComb = (RooRealVar*)fWorkSpace->var((string(nextObs->GetName()) + "_" + fParentMeasurement).c_str());
    // RooAbsArg* fromComb = (RooAbsArg*)fWorkSpace->var((string(nextObs->GetName())).c_str());
    if (!fromComb) {
      coutW(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") variable " << string(nextObs->GetName()) << " is not in workspace " << fWorkSpace->GetName() << "." << endl;
      continue;
    }

    if (!fPdf->dependsOn(*fromComb)) continue;

    fWorkSpace->import(*nextObs, Silence());
    tmpObservables->add(*(RooAbsArg*)fWorkSpace->obj(nextObs->GetName()));
    editStr << "," << nextObs->GetName() << "_" + fParentMeasurement + "=" << nextObs->GetName();
    nrEdits++;

    // fData->changeObservableName(nextObs->GetName(), fromComb->GetName());
    // cout << "adding " << fromComb->GetName() << endl;
    // tmpObservables->add(*fromComb);
  }
  delete obsSetItr;
  SetObservables(tmpObservables);

  // that should be all to edit for now
  editStr << ")";
  coutP(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") editing command so far " << editStr.str() << endl;

  // perform the editing
  if (nrEdits != 0) {
    std::string originalName = fPdf->GetName();
    fWorkSpace->factory(editStr.str().c_str());
    fPdf = fWorkSpace->pdf(originalName.c_str());
    if (!fPdf) {
      coutF(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") something went wrong when editing the Pdf for channel " << fName << "." << endl;
      exit(-1);
    }
  }

  // taking care of the global observables
  coutP(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") adding global observables for model " << fParentMeasurement << endl;

  TIterator* obsItr = fGlobalObservables->createIterator();
  RooRealVar* nextGlob;
  while ((nextGlob = (RooRealVar*)obsItr->Next())) {
    stringstream argName;
    if (thisRenamedGlobs[string(nextGlob->GetName())] == "")  {
      argName << nextGlob->GetName() << "_" + fParentMeasurement;
    } else {
      argName << thisRenamedGlobs[string(nextGlob->GetName())];
    }

    RooRealVar* nextGlob_store = nextGlob;

    // grab it from the ws
    nextGlob = (RooRealVar*)fWorkSpace->obj(argName.str().c_str());
    if (!nextGlob) {
      if (fPdf->dependsOn(*nextGlob_store)) {
        coutF(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") couldn't find global observable: " << argName.str() << " in workspace for model " << fParentMeasurement << endl;
        exit(-1);
      } else {
        coutW(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") global observable in list of global observables, but pdf has no dependence. Ignoring it." << endl;
        continue;
      }
    }

    coutI(ObjectHandling) << "Channel::RegulariseChannel(" << fName << ") adding global observable " << nextGlob->GetName() << " from model " << fParentMeasurement << endl;
    nextGlob->setConstant(true);
    nextGlob->setAttribute("GLOBAL_OBSERVABLE");
    tmpGlobalObservables->add(*nextGlob);
  }
  delete obsItr;
  SetGlobalObservables(tmpGlobalObservables);
  fGlobalObservables->setAttribAll("GLOBAL_OBSERVABLE");
}
