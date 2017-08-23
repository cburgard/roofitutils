#include "RooFitUtils/RenamingMap.h"
#include "RooFitUtils/CorrelationScheme.h"
#include "RooFitUtils/AbsMeasurement.h"

#include "RooGaussian.h"

const std::string RenamingMap::AttributeNames[] = { "Observable", "ObservableRange", "GlobalObservable", "GlobalObservableRange", "Sigma", "SigmaRange", "Constraint", "Type" };
const std::string RenamingMap::ConstraintTypeNames[] = { "automatic", "unconstrained", "Gaussian", "Poisson", "Lognormal", "unknown" };

// ____________________________________________________________________________|__________
// Constqructor
RenamingMap::RenamingMap( std::string MapName )
  :
  TNamed( MapName.c_str(), MapName.c_str() )
{
  coutP(InputArguments) << "RenamingMap::RenamingMap(" << fName <<") created" << endl;
}

// ____________________________________________________________________________|__________
// Destructor
RenamingMap::~RenamingMap()
{

}

// ____________________________________________________________________________|__________
// Generic interface for specifying a parameter that should be renamed
void RenamingMap::RenameParameter( std::string OldParameterName, std::string NewParameterName )
{
  // Check if the parameter exists already in the renaming map
  if (fRenamingMap.find(OldParameterName) != fRenamingMap.end()) {
    coutF(InputArguments) << "RenamingMap::AddRenamedParameter(" << fName <<") " << OldParameterName << " already exists! [Hint: renamed to " << fRenamingMap[OldParameterName] << "]" << endl;
    exit(-1);
  }

  // Add to renaming map
  fRenamingMap[OldParameterName] = NewParameterName;
  SetAttribute(OldParameterName, Observable, OldParameterName, individual);
  SetAttribute(NewParameterName, Observable, NewParameterName, combined);

  coutI(InputArguments) << "RenamingMap::RenameParameter(" << fName <<") " << OldParameterName << " will be renamed to " << NewParameterName << endl;
}

// ____________________________________________________________________________|__________
// Interface to set an attribute for a given parameter
void RenamingMap::SetAttribute( std::string ParameterName, Attribute thisAttribute, std::string AttributeValue, MeasurementType thisMeasurementType )
{
  // std::string thisAttributeName = AttributeNames[thisAttribute];

  switch (thisMeasurementType) {
    case individual: {
      fOldParameterMap[ParameterName][thisAttribute] = AttributeValue;
      ccoutD(InputArguments) << "RenamingMap::SetAttribute(" << fName <<") " << AttributeNames[thisAttribute] << " = " << AttributeValue << " for parameter " << ParameterName << " for individual measurement" << endl;
      break;
    }
    case combined: {
      fNewParameterMap[ParameterName][thisAttribute] = AttributeValue;
      ccoutD(InputArguments) << "RenamingMap::SetAttribute(" << fName <<") " << AttributeNames[thisAttribute] << " = " << AttributeValue << " for parameter " << ParameterName << " for combined measurement" << endl;
      break;
    }
    default: {
      coutF(InputArguments) << "RenamingMap::SetAttribute(" << fName <<") " << thisMeasurementType << " not recognised. Ignoring attribute " << AttributeNames[thisAttribute] << " = " << AttributeValue << " for parameter" << ParameterName << endl;
      exit(-1);
    }
  }
}

// ____________________________________________________________________________|__________
// Interface to retrieve an attribute for a given parameter
std::string RenamingMap::GetAttribute( std::string ParameterName, Attribute thisAttribute, MeasurementType thisMeasurementType )
{
  std::string AttributeValue = "";

  switch (thisMeasurementType) {
    case individual: {
      AttributeValue = fOldParameterMap[ParameterName][thisAttribute];
      break;
    }
    case combined: {
      AttributeValue = fNewParameterMap[ParameterName][thisAttribute];
      break;
    }
    default: {
      coutF(InputArguments) << "RenamingMap::GetAttribute(" << fName <<") " << thisMeasurementType << " not recognised. Ignoring attribute " << AttributeNames[thisAttribute] << " for parameter" << ParameterName << endl;
      exit(-1);
    }
  }

  if (AttributeValue == "") {
    coutW(InputArguments) << "RenamingMap::GetAttribute(" << fName <<") " << AttributeNames[thisAttribute] << " not found for parameter " << ParameterName << endl;
  }

  return AttributeValue;
}

// ____________________________________________________________________________|__________
// Automatically add missing attributes to the renaming map, overwrite
// explicitely specified old attributes, like constraint names and associated
// global observables
void RenamingMap::AddAttributes( ModelConfig* tmpModelConfig )
{
  // Load information needed to determine attributes from ModelConfig
  RooWorkspace* tmpWorkspace = tmpModelConfig->GetWorkspace();
  RooAbsPdf* tmpPdf = (RooAbsPdf*)tmpModelConfig->GetPdf();
  RooArgSet* tmpAllNuisanceParameters = (RooArgSet*) tmpModelConfig->GetNuisanceParameters();
  RooArgSet* tmpAllObservables = (RooArgSet*) tmpModelConfig->GetObservables();
  RooArgSet* tmpAllGlobalObservables = (RooArgSet*) tmpModelConfig->GetGlobalObservables();

  // Copies, to keep original sets in place after getAllconstraints call
  RooArgSet tmpAllNuisanceParameters2 = *tmpAllNuisanceParameters;
  RooArgSet tmpAllObservables2 = *tmpAllObservables;
  RooArgSet* AllConstraints = tmpPdf->getAllConstraints(tmpAllObservables2, tmpAllNuisanceParameters2, kFALSE);

  // Take care of the case where we have a product of constraint terms
  TIterator* ConstraintItr = AllConstraints->createIterator();
  RooAbsArg* nextConstraint;
  RooArgSet* tmpAllConstraints = new RooArgSet(AllConstraints->GetName());
  while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next())) {
    if (nextConstraint->IsA() == RooProdPdf::Class()) {
      RooArgSet thisComponents;
      FindUniqueProdComponents((RooProdPdf*)nextConstraint, thisComponents);
      tmpAllConstraints->add(thisComponents);
    } else {
      tmpAllConstraints->add(*nextConstraint);
    }
  }

  // AllConstraints->Print();
  // tmpAllConstraints->Print();

  std::set< std::string > notExistingParameters;

  for (std::map< std::string, std::string >::iterator it = fRenamingMap.begin(); it != fRenamingMap.end(); ++it) {
    std::string thisOldObservableName = it->first;
    std::string thisNewObservableName = it->second;

    // old attributes
    std::string thisOldConstraintName        = GetAttribute(thisOldObservableName, RenamingMap::Constraint, RenamingMap::individual);
    std::string thisOldObservableRange       = GetAttribute(thisOldObservableName, RenamingMap::ObservableRange, RenamingMap::individual);
    std::string thisOldGlobalObservableName  = GetAttribute(thisOldObservableName, RenamingMap::GlobalObservable, RenamingMap::individual);
    std::string thisOldGlobalObservableRange = GetAttribute(thisOldObservableName, RenamingMap::GlobalObservableRange, RenamingMap::individual);
    // std::string thisOldConstraintType        = GetAttribute(thisOldObservableName, RenamingMap::Type, RenamingMap::individual);

    // new attributes
    std::string thisNewConstraintName        = GetAttribute(thisNewObservableName, RenamingMap::Constraint, RenamingMap::combined);
    std::string thisNewObservableRange       = GetAttribute(thisNewObservableName, RenamingMap::ObservableRange, RenamingMap::combined);
    std::string thisNewGlobalObservableName  = GetAttribute(thisNewObservableName, RenamingMap::GlobalObservable, RenamingMap::combined);
    std::string thisNewGlobalObservableRange = GetAttribute(thisNewObservableName, RenamingMap::GlobalObservableRange, RenamingMap::combined);
    std::string thisNewConstraintType        = GetAttribute(thisNewObservableName, RenamingMap::Type, RenamingMap::combined);
    std::string thisNewSigmaName             = GetAttribute(thisNewObservableName, RenamingMap::Sigma, RenamingMap::combined);
    std::string thisNewSigmaRange            = GetAttribute(thisNewObservableName, RenamingMap::SigmaRange, RenamingMap::combined);

    // auto attributes
    std::string tmpConstraintName        = "";
    std::string tmpObservableRange       = "";
    std::string tmpGlobalObservableName  = "";
    std::string tmpGlobalObservableRange = "";
    std::string tmpSigmaRange            = "";
    std::string tmpConstraintType        = "unconstrained";

    // Get parameter from workspace
    RooRealVar* nextNuisanceParameter = tmpWorkspace->var(thisOldObservableName.c_str());
    if (!nextNuisanceParameter) {
      coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") observable " << thisOldObservableName <<" not found. Skipping." << endl;
      notExistingParameters.insert(thisOldObservableName);
      continue;
    }

    // find range of nuisance parameter in workspace
    double oldNuisVal = nextNuisanceParameter->getVal();
    double oldNuisMin = nextNuisanceParameter->getMin();
    double oldNuisMax = nextNuisanceParameter->getMax();

    stringstream tmpNuisSS;
    tmpNuisSS << setprecision(30) << "[" << oldNuisVal << "," << oldNuisMin <<"," << oldNuisMax << "]";
    tmpObservableRange = tmpNuisSS.str();

    // loop over all constraint of pdf to determine
    // constraint type of nuisance parameter
    TIterator* ConstraintItr = tmpAllConstraints->createIterator();
    RooAbsArg* nextConstraint;
    bool foundConstraint = kFALSE;
    while ((nextConstraint = (RooAbsArg*)ConstraintItr->Next()) && !foundConstraint) {
      if (nextConstraint->dependsOn(*nextNuisanceParameter)) {
        foundConstraint = kTRUE;

        // find name of constraint
        tmpConstraintName = string(nextConstraint->GetName());
        TString thisConstraintType = nextConstraint->ClassName();
        tmpConstraintType = (thisConstraintType.ReplaceAll("Roo", "")).Data();

        // loop over global observables to match nuisance parameter and
        // global observable in case of a constrained nuisnace parameter
        TIterator* GlobsItr = tmpAllGlobalObservables->createIterator();
        RooRealVar* nextGlobalObservable;
        bool foundGlobalObservable = kFALSE;
        while ((nextGlobalObservable = (RooRealVar*)GlobsItr->Next()) && !foundGlobalObservable) {
          if (nextConstraint->dependsOn(*nextGlobalObservable)) {
            foundGlobalObservable = kTRUE;

            // find name of globale observable
            tmpGlobalObservableName = string(nextGlobalObservable->GetName());

            // find range of global observable
            double oldGlobVal = nextGlobalObservable->getVal();
            // double oldGlobMin = nextGlobalObservable->getMin();
            // double oldGlobMax = nextGlobalObservable->getMax();

            stringstream tmpGlobSS;
            tmpGlobSS << setprecision(30) << "[" << oldGlobVal /* << "," << oldGlobMin <<"," << oldGlobMax */ << "]";
            tmpGlobalObservableRange = tmpGlobSS.str();

            // find constraint width in case of a Gaussian
            if (nextConstraint->IsA() == RooGaussian::Class()) {
              double oldSigmaVal = 1.0;
              TIterator* ServerItr = nextConstraint->serverIterator();
              RooRealVar* nextServer;
              bool foundSigma = kFALSE;
              while ((nextServer = (RooRealVar*)ServerItr->Next()) && !foundSigma) {
                if (nextServer != nextGlobalObservable && nextServer != nextNuisanceParameter) {
                  oldSigmaVal = nextServer->getVal();
                  foundSigma = kTRUE;
                }
              }

              stringstream tmpSigmaSS;
              tmpSigmaSS << setprecision(30) << oldSigmaVal;
              tmpSigmaRange = tmpSigmaSS.str();
              if (AbsMeasurement::AlmostEqualUlpsAndAbs(oldSigmaVal, 1.0, 0.001, 4)) {
                tmpSigmaRange = "1.0";
              }

              if (!foundSigma) {
                coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") sigma for pdf " << nextConstraint->GetName() << " not found. Using 1.0." << endl;
              } else {
                coutI(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") using " << tmpSigmaRange << " for sigma of pdf " << nextConstraint->GetName() << endl;
              }
            }
          }
        }
        delete GlobsItr;
      }
    }
    delete ConstraintItr;

    // individual observable range
    //   can be specified when defining the correlation scheme
    //   not explicitely needed for regularisation
    //   use automatic value in any case, but print warning if it does not match the specified value
    if (thisOldObservableRange != "" && thisOldObservableRange != tmpObservableRange) {
      coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") specified old observable range (" << thisOldObservableRange <<") does not match actual range (" << tmpObservableRange << "). Use that instead." << endl;
    }
    SetAttribute(thisOldObservableName, RenamingMap::ObservableRange, tmpObservableRange, RenamingMap::individual);

    // individual global observable name
    //   can be specified when defining the correlation scheme
    //   set only if constraint type is not unconstrained
    //   use automatic value in any case, but print warning if it does not match the specified value
    if (tmpConstraintType != "unconstrained") {
      if (thisOldGlobalObservableName != "" && thisOldGlobalObservableName != tmpGlobalObservableName) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") specified old global observable name (" << thisOldGlobalObservableName <<") does not match actual name (" << tmpGlobalObservableName << "). Use that instead." << endl;
      }
      SetAttribute(thisOldObservableName, RenamingMap::GlobalObservable, tmpGlobalObservableName, RenamingMap::individual);
    }

    // individual global observable range
    //   can be specified when defining the correlation scheme
    //   not explicitely needed for regularisation
    //   set only if constraint type is not uncosntrained
    //   use automatic value in any case, but print warning if it does not match the specified value
    if (tmpConstraintType != "unconstrained") {
      if (thisOldGlobalObservableRange != "" && thisOldGlobalObservableRange != tmpGlobalObservableRange) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") specified old global observable range (" << thisOldGlobalObservableRange <<") does not match actual range (" << tmpGlobalObservableRange << "). Use that instead." << endl;
      }
      SetAttribute(thisOldObservableName, RenamingMap::GlobalObservableRange, tmpGlobalObservableRange, RenamingMap::individual);
    }

    // individual constraint name
    //   can be specified when defining the correlation scheme
    //   set only if constraint type is not uncosntrained
    //   use automatic value in any case, but print warning if it does not match the specified value
    if (tmpConstraintType != "unconstrained") {
      if (thisOldConstraintName != "" && thisOldConstraintName != tmpConstraintName) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") specified old constraint name (" << thisOldConstraintName <<") does not match actual name (" << tmpConstraintName << "). Use that instead." << endl;
      }
      SetAttribute(thisOldObservableName, RenamingMap::Constraint, tmpConstraintName, RenamingMap::individual);
    }

    // individual constraint type
    //   no way of specifying by hand
    //   not explicitely needed for regularisation
    //   use automatic value in any case
    SetAttribute(thisOldObservableName, RenamingMap::Type, tmpConstraintType, RenamingMap::individual);

    // combined observable range
    //   can be specified when defining the correlation scheme
    //   if not defined, set it according to constraint type
    std::string tmpNewObservableRange = "";
    if (thisNewObservableRange == "") {
      if (thisNewConstraintType == "Gaussian" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Gaussian")) {
        tmpNewObservableRange = "[0.0,-5.0,5.0]";
      } else if (thisNewConstraintType == "Poisson" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Poisson")) {
        // tmpNewObservableRange = "[1.0,0.0,10.0]";
        tmpNewObservableRange = tmpObservableRange; // Potentially problematic in case the parameter is appearing in multiple measurements, but with different values. In that case need to write the value down explicitly
      } else if (thisNewConstraintType == "Lognormal" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Lognormal")) {
        tmpNewObservableRange = "[1.0,0.0,10.0]";
      } else if (thisNewConstraintType == "unconstrained" || (thisNewConstraintType == "automatic" && tmpConstraintType == "unconstrained")) {
        tmpNewObservableRange = tmpObservableRange;
      }
      SetAttribute(thisNewObservableName, RenamingMap::ObservableRange, tmpNewObservableRange, RenamingMap::combined);
    } else {
      if (thisNewObservableRange != tmpObservableRange) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") observable range of " << thisNewObservableName << " changed from " << tmpObservableRange << " to " << thisNewObservableRange << endl;
      }
    }

    // combined global observable name
    //   can be specified when defining the correlation scheme
    //   prefix nuisance parameter with specified string
    if (thisNewGlobalObservableName == "") {
      if (thisNewConstraintType != "unconstrained" && (thisNewConstraintType == "automatic" && tmpConstraintType != "unconstrained")) {
        std::string tmpNewGlobalObservableName = "nom_" + thisNewObservableName;
        SetAttribute(thisNewObservableName, RenamingMap::GlobalObservable, tmpNewGlobalObservableName, RenamingMap::combined);
      }
    } else {
      if (thisNewGlobalObservableName != tmpGlobalObservableName) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") global observable name of " << thisNewObservableName << " changed from " << tmpGlobalObservableName << " to " << thisNewGlobalObservableName << endl;
      }
    }

    // combined global observable range
    //   can be specified when defining the correlation scheme
    //   if not defined, set it according to the constraint type
    if (thisNewGlobalObservableRange == "") {
      std::string tmpNewGlobalObservableRange = "";
      if (thisNewConstraintType == "Gaussian" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Gaussian")) {
        tmpNewGlobalObservableRange = "[0.0]";
      } else if (thisNewConstraintType == "Poisson" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Poisson")) {
        // tmpNewGlobalObservableRange = "[1.0]";
        tmpNewGlobalObservableRange = tmpGlobalObservableRange; // Potentially problematic in case the parameter is appearing in multiple measurements, but with different values. In that case need to write the value down explicitly
      } else if (thisNewConstraintType == "Lognormal" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Lognormal")) {
        tmpNewGlobalObservableRange = "[1.0]";
      }
      if (thisNewConstraintType != "unconstrained" && (thisNewConstraintType == "automatic" && tmpConstraintType != "unconstrained")) {
        SetAttribute(thisNewObservableName, RenamingMap::GlobalObservableRange, tmpNewGlobalObservableRange, RenamingMap::combined);
      }
    } else {
      if (thisNewGlobalObservableRange != tmpGlobalObservableRange) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") global observable range of " << thisNewObservableName << " changed from " << tmpGlobalObservableRange << " to " << thisNewGlobalObservableRange << endl;
      }
    }

    // combined constraint name
    //   can be specified when defining the correlation scheme
    //   postfix nuisance parameter with specified string
    if (thisNewConstraintName == "") {
      if (thisNewConstraintType != "unconstrained" && (thisNewConstraintType == "automatic" && tmpConstraintType != "unconstrained")) {
        std::string tmpNewConstraintName = thisNewObservableName + "Constraint";
        SetAttribute(thisNewObservableName, RenamingMap::Constraint, tmpNewConstraintName, RenamingMap::combined);
      }
    } else {
      if (thisNewConstraintName != tmpConstraintName) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") constraint name of " << thisNewObservableName << " changed from " << tmpConstraintName << " to " << thisNewConstraintName << endl;
      }
    }

    // combined constraint type
    //   can be specified when defining the correlation scheme, default automatic
    //   if automatic, use the mapped individual constraint type
    //   otherwise, use the explicit type, but print warning that it is different
    if (thisNewConstraintType == "automatic") {
      SetAttribute(thisNewObservableName, RenamingMap::Type, tmpConstraintType, RenamingMap::combined);
    } else {
      if (thisNewConstraintType != tmpConstraintType) {
        coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") casting constraint type of " << thisNewObservableName << " from " << tmpConstraintType << " to " << thisNewConstraintType << endl;
      }
    }

    // combined sigma name
    //   only used for Gaussian constraints
    //   can only be specified explicitely, otherwise default values ("") will be used
    if (thisNewSigmaName == "") {
      if (thisNewConstraintType == "Gaussian" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Gaussian")) {
        std::string tmpNewSigmaName = "";
        SetAttribute(thisNewObservableName, RenamingMap::Sigma, tmpNewSigmaName, RenamingMap::combined);
      }
    }

    // combined sigma range
    //   only used for Gaussian constraints
    //   can only be specified explicitely, otherwise default value ("[1.0]") will be used
    if (thisNewSigmaRange == "") {
      if (thisNewConstraintType == "Gaussian" || (thisNewConstraintType == "automatic" && tmpConstraintType == "Gaussian")) {
        std::string tmpNewSigmaRange = tmpSigmaRange;
        SetAttribute(thisNewObservableName, RenamingMap::SigmaRange, tmpNewSigmaRange, RenamingMap::combined);
      }
    }
  }

  // Clean up spurious parameters
  for (std::set< std::string >::iterator it = notExistingParameters.begin(); it != notExistingParameters.end(); ++it) {
    std::string thisOldObservableName = *it;
    coutW(InputArguments) << "RenamingMap::AddAttributes(" << fName << ") removing spurious parameter " << thisOldObservableName << " from RenamingMap and CorrelationScheme" << endl;
    fOldParameterMap.erase(thisOldObservableName);
    fNewParameterMap.erase(fRenamingMap[thisOldObservableName]);
    fRenamingMap.erase(thisOldObservableName);
  }
}

// ____________________________________________________________________________|__________
// Format constraint terms in factory language
std::string RenamingMap::FactoryExpression(std::string ParameterName, RenamingMap::MeasurementType thisMeasurementType)
{
  std::string thisObservableName        = GetAttribute(ParameterName, Observable, thisMeasurementType);
  std::string thisObservableRange       = GetAttribute(ParameterName, ObservableRange, thisMeasurementType);
  std::string thisGlobalObservableName  = GetAttribute(ParameterName, GlobalObservable, thisMeasurementType);
  std::string thisGlobalObservableRange = GetAttribute(ParameterName, GlobalObservableRange, thisMeasurementType);
  std::string thisConstraintName        = GetAttribute(ParameterName, Constraint, thisMeasurementType);
  std::string thisConstraintType        = GetAttribute(ParameterName, Type, thisMeasurementType);
  std::string thisSigmaName             = GetAttribute(ParameterName, Sigma, thisMeasurementType);
  std::string thisSigmaRange            = GetAttribute(ParameterName, SigmaRange, thisMeasurementType);

  // strip "[" and "]" from sigma range, if not variable name given
  if (thisSigmaName == "") {
    TString tmpSigmaRange = thisSigmaRange.c_str();
    tmpSigmaRange.ReplaceAll("[", "");
    tmpSigmaRange.ReplaceAll("]", "");
    thisSigmaRange = tmpSigmaRange.Data();
  }

  std::string thisExpression = "";

  if (thisConstraintType == "unconstrained") {
    thisExpression = TString::Format("%s%s", thisObservableName.c_str(), thisObservableRange.c_str()).Data();
  } else if (thisConstraintType == "Gaussian") {
    thisExpression = TString::Format("Gaussian::%s(%s%s,%s%s,%s%s)", thisConstraintName.c_str(), thisObservableName.c_str(), thisObservableRange.c_str(), thisGlobalObservableName.c_str(), thisGlobalObservableRange.c_str(), thisSigmaName.c_str(), thisSigmaRange.c_str()).Data();
  } else if (thisConstraintType == "Poisson") {
    // no rounding as default implemented
    thisExpression = TString::Format("Poisson::%s(%s%s,Product::%spoisMean({%s%s,%stau%s}),1)", thisConstraintName.c_str(), thisGlobalObservableName.c_str(), thisGlobalObservableRange.c_str(), thisObservableName.c_str(), thisObservableName.c_str(), thisObservableRange.c_str(), thisObservableName.c_str(), thisGlobalObservableRange.c_str()).Data();
  } else if (thisConstraintType == "Lognormal") {
    coutW(InputArguments) << "RenamingMap::FactoryExpression(" << fName <<") " << thisConstraintType << " not implemented yet." << endl;
  } else {
    coutF(InputArguments) << "RenamingMap::FactoryExpression(" << fName <<") " << thisConstraintType << " not recognised for parameter " << ParameterName << endl;
    exit(-1);
  }

  if (thisExpression == "") {
    coutE(InputArguments) << "RenamingMap::FactoryExpression(" << fName <<") failed writing factory expression for parameter " << ParameterName << endl;
  }

  return thisExpression;
}

// ____________________________________________________________________________|__________
// Find all unique components of a RooProdPdf
void RenamingMap::FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components )
{
  static int counter = 0;
  counter++;

  if (counter > 50) {
    coutE(ObjectHandling) << "RenamingMap::FindUniqueProdComponents(" << fName << ") detected infinite loop. Please check." << endl;
    exit(1);
  }

  RooArgList pdfList = Pdf->pdfList();
  if (pdfList.getSize() == 1) {
    coutI(ObjectHandling) << "RenamingMap::FindUniqueProdComponents(" << fName << ") " << pdfList.at(0)->GetName() << " is fundamental." << endl;
    Components.add(pdfList);
  } else {
    TIterator* pdfItr = pdfList.createIterator();
    RooAbsArg* nextArg;
    while ((nextArg = (RooAbsArg*)pdfItr->Next())) {
      RooProdPdf* Pdf = (RooProdPdf*)nextArg;
      if (string(Pdf->ClassName()) != "RooProdPdf") {
        coutI(ObjectHandling) << "RenamingMap::FindUniqueProdComponents(" << fName << ") " << Pdf->GetName() << " is no RooProdPdf. Adding it." << endl;
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
// Print a RenamingMap
void RenamingMap::Print()
{
  int nrColumns = 1;
  string* header = new string[nrColumns];
  header[0] = "New name";

  int nrEntries = fRenamingMap.size();
  if (nrEntries == 0) {
    coutW(ObjectHandling) << "RenamingMap::Print(" << fName << ") " << fName << " is empty" << endl;
    delete[] header;
    return;
  }
  string* firstCol = new string[nrEntries+1];
  firstCol[0] = "Old name";
  std::string** matrix = new string*[nrEntries];
  for (int in = 0; in < nrEntries; in++) {
    matrix[in] = new std::string[nrColumns];
    for (int i = 0; i < nrColumns; i++) {
      matrix[in][i] = "";
    }
  }

  int irow = 1;
  for (std::map<string,string>::iterator it = fRenamingMap.begin(); it != fRenamingMap.end(); ++it) {
    std::string oldObservableName = it->first;
    firstCol[irow] = oldObservableName;

    std::string newObservableName = it->second;
    int icol = 0;
    matrix[irow-1][icol] = newObservableName;
    irow++;
  }
  CorrelationScheme::PrintTable(firstCol, matrix, NULL, header, nrEntries, nrColumns, 0, cout, "    ");

  delete[] header;
  delete[] firstCol;

  for (int in = 0; in < nrEntries; in++) {
    delete[] matrix[in];
  }

  delete[] matrix;
}
