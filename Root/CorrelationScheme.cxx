#include "RooFitUtils/CorrelationScheme.h"

// ____________________________________________________________________________|__________
// Constructor
RooFitUtils::CorrelationScheme::CorrelationScheme( std::string SchemeName, std::string ParametersOfInterest, bool AutoCorrelation )
  :
  TNamed(SchemeName.c_str(), SchemeName.c_str()),
  fParametersOfInterest(ParametersOfInterest),
  fAutoCorrelation(AutoCorrelation)
{
  coutP(InputArguments) << "CorrelationScheme::CorrelationScheme(" << fName <<") created" << endl;
  if (fAutoCorrelation) {
    cout << "\n****************************************************************************" << endl;
    cout << "* WARNING: Parameters with the same name will be correlated automatically! *" << endl;
    cout << "*          It is the responsibility of the user to make sure, that those   *" << endl;
    cout << "*          parameters describe the same effect in the input workspaces!    *" << endl;
    cout << "*                                                                          *" << endl;
    cout << "*          It is recommended to explicitly enter the correlation scheme!   *" << endl;
    cout << "****************************************************************************\n" << endl;
  }
}

// ____________________________________________________________________________|__________
// Destructor
RooFitUtils::CorrelationScheme::~CorrelationScheme()
{
  {}
}

// ____________________________________________________________________________|__________
// Enable automatic correlation, print warning, as potentially dangerous
void RooFitUtils::CorrelationScheme::SetAutoCorrelation( bool setting )
{
  fAutoCorrelation = setting;

  if (fAutoCorrelation) {
    cout << "\n****************************************************************************" << endl;
    cout << "* WARNING: Parameters with the same name will be correlated automatically! *" << endl;
    cout << "*          It is the responsibility of the user to make sure, that those   *" << endl;
    cout << "*          parameters describe the same effect in the input workspaces!    *" << endl;
    cout << "*                                                                          *" << endl;
    cout << "*          It is recommended to explicitly enter the correlation scheme!   *" << endl;
    cout << "****************************************************************************\n" << endl;
  }
}

// ____________________________________________________________________________|__________
// Interface to add fully correlated parameters among different channels,
// constraint terms will be determined automatically if not specified otherwise
void RooFitUtils::CorrelationScheme::CorrelateParameter( const char* OldParameterNamePlusMeasurement, const char* NewParameterName, RenamingMap::ConstraintType thisConstraintType )
{
  TString allOldParameterNamePlusMeasurement = OldParameterNamePlusMeasurement;
  TObjArray* allOldParameterNamePlusMeasurementArray = allOldParameterNamePlusMeasurement.Tokenize(",");
  unsigned int numCorrPars = allOldParameterNamePlusMeasurementArray->GetEntries();

  for (unsigned int itrPars = 0; itrPars < numCorrPars; ++itrPars) {
    TString thisOldParameterNamePlusMeasurement = ((TObjString*)allOldParameterNamePlusMeasurementArray->At(itrPars))->GetString();

    unsigned int offset = 0;
    while (thisOldParameterNamePlusMeasurement.Contains("(") && !thisOldParameterNamePlusMeasurement.Contains(")")) {
      offset++;
      thisOldParameterNamePlusMeasurement = thisOldParameterNamePlusMeasurement + "," + ((TObjString*)allOldParameterNamePlusMeasurementArray->At(itrPars+offset))->GetString();
    }

    if (thisOldParameterNamePlusMeasurement.Contains("(")) {
      assert(thisOldParameterNamePlusMeasurement.Contains(")"));
    }

    TObjArray* thisOldParameterNamePlusMeasurementArray = thisOldParameterNamePlusMeasurement.Tokenize("::");
    TString thisMeasurement = ((TObjString*)thisOldParameterNamePlusMeasurementArray->At(0))->GetString();
    TString thisOldParameter = ((TObjString*)thisOldParameterNamePlusMeasurementArray->At(1))->GetString();

    RenameParameter(thisMeasurement, thisOldParameter, NewParameterName, thisConstraintType);

    itrPars += offset;

    thisOldParameterNamePlusMeasurementArray->Delete();
    delete thisOldParameterNamePlusMeasurementArray;
  }

  allOldParameterNamePlusMeasurementArray->Delete();
  delete allOldParameterNamePlusMeasurementArray;
}

// ____________________________________________________________________________|__________
// Rename a parameter without correlating it
void RooFitUtils::CorrelationScheme::RenameParameter( const char* MeasurementName, const char* OldParameterName, const char* NewParameterName, RenamingMap::ConstraintType thisConstraintType )
{
  std::string thisOldParameterName         = string(OldParameterName);
  std::string thisOldConstraintName        = "";
  std::string thisOldObservableName        = "";
  std::string thisOldObservableRange       = "";
  std::string thisOldGlobalObservableName  = "";
  std::string thisOldGlobalObservableRange = "";
  std::string thisOldSigmaName             = "";
  std::string thisOldSigmaRange            = "";
  ParseInputs(thisOldParameterName, thisOldConstraintName, thisOldObservableName, thisOldObservableRange, thisOldGlobalObservableName, thisOldGlobalObservableRange, thisOldSigmaName, thisOldSigmaRange);

  std::string thisNewParameterName         = string(NewParameterName);
  std::string thisNewConstraintName        = "";
  std::string thisNewObservableName        = "";
  std::string thisNewObservableRange       = "";
  std::string thisNewGlobalObservableName  = "";
  std::string thisNewGlobalObservableRange = "";
  std::string thisNewSigmaName             = "";
  std::string thisNewSigmaRange            = "";
  ParseInputs(thisNewParameterName, thisNewConstraintName, thisNewObservableName, thisNewObservableRange, thisNewGlobalObservableName, thisNewGlobalObservableRange, thisNewSigmaName, thisNewSigmaRange);

  if (string(fCorrelationMap[MeasurementName].GetName()) == "") {
    fCorrelationMap[MeasurementName].SetName(MeasurementName);
  }

  fCorrelationMap[MeasurementName].RenameParameter(thisOldObservableName, thisNewObservableName);
  fCorrelationMap[MeasurementName].SetAttribute(thisNewObservableName, RenamingMap::Type, RenamingMap::ConstraintTypeNames[thisConstraintType], RenamingMap::combined);

  if (thisOldConstraintName != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisOldObservableName, RenamingMap::Constraint, thisOldConstraintName, RenamingMap::individual);
  }

  if (thisOldObservableRange != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisOldObservableName, RenamingMap::ObservableRange, thisOldObservableRange, RenamingMap::individual);
  }

  if (thisOldGlobalObservableName != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisOldObservableName, RenamingMap::GlobalObservable, thisOldGlobalObservableName, RenamingMap::individual);
  }

  if (thisOldGlobalObservableRange != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisOldObservableName, RenamingMap::GlobalObservableRange, thisOldGlobalObservableRange, RenamingMap::individual);
  }

  if (thisOldSigmaName != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisOldObservableName, RenamingMap::Sigma, thisOldSigmaName, RenamingMap::individual);
  }

  if (thisOldSigmaRange != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisOldObservableName, RenamingMap::SigmaRange, thisOldSigmaRange, RenamingMap::individual);
  }

  if (thisNewConstraintName != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisNewObservableName, RenamingMap::Constraint, thisNewConstraintName, RenamingMap::combined);
  }

  if (thisNewObservableRange != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisNewObservableName, RenamingMap::ObservableRange, thisNewObservableRange, RenamingMap::combined);
  }

  if (thisNewGlobalObservableName != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisNewObservableName, RenamingMap::GlobalObservable, thisNewGlobalObservableName, RenamingMap::combined);
  }

  if (thisNewGlobalObservableRange != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisNewObservableName, RenamingMap::GlobalObservableRange, thisNewGlobalObservableRange, RenamingMap::combined);
  }

  if (thisNewSigmaName != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisNewObservableName, RenamingMap::Sigma, thisNewSigmaName, RenamingMap::combined);
  }

  if (thisNewSigmaRange != "") {
    fCorrelationMap[string(MeasurementName)].SetAttribute(thisNewObservableName, RenamingMap::SigmaRange, thisNewSigmaRange, RenamingMap::combined);
  }
}

// ____________________________________________________________________________|__________
// Determine explicit names of constraint terms, observables and global observables
// as well as their ranges. The input format specified by the user should be
// ConstraintName(Observable[ObservableRange],GlobalObservable[GlobalObservableRange])
void RooFitUtils::CorrelationScheme::ParseInputs( std::string& InputName, std::string& InputConstraintName, std::string& InputObservableName, std::string& InputObservableRange, std::string& InputGlobalObservableName, std::string& InputGlobalObservableRange, std::string& InputSigmaName, std::string& InputSigmaRange )
{
  TString thisInputName = InputName.c_str();

  if (thisInputName.Contains("(")) {
    assert(thisInputName.Contains(")"));

    TObjArray* thisInputNameArray = thisInputName.Tokenize("(");

    // constraint name is part before the round opening bracket
    InputConstraintName = ((TObjString*)thisInputNameArray->At(0))->GetString();

    // arguments
    TString thisArguments = ((TObjString*)thisInputNameArray->At(1))->GetString();
    thisArguments.ReplaceAll(")", "");
    TObjArray* thisArgumentsArray;
    if (thisArguments.Contains("],")) {
      thisArguments.ReplaceAll("],", "|");
      thisArgumentsArray = thisArguments.Tokenize("|");
    } else {
      thisArgumentsArray = thisArguments.Tokenize(",");
    }

    int nrArguments = thisArgumentsArray->GetEntries();

    // find the observable including its range if specified
    if (nrArguments > 0) {
      TString thisInputObservableName = ((TObjString*)thisArgumentsArray->At(0))->GetString();
      std::string tmp = thisInputObservableName.Data();
      DecomposeVariable(tmp, InputObservableName, InputObservableRange);
    }

    // find the global observable including its range if specified
    if (nrArguments > 1) {
      TString thisInputGlobalObservableName = ((TObjString*)thisArgumentsArray->At(1))->GetString();
      std::string tmp = thisInputGlobalObservableName.Data();
      DecomposeVariable(tmp, InputGlobalObservableName, InputGlobalObservableRange);
    }

    // find the width/sigma including its range if specified, not used at the moment
    if (nrArguments > 2) {
      TString thisInputSigmaName = ((TObjString*)thisArgumentsArray->At(2))->GetString();
      std::string tmp = thisInputSigmaName.Data();
      DecomposeVariable(tmp, InputSigmaName, InputSigmaRange);
    }

    delete thisArgumentsArray;
    delete thisInputNameArray;
  } else {
    std::string tmp = InputName;
    DecomposeVariable(tmp, InputObservableName, InputObservableRange);
  }
}

// ____________________________________________________________________________|__________
// Decompose a factory like variable term into its name and range if specified
void RooFitUtils::CorrelationScheme::DecomposeVariable( std::string& InputName, std::string& InputVariableName, std::string& InputVariableRange )
{
  TString thisInputName = InputName.c_str();

  if (thisInputName.BeginsWith("[")) {
    assert(thisInputName.EndsWith("]"));
    InputVariableName = "";
    InputVariableRange = thisInputName;
  } else {
    TObjArray* thisInputNameArray = thisInputName.Tokenize("[");
    int nrTerms = thisInputNameArray->GetEntries();

    InputVariableName = ((TObjString*)thisInputNameArray->At(0))->GetString();

    if (nrTerms > 1) {
      TString thisInputVariableRange = ((TObjString*)thisInputNameArray->At(1))->GetString();
      InputVariableRange = "["+thisInputVariableRange;
      if (!thisInputVariableRange.EndsWith("]")) {
        InputVariableRange += "]";
      }
    } else {
      InputVariableRange = "";
    }

    delete thisInputNameArray;
  }
}








// ____________________________________________________________________________|__________
// Interface to add partly correlated parameters among 2 channels
void RooFitUtils::CorrelationScheme::CorrelateParameter( const char* OldParameterNamePlusMeasurement, const char* NewParameterName, double rho )
{
  // if parameters are fully correlated, use existing function
  if (rho == 1.0) {
    CorrelateParameter(OldParameterNamePlusMeasurement, NewParameterName);
  }

  TString allOldParameterNamePlusMeasurement = OldParameterNamePlusMeasurement;
  TObjArray* allOldParameterNamePlusMeasurementArray = allOldParameterNamePlusMeasurement.Tokenize(",");
  unsigned int numCorrPars = allOldParameterNamePlusMeasurementArray->GetEntries();

  // this interface takes correlates only 2 parameters
  if (numCorrPars != 2) {
    coutF(InputArguments) << "CorrelationScheme::CorrelateParameter(" << fName <<") " << numCorrPars << " specified. Can be 2 only." << endl;
    exit(-1);
  }

  TMatrixDSym cov(numCorrPars);
  for (unsigned int i = 0; i < numCorrPars; i++) {
    for (unsigned int j = 0; j < numCorrPars; j++) {
      if (i == j) cov(i,j) = 1.0;
      else        cov(i,j) = rho;
    }
  }

  CorrelateParameter(OldParameterNamePlusMeasurement, NewParameterName, cov);

  allOldParameterNamePlusMeasurementArray->Delete();
  delete allOldParameterNamePlusMeasurementArray;
}

// ____________________________________________________________________________|__________
// Interface to add partly correlated parameters among multiple channels
void RooFitUtils::CorrelationScheme::CorrelateParameter( const char* OldParameterNamePlusMeasurement, const char* NewParameterName, TMatrixDSym cov )
{
  // check dimensions
  if (cov.GetNcols() != cov.GetNrows()) {
    coutF(InputArguments) << "CorrelationScheme::CorrelateParameter(" << fName <<") expecting N x N matrix. Got " << cov.GetNcols() << " x " << cov.GetNrows() << endl;
    exit(-1);
  }

  // if parameters are fully correlated, use existing function
  bool fullycorrelated = true;
  for (int i = 0; i < cov.GetNcols(); i++) {
    for (int j = 0; j < cov.GetNrows(); j++) {
      cout << cov(i,j) << endl;
      if (cov(i,j) != 1.0) {
        fullycorrelated = false;
      }
    }
  }
  if (fullycorrelated) {
    coutW(InputArguments) << "CorrelationScheme::CorrelateParameter(" << fName <<") all parameters are fully correlated. Not using multivariate gaussian." << endl;
    CorrelateParameter(OldParameterNamePlusMeasurement, NewParameterName);
  }

  TString allOldParameterNamePlusMeasurement = OldParameterNamePlusMeasurement;
  TObjArray* allOldParameterNamePlusMeasurementArray = allOldParameterNamePlusMeasurement.Tokenize(",");
  int numCorrPars = allOldParameterNamePlusMeasurementArray->GetEntries();

  // match parameters to dimension of covariance matrix
  if (numCorrPars != cov.GetNcols()) {
    coutF(InputArguments) << "CorrelationScheme::CorrelateParameter(" << fName <<") " << numCorrPars << " specified but covariance matrix has " << cov.GetNcols() << " dimensions" << endl;
    exit(-1);
  }

  std::vector<TString> AllMeasurements;
  std::vector<TString> AllOldParameters;
  std::vector<TString> AllNewParameters;

  for (int itrPars = 0; itrPars < numCorrPars; ++itrPars) {
    TString thisOldParameterNamePlusMeasurement = ((TObjString*)allOldParameterNamePlusMeasurementArray->At(itrPars))->GetString();
    TObjArray* thisOldParameterNamePlusMeasurementArray = thisOldParameterNamePlusMeasurement.Tokenize("::");
    TString thisMeasurement = ((TObjString*)thisOldParameterNamePlusMeasurementArray->At(0))->GetString();
    TString thisOldParameter = ((TObjString*)thisOldParameterNamePlusMeasurementArray->At(1))->GetString();
    TString thisNewParamter = (string(NewParameterName)+"_"+string(thisMeasurement)).c_str();

    AllMeasurements.push_back(thisMeasurement);
    AllOldParameters.push_back(thisOldParameter);
    AllNewParameters.push_back(thisNewParamter);

    thisOldParameterNamePlusMeasurementArray->Delete();
    delete thisOldParameterNamePlusMeasurementArray;
  }

  for (unsigned int itrPars = 0; itrPars < AllOldParameters.size(); ++itrPars) {
    cout << AllMeasurements[itrPars] << endl;
    RenameParameter(AllMeasurements[itrPars], AllOldParameters[itrPars], AllNewParameters[itrPars]);
    IntroduceCorrelation(AllMeasurements[itrPars], AllNewParameters, cov);
  }

  allOldParameterNamePlusMeasurementArray->Delete();
  delete allOldParameterNamePlusMeasurementArray;
}

// ____________________________________________________________________________|__________
// Interface for specifying a parameter that should be renamed and correlated partially to others
void RooFitUtils::CorrelationScheme::IntroduceCorrelation( const char* MeasurementName, std::vector<TString> NewParameterNames, TMatrixDSym cov )
{
  unsigned int nDim = NewParameterNames.size();

  // get common new parameter name
  TString thisCommonName;
  TString thisMeasurementName = MeasurementName;
  for (unsigned int i = 0; i < nDim; ++i) {
    if (NewParameterNames[i].EndsWith(thisMeasurementName)) {
      thisCommonName = NewParameterNames[i];
      break;
    }
  }
  thisCommonName.Remove(thisCommonName.Sizeof()-1-thisMeasurementName.Sizeof());

  // build list of observables and means
  TString thisMean;
  TString thisObs;
  for (unsigned int i = 0; i < nDim-1; ++i) {
    thisObs += NewParameterNames[i] + "[0.0,-5.0,5.0],";
    thisMean += "nom_" + NewParameterNames[i] + "[0.0],";
  }
  thisObs += NewParameterNames[nDim-1] + "[0.0,-5.0,5.0]";
  thisMean += "nom_" + NewParameterNames[nDim-1] + "[0.0]";

  // Give a name to the correlation matrix and add it to the map with the pdf
  TString thisCorrName = thisCommonName + "_corr";
  TString OutputParameterName = TString::Format("MultiVarGaussian::%sCorr({%s},{%s},%s)", thisCommonName.Data(), thisObs.Data(), thisMean.Data(), thisCorrName.Data() ).Data();

  cout << OutputParameterName << endl;

  fCorrelationFactors[MeasurementName][thisCommonName.Data()].second.ResizeTo(cov);
  fCorrelationFactors[MeasurementName][thisCommonName.Data()].first = OutputParameterName;
  fCorrelationFactors[MeasurementName][thisCommonName.Data()].second = cov;
}
















// ____________________________________________________________________________|__________
// Print the correlation scheme as table
void CorrelationScheme::Print()
{
  std::set<std::string> allMeasurements;
  typedef std::map<std::string, RenamingMap>::iterator it_type;
  for (it_type corrItr = fCorrelationMap.begin(); corrItr != fCorrelationMap.end(); ++corrItr) {
    allMeasurements.insert(corrItr->first);
  }
  Print(allMeasurements);
}

// ____________________________________________________________________________|__________
// Print the correlation scheme as table for given numbers of measurements
void CorrelationScheme::Print( std::set<std::string> thisMeasurements )
{
  coutI(ObjectHandling) << "CorrelationScheme::Print(" << fName << ") printing correlations scheme" << endl;

  int nrMeasurements = thisMeasurements.size();
  if (nrMeasurements == 0) {
    coutW(ObjectHandling) << "CorrelationScheme::Print(" << fName << ") no measurement added to correlation map" << endl;
    return;
  }

  std::map<std::string, std::vector<std::string> > correlationMap;
  string* header = new string[nrMeasurements];

  typedef std::map<std::string, RenamingMap>::iterator it_type1;
  typedef std::map<std::string, std::string>::iterator it_type2;
  int iMeas = 0;

  for (it_type1 corrItr = fCorrelationMap.begin(); corrItr != fCorrelationMap.end(); ++corrItr) {
    if (thisMeasurements.find(corrItr->first) == thisMeasurements.end()) {
      continue;
    }

    header[iMeas] = corrItr->first;
    iMeas++;
    RenamingMap thisRenamingMap = corrItr->second;
    std::map< std::string, std::string > thisRenamingMapMap = thisRenamingMap.GetRenamingMap();
    for (it_type2 parItr = thisRenamingMapMap.begin(); parItr != thisRenamingMapMap.end(); ++parItr) {
      std::string thisNewObservableName = parItr->second;
      correlationMap[thisNewObservableName].push_back(corrItr->first);
    }
  }

  int nrNuis = correlationMap.size();
  string* firstCol = new string[nrNuis+1];
  firstCol[0] = "Parameter";
  std::string** matrix = new string*[nrNuis];
  for (int in = 0; in < nrNuis; in++) {
    matrix[in] = new std::string[nrMeasurements];
    for (int i = 0; i < nrMeasurements; i++) {
      matrix[in][i] = "";
    }
  }

  typedef std::map<std::string, std::vector<std::string> >::iterator it_type3;
  int irow = 1;
  for (it_type3 iterator = correlationMap.begin(); iterator != correlationMap.end(); ++iterator) {
    firstCol[irow] = iterator->first;
    int icol = 0;
    for (it_type1 corrItr = fCorrelationMap.begin(); corrItr != fCorrelationMap.end(); ++corrItr) {
      if (thisMeasurements.find(corrItr->first) == thisMeasurements.end()) {
        continue;
      }

      if (std::find(iterator->second.begin(), iterator->second.end(), corrItr->first) != iterator->second.end()) {
        matrix[irow-1][icol] = "x";
      }
      icol++;
    }
    irow++;
  }

  coutI(ObjectHandling) << "CorrelationScheme::Print(" << fName << ") the following correlation scheme is used" << endl;

  cout << "\\begin{tabular}{l";
  for (int i = 0; i < nrMeasurements; i++) cout << "|c";
  cout << "}\n";
  PrintTable(firstCol, matrix, NULL, header, nrNuis, nrMeasurements, 0, cout);
  cout << "\\end{tabular}\n";

  delete[] header;
  delete[] firstCol;

  for (int in = 0; in < nrNuis; in++) {
    delete[] matrix[in];
  }

  delete[] matrix;
}

// ____________________________________________________________________________|__________
// Helper function to print latex tables nicely
void CorrelationScheme::PrintTable(string* firstCol, std::string** matrix, std::string** matrixErr, string* header, int nrRows, int nrCols, int nSigFig, ostream& ost, string indent, string delim, string ending)
{
  int** lmatrix = new int*[nrRows];
  int** lmatrix_pm = new int*[nrRows];

  for (int i = 0; i < nrRows; i++) {
    lmatrix[i] = new int[nrCols];
    lmatrix_pm[i] = new int[nrCols];
  }

  int maxLfirst = 0;
  for (int i = 0; i < nrRows+1; i++)  {
    maxLfirst = int(max(double(maxLfirst), double(firstCol[i].size())));
  }

  string pm = " $\\pm$ ";

  int* maxLength = new int[nrCols];
  int* maxLength_pm = new int[nrCols];

  for (int j = 0; j < nrCols; j++) {
    if (header) maxLength[j] = header[j].size();
    else maxLength[j] = 0;
    maxLength_pm[j] = 0;
  }

  for (int i = 0; i < nrRows; i++) {
    for (int j = 0; j < nrCols;j++) {
      stringstream str;
      str << setprecision(nSigFig);
      str << matrix[i][j];

      lmatrix[i][j] = str.str().size();
      maxLength[j] = int(max(double(maxLength[j]), double(lmatrix[i][j])));

      stringstream str2;
      str2 << setprecision(nSigFig);
      if (/*j != 0 && */matrixErr) str2 << pm << matrixErr[i][j];

      lmatrix_pm[i][j] = str2.str().size();
      maxLength_pm[j] = int(max(double(maxLength_pm[j]), double(lmatrix_pm[i][j])));
    }
  }

  if (header) {
    ost << indent;
    ost << firstCol[0];
    for (int j = firstCol[0].size(); j < maxLfirst; j++) ost << " ";
    ost << " & ";
    for (int i = 0; i < nrCols; i++) {
      ost << header[i];

      for (int k = header[i].size(); k < maxLength[i] + maxLength_pm[i]; k++) {
        ost << " ";
      }
      if (i < nrCols-1) {
        if (delim != "") ost << delim;
      } else  {
        if (ending != "") ost << ending;
        ost << "\n" << indent << "\\hline\n";
      }
    }
  }

  for (int i = 0; i < nrRows; i++) {
    ost << indent;
    ost << firstCol[i+1];
    for (int j = firstCol[i+1].size(); j < maxLfirst; j++) ost << " ";
    ost << " & ";
    for (int j=0;j<nrCols;j++) {
      ost << setprecision(nSigFig);
      ost << matrix[i][j];
      for (int k = lmatrix[i][j]; k < maxLength[j]; k++) {
        ost << " ";
      }

      if (/*j != 0 && */matrixErr) {
        ost << pm << matrixErr[i][j];
        for (int k = lmatrix_pm[i][j]; k < maxLength_pm[j]; k++) {
          ost << " ";
        }
      }

      if (j < nrCols-1) {
        if (delim != "") ost << delim;
      } else {
        if (ending != "") ost << ending;
        ost << "\n";
      }
    }
  }

  delete[] maxLength;
  delete[] maxLength_pm;

  for (int i = 0; i < nrRows;i++) {
    delete[] lmatrix[i];
    delete[] lmatrix_pm[i];
  }

  delete[] lmatrix;
  delete[] lmatrix_pm;
}
