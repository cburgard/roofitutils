#ifndef ABSMEASUREMENT
#define ABSMEASUREMENT

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TNamed.h"
#include "TFile.h"
#include "TList.h"
#include "TIterator.h"
#include "Math/MinimizerOptions.h"

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooAbsArg.h"
#include "RooLinkedListIter.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooProdPdf.h"
#include "RooArgList.h"
#include "RooMsgService.h"
#include "RooRealSumPdf.h"

#include "RooStats/ModelConfig.h"

#include "RooFitUtils/WildcardList.h"

struct TOwnedList : public TList {
  // A collection class for keeping TObjects for deletion.
  // TOwnedList is like TList with SetOwner(), but really deletes all objects, whether or not on heap.
  // This is a horrible hack to work round the fact that RooArgSet and RooDataSet objects have have IsOnHeap() false.
  TOwnedList();
  virtual ~TOwnedList();
  virtual void Clear (Option_t* option="") override;
  ClassDefOverride(TOwnedList,0)
};


class AbsMeasurement : public TNamed {

// ____________________________________________________________________________|__________
public:

  // Constructor and destructor
  AbsMeasurement( const std::string& MeasurementName, const std::string& FileName, const std::string& WorkspaceName, const std::string& ModelConfigName, const std::string& DataName );
  AbsMeasurement( const std::string& MeasurementName, const std::string& WorkspaceName, const std::string& ModelConfigName, const std::string& DataName );
  AbsMeasurement( const std::string& MeasurementName, RooWorkspace* ws, const std::string& ModelConfigName, const std::string& DataName );
  virtual ~AbsMeasurement();

  // Accessors
  void SetFileName( const std::string& FileName ) { fFileName = FileName; }
  std::string GetFileName() { return fFileName; }

  void SetWorkspaceName( const std::string& WorkspaceName ) { fWorkspaceName = WorkspaceName; }
  std::string GetWorkspaceName() { return fWorkspaceName; }

  void SetModelConfigName( const std::string& ModelConfigName ) { fModelConfigName = ModelConfigName; }
  std::string GetModelConfigName() { return fModelConfigName; }

  void SetDataName( const std::string& DataName ) { fDataName = DataName; }
  std::string GetDataName() { return fDataName; }

  void SetNamePdf( std::string name ) { fPdfName = name; }
  std::string GetNamePdf() { return fPdfName; }

  void SetNameCategory( std::string name ) { fCategoryName = name; }
  std::string GetNameCategory() { return fCategoryName; }

  void SetPdf( RooAbsPdf* Pdf ) { fPdf = Pdf; }
  RooAbsPdf* GetPdf() { return fPdf; }

  void SetData( RooAbsData* Data ) { fData = Data; }
  RooAbsData* GetData() { return fData; }

  void SetNuisanceParameters( RooArgSet* NuisanceParameters ) { fNuisanceParameters = NuisanceParameters; }
  RooArgSet* GetNuisanceParameters() { return fNuisanceParameters; }

  void SetObservables( RooArgSet* Observables ) { fObservables = Observables; }
  RooArgSet* GetObservables() { return fObservables; }

  void SetGlobalObservables( RooArgSet* GlobalObservables ) { fGlobalObservables = GlobalObservables; }
  RooArgSet* GetGlobalObservables() { return fGlobalObservables; }

  void SetMinimizerType( std::string name ) { fMinimizerType = name; ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str()); }
  std::string GetMinimizerType() { return fMinimizerType; }

  void SetMinimizerAlgorithm( std::string name ) { fMinimizerAlgo = name; ROOT::Math::MinimizerOptions::SetDefaultMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), fMinimizerAlgo.c_str()); }
  std::string GetMinimizerAlgorithm() { return fMinimizerAlgo; }

  void SetMinimizerStrategy( int number ) { fDefaultStrategy = number; ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy); }
  int GetMinimizerStrategy() { return fDefaultStrategy; }

  void SetMinimizerNumCPU( int number ) { fNumCPU = number; }
  int GetMinimizerNumCPU() { return fNumCPU; }

  void enablePruning( const std::string& poi = "mu_ggF,mu_VBF,mu_WH,mu_ZH,mu_ttH", const std::string& filter = ".*", const std::string& weight = "19.12,1.573,0.6951,0.4102,0.1277", const std::string& threshold = "auto:3", int additionalDigit = 1 );
  void disablePruning() { fIsPrunable = kFALSE; }
  bool isPrunable() { return fIsPrunable; }

  RooWorkspace* GetWorkspace(){
    return fWorkSpace;
  }

  void SetPrunedNuisanceParameters( std::string parameters );
  void SetPrunedNuisanceParameters( std::list< std::string > parameters ) { fPrunedNuisanceParameters = parameters; }
  std::list< std::string > GetPrunedNuisanceParameters() { return fPrunedNuisanceParameters; }

  void SetDatasetBinning ( Int_t setNbins = 500, const char* generateBinnedTag = 0, const char* binnedCategories = 0, const char* unbinnedCategories = 0, const char* weightVarName = "weightVar" );

  // Steering
  virtual void initialise() = 0;
  void writeToFile();
  void writeToFile( const char *fileName );
  void PruneNuisanceParameters();

  static bool AlmostEqualUlpsAndAbs( float A, float B, float maxDiff, int maxUlpsDiff );

// ____________________________________________________________________________|__________
protected:

  std::list< std::string > PruneNuisanceParameters(const TMatrixDSym chesse, RooFitResult* fitresult, const std::string& poi = "mu_ggF,mu_VBF,mu_WH,mu_ZH,mu_ttH", const std::string& filter = ".*", const std::string& weight = "19.12,1.573,0.6951,0.4102,0.1277", const std::string& threshold = "auto:3", int additionalDigit = 1, std::list< std::string > prePrunedParameters = std::list< std::string >());
  void RemoveParameter( TMatrixDSym& hes, RooArgList& pars, std::list< std::string > names );
  void PrintRanking( std::set< std::pair< double, std::string > > uncerts, double initTotalError );
  std::pair< double, double > PDGrounding( double value, double error, int additionalDigit = 1 );
  int GetThreeDigits( double error );
  int GetNSigDigits( int threeDigits );
  double frexp10( double x, int* exp );
  double FormatValue( double value, int exponent, int nDigits, int extraRound = 0 );
  RooDataSet* SetDatasetBinning ( const RooAbsPdf* pdf, const RooAbsData* data, Int_t setNbins = 500, const char* generateBinnedTag = 0, const char* binnedCategories = 0, const char* unbinnedCategories = 0, const char* weightVarName = "weightVar" );

// ____________________________________________________________________________|__________
protected:

  friend class Measurement;
  friend class CombinedMeasurement;

// ____________________________________________________________________________|__________
private:

  bool fIsInitialised=false;
  bool fIsPrunable   =false;
  bool fSetBinning   =false;

  std::string fFileName;
  std::string fWorkspaceName;
  std::string fModelConfigName;
  std::string fDataName;
  std::string fPdfName;
  std::string fCategoryName;

  std::string fPruningPoi;
  std::string fPruningFilter;
  std::string fPruningWeight;
  std::string fPruningThreshold;
  int fPruningAdditionalDigit;
  std::list< std::string > fPrunedNuisanceParameters;

  Int_t fSetNbins;
  const char* fGenerateBinnedTag;
  const char* fBinnedCategories;
  const char* fUnbinnedCategories;
  const char* fWeightVarName;

  TFile* fFile = NULL;
  RooWorkspace* fWorkSpace = NULL;
  RooStats::ModelConfig* fModelConfig = NULL;
  RooAbsPdf* fPdf = NULL;
  RooAbsData* fData = NULL;
  RooAbsData* fAsimovData = NULL;
  RooArgSet* fParametersOfInterest = NULL;
  RooArgSet* fNuisanceParameters = NULL;
  RooArgSet* fObservables = NULL;
  RooArgSet* fGlobalObservables = NULL;

  std::string fMinimizerType = "Minuit2";
  std::string fMinimizerAlgo = "Migrad";
  int fDefaultStrategy = 0;
  int fNumCPU = 1;

// ____________________________________________________________________________|__________
protected:

  ClassDefOverride(AbsMeasurement, 1)

};

#endif
