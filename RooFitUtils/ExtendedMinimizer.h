// this file looks like plain C, but it's actually -*- c++ -*-
#ifndef EXTENDEDMINIMIZER
#define EXTENDEDMINIMIZER

#include <string>
#include <vector>

#include "TNamed.h"

#include "RooMinimizer.h"
class TGraph;
#include "TMatrixD.h"
#include "TVectorD.h"

#include <map>

namespace RooFitUtils {
  class ExtendedModel;
  class ExtendedMinimizer : public TNamed {
  public:
  typedef std::pair<TGraph *, TGraph *> GraphPair;
  typedef std::map<double, double> ValueMap;
  class Result {
  public:
    Result();
    class Scan {
    public:
      Scan(const std::string& name, const std::vector<std::string> &parnames);
      Scan(const std::string& name, const std::vector<std::string> &parnames, const std::vector<std::string> &extranames);      
      void add(const std::vector<double> &parvals, int fitStatus, double nllval);
      void add(const std::vector<double> &parvals, int fitStatus, double nllval, const std::vector<double>& extravals);      
      void printTable();
      std::string name;
      std::vector<std::string> parNames;
      std::vector<std::vector<double>> parValues;
      std::vector<std::string> extraParNames;      
      std::vector<std::vector<double>> extraParValues;      
      std::vector<int> fitStatus;
      std::vector<double> nllValues;
    };
    class Parameter {
    public:
      Parameter(const std::string &n, double v, double eH, double eL);
      std::string name;
      double value;
      double errHi;
      double errLo;
    };
    class Eigen {
    public:
      Eigen(const TVectorD &values, const TMatrixD &vectors);
      TVectorD values;
      TMatrixD vectors;
    };
    class Minimization {
    public:
      Minimization();
      bool ok() const;
      static bool ok(int status);
      int status = -1;
      int strategy = -1;
      double nll = 0;
      double nllOffset = 0;
      int constOpt = 0;            
      int ndim = 0;
      ROOT::Fit::FitConfig config;
      TMatrixDSym *hesse = NULL;
      TMatrixDSym *cov = NULL;      
      int covqual = -1;
      std::vector<Parameter> parameters;
      RooFitResult *fit = NULL;
    };
    std::vector<Scan> scans;
    Eigen *eigen;
    Minimization min;
    operator bool() const { return this->min.ok(); }
  };

  // Constructor and destructor
  ExtendedMinimizer(const char *minimizerName, RooFitUtils::ExtendedModel* model, const RooLinkedList &args);
  ExtendedMinimizer(const char *minimizerName, RooFitUtils::ExtendedModel* model);
  ExtendedMinimizer(const char *minimizerName, RooAbsPdf *pdf, RooAbsData *data,
                    RooWorkspace *workspace, const RooLinkedList &args);
  ExtendedMinimizer(const char *minimizerName, RooAbsPdf *pdf = NULL,
                    RooAbsData *data = NULL, RooWorkspace *workspace = NULL);
  ExtendedMinimizer() = default;
  virtual ~ExtendedMinimizer();

  void SetPdf(RooAbsPdf *Pdf) { fPdf = Pdf; }
  RooAbsPdf *GetPdf() { return fPdf; }

  void SetData(RooAbsData *Data) { fData = Data; }
  RooAbsData *GetData() { return fData; }

  void SetNllDirty();
  RooAbsReal *GetNll() {
    if (!fNll)
      this->setup();
    return fNll;
  }

  Result *getResult(bool make = false);
  RooFitResult *GetFitResult();
  TMatrixDSym GetHesseMatrix();
  double GetMinNll();

  inline void setPrintLevel(int plevel) { this->fPrintLevel = plevel; }
  inline int getPrintLevel() { return this->fPrintLevel; }

  // Steering
  Result minimize(const RooCmdArg &arg1, const RooCmdArg &arg2 = RooCmdArg::none(),
               const RooCmdArg &arg3 = RooCmdArg::none(),
               const RooCmdArg &arg4 = RooCmdArg::none(),
               const RooCmdArg &arg5 = RooCmdArg::none(),
               const RooCmdArg &arg6 = RooCmdArg::none(),
               const RooCmdArg &arg7 = RooCmdArg::none(),
               const RooCmdArg &arg8 = RooCmdArg::none(),
               const RooCmdArg &arg9 = RooCmdArg::none(),
               const RooCmdArg &arg10 = RooCmdArg::none(),
               const RooCmdArg &arg11 = RooCmdArg::none(),
               const RooCmdArg &arg12 = RooCmdArg::none());
  Result minimize(const RooLinkedList &cmdList);
  Result minimize();

  void findSigma(const std::vector<std::string> &pars);
  void findSigma();

  void scan(const std::map<const std::string, std::vector<double>> &params);
  void scan(const std::vector<std::string> &parnames,
            const std::vector<std::vector<double>> &points);

  ExtendedMinimizer::GraphPair
  prepareProfile(ExtendedMinimizer::ValueMap map_poi2nll);
  ExtendedMinimizer::GraphPair createProfile(RooRealVar *var, double lo,
                                             double hi, int nbins);
  ExtendedMinimizer::ValueMap createProfileValues(RooRealVar *var, double lo,
                                                  double hi, int nbins);
  // ____________________________________________________________________________|__________
public:
  static RooCmdArg Minimize(Bool_t flag = kTRUE) {
    return RooCmdArg("Minimize", flag, 0, 0, 0, 0, 0, 0, 0);
  }  
  static RooCmdArg Eigen(Bool_t flag = kTRUE) {
    return RooCmdArg("Eigen", flag, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg NumRetryFit(Int_t retry) {
    return RooCmdArg("NumRetryFit", retry, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg Verbose(Int_t verbose) {
    return RooCmdArg("verbose", verbose, 0, 0, 0, 0, 0, 0, 0);
  }  
  static RooCmdArg MaxCalls(Int_t maxcalls) {
    return RooCmdArg("MaxCalls", maxcalls, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg MaxIterations(Int_t maxiterations) {
    return RooCmdArg("MaxIterations", maxiterations, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg Eps(double eps) {
    return RooCmdArg("Eps", 0, 0, eps, 0, 0, 0, 0, 0);
  }
  static RooCmdArg NSigma(double nsigma) {
    return RooCmdArg("NSigma", 0, 0, nsigma, 0, 0, 0, 0, 0);
  }
  static RooCmdArg FindSigma(Bool_t flag = kTRUE, double nSigma=1, int maxIter = 25, double precision=0.005) {
    return RooCmdArg("FindSigma", flag, maxIter, nSigma, precision, 0, 0, 0, 0);
  }
  static RooCmdArg FindSigma(const RooArgSet &scanArgs, double nSigma=1, int maxIter = 25, double precision=0.005) {
    return RooCmdArg("FindSigma", kTRUE, maxIter, nSigma, precision, 0, 0, &scanArgs, 0);
  }
  static RooCmdArg Cond(const RooArgSet &condArgs) {
    return RooCmdArg("Cond", kTRUE, 0, 0, 0, 0, 0, &condArgs, 0);
  }
  static RooCmdArg ReuseMinimizer(Bool_t flag = kFALSE) {
    return RooCmdArg("ReuseMinimizer", flag, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg ReuseNLL(Bool_t flag = kTRUE) {
    return RooCmdArg("ReuseNLL", flag, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg UseChi2(Bool_t flag = kTRUE) {
    return RooCmdArg("UseChi2", flag, 0, 0, 0, 0, 0, 0, 0);
  }
  template <class A> int parseFitConfig(const A &cmdList);
  template <class A> int parseNllConfig(const A &cmdList);

  // ____________________________________________________________________________|__________
protected:
  void writeReproducer(const std::string& name);
  void initialize();
  Result *run();
  int runHesse(Result::Minimization& mini);  
  void setup();
  Result::Eigen *eigenAnalysis(const TMatrixDSym &hesse);
  void findSigma(Result* r);  
  void findSigma(Result *result, const RooAbsCollection &pois);
  double findSigma(Result *result, const Result::Minimization& min,
                   const double guessval, const double val_mle, RooRealVar *par, const double nsigma = +1.0);
  void scan(Result *r,
            const std::map<const std::string, std::vector<double>> &params);
  void scan(Result *r, const std::vector<std::string> &parnames,
            const std::vector<std::vector<double>> &points);
  Result::Minimization robustMinimize();

  // ____________________________________________________________________________|__________
protected:
  Result *fResult = NULL;

  ExtendedModel *fModel = NULL;
  RooWorkspace *fWorkspace = NULL;
  RooAbsPdf *fPdf = NULL;
  RooAbsData *fData = NULL;
  RooAbsReal *fNll = NULL;
  RooMinimizer *fMinimizer = NULL;

  std::vector<RooCmdArg> fOwnedArgs;
  RooLinkedList fNllCmdList;
  RooLinkedList fFitCmdList;

  Int_t fOffset;
  Int_t fOptConst;
  Int_t fVerbose;
  Int_t fSave;
  Int_t fTimer;
  Int_t fPrintLevel;
  Int_t fDefaultStrategy;
  Int_t fHesse;
  Int_t fMinimize;
  Int_t fMinos;
  Int_t fNumee;
  Int_t fDoEEWall;
  Int_t fRetry;
  Int_t fEigen;
  Int_t fReuseMinimizer;
  Int_t fReuseNLL;
  Int_t fChi2;
  Int_t fMaxCalls;
  Int_t fMaxIterations;  
  Double_t fEps;
  const RooArgSet *fPenaltyMini  = NULL;
  const RooArgSet *fPoiSet  = NULL;
  const RooArgSet *fMinosSet = NULL;
  const RooArgSet *fCondSet  = NULL;
  std::string fMinimizerType;
  std::string fMinimizerAlgo;

  Int_t fFindSigma = 0;
  Double_t fFindSigmaN = 1;
  Double_t fFindSigmaIter = 25;  
  Double_t fFindSigmaPrecision = 0.005;
  RooArgSet *fFindSigmaSet  = NULL;

  
  // ____________________________________________________________________________|__________
protected:
  ClassDefOverride(ExtendedMinimizer, 0)
};
}

#endif
