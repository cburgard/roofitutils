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
      Scan(const std::vector<std::string> &parnames);
      void add(const std::vector<double> &parvals, int fitStatus, double nllval);
      void printTable();
      std::vector<std::string> parNames;
      std::vector<std::vector<double>> parValues;
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
      bool ok();
      static bool ok(int status);
      int status;
      int strategy;
      double nll;
      int ndim;
      ROOT::Fit::FitConfig config;
    };
    std::vector<Scan> scans;
    std::vector<Parameter> parameters;
    Eigen *eigen;
    Minimization min;
    RooFitResult *fit;
    TMatrixDSym *hesse;
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
  int minimize(const RooCmdArg &arg1, const RooCmdArg &arg2 = RooCmdArg::none(),
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
  int minimize(const RooLinkedList &cmdList);
  int minimize();

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
  static RooCmdArg Eigen(Bool_t flag = kTRUE) {
    return RooCmdArg("Eigen", flag, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg NumRetryFit(Int_t retry) {
    return RooCmdArg("NumRetryFit", retry, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg MaxCalls(Int_t maxcalls) {
    return RooCmdArg("MaxCalls", maxcalls, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg Eps(double eps) {
    return RooCmdArg("Eps", 0, 0, eps, 0, 0, 0, 0, 0);
  }
  static RooCmdArg NSigma(double nsigma) {
    return RooCmdArg("NSigma", 0, 0, nsigma, 0, 0, 0, 0, 0);
  }
  static RooCmdArg Scan(Bool_t flag = kTRUE) {
    return RooCmdArg("Scan", flag, 0, 0, 0, 0, 0, 0, 0);
  }
  static RooCmdArg Scan(const RooArgSet &scanArgs) {
    return RooCmdArg("Scan", kTRUE, 0, 0, 0, 0, 0, &scanArgs, 0);
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

  template <class A> int parseFitConfig(const A &cmdList);
  template <class A> int parseNllConfig(const A &cmdList);

  // ____________________________________________________________________________|__________
protected:
  void initialize();
  Result *run();
  void setup();
  Result::Eigen *eigenAnalysis(const TMatrixDSym &hesse);
  void findSigma(Result *result, const RooAbsCollection &pois);
  virtual double findSigma(Result *result, const double guessval,
                           const double val_mle, RooRealVar *par,
                           const double nsigma = +1.0, const int maxiter = 25);
  void scan(Result *r,
            const std::map<const std::string, std::vector<double>> &params);
  void scan(Result *r, const std::vector<std::string> &parnames,
            const std::vector<std::vector<double>> &points);
  virtual Result::Minimization robustMinimize();

  // ____________________________________________________________________________|__________
protected:
  Result *fResult = NULL;

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
  Int_t fMinos;
  Int_t fScan;
  Int_t fNumee;
  Int_t fDoEEWall;
  Int_t fRetry;
  Int_t fEigen;
  Int_t fReuseMinimizer;
  Int_t fReuseNLL;
  Int_t fMaxIterations;
  Double_t fEps;
  Double_t fNsigma;
  Double_t fPrecision;
  const RooArgSet *fMinosSet = NULL;
  const RooArgSet *fCondSet  = NULL;
  const RooArgSet *fScanSet  = NULL;
  std::string fMinimizerType;
  std::string fMinimizerAlgo;

  // ____________________________________________________________________________|__________
protected:
  ClassDefOverride(ExtendedMinimizer, 0)
};
}

#endif
