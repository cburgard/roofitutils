#include "RooFitUtils/ExtendedMinimizer.h"

#include "TFile.h"
#include "TMath.h"
#include "TMatrixDSymEigen.h"
#include <math.h>
// unfortunately, this is needed
#define private public
#define protected public
#include "RooFitResult.h"
#undef private
#undef protected

#include "RooCmdConfig.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"
#include "RooRealVar.h"

#include "RooFitUtils/Utils.h"

#include "RooStats/RooStatsUtils.h"

#include <limits>

#define nan std::numeric_limits<double>::quiet_NaN()
#define inf std::numeric_limits<double>::infinity()

ClassImp(RooFitUtils::ExtendedMinimizer)

    // ____________________________________________________________________________|__________

    RooFitUtils::ExtendedMinimizer::Result::Result()
    : eigen(NULL), fit(NULL), hesse(NULL) {
  // nothing here
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Minimization::Minimization()
    : status(-1), strategy(-1), nll(nan) {
  // nothing here
}

// ____________________________________________________________________________|__________

bool RooFitUtils::ExtendedMinimizer::Result::Minimization::ok(int status) {
  // return true if status is 0 or 1, false otherwise
  return (status == 0 || status == 1);
}

// ____________________________________________________________________________|__________

bool RooFitUtils::ExtendedMinimizer::Result::Minimization::ok() {
  // return status of minization
  return ExtendedMinimizer::Result::Minimization::ok(this->status);
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Parameter::Parameter(
    const std::string &n, double v, double eH, double eL)
    : name(n), value(v), errHi(eH), errLo(eL) {
  // nothing here
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Eigen::Eigen(const TVectorD &vals,
                                                     const TMatrixD &vecs)
    : values(vals), vectors(vecs) {
  // nothing here
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Scan::Scan(
    const std::vector<std::string> &parnames)
    : parNames(parnames) {
  // nothing here
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::Result::Scan::add(
    const std::vector<double> &parvals, double nllval) {
  // add a new entry to a scan
  if (parvals.size() != this->parNames.size()) {
    throw std::runtime_error("cannot add parameter list with wrong length!");
  }
  this->parValues.push_back(parvals);
  this->nllValues.push_back(nllval);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::Result::Scan::printTable() {
  // print the scan values as a table
  auto prec = std::cout.precision();
  std::cout.precision(15);
  for (size_t i = 0; i < this->parNames.size(); ++i) {
    std::cout << parNames[i] << "\t";
  }
  std::cout << "nll" << std::endl;
  for (size_t j = 0; j < this->nllValues.size(); ++j) {
    for (size_t i = 0; i < this->parNames.size(); ++i) {
      std::cout << parValues[j][i] << "\t";
    }
    std::cout << nllValues[j] << std::endl;
  }
  std::cout.precision(prec);
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result *
RooFitUtils::ExtendedMinimizer::getResult(bool make) {
  // obtain the result from the minimizer
  if (!this->fResult && make) {
    this->fResult = new ExtendedMinimizer::Result();
  }
  return this->fResult;
}

// ____________________________________________________________________________|__________

RooFitResult *RooFitUtils::ExtendedMinimizer::GetFitResult() {
  // obtain the fit result from the minimizer
  if (this->fResult)
    return this->fResult->fit;
  return NULL;
}

// ____________________________________________________________________________|__________

TMatrixDSym RooFitUtils::ExtendedMinimizer::GetHesseMatrix() {
  // obtain the hesse matrix from the minimizer
  if (this->fResult && this->fResult->hesse)
    return *(this->fResult->hesse);
  return TMatrixDSym();
}

// ____________________________________________________________________________|__________

double RooFitUtils::ExtendedMinimizer::GetMinNll() {
  // obtain the Nll minimum from the minimizer
  if (this->fResult)
    return this->fResult->min.nll;
  return nan;
}

// ____________________________________________________________________________|__________

namespace {

// ____________________________________________________________________________|__________

RooLinkedList *makeList(const RooCmdArg &arg1 = RooCmdArg::none(),
                        const RooCmdArg &arg2 = RooCmdArg::none(),
                        const RooCmdArg &arg3 = RooCmdArg::none(),
                        const RooCmdArg &arg4 = RooCmdArg::none(),
                        const RooCmdArg &arg5 = RooCmdArg::none(),
                        const RooCmdArg &arg6 = RooCmdArg::none(),
                        const RooCmdArg &arg7 = RooCmdArg::none(),
                        const RooCmdArg &arg8 = RooCmdArg::none(),
                        const RooCmdArg &arg9 = RooCmdArg::none(),
                        const RooCmdArg &arg10 = RooCmdArg::none(),
                        const RooCmdArg &arg11 = RooCmdArg::none(),
                        const RooCmdArg &arg12 = RooCmdArg::none(),
                        const RooCmdArg &arg13 = RooCmdArg::none(),
                        const RooCmdArg &arg14 = RooCmdArg::none(),
                        const RooCmdArg &arg15 = RooCmdArg::none(),
                        const RooCmdArg &arg16 = RooCmdArg::none(),
                        const RooCmdArg &arg17 = RooCmdArg::none(),
                        const RooCmdArg &arg18 = RooCmdArg::none(),
                        const RooCmdArg &arg19 = RooCmdArg::none(),
                        const RooCmdArg &arg20 = RooCmdArg::none()) {
  // helper function
  RooLinkedList *l = new RooLinkedList();
  l->SetName("CmdList");
  l->Add(arg1.Clone());
  l->Add(arg2.Clone());
  l->Add(arg3.Clone());
  l->Add(arg4.Clone());
  l->Add(arg5.Clone());
  l->Add(arg6.Clone());
  l->Add(arg7.Clone());
  l->Add(arg8.Clone());
  l->Add(arg9.Clone());
  l->Add(arg10.Clone());
  l->Add(arg11.Clone());
  l->Add(arg12.Clone());
  l->Add(arg13.Clone());
  l->Add(arg14.Clone());
  l->Add(arg15.Clone());
  l->Add(arg16.Clone());
  l->Add(arg17.Clone());
  l->Add(arg18.Clone());
  l->Add(arg19.Clone());
  l->Add(arg20.Clone());
  return l;
}

inline const double *getAry(const std::vector<double> &numbers) {
  return &numbers[0];
}

inline void
addAllPoints(std::vector<std::vector<double>> &points,
             const std::vector<std::string> &parnames,
             const std::map<const std::string, std::vector<double>> &params,
             std::vector<double> &currentvals, size_t idx) {
  if (idx < parnames.size()) {
    for (auto val : params.at(parnames[idx])) {
      currentvals[idx] = val;
      addAllPoints(points, parnames, params, currentvals, idx + 1);
    }
  } else {
    points.push_back(currentvals);
  }
}

inline int addAllArgs(const RooLinkedList &orig, RooLinkedList &target) {
  int n = 0;
  RooLinkedListIter it = orig.iterator();
  while (true) {
    TObject *obj = it.Next();
    if (!obj)
      break;
    RooCmdArg *v = dynamic_cast<RooCmdArg *>(obj);
    if (!v)
      break;
    if (v != &RooCmdArg::none()) {
      target.Add(v);
      n++;
    }
  }
  return n;
}

void inverseFilterCmdList(RooLinkedList &cmdInList, const char *cmdNameList) {
  RooLinkedList filterList;
  if (!cmdNameList) {
    cmdInList.Clear();
    return;
  }

  // Copy command list for parsing
  char buf[1024];
  strlcpy(buf, cmdNameList, 1024);

  char *name = strtok(buf, ",");
  while (name) {
    TObject *cmd = cmdInList.FindObject(name);
    if (cmd) {
      filterList.Add(cmd);
    }
    name = strtok(0, ",");
  }
  cmdInList.Clear();
  addAllArgs(filterList, cmdInList);
}

void setVals(const RooAbsCollection &vars, const RooAbsCollection *snap,
             bool setConstant = false) {
  if (!snap)
    return;
  RooAbsArg *obj;
  RooFIter itr(vars.fwdIterator());
  while ((obj = itr.next())) {
    RooRealVar *v = dynamic_cast<RooRealVar *>(obj);
    if (!v)
      continue;
    RooAbsReal *sv = dynamic_cast<RooAbsReal *>(snap->find(v->GetName()));
    if (!sv)
      continue;
    v->setVal(sv->getVal());
    if (setConstant) {
      v->setConstant(true);
    }
  }
}

Double_t useLimits(const RooRealVar *par, Double_t val) {
  if (val < par->getMin()) {
    return par->getMin();
  } else if (val > par->getMax()) {
    return par->getMax();
  }
  return val;
}
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::ExtendedMinimizer(const char *minimizerName,
                                                  RooAbsPdf *pdf,
                                                  RooAbsData *data,
                                                  RooWorkspace *workspace,
                                                  const RooLinkedList &argList)
    : ExtendedMinimizer(minimizerName, pdf, data, workspace) {
  // Constructor
  parseConfig(argList);
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::ExtendedMinimizer(const char *minimizerName,
                                                  RooAbsPdf *pdf,
                                                  RooAbsData *data,
                                                  RooWorkspace *workspace)
    : // Constructor
      TNamed(minimizerName, minimizerName),
      fWorkspace(workspace), fPdf(pdf), fData(data),

      fOffset(0), fOptConst(2), fVerbose(0), fSave(0), fTimer(1),
      fPrintLevel(1), fDefaultStrategy(0), fHesse(0), fMinos(0), fScan(0),
      fNumee(5), fDoEEWall(1), fRetry(0), fEigen(0), fReuseMinimizer(0),
      fReuseNLL(0), fEps(1.0), fNsigma(1), // 1sigma 1dof
      fPrecision(0.005),

      fMinimizerType("Minuit2"), fMinimizerAlgo("Migrad") {
  // Constructor
  fMinosSet = NULL;
  fCondSet = NULL;
  fScanSet = NULL;

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(),
                                                    fMinimizerAlgo.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(fPrintLevel);

  fNllCmdList.SetName("NllCmdList");
  fFitCmdList.SetName("FitCmdList");

  coutP(InputArguments) << "ExtendedMinimizer::ExtendedMinimizer(" << fName
                        << ") created" << std::endl;
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::~ExtendedMinimizer() {
  // Destructor
}

// ____________________________________________________________________________|__________

int RooFitUtils::ExtendedMinimizer::minimize(
    const RooCmdArg &arg1, const RooCmdArg &arg2, const RooCmdArg &arg3,
    const RooCmdArg &arg4, const RooCmdArg &arg5, const RooCmdArg &arg6,
    const RooCmdArg &arg7, const RooCmdArg &arg8, const RooCmdArg &arg9,
    const RooCmdArg &arg10, const RooCmdArg &arg11, const RooCmdArg &arg12) {
  // Minimize function with iterative retry strategy adopted, simplified and
  // extended from RooAbsPdf::fitTo()
  RooLinkedList *l = makeList(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                              arg9, arg10, arg11, arg12);
  int status = minimize(*l);
  delete l;
  return status;
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::setup() {
  // initialization
  if (!fPdf) {
    throw std::runtime_error("[Error] ExtendedMinimizer::setup: No Pdf set!");
  }

  if (!fData) {
    throw std::runtime_error("[Error] ExtendedMinimizer::setup: No Data set!");
  }

  if (!fReuseNLL || !fNll) {
    if (fReuseNLL) {
      coutW(InputArguments)
          << "ExtendedMinimizer::createMinimizer: NLL is supposed to be "
             "re-used but not present, creating new NLL!"
          << std::endl;
    }
    coutI(InputArguments) << "Creating new Nll" << std::endl;
    if (fNll) {
      coutW(InputArguments)
          << "ExtendedMinimizer::createMinimizer: deleting previous NLL!"
          << std::endl;
      delete fNll;
    }

    if (fWorkspace) {
      if (RooFitUtils::RooStarMomentMorphFix) {
        int n = fixRooStarMomentMorph(fWorkspace);
        if (n > 0)
          coutP(InputArguments)
              << "Fixed cache of " << n << " instances of RooStarMomentMorph"
              << std::endl;
      }
    }

    fNll = fPdf->createNLL(*fData, fNllCmdList);

  } else {
    coutI(InputArguments) << "Using existing Nll" << std::endl;
  }

  if (fWorkspace) {
    if (RooFitUtils::RooStarMomentMorphFix) {
      int n = fixRooStarMomentMorph(fWorkspace);
      if (n > 0)
        coutP(InputArguments)
            << "Re-Fixed cache of " << n << " instances of RooStarMomentMorph"
            << std::endl;
    }
  }

  if (!fNll) {
    throw std::runtime_error(
        "[Error] ExtendedMinimizer::setup: Failed to obtain NLL");
  }

  if (!fMinimizer || !fReuseMinimizer) {
    if (fMinimizer)
      delete fMinimizer;

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(),
                                                      fMinimizerAlgo.c_str());
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(fPrintLevel);

    fMinimizer = new RooMinimizer(*fNll);
  }
}

// ____________________________________________________________________________|__________

template <class A>
int RooFitUtils::ExtendedMinimizer::parseConfig(const A &cmdList) {
  // parse the configuration
  RooCmdConfig pc(Form("ExtendedMinimizer::parseConfig(%s)", GetName()));
  pc.allowUndefined();

  fFitCmdList.Clear();
  addAllArgs(cmdList, fFitCmdList);
  fNllCmdList.Clear();
  addAllArgs(cmdList, fNllCmdList);
  inverseFilterCmdList(fNllCmdList, "NumCPU,Constrained,Constrain,CloneData,"
                                    "GlobalObservables,GlobalObservablesTag,"
                                    "OffsetLikelihood");

  pc.defineInt("optConst", "Optimize", 0, fOptConst);
  pc.defineInt("doOffset", "OffsetLikelihood", 0, fOffset);
  pc.defineInt("verbose", "Verbose", 0, fVerbose);
  pc.defineInt("doSave", "Save", 0, fSave);
  pc.defineInt("doTimer", "Timer", 0, fTimer);
  pc.defineInt("plevel", "PrintLevel", 0, fPrintLevel);
  pc.defineInt("strat", "Strategy", 0, fDefaultStrategy);
  pc.defineInt("hesse", "Hesse", 0, fHesse);
  pc.defineInt("minos", "Minos", 0, fMinos);
  pc.defineInt("scan", "Scan", 0, fScan);
  pc.defineInt("numee", "PrintEvalErrors", 0, fNumee);
  pc.defineInt("doEEWall", "EvalErrorWall", 0, fDoEEWall);
  pc.defineInt("retry", "NumRetryFit", 0, fRetry);
  pc.defineInt("eigen", "Eigen", 0, fEigen);
  pc.defineInt("reminim", "ReuseMinimizer", 0, fReuseMinimizer);
  pc.defineInt("renll", "ReuseNLL", 0, fReuseNLL);
  pc.defineDouble("eps", "Eps", 0, fEps);
  pc.defineDouble("nsigma", "NSigma", 0, fNsigma);
  pc.defineDouble("precision", "Precision", 0, fPrecision);
  pc.defineString("mintype", "Minimizer", 0, fMinimizerType.c_str());
  pc.defineString("minalg", "Minimizer", 1, fMinimizerAlgo.c_str());
  pc.defineObject("minosSet", "Minos", 0, 0);
  pc.defineObject("condSet", "Cond", 0, 0);
  pc.defineObject("scanSet", "Scan", 0, 0);
  pc.defineSet("cPars", "Constrain", 0, 0);
  pc.defineMutex("Scan", "Minos");

  pc.process(fFitCmdList);
  if (!pc.ok(kTRUE)) {
    throw std::runtime_error("unable to parse cmd list!");
  }

  fOptConst = pc.getInt("optConst");
  fOffset = pc.getInt("doOffset");
  fVerbose = pc.getInt("verbose");
  fSave = pc.getInt("doSave");
  fTimer = pc.getInt("doTimer");
  fPrintLevel = pc.getInt("plevel");
  fDefaultStrategy = pc.getInt("strat");
  fHesse = pc.getInt("hesse");
  fMinos = pc.getInt("minos");
  fScan = pc.getInt("scan");
  fNumee = pc.getInt("numee");
  fDoEEWall = pc.getInt("doEEWall");
  fRetry = pc.getInt("retry");
  fEigen = pc.getInt("eigen");
  fReuseMinimizer = pc.getInt("reminim");
  fReuseNLL = pc.getInt("renll");
  fEps = pc.getDouble("eps");
  fNsigma = pc.getDouble("nsigma");
  fPrecision = pc.getDouble("precision");
  fMinosSet = static_cast<RooArgSet *>(pc.getObject("minosSet"));
  fCondSet = static_cast<RooArgSet *>(pc.getObject("condSet"));
  fScanSet = static_cast<RooArgSet *>(pc.getObject("scanSet"));
  fMinimizerType = std::string(pc.getString("mintype", "Minuit2"));
  fMinimizerAlgo = std::string(pc.getString("minalg", "Migrad"));

  return 1;
}

// ____________________________________________________________________________|__________
RooFitUtils::ExtendedMinimizer::Result::Eigen *
RooFitUtils::ExtendedMinimizer::eigenAnalysis(const TMatrixDSym &hesse) {
  // perform an eigenvector decomposition of the hesse matrix
  int n = hesse.GetNrows();

  // Construct reduced Hessian matrix
  TMatrixDSym Gred(n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Gred(i, j) = hesse(i, j) / sqrt(hesse(i, i) * hesse(j, j));
    }
  }

  // Perform eigenvalue analysis using ROOT standard tools
  TMatrixDSymEigen Geigen(Gred);
  ExtendedMinimizer::Result::Eigen *result =
      new ExtendedMinimizer::Result::Eigen(Geigen.GetEigenValues(),
                                           Geigen.GetEigenVectors());
  return result;
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Minimization
RooFitUtils::ExtendedMinimizer::robustMinimize() {
  // Robust minimization, using an iterative retry strategy
  if (!fMinimizer)
    throw std::runtime_error("no minimizer set!");
  int strategy = fDefaultStrategy;
  int retry = fRetry;
  int status = -1;

  fMinimizer->setPrintLevel(fPrintLevel);
  fMinimizer->optimizeConst(fOptConst);
  fMinimizer->setMinimizerType(fMinimizerType.c_str());
  fMinimizer->setEvalErrorWall(fDoEEWall);
  fMinimizer->setOffsetting(fOffset);
  fMinimizer->setPrintEvalErrors(fNumee);
  fMinimizer->setVerbose(fVerbose);
  fMinimizer->setProfile(fTimer);
  fMinimizer->setStrategy(fDefaultStrategy);
  fMinimizer->setEps(fEps);

  while (true) {
    fMinimizer->setStrategy(strategy);
    coutP(ObjectHandling) << "ExtendedMinimizer::robustMinimize(" << fName
                          << "): starting minimization with strategy "
                          << strategy << std::endl;
    status =
        fMinimizer->minimize(fMinimizerType.c_str(), fMinimizerAlgo.c_str());
    const double nllval = fNll->getVal();

    if ((std::isnan(nllval) || std::isinf(nllval) || (status != 0 && status != 1)) &&
        strategy < 2 && retry > 0) {
      coutW(ObjectHandling)
          << "ExtendedMinimizer::robustMinimize(" << fName
          << ") fit failed with status " << status
          << ". Retrying with strategy " << strategy << std::endl;
      strategy++;
      retry--;
    } else {
      break;
    }
  }

  Result::Minimization mini;
  mini.status = status;
  mini.strategy = strategy;
  mini.nll = fNll->getVal();

  if (!mini.ok()) {
    coutE(ObjectHandling) << "ExtendedMinimizer::robustMinimize(" << fName
                          << ") fit failed with status " << status << std::endl;
  } else {
    coutP(ObjectHandling) << "ExtendedMinimizer::robustMinimize(" << fName
                          << ") fit completed with status " << status
                          << std::endl;
  }
  coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName
                        << "): Evaluating Nll" << std::endl;

  return mini;
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::initialize() {
  // apply all pre-minimization settings
  if (fCondSet) {
    setVals(*fNll->getVariables(), fCondSet, true);
  }
}

// ____________________________________________________________________________|__________

int RooFitUtils::ExtendedMinimizer::minimize(const RooLinkedList &cmdList) {
  // Minimize function  adopted, simplified and extended from RooAbsPdf::fitTo()
  parseConfig(cmdList);

  return minimize();
}

// ____________________________________________________________________________|__________

int RooFitUtils::ExtendedMinimizer::minimize() {
  // Minimize function  adopted, simplified and extended from RooAbsPdf::fitTo()
  initialize();

  setup();

  this->fResult = run();

  return this->fResult->min.status;
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result *RooFitUtils::ExtendedMinimizer::run() {
  // run the actual minimization
  Result *r = new Result();

  r->min = robustMinimize();

  // Evaluate errors with Hesse
  if (fHesse) {
    const int covqual = fMinimizer->lastMinuitFit()->covQual();
    //    TMatrixDSym test = fMinimizer->lastMinuitFit()->covarianceMatrix();
    if (covqual == 0) {
      coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName
                            << "): Covariance quality is " << covqual
                            << ". Running Hesse... " << std::endl;
      fMinimizer->hesse();
    }

    // Obtain Hessian matrix either from patched Minuit or after inversion
    // TMatrixDSym G = Minuit2::MnHesse::lastHessian();
    Double_t determ = 0;
    coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName
                          << "): attempting to invert covariance matrix... "
                          << std::endl;
    TMatrixDSym origG = fMinimizer->lastMinuitFit()->covarianceMatrix();
    TMatrixDSym G = origG.Invert(&determ);
    if (determ == 0 || std::isnan(determ)) {
      coutE(ObjectHandling)
          << "ExtendedMinimizer::minimize(" << fName
          << "): covariance matrix is singular! " << std::endl;
    } else {
      r->hesse = new TMatrixDSym(G);
      // Eigenvalue and eigenvector analysis
      if (fEigen) {
        coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName
                              << "): starting eigen analysis... " << std::endl;
        r->eigen = eigenAnalysis(G);
      }
    }
  }

  if (r->min.status >= 0) {
    // Evaluate errors with Minos
    coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName
                          << "): Running Minos" << std::endl;
    if (fMinos) {
      if (fMinosSet) {
        fMinimizer->minos(*fMinosSet);
      } else {
        fMinimizer->minos();
      }
    }
  }

  if (fScan) {
    coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName
                          << "): Running Scan" << std::endl;
    findSigma(r, *fScanSet);
  } else {
    RooArgSet *vars = fNll->getVariables();
    for (RooLinkedListIter it = vars->iterator();
         RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());) {
      Result::Parameter poi(v->GetName(), v->getVal(), v->getErrorHi(),
                            v->getErrorLo());
      r->parameters.push_back(poi);
    }
  }

  // Return fit result
  if (fSave) {
    std::string name = Form("fitresult_%s_%s", GetName(), fData->GetName());
    std::string title = Form("Result of fit of p.d.f. %s to dataset %s",
                             GetName(), fData->GetName());
    coutP(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName
                          << ") saving results as " << name << std::endl;
    r->fit = fMinimizer->save(name.c_str(), title.c_str());
    if (!r->hesse)
      r->fit->setCovQual(-1);
  }

  if (fCondSet) {
    coutP(ObjectHandling) << "Editing conditional set" << std::endl;
    RooArgSet *attachedSet = fNll->getVariables();
    for (RooLinkedListIter it = fCondSet->iterator();
         RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());) {
      if (RooRealVar *var =
              dynamic_cast<RooRealVar *>(attachedSet->find(v->GetName()))) {
        var->setVal(v->getVal());
        var->setConstant(v->isConstant());
      }
    }
  }

  if (!fReuseMinimizer) {
    delete fMinimizer;
    fMinimizer = NULL;
  }

  return r;
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::scan(
    const std::map<const std::string, std::vector<double>> &params) {
  // perform a scan over the given set of points
  this->scan(this->getResult(true), params);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::scan(
    const std::vector<std::string> &parnames,
    const std::vector<std::vector<double>> &points) {
  // perform a scan over the given set of points
  this->scan(this->getResult(true), parnames, points);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::scan(
    Result *r, const std::map<const std::string, std::vector<double>> &params) {
  // perform a scan over the given set of points
  std::vector<std::string> parnames;
  std::vector<std::vector<double>> parvalues;
  for (auto parit : params) {
    parnames.push_back(parit.first);
  }
  std::vector<std::vector<double>> points;
  std::vector<double> vals(parnames.size());
  addAllPoints(points, parnames, params, vals, 0);
  scan(r, parnames, points);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::scan(
    Result *r, const std::vector<std::string> &parnames,
    const std::vector<std::vector<double>> &points) {
  // perform a scan over the given set of points
  this->setup();
  RooArgSet *attachedSet = fNll->getVariables();

  std::vector<RooRealVar *> params;
  for (auto pname : parnames) {
    RooRealVar *v =
        dynamic_cast<RooRealVar *>(attachedSet->find(pname.c_str()));
    if (!v) {
      throw std::runtime_error(
          TString::Format("unknown parameter name: %s", pname.c_str()).Data());
    }
    params.push_back(v);
  }

  ExtendedMinimizer::Result::Scan scan(parnames);
  for (const auto &point : points) {
    if (point.size() != parnames.size()) {
      throw std::runtime_error("inconsistent vector lengths in scan!");
    }
    for (size_t i = 0; i < params.size(); ++i) {
      params[i]->setVal(point[i]);
      params[i]->setConstant(true);
    }
    auto min = this->robustMinimize();
    if (min.ok()) {
      std::vector<double> vals(params.size());
      for (size_t i = 0; i < params.size(); ++i) {
        vals[i] = params[i]->getVal();
      }
      scan.add(vals, min.nll);
    }
  }
  if (scan.nllValues.size() > 0)
    r->scans.push_back(scan);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::findSigma() {
  // run an iterative algorithm to find the 1-sigma-band
  if (!fNll) {
    throw std::runtime_error("invalid Nll!");
  }
  if (!fScanSet || fScanSet->getSize() == 0) {
    RooArgSet *attachedSet = fNll->getVariables();
    RooStats::RemoveConstantParameters(attachedSet);
    findSigma(this->getResult(true), *attachedSet);
  } else {
    findSigma(this->getResult(true), *fScanSet);
  }
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::findSigma(
    const std::vector<std::string> &pars) {
  // run an iterative algorithm to find the 1-sigma-band
  if (!fNll) {
    throw std::runtime_error("invalid Nll!");
  }
  RooArgList pois;
  RooArgSet *vars = fNll->getVariables();
  for (const auto &p : pars) {
    RooAbsArg *poi = vars->find(p.c_str());
    if (poi) {
      pois.add(*poi);
    }
  }
  findSigma(this->getResult(true), pois);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::findSigma(
    Result *r, const RooAbsCollection &scanSet) {
  // run an iterative algorithm to find the 1-sigma-band
  if (!r->min.ok()) {
    coutW(ObjectHandling) << "ExtendedMinimizer::findSigma(): no previous "
                             "minimization detected, geneating new minimum "
                          << std::endl;
    r->min = robustMinimize();
  }

  for (RooLinkedListIter it = scanSet.iterator();
       RooRealVar *v = dynamic_cast<RooRealVar *>(it.Next());) {
    if (!v) {
      throw std::runtime_error("invalid variable!");
    }
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(): starting scan of "
                          << v->GetName() << std::endl;

    RooArgSet vars(*fPdf->getVariables());
    RooStats::RemoveConstantParameters(&vars);
    vars.add(*v, kTRUE);
    RooArgSet *snap = dynamic_cast<RooArgSet *>(vars.snapshot());
    if (!snap) {
      throw std::runtime_error("invalid snapshot!");
    }

    double nsigma = fabs(fNsigma);
    Double_t val = v->getVal();
    Double_t err = fNsigma * v->getError();
    int maxitr = 25;

    setVals(vars, snap, false);
    Double_t shi = findSigma(r, val + err, val, v, nsigma, maxitr);
    setVals(vars, snap, false);
    Double_t slo = findSigma(r, val - err, val, v, -nsigma, maxitr);
    setVals(vars, snap, false);

    Result::Parameter poi(v->GetName(), val, shi, slo);
    r->parameters.push_back(poi);

    v->setAsymError(std::isnan(slo) ? 1.0 : slo, std::isnan(shi) ? -1.0 : shi);
    coutI(ObjectHandling) << "ExtendedMinimizer::minimize(" << fName << ") "
                          << std::endl;
    delete snap;
  }
}

// ____________________________________________________________________________|__________

Double_t RooFitUtils::ExtendedMinimizer::findSigma(
    RooFitUtils::ExtendedMinimizer::Result *result, const double guessval,
    const double val_mle, RooRealVar *par, const double nsigma,
    const int maxiter) {
  // Find the value of sigma evaluated at a specified nsigma, assuming NLL
  // -2logL is roughly parabolic in par.
  // The precision is specified as a fraction of the error, or based on the
  // Minuit default tolerance.
  // Based on
  // https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/HiggsPhys/HSG3/WWDileptonAnalysisCode/HWWStatisticsCode/trunk/macros/findSigma.C
  // by Aaron Armbruster <aaron.james.armbruster@cern.ch> and adopted by Tim
  // Adye <T.J.Adye@rl.ac.uk>.

  if (!result)
    throw std::runtime_error("cannot scan without previous result!");
  double precision = fPrecision;
  const double fitTol = fEps;
  bool isConst(par->isConstant());

  double val_guess = useLimits(par, guessval);
  int direction = nsigma >= 0.0 ? +1 : -1;
  int nDamping = 1;
  double damping_factor = 1.0;
  ExtendedMinimizer::ValueMap guess_to_corr;
  double tmu = TMath::QuietNaN();

  if (precision <= 0.0) {
    // RooFit default tolerance is 1.0.
    double eps = 0.001 * (fitTol > 0.0 ? fitTol : 1.0);
    precision = 5.0 * eps / (nsigma * nsigma);
  }

  int iter = 0;
  ExtendedMinimizer::Result::Scan values({par->GetName()});
  const double nllmin = result->min.nll;
  values.add({val_mle}, nllmin);
  for (; iter < maxiter; iter++) {
    coutI(ObjectHandling)
        << "ExtendedMinimizer::findSigma(" << fName << ") "
        << Form("Parameter %s %+gsigma iteration %d: start %g (MLE%+g)",
                par->GetName(), nsigma, iter + 1, val_guess,
                val_guess - val_mle)
        << std::endl;
    double val_pre = val_guess;
    par->setVal(val_pre);
    par->setConstant(true);
    Result::Minimization mini = robustMinimize();
    double nll = mini.nll;
    double poival = par->getVal();
    values.add({poival}, nll);

    tmu = 2.0 * (nll - nllmin);
    double sigma_guess = fabs(val_guess - val_mle);
    if (tmu > 0.01)
      sigma_guess /= sqrt(tmu);
    else
      sigma_guess *=
          10.0; // protect against tmu<=0, and also don't move too far
    double corr = damping_factor * (val_pre - val_mle - nsigma * sigma_guess);
    for (ExtendedMinimizer::ValueMap::iterator iguess = guess_to_corr.begin();
         iguess != guess_to_corr.end(); ++iguess) {
      if (fabs(iguess->first - val_pre) < direction * val_pre * 0.02) {
        damping_factor *= 0.8;
        coutW(ObjectHandling)
            << "ExtendedMinimizer::findSigma(" << fName
            << ") Changing damping factor to " << damping_factor << std::endl;
        if (nDamping++ > 10) {
          nDamping = 1;
          damping_factor = 1.0;
        }
        corr *= damping_factor;
        break;
      }
    }
    // subtract off the difference in the new and damped correction
    val_guess -= corr;
    guess_to_corr[val_pre] = corr;
    val_guess = useLimits(par, val_guess);
    double relprecision = precision * fabs(val_guess - val_mle);
    double delta = val_guess - val_pre;

    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName << ") "
                          << Form("%s %.3f (MLE%+.3f) -> %.3f (MLE%+.3f), "
                                  "change %+.3f, precision %.3f, -2lnL %.4f, "
                                  "sigma(guess) %.3f",
                                  par->GetName(), val_pre, val_pre - val_mle,
                                  val_guess, val_guess - val_mle, delta,
                                  relprecision, tmu, sigma_guess)
                          << std::endl;

    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") NLL:                 " << nll << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") delta(NLL):          " << nll - nllmin
                          << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") nsigma*sigma(pre):   "
                          << fabs(val_pre - val_mle) << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") sigma(guess):        " << sigma_guess
                          << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") par(guess):          " << val_guess + corr
                          << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") best-fit val:        " << val_mle << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") tmu:                 " << tmu << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") Precision:           "
                          << direction * val_guess * precision << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") Correction:          " << -corr << std::endl;
    coutI(ObjectHandling) << "ExtendedMinimizer::findSigma(" << fName
                          << ") nsigma*sigma(guess): "
                          << fabs(val_guess - val_mle) << std::endl;

    if (fabs(delta) <= relprecision)
      break;
  }
  bool ok = (!(iter >= maxiter));

  par->setConstant(isConst);

  if (!ok) {
    std::cerr << "findSigma failed after " << iter << " iterations"
              << std::endl;
    return TMath::QuietNaN();
  }

  result->scans.push_back(values);

  double err = val_guess - val_mle;
  coutI(ObjectHandling)
      << "ExtendedMinimizer::findSigma(" << fName << ") "
      << Form("%s %+gsigma = %.3f at -2lnL = %.4f after %d iterations",
              par->GetName(), nsigma, err, tmu, iter + 1)
      << std::endl;

  return err;
}

// _____________________________________________________________________________

RooFitUtils::ExtendedMinimizer::ValueMap
RooFitUtils::ExtendedMinimizer::createProfileValues(RooRealVar *var, double lo,
                                                    double hi, int nbins) {
  // perform a scan over the given set of points
  ExtendedMinimizer::ValueMap map_poi2nll;

  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  RooArgSet *vars = fPdf->getVariables();
  RooStats::RemoveConstantParameters(vars);

  // Unconditional fit
  Result::Minimization mini = robustMinimize();
  map_poi2nll[var->getVal()] = 2 * mini.nll;
  RooArgSet *ucmles = dynamic_cast<RooArgSet *>(vars->snapshot());

  // Perform the scan
  double delta_x = (hi - lo) / nbins;
  for (int i = 0; i <= nbins; i++) {
    *vars = *ucmles;
    var->setVal(lo + i * delta_x);
    var->setConstant(1);

    Result::Minimization mini = robustMinimize();
    map_poi2nll[var->getVal()] = 2 * mini.nll;

    var->setConstant(0);
  }
  return map_poi2nll;
}

// _____________________________________________________________________________

RooFitUtils::ExtendedMinimizer::GraphPair
RooFitUtils::ExtendedMinimizer::createProfile(RooRealVar *var, double lo,
                                              double hi, int nbins) {
  // perform a scan over a single parameter within the given bounds
  ExtendedMinimizer::ValueMap map_poi2nll =
      this->createProfileValues(var, lo, hi, nbins);
  ExtendedMinimizer::GraphPair graphs = prepareProfile(map_poi2nll);
  return graphs;
}

// _____________________________________________________________________________

RooFitUtils::ExtendedMinimizer::GraphPair
RooFitUtils::ExtendedMinimizer::prepareProfile(
    ExtendedMinimizer::ValueMap map_poi2nll) {
  // plot a 1d profile likelhood
  std::vector<double> x, y;
  int nbins = map_poi2nll.size() - 1;

  double xlo = std::numeric_limits<double>::infinity();
  double xhi = -std::numeric_limits<double>::infinity();

  for (ExtendedMinimizer::ValueMap::iterator it_poi = map_poi2nll.begin();
       it_poi != map_poi2nll.end(); ++it_poi) {
    double nll = it_poi->second;
    if (nll == nll && fabs(nll) < pow(10, 20)) {
      x.push_back(it_poi->first);
      y.push_back(it_poi->second);
    }
  }

  int nrPoints = x.size();

  for (int i = 0; i < nrPoints - 1; i++) {
    for (int j = 0; j < nrPoints - 1 - i; j++) {
      if (x[j] > x[j + 1]) {
        std::swap(x[j], x[j + 1]);
        std::swap(y[j], y[j + 1]);
      }
    }
  }

  if (x[0] < xlo)
    xlo = x[0];
  if (x[nrPoints - 1] > xhi)
    xhi = x[nrPoints - 1];

  TGraph *g = new TGraph(nrPoints, getAry(x), getAry(y));

  double minNll = TMath::Infinity();

  for (int i_point = 0; i_point < g->GetN(); ++i_point) {
    double xi, yi = 0;
    g->GetPoint(i_point, xi, yi);
    if (yi < minNll) {
      minNll = yi;
    }
  }

  for (int i_point = 0; i_point < g->GetN(); ++i_point) {
    double xi, yi = 0;
    g->GetPoint(i_point, xi, yi);
    yi -= minNll;
    g->SetPoint(i_point, xi, yi);
  }

  minNll = TMath::Infinity();

  // Make smooth interpolated graph for every folder in poi range, find minimum
  // nll
  std::vector<double> x_interpolated_coarse, y_interpolated_coarse;

  double stepsize_coarse = fabs(xhi - xlo) / nbins;
  for (double thisX = xlo; thisX <= xhi; thisX += stepsize_coarse) {
    double thisY = g->Eval(thisX, 0);
    x_interpolated_coarse.push_back(thisX);
    y_interpolated_coarse.push_back(thisY);
  }

  int nrPoints_interpolated_coarse = x_interpolated_coarse.size();
  TGraph *g_interpolated_coarse =
      new TGraph(nrPoints_interpolated_coarse, getAry(x_interpolated_coarse),
                 getAry(y_interpolated_coarse));

  std::vector<double> x_interpolated, y_interpolated;
  bool twoStepInterpolation = false;

  double stepsize = fabs(xhi - xlo) / (10 * nbins);
  for (double thisX = xlo; thisX <= xhi; thisX += stepsize) {
    double thisY = 0.0;
    if (twoStepInterpolation)
      thisY = g_interpolated_coarse->Eval(thisX, 0, "S");
    else
      thisY = g->Eval(thisX, 0, "S");
    x_interpolated.push_back(thisX);
    y_interpolated.push_back(thisY);
  }

  int nrPoints_interpolated = x_interpolated.size();
  TGraph *g_interpolated = new TGraph(
      nrPoints_interpolated, getAry(x_interpolated), getAry(y_interpolated));

  for (int i_point = 0; i_point < g_interpolated->GetN(); ++i_point) {
    double xi, yi = 0;
    g_interpolated->GetPoint(i_point, xi, yi);
    if (yi < minNll) {
      minNll = yi;
    }
  }

  for (int i_point = 0; i_point < g->GetN(); ++i_point) {
    double xi, yi = 0;
    g->GetPoint(i_point, xi, yi);
    yi -= minNll;
    g->SetPoint(i_point, xi, yi);
  }

  for (int i_point = 0; i_point < g_interpolated->GetN(); ++i_point) {
    double xi, yi = 0;
    g_interpolated->GetPoint(i_point, xi, yi);
    yi -= minNll;
    g_interpolated->SetPoint(i_point, xi, yi);
  }

  g->SetLineWidth(2);
  g->SetMarkerStyle(20);

  g_interpolated->SetLineWidth(2);
  g_interpolated->SetMarkerStyle(20);

  return std::make_pair(g, g_interpolated);
}
