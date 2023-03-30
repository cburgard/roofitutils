#include <limits>
#include <math.h>
#include <cmath>
#include <fstream>  
#include "RooFitUtils/ExtendedMinimizer.h"
#include "Math/Minimizer.h"

#include <iostream>
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TMatrixDSymEigen.h"

#include "RooCmdConfig.h"
#include "RooFitResult.h"
#include "RooAddition.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"
#include "RooRealVar.h"

#include "RooDataHist.h"
#include "RooDataSet.h"
#if ROOT_VERSION_CODE < ROOT_VERSION(6,25,0)
#include "RooMinimizerFcn.h"
#endif

#include "RooStats/RooStatsUtils.h"

#include "RooFitUtils/Utils.h"
#include "RooFitUtils/ExtendedModel.h"

#define nan std::numeric_limits<double>::quiet_NaN()
#define inf std::numeric_limits<double>::infinity()

ClassImp(RooFitUtils::ExtendedMinimizer)

#include "TMinuitMinimizer.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMinuit.h"

namespace {
  //somewhat complex but apparently standard conform hack to access RooMinimizer::getNPar.
  template <typename RooMinimizerTag>
  struct RooMinimizerHackResult {
    typedef typename RooMinimizerTag::type type;
    static type ptr;
  };

  template <typename RooMinimizerTag>
  typename RooMinimizerHackResult<RooMinimizerTag>::type RooMinimizerHackResult<RooMinimizerTag>::ptr;

  template<typename RooMinimizerTag, typename RooMinimizerTag::type p>
  struct RooMinimizerRob : RooMinimizerHackResult<RooMinimizerTag> {
    struct RooMinimizerFiller {
      RooMinimizerFiller() {RooMinimizerHackResult<RooMinimizerTag>::ptr = p;}
    };
    static RooMinimizerFiller RooMinimizerfiller_obj;
  };

  template<typename RooMinimizerTag, typename RooMinimizerTag::type p>
  typename RooMinimizerRob<RooMinimizerTag, p>::RooMinimizerFiller RooMinimizerRob<RooMinimizerTag, p>::RooMinimizerfiller_obj;

  //now expose some members of RooMinimizer that we need to access
  struct RooMinimizergetNPar { typedef Int_t(RooMinimizer::*type)() const; };
  template class RooMinimizerRob<RooMinimizergetNPar, &RooMinimizer::getNPar>;

  // struct RooMinimizerfitterFcn { typedef RooMinimizerFcn*(RooMinimizer::*type)(); };
  // template class RooMinimizerRob<RooMinimizerfitterFcn, &RooMinimizer::fitterFcn>;

}


namespace {
  //somewhat complex but apparently standard conform hack to access Fitter functions.
  template <typename RooFitterTag>
  struct RooFitterHackResult {
    typedef typename RooFitterTag::type type;
    static type ptr;
  };

  template <typename RooFitterTag>
  typename RooFitterHackResult<RooFitterTag>::type RooFitterHackResult<RooFitterTag>::ptr;

  template<typename RooFitterTag, typename RooFitterTag::type p>
  struct RooFitterRob : RooFitterHackResult<RooFitterTag> {
    struct RooFitterFiller {
      RooFitterFiller() {RooFitterHackResult<RooFitterTag>::ptr = p;}
    };
    static RooFitterFiller RooFitterfiller_obj;
  };

  template<typename RooFitterTag, typename RooFitterTag::type p>
  typename RooFitterRob<RooFitterTag, p>::RooFitterFiller RooFitterRob<RooFitterTag, p>::RooFitterfiller_obj;

  //now expose some members of RooFitter that we need to access
  struct RooFitterInit { typedef bool(ROOT::Fit::Fitter::*type)(); };
  template class RooFitterRob<RooFitterInit, &ROOT::Fit::Fitter::DoInitMinimizer>;
}


namespace {
  //somewhat complex but apparently standard conform hack to access MinuitMinimizer functions.
  template <typename RooMinuitMiniTag>
  struct RooMinuitMiniHackResult {
    typedef typename RooMinuitMiniTag::type type;
    static type ptr;
  };

  template <typename RooMinuitMiniTag>
  typename RooMinuitMiniHackResult<RooMinuitMiniTag>::type RooMinuitMiniHackResult<RooMinuitMiniTag>::ptr;

  template<typename RooMinuitMiniTag, typename RooMinuitMiniTag::type p>
  struct RooMinuitMiniRob : RooMinuitMiniHackResult<RooMinuitMiniTag> {
    struct RooMinuitMiniFiller {
      RooMinuitMiniFiller() {RooMinuitMiniHackResult<RooMinuitMiniTag>::ptr = p;}
    };
    static RooMinuitMiniFiller RooMinuitMinifiller_obj;
  };

  template<typename RooMinuitMiniTag, typename RooMinuitMiniTag::type p>
  typename RooMinuitMiniRob<RooMinuitMiniTag, p>::RooMinuitMiniFiller RooMinuitMiniRob<RooMinuitMiniTag, p>::RooMinuitMinifiller_obj;

  //now expose some members of RooMinuitMini that we need to access
  struct RooMinuitMiniMinuit { typedef TMinuit*(TMinuitMinimizer::*type); };
  template class RooMinuitMiniRob<RooMinuitMiniMinuit, &TMinuitMinimizer::fMinuit>;
}




namespace {
  //somewhat complex but apparently standard conform hack to access RooFitResult::setCovQual.
  template <typename RooFitResultTag>
  struct RooFitResultHackResult {
    typedef typename RooFitResultTag::type type;
    static type ptr;
  };

  template <typename RooFitResultTag>
  typename RooFitResultHackResult<RooFitResultTag>::type RooFitResultHackResult<RooFitResultTag>::ptr;

  template<typename RooFitResultTag, typename RooFitResultTag::type p>
  struct RooFitResultRob : RooFitResultHackResult<RooFitResultTag> {
    struct RooFitResultFiller {
      RooFitResultFiller() {RooFitResultHackResult<RooFitResultTag>::ptr = p;}
    };
    static RooFitResultFiller RooFitResultfiller_obj;
  };

  template<typename RooFitResultTag, typename RooFitResultTag::type p>
  typename RooFitResultRob<RooFitResultTag, p>::RooFitResultFiller RooFitResultRob<RooFitResultTag, p>::RooFitResultfiller_obj;

  //now expose some members of RooFitResult that we need to access
  struct RooFitResultsetCovQual { typedef void(RooFitResult::*type)(Int_t); };
  template class RooFitResultRob<RooFitResultsetCovQual, &RooFitResult::setCovQual>;
}


namespace {
  //somewhat complex but apparently standard conform hack to access RooAbsPdf::_norm.
  template <typename RooAbsPdfTag>
  struct RooAbsPdfHackResult {
    typedef typename RooAbsPdfTag::type type;
    static type ptr;
  };

  template <typename RooAbsPdfTag>
  typename RooAbsPdfHackResult<RooAbsPdfTag>::type RooAbsPdfHackResult<RooAbsPdfTag>::ptr;

  template<typename RooAbsPdfTag, typename RooAbsPdfTag::type p>
  struct RooAbsPdfRob : RooAbsPdfHackResult<RooAbsPdfTag> {
    struct RooAbsPdfFiller {
      RooAbsPdfFiller() {RooAbsPdfHackResult<RooAbsPdfTag>::ptr = p;}
    };
    static RooAbsPdfFiller RooAbsPdffiller_obj;
  };

  template<typename RooAbsPdfTag, typename RooAbsPdfTag::type p>
  typename RooAbsPdfRob<RooAbsPdfTag, p>::RooAbsPdfFiller RooAbsPdfRob<RooAbsPdfTag, p>::RooAbsPdffiller_obj;

  //now expose some members of RooAbsPdf that we need to access
  struct RooAbsPdf_norm { typedef RooAbsReal*(RooAbsPdf::*type); };
  template class RooAbsPdfRob<RooAbsPdf_norm, &RooAbsPdf::_norm>;
}


#if ROOT_VERSION_CODE < ROOT_VERSION(6,10,0)
#define FitterFcn RooMinimizerFcn
#else
#define FitterFcn auto
#endif


#if ROOT_VERSION_CODE < ROOT_VERSION(6,25,0)
namespace {
  class RooMinimizerHack : public RooMinimizer{
  public:
    FitterFcn* getFitterFcn(){ return this->fitterFcn(); }
  };
}
#endif


namespace {
  int countFloatParams(RooArgSet* args){
    int n = 0;
    RooAbsArg* obj;
    RooFIter itr(args->fwdIterator());
    while((obj = itr.next())){
      RooRealVar* v = dynamic_cast<RooRealVar*>(obj);
      if(!v) continue;
      if(!v->isConstant()) ++n;
    }
    return n;
  }
  int countFloatParams(RooMinimizer* minimizer){
    return (minimizer->*RooMinimizerHackResult<RooMinimizergetNPar>::ptr)();
  }
  FitterFcn* fitterFcn(RooMinimizer* minimizer){
#if ROOT_VERSION_CODE < ROOT_VERSION(6,25,0)
    return ((::RooMinimizerHack*) minimizer)->getFitterFcn();
#else
    return ((::RooMinimizer*) minimizer)->getMultiGenFcn();
#endif
  }
  bool init(ROOT::Fit::Fitter* fitter){
    return (fitter->*RooFitterHackResult<RooFitterInit>::ptr)();    
  }
  TMinuit* getMinuit(ROOT::Math::Minimizer* mini){
    TMinuitMinimizer* minuit = dynamic_cast<TMinuitMinimizer*>(mini);
    if(minuit){
      return (minuit->*RooMinuitMiniHackResult<RooMinuitMiniMinuit>::ptr);
    }
    return NULL;
  }
  template<class State> TMatrixDSym getHessian(const State& fState, int fDim) {
    // get value of Hessian matrix
    TMatrixDSym hess(fDim);
    // this is the second derivative matrices
    auto hessian = fState.Hessian();
    if (!fState.HasCovariance())
      throw std::runtime_error("no hessian available, minimization failed");
    for (unsigned int i = 0; i < fDim; ++i) {
      if (fState.Parameter(i).IsFixed() || fState.Parameter(i).IsConst()) {
        for (unsigned int j = 0; j < fDim; ++j) {
          hess[i][j] = 0;
        }
      } else {
        unsigned int l = fState.IntOfExt(i);
        for (unsigned int j = 0; j < fDim; ++j) {
          // could probably speed up this loop (if needed)
          if (fState.Parameter(j).IsFixed() || fState.Parameter(j).IsConst())
            hess[i][j] = 0;
          else {
            // need to transform from external to internal indices)
            // for taking care of the removed fixed row/columns in the Minuit2 representation
            unsigned int m = fState.IntOfExt(j);
            hess[i][j] = hessian(l, m);
          }
        }
      }
    }
    return hess;
  }
}

// some toggles
#define USE_ROOFITRESULT_NLL
#define FINDSIGMA_NO_CONSTOPT    

int RooFitUtils::ExtendedMinimizer::runHesse(RooFitUtils::ExtendedMinimizer::Result::Minimization& mini){
  // run the HESSE algorithm
  this->SetNllDirty();
  int ndim = mini.ndim;
  int hesseStatus = 0;
  auto* fitter = fMinimizer->fitter();        
  auto* minimizer = fitter->GetMinimizer();
  if(minimizer){
    // the simple case - just run hesse after minimization
    std::cout << "ExtendedMinimizer::runHesse(" << fName  << ") running after minimization (this might take a while) ... "<<std::endl;
    hesseStatus = fMinimizer->hesse();
  } else {
    // the more complicated case - not at minimum, need to do some gynmastics to get access to hesse matrix
    std::cout << "ExtendedMinimizer::runHesse(" << fName << ") running standalone (this might take a while) ... "<<std::endl;        
    auto* fcn = fitterFcn(fMinimizer);
    #if ROOT_VERSION_CODE < ROOT_VERSION(6,25,0)
    fcn->Synchronize(fitter->Config().ParamsSettings(),fOptConst,0);
    #endif
    fitter->SetFCN(*fcn);
    fitter->EvalFCN();          
    init(fitter);
    minimizer = fitter->GetMinimizer();
    minimizer->Hesse();
  }
  mini.covqual = minimizer->CovMatrixStatus();
  std::cout << "ExtendedMinimizer::runHesse(" << fName << ") finished with status "<< hesseStatus<< " (covqual=" << mini.covqual << ")" << std::endl;
  double det;
  TMinuit* minuit = getMinuit(minimizer);
  ROOT::Minuit2::Minuit2Minimizer* minuit2 = dynamic_cast<ROOT::Minuit2::Minuit2Minimizer*>(minimizer);
  if(minuit){
    mini.cov = new TMatrixDSym(ndim);
    minuit->mnemat(mini.cov->GetMatrixArray(),ndim);
    mini.hesse = new TMatrixDSym(mini.cov->Invert(&det));
  } else if(minuit2){
    mini.cov = new TMatrixDSym(ndim);
    fitter->Result().GetCovarianceMatrix(*mini.cov);
    mini.hesse = new TMatrixDSym(mini.cov->Invert(&det));    
  } else {
    mini.hesse = new TMatrixDSym(ndim);
    minimizer->GetHessianMatrix(mini.hesse->GetMatrixArray());
    mini.cov = new TMatrixDSym(mini.hesse->Invert(&det));
  }
  if(!(det > 0)){
    TMatrixDSymEigen m(*mini.hesse);
    double newdet = 1;
    for(int i=0; i<m.GetEigenValues().GetNrows(); ++i){
      newdet *= m.GetEigenValues()[i];
      if(m.GetEigenValues()[i] == 0.){
        std::cout << "ExtendedMinimizer::runHesse(" << fName << ") identified row "<<i<< " (parameter " << mini.parameters[i].name << ") as problematic" << std::endl;
      }
    }
    if(newdet == det){
      std::cout << "ExtendedMinimizer::runHesse(" << fName << ") Hesse matrix found invalid (D=" << det << "), deleting" << std::endl;
      delete mini.hesse;
      delete mini.cov;
      mini.hesse = NULL;
      mini.cov = NULL;
      hesseStatus = -1;
    } else {
      double invdet;
      TMatrixDSym closure(mini.cov->Invert(&invdet));
      bool ok = true;
      for(int i=0; i<closure.GetNcols(); ++i){
        for(int j=0; j<closure.GetNcols(); ++j){
          if(closure(i,j) != (*mini.hesse)(i,j)){
            ok = false;
          }
        }
      }
      if(ok) {
        std::cout << "ExtendedMinimizer::runHesse(" << fName << ") Hesse matrix determinant subject to numerical fluctuations (D=" << det << " or " << newdet << "), but closure seems fine, continuing!" << std::endl;
      } else {
        std::cout << "ExtendedMinimizer::runHesse(" << fName << ") Hesse matrix determinant subject to numerical fluctuations (D=" << det << " or " << newdet << ") and closure failed, reconstructing from eigenvectors!" << std::endl;
        TMatrixDSym evals(mini.hesse->GetNcols());
        TMatrixDSym invevals(mini.hesse->GetNcols());        
        TMatrixD evecs(m.GetEigenVectors());
        for(int i=0; i<mini.hesse->GetNcols(); ++i){
          evals(i,i) = m.GetEigenValues()[i];
          invevals(i,i) = 1./m.GetEigenValues()[i];
        }
        delete mini.cov;
        TMatrixD cov(evecs.T() * invevals * evecs);
        mini.cov = new TMatrixDSym(cov.GetNcols());
        for(int i=0; i<cov.GetNcols(); ++i){
          for(int j=0; j<cov.GetNcols(); ++j){
            (*mini.cov)(i,j) = cov(i,j);
          }
        }
      }
    }
  }
  
  return hesseStatus;
}


// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Result()
  : eigen(NULL) {
  // nothing here
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Minimization::Minimization(){
  // nothing here
}

// ____________________________________________________________________________|__________

bool RooFitUtils::ExtendedMinimizer::Result::Minimization::ok(int status) {
  // return true if status is 0 or 1, false otherwise
  return (status == 0 || status == 1);
}

// ____________________________________________________________________________|__________



// ____________________________________________________________________________|__________

bool RooFitUtils::ExtendedMinimizer::Result::Minimization::ok() const {
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

RooFitUtils::ExtendedMinimizer::Result::Scan::Scan(const std::string& n,
    const std::vector<std::string> &parnames)
     : name(n), parNames(parnames) {
  // nothing here
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result::Scan::Scan(const std::string& n,
    const std::vector<std::string> &parnames, const std::vector<std::string> &extraparnames)
      : name(n), parNames(parnames), extraParNames(extraparnames) {
  // nothing here
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::Result::Scan::add(const std::vector<double> &parvals, int fitstatus, double nllval) {
  // add a new entry to a scan
  if (parvals.size() != this->parNames.size()) {
    throw std::runtime_error("cannot add parameter list with wrong length!");
  }
  this->parValues.push_back(parvals);
  this->fitStatus.push_back(fitstatus);
  this->nllValues.push_back(nllval);
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::Result::Scan::add(const std::vector<double> &parvals, int fitstatus, double nllval, const std::vector<double>& extravals) {
  // add a new entry to a scan
  if (parvals.size() != this->parNames.size()) {
    throw std::runtime_error("cannot add parameter list with wrong length!");
  }
  if (extravals.size() != this->extraParNames.size()) {
    throw std::runtime_error("cannot add extra parameter list with wrong length!");
  }
  this->parValues.push_back(parvals);
  this->extraParValues.push_back(extravals);  
  this->fitStatus.push_back(fitstatus);
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
  for (size_t i = 0; i < this->extraParNames.size(); ++i) {
    std::cout << extraParNames[i] << "\t";
  }  
  std::cout << "nll" << std::endl;
  for (size_t j = 0; j < this->nllValues.size(); ++j) {
    for (size_t i = 0; i < this->parNames.size(); ++i) {
      std::cout << parValues[j][i] << "\t";
    }
    for (size_t i = 0; i < this->extraParNames.size(); ++i) {
      std::cout << extraParValues[j][i] << "\t";
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
    return this->fResult->min.fit;
  return NULL;
}

// ____________________________________________________________________________|__________

TMatrixDSym RooFitUtils::ExtendedMinimizer::GetHesseMatrix() {
  // obtain the hesse matrix from the minimizer
  if (this->fResult && this->fResult->min.hesse)
    return *(this->fResult->min.hesse);
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

  RooLinkedList *makeList(bool owned,
			  const RooCmdArg &arg1 = RooCmdArg::none(),
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
    l->Add(!owned?arg1 .Clone():const_cast<RooCmdArg*>(&arg1 ));
    l->Add(!owned?arg2 .Clone():const_cast<RooCmdArg*>(&arg2 ));
    l->Add(!owned?arg3 .Clone():const_cast<RooCmdArg*>(&arg3 ));
    l->Add(!owned?arg4 .Clone():const_cast<RooCmdArg*>(&arg4 ));
    l->Add(!owned?arg5 .Clone():const_cast<RooCmdArg*>(&arg5 ));
    l->Add(!owned?arg6 .Clone():const_cast<RooCmdArg*>(&arg6 ));
    l->Add(!owned?arg7 .Clone():const_cast<RooCmdArg*>(&arg7 ));
    l->Add(!owned?arg8 .Clone():const_cast<RooCmdArg*>(&arg8 ));
    l->Add(!owned?arg9 .Clone():const_cast<RooCmdArg*>(&arg9 ));
    l->Add(!owned?arg10.Clone():const_cast<RooCmdArg*>(&arg10));
    l->Add(!owned?arg11.Clone():const_cast<RooCmdArg*>(&arg11));
    l->Add(!owned?arg12.Clone():const_cast<RooCmdArg*>(&arg12));
    l->Add(!owned?arg13.Clone():const_cast<RooCmdArg*>(&arg13));
    l->Add(!owned?arg14.Clone():const_cast<RooCmdArg*>(&arg14));
    l->Add(!owned?arg15.Clone():const_cast<RooCmdArg*>(&arg15));
    l->Add(!owned?arg16.Clone():const_cast<RooCmdArg*>(&arg16));
    l->Add(!owned?arg17.Clone():const_cast<RooCmdArg*>(&arg17));
    l->Add(!owned?arg18.Clone():const_cast<RooCmdArg*>(&arg18));
    l->Add(!owned?arg19.Clone():const_cast<RooCmdArg*>(&arg19));
    l->Add(!owned?arg20.Clone():const_cast<RooCmdArg*>(&arg20));
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

  inline void clearContents(RooLinkedList& l, bool owned){
    // clear the contents of the list, deleting the objects if rquested
    if(owned){
      RooLinkedListIter it = l.iterator();
      while (true) {
	TObject *obj = it.Next();
	if (!obj)
	  break;
	delete obj;
      }
    }
    l.Clear();
  }

  inline int addAllArgs(const RooLinkedList &orig, RooLinkedList &target, bool clone) {
    // add all arguments from the original list to the target list, making clones if requested
    int n = 0;
    RooLinkedListIter it = orig.iterator();
    while (true) {
      TObject *obj = it.Next();
      if (!obj)
	break;
      RooCmdArg *v = dynamic_cast<RooCmdArg *>(obj);
      if (!v)
	break;
      if (strlen(v->GetName()) != 0){
	if(clone)
	  target.Add(v->Clone());
	else
	  target.Add(v);
	n++;
      }
    }
    return n;
  }

  void inverseFilterCmdList(const RooLinkedList &cmdInList, RooLinkedList& filteredList, const char *cmdNameList, bool clone) {
    // add all arguments passing the filter from the original list to the target list, making clones if requested
    char buf[1024];
    strlcpy(buf, cmdNameList, 1024);

    char *name = strtok(buf, ",");
    while (name) {
      RooCmdArg *cmd = (RooCmdArg*)(cmdInList.FindObject(name));
      if (cmd) {
        RooCmdArg* cl = (RooCmdArg*)(clone ? cmd->Clone() : cmd);
        filteredList.Add(cl);
      }
      name = strtok(0, ",");
    }
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

RooFitUtils::ExtendedMinimizer::ExtendedMinimizer(const char* minimizerName, RooFitUtils::ExtendedModel* model,
                                                  const RooLinkedList &argList)
	: ExtendedMinimizer(minimizerName,
                            model->GetPdf(),
                            model->GetData(),
                            model->GetWorkspace()){
  fModel = model;
  // Constructor
  RooLinkedList newargList(argList);
  if(model->GetGlobalObservables()){
    fOwnedArgs.push_back(RooFit::GlobalObservables(*(model->GetGlobalObservables())));
    newargList.Add(&fOwnedArgs.at(fOwnedArgs.size()-1));
  }
  fPenaltyMini = model->GetPenalty();
  parseNllConfig(newargList);
  parseFitConfig(newargList);
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::ExtendedMinimizer(const char* minimizerName, RooFitUtils::ExtendedModel* model)
	: ExtendedMinimizer(minimizerName,
                            model->GetPdf(),
                            model->GetData(),
                            model->GetWorkspace()){
  fModel = model;
  // Constructor
  RooLinkedList argList;
  if(model->GetGlobalObservables()){  
    fOwnedArgs.push_back(RooFit::GlobalObservables(*(model->GetGlobalObservables())));
    argList.Add(&fOwnedArgs.at(fOwnedArgs.size()-1));
  }
  fPenaltyMini = model->GetPenalty();
  parseNllConfig(argList);
}


// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::ExtendedMinimizer(const char *minimizerName,
                                                  RooAbsPdf *pdf,
                                                  RooAbsData *data,
                                                  RooWorkspace *workspace,
                                                  const RooLinkedList &argList)
    : ExtendedMinimizer(minimizerName, pdf, data, workspace) {
  // Constructor

  parseNllConfig(argList);
  parseFitConfig(argList);
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::ExtendedMinimizer(const char *minimizerName,
                                                  RooAbsPdf *pdf,
                                                  RooAbsData *data,
                                                  RooWorkspace *workspace)
    : TNamed(minimizerName, minimizerName),
      fWorkspace(workspace), fPdf(pdf), fData(data),
      fOffset(0), fOptConst(2), fVerbose(0), fSave(0), fTimer(1),
      fPrintLevel(1), fDefaultStrategy(0), fHesse(0), fMinimize(1), fMinos(0),
      fNumee(5), fDoEEWall(1), fRetry(0), fEigen(0), fReuseMinimizer(0),
      fReuseNLL(0), fChi2(0), fMaxCalls(10000), fMaxIterations(10000), fEps(1.0), // 1sigma 1dof
      fMinimizerType("Minuit2"), fMinimizerAlgo("Migrad") {
  // Constructor

  fNllCmdList.SetName("NllCmdList");
  fFitCmdList.SetName("FitCmdList");

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(),
                                                    fMinimizerAlgo.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(fPrintLevel);

  coutP(InputArguments) << "ExtendedMinimizer::ExtendedMinimizer(" << fName
                        << ") created" << std::endl;
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::~ExtendedMinimizer() {
  // Destructor
//  clearContents(this->fNllCmdList,true);
//  clearContents(this->fFitCmdList,true);
//  if(fNll) delete fNll;
//  if(fMinimizer) delete fMinimizer;
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result RooFitUtils::ExtendedMinimizer::minimize(
    const RooCmdArg &arg1, const RooCmdArg &arg2, const RooCmdArg &arg3,
    const RooCmdArg &arg4, const RooCmdArg &arg5, const RooCmdArg &arg6,
    const RooCmdArg &arg7, const RooCmdArg &arg8, const RooCmdArg &arg9,
    const RooCmdArg &arg10, const RooCmdArg &arg11, const RooCmdArg &arg12) {
  // Minimize function with iterative retry strategy adopted, simplified and
  // extended from RooAbsPdf::fitTo()
  RooLinkedList *l = makeList(false, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8,
                              arg9, arg10, arg11, arg12);
  auto result = minimize(*l);
  delete l;
  return result;
}

// ____________________________________________________________________________|__________

namespace {
  void printSet(std::ostream& out, const RooArgSet* set, const char* wsname){
    out << "RooArgSet(";
    bool first = true;
    for(const auto& e:*set){
      if(!first) out << ",";
      out << "*" << wsname << "->arg(\"" << e->GetName() << "\")";
      first = false;
    }
    out << ")";
  }
		
  void printCmd(std::ostream& out, const RooCmdArg* arg, const char* wsname){
    out << "RooCmdArg(\"" << arg->GetName() << "\", " <<
      arg->getInt(0) << ", " << arg->getInt(1) << ", " <<
      arg->getDouble(0) << ", " << arg->getDouble(1) << ", ";
    if(arg->getString(0)) out << '"' << arg->getString(0) << '"' << ", ";
    else out << "\"\", ";
    if(arg->getString(1)) out << '"' << arg->getString(1) << '"' << ", ";    
    else out << "\"\", ";
    if(arg->getObject(0)) const_cast<TObject*>(arg->getObject(0))->SavePrimitive(out);
    else out << " 0";
    out << ", ";
    if(arg->getObject(1)) const_cast<TObject*>(arg->getObject(1))->SavePrimitive(out);
    else out << " 0";
    out << ", ";
    out << " 0, ";
    if(arg->getString(2)) out << '"' << arg->getString(2) << '"' << ", ";    
    else out << "\"\", ";    
    if(arg->getSet(0)){ out << "new "; ::printSet(out,arg->getSet(0),wsname); }
    else out << " 0";
    out << ", ";
    if(arg->getSet(1)){ out << "new "; ::printSet(out,arg->getSet(1),wsname); }
    else out << " 0";        
    out << ")";
  }
}

void RooFitUtils::ExtendedMinimizer::writeReproducer(const std::string& name) {
  if(!this->fWorkspace){
    return;
  }
  this->fWorkspace->writeToFile((name+".root").c_str());  
  std::ofstream out(name+"_fit.C");
  out << "{\n";
  out << "  TFile* infile = TFile::Open(\"" << name << ".root\",\"READ\");\n";
  out << "  RooWorkspace* workspace = (RooWorkspace*)infile->Get(\""<<fWorkspace->GetName()<<"\");\n";
  out << "  RooAbsPdf* pdf = workspace->pdf(\""<< fPdf->GetName() << "\");\n";
  out << "  RooAbsData* data = workspace->data(\""<< fData->GetName() << "\");\n";
  out << "  RooAbsReal* nll = pdf->createNLL(*data";
  for(const auto& arg:fNllCmdList){
    out << ",\n    ";
    ::printCmd(out,dynamic_cast<RooCmdArg*>(arg),"workspace");
  }
  out << "  );\n";
  out << "  RooMinimizer mini(*nll);\n";
  out << "  mini.setMaxFunctionCalls("<<fMaxCalls<<");\n";
  out << "  mini.setMaxIterations("<<fMaxIterations<<");\n";
  out << "  mini.setPrintLevel("<<fPrintLevel<<");\n";
  out << "  mini.optimizeConst("<<fOptConst<<");\n";
  out << "  mini.setEvalErrorWall("<<fDoEEWall<<");\n";
  out << "  mini.setOffsetting("<<fOffset<<");\n";
  out << "  mini.setPrintEvalErrors("<<fNumee<<");\n";
  out << "  mini.setVerbose("<<fVerbose<<");\n";
  out << "  mini.setProfile("<<fTimer<<");\n";
  out << "  mini.setStrategy("<<fDefaultStrategy<<");\n";
  out << "  mini.setEps("<<fEps<<");\n";
  out << "  mini.minimize(\"" << fMinimizerType << "\",\"" << fMinimizerAlgo << "\");";
  out << "}\n";
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

  if (!fReuseNLL) {
    if(fNll){
      coutW(InputArguments) << "deleting previous NLL!" << std::endl;
      delete fNll;
    }
    fNll = NULL;
  }
  int ndim = 0;

  if(!fNll){
    coutI(InputArguments) << "Creating new Nll" << std::endl;

    if (fWorkspace) {
      if (RooFitUtils::RooStarMomentMorphFix) {
        int n = fixRooStarMomentMorph(fWorkspace);
        if (n > 0)
          coutP(InputArguments)
            << "Fixed cache of " << n << " instances of RooStarMomentMorph"
            << std::endl;
      }
    }
    double nllval = 0.;
    try {
      if(fChi2){
	if(fData->InheritsFrom(RooDataSet::Class())){
	  fNll = fPdf->createChi2(*(RooDataSet*)fData, fNllCmdList);
	} else if(fData->InheritsFrom(RooDataHist::Class())){
	  fNll = fPdf->createChi2(*(RooDataHist*)fData, fNllCmdList);
	}
      } else {
	fNll = fPdf->createNLL(*fData, fNllCmdList);
      }
      nllval = fNll->getVal();
    } catch (std::exception& ex){
      throw ex;
    } catch (std::string& s){
      throw std::runtime_error(s);
    } catch (const char*& s){
      throw std::runtime_error(s);
    } catch (int& code){
      throw std::runtime_error(TString::Format("encountered error code %d in RooAbsPdf::createNLL",code).Data());
    } catch (...){
      throw std::runtime_error("encountered unknown exception type in RooAbsPdf::createNLL");
    }

    if(std::isinf(nllval)){
      throw std::runtime_error("starting value of nll is inf!");
    }
    RooArgSet* args = fNll->getVariables();
    ndim = ::countFloatParams(args);
    delete args;
    coutI(InputArguments) << "Created new Nll with " << ndim << " parameters, starting value " << nllval << std::endl;
  } else {
    RooArgSet* args = fNll->getVariables();
    ndim = ::countFloatParams(args);
    delete args;
    coutI(InputArguments) << "Using existing Nll with " << ndim << " parameters" << std::endl;
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
    throw std::runtime_error("ExtendedMinimizer::setup: Failed to obtain NLL");
  }

  if(!fReuseMinimizer){
    if(fMinimizer){
      coutW(InputArguments) << "deleting previous Minimizer!" << std::endl;
      delete fMinimizer;
    }
    fMinimizer=NULL;
  }

  if (!fMinimizer) {
    coutI(InputArguments) << "Creating new Minimizer" << std::endl;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(fMinimizerType.c_str(),
                                                      fMinimizerAlgo.c_str());
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(fDefaultStrategy);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(fPrintLevel);

    fMinimizer = new RooMinimizer(*fNll);
    int npar = countFloatParams(fMinimizer);
    if(ndim > npar){
      throw std::runtime_error(TString::Format("construction of minimizer failed, number of parameters does not match between Nll (npar=%d) and minimizer (npar=%d)!",ndim,npar).Data());
    }
  } else {
    int npar = countFloatParams(fMinimizer);
    if(ndim > npar){
      throw std::runtime_error(TString::Format("minimizer seems to have wrong Nll set, number of parameters does not match between Nll (npar=%d) and minimizer (npar=%d)!",ndim,npar).Data());
    }
    coutI(InputArguments) << "Using existing Minimizer" << std::endl;
  }
}

// ____________________________________________________________________________|__________

template <class A>
int RooFitUtils::ExtendedMinimizer::parseNllConfig(const A &cmdList) {
  if(!this->fNll){
    clearContents(fNllCmdList,true);
    inverseFilterCmdList(cmdList, fNllCmdList,
			 "NumCPU,BatchMode,Constrained,Constrain,CloneData,"
			 "GlobalObservables,GlobalObservablesTag,"
			 "OffsetLikelihood",true);
  } else {
    std::cout << "cannot change Nll config with preexisting Nll!" << std::endl;
  }
  if(fPenaltyMini){
    fOwnedArgs.push_back(RooFit::ExternalConstraints(*fPenaltyMini));
    fNllCmdList.Add(&(fOwnedArgs[fOwnedArgs.size()-1]));
  }
  return fNllCmdList.GetSize();
}
template
int RooFitUtils::ExtendedMinimizer::parseNllConfig<RooLinkedList>(const RooLinkedList &cmdList);

// ____________________________________________________________________________|__________

template <class A>
int RooFitUtils::ExtendedMinimizer::parseFitConfig(const A &cmdList) {
  // parse the configuration
  RooCmdConfig pc(Form("ExtendedMinimizer::parseFitConfig(%s)", GetName()));
  pc.allowUndefined();

  clearContents(fFitCmdList,true);
  addAllArgs(cmdList, fFitCmdList, true);

  pc.defineInt("optConst", "Optimize", 0, fOptConst);
  pc.defineInt("doOffset", "OffsetLikelihood", 0, fOffset);
  pc.defineInt("verbose", "Verbose", 0, fVerbose);
  pc.defineInt("doSave", "Save", 0, fSave);
  pc.defineInt("doTimer", "Timer", 0, fTimer);
  pc.defineInt("plevel", "PrintLevel", 0, fPrintLevel);
  pc.defineInt("strat", "Strategy", 0, fDefaultStrategy);
  pc.defineInt("hesse", "Hesse", 0, fHesse);
  pc.defineInt("minimize", "Minimize", 0, fMinimize);  
  pc.defineInt("minos", "Minos", 0, fMinos);
  pc.defineInt("numee", "PrintEvalErrors", 0, fNumee);
  pc.defineInt("doEEWall", "EvalErrorWall", 0, fDoEEWall);
  pc.defineInt("retry", "NumRetryFit", 0, fRetry);
  pc.defineInt("maxcalls", "MaxCalls", 0, fMaxCalls);
  pc.defineInt("maxiterations", "MaxIterations", 0, fMaxIterations);
  pc.defineInt("eigen", "Eigen", 0, fEigen);
  pc.defineInt("reminim", "ReuseMinimizer", 0, fReuseMinimizer);
  pc.defineInt("renll", "ReuseNLL", 0, fReuseNLL);
  pc.defineInt("usechi2", "UseChi2", 0, fChi2);  
  pc.defineDouble("eps", "Eps", 0, fEps);
  pc.defineString("mintype", "Minimizer", 0, fMinimizerType.c_str());
  pc.defineString("minalg", "Minimizer", 1, fMinimizerAlgo.c_str());
  pc.defineObject("minosObjSet", "Minos", 0, 0);
  pc.defineSet("minosSet", "Minos", 0, 0);  
  pc.defineObject("condSet", "Cond", 0, 0);
  pc.defineSet("cPars", "Constrain", 0, 0);

  pc.defineInt("findSigma", "FindSigma", 0, fFindSigma);
  pc.defineInt("findSigmaIter", "FindSigma", 1, fFindSigmaIter);    
  pc.defineDouble("findSigmaN", "FindSigma", 0, fFindSigmaN);
  pc.defineDouble("findSigmaPrecision", "FindSigma", 1, fFindSigmaPrecision);  
  pc.defineObject("findSigmaSet", "FindSigma", 0, 0);
  
  pc.defineMutex("FindSigma", "Minos");

  pc.process(fFitCmdList);
  if (!pc.ok(kTRUE)) {
    throw std::runtime_error("unable to parse cmd list!");
  }

  fFindSigma = pc.getInt("findSigma");
  fFindSigmaN = fabs(pc.getDouble("findSigmaN"));
  fFindSigmaIter = pc.getInt("findSigmaIter");
  fFindSigmaSet = static_cast<RooArgSet *>(pc.getObject("findSigmaSet"));
  fFindSigmaPrecision = pc.getDouble("findSigmaPrecision");
  
  
  fOptConst = pc.getInt("optConst");
  fOffset = pc.getInt("doOffset");
  fVerbose = pc.getInt("verbose");
  fSave = pc.getInt("doSave");
  fTimer = pc.getInt("doTimer");
  fPrintLevel = pc.getInt("plevel");
  fDefaultStrategy = pc.getInt("strat");
  fHesse = pc.getInt("hesse");
  fMinimize = pc.getInt("minimize");  
  fMinos = pc.getInt("minos");
  fNumee = pc.getInt("numee");
  fDoEEWall = pc.getInt("doEEWall");
  fRetry = pc.getInt("retry");
  fMaxCalls = pc.getInt("maxcalls");
  fMaxIterations = pc.getInt("maxiterations");  
  fEigen = pc.getInt("eigen");
  fReuseMinimizer = pc.getInt("reminim");
  fReuseNLL = pc.getInt("renll");
  fChi2 = pc.getInt("usechi2");  
  fEps = pc.getDouble("eps");
  RooArgSet* minosSet = NULL;
  minosSet = static_cast<RooArgSet *>(pc.getObject("minosObjSet"));
  minosSet = pc.getSet("minosSet",minosSet);
  if(minosSet) fMinosSet = minosSet;
  fCondSet = static_cast<RooArgSet *>(pc.getObject("condSet"));
  fMinimizerType = std::string(pc.getString("mintype", "Minuit2"));
  fMinimizerAlgo = std::string(pc.getString("minalg", "Migrad"));

  return fFitCmdList.GetSize();
}
template
int RooFitUtils::ExtendedMinimizer::parseFitConfig<RooLinkedList>(const RooLinkedList &cmdList);

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

  try {
    int strategy = fDefaultStrategy;
    int retry = fRetry;
    int status = -1;

    RooArgSet* args = fNll->getVariables();
    int ndim = ::countFloatParams(args);
    delete args;

    fMinimizer->setMaxFunctionCalls(fMaxCalls);
    fMinimizer->setMaxIterations(fMaxIterations);
    fMinimizer->setPrintLevel(fPrintLevel);
    if(ndim > 0){
      fMinimizer->optimizeConst(fOptConst);
    }
    fMinimizer->setMinimizerType(fMinimizerType.c_str());
    fMinimizer->setEvalErrorWall(fDoEEWall);
    fMinimizer->setOffsetting(fOffset);
    fMinimizer->setPrintEvalErrors(fNumee);
    fMinimizer->setVerbose(fVerbose);
    fMinimizer->setProfile(fTimer);
    fMinimizer->setStrategy(fDefaultStrategy);
    fMinimizer->setEps(fEps);

    double nllval = 0.;
    RooFitResult* fitres = 0;
    bool minimize = fMinimize;

    
    //fNll->printTree(std::cout);
    bool ok = false;    
    while (!ok) {
      if(fitres) delete fitres;
      fitres = 0;
      this->SetNllDirty();      
      fMinimizer->setStrategy(strategy);
      if(minimize){
        // the following line is nothing but
        // int ndim = fMinimizer->getNPar();
        if(ndim > 0){
          std::cout << "ExtendedMinimizer::robustMinimize(" << fName
                    << "): starting minimization with strategy "
                    << strategy << std::endl;
          status = fMinimizer->minimize(fMinimizerType.c_str(), fMinimizerAlgo.c_str());
          if(status==0 || status==1) ok=true;
        } else {
          std::cout << "ExtendedMinimizer::robustMinimize(" << fName
                    << "): skipping minimization, no free parameters given!" << std::endl;
          minimize = false;
          ok = true;
          status = 0;
        }
      } else {
        std::cout << "ExtendedMinimizer::robustMinimize(" << fName
                  << "): skipping minimization" << std::endl;
        ok = true;
        status = 0;
      }      

      if(ok){
        std::string name = Form("fitresult_%s_%s", GetName(), fData->GetName());
        std::string title = Form("Result of fit of p.d.f. %s to dataset %s",
                                 GetName(), fData->GetName());
        if(minimize){
          fitres = fMinimizer->save(name.c_str(), title.c_str());
        }
        
#ifdef USE_ROOFITRESULT_NLL
        nllval = fitres ? fitres->minNll() : nan;
        if(!minimize){
          this->SetNllDirty();          
          nllval = fNll->getVal();
        }
#else
        this->SetNllDirty();
        nllval = fNll->getVal();
#endif
        std::cout
          << "ExtendedMinimizer::robustMinimize(" << fName
          << ") fit succeeded with status " << status << ", NLL=" << nllval;
        if(fNll->isOffsetting()) std::cout << " (offset=" << fNll->offset() << ")";
        std::cout << std::endl;
      }
      
      
      if (!ok || std::isnan(nllval) || std::isinf(nllval)){
        if(strategy < 2 && retry > 0) {
          strategy++;
          std::cout
            << "ExtendedMinimizer::robustMinimize(" << fName
            << ") fit failed with status " << status
            << ". Retrying with strategy " << strategy << std::endl;
          retry--;
        } else {
          std::cout
            << "ExtendedMinimizer::robustMinimize(" << fName
            << ") fit failed with status " << status << ", giving up." << std::endl;
          break;
        }
      } else {
        break;
      }
    }
    
    Result::Minimization mini;
    mini.status = status;
    mini.strategy = strategy;
    mini.ndim = ndim;
    mini.constOpt = fOptConst;
    mini.config = fMinimizer->fitter()->Config();
    mini.fit = fitres;
    
    if (ok && fHesse) {
      runHesse(mini);
    }
    
    //    if(ndim != ::countFloatParams(fMinimizer)){
    //      //      throw std::runtime_error(TString::Format("dimensionality inconsistency detected between minimizer (ndim=%d) and Nll (ndim=%d)!",::countFloatParams(fMinimizer),ndim).Data());
    //    }

    if (!mini.ok()) {
      std::cout << "ExtendedMinimizer::robustMinimize(" << fName
                            << ") fit failed with status " << status << std::endl;
    } else {
      mini.nll = nllval;
      mini.nllOffset = fNll->offset();
    }

    return mini;

  } catch (std::string& s){
    throw std::runtime_error(s);
  }

}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::initialize() {
  // apply all pre-minimization settings
  if (fCondSet && fNll) {
    setVals(*fNll->getVariables(), fCondSet, true);
  }
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result RooFitUtils::ExtendedMinimizer::minimize(const RooLinkedList &cmdList) {
  // Minimize function  adopted, simplified and extended from RooAbsPdf::fitTo()
  parseFitConfig(cmdList);

  return minimize();
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result RooFitUtils::ExtendedMinimizer::minimize() {
  // Minimize function  adopted, simplified and extended from RooAbsPdf::fitTo()
  initialize();

  setup();

  this->fResult = run();

  return *this->fResult;
}

namespace {
  TMatrixDSym reduce(const TMatrixDSym& mat){
    std::vector<int> keepRows;
    for(int i=0; i<mat.GetNcols(); ++i){
      if(mat(i,i) != 0) keepRows.push_back(i);
    }
    TMatrixDSym reduced(keepRows.size());
    for(size_t i=0; i<keepRows.size(); ++i){
      for(size_t j=0; j<keepRows.size(); ++j){
        reduced(i,j) = mat(keepRows[i],keepRows[j]);
      }
    }
    return reduced;
  }
}

namespace {
  bool find(const std::vector<RooFitUtils::ExtendedMinimizer::Result::Parameter>& pars,const char* name){
    for(const auto& p:pars){
      if(p.name == name) return true;
    }
    return false;
  }
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::SetNllDirty(){
  // set the value of the Nll to "dirty" and thus enforce the recalculation on next evaluation
  if(fNll){
    RooArgSet params(*fNll->getParameters((RooArgSet*)0));
    for(auto p:params){
      p->setValueDirty();
    }
    fNll->setValueDirty();
  }
}

// ____________________________________________________________________________|__________

RooFitUtils::ExtendedMinimizer::Result *RooFitUtils::ExtendedMinimizer::run() {
  // run the actual minimization
  Result *r = new Result();
  r->min = robustMinimize();

  if(!r->min.ok()){
    return r;
  }
  
  if(r->min.fit  && r->min.ndim > 0 && r->min.ndim != r->min.fit->floatParsFinal().getSize()){
    throw std::runtime_error("dimensionality inconsistency detected between minimizer and final floating parameter list!");
  }

  if(r->min.fit && !r->min.hesse){
    // if we don't have a hesse matrix yet, evaluate errors with Hesse
    r->min.covqual = r->min.fit->covQual();
    //if (covqual != -1) {
    /*std::cout*/ std::cout << "ExtendedMinimizer::minimize(" << fName
                                        << "): Covariance quality is " << r->min.covqual
                                        << ". " << std::endl;
    
    // Obtain Hessian matrix either from patched Minuit or after inversion
    // TMatrixDSym G = Minuit2::MnHesse::lastHessian();
    Double_t determ = 0;
    std::cout << "ExtendedMinimizer::minimize(" << fName
              << "): attempting to invert covariance matrix... "
              << std::endl;
    if (r->min.covqual != 0 && r->min.ndim>1 && r->min.fit) {
      const TMatrixDSym origCov = r->min.fit->covarianceMatrix();
      
      if(origCov.GetNcols() != r->min.ndim){
        throw std::runtime_error(TString::Format("inconsistency detected: correlation matrix size %d inconsistent with number float parameters %d!",origCov.GetNcols(),r->min.ndim).Data());
      }
      
      TMatrixDSym origG(::reduce(origCov));
      
      if(origG.GetNcols() != origCov.GetNcols()){
        std::cout   << "ExtendedMinimizer::minimize(" << fName
                    << "): removed " << origCov.GetNcols() - origG.GetNcols() << " dimensions from covariance matrix due to degeneration! " << std::endl;
      }
      
      TMatrixDSym G = origG.Invert(&determ);
      if (determ == 0 || std::isnan(determ)) {
        std::cout   << "ExtendedMinimizer::minimize(" << fName
                    << "): covariance matrix is singular! " << std::endl;
      } else {
        r->min.hesse = new TMatrixDSym(G);
        // Eigenvalue and eigenvector analysis
      }
    }
  }

  if(r->min.hesse && fEigen) {
    std::cout << "ExtendedMinimizer::minimize(" << fName
                          << "): starting eigen analysis... " << std::endl;
    r->eigen = eigenAnalysis(*(r->min.hesse));
  }
    
  
  if (r->min.ok()) {
    // Evaluate errors with Minos
    if (fMinos) {
      std::cout << "ExtendedMinimizer::minimize(" << fName
                            << "): Running Minos" << std::endl;
      if (fMinosSet) {
        fMinimizer->minos(*fMinosSet);
      } else {
        fMinimizer->minos();
      }

      delete r->min.fit;
      std::string name = Form("fitresult_%s_%s", GetName(), fData->GetName());
      std::string title = Form("Result of fit of p.d.f. %s to dataset %s",
                               GetName(), fData->GetName());
      r->min.fit = fMinimizer->save(name.c_str(), title.c_str());
    }
  }
  

  if (fFindSigma) {
    std::cout << "ExtendedMinimizer::minimize(" << fName
                          << "): Running findSigma" << std::endl;
    findSigma(r);
  }


  RooArgSet *vars = fNll->getVariables();
  for (auto arg:*vars){
    RooRealVar* v = dynamic_cast<RooRealVar*>(arg);
    if(!::find(r->min.parameters,v->GetName()) && !v->isConstant()){
      Result::Parameter poi(v->GetName(), v->getVal(), v->getErrorHi(),
                            v->getErrorLo());
      r->min.parameters.push_back(poi);
    }
  }

  if (fCondSet) {
    std::cout << "Editing conditional set" << std::endl;
    RooArgSet *attachedSet = fNll->getVariables();
    for (auto arg:*fCondSet){
      RooRealVar* v = dynamic_cast<RooRealVar *>(arg);
      RooRealVar *var = dynamic_cast<RooRealVar *>(attachedSet->find(v->GetName()));
      if(var){
        var->setVal(v->getVal());
        var->setConstant(v->isConstant());
      }
    }
  }
  // Return fit result
  if (r->min.ndim > 0) {
    if(r->min.fit->correlationMatrix().GetNcols() < r->min.fit->floatParsFinal().getSize()) throw std::runtime_error(TString::Format("fit result size %d is inconsistent with correlation matrix size %d!",r->min.fit->floatParsFinal().getSize(),r->min.fit->correlationMatrix().GetNcols()).Data());
  }
  
  if(!fReuseMinimizer) {
    if(fMinimizer)
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
  if(!attachedSet) return;
  
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

  std::vector<RooRealVar *> extraparams;
  std::vector<std::string> extraparnames;  
  if(fModel && fModel->GetParametersOfInterest()) {
    for (auto par : *fModel->GetParametersOfInterest()) {
      RooRealVar *v = dynamic_cast<RooRealVar *>(attachedSet->find(par->GetName()));
      if (!v) {
        throw std::runtime_error(
                                 TString::Format("unknown parameter name: %s", par->GetName()).Data());
      }
    }
  }

  bool hesse = fHesse;
  ExtendedMinimizer::Result::Scan scan(RooFitUtils::concat(parnames),parnames,extraparnames);

  for (const auto &point : points) {
    if (point.size() != parnames.size()) {
      throw std::runtime_error("inconsistent vector lengths in scan!");
    }
    for (size_t i = 0; i < params.size(); ++i) {
      params[i]->setVal(point[i]);
      params[i]->setConstant(true);
    }
    fHesse = false;
    auto min = this->robustMinimize();
    std::vector<double> vals(params.size());
    for (size_t i = 0; i < params.size(); ++i) {
      vals[i] = params[i]->getVal();
    }
    std::vector<double> extravals(extraparams.size());
    for (size_t i = 0; i < extraparams.size(); ++i) {
      extravals[i] = extraparams[i]->getVal();
    }
    scan.add(vals, min.status, min.nll, extravals);
  }

  if (scan.nllValues.size() > 0)
    r->scans.push_back(scan);
  fHesse = hesse;
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::findSigma() {
  // run an iterative algorithm to find the 1-sigma-band
  if (!fNll) this->setup();
  findSigma(this->getResult(true));
}

// ____________________________________________________________________________|__________

void RooFitUtils::ExtendedMinimizer::findSigma(ExtendedMinimizer::Result* r) {
  // run an iterative algorithm to find the 1-sigma-band
  if (!fNll) {
    throw std::runtime_error("invalid Nll!");
  }
  if (!fFindSigmaSet || fFindSigmaSet->getSize() == 0) {
    RooArgSet *attachedSet = fNll->getVariables();
    RooStats::RemoveConstantParameters(attachedSet);
    findSigma(r, *attachedSet);
  } else {
    findSigma(r, *fFindSigmaSet);
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

void RooFitUtils::ExtendedMinimizer::findSigma(Result *r, const RooAbsCollection &scanSet) {
  // run an iterative algorithm to find the 1-sigma-band
  int optConst = fOptConst;
#ifdef FINDSIGMA_NO_CONSTOPT      
  fOptConst = 0;
#endif
  Result::Minimization min(r->min);
  if (!min.ok()) {
    std::cout << "ExtendedMinimizer::findSigma(): no previous "
                             "minimization detected, geneating new minimum "
                          << std::endl;
    min = robustMinimize();
    if(!min.ok()){
      throw std::runtime_error("unable to findSigma without minimum!");
    }
#ifdef FINDSIGMA_NO_CONSTOPT    
  } else if (r->min.constOpt != 0){
    std::cout << "ExtendedMinimizer::findSigma(): previous minimization was using const optimization, generating new minimum "  << std::endl;
    min = robustMinimize();
    if(!min.ok()){
      throw std::runtime_error("unable to findSigma without minimum!");
    } else if (!RooFitUtils::compare(min.fit,r->min.fit)){
      throw std::runtime_error("rerunning minimization yielded inconsistent minimum!");      
    }
#endif
  }
  for (const auto& it:scanSet){
    RooRealVar *v = dynamic_cast<RooRealVar *>(it);
    if (!v) {
      throw std::runtime_error("invalid variable!");
    }
    std::cout << "ExtendedMinimizer::findSigma(): starting scan of "
                          << v->GetName() << std::endl;
    
    RooArgSet vars(*fPdf->getVariables());
    RooStats::RemoveConstantParameters(&vars);
    vars.add(*v, kTRUE);
    RooArgSet *snap = dynamic_cast<RooArgSet *>(vars.snapshot());
    if (!snap) {
      throw std::runtime_error("invalid snapshot!");
    }

    Double_t val = v->getVal();
    Double_t err = fFindSigmaN * v->getError();
    if(err < 1e-9) err = 0.1*val;

    setVals(vars, snap, false);
    Double_t shi = findSigma(r,min, val + err, val, v, fFindSigmaN);
    setVals(vars, snap, false);
    Double_t slo = findSigma(r,min, val - err, val, v, -fFindSigmaN);
    setVals(vars, snap, false);
    
    Result::Parameter poi(v->GetName(), val, shi, slo);
    r->min.parameters.push_back(poi);

    v->setAsymError(std::isnan(slo) ? val+err : slo, std::isnan(shi) ? val-err : shi);
    std::cout << "ExtendedMinimizer::minimize(" << fName << ") "
                          << std::endl;
    delete snap;
  }
  fOptConst = optConst;
}

// ____________________________________________________________________________|__________

double RooFitUtils::ExtendedMinimizer::findSigma(RooFitUtils::ExtendedMinimizer::Result* result,
                                                   const RooFitUtils::ExtendedMinimizer::Result::Minimization& min,
                                                   const double guessval, const double val_mle, RooRealVar *par, const double nsigma) {
  // Find the value of sigma evaluated at a specified nsigma, assuming NLL
  // -2logL is roughly parabolic in par.
  // The precision is specified as a fraction of the error, or based on the
  // Minuit default tolerance.
  // Based on
  // https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/HiggsPhys/HSG3/WWDileptonAnalysisCode/HWWStatisticsCode/trunk/macros/findSigma.C
  // by Aaron Armbruster <aaron.james.armbruster@cern.ch> and adopted by Tim
  // Adye <T.J.Adye@rl.ac.uk>.

  double precision = fFindSigmaPrecision;
  int maxiter = fFindSigmaIter;  
  const double fitTol = fEps;
  bool isConst(par->isConstant());

  TRandom rnd;
  
  double val_guess = useLimits(par, guessval);
  int direction = nsigma >= 0.0 ? +1 : -1;
  int nDamping = 1;
  double damping_factor = 1.0;
  ExtendedMinimizer::ValueMap guess_to_corr;
  double tmu = inf;

  if (precision <= 0.0) {
    // RooFit default tolerance is 1.0.
    double eps = 0.001 * (fitTol > 0.0 ? fitTol : 1.0);
    precision = 5.0 * eps / (nsigma * nsigma);
  }

  bool hesse = fHesse;
  int optConst = fOptConst;
  fHesse = false;
#ifdef FINDSIGMA_NO_CONSTOPT      
  fOptConst = 0;
#endif
  ExtendedMinimizer::Result::Scan values(TString::Format("findSigma_%s_%g",par->GetName(),val_guess).Data(),{par->GetName()});
  const double nllmin = min.nll;
  values.add({val_mle}, min.status, nllmin);
  int iter = 0;
  while(iter < maxiter){
    iter++;
    std::cout
        << "ExtendedMinimizer::findSigma(" << fName << ") "
        << Form("Parameter %s %+gsigma iteration %d: start %g (MLE%+g)",
                par->GetName(), nsigma, iter + 1, val_guess,
                val_guess - val_mle)
        << std::endl;
    double lastval = par->getVal();
    double val_pre = val_guess;
    par->setVal(val_pre);
    par->setConstant(true);
    Result::Minimization mini = robustMinimize();

    double nll = mini.nll;
    double tmu_new = 2.0 * (nll - nllmin);
    if(!mini.ok()){
      //      || (fabs(tmu_new - 1) > fabs(tmu - 1)+0.1)){
      val_guess = rnd.Gaus(lastval,0.1*useLimits(par, guessval));
      std::cout << "ExtendedMinimizer::findSigma(" << fName <<") failed to converge at " << par->getVal() << ", tmu = " << tmu << " => " << tmu_new << ", backtracking to " << val_guess << std::endl;
      continue;
    } else {
      tmu = tmu_new;
      values.add({val_pre}, mini.status, nll);
      double sigma_guess = fabs(val_guess - val_mle);
      if (tmu > 0.01)   sigma_guess /= sqrt(tmu);
      else              sigma_guess *=      10.0; // protect against tmu<=0, and also don't move too far

      double corr = damping_factor * (val_pre - val_mle - nsigma * sigma_guess);
      for (const auto& iguess:guess_to_corr){
        if (fabs(iguess.first - val_pre) < direction * val_pre * 0.02) {
          damping_factor *= 0.8;
          std::cout
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
  
      std::cout << "ExtendedMinimizer::findSigma(" << fName << ") "
                            << Form("%s %.3f (MLE%+.3f) -> %.3f (MLE%+.3f), "
                                    "change %+.3f, precision %.3f, -2lnL %.4f, "
                                    "sigma(guess) %.3f",
                                    par->GetName(), val_pre, val_pre - val_mle,
                                    val_guess, val_guess - val_mle, delta,
                                    relprecision, tmu, sigma_guess)
                            << std::endl;
  
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") NLL:                 " << nll << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") delta(NLL):          " << nll - nllmin
                            << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") nsigma*sigma(pre):   "
                            << fabs(val_pre - val_mle) << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") sigma(guess):        " << sigma_guess
                            << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") par(guess):          " << val_guess + corr
                            << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") best-fit val:        " << val_mle << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") tmu:                 " << tmu << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") Precision:           "
                            << direction * val_guess * precision << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") Correction:          " << -corr << std::endl;
      std::cout << "ExtendedMinimizer::findSigma(" << fName
                            << ") nsigma*sigma(guess): "
                            << fabs(val_guess - val_mle) << std::endl;
      if (fabs(delta) <= relprecision)
        break;
    }
  }
  bool ok = (!(iter >= maxiter));
  fHesse = hesse;
  fOptConst = optConst;  

  par->setConstant(isConst);
  if (!ok) {
    std::cerr << "findSigma failed after " << iter << " iterations"
              << std::endl;
    return TMath::QuietNaN();
  }

  result->scans.push_back(values);

  double err = val_guess - val_mle;
  std::cout
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

  for (auto it_poi:map_poi2nll){
    double nll = it_poi.second;
    if (nll == nll && fabs(nll) < pow(10, 20)) {
      x.push_back(it_poi.first);
      y.push_back(it_poi.second);
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
