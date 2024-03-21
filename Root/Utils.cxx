#include <chrono>
#include <iomanip>

#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTime.h"

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"

#include "RooFitUtils/Log.h"

#include "RooFitUtils/Utils.h"
namespace RooFitUtils {
bool RooStarMomentMorphFix = true;
bool RooMultiPdfFix = true;
}

namespace {
union MyFloat_t {
  MyFloat_t(float num = 0.0f) : f(num) {}
  // Portable extraction of components.
  bool Negative() const { return (i >> 31) != 0; }
  Int_t RawMantissa() const { return i & ((1 << 23) - 1); }
  Int_t RawExponent() const { return (i >> 23) & 0xFF; }

  Int_t i;
  float f;
};
}

#ifdef HAS_ROOSTARMOMENTMORPH
#include "RooStarMomentMorph.h"
#endif


#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
RooAbsCollection_IteratorHelper::RooAbsCollection_IteratorHelper(const RooAbsCollection& c, bool end) : itr(c.fwdIterator()), nextItem(end ? NULL : itr.next()) {}
RooAbsArg* RooAbsCollection_IteratorHelper::operator++(){ nextItem = itr.next(); return nextItem; }
bool RooAbsCollection_IteratorHelper::operator!=(const RooAbsCollection_IteratorHelper& other){ return this->nextItem != other.nextItem; }
bool RooAbsCollection_IteratorHelper::operator!=(const RooAbsArg* other){  return this->nextItem != other; }
RooAbsArg* RooAbsCollection_IteratorHelper::operator*(){ return nextItem; }
RooAbsCollection_IteratorHelper begin(const RooAbsCollection& c){
  return RooAbsCollection_IteratorHelper(c,false);
}
RooAbsCollection_IteratorHelper end(const RooAbsCollection& c){
  return RooAbsCollection_IteratorHelper(c,true);
}
#endif
RooLinkedList_IteratorHelper::RooLinkedList_IteratorHelper(const RooLinkedList& c, bool end) : itr(c.fwdIterator()), nextItem(end ? NULL : itr.next()) {}
RooAbsArg* RooLinkedList_IteratorHelper::operator++(){ nextItem = itr.next(); return nextItem; }
bool RooLinkedList_IteratorHelper::operator!=(const RooLinkedList_IteratorHelper& other){ return this->nextItem != other.nextItem; }
bool RooLinkedList_IteratorHelper::operator!=(const RooAbsArg* other){  return this->nextItem != other; }
RooAbsArg* RooLinkedList_IteratorHelper::operator*(){ return nextItem; }
RooLinkedList_IteratorHelper begin(const RooLinkedList& c){
  return RooLinkedList_IteratorHelper(c,false);
}
RooLinkedList_IteratorHelper end(const RooLinkedList& c){
  return RooLinkedList_IteratorHelper(c,true);
}
RooFIter_IteratorHelper::RooFIter_IteratorHelper(RooFIter& it, bool end) : itr(&it), nextItem(end ? NULL : itr->next()) {}
RooAbsArg* RooFIter_IteratorHelper::operator++(){ nextItem = itr->next(); return nextItem; }
bool RooFIter_IteratorHelper::operator!=(const RooFIter_IteratorHelper& other){ return this->nextItem != other.nextItem; }
bool RooFIter_IteratorHelper::operator!=(const RooAbsArg* other){  return this->nextItem != other; }
RooAbsArg* RooFIter_IteratorHelper::operator*(){ return nextItem; }
RooFIter_IteratorHelper begin(RooFIter& it){
  return RooFIter_IteratorHelper(it,false);
}
RooFIter_IteratorHelper end(RooFIter& it){
  return RooFIter_IteratorHelper(it,true);
}
TIterator_IteratorHelper::TIterator_IteratorHelper(TIterator* it, bool end) : itr(it), nextItem(end ? NULL : itr->Next()) {}
TObject* TIterator_IteratorHelper::operator++(){ nextItem = itr->Next(); return nextItem; }
bool TIterator_IteratorHelper::operator!=(const TIterator_IteratorHelper& other){ return this->nextItem != other.nextItem; }
bool TIterator_IteratorHelper::operator!=(const TObject* other){  return this->nextItem != other; }
TObject* TIterator_IteratorHelper::operator*(){ return nextItem; }
TIterator_IteratorHelper begin(TIterator* it){
  return TIterator_IteratorHelper(it,false);
}
TIterator_IteratorHelper end(TIterator* it){
  return TIterator_IteratorHelper(it,true);
}


// _____________________________________________________________________________

int RooFitUtils::fixRooStarMomentMorph(RooWorkspace *workspace) {
  int retval = 0;
#ifdef HAS_ROOSTARMOMENTMORPH
  RooFIter iter(workspace->components().fwdIterator());
  RooAbsArg *arg;
  while ((arg = iter.next())) {
    RooStarMomentMorph *pdf = dynamic_cast<RooStarMomentMorph *>(arg);
    if (pdf) {
      retval++;
      pdf->fixCache();
    }
  }
#else
  (void)workspace; // Remove unused parameter warning. We do nothing.
#endif
  return retval;
}
// ____________________________________________________________________________

void RooFitUtils::addArgSet(RooArgSet* args, const RooArgSet* addArgs){
  // convert const RooArgSet* to RooArgSet*
  if(addArgs) args->add(*addArgs);
}
// _____________________________________________________________________________

bool RooFitUtils::AlmostEqualUlpsAndAbs(float A, float B, float maxDiff,
                                        int maxUlpsDiff) {
  // Check if the numbers are really close -- needed  when comparing numbers
  // near zero.
  float absDiff = fabs(A - B);
  if (absDiff <= maxDiff)
    return true;

  MyFloat_t uA(A);
  MyFloat_t uB(B);

  // Different signs means they do not match.
  if (uA.Negative() != uB.Negative())
    return false;

  // Find the difference in ULPs.
  int ulpsDiff = abs(uA.i - uB.i);
  if (ulpsDiff <= maxUlpsDiff)
    return true;

  return false;
}

// _____________________________________________________________________________

void RooFitUtils::PrintResourcesUsed(const TTime &progStart) {
  // Print used resources
  // Courtesy of Tim Adye <T.J.Adye@rl.ac.uk>.
  ProcInfo_t info;
  if (gSystem->GetProcInfo(&info) < 0)
    return;
  Long_t cput = TMath::CeilNint(info.fCpuUser);
  Long_t wall =
      Long64_t(gSystem->Now() - progStart + TTime(500)) / Long64_t(1000);
  LOG(RooFitUtils::logINFO)
      << Form("resources used: cput=%02ld:%02ld:%02ld, mem=%ldkb, vmem=%ldkb, "
              "walltime=%02ld:%02ld:%02ld",
              cput / 3600, (cput / 60) % 60, cput % 60, info.fMemResident,
              info.fMemVirtual, wall / 3600, (wall / 60) % 60, wall % 60);
}

// _____________________________________________________________________________

std::vector<std::string> RooFitUtils::parseString(const std::string &str,
                                                  const std::string &sep) {
  // Split strings according to separator
  std::vector<std::string> parsed;
  int pos = 0;
  bool first = true;
  if (str.size() == 0)
    return parsed;
  if (str.find(sep) == std::string::npos) {
    parsed.push_back(str);
    return parsed;
  }

  while (true) {
    int newPos = str.find(sep, pos);
    if (str.find(sep, pos) == std::string::npos) {
      if (!first)
        parsed.push_back(str.substr(pos, newPos - pos));
      break;
    }

    std::string sub = str.substr(pos, newPos - pos);
    parsed.push_back(sub);
    pos = newPos + 1;
    first = false;
  }

  return parsed;
}

// _____________________________________________________________________________

void RooFitUtils::FindUniqueProdComponents(RooProdPdf *Pdf,
                                           RooArgSet &Components) {
  // Split a RooProdPdf into its components
  static int counter = 0;
  counter++;

  if (counter > 50) {
    LOG(RooFitUtils::logERROR)
        << "FindUniqueProdComponents detected infinite loop. Please check.";
    exit(1);
  }

  RooArgList pdfList = Pdf->pdfList();
  if (pdfList.getSize() == 1) {
    LOG(RooFitUtils::logINFO) << "FindUniqueProdComponents "
                              << pdfList.at(0)->GetName() << " is fundamental.";
    Components.add(pdfList);
  } else {
    for(auto* nextArg : pdfList){
      RooProdPdf *Pdf = (RooProdPdf *)nextArg;
      if (std::string(Pdf->ClassName()) != "RooProdPdf") {
        LOG(RooFitUtils::logINFO)
            << "FindUniqueProdComponents " << Pdf->GetName()
            << " is no RooProdPdf. Adding it.";
        Components.add(*Pdf);
        continue;
      }
      FindUniqueProdComponents(Pdf, Components);
    }
  }
  counter = 0;
}

// _____________________________________________________________________________

bool RooFitUtils::ensureDirectory(const TString &path) {
  // ensure that the directory with the given path exists
  // check if directory <path> exists
  Long_t flags = 0;
  gSystem->GetPathInfo(path.Data(), (Long_t *)0, (Long_t *)0, &flags,
                       (Long_t *)0);
  if (flags & 2) {
    // directory exists
    return true;
  }
  // create directory
  if (0 == gSystem->mkdir(path.Data(), true))
    return true;
  else
    return false;
}

// _____________________________________________________________________________

bool RooFitUtils::ensureDirectoryForFile(const TString &file) {
  // ensure that the directory for the given file exists
  Ssiz_t pos = file.Last('/');
  if (pos == kNPOS)
    return false;
  return ensureDirectory(file(0, pos));
}

// _____________________________________________________________________________

void RooFitUtils::PrintTable(std::string *firstCol, std::string **matrix,
                             std::string **matrixErr, std::string *header,
                             int nrRows, int nrCols, int nSigFig,
                             std::ostream &ost, std::string indent,
                             std::string delim, std::string ending) {
  // Helper function to print latex tables nicely
  int **lmatrix = new int *[nrRows];
  int **lmatrix_pm = new int *[nrRows];

  for (int i = 0; i < nrRows; i++) {
    lmatrix[i] = new int[nrCols];
    lmatrix_pm[i] = new int[nrCols];
  }

  int maxLfirst = 0;
  for (int i = 0; i < nrRows + 1; i++) {
    maxLfirst = int(std::max(double(maxLfirst), double(firstCol[i].size())));
  }

  std::string pm = " $\\pm$ ";

  int *maxLength = new int[nrCols];
  int *maxLength_pm = new int[nrCols];

  for (int j = 0; j < nrCols; j++) {
    if (header)
      maxLength[j] = header[j].size();
    else
      maxLength[j] = 0;
    maxLength_pm[j] = 0;
  }

  for (int i = 0; i < nrRows; i++) {
    for (int j = 0; j < nrCols; j++) {
      std::stringstream str;
      str << std::setprecision(nSigFig);
      str << matrix[i][j];

      lmatrix[i][j] = str.str().size();
      maxLength[j] = int(std::max(double(maxLength[j]), double(lmatrix[i][j])));

      std::stringstream str2;
      str2 << std::setprecision(nSigFig);
      if (/*j != 0 && */ matrixErr)
        str2 << pm << matrixErr[i][j];

      lmatrix_pm[i][j] = str2.str().size();
      maxLength_pm[j] =
          int(std::max(double(maxLength_pm[j]), double(lmatrix_pm[i][j])));
    }
  }

  if (header) {
    ost << indent;
    ost << firstCol[0];
    for (int j = firstCol[0].size(); j < maxLfirst; j++)
      ost << " ";
    ost << " & ";
    for (int i = 0; i < nrCols; i++) {
      ost << header[i];

      for (int k = header[i].size(); k < maxLength[i] + maxLength_pm[i]; k++) {
        ost << " ";
      }
      if (i < nrCols - 1) {
        if (delim != "")
          ost << delim;
      } else {
        if (ending != "")
          ost << ending;
        ost << "\n" << indent << "\\hline\n";
      }
    }
  }

  for (int i = 0; i < nrRows; i++) {
    ost << indent;
    ost << firstCol[i + 1];
    for (int j = firstCol[i + 1].size(); j < maxLfirst; j++)
      ost << " ";
    ost << " & ";
    for (int j = 0; j < nrCols; j++) {
      ost << std::setprecision(nSigFig);
      ost << matrix[i][j];
      for (int k = lmatrix[i][j]; k < maxLength[j]; k++) {
        ost << " ";
      }

      if (/*j != 0 && */ matrixErr) {
        ost << pm << matrixErr[i][j];
        for (int k = lmatrix_pm[i][j]; k < maxLength_pm[j]; k++) {
          ost << " ";
        }
      }

      if (j < nrCols - 1) {
        if (delim != "")
          ost << delim;
      } else {
        if (ending != "")
          ost << ending;
        ost << "\n";
      }
    }
  }

  delete[] maxLength;
  delete[] maxLength_pm;

  for (int i = 0; i < nrRows; i++) {
    delete[] lmatrix[i];
    delete[] lmatrix_pm[i];
  }

  delete[] lmatrix;
  delete[] lmatrix_pm;
}

//__________________________________________________________________________________|___________

bool RooFitUtils::matches(const std::string& text, const std::string& pattern) {
  // Performs a string match between the input string <text> and the string pattern
  // <pattern> and returns true in case of a match and false otherwise. The string
  // pattern may use wildcards "?" (matching exactly one character) and "*"
  // (matching any string sequence).
  //
  // Examples:
  //
  // - matches("hello", "h?llo") returns true
  // - matches("hallo", "h?llo") returns true
  // - matches("hello", "h*") returns true
  // - matches("hello", "h*a") returns false
 
  if (text.size() > 0 && pattern.size() > 0) {
    // direct and "?" match ?
    if (text[0] == pattern[0] || pattern[0] == '?') {
      // TODO: avoid recursion
      return matches(text.substr(1, text.size()), pattern.substr(1, pattern.size()));
    }
    // "*" match ?
    if (pattern[0] == '*') {
      // eating leading "*" in pattern ...
      return matches(text, pattern.substr(1, pattern.size()))
        // ... or matching leading character in text ?
        || matches(text.substr(1, text.size()), pattern);
    }
    // no match
    return false;
  } else {
    // empty text and/or pattern
    return (text.size() == 0 && (pattern == "*" || pattern.size() == 0));
  }
}

//__________________________________________________________________________________|___________

namespace {
    //freeing members from RooWorkspace:
  template <typename RooWorkspaceTag>
  struct RooWorkspaceHackResult {
    typedef typename RooWorkspaceTag::type type;
    static type ptr;
  };
  
  template <typename RooWorkspaceTag>
  typename RooWorkspaceHackResult<RooWorkspaceTag>::type RooWorkspaceHackResult<RooWorkspaceTag>::ptr;
  
  template<typename RooWorkspaceTag, typename RooWorkspaceTag::type p>
  struct RooWorkspaceRob : RooWorkspaceHackResult<RooWorkspaceTag> {
    struct RooWorkspaceFiller {
      RooWorkspaceFiller() {RooWorkspaceHackResult<RooWorkspaceTag>::ptr = p;}
    };
    static RooWorkspaceFiller rooworkspacefiller_obj;
  };
  
  template<typename RooWorkspaceTag, typename RooWorkspaceTag::type p>
  typename RooWorkspaceRob<RooWorkspaceTag, p>::RooWorkspaceFiller RooWorkspaceRob<RooWorkspaceTag, p>::rooworkspacefiller_obj;
  
  //now expose some members of RooWorkspace that we need to access
  struct RooWorkspace_snapshots { typedef RooLinkedList(RooWorkspace::*type); };
  template class RooWorkspaceRob<RooWorkspace_snapshots, &RooWorkspace::_snapshots>;
}


namespace {
  template<class listT, class stringT> void getParameterNames(const listT* l,std::vector<stringT>& names){
    // extract the parameter names from a list
    if(!l) return;
    for(auto* obj : *l){
      names.push_back(obj->GetName());
    }
  }
  void getArgs(RooWorkspace* ws, const std::vector<TString> names, RooArgSet& args){
    for(const auto& p:names){
      RooRealVar* v = ws->var(p.Data());
      if(v){
        args.add(*v);
      }
    }
  }
}



RooWorkspace* RooFitUtils::makeCleanWorkspace(RooWorkspace* oldWS, const char* newName, bool copySnapshots, const char* mcname, const char* newmcname){
  // clone a workspace, copying all needed components and discarding all others
	
  // butcher the old workspace
  auto objects = oldWS->allGenericObjects();
  RooStats::ModelConfig* oldMC = dynamic_cast<RooStats::ModelConfig*>(oldWS->obj(mcname));
  auto data = oldWS->allData();
  for(auto it:objects){
    if(!oldMC){
      oldMC = dynamic_cast<RooStats::ModelConfig*>(it);
    }
  }
  if(!oldMC) throw std::runtime_error("unable to retrieve ModelConfig");
  
  RooAbsPdf* origPdf = oldMC->GetPdf();
  
  
  // butcher the old modelconfig
  std::vector<TString> poilist;
  std::vector<TString> nplist;
  std::vector<TString> obslist;
  std::vector<TString> globobslist;
  RooAbsPdf* pdf = NULL;
  if(oldMC){
    pdf = oldMC->GetPdf();
    ::getParameterNames(oldMC->GetParametersOfInterest(),poilist);
    ::getParameterNames(oldMC->GetNuisanceParameters(),nplist);
    ::getParameterNames(oldMC->GetObservables(),obslist);
    ::getParameterNames(oldMC->GetGlobalObservables(),globobslist);
  } 
  if(!pdf){
    if(origPdf) pdf=origPdf;
  }
  if(!pdf){
    return NULL;
  }

  // create them anew
  RooWorkspace* newWS = new RooWorkspace(newName ? newName : oldWS->GetName());
  newWS->autoImportClassCode(true);
  RooStats::ModelConfig* newMC = new RooStats::ModelConfig(newmcname, newWS);

  // Copy snapshots
  // Fancy ways to avoid public-private hack used in the following, simplified version in comments above the respective lines
  if(copySnapshots){
    RooFIter itr( ((*oldWS).*::RooWorkspaceHackResult<RooWorkspace_snapshots>::ptr).fwdIterator());
    RooArgSet* snap;
    while((snap = (RooArgSet*)itr.next())){
      RooArgSet* snapClone = (RooArgSet*) snap->snapshot() ;
      snapClone->setName(snap->GetName()) ;
      //newWS->_snapshots.Add(snapClone) ;
      ((*newWS).*::RooWorkspaceHackResult<RooWorkspace_snapshots>::ptr).Add(snapClone) ;
    }
  }
  
  newWS->import(*pdf, RooFit::RecycleConflictNodes());
  RooAbsPdf* newPdf = newWS->pdf(pdf->GetName());
  newMC->SetPdf(*newPdf);
  
  for(auto d:data){
    newWS->import(*d);
  }
  
  RooArgSet poiset; ::getArgs(newWS,poilist,poiset);
  RooArgSet npset; ::getArgs(newWS,nplist,npset);
  RooArgSet obsset; ::getArgs(newWS,obslist,obsset);
  RooArgSet globobsset; ::getArgs(newWS,globobslist,globobsset);
  
  newMC->SetParametersOfInterest(poiset);
  newMC->SetNuisanceParameters  (npset);
  newMC->SetObservables         (obsset);
  newMC->SetGlobalObservables   (globobsset);
  newWS->import(*newMC);
  
  return newWS;
}

std::string RooFitUtils::concat(const std::vector<std::string>& text, const std::string& joint){
  if(text.size()==0) return "";
  std::string retval = text[0];
  for(size_t i=1; i<text.size(); ++i){
    retval += joint;
    retval += text[i];
  }
  return retval;
}
TMatrixDSym RooFitUtils::convertToCorrelationMatrix(const TMatrixDSym& _VM){
  // stolen from RooFitResult
  TMatrixDSym _CM(_VM);
  for (Int_t i=0 ; i<_CM.GetNrows() ; i++) {
    for (Int_t j=0 ; j<_CM.GetNcols() ; j++) {
      if (i!=j) {
        (_CM)(i,j) = (_CM)(i,j) / sqrt(_CM(i,i)*_CM(j,j)) ;
      }
    }
  }
  for (Int_t i=0 ; i<_CM.GetNrows() ; i++) {
    (_CM)(i,i) = 1.0 ;
  }
  return _CM;
}


bool RooFitUtils::compare(const RooFitResult* r1, const RooFitResult* r2,double tol_rel) {
  if(!r1 || !r2) return false;
  for(auto ai:r1->floatParsFinal()){
    RooRealVar* pi = dynamic_cast<RooRealVar*>(ai);
    if(!pi) continue;
    double re1 = pi->getError()/pi->getVal();
    for(auto aj:r2->floatParsFinal()){
      RooRealVar* pj = dynamic_cast<RooRealVar*>(aj);
      if(!pj) continue;      
      if(strcmp(pi->GetName(),pj->GetName()) != 0) continue;
      double re2 = pj->getError()/pj->getVal();
      double tol = std::min(tol_rel,std::max(re1,re2));
      if(!TMath::AreEqualRel(pi->getVal(),pj->getVal(),tol)) return false;
    }
  }
  exit(0);
  return true;
}

#if ROOT_VERSION_CODE < ROOT_VERSION(6,10,0)
TIterator* clients(const RooAbsArg* arg) { return arg->clientIterator(); }
TIterator* servers(const RooAbsArg* arg) { return arg->serverIterator(); }
#else
const RooAbsArg::RefCountList_t & clients(const RooAbsArg* arg) { return arg->clients(); }
const RooAbsArg::RefCountList_t & servers(const RooAbsArg* arg) { return arg->servers(); }
#endif
