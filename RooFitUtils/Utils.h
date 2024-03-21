// this file looks like plain C, but it's actually -*- c++ -*-
#ifndef _UTILS_
#define _UTILS_

#include <string>
#include <vector>

#include "TTime.h"
#include "TMatrixDSym.h"

class RooProdPdf;
class RooWorkspace;

namespace RooFitUtils {
  extern bool RooStarMomentMorphFix;
  extern bool RooMultiPdfFix;
  
  int fixRooStarMomentMorph(RooWorkspace *workspace);

  void addArgSet(RooArgSet* args, const RooArgSet* addArgs);

  bool AlmostEqualUlpsAndAbs(float A, float B, float maxDiff, int maxUlpsDiff);
  void PrintResourcesUsed(const TTime &progStart);
  std::vector<std::string> parseString(const std::string &str,
                                       const std::string &sep);
  void FindUniqueProdComponents(RooProdPdf *Pdf, RooArgSet &Components);
  
  bool ensureDirectory(const TString &path);
  bool ensureDirectoryForFile(const TString &file);
  RooWorkspace* makeCleanWorkspace(RooWorkspace* oldWS, const char* newName = NULL, bool copySnapshots=true, const char* mcname = "ModelConfig", const char* newmcname = "ModelConfig");
  
  void PrintTable(std::string *firstCol, std::string **matrix,
                  std::string **matrixErr, std::string *header, int nrRows,
                  int nrCols, int nSigFig, std::ostream &ost,
                  std::string indent = "", std::string delim = " & ",
                  std::string ending = " \\\\");
	
  bool matches(const std::string& text, const std::string& pattern);
  std::string concat(const std::vector<std::string>& text, const std::string& joint = ",");  
  TMatrixDSym convertToCorrelationMatrix(const TMatrixDSym& covmat);
  bool compare(const RooFitResult* r1, const RooFitResult* r2, double tol_rel=1e-6);
}

#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
struct RooAbsCollection_IteratorHelper { //EXCLUDE
  RooFIter itr;
  RooAbsArg* nextItem;
  RooAbsCollection_IteratorHelper(const RooAbsCollection& c, bool end);
  RooAbsArg* operator++();
  bool operator!=(const RooAbsCollection_IteratorHelper& other);
  bool operator!=(const RooAbsArg* other);
  RooAbsArg* operator*();
};
RooAbsCollection_IteratorHelper begin(const RooAbsCollection& c);
RooAbsCollection_IteratorHelper end(const RooAbsCollection& c);
#endif
struct RooLinkedList_IteratorHelper { //EXCLUDE
  RooFIter itr;
  RooAbsArg* nextItem;
  RooLinkedList_IteratorHelper(const RooLinkedList& c, bool end);
  RooAbsArg* operator++();
  bool operator!=(const RooLinkedList_IteratorHelper& other);
  bool operator!=(const RooAbsArg* other);
  RooAbsArg* operator*();
};
RooLinkedList_IteratorHelper begin(const RooLinkedList& c);
RooLinkedList_IteratorHelper end(const RooLinkedList& c);
struct RooFIter_IteratorHelper { //EXCLUDE
  RooFIter* itr;
  RooAbsArg* nextItem;
  RooFIter_IteratorHelper(RooFIter& it, bool end);
  RooAbsArg* operator++();
  bool operator!=(const RooFIter_IteratorHelper& other);
  bool operator!=(const RooAbsArg* other);
  RooAbsArg* operator*();
};
RooFIter_IteratorHelper begin(RooFIter& it);
RooFIter_IteratorHelper end(RooFIter& it);
struct TIterator_IteratorHelper { //EXCLUDE
  TIterator* itr;
  TObject* nextItem;
  TIterator_IteratorHelper(TIterator* it, bool end);
  TObject* operator++();
  bool operator!=(const TIterator_IteratorHelper& other);
  bool operator!=(const TObject* other);
  TObject* operator*();
};
TIterator_IteratorHelper begin(TIterator* it);
TIterator_IteratorHelper end(TIterator* it);


#if ROOT_VERSION_CODE < ROOT_VERSION(6,10,0)
TIterator* clients(const RooAbsArg* arg);
TIterator* servers(const RooAbsArg* arg);
#else
const RooAbsArg::RefCountList_t & clients(const RooAbsArg* arg);
const RooAbsArg::RefCountList_t & servers(const RooAbsArg* arg);
#endif


#endif
