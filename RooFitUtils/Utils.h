// this file looks like plain C, but it's actually -*- c++ -*-
#ifndef _UTILS_
#define _UTILS_

#include <string>
#include <vector>

#include "TTime.h"

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
  RooWorkspace* makeCleanWorkspace(RooWorkspace* oldWS, const char* newName = NULL);
  
  void PrintTable(std::string *firstCol, std::string **matrix,
                  std::string **matrixErr, std::string *header, int nrRows,
                  int nrCols, int nSigFig, std::ostream &ost,
                  std::string indent = "", std::string delim = " & ",
                  std::string ending = " \\\\");
	
  bool matches(const std::string& text, const std::string& pattern);
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
TIterator_IteratorHelper begin(TObject* it);
TIterator_IteratorHelper end(TObject* it);

#endif
