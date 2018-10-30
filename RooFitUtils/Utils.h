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

#endif
