//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef _UTILS_
#define _UTILS_

#include <string>
#include <vector>

#include "TTime.h"

class RooProdPdf;

namespace RooFitUtils {
  extern bool RooStarMomentMorphFix;
  extern bool RooMultiPdfFix;

  int fixRooStarMomentMorph(RooWorkspace* workspace);
  
  bool AlmostEqualUlpsAndAbs(float A, float B, float maxDiff, int maxUlpsDiff);
  void PrintResourcesUsed(const TTime& progStart);
  std::vector<std::string> parseString(const std::string& str, const std::string& sep);
  void FindUniqueProdComponents( RooProdPdf* Pdf, RooArgSet& Components );
  
  bool ensureDirectory(const TString& path);
  bool ensureDirectoryForFile(const TString& file);
  
  void PrintTable(std::string* firstCol, std::string** matrix, std::string** matrixErr, std::string* header, int nrRows, int nrCols, int nSigFig, std::ostream& ost, std::string indent = "", std::string delim = " & ", std::string ending = " \\\\");
}

#endif
