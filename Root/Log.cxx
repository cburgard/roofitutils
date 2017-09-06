#include "RooFitUtils/Log.h"
#include <iostream>
#include <sys/time.h>

namespace {
inline std::string GetTime() {
  struct timeval tv;
  gettimeofday(&tv, 0);
  char buffer[100];
  tm r;
  strftime(buffer, sizeof(buffer), "%X", localtime_r(&tv.tv_sec, &r));
  char result[100];
  sprintf(result, "%s.%06ld", buffer, (long)tv.tv_usec);
  return result;
}
}

// _____________________________________________________________________________

RooFitUtils::Log::Log(RooFitUtils::LogLevel _loglevel) {
  _os << "- " << GetTime();
  _os << " " << ToString(_loglevel) << ": ";
  _os << std::string(_loglevel > logDEBUG ? (_loglevel - logDEBUG) : 0, '\t');
}

// _____________________________________________________________________________

RooFitUtils::Log::~Log() {
  _os << std::endl;
  std::cout << _os.str();
}

static RooFitUtils::LogLevel reportingLevel = RooFitUtils::logDEBUG;

// _____________________________________________________________________________
RooFitUtils::LogLevel &RooFitUtils::Log::ReportingLevel() {
  return reportingLevel;
}

// _____________________________________________________________________________
void RooFitUtils::Log::SetReportingLevel(RooFitUtils::LogLevel lvl) {
  reportingLevel = lvl;
}

// _____________________________________________________________________________
std::string RooFitUtils::Log::ToString(RooFitUtils::LogLevel _loglevel) {
  static const char *const buffer[] = {"ERROR", "WARNING", "INFO", "DEBUG"};
  return buffer[_loglevel];
}

// _____________________________________________________________________________
RooFitUtils::LogLevel
RooFitUtils::Log::FromString(const std::string &_loglevel) {
  if (_loglevel == "DEBUG")
    return logDEBUG;
  if (_loglevel == "INFO")
    return logINFO;
  if (_loglevel == "WARNING")
    return logWARNING;
  if (_loglevel == "ERROR")
    return logERROR;
  Log(logWARNING) << "Unknown logging level '" << _loglevel
                  << "'. Using INFO level as default.";
  return logINFO;
}
