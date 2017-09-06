//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef _ROOFITUTILS_LOG_
#define _ROOFITUTILS_LOG_
#include <string>
#include <sstream>

namespace RooFitUtils {

enum LogLevel {logERROR, logWARNING, logINFO, logDEBUG};

class Log {
public:
  Log(LogLevel _loglevel);
  virtual ~Log();
public:
  static void SetReportingLevel(LogLevel _loglevel);
  static LogLevel& ReportingLevel();
  static std::string ToString(LogLevel _loglevel);
  static LogLevel FromString(const std::string& _loglevel);
  template <typename T>  Log& operator <<(T const& value){ _os << value; return *this; }
protected:
  std::ostringstream _os;
};

}

#define LOG(_loglevel)  if (_loglevel <= RooFitUtils::Log::ReportingLevel()) RooFitUtils::Log(_loglevel)
#endif
