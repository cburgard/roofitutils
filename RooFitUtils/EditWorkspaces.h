//this file looks like plain C, but it's actually -*- c++ -*-
#ifndef EDITWORKSPACE
#define EDITWORKSPACE

#include <RooWorkspace.h>
#include <TDirectory.h>
#include <vector>
#include <string>

int editws (const std::vector<std::string>& lines,RooWorkspace*w,const char* mcname=0);
int editws (const RooWorkspace* w, RooWorkspace& wout, const char* mcname=0, const char* dataname=0);
int editws (const char* edit, RooWorkspace* w, const char* mcname=0);
int editws (RooWorkspace* w, TDirectory* fout, const char* edit=0, const char* mcname=0, const char* dataname=0);
int editws (RooWorkspace* w, const char* outfile, const char* edit=0, const char* mcname=0, const char* dataname=0);
int editws (TDirectory* f, TDirectory* fout, const char* edit=0, const char* wsname=0, const char* mcname=0, const char* dataname=0);
int editws (TDirectory* f, const char* outfile, const char* edit=0, const char* wsname=0, const char* mcname=0, const char* dataname=0);
int editws (const char* outfile, const char* edit=0, const char* wsfile=0, const char* wsname=0, const char* mcname=0, const char* dataname=0);
int editws (RooWorkspace* w, RooWorkspace* wout, const char* edit, const char* mcname=0, const char* dataname=0);
int editws (RooWorkspace* w, RooWorkspace* wout, const std::vector<std::string>& edits, const char* mcname=0, const char* dataname=0);

#endif
