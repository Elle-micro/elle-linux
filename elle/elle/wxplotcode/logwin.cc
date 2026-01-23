// Written in 2003
// Author: Dr. J.K. Becker
// Copyright: Dr. J.K. Becker (becker@jkbecker.de)

#include "logwin.h"

LogWin::LogWin(wxWindow *parent) {
  Create();
  logframe = new wxLogWindow(parent, wxT("Logging window"), false, true);
  Run();
  logframe->Show(true);
}
