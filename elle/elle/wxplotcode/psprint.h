// Written in 2003
// Author: Dr. J.K. Becker
// Copyright: Dr. J.K. Becker (becker@jkbecker.de)
#ifndef _E_printps_h
#define _E_printps_h

#include "dsettings.h"
#include "file.h"
#include "general.h"
#include "init.h"
#include "interface.h"
#include "nodes.h"
#include "settings.h"
#include "unodes.h"
#include "wx/cmndata.h"
#include "wx/dcps.h"
#include "wx/generic/printps.h"
#include "wx/generic/prntdlgg.h"
#include "wx/print.h"
#include "wx/wx.h"

class PSPrint : public wxPrintout {
  DECLARE_CLASS(PSPrint)
public:
  PSPrint(wxWindow *parent);
  bool OnPrintPage(int n);

private:
  void DrawBNodes();
  void DrawLines();
  void CalcPolygon();
  void RWCoords(double x, double y, int *rwx, int *rwy);
  Settings *settings;
};
#endif //_E_printps_h
