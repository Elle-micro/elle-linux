// Minimal stub of the IPhreeqc C++ API to allow building without the real
// library. This is ONLY for compilation; it does not perform any geochemical
// calculations. If the real IPhreeqc headers are present on the system, they
// will be preferred.

#pragma once

#include <string>

// Value/result types mimicking the real API just enough for compilation
enum VARTYPE { TT_LONG = 0, TT_DOUBLE = 1, TT_STRING = 2 };

struct VAR {
  VARTYPE type{TT_LONG};
  long lVal{0};
  double dVal{0.0};
  const char *sVal{nullptr};
};

// Return code used by GetSelectedOutputValue etc.
static constexpr int VR_OK = 0;

class IPhreeqc {
public:
  IPhreeqc() = default;
  ~IPhreeqc() = default;

  // Database and error output
  int LoadDatabase(const char *) { return 0; }
  void OutputErrorString() {}
  void OutputAccumulatedLines() {}
  void ClearAccumulatedLines() {}

  // Output file controls (no-op)
  void SetOutputFileOn(bool) {}
  void SetOutputFileName(const char *) {}

  // Command accumulation and execution (no-op)
  void AccumulateLine(const char *) {}
  int RunAccumulated() { return 0; }
  int RunString(const char *) { return 0; }

  // Version string (dummy)
  const char *GetVersionString() const { return "stub"; }

  // Selected output accessors (empty table)
  int GetSelectedOutputRowCount() const { return 0; }
  int GetSelectedOutputColumnCount() const { return 0; }
  int GetSelectedOutputValue(int, int, VAR *v) const {
    if (v) {
      v->type = TT_STRING;
      v->sVal = "";
    }
    return VR_OK;
  }
};

// C-style helpers present in the real C API, used by legacy code paths
inline void VarInit(VAR *v) {
  if (!v)
    return;
  v->type = TT_LONG;
  v->lVal = 0;
  v->dVal = 0.0;
  v->sVal = nullptr;
}

inline void VarClear(VAR *v) {
  if (!v)
    return;
  // Nothing to release in stub; real API would free string buffers when needed
  v->type = TT_LONG;
  v->lVal = 0;
  v->dVal = 0.0;
  v->sVal = nullptr;
}
