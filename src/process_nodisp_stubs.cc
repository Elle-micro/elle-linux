#include <cstdio>

#include "setup.h"

extern "C" {

// Shared headless stubs used by legacy process mains that still call
// SetupApp/Run_App. Individual processes (e.g. latte, latte-clean-old)
// may provide their own nodisp_stubs.cc; if so, we don't add these to
// avoid duplicate symbol definitions.

// Mark as weak so that if any linked object file (process-specific)
// defines a real SetupApp or Run_App, it overrides these defaults
// without causing multiple definition link errors.
#if defined(__GNUC__)
__attribute__((weak))
#endif
int SetupApp(int argc, char** argv) {
  (void)argc;
  (void)argv; // suppress unused warnings
  return 0;
}

#if defined(__GNUC__)
__attribute__((weak))
#endif
int Run_App(FILE* fp) {
  (void)fp; // headless path directly calls init/run functions
  return 0;
}

} // extern "C"
