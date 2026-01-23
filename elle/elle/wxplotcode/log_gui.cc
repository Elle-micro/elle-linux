// Minimal GUI-side logging implementation to satisfy wx sources
// Matches declaration in wxplotcode/settings.h

#include <cstdio>

extern "C" void Log(int /*level*/, const char *msg) {
  if (msg)
    std::fprintf(stderr, "%s\n", msg);
}
