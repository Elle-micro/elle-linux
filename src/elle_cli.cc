// minimal headless CLI for link/smoke tests (behavior-neutral)
#include <cstring>
#include <iostream>
#include <string>

// Pull a safe, side-effect-free query from elle_core
#include "version.h" // provides ElleGetLibVersionString()

static void print_help(const char *prog) {
  std::cout << "Usage: " << prog << " [--version] [--selftest]\n";
}

int main(int argc, char **argv) {
  const char *prog = (argc > 0 && argv[0]) ? argv[0] : "elle_cli";

  if (argc <= 1) {
    print_help(prog);
    return 0; // modernize: show help by default (analysis fast-track,
              // behavior-neutral)
  }

  if (std::strcmp(argv[1], "--version") == 0) {
    // Prefer the build's PROJECT_VERSION injected via -DELLE_VERSION (Task A20)
    // modernize: prefer ELLE_VERSION define with legacy fallback (analysis
    // fast-track)
#ifdef ELLE_VERSION
    std::cout << "elle_cli (elle_core) version " << ELLE_VERSION << "\n";
#else
    std::string ver = ElleGetLibVersionString();
    std::cout << "elle_cli (elle_core) version " << ver << "\n";
#endif
    return 0; // modernize: simple version print (analysis fast-track,
              // behavior-neutral)
  }

  if (std::strcmp(argv[1], "--selftest") == 0) {
    // Selftest: call a harmless function from elle_core; verify we can link and
    // return success
    volatile std::string ver = ElleGetLibVersionString();
    (void)ver; // suppress unused warning
    std::cout << "selftest: OK\n";
    return 0; // modernize: link smoke test only (analysis fast-track,
              // behavior-neutral)
  }

  print_help(prog);
  return 1;
}
