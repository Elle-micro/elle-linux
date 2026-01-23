// Wrapper added to expose InitThisProcess symbol under the stem 'plot_ugrid'
// so automated CMake glob (matching stem*.cc) picks it up alongside plot_ugrid.main.cc.
// The real implementation lives in plot_misorient_unodebound.elle.cc.

extern int InitThisProcess(); // implemented in plot_misorient_unodebound.elle.cc

// Nothing else needed; presence of this TU ensures the linker sees InitThisProcess
// when building the plot_ugrid target without renaming legacy file.