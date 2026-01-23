#include "timefn.h"
#include <iostream.h>
#include <string>

int main(int argc, char **argv) {
  char *t;

  string x = "time is ";
  x += GetLocalTime();
  cout << x;
}
