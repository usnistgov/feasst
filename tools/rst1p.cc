/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *
 * This example restarts a single processor simulation
 * 
 */

#include "mc_wltmmc.h"

int main() {

  long long int npr = 1e6;

  // read checkpoint files
  WLTMMC mc("tmp/rst");
  
  // run simulation
  mc.runNumTrials(npr);
}
