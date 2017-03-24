/**
 * \mainpage
 *
 * Developed by Harold Wickes Hatch, 12/13/2013, hhatch.com, harold@hhatch.com
 *aa
 * 
 */

#include "mc_wltmmc.h"

int main() {

  // set input variables
  std::ostringstream rstFileName("tmp/rst");
  
  // read restart file
  WLTMMC mc(rstFileName.str().c_str());

  // run sweeps
  mc.runNumSweepsRestart(100, rstFileName.str().c_str());
}
