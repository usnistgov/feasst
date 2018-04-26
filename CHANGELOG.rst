*************************************
Change log
*************************************

Depreciated from v0.3.1 to v0.4
####################################

Log file printing is now almost completely controlled by the Trial class instead of MC.

Renamed
***************

Space::l() -> boxLength()

Space::l(const int)  -> boxLength(const int)

Space::vol() -> volume()

Accumulator::stdev() -> std()

Space::readxyz2 -> readXYZ

Space/Pair::readxyz -> readXYZ

Space::randDispMulti -> randDispNoWrap

Space::printxyz -> printXYZ

Pair*::printxyz -> printXYZ

Space::contact2clusterAlt > contact2cluster

Space::floodFillContactAlt -> floodFillContact

Removed
***************

PairLJCoulEwald::delPart(const int)

Functions.vecSPCE

Space::lset(double, int)

Space::pbc(double, int)

Space::ranDisp(int, double)

Space::checkBond(char*, double)

Space::contact2cluster

Space::floodFillContact

Space::nRadialHist/printRadial

Space::floppyBox

Analyze::initPrintFreq

Criteria::printBeta

Criteria::printPressure

WLTMMC::nMolResizeWindow

WLTMMC::initGR, nFreqGR, GRFileName, grt, gr

Histogram::iType/jType

Random::int32

