%module pair

%{
#include "pair.h"
#include "pair_lj.h"
#include "pair_lj_multi.h"
#include "pair_patch_kf.h"
#include "pair_lj_coul_ewald.h"
#include "pair_ideal.h"
#include "pair_wall.h"
#include "pair_round_square.h"
#include "pair_hard_circle.h"
#include "pair_squarewell.h"
#include "pair_hybrid.h"
#include "pair_hs.h"
#include "pair_tabular.h"
#include "pair_tabular_1d.h"
%}

%pythonnondynamic;

%include pair.h
%include pair_lj.h
%include pair_lj_multi.h
%include pair_patch_kf.h
%include pair_lj_coul_ewald.h
%include pair_ideal.h
%include pair_wall.h
%include pair_round_square.h
%include pair_hard_circle.h
%include pair_squarewell.h
%include pair_hybrid.h
%include pair_hs.h
%include pair_tabular.h
%include pair_tabular_1d.h

