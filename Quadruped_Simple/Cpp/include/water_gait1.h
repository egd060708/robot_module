#pragma once
#include <iostream>


#define WG_ROWS 162
#define WG_COLS 2

extern double forelimb[WG_ROWS][WG_COLS];
extern double forelimb_cont[WG_ROWS][WG_COLS];
extern double hindlimb[WG_ROWS][WG_COLS];
extern double hindlimb_cont[WG_ROWS][WG_COLS];

void getGait(double _dst[4][2], double _Ts, double _t, bool _is_gait);
void slopeGait(double _dst[4][2], double _tar[4][2], double _cur[4][2], double slope);