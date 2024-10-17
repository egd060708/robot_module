#pragma once
#include <iostream>

#define ROWS 162
#define COLS 2

extern double forelimb[ROWS][COLS];
extern double forelimb_cont[ROWS][COLS];
extern double hindlimb[ROWS][COLS];
extern double hindlimb_cont[ROWS][COLS];

void getGait(double _dst[4][2], double _Ts, double _t, double _start_t);