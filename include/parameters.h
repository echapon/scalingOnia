#ifndef parameters_h
#define parameters_h

#include "utils.h"

// global settings
const Escaling       gscaling = pt2; // pt2 or mtpt
const Einterpolation ginterpolation = logcspline; // lin, cspline, loglin, logcspline
float                gTextSize = 0.04;
// should we use the Lafferty & Wyatt prescription to change the x position of the points? (assuming a locally exponentially falling spectrum)
bool                 doxLW = true;//true;
Elwmode              lwmode = powlaw; // expo or powlaw
bool                 plotxt = true;//true;

#endif // #ifndef parameters_h
