#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"

/* compute line-of-sight integrated quantities (for e.g. Lyman-alpha forest studies) */

/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). It is here in GIZMO
 * as legacy code at the moment, and needs to be re-written or removed.
 */


