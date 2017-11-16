#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "../allvars.h"
#include "../proto.h"

/*! \file twopoint.c
 *  \brief computes the two-point mass correlation function on the fly
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel (volker.springel@h-its.org). It is here in GIZMO
 * as legacy code at the moment, and needs to be re-written or removed.
 */

/* Note: This routine will only work correctly for particles of equal mass ! */


