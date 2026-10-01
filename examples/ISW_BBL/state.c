/*
 * File: state.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains functions that define and implement the equation of state for
 * the density from the salinity, temperature, and pressure.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "state.h"

/*
 * Function: StateEquation
 * Usage: rho = StateEquation(prop,s,T,p);
 * ---------------------------------------
 * Returns the density as a function of temperature, salinity, and
 * pressure, where pressure is the hydrostatic pressure p=RHO0*GRAV*z,
 * and RHO0 and GRAV are defined in suntans.h.  Note that rho should
 * always normalized by RHO0 so that this function returns a dimensionless
 * quantity.
 *
 */
REAL StateEquation(const propT *prop, const REAL s, const REAL T, const REAL p) {
  return prop->beta*s/1e5; //1e5 because of the rigid lid gravity
  //REAL rho_s = 2650.0; 
  //return prop->beta*s/1e5*(1-(T/rho_s)) + rho_s*(T/rho_s)/RHO0/1e5;
  //return prop->beta*s;
}


