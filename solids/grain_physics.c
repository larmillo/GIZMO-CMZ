#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*
 
 This module contains the self-contained sub-routines needed for
 grain-specific physics in proto-planetary/proto-stellar/planetary cases.
 It's also potentially use-able for GMC and ISM scales, and terrestrial
 turbulence. Anything where aerodynamic particles are interesting
 
 
 This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 
 */



#ifdef GRAIN_FLUID


/* function to apply the drag on the grains from surrounding gas properties */
void apply_grain_dragforce(void)
{
    
    CPU_Step[CPU_MISC] += measure_time();
    
    int i, k;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if((P[i].Type != 0)&&(P[i].Type != 4))
        {
            if(P[i].Gas_Density > 0)
            {
                double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
                if(dt > 0)
                {
                    double cs = sqrt( GAMMA * GAMMA_MINUS1 * P[i].Gas_InternalEnergy);
                    double R_grain_cgs = P[i].Grain_Size;
                    double R_grain_code = R_grain_cgs / (All.UnitLength_in_cm / All.HubbleParam);
                    double rho_gas = P[i].Gas_Density * All.cf_a3inv;
                    double rho_grain_physical = All.Grain_Internal_Density; // cgs units //
                    double rho_grain_code = rho_grain_physical / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam); // code units //
                    double vgas_mag = 0.0;
                    for(k=0;k<3;k++) {vgas_mag+=(P[i].Gas_Velocity[k]-P[i].Vel[k])*(P[i].Gas_Velocity[k]-P[i].Vel[k]);}
                    
                    
                    if(vgas_mag > 0)
                    {
                        vgas_mag = sqrt(vgas_mag) / All.cf_atime;
                        double x0 = 0.469993*sqrt(GAMMA) * vgas_mag/cs; // (3/8)*sqrt[pi/2]*|vgas-vgrain|/cs //
                        double tstop_inv = 1.59577/sqrt(GAMMA) * rho_gas * cs / (R_grain_code * rho_grain_code); // 2*sqrt[2/pi] * 1/tstop //
#ifdef GRAIN_EPSTEIN
                        double mu = 2.3 * PROTONMASS;
                        double temperature = mu * (P[i].Gas_InternalEnergy*All.UnitEnergy_in_cgs*All.HubbleParam/All.UnitMass_in_g) / BOLTZMANN;
                        double cross_section = GRAIN_EPSTEIN * 2.0e-15 * (1. + 70./temperature);
                        cross_section /= (All.UnitLength_in_cm * All.UnitLength_in_cm / (All.HubbleParam*All.HubbleParam));
                        double n_mol = rho_gas / (mu * All.HubbleParam/All.UnitMass_in_g);
                        double mean_free_path = 1 / (n_mol * cross_section); // should be in code units now //
                        double corr_mfp = R_grain_code / ((9./4.) * mean_free_path);
                        if(corr_mfp > 1) {tstop_inv /= corr_mfp;}
#endif
                        double C1 = (-1-sqrt(1+x0*x0)) / x0;
                        double xf = 0.0;
                        double dt_tinv = dt * tstop_inv;
                        if(dt_tinv < 100.)
                        {
                            double C2 = C1 * exp( dt_tinv );
                            xf = -2 * C2 / (C2*C2 -1);
                        }
                        double slow_fac = 1 - xf / x0;
                        // note that, with an external (gravitational) acceleration, we can still solve this equation for the relevant update //

                        double external_forcing[3];
                        for(k=0;k<3;k++) {external_forcing[k] = 0;}
                        /* this external_forcing parameter includes additional grain-specific forces. note that -anything- which imparts an 
                            identical acceleration onto gas and dust will cancel in the terms in t_stop, and just act like a 'normal' acceleration
                            on the dust. for this reason the gravitational acceleration doesn't need to enter our 'external_forcing' parameter */
#ifdef GRAIN_LORENTZFORCE
                        /* Lorentz force on a grain = Z*e/c * ([v_grain-v_gas] x B) */
                        double v_cross_B[3];
                        v_cross_B[0] = (P[i].Vel[1]-P[i].Gas_Velocity[1])*P[i].Gas_B[2] - (P[i].Vel[2]-P[i].Gas_Velocity[2])*P[i].Gas_B[1];
                        v_cross_B[1] = (P[i].Vel[2]-P[i].Gas_Velocity[2])*P[i].Gas_B[0] - (P[i].Vel[0]-P[i].Gas_Velocity[0])*P[i].Gas_B[2];
                        v_cross_B[2] = (P[i].Vel[0]-P[i].Gas_Velocity[0])*P[i].Gas_B[1] - (P[i].Vel[1]-P[i].Gas_Velocity[1])*P[i].Gas_B[0];

                        double grain_mass = (4.*M_PI/3.) * R_grain_code*R_grain_code*R_grain_code * rho_grain_code; // code units
                        double lorentz_units = sqrt(4.*M_PI*All.UnitPressure_in_cgs*All.HubbleParam*All.HubbleParam); // code B to Gauss
                        lorentz_units *= (ELECTRONCHARGE/C) * All.UnitVelocity_in_cm_per_s / (All.UnitMass_in_g / All.HubbleParam); // converts acceleration to cgs
                        lorentz_units /= All.UnitVelocity_in_cm_per_s / (All.UnitTime_in_s / All.HubbleParam); // converts it to code-units acceleration

                        /* calculate the grain charge following Draine & Sutin */
                        double cs_cgs = cs * All.UnitVelocity_in_cm_per_s;
                        double tau_draine_sutin = R_grain_cgs * (2.3*PROTONMASS) * (cs_cgs*cs_cgs) / (ELECTRONCHARGE*ELECTRONCHARGE);
                        double Z_grain = -DMAX( 1./(1. + sqrt(1.0e-3/tau_draine_sutin)) , 2.5*tau_draine_sutin );
                        if(isnan(Z_grain)||(Z_grain>=0)) {Z_grain=0;}
                        
                        /* define unit vectors and B for evolving the lorentz force */
                        double bhat[3]={0}, bmag=0, dv[3]={0}, efield[3]={0}, efield_coeff=0;
                        for(k=0;k<3;k++) {bhat[k]=P[i].Gas_B[k]; bmag[k]+=bhat[k]*bhat[k]; dv[k]=P[i].Vel[k]-P[i].Gas_Velocity[k];}
                        if(bmag>0) {bmag=sqrt(bmag); for(k=0;k<3;k++) {bhat[k]/=bmag;}} else {bmag=0;}
                        double lorentz_coeff = (0.5*dt) * bmag * Z_grain / grain_mass * lorentz_units; // multiply in full timestep //
                        
                        /* now apply the boris integrator */
                        double v_m[3]={0}, v_t[3]={0}, v_p[3]={0}, vcrosst[3]={0};
                        for(k=0;k<3;k++) {v_m[k] = dv[k] + 0.5*efield_coeff*efield[k];} // half-step from E-field
                        /* cross-product for rotation */
                        vcrosst[0] = v_m[1]*bhat[2] - v_m[2]*bhat[1]; vcrosst[1] = v_m[2]*bhat[0] - v_m[0]*bhat[2]; vcrosst[2] = v_m[0]*bhat[1] - v_m[1]*bhat[0];
                        for(k=0;k<3;k++) {v_t[k] = v_m[k] + lorentz_coeff * vcrosst[k];} // first half-rotation
                        vcrosst[0] = v_t[1]*bhat[2] - v_t[2]*bhat[1]; vcrosst[1] = v_t[2]*bhat[0] - v_t[0]*bhat[2]; vcrosst[2] = v_t[0]*bhat[1] - v_t[1]*bhat[0];
                        for(k=0;k<3;k++) {v_p[k] = v_m[k] + (2.*lorentz_coeff/(1.+lorentz_coeff*lorentz_coeff)) * vcrosst[k];} // second half-rotation
                        for(k=0;k<3;k++) {v_p[k] += 0.5*efield_coeff*efield[k];} // half-step from E-field
                        /* calculate effective acceleration from discrete step in velocity */
                        for(k=0;k<3;k++) {external_forcing[k] += (v_p[k] - dv[k]) / dt;}
                        /* note: if grains moving super-sonically with respect to gas, and charge equilibration time is much shorter than the 
                            streaming/dynamical timescales, then the charge is slightly reduced, because the ion collision rate is increased while the 
                            electron collision rate is increased less (since electrons are moving much faster, we assume the grain is still sub-sonic 
                            relative to the electron sound speed. in this case, for the large-grain limit, the Draine & Sutin results can be generalized; 
                            the full expressions are messy but can be -approximated- fairly well for Mach numbers ~3-30 by simply 
                            suppressing the equilibrium grain charge by a power ~exp[-0.04*mach]  (weak effect, though can be significant for mach>10) */
#endif
                        
                        double delta_egy = 0;
                        double delta_mom[3];
                        double dv[3];
                        for(k=0; k<3; k++)
                        {
                            /* measure the imparted energy and momentum as if there were no external acceleration */
                            double v_init = P[i].Vel[k];
                            double vel_new = v_init + slow_fac * (P[i].Gas_Velocity[k]-v_init);
                            delta_mom[k] = P[i].Mass * (vel_new - v_init);
                            delta_egy += 0.5*P[i].Mass * (vel_new*vel_new - v_init*v_init);
                            /* now calculate the updated velocity accounting for any external, non-standard accelerations */
                            double vdrift = 0;
                            if(tstop_inv > 0) {vdrift = external_forcing[k] / (tstop_inv * sqrt(1+x0*x0));}
                            dv[k] = slow_fac * (P[i].Gas_Velocity[k] - v_init + vdrift);
                            if(isnan(vdrift)||isnan(slow_fac)) {dv[k] = 0;}
                            /* note, we can directly apply this by taking P[i].Vel[k] += dv[k]; but this is not as accurate as our
                                normal leapfrog integration scheme.
                                we can also account for the -gas- acceleration, by including it like vdrift;
                                for a constant t_stop, the gas acceleration term appears as 
                                P[i].Vel[l] += Gas_Accel[k] * dt + slow_fac * (Gas-Accel[k] / tstop_inv) */
                            /* note that we solve the equations with an external acceleration already (external_forcing above): therefore add to forces
                             like gravity that are acting on the gas and dust in the same manner (in terms of acceleration) */
                            P[i].GravAccel[k] += dv[k] / dt;
                            //P[i].Vel[k] += dv[k];
                        }

                    

                    } // closes check for if(v_mag > 0)
                } // closes check for if(dt > 0)
            } // closes check for if(P[i].Gas_Density > 0)
        } // closes check for if(P[i].Type != 0)
    } // closes main particle loop
    
    CPU_Step[CPU_DRAGFORCE] += measure_time();
    
}











#endif


