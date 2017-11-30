#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <ctype.h>

#include "allvars.h"
#include "proto.h"
#include "kernel.h"

void gas_to_star()
{
	double tff, sfr, mstar, Prandom, Jeans_mass, stars_converted, tot_stars_converted;
	double cspeed, RJ, starfraction, star_probability, dtime;
	double alpha, fsh, tau_a, phi_a, theta_a, ModRhoGrad, sfrate = 0;
	double sfrrate, totsfrrate, sum_sm, total_sm, sum_mass_stars, total_sum_mass_stars, rate, rate_in_msunperyear;
	sum_sm = sum_mass_stars = 0;
	int i,j,k;
	Stars_converted = stars_converted = 0;
	rearrange_particle_sequence();
	
	hydro_gradient_calc();
	
	//TESTARE!!!! 
	//for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinActive[bin]) {TimeBinSfr[bin] = 0;}}
		
	/* My implementation for SPH code */ 
	for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{	
	if((P[i].Type == 0)&&(P[i].Mass>0)&&(dtime>0)&&(P[i].TimeBin)) //TESTARE ULTIMA CONDIZIONE!!!!
	{	
		//dtime = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; //TESTARE!!!!
		//SphP[i].Sfr = 0;
			
		SphP[i].Pressure = get_pressure(i);
		cspeed = sqrt(GAMMA * SphP[i].Pressure / Particle_density_for_energy_i(i));
		VelGradTens = 0.;
		for(k = 0; k < 3; k++)
		{
			for(j = 0; j < 3; j++) VelGradTens += SphP[i].Gradients.Velocity[j][k]*SphP[i].Gradients.Velocity[j][k];
		}
		alpha = (VelGradTens + pow(cspeed/Get_Particle_Size(i),2))/(8*PI_VAL*All.G*SphP[i].density);
		/* Self-gravitating criterion */	
		if(alpha < 1)
		{
			//sfrate = 1;
			/* Self-shielding criterion */
			ModRhoGrad = 0.;
			for (k = 0; k < 3; k++) ModRhoGrad += SphP[i].Gradients.Density[k]*SphP[i].Gradients.Density[k];
			ModRhoGrad = sqrt(ModRhoGrad);
			theta_a = 0.756*pow(1.+3.1*P[i].metallicity,0.365);
			tau_a = 434.8 * SphP[i].density * (Get_Particle_Size(i)+SphP[i].density/ModRhoGrad); //cgs units - trasformare in code units cm**-2 g**-1
			phi_a = 0.6*tau_a*(0.01+P[i].metallicity)/(log(1.+0.6*theta_a+0.01*theta_a*theta_a));
			fsh = 1.-3./(1.+4./phi_a);
			
			if(fsh > 0)
			{	 
				//sfrate *= fsh;			
				// Dovremmo calcolare la temperatura e quindi la densità numerica
				/* Sufficiently-Dense */
				if (SphP[i].Density > All.CritPhysDensity)
				{
					RJ = sqrt(PI_VAL*cspeed*cspeed/All.G/SphP[i].density);
					Jeans_mass = 4./3.*PI_VAL*SphP[i].density*pow(RJ,3.);
					/* Jeans unstable criterion */
					if (P[i].Mass > Jeans_Mass) 
					{
						t_ff = sqrt(3*PI_VAL/(32*All.G*SphP[i].Density));
						
			            //SphP[i].Sfr = sm / dtime * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
			  	      	//if(dtime>0) TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
					  
						star_probability = All.SfEffPerFreeFall * All.time / tff;
						Prandom = get_random_number(ThisTask);
						if (star_probability > Pramdom)
						{
							P[i].Type = 4;
			  		      	//TimeBinCountSph[P[i].TimeBin]--;
			  		      	//TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;
						  
							Stars_converted+=1;
							stars_converted+=1;
							P[i].StellarAge = All.Time;
							sum_mass_stars += P[i].Mass;
							
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
                			P[i].DensAroundStar = SphP[i].Density;
#endif
							
#ifdef SLUG
							//Slug part
							double particle_mass = P[i].Mass;
							P[i].SlugOb = slug_object_new();
							slug_construct_cluster((void*)(P[i].SlugOb), double particle_mass);
#endif						
						}		
					}
				}
			}	
		}
	}} /* end of main loop over active particles*/
	MPI_Allreduce(&stars_converted, &total_stars_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	All.TotN_gas -= tot_converted;
	rearrange_particle_sequence();
	
    //Write File
	
    //for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
      //if(TimeBinCount[bin])
        //sfrrate += TimeBinSfr[bin];
	
	MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        if(All.TimeStep > 0)
            rate = total_sm / (All.TimeStep / (All.cf_atime*All.cf_hubble_a));
        else
            rate = 0;
        /* convert to solar masses per yr */
        rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
        fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars);
        fflush(FdSfr); // can flush it, because only occuring on master steps anyways
    } // thistask==0
		
	CPU_Step[CPU_COOLINGSFR] += measure_time();		
}/* end of gas_to_star routine*/

