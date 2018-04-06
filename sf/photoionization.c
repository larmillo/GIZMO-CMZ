#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#ifdef PHOTOIONIZATION

void HII_region(void)
{
    double *Distance = (double*)malloc(N_gas*sizeof(double*));
	double *IonRate = (double*)malloc(N_gas*sizeof(double*));
	double *Tini = (double*)malloc(N_gas*sizeof(double*));
	int *ParticleNum = (int*)malloc(N_gas*sizeof(int*));
	//int *Tag_HIIregion = (int*)malloc(N_gas*sizeof(int*));
	double ord1,ord2, ord5;
	int ord3,ord4;
	double fbtime = 1e20, tcooling;
	double Tfin = 1e4;
	double molw_i = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* assuming full ionization */
	
	//for(int j = 0; j < N_gas; j++) SphP[j].HIIregion=0; 
		
    for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 4) continue;
        if(P[i].Mass <= 0) continue;
		P[i].Feedback_timestep = fbtime;
		
		//get number of photons per second
		P[i].N_photons = slug_get_photometry_QH0(P[i].SlugOb);
		if(P[i].N_photons <= 0) continue;
		//printf("Initial number of photons = %e \n", P[i].N_photons);
		int count = 0;
		for(int j = 0; j < N_gas; j++) /* loop over the gas block */
		{
			if (SphP[j].HIIregion==1) continue; // The particle belongs to another HII region
			double distx = P[i].Pos[0] - P[j].Pos[0];
			double disty = P[i].Pos[1] - P[j].Pos[1];
			double distz = P[i].Pos[2] - P[j].Pos[2];
			double dist = sqrt(distx*distx+disty*disty+distz*distz);
			if (dist>300.) continue;
			Distance[count] = sqrt(distx*distx+disty*disty+distz*distz);
			ParticleNum[count] = j;
			double Rhob = SphP[j].Density * All.UnitDensity_in_cgs;
			double Mb   = P[j].Mass * All.UnitMass_in_g;
			Tini[count] = CallGrackle(SphP[j].InternalEnergy,SphP[j].Density,0,&(SphP[j].Ne),j,2);
			double molw_n = Tini[count]*BOLTZMANN/(GAMMA-1)/(SphP[j].InternalEnergy*All.UnitEnergy_in_cgs/All.UnitMass_in_g)/PROTONMASS;
			double beta = 3e-13; //cm**3s*-1, cgs units
			IonRate[count] = HYDROGEN_MASSFRAC*beta*Tini[j]*Rhob*Mb/(2*PROTONMASS*PROTONMASS*molw_n*molw_i*Tfin); 
			count += 1;
			//printf("IonRate Tini Mb Rhob %e %e %e %e \n", IonRate[j], Tini[j], Mb, molw);
		}	
		printf("Number of particles to be photionized %d %d \n", P[i].ID, count+1);
		//Cycle to sort particles by increasing distance
        for (int j = 0; j < (count - 1); j++)
		{
		    for (int k = 0; k < count - 1 - j; k++)
		    {
		        if (Distance[k] > Distance[k+1])
		        {
		            ord1 = Distance[k+1];
		            Distance[k+1] = Distance[k];
		            Distance[k] = ord1;
					
		            ord2 = IonRate[k+1];
		            IonRate[k+1] = IonRate[k];
		            IonRate[k] = ord2;
					
		            ord4 = Tini[k+1];
		            Tini[k+1] = Tini[k];
		            Tini[k] = ord4;
					
		            ord3 = ParticleNum[k+1];
		            ParticleNum[k+1] = ParticleNum[k];
		            ParticleNum[k] = ord3;
		        }
		     }
		}
		
		
		for(int j = 0; j < count; j++) /* loop over the gas block */
		{
			if (Tini[j] >= Tfin) continue; //Particle already ionized
			if (IonRate[j] <= P[i].N_photons) 
			{
				SphP[ParticleNum[j]].InternalEnergy = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
				SphP[ParticleNum[j]].InternalEnergyPred = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
				SphP[ParticleNum[j]].HIIregion = 1;
				//double T = CallGrackle(SphP[ParticleNum[j]].InternalEnergy,SphP[ParticleNum[j]].Density,0,&(SphP[ParticleNum[j]].Ne),ParticleNum[j],2);
				tcooling= CallGrackle(SphP[ParticleNum[j]].InternalEnergy,SphP[ParticleNum[j]].Density,0,&(SphP[ParticleNum[j]].Ne),ParticleNum[j],1);
				//printf("IonRate Tini Nphotons Prandom mol %e %e %e %e \n", IonRate[j], T, Tfin, tcooling*All.UnitTime_in_s/SEC_PER_YEAR);
				P[i].N_photons -= IonRate[j];
				/*
				fbtime = (1 << P[ParticleNum[j]].TimeBin) * All.Timebase_interval / All.cf_hubble_a;
				fbtime = DMAX(fbtime,fabs(tcooling));
				printf("IonRate Tini Nphotons Prandom mol %e %e \n", fbtime, tcooling);*/
			}	
			else 
			{
				double Prandom = get_random_number(ThisTask);
				if (IonRate[j]/P[i].N_photons > Prandom)
				{
					SphP[ParticleNum[j]].InternalEnergy = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
					SphP[ParticleNum[j]].InternalEnergyPred = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
					SphP[ParticleNum[j]].HIIregion = 1;
					P[i].N_photons -= IonRate[j];
					/*tcooling= CallGrackle(SphP[ParticleNum[j]].InternalEnergy,SphP[ParticleNum[j]].Density,0,&(SphP[ParticleNum[j]].Ne),ParticleNum[j],1);
					fbtime = (P[ParticleNum[j]].TimeBin ? (1 << P[ParticleNum[j]].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
					fbtime = DMAX(fbtime,fabs(tcooling));*/
				}	
			}
			P[i].Feedback_timestep = DMIN(P[i].Feedback_timestep,fbtime);
			if (P[i].N_photons <= 0) break;
		}
		
		//printf("Final number of photons = %e \n", P[i].N_photons);
		//exit(0);
	}	
	
    free(IonRate);
    IonRate=NULL;
    free(Distance);
    Distance=NULL;
    free(Tini);
    Tini=NULL;
    free(ParticleNum);
    ParticleNum=NULL;
    //free(Tag_HIIregion);
    //Tag_HIIregion=NULL;

}	
	

#endif //photionization