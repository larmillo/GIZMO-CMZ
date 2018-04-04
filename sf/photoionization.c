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
	int *Tag_HIIregion = (int*)malloc(N_gas*sizeof(int*));
	double ord1,ord2, ord5;
	int ord3,ord4;
	double Tfin = 1e4;
	double molw_i = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* assuming full ionization */
	
	//for(int j = 0; j < N_gas; j++) SphP[j].HIIregion=0; 
		
    for(int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 4) continue;
        if(P[i].Mass <= 0) continue;
		
		//get number of photons per second
		P[i].N_photons = slug_get_photometry_QH0(P[i].SlugOb);
		if(P[i].N_photons <= 0) continue;
		//printf("Initial number of photons = %e \n", P[i].N_photons);
		
		for(int j = 0; j < N_gas; j++) /* loop over the gas block */
		{
			ParticleNum[j] = j;
			Tag_HIIregion[j] = SphP[j].HIIregion; 
			double distx = P[i].Pos[0] - P[j].Pos[0];
			double disty = P[i].Pos[1] - P[j].Pos[1];
			double distz = P[i].Pos[2] - P[j].Pos[2];
			Distance[j] = sqrt(distx*distx+disty*disty+distz*distz);
			double Rhob = SphP[j].Density * All.UnitDensity_in_cgs;
			double Mb   = P[j].Mass * All.UnitMass_in_g;
			Tini[j] = CallGrackle(SphP[j].InternalEnergy,SphP[j].Density,0,&(SphP[j].Ne),j,2);
			double molw_n = Tini[j]*BOLTZMANN/(GAMMA-1)/(SphP[j].InternalEnergy*All.UnitEnergy_in_cgs/All.UnitMass_in_g)/PROTONMASS;
			double beta = 3e-13; //cm**3s*-1, cgs units
			IonRate[j] = HYDROGEN_MASSFRAC*beta*Tini[j]*Rhob*Mb/(2*PROTONMASS*PROTONMASS*molw_n*molw_i*Tfin); 
			//printf("IonRate Tini Mb Rhob %e %e %e %e \n", IonRate[j], Tini[j], Mb, molw);
		}	
		
		//Cycle to sort particles by increasing distance
        for (int j = 0; j < (N_gas - 1); j++)
		{
		    for (int k = 0; k < N_gas - 1 - j; k++)
		    {
		        if (Distance[k] > Distance[k+1])
		        {
		            ord1 = Distance[k+1];
		            Distance[k+1] = Distance[k];
		            Distance[k] = ord1;
					
		            ord2 = IonRate[k+1];
		            IonRate[k+1] = IonRate[k];
		            IonRate[k] = ord2;
					
		            ord5 = Tini[k+1];
		            Tini[k+1] = Tini[k];
		            Tini[k] = ord5;
					
		            ord3 = ParticleNum[k+1];
		            ParticleNum[k+1] = ParticleNum[k];
		            ParticleNum[k] = ord3;
					
		            ord4 = Tag_HIIregion[k+1];
		            Tag_HIIregion[k+1] = Tag_HIIregion[k];
		            Tag_HIIregion[k] = ord4;
		        }
		     }
		}
		
		
		for(int j = 0; j < N_gas; j++) /* loop over the gas block */
		{
			
			if (Tag_HIIregion[j] == 1) continue; // The particle belongs to another HII region
			if (Tini[j] >= Tfin) continue; //Particle already ionized
			if (IonRate[j] <= P[i].N_photons) 
			{
				SphP[ParticleNum[j]].InternalEnergy = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
				SphP[ParticleNum[j]].InternalEnergyPred = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
				SphP[ParticleNum[j]].HIIregion = 1;
				//double T = CallGrackle(SphP[ParticleNum[j]].InternalEnergy,SphP[ParticleNum[j]].Density,0,&(SphP[ParticleNum[j]].Ne),ParticleNum[j],2);
				//double tcooling= CallGrackle(SphP[ParticleNum[j]].InternalEnergy,SphP[ParticleNum[j]].Density,0,&(SphP[ParticleNum[j]].Ne),ParticleNum[j],1);
				//printf("IonRate Tini Nphotons Prandom mol %e %e %e %e \n", IonRate[j], T, Tfin, tcooling*All.UnitTime_in_s/SEC_PER_YEAR);
				P[i].N_photons -= IonRate[j];
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
				}	
			}
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
    free(Tag_HIIregion);
    Tag_HIIregion=NULL;

}	
	

#endif //photionization