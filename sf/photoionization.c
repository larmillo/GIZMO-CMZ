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
static inline void downheap2 (double* data1, double* data2, double* data3, int* data4, int* data5, const size_t N, size_t k) {
    
    double v1 = data1[k];
    double v2 = data2[k];
	double v3 = data3[k];
	int v4 = data4[k];
	int v5 = data5[k];

    while (k<=N/2) {
        size_t j = 2 * k;
        if (j < N && data1[j] < data1[(j + 1)]) j++;
        if (!(v1 < data1[j]))  break;
        data1[k] = data1[j];
        data2[k] = data2[j];
		data3[k] = data3[j];
		data4[k] = data4[j];
		data5[k] = data5[j];
        k = j;
    }
    data1[k] = v1;
    data2[k] = v2;
    data3[k] = v3;
    data4[k] = v4;
    data5[k] = v5;
}


void sort (double *data1, double* data2, double* data3, int* data4, int* data5, const size_t n) {
    
    if (n==0) return;
    size_t N = n - 1;
    size_t k = N / 2;
    k++;     
    do {
        k--;
        downheap2(data1, data2, data3, data4, data5, N, k);
    }
    while (k > 0);

    while (N > 0) {
      /* first swap the elements */
      double tmp;
      
      tmp = data1[0];
      data1[0] = data1[N];
      data1[N] = tmp;

      tmp = data2[0];
      data2[0] = data2[N];
      data2[N] = tmp;
	  
      tmp = data3[0];
      data3[0] = data3[N];
      data3[N] = tmp;
	  
	  int tmp2;

      tmp2 = data4[0];
      data4[0] = data4[N];
      data4[N] = tmp2;
	  
      tmp2 = data5[0];
      data5[0] = data5[N];
      data5[N] = tmp2;

      N--;
      downheap2(data1, data2, data3, data4, data5, N, 0);
    }
}

void HII_region(void)
{
    double *Distance = (double*)malloc(N_gas*sizeof(double*));
	double *IonRate = (double*)malloc(N_gas*sizeof(double*));
	double *Tini = (double*)malloc(N_gas*sizeof(double*));
	int *ParticleNum = (int*)malloc(N_gas*sizeof(int*));
	int *Tag_HIIregion = (int*)malloc(N_gas*sizeof(int*));
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
		sort(Distance, IonRate, Tini, ParticleNum, Tag_HIIregion, N_gas);
		
		
		for(int j = 0; j < N_gas; j++) /* loop over the gas block */
		{
			
			if (Tag_HIIregion[j] == 1) continue; // The particle belongs to another HII region
			if (Tini[j] >= Tfin) continue; //Particle already ionized
			if (IonRate[j] <= P[i].N_photons) 
			{
				SphP[ParticleNum[j]].InternalEnergy = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
				SphP[ParticleNum[j]].InternalEnergyPred = BOLTZMANN*Tfin/((GAMMA-1)*molw_i*PROTONMASS)*All.UnitMass_in_g / All.UnitEnergy_in_cgs;
				SphP[ParticleNum[j]].HIIregion = 1;
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
    free(Tag_HIIregion);
    Tag_HIIregion=NULL;

}	
	

#endif //photionization