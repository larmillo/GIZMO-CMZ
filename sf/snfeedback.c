#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif
#ifdef PTHREADS_NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

#ifdef SN_FEEDBACK

void SNproduction(void)
{
	double time_cluster, Cur_stellar_mass;
	int Num_supernovae_tot;
	
    int i;
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].Type != 4) continue;
        if(P[i].Mass <= 0) continue;
		
		//advance slug object in time
		time_cluster = All.Time - P[i].StellarAge;
		time_cluster *= (All.UnitTime_in_s / SEC_PER_YEAR);
		slug_advance(P[i].SlugOb,time_cluster);
	
		//get number of SNe
		P[i].Nsn_timestep = 0;
		Num_supernovae_tot = slug_get_stoch_sn(P[i].SlugOb);
		P[i].Nsn_timestep = Num_supernovae_tot - P[i].Nsn_tot; //Nsn_tot(current_time) - Nsn_tot(current_time-dt)
		if(P[i].Nsn_timestep > 0) printf("SN explosion %d %d %e %e \n", P[i].ID, P[i].Nsn_timestep, P[i].DensAroundStar* All.UnitDensity_in_cgs, P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]);
		P[i].Nsn_tot = Num_supernovae_tot; //update of cumulative number of SNe
	
		//get mass of ejecta
		P[i].Mej = 0.;
		Cur_stellar_mass = slug_get_stellar_mass(P[i].SlugOb);
		Cur_stellar_mass /= (All.UnitMass_in_g / SOLAR_MASS); //solar masses
		P[i].Mej = P[i].SlugMass - Cur_stellar_mass;
		if (P[i].Mej < 0.0 && fabs(P[i].Mej) > Cur_stellar_mass*1e-4)
			printf("Warning: Large fluctuation! relatTime = %e, |dm|/m = %e, prev_stellar_mass = %e, curr_stellar_mass = %e, Num = %d\n", time_cluster, fabs(P[i].Mej)/Cur_stellar_mass, P[i].SlugMass, Cur_stellar_mass, P[i].ID);
    	if (P[i].Mej < 0.0) P[i].Mej = 0.0;
		P[i].SlugMass = Cur_stellar_mass;
		
		//get number of photons per second
	}	
}

void Check_conservation(void)
{
	int i;

	
	for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
		double Ej = 1.e51 * P[i].Nsn_timestep / All.UnitEnergy_in_cgs; // energy per event
		double vj = sqrt(2*Ej*P[i].Mej);
		
        if(P[i].Type != 4) continue;
        if(P[i].Mass <= 0) continue;
	    if((P[i].Nsn_timestep<=0)||(P[i].DensAroundStar<=0)) continue;
		
		if (P[i].tocons[0]==0 || P[i].tocons[0]<0 || fabs(P[i].tocons[0]-P[i].Mej)>1e-8) 
		{
		printf("ERROR IN MASS CONSERVATION!!! %e %e %e\n", P[i].tocons[0], P[i].Mej, fabs(P[i].tocons[0]-P[i].Mej));
		exit(0);
		}
		if (P[i].tocons[1]==0 || P[i].tocons[1]<0 || fabs(P[i].tocons[1]-vj)>1e-8) 
		{
		printf("ERROR IN MOMENTUM CONSERVATION!!! %e %e %e\n", P[i].tocons[1], vj, fabs(P[i].tocons[1]-vj));
		exit(0);
		}
		if (fabs(P[i].tocons[2])/P[i].tocons[1] > 1e-8) 
		{
		printf("ERROR IN TOTAL MOMENTUM CONSERVATION!!! %e %e \n", P[i].tocons[2], fabs(P[i].tocons[2])/P[i].tocons[1]);
		exit(0);
		}
		if (P[i].tocons[3]==0 || P[i].tocons[3]<0 || fabs(P[i].tocons[3]-Ej)>1e-8) 
		{
		printf("ERROR IN ENERGY CONSERVATION!!! %e %e %e\n", P[i].tocons[3], Ej, fabs(P[i].tocons[3]-Ej));
		exit(0);
		}
	}	
}		

struct FBdata_in
{
    MyDouble Pos[3];
    MyDouble Hsml;
	MyDouble wb_tot;
	MyDouble Vel[3];
	int Nsn;
    MyDouble Mej;
	MyDouble Eej;
	MyDouble Vej;
	
    int NodeList[NODELISTLENGTH];
}
*FBDataIn, *FBDataGet;
           
/* define properties to be injected */
void particle2in_FB(struct FBdata_in *in, int i);
void particle2in_FB(struct FBdata_in *in, int i)
{
    if((P[i].Nsn_timestep<=0)||(P[i].DensAroundStar<=0)) {in->Mej=0; return;} // trap for no sne
    int k;
    for(k=0; k<3; k++) {in->Pos[k] = P[i].Pos[k];}
	for(k=0; k<3; k++) {in->Vel[k] = P[i].Vel[k];}
	in->wb_tot = P[i].wb_tot;
    in->Hsml = PPP[i].Hsml;
	in->Nsn = P[i].Nsn_timestep;
    in->Mej = P[i].Mej;
    in->Eej = 1.e51 * P[i].Nsn_timestep / All.UnitEnergy_in_cgs; // energy per event
	in->Vej = sqrt(2*in->Eej/in->Mej);
}



struct FBdata_out
{
    MyDouble M_coupled;
	MyDouble p_coupled;
	MyDouble pvec_coupled;
	MyDouble E_coupled;
}
*FBDataResult, *FBDataOut;

void out2particle_FB(struct FBdata_out *out, int i, int mode);
void out2particle_FB(struct FBdata_out *out, int i, int mode)
{
    P[i].Mass -= out->M_coupled;

    if((P[i].Mass<0)||(isnan(P[i].Mass))) {P[i].Mass=0;}
	
	//quantity to verify conservation in the star frame
	if (mode == 0) 
	{
		P[i].tocons[0] = out->M_coupled; //mass of ejecta
		P[i].tocons[1] = out->p_coupled; //momentum of ejecta
		P[i].tocons[2] = out->pvec_coupled; //total momentum = 0
		P[i].tocons[3] = out->E_coupled; //energy of ejecta
	}	
	
	if (mode == 1) 
	{
		P[i].tocons[0] += out->M_coupled; //mass of ejecta
		P[i].tocons[1] += out->p_coupled; //momentum of ejecta
		P[i].tocons[2] += out->pvec_coupled; //total momentum = 0
		P[i].tocons[3] += out->E_coupled; //energy of ejecta
	}		
}


struct kernel_FB
{
    double	dp[3],r,wk_i, wk_j, dwk_i, dwk_j,h_i;
};


int FB_evaluate_active_check(int i)
{
    if(P[i].Type != 4) return 0;
    if(P[i].Mass <= 0) return 0;
    if(PPP[i].Hsml <= 0) return 0;
    if(PPP[i].NumNgb <= 0) return 0;
    if(P[i].Nsn_timestep>0) {return 1;}
    return 0;
}


void *FB_evaluate_primary(void *p)
{
#define CONDITION_FOR_EVALUATION if(FB_evaluate_active_check(i)==1)
#define EVALUATION_CALL FB_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist)
#include "../system/code_block_primary_loop_evaluation.h"
#undef CONDITION_FOR_EVALUATION
#undef EVALUATION_CALL
}
void *FB_evaluate_secondary(void *p)
{
#define EVALUATION_CALL FB_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
#include "../system/code_block_secondary_loop_evaluation.h"
#undef EVALUATION_CALL
}


void sn_feedback_calc(void)
{
    int j, k, ngrp, ndone, ndone_flag;
    int recvTask, place;
    int save_NextParticle;
    long long n_exported = 0;
    
    /* allocate buffers to arrange communication */
    long long NTaskTimesNumPart;
    NTaskTimesNumPart = maxThreads * NumPart;
    Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
    size_t MyBufferSize = All.BufferSize;
    All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct FBdata_in) + sizeof(struct FBdata_out) + sizemax(sizeof(struct FBdata_in),sizeof(struct FBdata_out))));
    DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
    DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
    
    NextParticle = FirstActiveParticle;	/* begin with this index */
    do
    {
        
        BufferFullFlag = 0;
        Nexport = 0;
        save_NextParticle = NextParticle;
        
        for(j = 0; j < NTask; j++)
        {
            Send_count[j] = 0;
            Exportflag[j] = -1;
        }
        
        /* do local particles and prepare export list */
#ifdef PTHREADS_NUM_THREADS
        pthread_t mythreads[PTHREADS_NUM_THREADS - 1];
        int threadid[PTHREADS_NUM_THREADS - 1];
        pthread_attr_t attr;
        
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        pthread_mutex_init(&mutex_nexport, NULL);
        pthread_mutex_init(&mutex_partnodedrift, NULL);
        
        TimerFlag = 0;
        
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
        {
            threadid[j] = j + 1;
            pthread_create(&mythreads[j], &attr, FB_evaluate_primary, &threadid[j]);
        }
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            FB_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
        }
        
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_join(mythreads[j], NULL);
#endif
        
        
        if(BufferFullFlag)
        {
            int last_nextparticle = NextParticle;
            
            NextParticle = save_NextParticle;
            
            while(NextParticle >= 0)
            {
                if(NextParticle == last_nextparticle)
                    break;
                
                if(ProcessedFlag[NextParticle] != 1)
                    break;
                
                ProcessedFlag[NextParticle] = 2;
                
                NextParticle = NextActiveParticle[NextParticle];
            }
            
            if(NextParticle == save_NextParticle)
            {
                /* in this case, the buffer is too small to process even a single particle */
                endrun(116608);
            }
            
            int new_export = 0;
            
            for(j = 0, k = 0; j < Nexport; j++)
                if(ProcessedFlag[DataIndexTable[j].Index] != 2)
                {
                    if(k < j + 1)
                        k = j + 1;
                    
                    for(; k < Nexport; k++)
                        if(ProcessedFlag[DataIndexTable[k].Index] == 2)
                        {
                            int old_index = DataIndexTable[j].Index;
                            
                            DataIndexTable[j] = DataIndexTable[k];
                            DataNodeList[j] = DataNodeList[k];
                            DataIndexTable[j].IndexGet = j;
                            new_export++;
                            
                            DataIndexTable[k].Index = old_index;
                            k++;
                            break;
                        }
                }
                else
                    new_export++;
            
            Nexport = new_export;
            
        }
        
        n_exported += Nexport;
        
        for(j = 0; j < NTask; j++)
            Send_count[j] = 0;
        for(j = 0; j < Nexport; j++)
            Send_count[DataIndexTable[j].Task]++;
        
        MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);
        MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        
        for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
            Nimport += Recv_count[j];
            
            if(j > 0)
            {
                Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }
        
        FBDataGet = (struct FBdata_in *) mymalloc("FBDataGet", Nimport * sizeof(struct FBdata_in));
        FBDataIn = (struct FBdata_in *) mymalloc("FBDataIn", Nexport * sizeof(struct FBdata_in));
        
        /* prepare particle data for export */
        
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            particle2in_FB(&FBDataIn[j], place);
            memcpy(FBDataIn[j].NodeList,
                   DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }
        
        /* exchange particle data */
        int TAG_TO_USE = TAG_FBLOOP_1A;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* get the particles */
                    MPI_Sendrecv(&FBDataIn[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct FBdata_in), MPI_BYTE,
                                 recvTask, TAG_TO_USE,
                                 &FBDataGet[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct FBdata_in), MPI_BYTE,
                                 recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
        
        myfree(FBDataIn);
        FBDataResult = (struct FBdata_out *) mymalloc("FBDataResult", Nimport * sizeof(struct FBdata_out));
        FBDataOut = (struct FBdata_out *) mymalloc("FBDataOut", Nexport * sizeof(struct FBdata_out));
        
        /* now do the particles that were sent to us */
        NextJ = 0;
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_create(&mythreads[j], &attr, FB_evaluate_secondary, &threadid[j]);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int mainthreadid = omp_get_thread_num();
#else
            int mainthreadid = 0;
#endif
            FB_evaluate_secondary(&mainthreadid);
        }
        
#ifdef PTHREADS_NUM_THREADS
        for(j = 0; j < PTHREADS_NUM_THREADS - 1; j++)
            pthread_join(mythreads[j], NULL);
        
        pthread_mutex_destroy(&mutex_partnodedrift);
        pthread_mutex_destroy(&mutex_nexport);
        pthread_attr_destroy(&attr);
#endif
        
        if(NextParticle < 0)
            ndone_flag = 1;
        else
            ndone_flag = 0;
        
        MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        /* get the result */
        TAG_TO_USE = TAG_FBLOOP_1B;
        for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
            recvTask = ThisTask ^ ngrp;
            if(recvTask < NTask)
            {
                if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                    /* send the results */
                    MPI_Sendrecv(&FBDataResult[Recv_offset[recvTask]],
                                 Recv_count[recvTask] * sizeof(struct FBdata_out),
                                 MPI_BYTE, recvTask, TAG_TO_USE,
                                 &FBDataOut[Send_offset[recvTask]],
                                 Send_count[recvTask] * sizeof(struct FBdata_out),
                                 MPI_BYTE, recvTask, TAG_TO_USE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
		
		
        
        /* add the result to the local particles */
        for(j = 0; j < Nexport; j++)
        {
            place = DataIndexTable[j].Index;
            out2particle_FB(&FBDataOut[j], place, 1);
        }
        myfree(FBDataOut);
        myfree(FBDataResult);
        myfree(FBDataGet);
    }
    while(ndone < NTask);
    
    myfree(DataNodeList);
    myfree(DataIndexTable);
    myfree(Ngblist);
}


int FB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
	double r2;
    struct kernel_FB kernel;
    struct FBdata_in local;
    struct FBdata_out out;
    memset(&out, 0, sizeof(struct FBdata_out));
    //kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1);
	
    out.M_coupled = out.p_coupled = out.pvec_coupled = out.E_coupled = 0.;
	
	/* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_FB(&local, target);} else {local = FBDataGet[target];}
    if(local.Nsn<=0) return 0; // no SNe for the master particle! nothing to do here //
    if(local.Hsml<=0) return 0; // zero-extent kernel, no particles //
    //h2 = local.Hsml*local.Hsml;
    //kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    kernel.h_i = local.Hsml;
    double h2_i = kernel.h_i*kernel.h_i;
		 
    /* Now start the weights computation for this particle */
    if(mode == 0)
    {
        startnode = All.MaxPart;	/* root node */
    }
    else
    {
        startnode = FBDataGet[target].NodeList[0];
        startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }
    while(startnode >= 0)
    {
        while(startnode >= 0)
        { 
			//find neighbors that can -mutually- see one another, not just single-directional searching here
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) return -1;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n];
                if(P[j].Type != 0) continue; // require a gas particle //
                if(P[j].Mass <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
#ifdef PERIODIC
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
#endif
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                double h_j = PPP[j].Hsml;
                if(r2 <= 0) continue; // same particle //
                if((r2 >= h2_i) && (r2 >= h_j * h_j)) continue; // outside kernel //
				if(r2 > 2 * All.UnitLength_in_cm * CM_PER_KPC) continue; //search radius not larger than 2 kpc //
                kernel.r = sqrt(r2);
				
				double wb_mod = 0;
				
				for(k=0;k<3;k++) 
				{
					SphP[j].wb[k] /= local.wb_tot; //real weight// 
					wb_mod += SphP[j].wb[k] * SphP[j].wb[k];
				}
				
				double dmass, dmom[3], mom_tot[3], dmom_tocons = 0, dmom_mod = 0, dmomrf_mod = 0;
				double dtot_en, tot_en, v2 = 0.;
				
				for(k=0;k<3;k++) {mom_tot[k] = P[j].Vel[k]*P[j].Mass;}
				for(k=0;k<3;k++) {v2 += P[j].Vel[k]*P[j].Vel[k];}
				tot_en = P[j].Mass * (SphP[j].InternalEnergy + 0.5 * v2);
				
				//calculation of terminal momentum
				double pt, nb, fZ, mu, Rcool;
				mu = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);	/* note: assuming NEUTRAL GAS */
				nb = SphP[j].Density * All.UnitDensity_in_cgs / (mu * 1.67e-24);
				if(P[j].Metallicity[0] < 0.01) fZ = 2;
				else fZ = pow(P[j].Metallicity[0],-0.14);	
				pt = 4.8e5 * pow(P[j].Nsn_timestep,13./14.) * pow(nb,-1./7.) * pow(fZ, 3./2.); //terminal momentum in solar mass * kms^-1
				pt *= SOLAR_MASS / All.UnitMass_in_g;
				pt *= 1e5 / All.UnitVelocity_in_cm_per_s;
				Rcool = 28.4 * pow(P[j].Nsn_timestep,2./7.) * pow(nb,-3./7.) * fZ; //pc
				Rcool *= CM_PER_PC / All.UnitLength_in_cm;
				//////////////////////////////////
				
				//INJECTION OF MASS//
				dmass = local.Mej * sqrt(wb_mod);
				SphP[j].Density *= (1 + dmass/P[j].Mass);
				P[j].Mass += dmass; 
				out.M_coupled += local.Mej * sqrt(wb_mod); 
				//END OF INJECTION OF MASS//
				
				//INJECTION OF MOMENTUM//
				for(k=0;k<3;k++)
				{
					dmom[k] = SphP[j].wb[k] * local.Vej * local.Mej;
					//printf("Check segno %e %e \n",P[j].Pos[k]-local.Pos[k], dmom[k]); 
					out.pvec_coupled += dmom[k];
					dmom_tocons += dmom[k]*dmom[k]; //to check conservation
					dmom[k] *= DMIN(sqrt(1+(P[j].Mass-dmass)/dmass), pt/(local.Vej * local.Mej)); //accounting for effect of not-solved Sedov-Taylor phase
					dmomrf_mod += dmom[k]*dmom[k]; //rest frame
					dmom[k] += dmass * local.Vel[k];
					mom_tot[k] += dmom[k];
					P[j].Vel[k] = mom_tot[k]/P[j].Mass;
					dmom_mod +=  dmom[k]*dmom[k]; //for energy calculation
				}
				out.p_coupled += sqrt(dmom_tocons);
				//END OF INJECTION OF MOMENTUM//
				
				//INJECTION OF ENERGY//
				dtot_en = local.Eej * sqrt(wb_mod) + 0.5/dmass * (dmom_mod - dmomrf_mod);
				out.E_coupled += local.Eej * sqrt(wb_mod);
				tot_en += dtot_en;
				v2 = 0;
				for(k=0;k<3;k++) {v2 += P[j].Vel[k]*P[j].Vel[k];}
				tot_en /= P[j].Mass;
				SphP[j].InternalEnergy = tot_en - 0.5 * v2;
				if (Rcool < kernel.r) SphP[j].InternalEnergy *= pow(kernel.r/Rcool,-6.5);
				SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;	
				//END OF INJECTION OF ENERGY//
						
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1)
        {
            listindex++;
            if(listindex < NODELISTLENGTH)
            {
                startnode = FBDataGet[target].NodeList[listindex];
                if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}	/* open it */
            }
        } // if(mode == 1)
    } // while(startnode >= 0)
    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_FB(&out, target, 0);} else {FBDataResult[target] = out;}
    return 0;
} // int FB_evaluate


#endif