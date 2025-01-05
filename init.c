#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "halo_props.h"
#include "allvars.h"
#include "proto.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/*! \file init.c
 *  \brief Code for initialisation of a simulation from initial conditions
 */

double MS_UM(double z) {

  FILE *myFile1, *myFile2;
  int HALO_NO = HALO_ID;
  char buffer[1024];
  if ((HALO_NO == 9749) || (HALO_NO == 9829)){
      sprintf(buffer, "%04d", HALO_NO);
  }
  else{
      sprintf(buffer, "%03d", HALO_NO);
  }
  char* halo_no_str = buffer;
  char fname1[1024];
  char fname2[1024];
  char* file_head = "/sdf/home/y/ycwang/Gadget2_disk/um_sfh/Halo";
  char* file_tail1 = "8_Z.txt";
  char* file_tail2 = "8_Mstar.txt";
  snprintf(fname1, sizeof(fname1), "%s%s", file_head, halo_no_str);
  snprintf(fname2,sizeof(fname2), "%s%s", file_head, halo_no_str);
  char* file_name1 = fname1;
  char* file_name2 = fname2;
  strcat(file_name1, file_tail1);
  strcat(file_name2, file_tail2);
  myFile1 = fopen(file_name1, "r");
  myFile2 = fopen(file_name2, "r");
  int snap = TOT_SNAP;
  
  double red[snap];
  double mstar[snap];
  int i; 
  for (i = 0; i < snap; i++){
      fscanf(myFile1, "%lf", &red[i]);
      fscanf(myFile2, "%lf", &mstar[i]);
  }

  double lgms_z;
  {
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_linear, snap);

    gsl_spline_init (spline, red, mstar, snap);
    lgms_z = gsl_spline_eval (spline, z, acc);

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

  fclose(myFile1);
  fclose(myFile2);
  return pow(10., lgms_z);
}

double FS_Guo(double z) {

  FILE *myFile1, *myFile2;
  char* file1 = "/sdf/home/y/ycwang/Gadget2_disk/um_sfh/Guo_Z.txt";
  char* file2 = "/sdf/home/y/ycwang/Gadget2_disk/um_sfh/Guo_Fstar.txt";
  myFile1 = fopen(file1, "r");
  myFile2 = fopen(file2, "r");
  int snap = 14;

  double red[snap];
  double fstar[snap];
  int i;
  for (i = 0; i < snap; i++){
      fscanf(myFile1, "%lf", &red[i]);
      fscanf(myFile2, "%lf", &fstar[i]);
  }

  double fs_z;
  {
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_linear, snap);

    gsl_spline_init (spline, red, fstar, snap);
    fs_z = gsl_spline_eval (spline, z, acc);

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

  fclose(myFile1);
  fclose(myFile2);
  return fs_z;
}

/*! This function reads the initial conditions, and allocates storage for the
 *  tree. Various variables of the particle data are initialised and An intial
 *  domain decomposition is performed. If SPH particles are present, the inial
 *  SPH smoothing lengths are determined.
 */
void init(void)
{
  int i, j;
  double a3;

  All.Time = All.TimeBegin;

  switch (All.ICFormat)
    {
    case 1:
#if (MAKEGLASS > 1)
      seed_glass();
#else
      read_ic(All.InitCondFile);
#endif
      break;
    case 2:
    case 3:
      read_ic(All.InitCondFile);
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;
  All.Ti_Current = 0;

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      a3 = All.Time * All.Time * All.Time;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      a3 = 1;
    }

  set_softenings();

// DY==============================================================
  //RW: Move PartType4 that's closest to halo center at A_START to halo center, convert it to Type 6
  id0=SINK_ID; 
  int i0=-1;
  for(i = 0; i < NumPart; i++){
     if(P[i].ID==id0){
        if(All.Time <= START_A){ // Do this only once !!!! 
           cenx=HALO_CEN_X;
           ceny=HALO_CEN_Y;
           cenz=HALO_CEN_Z;
           P[i].Pos[0]=cenx;
           P[i].Pos[1]=ceny;
           P[i].Pos[2]=cenz;
           double zz = 1./All.Time - 1.;
           Mdisk = MS_UM(zz); //Assign UM disk mass
	   Fdisk = FS_Guo(zz); //Assign Guo stellar fraction 
           // 
	   P[i].Vel[0]=HALO_VX/sqrt(All.Time);
	   P[i].Vel[1]=HALO_VY/sqrt(All.Time);
	   P[i].Vel[2]=HALO_VZ/sqrt(All.Time);
           P[i].Mass=0.015; //
           P[i].Type=6; //
        }
	else{
           cenx=P[i].Pos[0];
           ceny=P[i].Pos[1];
           cenz=P[i].Pos[2];
           P[i].Mass = 0.015; //+mdisk/1e10; //
	   P[i].Type = 6; // DY convert the type from snapshots
	}
        i0=i;
     }
  }
  if(i0>0){
     printf("Task: %d, x=%f, y=%f, z=%f\n",ThisTask,P[i0].Pos[0],P[i0].Pos[1],P[i0].Pos[2]);
     if(P[i0].Type!=6 || fabs(P[i0].Mass-0.015)>0.001){
         printf("CANNOT FIND THE SINK PARTICLE CORRECTLY!!!! \n");
         endrun(888); 
         return;
     }
  }
// DY==============================================================

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

  All.TotNumOfForces = 0;
  All.NumForcesSinceLastDomainDecomp = 0;

  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  if(All.ComovingIntegrationOn)	/*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
	for(j = 0; j < 3; j++)
	  P[i].Vel[j] *= sqrt(All.Time) * All.Time;
    }

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;
#ifdef PMGRID
      for(j = 0; j < 3; j++)
	P[i].GravPM[j] = 0;
#endif
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
      P[i].Potential = 0;
    }

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef FLEXSTEPS
  All.PresentMinStep = TIMEBASE;
  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
    }
#endif


  for(i = 0; i < N_gas; i++)	/* initialize sph_properties */
    {
      for(j = 0; j < 3; j++)
	{
	  SphP[i].VelPred[j] = P[i].Vel[j];
	  SphP[i].HydroAccel[j] = 0;
	}

      SphP[i].DtEntropy = 0;

      if(RestartFlag == 0)
	{
	  SphP[i].Hsml = 0;
	  SphP[i].Density = -1;
	}
    }

  ngb_treeallocate(MAX_NGB);

  force_treeallocate(All.TreeAllocFactor * All.MaxPart, All.MaxPart);

  All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;

  Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */

  domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */

  ngb_treebuild();		/* will build tree */

  setup_smoothinglengths();

  TreeReconstructFlag = 1;

  /* at this point, the entropy variable normally contains the 
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly. 
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
  if(header.flag_entropy_instead_u == 0)
    for(i = 0; i < N_gas; i++)
      SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif
}


/*! This routine computes the mass content of the box and compares it to the
 *  specified value of Omega-matter.  If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass += P[i].Mass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}



/*! This function is used to find an initial smoothing length for each SPH
 *  particle. It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;

  if(RestartFlag == 0)
    {

      for(i = 0; i < N_gas; i++)
	{
	  no = Father[i];

	  while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
	    {
	      p = Nodes[no].u.d.father;

	      if(p < 0)
		break;

	      no = p;
	    }
#ifndef TWODIMS
	  SphP[i].Hsml =
	    pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
	  SphP[i].Hsml =
	    pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
	}
    }

  density();
}


/*! If the code is run in glass-making mode, this function populates the
 *  simulation box with a Poisson sample of particles.
 */
#if (MAKEGLASS > 1)
void seed_glass(void)
{
  int i, k, n_for_this_task;
  double Range[3], LowerBound[3];
  double drandom, partmass;
  long long IDstart;

  All.TotNumPart = MAKEGLASS;
  partmass = All.Omega0 * (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G))
    * (All.BoxSize * All.BoxSize * All.BoxSize) / All.TotNumPart;

  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */

  allocate_memory();

  header.npartTotal[1] = All.TotNumPart;
  header.mass[1] = partmass;

  if(ThisTask == 0)
    {
      printf("\nGlass initialising\nPartMass= %g\n", partmass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
    }

  /* set the number of particles assigned locally to this task */
  n_for_this_task = All.TotNumPart / NTask;

  if(ThisTask == NTask - 1)
    n_for_this_task = All.TotNumPart - (NTask - 1) * n_for_this_task;

  NumPart = 0;
  IDstart = 1 + (All.TotNumPart / NTask) * ThisTask;

  /* split the temporal domain into Ntask slabs in z-direction */

  Range[0] = Range[1] = All.BoxSize;
  Range[2] = All.BoxSize / NTask;
  LowerBound[0] = LowerBound[1] = 0;
  LowerBound[2] = ThisTask * Range[2];

  srand48(ThisTask);

  for(i = 0; i < n_for_this_task; i++)
    {
      for(k = 0; k < 3; k++)
	{
	  drandom = drand48();

	  P[i].Pos[k] = LowerBound[k] + Range[k] * drandom;
	  P[i].Vel[k] = 0;
	}

      P[i].Mass = partmass;
      P[i].Type = 1;
      P[i].ID = IDstart + i;

      NumPart++;
    }
}
#endif
