#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "proto.h"
#include "halo_props.h"

double Ms_UM(double z) {

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
  snprintf(fname2, sizeof(fname2), "%s%s", file_head, halo_no_str);
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

double Fs_Guo(double z) {

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

/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */

/*! This function drifts all particles from the current time to the future:
 *  time0 - > time1
 *
 *  If there is no explicit tree construction in the following timestep, the
 *  tree nodes are also drifted and updated accordingly. Note: For periodic
 *  boundary conditions, the mapping of coordinates onto the interval
 *  [0,All.BoxSize] is only done before the domain decomposition, or for
 *  outputs to snapshot files.  This simplifies dynamic tree updates, and
 *  allows the domain decomposition to be carried out only every once in a
 *  while.
 */
void move_particles(int time0, int time1)
{
  int i, j, task;
  double dt_drift, dt_gravkick, dt_hydrokick, dt_entr;
  double t0, t1;
  const int root = 0;

  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      dt_drift = get_drift_factor(time0, time1);
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
  else
    {
      dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
    }
 
      
  int itask=0;
  int otask=0;
  for(i = 0; i < NumPart; i++){
	// DY update disk center
	if(P[i].ID==id0){
          double z_now = 1.0/All.Time - 1.0;
          Mdisk = Ms_UM(z_now);
          Mdisk *= HALO_NORM;
	  Fdisk = Fs_Guo(z_now);
	  Fdisk *= HALO_NORM;

	  printf("Update disk center, a,m,h,pos,vphys, mstar, fstar:  "
			 "%f %f %f   %f %f %f   %f %f %f  %f\n",
			 All.Time,P[i].Mass,All.ForceSoftening[P[i].Type],
			 P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],
			 P[i].Vel[0]/sqrt(All.Time),
			 P[i].Vel[1]/sqrt(All.Time),
			 P[i].Vel[2]/sqrt(All.Time),
                         Mdisk/1e10*All.HubbleParam,
			 Fdisk);

          cenx=P[i].Pos[0];
          ceny=P[i].Pos[1];
          cenz=P[i].Pos[2];
	  itask=ThisTask;
	}
  }
  // Need to call on all tasks! 
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(&itask, &otask, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#ifdef DOUBLEPRECISION
  MPI_Bcast(&cenx, 1, MPI_DOUBLE, otask, MPI_COMM_WORLD);
  MPI_Bcast(&ceny, 1, MPI_DOUBLE, otask, MPI_COMM_WORLD);
  MPI_Bcast(&cenz, 1, MPI_DOUBLE, otask, MPI_COMM_WORLD);
  MPI_Bcast(&Mdisk, 1, MPI_DOUBLE, otask, MPI_COMM_WORLD);
  MPI_Bcast(&Fdisk, 1, MPI_DOUBLE, otask, MPI_COMM_WORLD);
#else
  MPI_Bcast(&cenx, 1, MPI_FLOAT, otask, MPI_COMM_WORLD);
  MPI_Bcast(&ceny, 1, MPI_FLOAT, otask, MPI_COMM_WORLD);
  MPI_Bcast(&cenz, 1, MPI_FLOAT, otask, MPI_COMM_WORLD);
  MPI_Bcast(&Mdisk, 1, MPI_FLOAT, otask, MPI_COMM_WORLD);
  MPI_Bcast(&Fdisk, 1, MPI_FLOAT, otask, MPI_COMM_WORLD);
#endif
  MPI_Barrier(MPI_COMM_WORLD); 

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	P[i].Pos[j] += P[i].Vel[j] * dt_drift;

      if(P[i].Type == 0)
	{
#ifdef PMGRID
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] +=
	      (P[i].GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick;
#else
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick;
#endif
	  SphP[i].Density *= exp(-SphP[i].DivVel * dt_drift);
	  SphP[i].Hsml *= exp(0.333333333333 * SphP[i].DivVel * dt_drift);

	  if(SphP[i].Hsml < All.MinGasHsml)
	    SphP[i].Hsml = All.MinGasHsml;

	  dt_entr = (time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

	  SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
	}
    }

  /* if domain-decomp and tree are not going to be reconstructed, update dynamically.  */
  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      for(i = 0; i < Numnodestree; i++)
	for(j = 0; j < 3; j++)
	  Nodes[All.MaxPart + i].u.d.s[j] += Extnodes[All.MaxPart + i].vs[j] * dt_drift;

      force_update_len();

      force_update_pseudoparticles();
    }

  t1 = second();

  All.CPU_Predict += timediff(t0, t1);
}

/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize[j];

	while(P[i].Pos[j] >= boxsize[j])
	  P[i].Pos[j] -= boxsize[j];
      }
}
#endif
