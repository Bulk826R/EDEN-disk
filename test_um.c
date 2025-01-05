#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "halo_props.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/*
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
  printf("halo_no = %s \n", buffer);
  //return 1.;
   
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

  printf("File1 %s \n", file_name1);
  printf("File2 %s \n", file_name2);
  printf("Snap = %d \n", snap);
  
  double red[snap];
  double mstar[snap];
  int i;
  for (i = 0; i < snap; i++){
      fscanf(myFile1, "%lf", &red[i]);
      fscanf(myFile2, "%lf", &mstar[i]);
  }
  //return 1.;

  //for (i=0; i<snap; i++){
  //  printf("%f \n", red[i]);
  //}

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
  printf("lg %f", lgms_z);
  fclose(myFile1);
  fclose(myFile2);
  return pow(10., lgms_z);
 
}
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
  int snap = 13;
  printf("snap = %d\n", snap);
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

double Fs_Guo(double z) {

  FILE *myFile1, *myFile2;
  char* file1 = "/sdf/home/y/ycwang/Gadget2_disk/um_sfh/Guo_Z.txt";
  char* file2 = "/sdf/home/y/ycwang/Gadget2_disk/um_sfh/Guo_Fstar.txt";
  myFile1 = fopen(file1, "r");
  myFile2 = fopen(file2, "r");
  int snap = 13;

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


int main(void){
  double MM = 0.;
  double fs = 0.;
  MM = MS_UM(0.2);
  fs = Fs_Guo(0.2);
  printf("Mstar, fstar = %f %f \n", MM, fs);
}
