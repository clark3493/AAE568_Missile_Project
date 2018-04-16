#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Simple_2D_DP.h"
#include "utils.h"

int main(void)
{
  /****************************
   *                          *
   *     USER INPUTS          *
   *                          *
   ****************************/
  
  /* initial conditions */
  const double x0 = 0;
  const double z0 = -2001;
  const double xdot0 = 100;
  const double zdot0 = 0;
  
  /* final conditions */
  const double x_target = 2000;
  const double z_target = -1;
  const double eta = 45;
  
  /* cost function weights */
  const double A = 1;
  const double B = 1;
  const double C = 1000;
  
  /* constant parameters */
  const double m = 100;
  const double T = 1000;
  const double g = 9.81;
  
  /* discretization parameters */
  const int Nx = 21;
  const int Nz = 21;
  const int Nxdot = 21;
  const int Nzdot = 21;
  
  const int Nu = 19;
  
  const int Ntau = 101;
  const int Ntf = 11;
  
  const double u_min = -90;
  const double u_max =  90;
  
  const double rmin_v = 0.75;
  const double rmax_v = 1.50;
  
  const double rmin_tf = 0.75;
  const double rmax_tf = 2.00;
  
  /* numerical solution parameters */
  const double md_tol = .001;
  const double MyInf = 999999999;
  
  /* Output parameters */
  const char *paramsfilename = "data\\params.csv";
  const char *finalcostfilename = "data\\finalcost.csv";
  const char *controlfilename = "data\\control.csv";
  const char *statefilename = "data\\state.csv";
  
  /****************************
   *                          *
   *     END USER INPUTS      *
   *                          *
   ****************************/
  /* declare various repeatedly needed variables */
  int ii,jj,kk,ll,mm,nn,tt,len1,len2;
  char tfile1[128];   /* temporary filenames */
  char tfile2[128];
  
  /* Create data storage file */
  FILE *fparams = fopen(paramsfilename,"w");
  fprintf(fparams,"INPUTS\n");
  fprintf(fparams,"x0,%f\n",x0);
  fprintf(fparams,"z0,%f\n",z0);
  fprintf(fparams,"xdot0,%f\n",xdot0);
  fprintf(fparams,"zdot0,%f\n",zdot0);
  fprintf(fparams,"x_target,%f\n",x_target);
  fprintf(fparams,"z_target,%f\n",z_target);
  fprintf(fparams,"eta,%f\n",eta);
  fprintf(fparams,"A,%f\n",A);
  fprintf(fparams,"B,%f\n",B);
  fprintf(fparams,"C,%f\n",C);
  fprintf(fparams,"m,%f\n",m);
  fprintf(fparams,"T,%f\n",T);
  fprintf(fparams,"g,%f\n",g);
  fprintf(fparams,"Nx,%d\n",Nx);
  fprintf(fparams,"Nz,%d\n",Nz);
  fprintf(fparams,"Nxdot,%d\n",Nxdot);
  fprintf(fparams,"Nzdot,%d\n",Nzdot);
  fprintf(fparams,"Nu,%d\n",Nu);
  fprintf(fparams,"Ntau,%d\n",Ntau);
  fprintf(fparams,"Ntf,%d\n",Ntf);
  fprintf(fparams,"u_min,%f\n",u_min);
  fprintf(fparams,"u_max,%f\n",u_max);
  fprintf(fparams,"rmin_v,%f\n",rmin_v);
  fprintf(fparams,"rmax_v,%f\n",rmax_v);
  fprintf(fparams,"rmin_tf,%f\n",rmin_tf);
  fprintf(fparams,"rmax_tf,%f\n",rmax_tf);
  fprintf(fparams,"md_tol,%f\n",md_tol);
  fprintf(fparams,"MyInf,%f\n\n",MyInf);
  
  
  /* ---- CALCULATE SOLUTION TO MINIMUM DISTANCE PROBLEM -------- */
  const double theta_md = calc_theta_md(m,T,g,x0,x_target,z0,z_target,-89,89,.000000001);
  if( theta_md == -999 ) {
    printf("Something went wrong\n");
  }
  
  const double V0 = sqrt(xdot0*xdot0 + zdot0*zdot0);
  const double eta_md = atand((z_target-z0)/(x_target-x0));
  
  double *tf_mdxv = quad_roots(0.5*T/m*cosd(theta_md),V0*cosd(eta_md),x0-x_target);
  const double tf_mdx = max(tf_mdxv,2);
  free(tf_mdxv);
  
  double *tf_mdzv = quad_roots(0.5*(g+T/m*sind(theta_md)),V0*sind(eta_md),z0-z_target);
  const double tf_mdz = max(tf_mdzv,2);
  free(tf_mdzv);
  
  double tf_md;  
  
  if( fabs(tf_mdx-tf_mdz) > md_tol ) {
    printf("Difference between tf solutions for x and z min distance problems is out of tolerance\n");
    printf("tf_mdz=%f,tf_mdx=%f\n",tf_mdz,tf_mdx);
    return 0;
  }
  else {
    tf_md = (tf_mdx + tf_mdz)/2;
  }
  
  const double xf_md = x0 + tf_md*V0*cosd(eta_md) + 0.5*tf_md*tf_md*T/m*cosd(theta_md);
  const double zf_md = z0 + tf_md*V0*sind(eta_md) + 0.5*tf_md*tf_md*(g+T/m*sind(theta_md));
  const double xdotf_md = V0*cosd(eta_md) + tf_md*T/m*cosd(theta_md);
  const double zdotf_md = V0*sind(eta_md) + tf_md*(g+T/m*sind(theta_md));
  
  if( fabs(xf_md-x_target) > md_tol ) {
    printf("Difference between xf and xtarget for min distance problem is out of tolerance\n");
    return 0;
  }
  else if( fabs(zf_md-z_target) > md_tol ) {
    printf("Difference between zf and ztarget for min distance problem is out of tolerance\n");
    return 0;
  }
  else if( fabs( atand(zdotf_md/xdotf_md)-eta_md ) > md_tol ) {
    printf("Difference between etaf_md/eta_md for min distance problem is out of tolerance\n");
    return 0;
  }
  
  fprintf(fparams,"MIN DISTANCE PROBLEM\n");
  fprintf(fparams,"theta_md,%f\n",theta_md);
  fprintf(fparams,"V0,%f\n",V0);
  fprintf(fparams,"eta_md,%f\n",eta_md);
  fprintf(fparams,"tf_mdx,%f\n",tf_mdx);
  fprintf(fparams,"tf_mdz,%f\n",tf_mdz);
  fprintf(fparams,"tf_md,%f\n",tf_md);
  fprintf(fparams,"xf_md,%f\n",xf_md);
  fprintf(fparams,"zf_md,%f\n",zf_md);
  fprintf(fparams,"xdotf_md,%f\n",xdotf_md);
  fprintf(fparams,"zdotf_md,%f\n\n",zdotf_md);
  
  /* ----- DETERMINE MIN AND MAX VALUES FOR DISCRETIZED RANGES ----- */
  double x_min = x0 - 100;
  double x_max = x_target + 100;
  
  double z_min = z0 - 100;
  double z_max = z_target + 100;
  
  double xdot_min = rmin_v * xdot0;
  double xdot_max = rmax_v * xdotf_md;
  
  double zdot_min = rmin_v * zdot0 - 100;
  double zdot_max = rmax_v * zdotf_md;
  
  double tf_min = rmin_tf * tf_md;
  double tf_max = rmax_tf * tf_md;
  
  /* ----- DISCRETIZE THE STATE, CONTROL, AND TIME VARIABLES ----- */
  double dx = (x_max - x_min) / (Nx-1);
  double dz = (z_max - z_min) / (Nz-1);
  double dxdot = (xdot_max - xdot_min) / (Nxdot-1);
  double dzdot = (zdot_max - zdot_min) / (Nzdot-1);
  double du = (u_max - u_min) / (Nu-1);
  double dtau = 1. / (Ntau-1);
  double dtf = (tf_max - tf_min) / (Ntf-1);
  
  double *x_vec = linspace(x_min,x_max,Nx);
  double *z_vec = linspace(z_min,z_max,Nz);
  double *xdot_vec = linspace(xdot_min,xdot_max,Nxdot);
  double *zdot_vec = linspace(zdot_min,zdot_max,Nzdot);
  double *u_vec = linspace(u_min,u_max,Nu);
  double *tau_vec = linspace(0,1,Ntau);
  double *tf_vec = linspace(tf_min,tf_max,Ntf);
  
  fprintf(fparams,"DISCRETIZED VARIABLES");
  fprintf(fparams,"\nx"); for(ii=0;ii<Nx;ii++) fprintf(fparams,",%f",x_vec[ii]);
  fprintf(fparams,"\nz"); for(ii=0;ii<Nz;ii++) fprintf(fparams,",%f",z_vec[ii]);
  fprintf(fparams,"\nxdot"); for(ii=0;ii<Nxdot;ii++) fprintf(fparams,",%f",xdot_vec[ii]);
  fprintf(fparams,"\nzdot"); for(ii=0;ii<Nzdot;ii++) fprintf(fparams,",%f",zdot_vec[ii]);
  fprintf(fparams,"\nu"); for(ii=0;ii<Nu;ii++) fprintf(fparams,",%f",u_vec[ii]);
  fprintf(fparams,"\ntau"); for(ii=0;ii<Ntau;ii++) fprintf(fparams,",%f",tau_vec[ii]);
  fprintf(fparams,"\ntf"); for(ii=0;ii<Ntf;ii++) fprintf(fparams,",%f",tf_vec[ii]);
  fprintf(fparams,"\n\n");
  fclose(fparams);
  
  /* ----- CALCULATE TERMINAL COSTS FOR EVERY STATE COMBINATION ----- */
  double cost;
  FILE *fcost = fopen(finalcostfilename,"w");
  fprintf(fcost,"x,z,xdot,zdot,cost\n");
  for( ii=0; ii<Nx; ii++ ) {
    for( jj=0; jj<Nz; jj++ ) {
      for( kk=0; kk<Nxdot; kk++ ) {
        for( ll=0; ll<Nzdot; ll++ ) {
          cost = 
              A * (atand(zdot_vec[ll]/xdot_vec[kk]) - eta) *
                  (atand(zdot_vec[ll]/xdot_vec[kk]) - eta) +
              C * sqrt( (x_vec[ii]-x_target)*(x_vec[ii]-x_target) +
                        (z_vec[jj]-z_target)*(z_vec[jj]-z_target) );
          fprintf(fcost,"%f,%f,%f,%f,%f\n",
              x_vec[ii],z_vec[jj],xdot_vec[kk],zdot_vec[ll],cost);
        }
      }
    }
  }
  fclose(fcost);
  
  double state[] = {-50,-2000,80,90};
  FILE *test = fopen(finalcostfilename,"r");
  InterpBound ib = futurecost_interp_vals(test,state);
  printf("x_low=%f,x_low_cost=%f\n",ib.low_pts[0],ib.low_vals[0]);
  printf("x_high=%f,x_high_cost=%f\n",ib.high_pts[0],ib.high_vals[0]);
  
  /* ----- CALCULATE COST-TO-GO FOR EVERY STATE, AT EVERY TIME STEP -----
     ----- FOR EVERY POSSIBLE VALUE OF Tf ------------------------------- */
  /*
  time_t tstart = clock(), tstart_local;
  double tf,dt,L;
  
  double *x_next = (double *)malloc(Nu*sizeof(double));
  double *z_next = (double *)malloc(Nu*sizeof(double));
  double *xdot_next = (double *)malloc(Nu*sizeof(double));
  double *zdot_next = (double *)malloc(Nu*sizeof(double));
  double *state = (double *)malloc(4*sizeof(double));
  
  FILE *ffuturecost;
  
  for( int itf=0; itf<Ntf; itf++ ) {
    tf = tf_vec[itf];
    dt = tf/((double)Ntau-1);
    L = B * dt;
    
    printf("Performing solution for tf = %f\n",tf);
    tstart_local = clock();
    
    for( ii=0; ii<Nx; ii++ ) {
      for( jj=0; jj<Nz; jj++ ) {
        for( kk=0; kk<Nxdot; kk++ ) {
          for( ll=0; ll<Nzdot; ll++ ) {
            /* get next state for every possible control */
            /*for( mm=0; mm<Nu; mm++ ) {
              state = f(x_vec[ii],z_vec[jj],xdot_vec[kk],zdot_vec[ll],u_vec[mm],
                          dt,m,T,g);
              x_next[mm] = state[0];
              z_next[mm] = state[1];
              xdot_next[mm] = state[2];
              zdot_next[mm] = state[3];
            }*/
            /* Get the optimal future cost for every time step */
            /*for( tt=Ntau-2; tt>0; tt-- ) {*/
              /* get future cost file for reading
                 get current cost file for writing */
              /*memset(tfile1,0,sizeof(tfile1));
              memset(tfile2,0,sizeof(tfile2));
              
              sprintf(tfile2,"data\\cost_tf%d_tau%d.csv",itf,tt);
              if( tt == Ntau - 1 ) {
                strcpy(tfile1,"data\\finalcost.csv");
              }
              else {
                sprintf(tfile1,"data\\cost_tf%d_tau%d.csv",itf,tt+1);
              }
              ffuturecost = fopen(tfile1,"r");
              
              if( ii == 0 && jj == 0 && kk == 0 & ll == 0 ) {
                fcost = fopen(tfile2,"w");
              }
              else {
                fcost = fopen(tfile2,"a");
              }
              
              
              fclose(ffuturecost);
              fclose(fcost);
            }
          }
        }
      }
    }
  }
  */
  
  return 1;  
}

double *f(double x, double z, double xdot, double zdot, double u, double dt,
            double m, double T, double g)
{
  double *xout = (double *)malloc(4*sizeof(double));
  
  xout[0] = x + dt * xdot;
  xout[1] = z + dt * zdot;
  xout[2] = xdot + dt * (  T/m * cosd(u) );
  xout[3] = zdot + dt * ( -T/m * sind(u) + g );
  
  return xout;
}

double theta_md_fun(double m,
                    double T,
                    double g,
                    double x0,
                    double x_target,
                    double z0,
                    double z_target,
                    double theta)
{
  return (m*g + T*sind(theta))/(T*cosd(theta)) - (z_target - z0)/(x_target - x0);
}

double calc_theta_md(double m,
                     double T,
                     double g,
                     double x0,
                     double x_target,
                     double z0,
                     double z_target,
                     double a0,
                     double b0,
                     double tol)
{
  double err;
  double ta0 = theta_md_fun(m,T,g,x0,x_target,z0,z_target,a0);
  double tb0 = theta_md_fun(m,T,g,x0,x_target,z0,z_target,b0);
  
  if ( (ta0>0 && tb0 >0) || (ta0<0 && tb0<0) ) {
    printf("Input limits to calc_theta_md have the same output sign\n");
    printf("  ta = %f\n",ta0);
    printf("  tb = %f\n",tb0);
    err = -999;
    return err;
  }
  
  int i = 0;
  double an = a0, bn = b0;
  double theta = (an+bn)/2;
  double r_an,r_theta;
  
  while( fabs(bn-an)/2 > tol ) {
    i = i + 1;
    theta = (an+bn)/2;
    r_an = theta_md_fun(m,T,g,x0,x_target,z0,z_target,an);
    r_theta = theta_md_fun(m,T,g,x0,x_target,z0,z_target,theta);
    if( r_an>0 && r_theta>0 ) {
      an = theta;
    }
    else if( r_an<0 && r_theta>0 ) {
      bn = theta;
    }
    else if( r_an>0 && r_theta<0 ) {
      bn = theta;
    }
    else if( r_an<0 && r_theta<0 ) {
      an = theta;
    }
    else if( r_an == 0 ) {
      return an;
    }
    else if( r_theta == 0 ) {
      return theta;
    }
    else {
      err = -999;
      return err;
    }
  }
  
  return theta;
}