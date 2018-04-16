#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

#define BUFF_SIZE 1024
#define MAX_FLD_SIZE 32
#define MAX_FLDS 10

#define PI 3.14159265359

InterpBound newInterpBound(int size)
{
  InterpBound ib;
  ib.n = size;
  ib.low_pts   = (double *)malloc(size*sizeof(double));
  ib.low_vals  = (double *)malloc(size*sizeof(double));
  ib.high_pts  = (double *)malloc(size*sizeof(double));
  ib.high_vals = (double *)malloc(size*sizeof(double));
  
  return ib;
}

double sind(double x)
{
  return sin(x * PI / 180 );
}

double cosd(double x)
{
  return cos(x * PI / 180 );
}

double atand(double x)
{
  return atan(x) * 180 / PI;
}

double *quad_roots(double a, double b, double c)
{
  double det;
  double *r = (double *)malloc(2*sizeof(double));
  
  det = b*b - 4*a*c;
  
  /* condition for real and different roots */
  if( det > 0 ) {
    r[0] = (-b+sqrt(det))/(2*a);
    r[1] = (-b-sqrt(det))/(2*a);
  }
  
  /* condition for real and equal roots */
  else if( det == 0 ) {
    r[0] = r[1] = -b/(2*a);
  }
  
  /* condition for non real roots */
  else {
    printf("Roots are imaginary\n");
    return NULL;
  }
  
  return r;  
}

double min(double *x, int n)
{
  double m = x[0];
  for( int i = 1; i < n; i++ ) {
    if( x[i] < x[i-1] )
      m = x[i];
  }
  return m;
}

double max(double *x, int n)
{
  double m = x[0];
  for( int i = 1; i < n; i++ ) {
    if( x[i] > x[i-1] )
      m = x[i];
  }
  return m;
}

double mean(double *x, int n)
{
  double sum = 0;
  for( int i = 0; i < n; i++ )
    sum  = sum + x[i];
  
  return sum / n;
}

double *linspace(double x0, double xf, int n)
{
  int i;
  
  double dx = (xf-x0)/((double)n-1);
  double *x = (double *)malloc(n*sizeof(double));
  x[0] = x0;
  for( i = 1; i < n; i++ ) {
    x[i] = x[i-1] + dx;
  }
    
  return x;
}

InterpBound futurecost_interp_vals(FILE *f, double *state)
{
  double l=-99999999, h=99999999;
  double x_low=l, x_low_cost, x_high=h, x_high_cost;
  double z_low=l, z_low_cost, z_high=h, z_high_cost;
  double xdot_low=l, xdot_low_cost, xdot_high=h, xdot_high_cost;
  double zdot_low=l, zdot_low_cost, zdot_high=h, zdot_high_cost;
  
  double x,z,xdot,zdot,cost;
  
  int init = 1,i,n=0;
  char buff[BUFF_SIZE];
  char **flds = (char **)malloc(MAX_FLDS*sizeof(char*));
  for(int c=0; c<MAX_FLDS; c++)
    flds[c] = (char *)malloc(MAX_FLD_SIZE*sizeof(char));
  char *fld;
  
  InterpBound ib = newInterpBound(4);
  
  while(!feof(f)) {
    i = 0; n++;
    if( n != 1 ) {
      fgets(buff,sizeof(buff),f);
      fld = strtok(buff,",");
      while( fld != NULL ) {
        memset(flds[i],'0',sizeof(flds[i]));
        strcpy(flds[i],fld);
        i++;
        fld = strtok(NULL,",");
      }
      
      x=atof(flds[0]);
      z=atof(flds[1]);
      xdot=atof(flds[2]);
      zdot=atof(flds[3]);
      cost=atof(flds[4]);
      
      if( x < x_high && x >= state[0] ) {
        x_high = x;
        x_high_cost = cost;
      }
      if( x > x_low && x <= state[0] ) {
        x_low = x;
        x_low_cost = cost;
      }
      if( z < z_high && z >= state[1] ) {
        z_high = z;
        z_high_cost = cost;
      }
      if( z > z_low && z <= state[1] ) {
        z_low = z;
        z_low_cost = cost;
      }
      if( xdot < xdot_high && xdot >= state[2] ) {
        xdot_high = xdot;
        xdot_high_cost = cost;
      }
      if( xdot > xdot_low && xdot <= state[2] ) {
        xdot_low = xdot;
        xdot_low_cost = cost;
      }
      if( zdot < zdot_high && zdot >= state[3] ) {
        zdot_high = zdot;
        zdot_high_cost = cost;
      }
      if( zdot > zdot_low && zdot <= state[3] ) {
        zdot_low = zdot;
        zdot_low_cost = cost;
      }
    }
  }
  
  ib.low_pts[0] = x_low;
  ib.low_vals[0] = x_low_cost;
  ib.high_pts[0] = x_high;
  ib.high_vals[0] = x_high_cost;
  
  ib.low_pts[1] = z_low;
  ib.low_vals[1] = z_low_cost;
  ib.high_pts[1] = z_high;
  ib.high_vals[1] = z_high_cost;
  
  ib.low_pts[2] = xdot_low;
  ib.low_vals[2] = xdot_low_cost;
  ib.high_pts[2] = xdot_high;
  ib.high_vals[2] = xdot_high_cost;
  
  ib.low_pts[3] = x_low;
  ib.low_vals[3] = x_low_cost;
  ib.high_pts[3] = x_high;
  ib.high_vals[3] = x_high_cost;
  
  for( int c=0; c< MAX_FLDS; c++ )
    free(flds[c]);
  free(flds);
  
  return ib;
}
  
  