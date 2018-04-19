#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csv.h"
#include "utils.h"

#define BUFF_SIZE 1024
#define MAX_FLD_SIZE 32
#define MAX_FLDS 10
#define epsilon .000001

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

double interp4(double **states, double *xq, int *m, double ****V)
{
  int i,j,k,l;
  int iup,jup,kup,lup;
  int idown,jdown,kdown,ldown;
  double t1,t2,t3,t4,N,cost,y=0;
  double den = 1;
  
  int **bound_inds = interp_bound_inds(states,xq,4,m);
  if( bound_inds == NULL )
    return -1;
  
  idown = bound_inds[0][0];
  iup   = bound_inds[0][1];
  
  jdown = bound_inds[1][0];
  jup   = bound_inds[1][1];
  
  kdown = bound_inds[2][0];
  kup   = bound_inds[2][1];
  
  ldown = bound_inds[3][0];
  lup   = bound_inds[3][1];
  
  for( i=0; i<4; i++ )
    den = den * (states[i][bound_inds[i][1]]-states[i][bound_inds[i][0]]);
  
  for( i=0; i<2; i++ ) {
    for( j=0; j<2; j++ ) {
      for( k=0; k<2; k++ ) {
        for( l=0; l<2; l++ ) {
          
          if( i == 0 ) { t1 = states[0][iup] - xq[0]; }
          else { t1 = xq[0] - states[0][idown]; }
          
          if( j == 0 ) { t2 = states[1][jup] - xq[1]; }
          else { t2 = xq[1] - states[1][jdown]; }
          
          if( k == 0 ) { t3 = states[2][kup] - xq[2]; }
          else { t3 = xq[2] - states[2][kdown]; }
          
          if( l == 0 ) { t4 = states[3][lup] - xq[3]; }
          else { t4 = xq[3] - states[3][ldown]; }
          
          N = t1*t2*t3*t4 / den;          
          cost = V[bound_inds[0][i]][bound_inds[1][j]][bound_inds[2][k]][bound_inds[3][l]];
          
          if( N > 1 ) {
            printf("  POSSIBLE INTERPOLATION ERROR\n");
            printf("    N=%f\n",N);
            printf("    i,j,k,l=%d,%d,%d,%d\n",i,j,k,l);
            printf("    den = %f\n",den);
            printf("    t1,t2,t3,t4=%f,%f,%f,%f\n",t1,t2,t3,t4);
            printf("    state bound = %f,%f,%f,%f\n",states[0][bound_inds[0][i]],states[1][bound_inds[1][j]],states[2][bound_inds[2][k]],states[3][bound_inds[3][l]]);
            printf("    state = %f,%f,%f,%f\n",xq[0],xq[1],xq[2],xq[3]);
            printf("    cost = %f\n",cost);
            printf("\n");
          }
          
          y = y + N*cost;
        }
      }
    }
  }
  for( i=0; i<4; i++ )
    free(bound_inds[i]);
  free(bound_inds);
  return y;
}

int **interp_bound_inds(double **states, double *xq, int n, int *m)
{
  int i,j, temp;
  int **bound_inds = (int **)malloc(n*sizeof(int *));
  if( bound_inds == NULL )
    return NULL;
  
  for( i=0; i<n; i++ ) {
    bound_inds[i] = (int *)malloc(2*sizeof(int));
    if( bound_inds[i] == NULL ) {
      for( j=0; j<i; j++ )
        free(bound_inds[j]);
      free(bound_inds);
      return NULL;
    }
    if( xq[i] <= states[i][0] ) {
      bound_inds[i][0] = 0;
      bound_inds[i][1] = 1;
    }
    else if( xq[i] >= states[i][m[i]-1] ) {
      bound_inds[i][0] = m[i]-2;
      bound_inds[i][1] = m[i]-1;
    }
    else {
      for( j=0; j<m[i]; j++ ) {
        if( states[i][j] > xq[i] ) {
          bound_inds[i][0] = j-1;
          bound_inds[i][1] = j;
          break;
        }
      }
    }
  }
  return bound_inds;
}

double **interp_bounds(double **states, double *xq, int n, int *m)
{
  /*
   *
   * states = n x m
   *    n = # variables
   *    m = # states for the i'th variable
   *
   * bounds = n x 2
   *    0th index is lower bound
   *    1st index is upper bound
   *
   */
  int i,j, temp;
  double **bounds = (double **)malloc(n*sizeof(double *));
  if( bounds == NULL )
    return NULL;
  
  for( i=0; i<n; i++ ) {
    bounds[i] = (double *)malloc(2*sizeof(double));
    if( bounds[i] == NULL ) {
      for( j=0; j<i; j++ )
        free(bounds[j]);
      free(bounds);
      return NULL;
    }
    if( xq[i] <= states[i][0] ) {
      bounds[i][0] = states[i][0];
      bounds[i][1] = states[i][1];
    }
    else if( xq[i] >= states[i][m[i]-1] ) {
      bounds[i][0] = states[i][m[i]-2];
      bounds[i][1] = states[i][m[i]-1];
    }
    else {
      for( j=0; j<m[i]; j++ ) {
        if( states[i][j] > xq[i] ) {
          bounds[i][0] = states[i][j-1];
          bounds[i][1] = states[i][j];
          break;
        }
      }
    }
  }
  
  return bounds;
}

void cost_from_file4(FILE *f, int *dims, double ****V) {
  int i,j,k,l;
  char buff[BUFF_SIZE];
  double cost;
  char **flds;
  
  fgets(buff,sizeof(buff),f); /* skip header line */
  
  for( i=0; i<dims[0]; i++ ) {
    for( j=0; j<dims[1]; j++ ) {
      for( k=0; k<dims[2]; k++ ) {
        for( l=0; l<dims[3]; l++ ) {
          fgets(buff,sizeof(buff),f);
          flds = parse_csv(buff);
          V[i][j][k][l] = atof(flds[4]);
          free_csv_line(flds);
        }
      }
    }
  }
}

/*
double cost_from_file(FILE *f, double *state, int n)
{
  int i, ln=0, found=0;
  char buff[BUFF_SIZE];
  double cost, tmp;
  char **flds;
  
  while( !feof(f) ) {
    fgets(buff,sizeof(buff),f);
    ln++;
    if( found )
      break;
    if( ln != 1 ) {
      flds = parse_csv(buff);
      found=1;
      for( i=0; i<n; i++ ) {
        tmp = atof(flds[i]);
        if( !fequal(tmp,state[i]) ) {
          found=0;
          break;
        }
      }
      if ( found ) {
        cost = atof(flds[4]);
      }
      free_csv_line(flds);
    }
  }
  if( !found ) {
    printf("state not found\n");
    printf("state: %f,%f,%f,%f\n",state[0],state[1],state[2],state[3]);
    return -1;
  }
  
  return cost;
}
*/

double ****p4(int *size) {
  int i,j,k;
  double ****p = (double ****)malloc(size[0]*sizeof(double***));
  if( p == NULL )
    return NULL;
  
  for( i=0; i<size[0]; i++ ) {
    p[i] = (double ***)malloc(size[1]*sizeof(double **));
    if( p[i] == NULL )
      return NULL;
    for( j=0; j<size[1]; j++ ) {
      p[i][j] = (double **)malloc(size[2]*sizeof(double *));
      if( p[i][j] == NULL )
        return NULL;
      for( k=0; k<size[2]; k++ ) {
        p[i][j][k] = (double *)malloc(size[3]*sizeof(double));
        if( p[i][j][k] == NULL )
          return NULL;
      }
    }
  }
  return p;
}

int fequal(double a, double b)
{
  if( fabs(a-b) < epsilon )
    return 1;
  return 0;
}
  
  