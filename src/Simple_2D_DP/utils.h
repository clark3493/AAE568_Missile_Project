#ifndef UTILS_H
#define UTILS_H

typedef struct
{
  double *low_pts;
  double *low_vals;
  double *high_pts;
  double *high_vals;
  int n;
} InterpBound;

InterpBound newInterpBound(int size);

double sind(double x);
double cosd(double x);
double atand(double x);

double *quad_roots(double a, double b, double c);

double min(double *x, int n);
double max(double *x, int n);
double mean(double *x, int n);

double *linspace(double x0, double xf, int n);

double interp4(FILE *f, double **states, double *xq, int *m);
double **interp_bounds(double **states, double *xq, int n, int *m);

double cost_from_file(FILE *f, double *state, int n);

int fequal(double a, double b);

InterpBound futurecost_interp_vals(FILE *f, double *state);

#endif