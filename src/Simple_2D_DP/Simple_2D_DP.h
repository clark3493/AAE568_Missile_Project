#ifndef SIMPLE_2D_DP_H
#define SIMPLE_2D_DP_H

double *f(double x, double z, double xdot, double zdot, double u, double dt,
            double m, double T, double g);

double theta_md_fun(double m,
                    double T,
                    double g,
                    double x0,
                    double x_target,
                    double z0,
                    double z_target,
                    double theta);
                    
double calc_theta_md(double m,
                     double T,
                     double g,
                     double x0,
                     double x_target,
                     double z0,
                     double z_target,
                     double a0,
                     double b0,
                     double tol);
                    
#endif