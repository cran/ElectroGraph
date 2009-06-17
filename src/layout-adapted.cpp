/*
######################################################################
#
# layout-adapted.cpp, adapted from the R package "network" plotting functions.
#
######################################################################
*/
 
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "layout-adapted.h"
//const int debug_mode_adapt = 1;

extern "C" {
	
void network_layout_adapted (double * true_dist,
	const int * nn_c,
	const int * niter_c,
	const double * coolexp_c,
	const int * dimension_c,

	const int * fruchterman_c,
	const double * egocenters,
	double * coords //nn-by-dimension
) {
	
	int nn = *nn_c, niter = *niter_c, dim = *dimension_c, fruchterman = *fruchterman_c;
	double coolexp = *coolexp_c, t_dist, maxmove, force;
	double maxdelta = (double)nn;

	int iteration, jj, kk, dd;
	int magic_power = 4;
	double *change, *diffs, *power_dist;
	change = (double *)R_alloc(nn * dim,sizeof(double)); //simultaneous.
	diffs = (double *)R_alloc(dim, sizeof(double));
	power_dist = (double *)R_alloc(nn*nn, sizeof(double));
	//true distances cubed. Guaranteed connected in this routine.
	for (jj = 0; jj < nn*nn; jj++) power_dist[jj] = pow(true_dist[jj],magic_power);

	for (iteration = niter; iteration >= 0; iteration--) {
		maxmove = maxdelta*pow(iteration/(double)niter,coolexp);		
	//	std::cout << "i:" << iteration << " t:" << maxmove << std::endl;

		for (jj = 0; jj<nn*dim; jj++) change[jj]=0; //clear changes.

		for (jj = 0; jj < (nn-1); jj++) for (kk = jj+1; kk < nn; kk++) {
			t_dist = 0; 
			for (dd = 0; dd < dim; dd++) {
				diffs[dd] = coords[jj + nn*dd] - coords[kk + nn*dd]; 
				t_dist += diffs[dd]*diffs[dd];
			}
			t_dist = sqrt(t_dist);

			for (dd = 0; dd < dim; dd++) diffs[dd] /= t_dist; //normalized components.
			if (fruchterman) {
				force = (pow(t_dist,magic_power-1)/power_dist[jj+nn*kk] - 1/t_dist); //all repulsive, attractive here. positive is net-attractive.
			} else { //kamada-kawai spring force
				force = (t_dist - true_dist[jj+nn*kk] - 1/t_dist)/true_dist[jj+nn*kk]/true_dist[jj+nn*kk];
			}
			force = force*egocenters[jj]*egocenters[kk];
			for (dd = 0; dd < dim; dd++) {
				change[jj+nn*dd] -= diffs[dd]*force;
				change[kk+nn*dd] += diffs[dd]*force;
			}
		}
		for (jj = 0; jj < nn; jj++) {
			t_dist = 0; for (dd = 0; dd < dim; dd++) t_dist += change[jj+nn*dd]*change[jj+nn*dd];
			t_dist = sqrt(t_dist);
			if (t_dist > maxmove) {
				t_dist = maxmove/t_dist;
				for (dd = 0; dd < dim; dd++) change[jj+dd*nn] *= t_dist;
			}
			for (dd = 0; dd < dim; dd++) coords[jj+dd*nn] += change[jj+dd*nn];
		}
	}
} //network_layout_fr_adapted


} //extern "C"



