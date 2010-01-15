/*
######################################################################
#
# layout-adapted-2.cpp, adapted from the R package "network" plotting functions. This routine shows the differential force method in order to reach an earlier stopping time.
# 1-12-10
#
######################################################################
*/
 
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "layout-adapted.h"
const int debug_mode_adapt = 0;

extern "C" {
	
void network_layout_adapted (double * true_dist,
//void network_layout_adapted_new (double * true_dist,
	const int * nn_c,
	const int * niter_c,
	const double * coolexp_c,
	const int * dimension_c,

	const int * fruchterman_c,
	const double * egocenters,
	double * coords //nn-by-dimension
) {
	
	int nn = *nn_c, niter = *niter_c, dim = *dimension_c, fruchterman = *fruchterman_c;
	double coolexp = *coolexp_c, t_dist, maxmove, force, total_force, total_force2;
	double maxdelta = (double)nn;

	int iteration, jj, kk, dd;
	int magic_power = 4;
	double *change, *diffs, *power_dist;

	double mindist = 1; for (jj = 1; jj<nn*nn; jj++) {if (true_dist[jj]>0) mindist = (mindist<true_dist[jj]?mindist:true_dist[jj]);}

	double interval = 1e-4*mindist;
	double threshold = 1e-2*mindist; //change this!
	double tiniest = 1e-10*mindist;
	double biggest_move;

	change = (double *)R_alloc(dim,sizeof(double)); //simultaneous.
	diffs = (double *)R_alloc(dim, sizeof(double));
	power_dist = (double *)R_alloc(nn*nn, sizeof(double));
	//true distances cubed. Guaranteed connected in this routine.
	for (jj = 0; jj < nn*nn; jj++) power_dist[jj] = pow(true_dist[jj],magic_power);

	if (debug_mode_adapt) for (jj=0; jj < nn; jj++) {
		for (dd = 0; dd < dim; dd++) {
			std::cout << coords[jj+dd*nn] << " ";
		}
		std::cout << std::endl;
	}

	for (iteration = niter; iteration >= 0; iteration--) {
		maxmove = maxdelta*pow(iteration/(double)niter,coolexp);		
		biggest_move = 0;
	//	if (debug_mode_adapt) std::cout << "i:" << iteration << " t:" << maxmove << " c1:" << coords[0] << std::endl;

		for (jj = 0; jj < nn; jj++) {
			for (kk = 0; kk < dim; kk++) {change[kk]=0;}
			for (kk = 0; kk < nn; kk++) if (kk != jj) {
				t_dist = 0; 
				for (dd = 0; dd < dim; dd++) {
					diffs[dd] = coords[jj + nn*dd] - coords[kk + nn*dd]; 
					t_dist += diffs[dd]*diffs[dd];
				}
				t_dist = sqrt(t_dist);
	
				for (dd = 0; dd < dim; dd++) diffs[dd] /= t_dist; //normalized components.
				if (fruchterman) {
					force = (pow(t_dist,magic_power-1)/power_dist[jj+nn*kk] - 1/t_dist); //all repulsive, 	attractive here. positive is net-attractive.
				} else { //kamada-kawai spring force
					force = (t_dist - true_dist[jj+nn*kk])/true_dist[jj+nn*kk]/true_dist[jj+nn*kk]; // - 1/t_dist
				}
			//	force = force*egocenters[jj]*egocenters[kk];

				for (dd = 0; dd < dim; dd++) {
					change[dd] -= diffs[dd]*force;
				}
			}

			// move it a little bit.
			t_dist = 0; for (dd = 0; dd < dim; dd++) t_dist += change[dd]*change[dd];
			t_dist = sqrt(t_dist);
			total_force = t_dist;

			if (debug_mode_adapt) {
				std::cout << "node:" << jj << ": total_force = " << total_force << ", ";
				for (dd = 0; dd < dim; dd++) std::cout << coords[jj+dd*nn] << " " << change[dd] << " ";
				std::cout << std::endl;
			}  

			if (total_force > tiniest) {
				for (dd = 0; dd < dim; dd++) coords[jj+dd*nn] += change[dd]/t_dist*interval;
			
				//second go-around.
				for (kk = 0; kk < dim; kk++) {change[kk]=0;}
				for (kk = 0; kk < nn; kk++) if (kk != jj) {
					t_dist = 0; 
					for (dd = 0; dd < dim; dd++) {	
						diffs[dd] = coords[jj + nn*dd] - coords[kk + nn*dd]; 
						t_dist += diffs[dd]*diffs[dd];
					}
					t_dist = sqrt(t_dist);
		
					for (dd = 0; dd < dim; dd++) diffs[dd] /= t_dist; //normalized components.
					if (fruchterman) {
						force = (pow(t_dist,magic_power-1)/power_dist[jj+nn*kk] - 1/t_dist); //all repulsive, 	attractive 	here. positive is net-attractive.
					} else { //kamada-kawai spring force
						force = (t_dist - true_dist[jj+nn*kk])/true_dist[jj+nn*kk]/true_dist[jj+nn*kk]; // - 1/t_dist
					}
				//	force = force*egocenters[jj]*egocenters[kk];

					for (dd = 0; dd < dim; dd++) {
						change[dd] -= diffs[dd]*force;
					}
				}

				t_dist = 0; for (dd = 0; dd < dim; dd++) t_dist += change[dd]*change[dd];
				t_dist = sqrt(t_dist);
				total_force2 = t_dist;
				for (dd = 0; dd < dim; dd++) change[dd] /= t_dist; // unit vector again.
				t_dist = interval*total_force2/(total_force-total_force2)*(total_force > total_force2);

				if (debug_mode_adapt) {
					std::cout << "after: node " << jj << ": total_force = " << total_force2 << ", t_dist = " << t_dist << "; ";
					for (dd = 0; dd < dim; dd++) std::cout << coords[jj+dd*nn] << " " << change[dd] << " ";
					std::cout << std::endl;
				}  


				biggest_move = (t_dist>biggest_move?t_dist:biggest_move);
				if (t_dist > maxmove) {t_dist = maxmove/t_dist;}
				
				for (dd = 0; dd < dim; dd++) coords[jj+dd*nn] += change[dd]*t_dist;
			} //if minimai force.


		}
		if (debug_mode_adapt) std::cout << ", bigmove: " << biggest_move << ", threshold: " << threshold << std::endl;
		if (biggest_move < threshold) break;//iteration = 0;

	}

	if (debug_mode_adapt) for (jj=0; jj < nn; jj++) {
		for (dd = 0; dd < dim; dd++) {
			std::cout << coords[jj+dd*nn] << " ";
		}
		std::cout << std::endl;
	}

} //network_layout_fr_adapted


} //extern "C"



