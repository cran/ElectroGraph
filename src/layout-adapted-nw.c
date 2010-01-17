/*
######################################################################
#
# layout-adapted-nw.cpp, adapted from the R package "network" plotting functions. 
# This routine shows the differential force method in order to reach an earlier stopping time.
# This routine also considers connectivity to be the prime determinant of force.
# January 18, 2010
#
######################################################################
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
const int debug_mode_adapt = 0;
	
void network_layout_by_connectivity (double * connectivity,
	const int * nn_c,
	const int * niter_c,
	const int * dimension_c,
	const int * fruchterman_c,
	const double * egocenters,
	double * coords,
	const int * verbose_c 
) {
	
	int nn = *nn_c, niter = *niter_c, dim = *dimension_c, fruchterman = *fruchterman_c;
	int verbose = *verbose_c;
	double t_dist, force, total_force, total_force2;
	//double maxdelta = (double)nn;

	int iteration, jj, kk, dd;
	double *change, *diffs, *con_power, *con_square;

	double mindist = 1; 
	for (jj = 1; jj<nn*nn; jj++) {
		mindist = (mindist<(1/connectivity[jj])?mindist:(1/connectivity[jj]));
	}
	double interval = 1e-4*mindist;
	double threshold = 1e-2*mindist; 
	double tiniest = 1e-10*mindist;
	double biggest_move;
	double maxmove = 1;

	change = (double *)R_alloc(dim, sizeof(double)); 
	diffs = (double *)R_alloc(dim, sizeof(double));
	con_power = (double *)R_alloc(nn*nn, sizeof(double));
	con_square = (double *)R_alloc(nn*nn, sizeof(double));


	/*int magic_power = 4;*/
	for (jj = 0; jj < nn*nn; jj++) {
		con_square[jj] = connectivity[jj]*connectivity[jj];
		con_power[jj] = con_square[jj]*connectivity[jj];
	}

	for (iteration = niter; iteration >= 0; iteration--) {
		biggest_move = 0;

		for (jj = 0; jj < nn; jj++) {
			for (kk = 0; kk < dim; kk++) {change[kk]=0;} 
			for (kk = 0; kk < nn; kk++) if (kk != jj) {
				t_dist = 0; 
				for (dd = 0; dd < dim; dd++) {
					diffs[dd] = coords[jj + nn*dd] - coords[kk + nn*dd]; 
					t_dist += diffs[dd]*diffs[dd];
				}
				t_dist = sqrt(t_dist);
	
				for (dd = 0; dd < dim; dd++) diffs[dd] /= t_dist; 
				if (fruchterman) {
					force = t_dist*t_dist*con_power[jj+nn*kk] - 1/t_dist; 
				} else { 
					force = (t_dist - 1/connectivity[jj+nn*kk])*con_square[jj+nn*kk];
				} 
		/*
				force = (t_dist*t_dist*con_power[jj+nn*kk] - 1/t_dist)*fruchterman +
 					((t_dist - 1/connectivity[jj+nn*kk])*con_square[jj+nn*kk])*(1-fruchterman);
 */
				for (dd = 0; dd < dim; dd++) {
					change[dd] -= diffs[dd]*force;
				}
			}

			t_dist = 0; for (dd = 0; dd < dim; dd++) t_dist += change[dd]*change[dd];
			t_dist = sqrt(t_dist);
			total_force = t_dist;


			if (total_force > tiniest) {
				for (dd = 0; dd < dim; dd++) coords[jj+dd*nn] += change[dd]/t_dist*interval;
			
				for (kk = 0; kk < dim; kk++) {change[kk]=0;} 
				for (kk = 0; kk < nn; kk++) if (kk != jj) {
					t_dist = 0; 
					for (dd = 0; dd < dim; dd++) {	
						diffs[dd] = coords[jj + nn*dd] - coords[kk + nn*dd]; 
						t_dist += diffs[dd]*diffs[dd];
					}
					t_dist = sqrt(t_dist);
		
					for (dd = 0; dd < dim; dd++) diffs[dd] /= t_dist; 

				if (fruchterman) {
					force = t_dist*t_dist*con_power[jj+nn*kk] - 1/t_dist; 
				} else { 
					force = (t_dist - 1/connectivity[jj+nn*kk])*con_square[jj+nn*kk];
				}  
			/*
				force = (t_dist*t_dist*con_power[jj+nn*kk] - 1/t_dist)*fruchterman +
 						((t_dist - 1/connectivity[jj+nn*kk])*con_square[jj+nn*kk])*(1-fruchterman);
*/
					for (dd = 0; dd < dim; dd++) {
						change[dd] -= diffs[dd]*force;
					}
				}

				t_dist = 0; for (dd = 0; dd < dim; dd++) t_dist += change[dd]*change[dd];
				t_dist = sqrt(t_dist);
				total_force2 = t_dist;
				for (dd = 0; dd < dim; dd++) change[dd] /= t_dist; 
				t_dist = interval*total_force2/(total_force-total_force2)*(total_force > total_force2);
				if (total_force < total_force2) t_dist = total_force2;

		/*		t_dist = interval*total_force2/abs(total_force-total_force2); */
				biggest_move = (t_dist>biggest_move?t_dist:biggest_move);
				if (t_dist > maxmove) {t_dist = maxmove/t_dist;}
				
				for (dd = 0; dd < dim; dd++) coords[jj+dd*nn] += change[dd]*t_dist;
			} 

		}
		if (biggest_move < threshold) break;
	}

} 



void network_components (const double * sociomatrix,
	const int * nn_c,
	int * components
) {
		
	/* components is preloaded to be 1:n. */
	int nn = *nn_c;
	int moves = 1, tries = 0;
	int jj,kk;
	while (moves) {
		moves = 0;
		for (jj = 0; jj < (nn-1); jj++) for (kk = (jj+1); kk < nn; kk++) {
			if ((sociomatrix[jj+nn*kk]>0 || sociomatrix[kk+nn*jj]>0) && (components[jj] != components[kk])) {
				moves++;
				/* make it the lower one. */
				components[kk]=(components[jj]<components[kk]?components[jj]:components[kk]);
				components[jj] = components[kk];
			}
		}
		tries++;
		if (tries>nn) {
			/* break if it doesn't do it within nn tries. */
			for (jj = 0; jj < (nn-1); jj++) components[jj] = -1;
			moves = 0;
		}
	}
}






