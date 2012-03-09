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
#include <sys/time.h>

const int debug_mode_adapt = 0;

double seconds () {
	struct timeval tv;
	gettimeofday (&tv, NULL);
	double seed = (double)tv.tv_sec + (double)tv.tv_usec *1e-6;
	return seed;
}


void network_projection_by_connectivity (double * connectivity,
	const int * nn_c,
	const int * niter_c,
	const int * dimension_c,
	const int * fruchterman_c,
	const double * egocenters,
	double * coords,
	const int * verbose_c,
	int * trial_count 
) {
	
	int nn = *nn_c, niter = *niter_c, dim = *dimension_c, fruchterman = *fruchterman_c;
	int verbose = *verbose_c;
	double t_dist, force, total_force, total_force2;
	double time_1, time_2;  time_1 = seconds();
	
	//double maxdelta = (double)nn;

	int jj, kk, dd, iteration;
	double *change, *diffs, *con_power, *con_square;

	double mindist = 1; 
	for (jj = 1; jj<nn*nn; jj++) {
		mindist = (mindist<(1/connectivity[jj])?mindist:(1/connectivity[jj]));
	}
	double interval = 1e-4*mindist;
	double threshold = 1e-1*mindist; 
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
		if (niter - iteration == 9) {
			time_2 = seconds(); time_1 = (time_2 - time_1)*(double)niter/10;
			if (verbose) Rprintf ("Maximum seconds remaining in plot: %f\n", time_1);
		}
	}

	*trial_count = iteration;

} 




/* Not currently in use; trial to see if it would be effective failed in graphical presentation.   */

void network_projection_by_connectivity_sparsish (double * connectivity,
	const int * nn_c,
	const int * niter_c,
	const int * dimension_c,
	const int * fruchterman_c,
	const double * egocenters,
	const double * min_closeness,  /* is it closer than this?  */
	double * coords,
	const int * verbose_c,
	int * trial_count 
) {
	
	int nn = *nn_c, niter = *niter_c, dim = *dimension_c, fruchterman = *fruchterman_c;
	int verbose = *verbose_c;
	double t_dist, force, total_force, total_force2;
	double time_1, time_2;  time_1 = seconds();
	
	//double maxdelta = (double)nn;

	int jj, kk, dd, ee, iteration;
	double *change, *diffs, *con_power, *con_square;

	int ** nonzeros = (int **)R_alloc(nn, sizeof(int *));;
	int * holds = (int *)R_alloc(nn, sizeof(int));
	int * degree_effective = (int *)R_alloc(nn, sizeof(int));

	for (jj = 0; jj<nn; jj++) {	
	/*	printf("%i ", jj);   */
		iteration = 0;		
		for (kk = 0; kk<nn; kk++) if (connectivity[nn*kk+jj] >= *min_closeness) {
			holds[iteration] = kk; iteration++;
		}
	/*	printf("%i;", jj);   */

		nonzeros[jj] = (int *)R_alloc(iteration, sizeof(int));
		degree_effective[jj] = iteration;
		for (kk = 0; kk<iteration; kk++) nonzeros[jj][kk] = holds[kk];
	/*	if (iteration == 0) printf ("Trouble: Essentially disconnected node.\n"); else printf ("%i ", iteration);  */
	}


	double mindist = 1; 
	for (jj = 1; jj<nn*nn; jj++) {
		mindist = (mindist<(1/connectivity[jj])?mindist:(1/connectivity[jj]));
	}
	double interval = 1e-4*mindist;
	double threshold = 1e-1*mindist; 
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

		/*	for (kk = 0; kk < nn; kk++) if (kk != jj) {  */
			for (ee = 0; ee < degree_effective[jj]; ee++) {kk = nonzeros[jj][ee];
				t_dist = 0; 
				for (dd = 0; dd < dim; dd++) {
					diffs[dd] = coords[jj + nn*dd] - coords[kk + nn*dd]; 
					t_dist += diffs[dd]*diffs[dd];
				}
				t_dist = sqrt(t_dist); /* distance (  */
	
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
			/*	for (kk = 0; kk < nn; kk++) if (kk != jj) {  */
				for (ee = 0; ee < degree_effective[jj]; ee++) {kk = nonzeros[jj][ee];
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
		if (niter - iteration == 9) {
			time_2 = seconds(); time_1 = (time_2 - time_1)*(double)niter/10;
			if (verbose) Rprintf ("Maximum seconds remaining in plot: %f\n", time_1);
		}
	}

	*trial_count = iteration;

} 



