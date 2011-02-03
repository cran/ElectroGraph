/*
####################################################################
# components.cpp 
#
# Andrew C. Thomas <act@acthomas.ca>
# Last modified: January 31, 2011
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Plot functions have been updated for enhanced visualization of weighted 
# edge graphs, egocentric views and potential three-dimensional plotting.
# (as yet unimplemented)
#
####################################################################
*/

#include <iostream>

extern "C" {

void network_components_edges (
			const int * edgelist,
			int * components,
			const int * pnn,
			const int * pedges
	) {
		
	/* components is preloaded to be 1:n. */
	int nn = *pnn; int ee = *pedges;
	int moves = 1, tries = 0;
	int jj,kk;
	//std::cout << "?" << tries << std::endl;
	while (moves) {
		moves = 0;
		for (jj=0; jj<ee; jj++) 
			if (components[edgelist[jj]] != components[edgelist[jj+ee]]) {
				moves++;
				components[edgelist[jj+ee]] = (components[edgelist[jj]]<components[edgelist[jj+ee]]?components[edgelist[jj]]:components[edgelist[jj+ee]]);
				components[edgelist[jj]] = components[edgelist[jj+ee]];
			}
		tries++;
		if (tries>nn) {
			/* break if it doesn't do it within nn tries. */
			for (jj = 0; jj < (nn-1); jj++) components[jj] = -1;
			moves = 0;
		}
		//std::cout << "?" << tries << std::endl;
	}
}


void network_components_socio (const double * sociomatrix,
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





}
