/* 
#####################################################################
# geodesic-routines.cpp 
#
# Andrew C. Thomas <act@acthomas.ca>
# Last modified: January 19, 2011
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Handles the calculation of geodesically derived quantities,.
#
#
#####################################################################
*/

#include <iostream>
#include <limits>
#include <R_ext/Print.h>

double abso (double arg) {return(arg*(arg>0)-arg*(arg<0));}

extern "C" {

	// Floyd-Warshall algorithm.
	void floyd_warshall (double * output, const int * pnn) { //output is distances.
		int nn = *pnn;
		int change = 1; 
		int jj,kk,ll, iter;
		double val,tt;
		for (iter=0; iter<nn; iter++) {
			change = 0;
			for (kk=0; kk<nn; kk++) for (jj=0; jj<nn; jj++) {
				val = output[kk+nn*jj];
				for (ll=0; ll<nn; ll++) {
					tt = output[kk+nn*ll]+output[ll+nn*jj]; 
					val = (val<tt?val:tt);
				}
				if (val < output[kk+nn*jj]) {
					change++;
					output[kk+nn*jj] = val;
				}
			}
			if (!change) break;
		}
	}

	void dijkstra_single (const int * edges, // input is ee-by-2
				const double * strengths, // ee long.
				double * spdists, // vector of distances from source.
				const int * pnn, const int * pedges, // edgelist
				const int * psrc, const int * pdest // if outside [0, nn-1], won't worry about it.
		) {

		int nn = *pnn; int ee = *pedges;	int kk, choice, dd; int condition = 1; int reps=0;
		double ddchoice;	int * done = new int[nn];	double maxval = std::numeric_limits<double>::max( ); 

		//for (kk = 0; kk<ee; kk++) {if (strengths[kk]>0) maxval += nn/strengths[kk];}
		for (kk = 0; kk<nn; kk++) {spdists[kk] = maxval*(kk != *psrc);	done[kk] = 0;}

		while (condition) {
			choice = 0;	ddchoice = maxval;		//initial seeding.
			for (kk = 0; kk<nn; kk++) { 			// pick the closest that isn't "done".
				if (!done[kk] & ddchoice>spdists[kk]) {choice = kk; ddchoice=spdists[kk];}}
			if (ddchoice < maxval) {		// any more reachable points?
			//	std::cout << "src " << *psrc << " cc " << choice << " " << spdists[choice];
				kk=0; while (edges[kk] != choice && kk < ee) kk++;
				while (kk < ee && edges[kk] == choice) { // edges.
				//	std::cout << "dst " << edges[kk+ee] << " kk " << kk << " strength " << strengths[kk];

					if (strengths[kk]>0 & spdists[edges[kk+ee]]>spdists[choice]+1/strengths[kk]) {
						spdists[edges[kk+ee]] = spdists[choice]+1/strengths[kk];}
					kk++; while (edges[kk] != choice && kk < ee) kk++;}
	
				done[choice] = 1; condition = 0;	for (kk=0; kk<nn; kk++) condition += (done[kk]==0);
				reps++; condition = 1*(reps<nn+1);  // safety.
				condition = (choice==(*pdest)?0:condition);
			//	std::cout << std::endl;
			} else condition = 0;
		}
		
		//catch this or problems happen!
		delete[] done;  //for (kk = 0; kk<nn; kk++) spdists[kk] = (spdists[kk]==maxval?-1:spdists[kk]);
	}

	void dijkstra_all (const int * edges, const double * strengths,
						double * output, // nn by nn
						const int * pnn, const int * pedges) {

		int kk, ll; int nn = *pnn;	double * temp_out = new double[nn];
		
		for (kk=0; kk < *pnn; kk++) {
			for (ll=0; ll< *pnn; ll++) temp_out[ll]=0;
			dijkstra_single (edges, strengths, temp_out, pnn, pedges, &kk, pnn);
			for (ll=0; ll<nn; ll++) output[ll+nn*kk] = temp_out[ll];
		}
		delete[] temp_out;

	}



	// calculates the shortest paths for all nodes beginning with "source", intelligently recursively.
	void betweenness_geodesic_single (
								const int * edges,
								const double * strengths, 

								double * path_count, // matrix, ee by nn.
								double * spdists,

								const int * pnn,
								const int * pedges,

								const int * psrc,
								const int * pdest) {

		int nn = *pnn; int src = *psrc; int ee = *pedges;
		int current_node, kk, ll, mm, keep_going, outpath, current_count;	
		int iters = 0;		double maxval = 0; 
		for (kk = 0; kk<ee; kk++) if (strengths[kk]>0) maxval += 1/strengths[kk];
		for (ll = 0; ll<ee*nn; ll++) path_count[ll] = 0;
		
		int * done = new int[nn]; 	int * node_count = new int[nn]; 		

		dijkstra_single (edges, strengths, spdists, pnn, pedges, psrc, pdest);
	
		//bottom up!
		outpath = 0; iters=0;
		// find the farthest node.
		for (kk=0; kk<nn; kk++) {
			//std::cout << spdists[kk] << " ";
			done[kk]=0; node_count[kk]=0; if (spdists[kk]>spdists[outpath]) outpath = kk;
		}
		done[src] = 1; node_count[src] = 1;

		keep_going = 1; while (keep_going) {
			current_node = outpath; //find the closest untouched node.
			for (ll=0; ll<nn; ll++) {if (spdists[ll] <= spdists[current_node] & !done[ll]) current_node = ll;}
			
			kk=0; while (edges[kk+ee] != current_node && kk < ee) kk++;
			while (kk < ee && edges[kk+ee] == current_node) { // edges.

				if (strengths[kk]>0 & abso(spdists[current_node]-spdists[edges[kk]]-1/strengths[kk])<1/maxval) {
					current_count = 0;//add all the paths that came from the previous steps, times the path count.
					for (ll=0; ll<ee; ll++) {
						path_count[ll+ee*current_node] += path_count[ll+ee*edges[kk]];  //source edge.
					}
					//add the direct path!
					path_count[kk+ee*current_node] += node_count[edges[kk]];
				}

				kk++; while (edges[kk+ee] != current_node && kk < ee) kk++;
			}
			kk=0; while (edges[kk+ee] != current_node && kk < ee) kk++;
			while (kk < ee && edges[kk+ee] == current_node) { // edges.
				if (strengths[kk]>0 & abso(spdists[current_node]-spdists[edges[kk]]-1/strengths[kk])<1/maxval) node_count[current_node] += node_count[edges[kk]];
				kk++; while (edges[kk+ee] != current_node && kk < ee) kk++;
			}

			done[current_node] = 1;	keep_going = 0; for (kk=0; kk<nn; kk++) keep_going += !done[kk];
			iters++; if (iters>2*nn) {REprintf("Something went wrong. Break 2"); break;}
			keep_going = (current_node==(*pdest)?0:keep_going);
		}
		delete[] done; delete[] node_count; //delete[] spdists;
	}





	void betweenness_geodesic_full (
									const int * edges, 
									const double * strengths, 

									double * outweight, // matrix, n by n.
									double * nodeweight,

									const int * pnn,
									const int * pedges,
									const int * distance_weight
		
	) {

		int sources, kk, ll, dest;
		int nn = *pnn; int ee = *pedges; 
		int ddw = *distance_weight;

		double totalweight = 0;
		double * path_count = new double[nn*ee]; // matrix, n^2 by n.
		double * spdists = new double[nn]; 
		double totalpaths;

		//int * psrc = new int[1];
		for (kk=0; kk<ee; kk++) outweight[kk] = 0;
		for (kk=0; kk<nn; kk++) nodeweight[kk] = 0;

		for (sources=0; sources<nn; sources++) {
			
			for (kk=0; kk<nn*ee; kk++) path_count[kk]=0;
			for (kk=0; kk<nn; kk++) spdists[kk]=0;

			betweenness_geodesic_single (edges, strengths, 
				path_count, spdists,	pnn, pedges, &sources, pnn); 

			for (dest=0; dest<nn; dest++) if (dest != sources) {
				totalpaths = 0; for (ll=0; ll<ee; ll++) {totalpaths += path_count[ll+ee*dest]*(edges[ll+ee]==dest);}
				totalpaths += 1*(totalpaths==0);
				for (ll=0; ll<ee; ll++) {
					outweight[ll] += path_count[ll + ee*dest]/totalpaths/(spdists[dest]*ddw+1*(1-ddw));
					nodeweight[edges[ll]] += path_count[ll + ee*dest]/totalpaths/(spdists[dest]*ddw+1*(1-ddw))*(edges[ll] != sources && edges[ll] != dest);
				}
				totalweight += 1/(spdists[dest]*ddw+1*(1-ddw));
			}
		//	std::cout << sources << std::endl;
		}

		for (kk=0; kk<ee; kk++) outweight[kk] /= totalweight;

		delete[] path_count; delete[] spdists;

	}





	void recourse_betweenness_single_sd (const int * edges,	double * strengths, 
											double * path_extra, // matrix, ee by 1.
											double * spdists,

											const int * pnn,	const int * pedges,
											const int * psrc,	const int * pdest,
											const double * ppenalty) {

		//std::cout << *pnn << " " << *pedges << " " << std::endl;

		int nn = *pnn; int src = *psrc; int ee = *pedges; int out_s; int total_paths = 0;
		int kk, ll, mm;	int iters = 0;		double maxval = 0; double stren_hold;

		for (kk = 0; kk<ee; kk++) {if (strengths[kk]>0) maxval += 1/strengths[kk]; path_extra[kk]=0;}
		double * path_count = new double[nn*ee];	for (ll = 0; ll<ee*nn; ll++) path_count[ll] = 0;
		double * hold_dists = new double[nn];

		betweenness_geodesic_single (edges, strengths, path_count, spdists, 
												pnn, pedges, psrc, pdest);
		for (kk=0; kk<ee; kk++) total_paths += path_count[kk+(*pdest)*ee]*(edges[kk+ee]==*pdest);

		for (kk=0; kk<ee; kk++) if (path_count[kk+(*pdest)*ee] != 0) {
			stren_hold = strengths[kk]; strengths[kk] = 1/(1/strengths[kk] + *ppenalty);
			if (spdists[edges[kk]]<spdists[edges[kk+ee]]) {out_s = edges[kk];} else 
				{out_s = edges[kk+ee];}

			dijkstra_single (edges, strengths, hold_dists, 
				pnn, pedges, &out_s, pdest);
			path_extra[kk] = (spdists[out_s] + hold_dists[*pdest] - spdists[*pdest])*path_count[kk+(*pdest)*ee]/total_paths;
			strengths[kk] = stren_hold;			
		}
		delete[] path_count; delete[] hold_dists;
	}


	void recourse_betweenness_one_source (const int * edges,	double * strengths, 
											double * path_extra, // matrix, ee by nn.
											double * spdists,

											const int * pnn,	const int * pedges,
											const int * psrc,	
											const double * ppenalty) {

		//std::cout << *pnn << " " << *pedges << " " << std::endl;

		int nn = *pnn; int src = *psrc; int ee = *pedges; int out_s; int total_paths = 0;
		int pdest = nn;

		int kk, dest, mm;	int iters = 0;		double maxval = 0; double stren_hold;

		for (kk = 0; kk<ee; kk++) {if (strengths[kk]>0) maxval += 1/strengths[kk]; path_extra[kk]=0;}
		double * path_count = new double[nn*ee];	for (kk = 0; kk<ee*nn; kk++) path_count[kk] = 0;
		double * hold_dists = new double[nn];

		betweenness_geodesic_single (edges, strengths, path_count, spdists, 
												pnn, pedges, psrc, &pdest);
	
		for (dest=0; dest<nn; dest++) if (dest != src) { //destinations.
			total_paths=0;	for (kk=0; kk<ee; kk++) total_paths += path_count[kk+dest*ee]*(edges[kk+ee]==dest);
			for (kk=0; kk<ee; kk++) if (path_count[kk+dest*ee] != 0) {
				stren_hold = strengths[kk]; strengths[kk] = 1/(1/strengths[kk] + *ppenalty);
				out_s = (spdists[edges[kk]]<spdists[edges[kk+ee]]?edges[kk]:edges[kk+ee]);

				dijkstra_single (edges, strengths, hold_dists, 
					pnn, pedges, &out_s, &dest);
				path_extra[kk+dest*ee] = (spdists[out_s] + hold_dists[dest] - spdists[dest])*path_count[kk+(dest)*ee]/total_paths;
				strengths[kk] = stren_hold;			
			}
		}

		delete[] path_count; delete[] hold_dists;
	}


	void recourse_betweenness_full (
				const int * edges, 	double * strengths, 
				double * outweight, // ee vector.

				const int * pnn,	const int * pedges,
				const int * distance_weight,	const double * ppenalty) {

		int sources, kk, ll, dest;
		int nn = *pnn; int ee = *pedges; 
		int ddw = *distance_weight;

	//	std::cout << "cbf: " << nn << " " << ee << std::endl;

		double totalweight = 0;
		double * path_extra = new double[nn*ee]; // matrix, n^2 by n.
		double * spdists = new double[nn]; 
		double totalpaths;

		//int * psrc = new int[1];
		for (kk=0; kk<ee; kk++) outweight[kk] = 0;

		for (sources=0; sources<nn; sources++) {
		//	std::cout << "source: " << sources << std::endl;
			for (kk=0; kk<nn*ee; kk++) path_extra[kk]=0;
			for (kk=0; kk<nn; kk++) spdists[kk]=0;
			recourse_betweenness_one_source (edges, strengths, 
								path_extra, spdists,	pnn, 
								pedges,	&sources,	ppenalty); 

			for (dest=0; dest<nn; dest++) if (dest != sources) {
				for (ll=0; ll<ee; ll++) {outweight[ll] += path_extra[ll + ee*dest]/(spdists[dest]*ddw+1*(1-ddw));}
			//	std::cout << path_extra[ll + ee*dest] << " " << 1/(spdists[dest]*ddw+1*(1-ddw)) << " " << std::endl;
				totalweight += 1/(spdists[dest]*ddw+1*(1-ddw));
			}
		}
		for (kk=0; kk<ee; kk++) outweight[kk] /= totalweight;
		delete[] path_extra; delete[] spdists;
	}



	void clustering_statistics_socio (
		const double * sociomatrix_c,
		const int * nn_row_c,
		double * transitives,
		double * cycles) {

		int nn_row = *nn_row_c;
		int ii,jj,kk;

		double sum_n, sum_d;		

		//for (kk=0; kk<nn_row; kk++) {transitives[kk]=0;	cycles[kk]=0;}
		for (kk=0; kk<nn_row; kk++) {
			sum_n=0; sum_d=0; //cycles.
			for (ii=0; ii<nn_row; ii++) for (jj=0; jj<nn_row; jj++) {
				sum_n += (sociomatrix_c[ii+nn_row*jj]*sociomatrix_c[jj+nn_row*kk]*
								sociomatrix_c[kk+nn_row*ii]);
				sum_d += (sociomatrix_c[jj+nn_row*kk]*sociomatrix_c[kk+nn_row*ii]);
			}
			if (sum_d==0) sum_d = 1;
			cycles[kk] = sum_n/sum_d;

			sum_n=0; sum_d=0; //transitives.
			for (ii=0; ii<nn_row; ii++) for (jj=0; jj<nn_row; jj++) {
				sum_n += (sociomatrix_c[ii+nn_row*jj]*sociomatrix_c[jj+nn_row*kk]*
								sociomatrix_c[ii+nn_row*kk]);
				sum_d += (sociomatrix_c[jj+nn_row*kk]*sociomatrix_c[ii+nn_row*kk]);
			}
			if (sum_d==0) sum_d = 1;
			transitives[kk] = sum_n/sum_d;
		}

	}


	void short_length_statistics (const int * edges, const double * values, //edges 0:n-1.
								const int * pnn, const int *pedge,
								double * indegree, double * outdegree,
								double * transitives, double * cycles) {


		int nn = *pnn; int ee = *pedge;
		double * sociomatrix = new double[nn*nn];
		int kk,ll;

		//std::cout << "sls " << nn << " " << ee << std::endl;

		for (kk=0; kk<nn; kk++) {indegree[kk]=0; outdegree[kk]=0; transitives[kk]=0; cycles[kk]=0;}
		for (kk=0; kk<nn*nn; kk++) sociomatrix[kk]=0;
		//for (kk=0; kk<nn; kk++) std::cout << outdegree[kk] << " " << indegree[kk] << " " << transitives << " " << cycles << std::endl;

		for (ll=0; ll<ee; ll++) {
			sociomatrix[edges[ll]+nn*edges[ll+ee]] = values[ll];
			outdegree[edges[ll]] += values[ll];
			indegree[edges[ll+ee]] += values[ll];
		}
		//std::cout << "sls2 " << nn << " " << ee << " "; for (kk=0; kk<nn; kk++) std::cout << outdegree[kk]; std::cout << std::endl;

		clustering_statistics_socio (sociomatrix, pnn, transitives, cycles);

		//std::cout << "sls3 " << nn << " " << ee << " ";
		//for (kk=0; kk<nn; kk++) std::cout << outdegree[kk] << " " << indegree[kk] << " " << transitives << " " << cycles << std::endl;

		delete[] sociomatrix;
	}



/*   /////////////////////////////////////////////////////////////////////////////
Legacy routines here.
/////////////////////////////////////////////////////////////////////////////   */


	void dijkstra_single_socio (const double * input, // input matrix of closenesses.
								double * output, // vector of distances from source.
								const int * pnn, 
								const int * sourcey) {
		// let the edges represent inverse distances.
		// assumes complete connectivity.

		int nn = *pnn;
		int kk, choice; int condition = 1; int reps=0;
		
		double ddchoice;
		int * done = new int[nn];		

		//initial loading.
		double maxval = 0; for (kk = 0; kk<nn*nn; kk++) {if (input[kk]>0) maxval += 1/input[kk];}
		for (kk = 0; kk<nn; kk++) {
			output[kk] = maxval*(kk != *sourcey);
			done[kk] = 0;
		}
 
	//	std::cout << maxval << " " << *sourcey << std::endl;

		while (condition) {
			//initial seeding.
			choice = 0;	ddchoice = maxval;
			for (kk = 0; kk<nn; kk++) { // pick the closest that isn't "done".
		//		std::cout << done[kk] << " " << ddchoice << " " << output[kk] << " ;";

				if (!done[kk] & ddchoice>output[kk]) {choice = kk; ddchoice=output[kk];}
			}
		//	std::cout << output[choice] << std::endl;
			for (kk = 0; kk<nn; kk++) { // check all edges from choice to kk.
				if (input[choice+nn*kk]>0 & output[kk]>output[choice]+1/input[choice+nn*kk]) {output[kk] = output[choice]+1/input[choice+nn*kk];}
			//	std::cout << done[kk] << " " << ddchoice << " " << output[kk] << " ;";
			}
		//	std::cout << output[choice] << std::endl;
			
			done[choice] = 1; condition = 0;
			for (kk=0; kk<nn; kk++) condition += (done[kk]==0);

		//	std::cout << choice << " " << output[choice] << std::endl;

			reps++; condition = 1*(reps<nn+1);
		}

		//std::cout << "???"; for (kk = 0; kk<nn; kk++) {std::cout << output[kk] << " ";} std::cout << std::endl;

		delete[] done;
	}



	// calculates the shortest paths for all nodes beginning with "source", intelligently recursively.
	void betweenness_geodesic_single_socio (const double * sociomatrix, 
										double * path_count, // matrix, n^2 by n.
										double * spdists,
										const int * pnn,
										const int * psrc) {

		int nn = *pnn; int src = *psrc;
		int current_node, kk, ll, mm, keep_going, outpath, current_count;	
		int iters = 0;
		double maxval = 0; 
		for (kk = 0; kk<nn*nn; kk++) {
			if (sociomatrix[kk]>0) maxval += 1/sociomatrix[kk];
			for (ll=0; ll<nn; ll++) path_count[kk+nn*nn*ll] = 0;
		}

		int * done = new int[nn]; 	int * node_count = new int[nn]; 		

		//std::cout << " " << nn << " " << spdists[nn-1] << " " << std::endl;
		dijkstra_single_socio (sociomatrix, spdists, pnn, psrc);

		//std::cout << " " << nn << " " << spdists[nn-1] << " " << std::endl;
		//for (kk=0; kk<nn; kk++) {std::cout << spdists[kk] << " ";} std::cout << std::endl;

		//bottom up!
		outpath = 0; iters=0;
		// find the farthest node.
		for (kk=0;kk<nn;kk++) {done[kk]=0; node_count[kk]=0; if (spdists[kk]>spdists[outpath]) outpath = kk;}
		done[src] = 1; node_count[src] = 1;
		keep_going = 1; while(keep_going) {
			current_node = outpath; 
			for (kk=0; kk<nn; kk++) {if (spdists[kk]<=spdists[current_node] & !done[kk]) current_node = kk;}
	
			//std::cout << " " << current_node << " " << spdists[current_node] << " " << std::endl;
			// what are the paths for all predecessors?			
			for (kk=0; kk<nn; kk++) {

				if (sociomatrix[kk+nn*current_node]>0 & abso(spdists[current_node]-spdists[kk]-1/sociomatrix[kk+nn*current_node])<1/maxval) {
					// then, there is a path here.
				//	std::cout << current_node << " " << kk << " ";

					current_count = 0;
					//add all the paths that came from the previous steps, times the path count.
					for (ll=0; ll<nn; ll++) for (mm=0; mm<nn; mm++) {
						path_count[mm+nn*ll+nn*nn*current_node] += path_count[mm+nn*ll+nn*nn*kk];
				//		std::cout << mm << " " << kk << ", ";
					}
				//	std::cout << std::endl;
					//add the direct path!
					path_count[kk+nn*current_node+nn*nn*current_node] += node_count[kk];
				}
			}
			for (kk=0; kk<nn; kk++) if (sociomatrix[kk+nn*current_node]>0 & abso(spdists[current_node]-spdists[kk]-1/sociomatrix[kk+nn*current_node])<1/maxval) {node_count[current_node] += node_count[kk];}


			done[current_node] = 1;
			keep_going = 0; for (kk=0; kk<nn; kk++) keep_going += !done[kk];
		//	for (kk=0; kk<nn; kk++) std::cout << done[kk];			std::cout << std::endl;

			iters++; if (iters>2*nn) {REprintf("Something went wrong. Break 2"); break;}
		}

//		for (kk=0; kk<nn; kk++) std::cout << kk << "," << node_count[kk] << "; ";	std::cout << std::endl;

		delete[] done; delete[] node_count; //delete[] spdists;

//		for (kk=0; kk<nn; kk++) std::cout << kk << ",; ";		std::cout << std::endl;

	}


	void betweenness_geodesic_full_socio (const double * sociomatrix, 
									double * outweight, // matrix, n by n.
									double * nodeweight,
									const int * pnn,
									const int * distance_weight) {

		int sources, kk, ll, dest;
		int nn = *pnn; int nn2 = nn*nn; int nn3 = nn2*nn;
		int ddw = *distance_weight;

		double totalweight = 0;
		double * path_count = new double[nn3]; // matrix, n^2 by n.
		double * spdists = new double[nn]; 
		double totalpaths;

		//int * psrc = new int[1];
		for (kk=0; kk<nn2; kk++) outweight[kk] = 0;
		for (kk=0; kk<nn; kk++) nodeweight[kk] = 0;

		for (sources=0; sources<nn; sources++) {
			
			for (kk=0; kk<nn3; kk++) path_count[kk]=0;
			for (kk=0; kk<nn; kk++) spdists[kk]=0;

			betweenness_geodesic_single_socio (sociomatrix, path_count, // matrix, n^2 by n.
											spdists,	pnn, &sources); 
			
			for (dest=0; dest<nn; dest++) if (dest != sources) {
				totalpaths = 0; for (ll=0; ll<nn; ll++) {totalpaths += path_count[ll + (nn+1)*nn*dest];} 
				totalpaths += 1*(totalpaths==0);
				for (ll=0; ll<nn2; ll++) {
					outweight[ll] += path_count[ll + nn2*dest]/totalpaths/(spdists[dest]*ddw+1*(1-ddw));
					nodeweight[ll % nn] += path_count[ll + nn2*dest]/totalpaths/(spdists[dest]*ddw+1*(1-ddw))*(ll % nn != sources && ll % nn != dest);
				}
				totalweight += (spdists[dest]*ddw+1*(1-ddw));
			}
		//	std::cout << sources << std::endl;
		}

		for (kk=0; kk<nn2; kk++) outweight[kk] /= totalweight;

		delete[] path_count; delete[] spdists;

	}


	void sym_test (const int * edges,
			double * strengths,
			const int * pnn, const int * pedge,
			int * sym) {

		int * knockout = new int[*pedge];
		int current, kk; for (kk=0; kk<(*pedge); kk++) knockout[kk]=1;
		int keep_going = 1;
		while (keep_going) {
			current=0; while (!knockout[current] && current < *pedge) current++;
		//	std::cout << *pnn << ":" << current << ":";
			if (current < *pedge) {
				kk = current + 1; while ((edges[current] != edges[kk+*pedge] || edges[kk] != edges[current+*pedge]) && kk < *pedge) kk++;
				
		//		std::cout << kk << std::endl;

				if (kk < *pedge && strengths[current] == strengths[kk]) {
					knockout[kk] = 0; knockout[current] = 0;	
				} else {
					keep_going = 0;	
				}
			} else keep_going = 0;
		}
		*sym = 0; for (kk=0; kk<(*pedge); kk++) *sym += knockout[kk];
	}

























}
