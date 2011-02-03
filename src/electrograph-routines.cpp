/* 
#####################################################################
# electrograph-routines.cpp 
#
# Andrew C. Thomas <act@acthomas.ca>
# Last modified: January 31, 2011
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Plot functions have been updated for enhanced visualization of weighted 
# edge graphs, egocentric views and potential three-dimensional plotting.
# (as yet unimplemented)
#
# Modification: Placed the option of saving voltages/currents.
#
#####################################################################
*/

#include "scythe-extended.h"
#include "electrograph-routines.h"

int debug_mode = 0;

inline double abs (double input) {return (input<0?-input:input);}

Matrix<int> true_order(Matrix<double> input) {
	Matrix<int> out (input.rows()*input.cols(),1,1,0);
	Matrix<int> hold = order(input);
	for (int ii=0; ii<(input.rows()*input.cols()); ii++) {out(hold(ii)) = ii;}
	return out;	
}

Matrix<int> rev_true_order(Matrix<double> input) {
	int l_out = input.rows()*input.cols();
	Matrix<int> out (l_out,1,1,0);
	Matrix<int> hold = order(input);
	for (int ii=0; ii<(l_out); ii++) {out(l_out-1-hold(ii)) = ii;}
	return out;	
}

inline double dst (Matrix<double> c1, Matrix<double> c2) {
	double out = 0;
	for (int dd=0; dd<c1.cols(); dd++) out += pow(c1(dd)-c2(dd),2);
	out = sqrt(out);
	return out;
} 




//given standard current, flip where necessary.
void process_for_enemies (
	Matrix<double> sociomatrix,
	Matrix<double> fidelities,
	Matrix<double> voltages,	
	Matrix<double> & average_current,	
	Matrix<double> & red_current,
	double & total_conductance_strength, 	
	double & total_conductance_fidelity 
) {

	int nn_row = voltages.rows();
	Matrix<int> t_ord = rev_true_order(voltages);
	int kk, cc, curr_pick;

	double total_current = sum(average_current(_,t_ord(nn_row-1)));
	double av_hl, sd_hl, pos_frac, t_curr;
	
	for (kk=0; kk<(nn_row-1); kk++) {
		curr_pick = t_ord[kk];
		av_hl = 0; sd_hl = 0; //black current, red current.
		if (kk>0) for (cc=0; cc<nn_row; cc++) {
			av_hl += average_current(cc,curr_pick);
			sd_hl += red_current(cc,curr_pick);
		} else {
			av_hl = total_current;
		}
		if (av_hl+sd_hl>0) {
			pos_frac = av_hl/(av_hl+sd_hl);
			for (cc=0; cc<nn_row; cc++) {
				t_curr = average_current(curr_pick,cc)+red_current(curr_pick,cc); //outbound total.

				average_current(curr_pick,cc) = pos_frac*t_curr*(1+fidelities(curr_pick,cc))/2 + 
					(1-pos_frac)*t_curr*(1-fidelities(curr_pick,cc))/2;
				red_current(curr_pick,cc) = pos_frac*t_curr*(1-fidelities(curr_pick,cc))/2 + 
					(1-pos_frac)*t_curr*(1+fidelities(curr_pick,cc))/2;

/*				average_current(curr_pick,cc) = pos_frac*t_curr*(sociomatrix(curr_pick,cc)>0)+
					(1-pos_frac)*t_curr*(sociomatrix(curr_pick,cc)<0);
				red_current(curr_pick,cc) = pos_frac*t_curr*(sociomatrix(curr_pick,cc)<0)+
					(1-pos_frac)*t_curr*(sociomatrix(curr_pick,cc)>0);  */
				

			}
		}
	}

	curr_pick = t_ord[nn_row-1];
	av_hl = 0; sd_hl = 0;
	for (cc=0; cc<nn_row; cc++) {
		av_hl += average_current(cc,curr_pick);
		sd_hl += red_current(cc,curr_pick);
	}	
	total_conductance_strength = av_hl+sd_hl;
	total_conductance_fidelity = (av_hl-sd_hl)/(av_hl+sd_hl); //1 to -1.

}

extern "C" {
	void process_for_enemies_c (
		const double * sociomatrix_c,
		const double * fidelities_c,
		const double * voltages_c,
		double * average_current_c,
		double * red_current_c,
		double * total_conductance_strength_c,
		double * total_conductance_fidelity_c,

		const int * nn_c
	) {
		int nn = *nn_c;
		int jj,kk;
		Matrix<double> sociomatrix(nn,nn,sociomatrix_c);
		Matrix<double> fidelities(nn,nn,fidelities_c);
		Matrix<double> voltages(nn,1,voltages_c);
		Matrix<double> average_current(nn,nn,1,0);
		Matrix<double> red_current(nn,nn,1,0);
		double total_conductance_strength = * total_conductance_strength_c;
		double total_conductance_fidelity = * total_conductance_fidelity_c;

		for (jj=0; jj<nn; jj++) for (kk=0; kk<nn; kk++) {
			average_current(jj,kk) = (voltages(jj)-voltages(kk))*(voltages(jj)>voltages(kk))*abs(sociomatrix(jj,kk));
		}

		process_for_enemies(sociomatrix,	fidelities,	voltages, 
									average_current,	red_current, 
									total_conductance_strength, total_conductance_fidelity);

		*total_conductance_strength_c = total_conductance_strength;
		*total_conductance_fidelity_c = total_conductance_fidelity;
		for (jj=0; jj<nn*nn; jj++) {
			average_current_c[jj] = average_current(jj);
			red_current_c[jj] = red_current(jj);
		}

	}

}



void process_without_enemies (
	Matrix<double> & average_current,	int snk, 
	double & total_conductance_strength, 	double & total_conductance_fidelity 
) {
	double av_hl = 0; int cc;
	for (cc=0; cc<average_current.rows(); cc++) {
		av_hl += average_current(cc,snk);
	}	
	if (debug_mode > 2) std::cout << "wo " << snk << ":" << av_hl << std::endl;
	total_conductance_strength = av_hl;
	total_conductance_fidelity = 1; //1 to -1.

}



void solve_votes_symmetric (
	Matrix<double> sociomatrix,		Matrix<double> fidelities,		
	int src,	int snk,
	double & equiv_resist,	double & fidelity,
	Matrix<double> & voltages,	Matrix<double> & currents, //both clean.
	Matrix<double> & holder, //(nn_row-2,nn_row-2,1,0);
	Matrix<double> & average_current,
	Matrix<double> & red_current
) {
	
	//if (debug_mode) std::cout << "in svs" << std::endl;
	// sociomatrix is checked symmetric.
	int nn_row = sociomatrix.rows();
	Matrix<double> bb((nn_row - 2), 1, 1, 0);
	int jj,kk,rr,cc,enemies;
	double running_sum, total_conductance=0;

	//	Preprocessing takes care of "enemies".
	cc = -1;
	for (jj=0; jj<nn_row; jj++) {
		if (jj != src && jj != snk) {
			cc++;
			rr = -1;
			running_sum = 0;
			for (kk=0; kk<nn_row; kk++) {
				running_sum += sociomatrix(jj,kk);
				if (kk != src && kk != snk) {
					rr++;
					holder(rr,cc) = sociomatrix(jj,kk);
				}				
			}
			holder(cc,cc) = (running_sum>0?-running_sum:-1); //should fix isolates.
		}
		if (jj == src) {
			rr = -1;
			for (kk=0; kk<nn_row; kk++) if (kk != src && kk != snk) {
				rr++; bb(rr) = -sociomatrix(jj,kk);
			}
		}
	}
	bb = lu_solve(holder,bb);

	rr = -1;
	for (kk=0; kk<nn_row; kk++) if (kk == src) {
		voltages(kk) = 1;} else if (kk == snk) {voltages(kk) = 0;
	} else {rr++; voltages(kk) = bb(rr);}


	// calculate currents through nodes.
	enemies = 0;
	for (jj=0; jj<nn_row; jj++) {
		running_sum = 0;
		for (kk = 0; kk<nn_row; kk++) {
			average_current(jj,kk) = (voltages(jj)-voltages(kk))*(voltages(jj)>voltages(kk))*sociomatrix(jj,kk);
			red_current(jj,kk) = 0;

			running_sum += average_current(jj,kk);
			enemies += 1*(fidelities(jj,kk)<1); 
		}
		currents(jj) = running_sum;
	}

	if (enemies) {
		process_for_enemies (sociomatrix, fidelities, 
						voltages, average_current,
						red_current, total_conductance, fidelity);
	} else {
		process_without_enemies (average_current,	snk, 
			total_conductance, fidelity);
	}
	equiv_resist = 1/total_conductance;

}

extern "C" {	
	void solve_volts_symmetric_c (
			const double * sociomatrix_c, 
			const double * fidelities_c,
			const int * nn_row_c,
			const int * src_c, const int * snk_c,
			double * equiv_resist_c, double * fidelity_c,
			double * voltages_c, double * currents_c) {

		int nn_row = *nn_row_c; int src = *src_c; int snk = *snk_c;
		int ii;
		Matrix<double> sociomatrix(nn_row,nn_row,sociomatrix_c);
		Matrix<double> fidelities(nn_row,nn_row,fidelities_c);

		double equiv_resist = *equiv_resist_c; double fidelity = *fidelity_c;
		Matrix<double> voltages(nn_row,1,voltages_c); 	Matrix<double> currents(nn_row,1,currents_c); 		
		Matrix<double> holder (nn_row-2,nn_row-2,1,0);
		Matrix<double> average_current_hold(nn_row,nn_row,1,0);
		Matrix<double> red_current_hold(nn_row,nn_row,1,0);
	
		//sociomatrix = abs(sociomatrix);

		solve_votes_symmetric (sociomatrix,	fidelities,	
				src, snk, equiv_resist,	fidelity, voltages, currents, 
				holder, average_current_hold, red_current_hold);	

		*equiv_resist_c = equiv_resist; *fidelity_c = fidelity;	
		for (ii=0; ii < nn_row; ii++) {
			voltages_c[ii] = voltages(ii);
			currents_c[ii] = currents(ii);	
		}

	}
}



void solve_votes_notsymmetric_raw (
		Matrix<double> sociomatrix,	Matrix<double> fidelities,	
		int src,	int snk,
		Matrix<int> asymmetric_entries, // k by 2.
		double & equiv_resist,	double & fidelity,
		Matrix<double> & voltages,	Matrix<double> & currents,
		Matrix<double> & holder, //(nn_row-2,nn_row-2,1,0);
		Matrix<double> & average_current, Matrix<double> & red_current 
) {

	int nn_row = sociomatrix.rows();
	//Matrix<double> current_flows(asymmetric_entries.rows(),1);
	Matrix<int> ascending(asymmetric_entries.rows(),1,1,1);
	int asy_count,kk; //,jj,kk;
	double current,flipper, flipper_fid;
	int checker_counter = 0;
	double min_nonzero = 0;
	for (kk=0; kk<nn_row*nn_row; kk++) 
		min_nonzero = (sociomatrix(kk)>min_nonzero?sociomatrix(kk):min_nonzero);
	for (kk=0; kk<nn_row*nn_row; kk++) 
		min_nonzero = (sociomatrix(kk)<min_nonzero && sociomatrix(kk)>0?sociomatrix(kk):min_nonzero);

	min_nonzero = min_nonzero*1e-9;

	//if (debug_mode) std::cout << "in sv-non" << std::endl;

	int done = 0;
	
	for (asy_count=0; asy_count<ascending.rows(); asy_count++) {
		if (sociomatrix(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1))==0) sociomatrix(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1)) = min_nonzero;
		if (sociomatrix(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0))==0) sociomatrix(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0)) = min_nonzero;

	//	if (sociomatrix(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1))==min_nonzero &&
	//		sociomatrix(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0))>0) ascending(asy_count) = 0; 
	}
	Matrix<double> socio_hold = sociomatrix;
	Matrix<double> fidel_hold = fidelities;
	

	while (!done) {

		checker_counter += 1;
		done = 1;
		for (asy_count = 0; asy_count < asymmetric_entries.rows(); asy_count++) {
			flipper = ascending(asy_count)*sociomatrix(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1))+
				(1-ascending(asy_count))*sociomatrix(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0));
			socio_hold(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1)) = flipper;
			socio_hold(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0)) = flipper;

			flipper_fid = ascending(asy_count)*fidelities(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1))+
				(1-ascending(asy_count))*fidelities(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0));
			fidel_hold(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1)) = flipper_fid;
			fidel_hold(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0)) = flipper_fid;


		}

		//	if (debug_mode) std::cout << symmetric(socio_hold) << std::endl << t(voltages) << std::endl;	

		solve_votes_symmetric(socio_hold, fidel_hold, src, snk, 
			equiv_resist, fidelity, voltages,	currents, 
			holder, average_current, red_current);

		//	if (debug_mode) std::cout << t(voltages) << std::endl;	

		for (asy_count = 0; asy_count < asymmetric_entries.rows(); asy_count++) {
			current = voltages(asymmetric_entries(asy_count,0))-voltages(asymmetric_entries(asy_count,1));
			if (ascending(asy_count)) {
				current *= sociomatrix(asymmetric_entries(asy_count,0),asymmetric_entries(asy_count,1));
			} else {
				current *= sociomatrix(asymmetric_entries(asy_count,1),asymmetric_entries(asy_count,0));
			}

			if ((ascending(asy_count) && current<0) || (!ascending(asy_count) && current>0)) {
				done = 0; ascending(asy_count) = 1 - ascending(asy_count);
			}
		}

		if (checker_counter > nn_row*sqrt(nn_row)) {
			done = 1;
		//	std::cout << "Convergence trouble: Solve-votes-asymmetric is not converging for node pair " << src+1 << ", " << snk+1 << "." << std::endl;
			for (asy_count=0; asy_count < voltages.rows(); asy_count++) {
				voltages(asy_count) = -1;
			}
		}
	}
}

extern "C" {
	void solve_volts_notsymmetric_c (
			const double * sociomatrix_c, 
			const double * fidelities_c, const int * nn_row_c,
			const int * src_c, const int * snk_c,
			double * equiv_resist_c, double * fidelity_c,
			double * voltages_c, double * currents_c) {

		//	if (debug_mode) std::cout << "svn_c" << std::endl;

		int nn_row = *nn_row_c;
		Matrix<double> sociomatrix(nn_row,nn_row,sociomatrix_c);
		Matrix<double> fidelities(nn_row,nn_row,fidelities_c);
		int src = *src_c; int snk = *snk_c;
		double equiv_resist = 420; double fidelity = *fidelity_c;
		Matrix<double> voltages(nn_row,1,voltages_c); 	Matrix<double> currents(nn_row,1,currents_c); 		
		int ii,rr,cc,ss;
		Matrix<double> holder (nn_row-2,nn_row-2,1,0);
		Matrix<double> average_current_hold(nn_row,nn_row,1,0);
		Matrix<double> red_current_hold(nn_row,nn_row,1,0);

		//	if (debug_mode) std::cout << t(voltages) << std::endl;

		int total_asymm = 0;
		for (rr=0; rr<(nn_row-1); rr++) for (cc=(rr+1); cc<nn_row; cc++) total_asymm += 1*(sociomatrix(rr,cc)!=sociomatrix(cc,rr));
		int t_symm = 1*(total_asymm == 0) + total_asymm;
		Matrix<int> asymmetric(t_symm,2,1,0);
		ss = 0;
		for (rr=0; rr<(nn_row-1); rr++) for (cc=(rr+1); cc<nn_row; cc++) if (sociomatrix(rr,cc)!=sociomatrix(cc,rr)) {
			asymmetric(ss,0) = rr; asymmetric(ss,1) = cc; ss++;
		}
	
		//	if (debug_mode) std::cout << asymmetric << std::endl;

		if (total_asymm) {
			solve_votes_notsymmetric_raw (sociomatrix, fidelities, src, snk, 
				asymmetric,	equiv_resist, fidelity, voltages, currents, 
				holder, average_current_hold, red_current_hold);		
		} else {
			solve_votes_symmetric (sociomatrix, fidelities, src, snk, 
				equiv_resist, fidelity, voltages, currents, holder, 
				average_current_hold, red_current_hold);
		}

		//	if (debug_mode) std::cout << "svn_c 3" << std::endl;

		*equiv_resist_c = equiv_resist; *fidelity_c = fidelity;	
		for (ii=0; ii < nn_row; ii++) {
			voltages_c[ii] = voltages(ii);
			currents_c[ii] = currents(ii);	
		}
	//	for (ii=0; ii < nn_row*nn_row; ii++) {
	//		average_current_c[ii] = average_current(ii);
	//	}

	}
}



int symmetric (Matrix<double> sym_obj) {
	int total_asymm; int nn_row = sym_obj.rows();
	for (int rr=0; rr<(nn_row-1); rr++) for (int cc=(rr+1); cc<nn_row; cc++) total_asymm += 1*(sym_obj(rr,cc)!=sym_obj(cc,rr));
	return 1*(total_asymm==0);
}


void get_resistances_fast_symmetric (
		Matrix<double> sociomatrix,	Matrix<double> fidelities,  
		Matrix<int> sources, Matrix<int> sinks,
		Matrix<double> & equiv_resist, Matrix<double> & fidelity,
		Matrix<double> & blk_curr_a, Matrix<double> & red_curr_a, //) {
		Matrix<double> & blk_curr_v, Matrix<double> & red_curr_v,
		Matrix<double> & blk_curr_p, Matrix<double> & red_curr_p,
		Matrix<double> & ll_save_1, Matrix<double> & ll_save_2, Matrix<double> & ll_save_3	
	) {
	
	int nn_row = sociomatrix.rows();
	int dyads = sources.rows();
	int dyad,jj,kk,rr,cc, t_source, t_sink;
	//output quantities are already defined.

	int enemies_present = 0;
	for (jj=0; jj<nn_row*nn_row; jj++) enemies_present += 1*(fidelities(jj)<1);

	Matrix<double> blk_curr_hold(nn_row,nn_row,1,0);
	Matrix<double> red_curr_hold(nn_row,nn_row,1,0);
	Matrix<double> voltages(nn_row,1,1,0);
	double resist_hold, fidelity_hold;

	//transfer contents of sociomatrix to laplacian.
	double rowsum;	
	for (jj=0; jj<nn_row; jj++) {
		rowsum=0;
		for (kk=0; kk<nn_row; kk++) {
			rowsum += sociomatrix(jj,kk);
			//ll_save		
			if (jj != 0 && kk != 0) {ll_save_1(jj-1,kk-1) = -sociomatrix(jj,kk);}
			if (jj != 1 && kk != 1) {
				rr = (jj<1?jj:jj-1); cc = (kk<1?kk:kk-1); ll_save_2(rr,cc) = -sociomatrix(jj,kk);
			}
			if (jj != 2 && kk != 2) {
				rr = (jj<2?jj:jj-1); cc = (kk<2?kk:kk-1); ll_save_3(rr,cc) = -sociomatrix(jj,kk);
			}
		}
		if (jj!=0) {ll_save_1(jj-1,jj-1) = rowsum;}
		if (jj!=1) {if (jj<1) ll_save_2(jj,jj) = rowsum; else ll_save_2(jj-1,jj-1) = rowsum;}
		if (jj!=2) {if (jj<2) ll_save_3(jj,jj) = rowsum; else ll_save_3(jj-1,jj-1) = rowsum;}
	}

//	if (debug_mode > 1) std::cout << ll_save_1 << ll_save_2 << ll_save_3 << std::endl;


	//inversions.
	ll_save_1 = inv(ll_save_1);
	ll_save_2 = inv(ll_save_2);
	ll_save_3 = inv(ll_save_3);

	if (debug_mode > 0) std::cout << "Inversions complete." << std::endl;

	double max_v, min_v;
	for (dyad=0; dyad < dyads; dyad++) {
		for (kk=0; kk < nn_row; kk++) voltages(kk)=0;
		for (kk=0; kk < nn_row*nn_row; kk++) blk_curr_hold(kk)=0;

		if (sources(dyad) != 0 && sinks(dyad) != 0) { //1+, 2+.
			for (kk=0; kk<nn_row-1; kk++) voltages(kk+1) = ll_save_1(kk,sources(dyad)-1) - ll_save_1(kk,sinks(dyad)-1);
		} else {
			if (sources(dyad) != 1 && sinks(dyad) != 1) { //0 and 2+.
				t_source = (sources(dyad)>0?sources(dyad)-1:0); 
				t_sink = (sinks(dyad)>0?sinks(dyad)-1:0);
				voltages(0) = ll_save_2(0,t_source) - ll_save_2(0,t_sink);
				for (kk=1; kk<nn_row-1; kk++) {
					voltages(kk+1) = ll_save_2(kk,t_source) - ll_save_2(kk,t_sink);
				}
			} else { //0,1. Omit 2
				t_sink = sinks(dyad); t_source = sources(dyad);
				for (kk=0; kk<2; kk++) voltages(kk) = ll_save_3(kk,t_source) - ll_save_3(kk,t_sink);
				for (kk=2; kk<nn_row-1; kk++) {
					voltages(kk+1) = ll_save_3(kk,t_source) - ll_save_3(kk,t_sink);
				}
			}
		}

		max_v=0; min_v=0; for (kk=0; kk<nn_row; kk++) {
			max_v = (max_v<voltages(kk)?voltages(kk):max_v);
			min_v = (min_v>voltages(kk)?voltages(kk):min_v);
		}
		//equiv_resist(jj)=max_v - min_v; 
		resist_hold = max_v-min_v; fidelity_hold=1;
		for (jj=0; jj<nn_row; jj++) for (kk=0; kk<nn_row; kk++) {
			blk_curr_hold(jj,kk) = (voltages(jj)-voltages(kk))*
				(voltages(jj)>voltages(kk))*abs(sociomatrix(jj,kk));
			red_curr_hold(jj,kk) = 0;
		}
		
	//	if (debug_mode > 2) std::cout << "Black current:" << std::endl << blk_curr_hold << std::endl;

		if (enemies_present) {
			process_for_enemies (sociomatrix, fidelities, 
				voltages, blk_curr_hold, 
				red_curr_hold, resist_hold, fidelity_hold);
		} 

		equiv_resist(dyad) = resist_hold; fidelity(dyad) = fidelity_hold;

		blk_curr_a += blk_curr_hold; red_curr_a += red_curr_hold;

		blk_curr_hold /= sqrt(resist_hold); red_curr_hold /= sqrt(resist_hold);
		blk_curr_p += blk_curr_hold; red_curr_p += red_curr_hold;
	//	if (debug_mode > 2) std::cout << "Black current 2 :" << std::endl << blk_curr_hold << std::endl;

		blk_curr_hold /= sqrt(resist_hold); red_curr_hold /= sqrt(resist_hold);
		blk_curr_v += blk_curr_hold; red_curr_v += red_curr_hold;
	//	if (debug_mode > 2) std::cout << "Black current 3 :" << std::endl << blk_curr_hold << std::endl;
	//	if (debug_mode > 2) std::cout << "Black current main :" << std::endl << blk_curr_a << std::endl;

		if (dyad % 1000 == 0) std::cout << "Done: " << dyad << " of " << dyads << std::endl;
	}

}

extern "C" {

	void get_resistances_fast_symmetric_c (
		const double * sociomatrix_c, 
		const double * fidelities_c,
		const int * nn_row_c,
		const int * src_c, const int * snk_c, const int * dyad_lengths_c, 

		double * equiv_resist_c, double * fidelity_c,
		double * blk_curr_a_c, double * red_curr_a_c,
		double * blk_curr_v_c, double * red_curr_v_c,
		double * blk_curr_p_c, double * red_curr_p_c,

		double * ll_save_1_c, double * ll_save_2_c, double * ll_save_3_c
	) {
		if (debug_mode) std::cout << "get_resistances_fast_symmetric_c entering." << std::endl;

		int nn_row = *nn_row_c, dyad_lengths = *dyad_lengths_c;
		int kk;
		Matrix<double> sociomatrix(nn_row, nn_row, sociomatrix_c);
		Matrix<double> fidelities(nn_row, nn_row, fidelities_c);
		Matrix<int> src(dyad_lengths, 1, src_c);
		Matrix<int> snk(dyad_lengths, 1, snk_c);

		Matrix<double> equiv_resist(dyad_lengths, 1, equiv_resist_c);
		Matrix<double> fidelity(dyad_lengths, 1, fidelity_c);

		Matrix<double> blk_curr_a(nn_row, nn_row, 1, 0);
		Matrix<double> red_curr_a(nn_row, nn_row, 1, 0);
		Matrix<double> blk_curr_v(nn_row, nn_row, 1, 0);
		Matrix<double> red_curr_v(nn_row, nn_row, 1, 0);
		Matrix<double> blk_curr_p(nn_row, nn_row, 1, 0);
		Matrix<double> red_curr_p(nn_row, nn_row, 1, 0);

		Matrix<double> ll_save_1(nn_row-1, nn_row-1, ll_save_1_c);
		Matrix<double> ll_save_2(nn_row-1, nn_row-1, ll_save_2_c);
		Matrix<double> ll_save_3(nn_row-1, nn_row-1, ll_save_3_c);
	
		if (debug_mode > 1) std::cout << "Fidelities:" << fidelities << std::endl;		

		if (debug_mode) std::cout << "get_resistances_fast_symmetric_c terms loaded." << std::endl;

		get_resistances_fast_symmetric (sociomatrix, fidelities, src, snk,
			equiv_resist, fidelity, blk_curr_a, red_curr_a, 
			blk_curr_v, red_curr_v, blk_curr_p, red_curr_p,
			ll_save_1, ll_save_2, ll_save_3);

		if (debug_mode) std::cout << "get_resistances_fast_symmetric complete." << std::endl;

		for (kk=0; kk<dyad_lengths; kk++) {
			equiv_resist_c[kk] = equiv_resist(kk); fidelity_c[kk] = fidelity(kk);
		}

		for (kk=0; kk<nn_row*nn_row; kk++) {
			blk_curr_a_c[kk] = blk_curr_a(kk);	blk_curr_v_c[kk] = blk_curr_v[kk];
			red_curr_a_c[kk] = red_curr_a(kk);	red_curr_v_c[kk] = red_curr_v[kk];
			blk_curr_p_c[kk] = blk_curr_p(kk);	red_curr_p_c[kk] = red_curr_p[kk];
		} 
		for (kk=0; kk<(nn_row-1)*(nn_row-1); kk++) {
			ll_save_1_c[kk] = ll_save_1(kk); ll_save_2_c[kk] = ll_save_2(kk); ll_save_3_c[kk] = ll_save_3(kk);
		}

	}
}




void get_resistances (
		Matrix<double> sociomatrix, Matrix<double> fidelities,  
		Matrix<int> sources, Matrix<int> sinks,
		Matrix<double> & equiv_resist, Matrix<double> & fidelity,
		Matrix<double> & blk_curr_a, Matrix<double> & red_curr_a, 
		Matrix<double> & blk_curr_v, Matrix<double> & red_curr_v,
		Matrix<double> & blk_curr_p, Matrix<double> & red_curr_p) {

	int nn_row = sociomatrix.rows();
	int rr,cc,ss, t_symm, t_src, t_snk; int total_asymm = 0;
	for (rr=0; rr<(nn_row-1); rr++) for (cc=(rr+1); cc<nn_row; cc++) {
		total_asymm += (1 * (sociomatrix(rr,cc)!=sociomatrix(cc,rr)) );
	}

	t_symm = 1*(total_asymm == 0) + total_asymm; Matrix<int> asymmetric(t_symm,2);
	Matrix<double> holder (nn_row-2,nn_row-2,1,0);

	Matrix<double> voltage_hold(nn_row,1);
	Matrix<double> current_hold(nn_row,1);
	double equiv_resist_hold, fidelity_hold;
	Matrix<double> average_current_hold(nn_row,nn_row,1,0);
	Matrix<double> red_current_hold(nn_row,nn_row,1,0);

	if (total_asymm > 0) {
		ss = 0;
		for (rr=0; rr<(nn_row-1); rr++) for (cc=(rr+1); cc<nn_row; cc++) if (sociomatrix(rr,cc)!=sociomatrix(cc,rr)) {
			asymmetric(ss,0) = rr; asymmetric(ss,1) = cc; ss++;
		}
	//	if (debug_mode) std::cout << asymmetric << std::endl;
	}

	for (rr=0; rr<equiv_resist.cols(); rr++) {
		t_src = sources(rr); t_snk = sinks(rr);
		if (total_asymm > 0) {
			solve_votes_notsymmetric_raw (sociomatrix, fidelities, t_src, t_snk,
				asymmetric, equiv_resist_hold, fidelity_hold, voltage_hold,	
				current_hold, holder, average_current_hold, red_current_hold);
		} else {
			solve_votes_symmetric (sociomatrix, fidelities, t_src, t_snk,
				equiv_resist_hold, fidelity_hold, voltage_hold,	current_hold, 
				holder, average_current_hold, red_current_hold);
		}
		if (debug_mode>2) std::cout << "g-res " << t_src << " " << t_snk << " " << equiv_resist_hold << std::endl;

		equiv_resist(rr) = equiv_resist_hold;
		fidelity(rr) = fidelity_hold;
		// results are in the form of constant voltage,

		blk_curr_v += average_current_hold; red_curr_v += red_current_hold; 

		//trouble if no current to start with...
		if (1/equiv_resist_hold > 0) {
			average_current_hold *= sqrt(equiv_resist_hold); 
			red_current_hold *= sqrt(equiv_resist_hold); 
			blk_curr_p += average_current_hold; red_curr_p += red_current_hold; 

			average_current_hold *= sqrt(equiv_resist_hold); 
			red_current_hold *= sqrt(equiv_resist_hold); 
			blk_curr_a += average_current_hold; red_curr_a += red_current_hold; 
		}

	//	if (voltage_hold(0)==-1) break;
		if ((rr+1) % 1000 == 0) std::cout << "Done: " << rr+1 << " " << nn_row << std::endl;

	}

}

extern "C" {

	void get_resistances_c (const double * sociomatrix_c, const double * fidelities_c, 
			const int * nn_row_c,
			const int * src_c, const int * snk_c, const int * dyad_lengths_c, 

			double * equiv_resist_c, double * fidelity_c,
		//	double * voltages_c, double * currents_c,

			double * blk_curr_a_c, double * red_curr_a_c,
			double * blk_curr_p_c, double * red_curr_p_c,
			double * blk_curr_v_c, double * red_curr_v_c) {

		int nn_row = *nn_row_c; //int saver = *save_results_c;
		int ii;
	//	int nn_t, dd_t; 
		int dyad_lengths = *dyad_lengths_c;

		Matrix<double> sociomatrix(nn_row,nn_row,sociomatrix_c);
		Matrix<double> fidelities(nn_row,nn_row,fidelities_c);
		Matrix<double> equiv_resist(1,dyad_lengths, 1, 0);
		Matrix<double> fidelity(1,dyad_lengths, 1, 0);

	//	nn_t = saver*nn_row + (1-saver);
	//	dd_t = saver*dyad_lengths + (1-saver);

	//	if (debug_mode > 0) std::cout << saver << " " << nn_t << " " << dd_t << std::endl;
		if (debug_mode > 0) std::cout << nn_row << " " << dyad_lengths << std::endl;
	
		Matrix<int> sources(1, dyad_lengths, src_c);
		Matrix<int> sinks(1, dyad_lengths, snk_c);
		Matrix<double> blk_curr_a(nn_row, nn_row, 1, 0);
		Matrix<double> red_curr_a(nn_row, nn_row, 1, 0);
		Matrix<double> blk_curr_p(nn_row, nn_row, 1, 0);
		Matrix<double> red_curr_p(nn_row, nn_row, 1, 0);
		Matrix<double> blk_curr_v(nn_row, nn_row, 1, 0);
		Matrix<double> red_curr_v(nn_row, nn_row, 1, 0);

		//	if (debug_mode) std::cout << "in g_r_c" << sociomatrix << std::endl;

		// get_resistances.
		get_resistances(sociomatrix, fidelities, sources, sinks, equiv_resist, fidelity, 
			blk_curr_a, red_curr_a, blk_curr_p, red_curr_p, blk_curr_v, red_curr_v);

	//	if (debug_mode > 0) std::cout << ":" << saver << " " << nn_t << " " << dd_t << std::endl;
	
		//	if (debug_mode) std::cout << equiv_resist << std::endl;

		// return.
		for (ii=0; ii<dyad_lengths; ii++) {
			equiv_resist_c[ii] = equiv_resist(ii); 
			fidelity_c[ii] = fidelity(ii);
		}	
		if (debug_mode > 0) std::cout << "post eq-fid" << std::endl;

/*		if (saver) {
			if (debug_mode) std::cout << "saving results" << std::endl;
			for (ii=0; ii < dyad_lengths*nn_row; ii++) {
				voltages_c[ii] = voltages(ii);
				currents_c[ii] = currents(ii);	
			}
		} else {
			voltages_c[0] = 0;
			currents_c[0] = 0;
		} */

		if (debug_mode > 0) std::cout << "post volt" << std::endl;
		for (ii=0; ii < nn_row*nn_row; ii++) {
			blk_curr_a_c[ii] = blk_curr_a(ii);
			blk_curr_p_c[ii] = blk_curr_p(ii);
			blk_curr_v_c[ii] = blk_curr_v(ii);
			red_curr_a_c[ii] = red_curr_a(ii);
			red_curr_p_c[ii] = red_curr_p(ii);
			red_curr_v_c[ii] = red_curr_v(ii);
		}
		if (debug_mode > 0) std::cout << "post-av-curr" << std::endl;


		//	if (debug_mode) std::cout << "done g_r_c" << std::endl;

	}




}



