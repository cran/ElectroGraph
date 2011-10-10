/* 
#####################################################################
# electrograph-routines.h 
#
# Andrew C. Thomas <acthomas@stat.cmu.edu>
# Last modified: April 19, 2011
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Headers for electrograph-routines.cpp
#
#####################################################################


void solve_votes_symmetric (
	Matrix<double> sociomatrix,	int src,	int snk,
	double & resist_eq,	double & ref_mean,
	Matrix<double> & voltages,	Matrix<double> & currents);

void solve_votes_notsymmetric_raw (
		Matrix<double> sociomatrix,	int src,	int snk,
		Matrix<int> asymmetric_entries, // k by 2.
		double & resist_eq,	double & ref_mean,
		Matrix<double> & voltages,	Matrix<double> & currents);

extern "C" {
	void dijkstra_position_valued (const double * sociomatrix_c,		
		const int * src_c,	const int * snk_c,
		const int * nn_row_c,	double * voltages);
}*/

int symmetric (Matrix<double> sym_obj);

