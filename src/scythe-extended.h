// Scythe-based matrix random number generation
// Andrew C. Thomas
// January 20, 2009

//#include "scythe-local.h"
#include <stdint.h>
#include "scythestat-act/distributions.h" //	Definitions for probability density functions (PDFs), cumulative distribution functions (CDFs), and some common functions (gamma, beta, etc) */
//#include "scythestat/ide.h" 
//#include "scythestat/la.h" // manipulations on Scythe Matrix objects
#include "scythestat-act/matrix.h" //	Definitions of Matrix and related classes and functions
//#include "scythestat/optimize.h" //	Definitions of functions for doing numerical optimization and related operations
#include "scythestat-act/rng.h" //	The definition of the random number generator base class
#include "scythestat-act/smath.h" //	Definitions for functions that perform common mathematical operations on every element of a Matrix
#include "scythestat-act/stat.h" //	Definitions for functions that perform common statistical operations on Scythe Matrix objects
#include "scythestat-act/rng-mersenne.h" //	The Mersenne Twister random number generator
//#include "scythestat/rng/rtmvnorm.h" //	A truncated multivariate normal random number generator

#include <iostream>
#include <sys/time.h>
#include <cmath>

using namespace scythe;
//namespace scythe {

	// expand matrices to fit dimensions when necessary.

	template <typename T1, typename T2>
	Matrix<T1> efts (Matrix<T1> obj_one, Matrix<T2> obj_two) { 
		//	expand first matrix to the dimensions of the second.
		//	doesn't need to be an even multiple of cells.

		Matrix<T1> obj_one_hold(obj_two.rows(),obj_two.cols());
		int index, size_one;

		size_one = obj_one.rows()*obj_one.cols();
		for (index = 0; index < obj_two.rows()*obj_two.cols(); index++) {
			obj_one_hold[index] = obj_one[index % size_one];
		}

		return obj_one_hold;
	}

	//assume that one matrix fits in another.
	template <typename T1, typename T2>
	void expand_to_fit (Matrix<T1> & obj_one, Matrix<T2> & obj_two) {
		//which is bigger?
		//int one_items = obj_one.rows()*obj_one.cols(); int two_items = obj_two.rows()*obj_two.cols();
		if (obj_one.rows()*obj_one.cols()>obj_two.rows()*obj_two.cols()) {
			obj_two = efts(obj_two,obj_one);} else {
			obj_one = efts(obj_one,obj_two);}
	}

	template <typename T1,typename T2,typename T3>
	void expand_to_fit (Matrix<T1> & obj_one, Matrix<T2> & obj_two, Matrix<T3> & obj_three) {
		int one_items = obj_one.rows()*obj_one.cols(); 
		int two_items = obj_two.rows()*obj_two.cols();
		int three_items = obj_three.rows()*obj_three.cols();

		if (one_items>=two_items && one_items>=three_items) {
			obj_two = efts(obj_two,obj_one);
			obj_three = efts(obj_three,obj_one);
		} else if (two_items>=one_items && two_items>=three_items) {
			obj_one = efts(obj_one,obj_two);
			obj_three = efts(obj_three,obj_two);
		} else {
			obj_one = efts(obj_one,obj_three);
			obj_two = efts(obj_two,obj_three);
		}
	}

	
	// RNG routines that aren't standard with Scythe, but are needed here.

	// cribbed from R.
	unsigned int seed_value () {
		struct timeval tv;
		gettimeofday (&tv, NULL);
		unsigned int seed = ((int64_t) tv.tv_usec << 16) ^ tv.tv_sec;
		return seed;
	}

	// Bernoulli generation with different p in each cell.
	template <typename RNGTYPE>
	Matrix<> Berns (Matrix<> ps, rng<RNGTYPE>& rands) {
		Matrix<> output(ps.rows(),ps.cols());
		for (size_t rr = 0; rr < ps.rows(); rr++) for (size_t cc = 0; cc < ps.cols(); cc++)	output(rr,cc) = rands.rbern(ps(rr,cc));

		return output;
	}

	template <typename RNGTYPE>
	Matrix<> rgamma_big (Matrix<> shape, Matrix<> rate, rng<RNGTYPE>& rands) {
		Matrix<> output (shape.rows(),shape.cols());
		expand_to_fit (shape,rate);
		std::cout << ".";
	//	std::cout << "etf " << shape.rows() << " " << shape.cols() << " " << rate.rows() << " " << rate.cols() << std::endl;

		for (size_t rr = 0; rr < shape.rows(); rr++) for (size_t cc = 0; cc < shape.cols(); cc++)	{
			if (shape(rr,cc) <= 0 || rate(rr,cc) <= 0) std::cout << rr << " " << cc << ": shape=" << shape(rr,cc) << ",rate=" << rate(rr,cc) << std::endl;
			output(rr,cc) = rands.rgamma(shape(rr,cc),rate(rr,cc));
		}
		return output;	
	}

// Generates a Matrix of truncated normals with sd=1.
	template <typename RNGTYPE>
	Matrix<> rtnorm_big (Matrix<> mu, Matrix<> below, rng<RNGTYPE>& rands) {
		Matrix<> output(mu.rows(),mu.cols());
		for (size_t rr = 0; rr < mu.rows(); rr++) for (size_t cc = 0; cc < mu.cols(); cc++)	{
			output(rr,cc) = (below(rr,cc)?rands.rtbnorm_combo(mu(rr,cc),1,0):rands.rtanorm_combo(mu(rr,cc),1,0));  //rands.rbern(ps(rr,cc));
		}
		return output;
	}


	Matrix<> lndbern_big (Matrix<> yy, Matrix<> pp) {
		expand_to_fit(yy,pp);
		Matrix<> output(pp.rows(),pp.cols());
		for (size_t rr = 0; rr < pp.rows(); rr++) for (size_t cc = 0; cc < pp.cols(); cc++)	
			output(rr,cc) = log(dbinom(yy(rr,cc),1,pp(rr,cc)));
		return output;
	}

	Matrix<> lndnorm_big (Matrix<> x, Matrix<> mu, Matrix<> sigma) {
		expand_to_fit(x,mu,sigma);
		Matrix<> output(x.rows(),x.cols());
		for (size_t rr = 0; rr < x.rows(); rr++) for (size_t cc = 0; cc < x.cols(); cc++)	
			output(rr,cc) = lndnorm(x(rr,cc),mu(rr,cc),sigma(rr,cc));
		return output;	
	}

	Matrix<> pnorm_big (Matrix<> x, Matrix<> mu, Matrix<> sigma) {
		expand_to_fit(x,mu,sigma);
		double zscore;
	//	std::cout << "pnorm" << x.rows() << "-" << x.cols() << "," << mu.rows() << "-" << mu.cols() << "," << sigma.rows() << "-" << sigma.cols() << std::endl;
		Matrix<> output(x.rows(),x.cols());	
		for (size_t rr = 0; rr < x.rows(); rr++) for (size_t cc = 0; cc < x.cols(); cc++) {	
			zscore = (x(rr,cc)-mu(rr,cc))/sigma(rr,cc);
			if (abs(zscore)<=7) output(rr,cc) = pnorm1(zscore,1,0); else output(rr,cc) = (zscore>0?pnorm1(7,1,0):pnorm1(-7,1,0));
		}
	//	std::cout << "pnorm2" << std::endl;
		return output;	
	}

	double dinvgamma (double x, double a, double b, int log_out=1) {
		double quant = a*log(b)-lngammafn(a)-(a+1)*log(x)-b/x;
		if (log_out) return quant; else return exp(quant);
	}

	Matrix<> ev_tnorm (Matrix<> mu, Matrix<int> below) {
		//check if dimensions of mu and below match? below == "truncated below at 0"
		expand_to_fit(mu,below);
		Matrix<> outcome(mu.rows(),mu.cols());
		//std::cout << "ev1" << std::endl;
		for (size_t rr = 0; rr < mu.rows(); rr++) for (size_t cc = 0; cc < mu.cols(); cc++) {
			//if (abs(mu(rr,cc)>=8)) std::cout << mu(rr,cc) << " ";
			outcome(rr,cc) = mu(rr,cc) + (below(rr,cc)==0?
							-exp(lndnorm(-mu(rr,cc),0,1)-pnorm1(-mu(rr,cc),1,1)):   //pnorm1: x-mu/sigma,lower tail, log?
							exp(lndnorm(-mu(rr,cc),0,1)-pnorm1(mu(rr,cc),1,1)));
		}
		//std::cout << "ev2" << std::endl;
		return outcome;
	}

//}


