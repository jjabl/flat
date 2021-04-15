#ifndef __METRICS_HPP_
#define __METRICS_HPP_

const double eps = 0.00001;

double __wasserstein1_norm(double *pos1, double *mass1, int size1, double scale1, double *pos2, double *mass2, int size2, double scale2);
double wasserstein1(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double wasserstein1_norm(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double radon(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double __flat(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double __flat_nlogn(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double __flat_n(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double flat(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double flat_nlogn(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double flat_n(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);
double wasserstein1_restricted(double *pos1, double *mass1, int size1, double *pos2, double *mass2, int size2);

class fentry {
	public:
		double x;
		double a;

		fentry() { x=a=0; }
		fentry(double ina) { x=0.; a=ina; }
		fentry(double inx, double ina) {x=inx;a=ina;}

		bool operator<(const fentry& r) const { return (a<r.a); }
};

#endif
