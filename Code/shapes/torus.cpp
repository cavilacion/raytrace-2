/* Finding intersections of a torus with a line in 3D requires solving
 * a fourth degree polynomial. We need complex numbers for this.
 */

#include "torus.h"

Torus::Torus(Point const &C, double R, double r)
:
  C (C),
	R (R),
	r (r)
{}

/* We use vector equations as derived in
 * https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
 */
Hit Torus::intersect(Ray const &ray)
{
	Vector x = ray.O - C;
	Vector s = ray.D.normalized();
	
	/* helper variables for polynomial coefficients */
	double beta = 2*s.dot(x);
	double gamma = x.dot(x) + R*R - r*r;
	
	/* polynomial coefficients */
	double *c = new double [4];
	c[0] = 2.0*beta;
	c[1] = beta*beta + 2.0*gamma  - 4.0*sqr(R)*(s.x*s.x+s.y*s.y);
	c[2] = 2.0*beta*gamma - 8.0*sqr(R)*(x.x*s.x+x.y*s.y);
	c[3] = sqr(gamma) - 4.0*sqr(R)*(sqr(x.x)+sqr(x.y));
	double t=0.0;
	int n = solveQuarticReal (t, c);
	delete c;
	if (n==0) { // negative or nonreal solutions
		return Hit::NO_HIT();
	}
	Point P = ray.at(t) - C;
	double norm=sqrt(sqr(P.x) + sqr(P.y));
	Vector dir (P.x/norm,P.y/norm,0.0);
	Point C2 = dir*R;
	Vector N = (P-C2).normalized();
  return Hit(t,N);
}

/* Function solves a quartic polynomial and returns 0 if equation has 
 * no real positive roots. 
 * The smallest positive real solution will will be in t.
 */
int solveQuarticReal (double &t, const double *c) {
	std::complex<double> *z;
	std::complex<double> *a=new std::complex<double> [4];
	for (int i=0; i<4; i++) {
		a[i] = std::complex<double>(c[i],0);
	}
	z = solveQuarticComplex (a);
	double *x= new double [4];
	int n = numRealSolutions (c);
	if (n==0) {
		delete x;
		delete z;
		delete a;
		return 0;
	}
	for (int i=0; i<n; i++) {
		x[i] = z[i].real();
	}
	double minimum=infty;
	int flag = 0;
	for (int i=0; i<n; i++) {
		if (x[i] < minimum && x[i] > 0.0) {
			flag = 1; // at least one is positive
			minimum = x[i];
		}
	}
	if (flag) {
		t = minimum;
	}
	else {
		// solutions were negative
		n = 0;
	}
	delete z;
	delete a;
	delete x;
	return n;
}		

/* Returns the number of real solutions of a quartic polynomial.
 * Can be zero, two, or four 
 */
int numRealSolutions (const double *a) {
	double b,c,d,e;
	b = a[0];
	c = a[1]; 
	d = a[2]; 
	e = a[3];
	/* Discriminant values and logic from
	 * https://en.wikipedia.org/wiki/Quartic_function
	 */
	double delta = 256.0*e*e*e-192.0*b*d*e*e-128.0*c*c*e*e+144.0*c*d*d*e-27.0*d*d*d*d 
	               +144.0*b*b*c*e*e-6.0*b*b*d*d*e-80.0*b*c*c*d*e+18.0*b*c*d*d*d+16.0*c*c*c*c*e 
	               -4.0*c*c*c*d*d-27.0*b*b*b*b*e*e+18.0*b*b*b*c*d*e-4.0*b*b*b*d*d*d-4.0*b*b*c*c*c*e+b*b*c*c*d*d;
	double P = 8.0*c - 3.0*b*b;
	double D = 64.0*e - 16.0*c*c + 16.0*b*b*c - 16.0*b*d - 3.0*b*b*b*b;
	double del0 = c*c - 3.0*b*d + 12.0*e;
	double R = b*b*b + 8.0*d - 4.0*b*c;
	if (delta < 0.0) {
		return 2;
	}
	if (delta > 0.0) {
		if (P < 0.0 && D < 0.0) return 4;
		if (P > 0.0 || D > 0.0) return 0;
	}
	/* very rare case detla = 0 */
	if (P < 0.0 && D < 0.0 && fabs(del0) > 0.0) return 4; 
	if (D > 0.0 || (P > 0.0 && (fabs(D)>0.0 || fabs(R)>0.0))) return 2; 
	if (fabs(del0) < 1.0E-10 && fabs(D) > 1.0E-10) return 4;
	if (fabs(D)>0.0) {
		if (P < 0.0) return 4;
		if (P > 0.0 && fabs(R) < 1.0E-10) return 0;
		if (fabs(del0) < 1.0E-10) return 4;
	}
	
	std::cout << "warning: uncovered case.\n";
	return 0;
}

/* Solving a fourth degree polynomial requires taking cubic roots of
 * complex numbers
 */
std::complex<double> cuberoot(std::complex<double> z) {
    if (z.real() < 0) {
        return -pow(-z, 1.0/3.0);
    } else {
        return pow(z, 1.0/3.0);
    }
}

std::complex<double> *solveQuadratic (const std::complex<double> *a) {
	std::complex<double> *x=new std::complex<double>[2];
	/* x^2 + a[0]*x + a[1] = 0 */
	std::complex<double> D = sqr(a[0])-4.0*a[1];
	x[0] = (-a[0]+sqrt(D))/2.0;
	x[1] = (-a[0]-sqrt(D))/2.0;
	return x;
}
	
std::complex<double> *solveQuarticComplex (const std::complex<double> *a) {
	std::complex<double> *x=new std::complex<double>[4];
	/* x^4 + a[0]x^3 + a[1]x^2 + a[2]x + a[3] = 0*/
	/* x = y - a[0]/4 */
	std::complex<double> b[3];
	b[0] = (a[1]-6.0*sqr(a[0])/16.0)/2.0;
	b[1] = cub(a[0])/8.0 - a[0]*a[1]/2.0 + a[2];
	b[2] = -3.0*qrt(a[0])/256.0 + a[1]*sqr(a[0])/16.0 - a[0]*a[2]/4.0 + a[3];
	/* y^4 + 2*b[0]*y^2 + b[1]*y + b[2] = 0 */
	std::complex<double> c[3];
	c[0] = 2.0*b[0];
	c[1] = sqr(b[0])-b[2];
	c[2] = -sqr(b[1])/8.0;
	/* u^3 + c[0]*u^2 + c[1]*u + c[2] = 0 */
	std::complex<double> d[2];
	d[0] = c[1] - sqr(c[0])/3.0;
	d[1] = 2.0*cub(c[0])/27.0 - c[1]*c[0]/3.0 + c[2];
	/* v^3 + d[0]*v + d[1] = 0 */
	std::complex<double> e[2];
	e[0] = d[1];
	e[1] = -cub(d[0])/27.0;
	/* w^2 + e[0]*w + e[1] = 0 */
	std::complex<double> *w=solveQuadratic (e);
	/* find one of three solutions for u */
	std::complex<double> u;
	u=cuberoot(w[0])+cuberoot(w[1])-c[0]/3.0;
	delete[] w;
	std::complex<double> f[2], g[2];
	f[0] = -sqrt(2.0*u);
	g[0] = sqrt(2.0*u);
	f[1] = b[0]+u+b[1]/sqrt(8.0*u);
	g[1] = b[0]+u-b[1]/sqrt(8.0*u);
	std::complex<double> *m=solveQuadratic (f);
	std::complex<double> *n=solveQuadratic (g);
	x[0] = m[0]; x[1] = m[1];
	x[2] = n[0]; x[3] = n[1];
	delete[] m;
	delete[] n;
	/* sort roots ascending on imaginary part */
	for (int i=0; i<4; i++) {
		for (int j=i+1; j<4; j++) {
			if (fabs(x[j].imag()) < fabs(x[i].imag())) {
				x[i] = x[i]+x[j];
				x[j] = x[i]-x[j];
				x[i] = x[i]-x[j];
			}
		}
		x[i]-=a[0]/4.0;
	}
	/* If some but not all solutions are real, they are in x[0] and x[1] */
	return x;
}



