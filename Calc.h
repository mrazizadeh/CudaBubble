#pragma once
#include "cufft.h"
#include "tool.h"
typedef double* DARR3D;



class Calc
{
public:
	Calc(int nx, int nz, double lenx);
	virtual ~Calc();

	void r2c_3d_f(DARR3D u, cufftDoubleComplex* s);
	void c2r_3d_b(cufftDoubleComplex* s, DARR3D u);

	void f(DARR3D f, DARR3D lap, DARR3D d, double eps);
	void g0(DARR3D g0, DARR3D d, DARR3D f, double eps);
	void intf(DARR3D f, DARR3D ret1, double* val, double addition, double multiply);
	void nld(DARR3D nldp, DARR3D d, DARR3D f, DARR3D tmp1, double eps, double kappa, double D);
	void lap(cufftDoubleComplex* fd, DARR3D out);
	void lap(DARR3D d, DARR3D out);
	void fftgradient(cufftDoubleComplex* uf, DARR3D udx, DARR3D udy, DARR3D udz);
	void s_term(DARR3D s, DARR3D d, DARR3D dx, DARR3D dy, DARR3D dz, double eps);

	void multi(DARR3D a, DARR3D b, DARR3D c);
	void plus(DARR3D u, DARR3D nld1, double t);

	void nld1(DARR3D nld, DARR3D f, DARR3D g, double lamb1, double lamb2, double eta);
	void nld2(DARR3D nld, DARR3D f, DARR3D g, double lamb1, double lamb2, double eta);
	void etd1(cufftDoubleComplex* d, cufftDoubleComplex * fd, cufftDoubleComplex* nld, DARR3D Het, double time_step);
	void etd2(cufftDoubleComplex* d, cufftDoubleComplex* fd, cufftDoubleComplex* nnld, cufftDoubleComplex* nld, DARR3D Het, double time_step);

	void etd41(DARR3D d, cufftDoubleComplex* fd, cufftDoubleComplex* nld, DARR3D Het, double time_step);
	void etd42(cufftDoubleComplex* b, cufftDoubleComplex* fd, cufftDoubleComplex* nld, DARR3D Het, double time_step);
	void etd43(cufftDoubleComplex* c, cufftDoubleComplex* fd, cufftDoubleComplex* nld, cufftDoubleComplex* nldb, DARR3D Het, double time_step);
	void etd44(DARR3D d, cufftDoubleComplex* fd, cufftDoubleComplex* nlda, cufftDoubleComplex* nldb, cufftDoubleComplex* nldc, DARR3D Het, double time_step);

	void init_fEIF(DARR3D Het, double eps);

	void fddx(DARR3D u, DARR3D s = NULL);
	void fddy(DARR3D u, DARR3D s = NULL);
	void fddx_b(DARR3D u, DARR3D s = NULL);
	void fddy_b(DARR3D u, DARR3D s = NULL);
	void fddxx(DARR3D u, DARR3D s = NULL);
	void fddyy(DARR3D u, DARR3D s = NULL);
	void fdlap(DARR3D u, DARR3D s = NULL);
	void fdgdsq(DARR3D u, DARR3D s = NULL);

	/*
	void fftdx(DARR3D u, DARR3D s, cufftDoubleComplex* ctmp);
	void fftdy(DARR3D u, DARR3D s);
	void fftdz(DARR3D u, DARR3D s);
	void z_fft_forward(DARR3D u, DARR3D s = NULL);
	void z_fft_backward(DARR3D u, DARR3D s = NULL);
	void z_fft_forward(ARR2D u, ARR2D s = NULL);
	void z_fft_backward(ARR2D u, ARR2D s = NULL);
	double Euler_Number1(DARR3D d);
	double Euler_Number(DARR3D d);
	double intf(fftw_complex* f);
	double intf(ARR1D f, ARR1D g = NULL);

	double L2norm_square(DARR3D u);
	void fftdz_b(DARR3D u, DARR3D s = NULL);
	void fftdzz(DARR3D u, DARR3D s = NULL);
	void fohmz2(double alp, DARR3D  u);
	void fftproj(DARR3D f, DARR3D g, DARR3D h);
	void fodz(DARR3D  u);
	void fody(DARR3D  u);
	void fodx(DARR3D u);
	void period(DARR3D u);
	void set_mesh(double* xx, double* yy, double* zz);
	void time_diff_solve_phi(DARR3D phi, double a, double b, double c, DARR3D f0, DARR3D f1, double t);
	void init_fEIF(DARR3D Het, double eps);
	void c2r(fftw_complex* c, DARR3D r);
	void c2i(fftw_complex* c, DARR3D r);
	void r2c(DARR3D r, fftw_complex* c);
	void fft_3d_b(fftw_complex* u, fftw_complex* s = NULL);
	void fft_3d_f(fftw_complex* u, fftw_complex* s = NULL);
	void expend(DARR3D d, DARR3D w);

	FdHelmhz* fd;
	*/

private:

	cufftHandle planf, planb;
	void laplace(cufftDoubleComplex* u);

	//cufftDoubleComplex* uf, * uf_copy;
	//cufftDoubleReal* udx, * udy, * udz;
	/*
	double euler(DARR3D f, DARR3D d, DARR3D w, double a, double b);
	void r2r_2d_f(ARR2D u, ARR2D s = NULL, int pind = 0); //only used for x-y plan
	void r2r_2d_b(ARR2D u, ARR2D s = NULL, int pind = 0); //only used for x-y plan

	//fftw_plan plan_f, plan_b;
	fftw_plan* plan_fx, * plan_bx, * plan_fy, * plan_by, * plan_fz, * plan_bz;
	fftw_plan* plan_fxy, * plan_bxy;

	fftw_plan plan_3d_f, plan_3d_b;
	fftw_complex* fftbuff_3d_in, * fftbuff_3d_out;


	ARR2D fftbuff;
	ARR2D fdbuff;
	DARR3D fftbuffxy;
	int nproc;


	DARR3D w1, w2, w3;   //tmp working array

	DARR3D bigw; //big working array
	*/
	int nx, ny, nz;
	double xlen, ylen, zlen;   //mesh size
	cufftDoubleComplex* ctmp1, * ctmp2, * ctmp3;
	//cufftDoubleComplex* fftbackcopy;  //special used for fftback
};

