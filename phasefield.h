// phasefield.h: interface for the phasefield class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "tool.h"	// Added by ClassView
#include <cufft.h>
#include "Calc.h"


class phasefield  
{
	friend class phasefield;
public:
	double *surface();
	double *vol();
	double* elastic_eng();
	double* LSerror(DARR3D u);
	DARR3D s();
	DARR3D d();
	cufftDoubleComplex* fd(bool calulate = true);
	DARR3D lap();
	DARR3D lapf();
	DARR3D f();
	DARR3D fc();
	DARR3D g0();
	DARR3D g();
	void gradient();
	DARR3D nld();
	void update_d();
	void update_fd(cufftDoubleComplex* fd);
	phasefield(ARR3D d, int nx1, int nz1, double xlen, double eps1, double xi1, ARR3D C, Calc* calcp);
	virtual ~phasefield();
	void do_fd(bool b) {dofft = !b;};
	void copy_array(DARR3D s, DARR3D d);
	void time_diff_solve_phi(DARR3D nld1, double t);

private:
	bool dofft;  //fft or fd
	int nx, ny, nz;
	double eps, h, xi;
	DARR3D m_C;
	cufftDoubleComplex* m_fd;
	cufftDoubleComplex* cptmp1;
	DARR3D m_d, m_f, m_fc, m_g, m_g0, m_lap, m_lapf, m_dx, m_dy, m_dz, m_s, m_nld;
	DARR3D tmp1, tmp2, tmp3;
	bool b_fd, b_f, b_fc, b_g0, b_g, b_q, b_lap, b_lapf, b_dx, b_dy, b_dz, b_s, b_t, b_tanh, b_tanhdiff, b_vf, b_nld;
	double * m_cur_vol, *m_cur_surf, * m_Q;  //on device
	double* m_elastic_eng;
	double* m_LSerror;
	bool b_cur_vol, b_cur_surf, b_Q, b_elastic_eng;
	void serial_update_d();

	void copy_array(cufftDoubleComplex* d, cufftDoubleComplex* s);

	Calc* calc;

};

