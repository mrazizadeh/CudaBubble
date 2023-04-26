// phasefield.cpp: implementation of the phasefield class.
//
//////////////////////////////////////////////////////////////////////

#include "phasefield.h"
#include <cuda_runtime.h>
#include <stdio.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//input example nx1 = 64, nz1 = 64.
phasefield::phasefield(ARR3D d, int nx1, int nz1, double xlen, double eps1, double xi1, ARR3D C, Calc* calcp)
{
	dofft = true;
	nx = nx1; ny = nx1; 
	nz = nz1;
	h = xlen / nx1;
	//nx = ny = nx1 + 1;
	calc = calcp;

	cudaMalloc(&m_d, nx * ny * nz * sizeof(double));
	if (d) cudaMemcpy(m_d, &d[0][0][0], nx * ny * nz * sizeof(double), cudaMemcpyHostToDevice);
	serial_update_d();

	cudaMalloc(&m_fd, (nx / 2 + 1) * ny * nz * sizeof(cufftDoubleComplex));

	m_f = m_fc = NULL;
	m_g0 = m_g = NULL;
	m_lap = NULL;
	m_lapf = NULL;
	m_dx = m_dy = m_dz = NULL;
	m_s = NULL;
	m_nld = NULL;

	cudaMalloc(&tmp1, nx * ny * nz * sizeof(double));
	cudaMalloc(&tmp2, nx * ny * nz * sizeof(double));
	cudaMalloc(&tmp3, nx * ny * nz * sizeof(double));

	cudaMalloc(&m_cur_vol, sizeof(double));
	cudaMalloc(&m_cur_surf, sizeof(double));
	cudaMalloc(&m_Q, sizeof(double));
	cudaMalloc(&m_elastic_eng, sizeof(double));
	cudaMalloc(&m_LSerror, sizeof(double));

	cudaMalloc(&cptmp1, (nx/2+1) * ny * nz * sizeof(cufftDoubleComplex));

	if (C) {
		cudaMalloc(&m_C, nx * ny * nz * sizeof(double));
		cudaMemcpy(m_C, &C[0][0][0], nx * ny * nz * sizeof(double), cudaMemcpyHostToDevice);
	}
	else {
		m_C = NULL;
	}
	eps = eps1;
	xi = xi1;
}

phasefield::~phasefield()
{
	if (m_C) cudaFree(m_C);
	if (tmp3) cudaFree(tmp3);
	if (tmp2) cudaFree(tmp2);
	if (tmp1) cudaFree(tmp1);
	if (m_nld) cudaFree(m_nld);
	if (m_s) cudaFree(m_s);
	if (m_dz) cudaFree(m_dz);
	if (m_dy) cudaFree(m_dy);
	if (m_dx) cudaFree(m_dx);
	if (m_lap) cudaFree(m_lap);
	if (m_lapf) cudaFree(m_lapf);
	if (m_g) cudaFree(m_g);
	if (m_g0) cudaFree(m_g0);
	if (m_fc) cudaFree(m_fc);
	if (m_f) cudaFree(m_f);
	if (m_fd) cudaFree(m_fd);
	if (m_d) cudaFree(m_d);

	cudaFree(cptmp1);

	cudaFree(m_cur_vol);
	cudaFree(m_cur_surf);
	cudaFree(m_Q);
	cudaFree(m_elastic_eng);
	cudaFree(m_LSerror);
}

void  phasefield::time_diff_solve_phi(DARR3D nld1, double t)
{
	calc->plus(m_d, nld1, t);
	update_d();
}

void phasefield::serial_update_d()
{
	b_fd = b_f = b_fc = b_g = b_g0 = b_q = b_lap = b_lapf = b_dx = b_dy = b_dz = b_s = b_t = b_tanh = b_tanhdiff = b_vf = b_nld = false;
	b_cur_surf = b_cur_vol = b_Q = b_elastic_eng = false;
}

void phasefield::update_d()
{
#pragma omp barrier
#pragma omp single
	serial_update_d();
}

void phasefield::copy_array(cufftDoubleComplex* d, cufftDoubleComplex* s)
{
	cudaMemcpy(d, s, nz * ny * (nx / 2 + 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToDevice);
}

void phasefield::copy_array(DARR3D d, DARR3D s)
{
	cudaMemcpy(d, s, nz * ny * nx * sizeof(double), cudaMemcpyDeviceToDevice);
}

//fd should be different with m_d before use. If they are the same, need copy to another one for safe.
void phasefield::update_fd(cufftDoubleComplex* fd)
{
	calc->c2r_3d_b(fd, m_d);
	serial_update_d();
}

DARR3D phasefield::d()
{
	return m_d;
}

cufftDoubleComplex* phasefield::fd(bool calulate)
{
	if (!calulate) return m_fd;

	if (b_fd) return m_fd;
	calc->r2c_3d_f(m_d, m_fd);
	b_fd = true;
	return m_fd;
}

DARR3D phasefield::lap()
{
	if (!m_lap) cudaMalloc(&m_lap, nx * ny * nz * sizeof(double));
	if (b_lap) return m_lap;
	if (dofft) {
		fd();
		copy_array(cptmp1, m_fd);
		calc->lap(cptmp1, m_lap);
	}
	else {
		calc->fdlap(m_d, m_lap);
	}

	b_lap = true;
	return m_lap;
}

DARR3D phasefield::lapf()
{
	if (!m_lapf) cudaMalloc(&m_lapf, nx * ny * nz * sizeof(double));
	if (b_lapf) return m_lapf;
	f();
	if (dofft) {
		calc->lap(m_f, m_lapf);
	}
	else {  //to do...
		//calc->fdlap(m_lapf);
	}

	b_lapf = true;
	return m_lapf;
}

DARR3D phasefield::f()
{
	if (!m_f) cudaMalloc(&m_f, nx * ny * nz * sizeof(double));
	if (b_f) return m_f;
	lap();
	calc->f(m_f, m_lap, m_d, eps);
	b_f = true;
	return m_f;
}

DARR3D phasefield::fc()
{
	if (m_C == NULL) return f();
	if (!m_fc) cudaMalloc(&m_fc, nx * ny * nz * sizeof(double));
	if (b_fc) return m_fc;
	f();
	/*  to do...
	//calc->fc(m_fc, m_f, m_d, m_C);
	int i, j, k;
#pragma omp for 
	for (i = 0; i <nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				m_fc[i][j][k] = m_f[i][j][k] + (m_d[i][j][k]*m_d[i][j][k] - 1.0)*m_C[i][j][k];
			}
		}
	}
	*/
	b_fc = true;
	return m_fc;
}

DARR3D phasefield::g0()
{
	if (!m_g0) cudaMalloc(&m_g0, nx * ny * nz * sizeof(double));
	if (b_g0) return m_g0;
	f();
	if (dofft) {
		calc->lap(m_f, m_g0);
	}
	else { // to do...
		//calc->fdlap(m_f, m_g0);
	}

	//m_g0[i][j][k] = -m_g0[i][j][k] + (3.0*m_d[i][j][k]*m_d[i][j][k] - 1.0)*m_f[i][j][k] / (eps*eps);
	calc->g0(m_g0, m_d, m_f, eps);
	b_g0 = true;
	return m_g0;
}

DARR3D phasefield::g()
{
	if (m_C == NULL) return g0();

	if (!m_g) cudaMalloc(&m_g, nx * ny * nz * sizeof(double));
	if (b_g) return m_g;
	if (dofft) {
		calc->lap(fc(), m_g);
	}
	else calc->fdlap(fc(), m_g);

	/*

	double tmp = eps*eps;
	int i, j, k;
#pragma omp for 
	for (i = 0; i <nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				m_g[i][j][k] = -m_g[i][j][k] + 
					(3.0*m_d[i][j][k]*m_d[i][j][k] - 1.0 + 2.0*eps*m_d[i][j][k]*m_C[i][j][k])*m_fc[i][j][k] / tmp;
			}
		}
	}
	*/
	b_g = true;
	return m_g;
}

//nld = -1/eps*Laplace(phi^3 - (1+k)phi) + 1/eps^3 *D*(phi^3 - (1+k)phi) - 1/eps(3phi^2-D-1)(Laplace phi - 1/eps^2(phi^3-phi))
DARR3D phasefield::nld()
{
	if (!m_nld) cudaMalloc(&m_nld, nx * ny * nz * sizeof(double));

	if (b_nld) return m_nld;

	double kappa = 2.0;
	double D = 1.5;

	f();
	calc->nld(m_nld, m_d, m_f, tmp1, eps, kappa, D);

	b_nld = true;
	return m_nld;
}


void phasefield::gradient()
{
	if (!m_dx) cudaMalloc(&m_dx, nx * ny * nz * sizeof(double));
	if (!m_dy) cudaMalloc(&m_dy, nx * ny * nz * sizeof(double));
	if (!m_dz) cudaMalloc(&m_dz, nx * ny * nz * sizeof(double));

	if (b_dx && b_dy && b_dz) return;
	fd();
	calc->fftgradient(m_fd, m_dx, m_dy, m_dz);
	b_dx = b_dy = b_dz = true;

}

double * phasefield::vol()
{
	if (b_cur_vol) return m_cur_vol;
	calc->intf(m_d, tmp1, tmp2, m_cur_vol, 1.0, 0.5);
	b_cur_vol = true;
	return m_cur_vol;
}



DARR3D phasefield::s()
{
	if (!m_s) cudaMalloc(&m_s, nx * ny * nz * sizeof(double));
	if (b_s) return m_s;

	if (dofft) {
		gradient();
		calc->s_term(m_s, m_d, m_dx, m_dy, m_dz, eps);

	}
	else {
	}

	b_s = true;
	return m_s;
}

double* phasefield::surface()
{
	if (b_cur_surf) return m_cur_surf;
	calc->intf(s(), tmp1, tmp2, m_cur_surf, 0.0, 3.0 / (2.0 * sqrt(2.0)));
	b_cur_surf = true;
	return m_cur_surf;
}

double* phasefield::elastic_eng()
{
	if (b_elastic_eng) return m_elastic_eng;
	calc->multi(fc(), fc(), tmp1);
	calc->intf(tmp1, tmp2, tmp3, m_elastic_eng, 0.0, 1.0 / (eps * 2));
	b_elastic_eng = true;
	return m_elastic_eng;
}

double* phasefield::LSerror(DARR3D u)
{
	calc->multi(u, u, tmp1);
	calc->intf(tmp1, tmp2, tmp3, m_LSerror, 0.0, 1.0);
	return m_LSerror;
}

