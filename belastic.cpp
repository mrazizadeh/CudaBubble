// belastic.cpp: implementation of the belastic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "belastic.h"
#include "mem_array.h"
#include <cuda_runtime.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

belastic::belastic(int folderindex):fluidbase(folderindex)
{
	calc = new Calc(nx, nz, xlen);
	init();
	m_pd = new phasefield(m_d, nx, nz, xlen, eps, xi, C, calc);
	m_pa = new phasefield(NULL, nx, nz, xlen, eps, xi, C, calc);
	output("Elastic Bending energy program!\n");
	cudaMalloc(&tmp, nx * ny * nz * sizeof(double));
}

belastic::~belastic()
{
	delete calc;
	delete m_pd;
	delete m_pa;

	cudaFree(tmp);
}

//integral of s*s
double belastic::intf(DARR3D s)
{
	mem_array mem(this);
	DARR3D d, ret1, ret2;
	d = mem.get_array();
	ret1 = mem.get_array();
	double ret;
	double* value;
	cudaMalloc(&value, sizeof(double));
	calc->multi(s, s, d);
	calc->intf(d, ret1, value, 0, 1);
	cudaMemcpy(&ret, value, sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(value);
	return ret;
}

//nld = g + lamb1 + lamb2 * f + lamb3 * q,.
void belastic::get_sum_fgqterm(phasefield * pd, DARR3D nld)
{
	double v, s;
	cudaMemcpy(&v, pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&s, pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
	double lamb1 = M1*(v - alpha);
	double lamb2 = M2*(s - beta);
	DARR3D f = pd->f();
	DARR3D g = pd->g();

	//nld[i][j][k] = lamb1 * 0.5 + lamb2 * f[i][j][k] + g[i][j][k] * eta;
	calc->nld1(nld, f, g, lamb1, lamb2, eta);
}


//static energy.
double belastic::getstaticenergy(phasefield * pd, bool constrained)
{
	double volerr, surerr;
	if (constrained) {
		double v, s;
		cudaMemcpy(&v, pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(&s, pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
		volerr = v - alpha;
		surerr = s - beta;
	}
	double eng1;
	cudaMemcpy(&eng1, pd->elastic_eng(), sizeof(double), cudaMemcpyDeviceToHost);
	eng1 *= eta;

	if (constrained) {
		return  eng1 + 0.5*M1 * volerr* volerr + 0.5*M2*surerr * surerr;   
	}
	else return eng1;
}

double belastic::outputenergy()
{
	double volerr, surerr;

	double v, s;
	cudaMemcpy(&v, m_pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&s, m_pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
	volerr = v - alpha;
	surerr = s - beta;
	double eng1;
	cudaMemcpy(&eng1, m_pd->elastic_eng(), sizeof(double), cudaMemcpyDeviceToHost);
	eng1 *= eta;

	double eng2 = 0;
	double eng3 = 0;
	double eng4 = 0.5*M1 * volerr* volerr;
	double eng5 = 0.5*M2*surerr * surerr;
	double eng6 = 0;
	double eng7 = eng1 + eng2 + eng3 + eng4 + eng5 + eng6;
	double eng8 = 0;
	double eng9 = eng7 + eng8;

	//output: 1, elstic bending energy, 2, weighted energy, 3, fix part energy, 4, volume fix energy, 5, surface fix energy, 
	//		6, tanh energy, 7, total phase field energy, 8, fluid energy, 9, total energy)
#pragma omp single
	output("energy", eng1, eng4, eng5, eng7); 
	return  eng9;  
}

double belastic::gettotalenergy()
{
	double volerr, surerr;

	double v, s;
	cudaMemcpy(&v, m_pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&s, m_pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
	volerr = v - alpha;
	surerr = s - beta;
	double eng1;
	cudaMemcpy(&eng1, m_pd->elastic_eng(), sizeof(double), cudaMemcpyDeviceToHost);
	eng1 *= eta;

	double eng4 = 0.5*M1 * volerr* volerr;
	double eng5 = 0.5*M2*surerr * surerr;

	return  eng1+eng4+eng5;  
}

double belastic::geterror(phasefield* pd)
{
	get_sum_fgqterm(pd, tmp);
	double* err = pd->LSerror(tmp);
	double ret;
	cudaMemcpy(&ret, err, sizeof(double), cudaMemcpyDeviceToHost);
	return ret;
}



void belastic::outputall()
{
	double vol1, sur1;
	cudaMemcpy(&vol1, m_pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&sur1, m_pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
	double lamb1 = M1*(vol1 - alpha);
	double lamb2 = M2*(sur1 - beta);
	double err = geterror(m_pd);
	double eng = outputenergy();
		output("vol_sur", vol1, sur1);
		output("lamb", lamb1, lamb2);
		output("error", err);
		output("time_step", time_step);
		output("alpha_beta", alpha, beta);

		output_info();
		cudaMemcpy(&m_d[0][0][0], m_pd->d(), nx * ny * nz * sizeof(double), cudaMemcpyDeviceToHost);
		output("d", m_d);
		output_graph("d", m_d, NULL, 0, 0, m_angle);
		
		sprintf(sout, "current time is %12.8f\n", time_current);
		output(sout);
		sprintf(sout, "vol: %12.7f surf: %12.7f eng: %12.10f Error: %12.10f Time_step: %e\n", vol1, sur1, eng, err, time_step);
		output(sout);
}

//run till tmax time to adjust the shape. Here is no fluid, no constains.
void belastic::shape_adjust()
{
	DARR3D od;

	cudaMalloc(&od, nx * ny * nz * sizeof(double));

	//the following are init parameters.
	double tcur = 0.0;

	double eng0, eng1, eng;
	int steps = 1;
	double old_time_step;
	old_time_step = time_step;
	time_step = 1e-7;

	//start parallel.
#pragma omp parallel 
	{
		//eng
		eng0 = getstaticenergy(m_pd, false);
		eng1 = eng0;
		eng = eng0;

#pragma omp single
		{
			output("Start Shape Adjusting......\n");
			sprintf(sout, "Start energy is %12.8f\n", eng0);
			output(sout);
		}

		for (;;) {
#pragma omp single
			{
				tcur += time_step;
				steps++;
			}
			//		if (steps % 10 == 0) {
			//			time_current = tcur;
			//			outputall();
			//		}

			if (time_step > 1e-4 || int(tcur / time_step) % int(1e-4 / time_step) == 0) {
				//			time_current = tcur;
				//			outputall();

				eng1 = getstaticenergy(m_pd, false);
#pragma omp single
				{
					sprintf(sout, "energy is %12.8f\n", eng1);
					output(sout);
				}
				if (eng1 < 100 * eta || fabs((eng - eng1) / eng) < 0.02) break;
#pragma omp barrier
				eng = eng1;
			}

#pragma omp barrier   //Important, please be care for previous judgement.
#pragma omp single
			{
				if (steps % 10 == 0) time_step *= 2;
			}
			m_pd->copy_array(od, m_pd->d());
			// ------------------START ------------------
			// 2. Update d
			DARR3D nld = m_pd->g0();  //to be safe, use copy_array.
			for (;;) {
				m_pd->time_diff_solve_phi(nld, -time_step * 3);
				eng1 = getstaticenergy(m_pd, false);

				if (eng1 < eng0) {
#pragma omp barrier
					eng0 = eng1;
					break;
				}
				else {
					m_pd->copy_array(m_pd->d(), od);
					m_pd->update_d();
#pragma omp single
					time_step /= 2;
				}
				if (time_step < 1e-13) break; //too small step, get exit.
			}

			if (time_step < 1e-13) {
#pragma omp single
				output("time_step in shape adjust < 1e-13, get exit.\n");
				break;
			}
			// 3. check for blow up!
			if (!_finite(m_d[nz / 2][ny / 2][nx / 2])) {
#pragma omp single
				{
					sprintf(sout, "Shape Adjust Blow up! Blow up at time %12.8f, %d steps\n", tcur, int(tcur / time_step));
					output(sout);
				}
				break;
			}

			// ------------------END ------------------
		}
		//end parallel.
	}

	time_step = old_time_step;
	output("Shape Adjust Ended successfully\n");

	//free all memory	

	cudaFree(od);
}


void belastic::run_static()
{
	//correct initial shape by trying run
	if (time_current == 0.0 && b_shape != 0) {
#ifndef _DEBUG
		shape_adjust();
#endif
	}

	DARR3D nld1; //non linear term -- all in spectral space.

	DARR3D od;

	cudaMalloc(&nld1, nx * ny * nz * sizeof(double));
	cudaMalloc(&od, nx * ny * nz * sizeof(double));

	//the following are init parameters.
	double eng0, eng1;
	int steps;
	time_t ltime1, ltime2, ltimec;

	//start parallel.
#pragma omp parallel 
	{

		//set alpha and beta
		/*
		int loop;
		for (loop = 0; loop < 100; loop++) {
			get_sum_fgqterm(m_pd, nld1);
			double a1 = calc->intf(nld1, nld1);
			double a2 = m_pd->surface();
			printf("alpha = %15.12f, beta = %15.12f\n", a1, a2);
			m_pd->update_d();
		}
	*/
		if (beta <= 0.0) {
			cudaMemcpy(&alpha, m_pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(&beta, m_pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
		}

#pragma omp single
		{
			sprintf(sout, "grid(xyz): %dx%dx%d, timestep = %e\n", nx, ny, nz, time_step);
			output(sout);
			if (fixed_time_step > 0) {
				time_step = min(fixed_time_step, time_step);
				sprintf(sout, "FIXED timestep = %e    ", time_step);
				output(sout);
			}
			else {
				output("Adjustment time_steps; ");
			}
			sprintf(sout, "eps = %9.6f, ga = %9.6f, eta = %9.6f, tmax = %9.6f\n", eps, ga, eta, tmax);
			output(sout);
			sprintf(sout, "alpha = %9.6f, beta = %9.6f, M1 = %8.2f, M2 = %8.2f, M3 = %8.2f, M4 = %8.2f, M5 = %8.2f, gweight = %8.2f\n", alpha, beta, M1, M2, M3, M4, M5, gweight);
			output(sout);

			output("********* Start STATIC Running ***********\n");
		}
		//get 0, 1 step values
		get_sum_fgqterm(m_pd, nld1);

		//eng
		eng0 = getstaticenergy(m_pd);
		eng1 = eng0;

		if (time_current == 0.0) {
#pragma omp single
			{
				sprintf(sout, "Preserving volume and surface. (Initial volume = %12.8f, surface = %12.8f\n", alpha, beta);
				output(sout);
			}
			outputall();
		}
#pragma omp single
		{
			time(&ltime1);
			ltimec = ltime1;
			steps = 0;
		}
		for (;;) {
			// ------------------START ------------------
			// 2. Update d
			//adjust time_step
			if (fixed_time_step <= 0) {
				eng0 = getstaticenergy(m_pd);
				m_pd->copy_array(od, m_pd->d());
#pragma omp single
				{
					if (steps % 4 == 0) {
						time_step = time_step * 1.4;
					}
				}
				for (;;) {
					if (time_step < 1e-23) break; //too small step, get exit.
					m_pd->time_diff_solve_phi(nld1, -time_step * ga);

					eng1 = getstaticenergy(m_pd);

					if (eng1 > eng0) {
						m_pd->copy_array(m_pd->d(), od);
						m_pd->update_d();
#pragma omp single
						time_step /= 2;
						continue;
					}
					else break;
				}
			}
			else {
				m_pd->time_diff_solve_phi(nld1, -time_step * ga);
			}

			//exit loop.
			if (time_step < 1e-23) {
#pragma omp single
				output("time_step < 1e-23, get exit.\n");
				outputall();
				break;
			}
			if (time_current > tmax) {
#pragma omp single
				output("Too long time running, get exit.\n");
				outputall();
				break;
			}
			double err;
			cudaMemcpy(&err, m_pd->LSerror(nld1), sizeof(double), cudaMemcpyDeviceToHost);
			if (err < 1e-25) {
#pragma omp single
				output("Successfully reached a local minimum! Exit.\n");
				outputall();
				break;
			}
			if (!_finite(m_d[nz / 2][ny / 2][nx / 2])) {  //blow up and exit.
#pragma omp single
				{
					sprintf(sout, "Static Director Blow up! Blow up at time %12.8f, %d steps\n", time_current, steps);
					output(sout);
				}
				break;
			}

#pragma omp single
			{
				time_current += time_step;
				steps++;
			}
			// 4. Finally update the nlterm for step2
			get_sum_fgqterm(m_pd, nld1);

			// ------------------END ------------------
			// Till here is the end, the following is for output

#pragma omp single
			time(&ltime2);
			int intv_sec = output_interval * 60;
			if ((ltimec - ltime1 < intv_sec / 50 && ltime2 - ltime1 >= intv_sec / 50) ||
				(ltimec - ltime1 < intv_sec / 10 && ltime2 - ltime1 >= intv_sec / 10) ||
				(ltimec - ltime1 < intv_sec / 2 && ltime2 - ltime1 >= intv_sec / 2) ||
				(ltimec - ltime1 < intv_sec && ltime2 - ltime1 >= intv_sec) ||
				(ltime2 - ltimec >= intv_sec)) {
				//		if ((ltime2 - ltimec >= 20)) {
#pragma omp barrier
				ltimec = ltime2;
				outputall();
			}
			bool chk_brk;
#pragma omp single copyprivate(chk_brk)
			chk_brk = (steps % 16 == 0) && check_folder();
			if (chk_brk) break;
		}

		//end parallel.
	}

	//free all memory	

	cudaFree(nld1);

	cudaFree(od);
	output("********* End Running ***********\n\n\n");
}

//FF in matlab
//nld = 1/eps*Laplace(phi^3 - (1+k)phi) - 1/eps^3 *D*(phi^3 - (1+k)phi) + 1/eps(3phi^2-D-1)(Laplace phi - 1/eps^2(phi^3-phi)) - lamb1*0.5 - lamb2 * f *3/(2*sqrt(2))
void belastic::get_nld_term(phasefield* pd, cufftDoubleComplex * nld)
{
	double vol1, sur1;
	cudaMemcpy(&vol1, pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(&sur1, pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
	double lamb1 = M1 * (vol1 - alpha);
	double lamb2 = M2 * (sur1 - beta);

	DARR3D f = pd->f();
	DARR3D nl = pd->nld();
	mem_array mem(this);
	DARR3D tmp = mem.get_array();
	calc->nld2(tmp, f, nl, lamb1, lamb2, eta);
	calc->r2c_3d_f(tmp, nld);
	//set_boundary(nld, 0.0);
}


void belastic::ETD2RK(DARR3D Het)
{
	mem_array mem(this);
	cufftDoubleComplex *a, * nld, *nnld, *tmp;

#pragma omp single copyprivate(a, nld, nnld, tmp)
	{
		a = mem.get_complex_array();
		nld = mem.get_complex_array();
		nnld = mem.get_complex_array();
		tmp = mem.get_complex_array();
	}
	m_pa->update_d();

	get_nld_term(m_pd, nld);

	cufftDoubleComplex* fd = m_pd->fd();

	calc->etd1(a, fd, nld, Het, time_step);
	m_pa->update_fd(a);

	get_nld_term(m_pa, nnld);

	fd = m_pa->fd();
	calc->etd2(tmp, fd, nnld, nld, Het, time_step);
	m_pd->update_fd(tmp);
}

void belastic::ETD4RK(DARR3D Het)
{
	mem_array mem(this);
	cufftDoubleComplex *b, *c, * nld, *nlda, * nldb, * nldc;

#pragma omp single copyprivate(a, b, c, nld, nlda, nldb, nldc, e1, e2)
	{
		b = mem.get_complex_array();
		c = mem.get_complex_array();
		nld = mem.get_complex_array();
		nlda = mem.get_complex_array();
		nldb = mem.get_complex_array();
		nldc = mem.get_complex_array();
	}
	m_pa->update_d();

	get_nld_term(m_pd, nld);
	calc->etd41(m_pa->d(), m_pd->fd(), nld, Het, time_step);
	m_pa->update_d();

	get_nld_term(m_pa, nlda);
	calc->etd42(b, m_pd->fd(), nlda, Het, time_step);

	cudaMemcpy(c, m_pa->fd(), (nx / 2 + 1) * ny * nz * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToDevice);
	m_pa->update_fd(b);

	get_nld_term(m_pa, nldb);
	calc->etd43(c, m_pd->fd(), nld, nldb, Het, time_step);
	m_pa->update_fd(c);

	get_nld_term(m_pa, nldc);
	calc->etd44(m_pd->d(), m_pd->fd(), nlda, nldb, nldc, Het, time_step);
	m_pd->update_d();
}

void belastic::run_static_fEIF()
{
	mem_array mem(this);
	//correct initial shape by trying run
	if (time_current == 0.0 && b_shape != 0) {
		shape_adjust();
	}

	DARR3D Het = mem.get_array();

	DARR3D od = mem.get_array();

	//the following are init parameters.
	double eng0, eng1, alpha0, beta0;
	int steps;
	time_t ltime1, ltime2, ltimec;

	//start parallel.
#pragma omp parallel 
	{

		//set alpha and beta
		/*
		int loop;
		for (loop = 0; loop < 100; loop++) {
			get_sum_fgqterm(m_pd, nld1);
			double a1 = calc->intf(nld1, nld1);
			double a2 = m_pd->surface();
			printf("alpha = %15.12f, beta = %15.12f\n", a1, a2);
			m_pd->update_d();
		}
	*/
		//set_boundary(m_d, -1.0);

		if (beta <= 0.0) {
			cudaMemcpy(&alpha, m_pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(&beta, m_pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
		}
		alpha0 = alpha;
		beta0 = beta;

#pragma omp single
		{
			sprintf(sout, "grid(xyz): %dx%dx%d, timestep = %e\n", nx, ny, nz, time_step);
			output(sout);
			if (fixed_time_step > 0) {
				time_step = min(fixed_time_step, time_step);
				sprintf(sout, "FIXED timestep = %e    ", time_step);
				output(sout);
			}
			else {
				output("Adjustment time_steps; ");
			}
			sprintf(sout, "eps = %9.6f, ga = %9.6f, eta = %9.6f, tmax = %9.6f\n", eps, ga, eta, tmax);
			output(sout);
			sprintf(sout, "alpha = %9.6f, beta = %9.6f, M1 = %8.2f, M2 = %8.2f, M3 = %8.2f, M4 = %8.2f, M5 = %8.2f, gweight = %8.2f\n", alpha, beta, M1, M2, M3, M4, M5, gweight);
			output(sout);

			output("********* Start STATIC Running ***********\n");
		}
		calc->init_fEIF(Het, eps);

		//eng
		eng0 = getstaticenergy(m_pd);
		eng1 = eng0;

		if (time_current == 0.0) {
#pragma omp single
			{
				sprintf(sout, "Preserving volume and surface. (Initial volume = %12.8f, surface = %12.8f\n", alpha, beta);
				output(sout);
			}
			outputall();
		}
#pragma omp single
		{
			time(&ltime1);
			ltimec = ltime1;
			steps = 0;
		}
		double c = 0;
		for (;;) {
			// ------------------START ------------------
			// 2. Update d
			//adjust time_step
			if (fixed_time_step <= 0) {
				eng0 = getstaticenergy(m_pd);
				m_pd->copy_array(od, m_pd->d());
#pragma omp single
				{
					if (steps % 4 == 0) {
						if (time_step * 1.4 < 0.001) {
							time_step = time_step * 2;
						}
						else {
							time_step = 0.001;
						}
					}
				}
				for (;;) {
					if (time_step < 1e-23) break; //too small step, get exit.

					ETD4RK(Het);

					eng1 = getstaticenergy(m_pd);

					if (eng1 > eng0) {
						m_pd->copy_array(m_pd->d(), od);
						m_pd->update_d();
#pragma omp single
						time_step /= 2;
						continue;
					}
					else break;
				}
			}
			else if (M1 != 0 || M2 != 0) {
				double a, b;
				a = alpha;
				b = beta;
				double v;
				double s;
				m_pd->copy_array(od, m_pd->d());
				ETD4RK(Het);
				cudaMemcpy(&v, m_pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(&s, m_pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
				alpha -= (v - alpha0);
				beta -= (s - beta0);
				m_pd->copy_array(m_pd->d(), od);
				m_pd->update_d();
				ETD4RK(Het);
				cudaMemcpy(&v, m_pd->vol(), sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(&s, m_pd->surface(), sizeof(double), cudaMemcpyDeviceToHost);
				alpha -= (v - alpha0);
				beta -= (s - beta0);
				//			eng1 = getstaticenergy(m_pd);
				//			eng0 = eng1;
			}
			else {
				ETD4RK(Het);
				/*
								copy_array(od, m_d);
								eng0 = getstaticenergy(m_pd);
								time_step = 1e-3;
								for (int i = 0; i < 20; i++) {
									time_step *= 2;
									ETD2RK(Het);
									eng1 = getstaticenergy(m_pd);
									output("time_step_energy", time_step, eng0 - eng1);
									copy_array(m_d, od);
									m_pd->update_d();
								}
								outputall();
								break;
				*/
				//outputall();
			}

			//exit loop.
			if (time_step < 1e-23) {
#pragma omp single
				output("time_step < 1e-23, get exit.\n");
				outputall();
				break;
			}
			/*
			if (time_current > tmax) {
	#pragma omp single
				output("Too long time running, get exit.\n");
				outputall();
				break;
			}
			*/
			if (eng0 < 1e-3) {
#pragma omp single
				output("Energy is too small < 1e-3, get exit.\n");
				outputall();
				break;
			}
			if (!_finite(m_d[nz / 2][ny / 2][nx / 2])) {  //blow up and exit.
#pragma omp single
				{
					sprintf(sout, "Static Director Blow up! Blow up at time %12.8f, %d steps\n", time_current, steps);
					output(sout);
				}
				break;
			}

#pragma omp single
			{
				time_current += time_step;
				steps++;
			}

			// ------------------END ------------------
			// Till here is the end, the following is for output

#pragma omp single
			time(&ltime2);
			int intv_sec = output_interval * 60;
			/*
					if ((ltimec - ltime1 < intv_sec/50 && ltime2 - ltime1 >= intv_sec/50) ||
						(ltimec - ltime1 < intv_sec/10 && ltime2 - ltime1 >= intv_sec/10) ||
						(ltimec - ltime1 < intv_sec/2 && ltime2 - ltime1 >= intv_sec/2) ||
						(ltimec - ltime1 < intv_sec && ltime2 - ltime1 >= intv_sec) ||
						(ltime2 - ltimec >= intv_sec)) {
			*/
					if ((ltime2 - ltimec >= 20)) {
			#pragma omp barrier
						outputall();
						time(&ltimec);
					}
			/*
			if (fabs(time_current - floor(time_current / 0.01 + 0.5) * 0.01) <= time_step / 2 || time_step >= 0.01) {
				//			printf("steps = %d, time = %d\n", steps, ltime2 - ltime1);
				outputall();

			}
			*/

			bool chk_brk;
#pragma omp single copyprivate(chk_brk)
			chk_brk = (steps % 16 == 0) && check_folder();
			if (chk_brk) break;
		}



		//end parallel.
	}

	output("********* End Running ***********\n\n\n");
}


void belastic::test1()
{
	DARR3D nld1; //non linear term -- all in spectral space.

	DARR3D od;
	/*
	ARR3D d;
	malloc_array(d);
	for (int i = 0; i < nx*ny*nz; i++) (&d[0][0][0])[i] = 1.0;
	cudaMemcpy(m_pd->d(), &d[0][0][0], nx * ny * nz * sizeof(double), cudaMemcpyHostToDevice);
	m_pd->update_d();
	*/
	printf("starting test...\n");
	cudaMalloc(&nld1, nx * ny * nz * sizeof(double));
	cudaMalloc(&od, nx * ny * nz * sizeof(double));


	double ret;
	ret = geterror(m_pd);
	//cudaMemcpy(&ret, geterror(m_pd), sizeof(double), cudaMemcpyDeviceToHost);

	printf("surface = %lf\n", ret);


	cudaFree(nld1);

	cudaFree(od);

}