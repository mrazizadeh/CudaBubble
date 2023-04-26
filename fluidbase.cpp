// fluidbase.cpp: implementation of the fluidbase class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "fluidbase.h"
#include "mem_array.h"
#include <cuda_runtime.h>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

fluidbase::fluidbase(int findex)
{
	sout = new char[1000];
	folderindex = findex;
	if (folderindex < 0) folderindex = 0;

	//the following are init parameters.
	int epssize = 2; //1 is small, 2 is large

	int ifold = 101;
	char inifile[100], inibakfile[100];

	einit = 0; m_k = 1.0; m_c = 0.0; m_alpha3 = 0.0;

	sprintf(inifile, "%sini%d.dat", RESULT_PATH, folderindex);
	FILE * f = fopen(inifile, "r");
	if (!f) {
		printf("Can not open the ini file [%s]\n", inifile);
		exit(0);
	}
	fscanf(f, "%d %d %d %d %lf %lf %lf %d %lf %lf %lf %d %d %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %d %lf %lf %lf %lf", 
		       &ifold, &nx, &ny, &nz, &xlen, &ylen, &zlen, &b_shape, &time_step, &fluid_time_step, &fixed_time_step, &ninit, &nflush, &surctrl, &epssize, &ncur, 
			   &m_angle, &re, &eta, &ga, &alpha, &beta, &M1, &M2, &M3, &M4, &tmax, &output_interval, &gweight, &einit, &m_k, &m_c, &m_linet, &m_alpha3);
	fclose(f);

	//ny, nz,  nx;  may changed in other classes' construction function
	if (nx <= 0 || ny <= 0 || nz <= 0)  {
		printf("Illegal ngrid number!\n");
		exit(0);
	}

	sprintf(folderpath, "%s%d/", RESULT_PATH, ifold);
	output(folderpath);
	output("\n");
	sprintf(inibakfile, "%sini%d.dat", folderpath, folderindex);
	copyfile(inifile, inibakfile);

	openmp_info();
	
	if (xlen <= 0) xlen = 1;
	if (ylen <= 0) ylen = 1;
	if (zlen <= 0) zlen = 1;
	xlen *= 2.0*PI;
	ylen *= 2.0*PI;
	zlen *= 2.0*PI;

/*	0: middle 2.0h
	1: small, 1.77h
	2: larger, 2.5h
	3: smaller, 1.5h
	4: more smaller, 1.25h
	5: much smaller, 1h,
	6: smallest, 0.7h,
	7: most smallest, 0.5h
	other: middle 2.0h*/
	if (epssize == 1) {
		eps = 2.5*(sqrt(2.)*xlen/2.0/nx);
		output("Smaller eps = 1.768 h\n");
	}
	else if (epssize == 2) {
		eps = 2.5*(2*xlen/2.0/nx);
		output("Larger eps = 2.5h\n");
	}
	else if (epssize == 3) {
		eps = 1.5*(2*xlen/2.0/nx);
		output("smaller eps = 1.5h\n");
	}
	else if (epssize == 4) {
		eps = 1.25*(2*xlen/2.0/nx);
		output("more smaller eps = 1.25h\n");
	}
	else if (epssize == 5) {
		eps = 1.0*(2*xlen/2.0/nx);
		output("much smaller eps = 1h\n");
	}
	else if (epssize == 6) {
		eps = 0.7*(2*xlen/2.0/nx);
		output("smallest eps = 0.7h\n");
	}
	else if (epssize == 7) {
		eps = 4.0*(2*xlen/2.0/nx);
		output("most smallest eps = 4.0h\n");
	}
	else if (epssize == 8) {
		eps = 2.5*(sqrt(2.)*xlen/2.0/nx)*1.22;
		output("most smallest eps = 2.16h\n");
	}
	else if (epssize > 10) {
		eps = epssize/100.0 * (2*xlen/2.0/nx);
		sprintf(sout, "eps = %4.2fh\n", epssize/100.0);
		output(sout);
	}
	else {
		eps = 2.0*(2*xlen/2.0/nx);
		output("Middle eps = 2h\n");
	}

	xi = eps;

	if (time_step <= 0) time_step = 1e-6;
	if (time_start < 0) time_start = 0;
	time_current = time_start;

	if (ninit < -1) ninit = 0;
	if (einit < -1) einit = 0;
	if (nflush < -1) nflush = 0;

	if (surctrl > 0.0) {
		sprintf(sout, "All Control at %04.2f; ", surctrl);
		output(sout);
	}
	else {
		output("No Surface Control; ");
	}

	if (ncur < 0) ncur = 0;

	m_angle = m_angle*PI/180.0;

	if (re < 0) re = 1.0;

	//eta: elastic bending rigidity. For two component vesicles, it is set for elastic relaxation.
	if (eta < 0) eta = 10.0;
	if (ga <= 0) ga = 1.0;
	if (tmax < 0) tmax = 4.0;
	if (output_interval <= 0) output_interval = 30;
//	if (gweight <= 0) gweight = 0;

	if (m_linet < 0) m_linet = 1.0;

	if (beta < 0) alpha = beta = 0;

	if (m_k < 0) m_k = 1.0;
	if (fabs(m_alpha3) > beta) m_alpha3 = 0;

	M5 = 0.0; //xwang.

	C = NULL;

	xx = new double[nx]; 
	yy = new double[ny]; 
	zz = new double[nz]; 
	int i;
	for (i = 0; i < nx; i++) 
		xx[i] = i * xlen / nx - xlen / 2.0;
	for (i = 0; i < ny; i++)
		yy[i] = i * ylen / ny - ylen / 2.0;
	for (i = 0; i < nz; i++)
		zz[i] = i * zlen / nz - zlen / 2.0;

	time_start = time_current = 0;

	load_info();   //try it. 
	
	tmp_array_len = tmp_complex_array_len = 0;

	if (nflush != -1) {
		malloc_array(m_u);
		malloc_array(m_v);
		malloc_array(m_w);
	}
	if (ninit != -1) malloc_array(m_d); 
	if (einit != -1) malloc_array(m_e); 
}

fluidbase::~fluidbase()
{
	for (int i = 0; i < tmp_array_len; i++)
		cudaFree(base_tmp_array[i]);
	for (int i = 0; i < tmp_complex_array_len; i++)
		cudaFree(base_tmp_complex_array[i]);

	if (C) free_array(C);
	free(zz); free(yy); free(xx);
	if (einit != -1) free_array(m_e);
	if (ninit != -1) free_array(m_d);
	if (nflush != -1) {
		free_array(m_w);
		free_array(m_v);
		free_array(m_u);
	}
	delete[] sout;
}

bool fluidbase::malloc_array(ARR3D  &x)
{
	return ::malloc_array(x, nx, ny, nz);
}



void fluidbase::output(const char *s)
{
	char fname[100];
	sprintf(fname, "%soutput%dx%dx%d.txt", folderpath, nx, ny, nz);
	FILE * f = fopen(fname, "a+");
	if (!f) {
		printf("Can not write file %s!\n", fname);
		return;
	}
#ifdef _OUTPUT_TIME_
    time_t ltime;
    time( &ltime );
	fprintf(f, "%s %s", ctime( &ltime ), s);
#else
	fprintf(f, "%s", s);
#endif
	printf("%s", s);
	fclose(f);
}

void fluidbase::output(char *name, ARR3D u)
{
	int i, j;
	char fname[100];
	sprintf(fname, "%s%s%dx%dx%d_%08d_xz.txt", folderpath, name, nx, ny, nz, roundx(time_current/TIME_STEP_OUTPUT));
	FILE * f = fopen(fname, "w+");
	if (!f) {
		printf("Can not write file %s!\n", fname);
		return;
	}
	//x-z plan
	for(i = 0; i < nz; i++){
		for(j = 0; j <nx; j++)
		{
			fprintf(f, "%7.5f ", u[i][ny/2][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);

	sprintf(fname, "%s%s%dx%dx%d_%08d_xy.txt", folderpath, name, nx, ny, nz, roundx(time_current/TIME_STEP_OUTPUT));
	f = fopen(fname, "w+");
	if (!f) {
		printf("Can not write file %s!\n", fname);
		return;
	}
	//x-y plan
	for(i = 0; i < ny; i++){
		for(j = 0; j <nx; j++)
		{
			fprintf(f, "%6.3f ", u[nz/2][i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);

	if (nx*ny*nz*8 > 100000) 
		sprintf(fname, "%s%s%dx%dx%d.out", folderpath, name, nx, ny, nz);
	else 
		sprintf(fname, "%s%s%dx%dx%d_%08d.out", folderpath, name, nx, ny, nz, roundx(time_current/TIME_STEP_OUTPUT));
	write_array(u, fname);
}

void fluidbase::input(char *name, ARR3D u)
{
	char fname[100];
	if (nx*ny*nz*8 > 100000) 
		sprintf(fname, "%s%s%dx%dx%d.out", folderpath, name, nx, ny, nz);
	else 
		sprintf(fname, "%s%s%dx%dx%d_%08d.out", folderpath, name, nx, ny, nz, roundx(time_current/TIME_STEP_OUTPUT));
	read_array(u, fname);
}

/*
void fluidbase::input(char *name, ARR3D u)
{
	char fname[100];
	sprintf(fname, "%s%s%dx%dx%d.out", folderpath, name, 65, 65, 64);
	read_array_intepolate(u, fname, 65, 65, 64, 45);
}
*/

void fluidbase::output(char *name, double x)
{
	char fname[100];
	sprintf(fname, "%s%s%dx%dx%d.txt", folderpath, name, nx, ny, nz);
	FILE * f = fopen(fname, "a+");
	if (!f) {
		printf("Can not write file %s!\n", fname);
		return;
	}
	fprintf(f, "%08d, %15.12f\n", roundx(time_current/TIME_STEP_OUTPUT), x);
	fclose(f);
}

void fluidbase::output(char *name, double x1, double x2)
{
	char fname[100];
	sprintf(fname, "%s%s%dx%dx%d.txt", folderpath, name, nx, ny, nz);
	FILE * f = fopen(fname, "a+");
	if (!f) {
		printf("Can not write file %s!\n", fname);
		return;
	}
	fprintf(f, "%08d, %15.12f, %15.12f\n", roundx(time_current/TIME_STEP_OUTPUT), x1, x2);
	fclose(f);
}

void fluidbase::output(char *name, double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9)
{
	char fname[100];
	sprintf(fname, "%s%s%dx%dx%d.txt", folderpath, name, nx, ny, nz);
	FILE * f = fopen(fname, "a+");
	if (!f) {
		printf("Can not write file %s!\n", fname);
		return;
	}
	fprintf(f, "%08d, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f, %8.5f\n", roundx(time_current/TIME_STEP_OUTPUT), x1, x2, x3, x4, x5, x6, x7, x8, x9);
	fclose(f);
}

void fluidbase::release_tmp_array(int i)
{
	tmp_array_used[i] = false;
}

void fluidbase::release_tmp_complex_array(int i)
{
	tmp_complex_array_used[i] = false;
}

DARR3D  fluidbase::get_tmp_array(int& ind)
{
	int i;
	for (i = 0; i < tmp_array_len; i++) {
		if (!tmp_array_used[i]) {
			tmp_array_used[i] = true;
			ind = i;
			return base_tmp_array[i];
		}
	}
	if (tmp_array_len == TMP_ARRAY_MAX) {
		output("Sorry, base_tmp_array all used!\n");
		return NULL;
	}
	else {
		tmp_array_len++;
		ind = tmp_array_len - 1;
		tmp_array_used[ind] = true;
		cudaMalloc(&base_tmp_array[ind], nx * ny * nz * sizeof(double));
		return base_tmp_array[ind];
	}
}

cufftDoubleComplex *  fluidbase::get_tmp_complex_array(int& ind)
{
	int i;
	for (i = 0; i < tmp_complex_array_len; i++) {
		if (!tmp_complex_array_used[i]) {
			tmp_complex_array_used[i] = true;
			ind = i;
			return base_tmp_complex_array[i];
		}
	}
	if (tmp_complex_array_len == TMP_ARRAY_MAX) {
		output("Sorry, base_tmp_array all used!\n");
		return NULL;
	}
	else {
		tmp_complex_array_len++;
		ind = tmp_complex_array_len - 1;
		tmp_complex_array_used[ind] = true;
		cudaMalloc(&base_tmp_complex_array[ind], (nx/2+1) * ny * nz * sizeof(cufftDoubleComplex));
		return base_tmp_complex_array[ind];
	}
}

//init for memorym, shape, display angle, fluid, spontaneous curvature.
//must be done after set_mesh.
void fluidbase::init()
{
	if (time_start != 0.0) {
		reload();
		output("Initial shape and fluid velocity from files Loaded.\n");
	}
	else {
		if (nflush == 0) {
			init_fluid_zero();
		}
		else if (nflush == 1) {
			init_fluid_flush();
		}
		else if (nflush == 2) {
			init_fluid_flushxz();
		}
		else if (nflush == 3) {
			init_fluid_flushp();
		}
		else if (nflush == 4) {
			init_fluid_flushp_slow();
		}
		else if (nflush == 5) {
			init_fluid_flush_slow();
		}
		else if (nflush == 6) {
			init_fluid_flush_cos();
		}
		else if (nflush == 7) {
			input("u", m_u);
			input("v", m_v);
			input("w", m_w);
			output("Loading initial fluid from file, multiply by 3; ");
			int i, j, k;
			for (i = 0; i < nz; i++) {
				for (j = 0; j < ny; j++) {
					for (k = 0; k < nx; k++) {
						m_u[i][j][k] *= 3;
						m_v[i][j][k] *= 3;
						m_w[i][j][k] *= 3;
					}
				}
			}
		}
		else if (nflush == 8) {
			init_fluid_flush_middle();
		}
		else if (nflush == 9) {
			init_fluid_flushxz2();
		}
		else if (nflush == 10) {
			init_fluid_curve();
		}
		else if (nflush == 11) {
			init_fluid_flush_cos_fast();
		}
		else if (nflush == 12) {
			init_fluid_tanh();
		}
		else if (nflush == 13) {
			init_fluid_flushp_middle();
		}

		if (einit == 0) {
			init_eta1();
		}
		else if (einit == 1) {
			init_eta2();
		}
		else if (einit == 2) {
			init_eta3();
		}
		else if (einit == 3) {
			init_eta4();
		}
		else if (einit == 4) {
			init_eta5();
		}
		else if (einit == 5) {
			init_eta6();
		}
		else if (einit == 7) {
			input("e", m_e);
			output("Loading initial eta phase from file; ");
		}
		else if (einit == 8) {
			init_eta8();
		}
		else if (einit == 9) {
			init_eta9();
		}
		else if (einit == 10) {
			init_eta10();
		}
		else if (einit == 11) {
			init_eta11();
		}
		else if (einit == 12) {
			init_eta12();
		}
		else if (einit == 13) {
			init_eta13();
		}
		else if (einit == 14) {
			init_etaleftsphere();
		}
		else if (einit == 35) {
//			init_disk();
//			copy_array(m_e, m_d);
//			shiftright(m_e, m_d, 7);
//			shiftright(m_d, m_e, -14);
//			shift(m_e, 8);
//			shift(m_d, -8);
			init_disk();
			shiftright(m_d, m_e, 7);
			copy_array(m_d, m_e);
			shift(m_e, 14);
			shift(m_d, -14);
			int i, j, k;
			for (i = 0; i < nz; i++) {
				for (j = 0; j < ny; j++) {
					for (k = 0; k < nx; k++) {
						m_e[i][j][k] += m_d[i][j][k] + 1.0;
					}
				}
			}
//			init_disk();
//			copy_array(m_u, m_d);
//			shiftright(m_u, m_d, -7);
			shiftright(m_e, m_d, -14);
			shift(m_d, -14);

		}
		else if (einit == 36) {
	double h = zlen/nz;
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i] - h/2;
		for (j = 0; j < ny; j++) {
			double y  = yy[j] - h/2;
			for (k = 0; k < nx; k++) {
				double x = xx[k] - h/2;
				double r;
				r = sqrt(x*x + y*y + z*z);
				m_e[i][j][k] = tanh((0.3*ylen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
			shift(m_e, 15);
			shift(m_d, -15);
		}
		else {
			init_eta1();
		}


		if (ninit == 0) {
			init_1sphere();
		}
		else if (ninit == 1) {
			init_2spheres();
		}
		else if (ninit == 2) {
			init_ellipse();
		}
		else if (ninit == 3) {
			init_triflower();
		}
		else if (ninit == 4) {
			init_gtorus();
		}
		else if (ninit == 5) {
			init_ellipse2();
		}
		else if (ninit == 6) {
			init_ellipse3();
		}
		else if (ninit == 7) {
			input("d", m_d);
			output("Loading initial shape from file; ");
			if (m_angle < 0) m_angle = PI/4;

//			changeshape(m_d, m_w);
//			copy_array(m_d, m_w);
		}
		else if (ninit == 8) {
			init_vtorus();
		}
		else if (ninit == 9) {
			init_htorus();
		}
		else if (ninit == 10) {
			init_1sphere_right();
		}
		else if (ninit == 11) {
			init_ellipse4();
		}
		else if (ninit == 12) {
			init_gtorus2();
		}
		else if (ninit == 13) {
			init_2ellipse();
		}
		else if (ninit == 14) {
			init_3spheres();
		}
		else if (ninit == 15) {
			init_hulu();
		}
		else if (ninit == 16) {
			init_cherry();
		}
		else if (ninit == 17) {
			init_3spheres_touch();
		}
		else if (ninit == 18) {
			init_3spheres_curv();
		}
		else if (ninit == 19) {
			init_3spheres_curv_touch();
		}
		else if (ninit == 20) {
			init_3spheres_vertical();
		}
		else if (ninit == 21) {
			init_ciga_vertical();
		}
		else if (ninit == 22) {
			init_ciga_horizontal();
		}
		else if (ninit == 23) {
			init_torus_small_ball();
		}
		else if (ninit == 24) {
			init_torus_big_ball();
		}
		else if (ninit == 25) {
			init_torus_bing();
		}
		else if (ninit == 26) {
			init_two_torus_ring();
		}
		else if (ninit == 27) {
			init_cherry_inverse();
		}
		else if (ninit == 28) {
			init_art_cherry();
		}
		else if (ninit == 29) {
			init_1sphere_small_right();
		}
		else if (ninit == 30) {
			init_2sphere_small();
		}
		else if (ninit == 31) {
			init_ciga_left();
		}
		else if (ninit == 32) {
			init_ciga_two();
		}
		else if (ninit == 33) {
			init_ellipse5();
		}
		else if (ninit == 34) {
			init_ciga_long_thin();
		}
		else if (ninit == 35) {
			init_disk();
		}
		else if (ninit == 36) {
			init_8spheres();
		}
		else if (ninit == 37) {
			init_ellipse6();
		}
		else if (ninit == 38) {
			init_ellipse7();
		}
		else if (ninit == 39) {
			init_ellipse8();
		}
		else if (ninit == 40) {
			init_25spheres();
		}
		else if (ninit == 41) {
			init_26spheres();
		}
		else if (ninit == 42) {
			init_5spheres();
		}
		else if (ninit == 43) {
			init_drop();
		}
		else if (ninit == 44) {
			init_cube();
		}
		else if (ninit == 45) {
			init_ellipse_left();
		}
		else if (ninit == 46) {
			init_rightsphere();
		}
		else if (ninit == 47) {
			init_4spheres_touch();
		}
		else if (ninit == 48) {
			init_9_cylinder();
		}
		else if (ninit == 49) {
			init_4_bings();
		}
		else if (ninit == 50) {
			init_tian();
		}
		else if (ninit == 51) {
			init_Gyroid();
		}
		else if (ninit == 52) {
			init_3cylinder();
		}
	}

	if (ncur == 0) {
		output("No spontenuous curvature.\n");
		C = NULL;
	}
	else if (ncur == 1) {
		init_curvature1();
	}
	else if (ncur == 2) {
		init_curvature2();
	}
	else if (ncur == 3) {
		init_curvature3();
	}
}

void fluidbase::init_eta1()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = fabs(z - ylen/2*0.4);
				m_e[i][j][k] = tanh((0.4*ylen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("eta phase field half by half the domain.\n");
}

void fluidbase::init_eta2()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = fabs(z);
				//m_e[i][j][k] = tanh((0.15*ylen/2.0-r)/(sqrt(2.0)*eps));
				m_e[i][j][k] = tanh((0.375*ylen/2.0-r)/(sqrt(2.0)*eps));

//				if (fabs(x) > ylen/2.0*0.8 || fabs(y) > ylen/2.0*0.8) m_e[i][j][k] = -1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("eta phase field small middle the domain.\n");
}

void fluidbase::init_eta3()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = fabs(z - ylen/2*0.6);
				m_e[i][j][k] = tanh((0.25*ylen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("eta phase field small half by large half the domain.\n");
}

void fluidbase::init_eta4()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = fabs(z - ylen/2*0.6);
				m_e[i][j][k] = -tanh((0.25*ylen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("eta phase field small half by large half the domain -- inverse.\n");
}

void fluidbase::init_eta5()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double a, b, c, t;

				a = 1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d1 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = -1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d2 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = -1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d3 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = 1; c = -1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d4 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				if (d1 < 0.2*ylen/2.0 || d2 < 0.2*ylen/2.0 || d3 < 0.2*ylen/2.0 || d4 < 0.2*ylen/2.0) m_e[i][j][k] = 1.0;
				else m_e[i][j][k] = -1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("eta phase field Eight corner bumps.\n");
}

void fluidbase::init_eta6()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double a, b, c, t;

				a = 1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d1 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = -1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d2 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = -1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d3 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = 1; c = -1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d4 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				double d5 = sqrt(x*x + y*y);
				double d6 = sqrt(x*x + z*z);
				double d7 = sqrt(y*y + z*z);

				if (d1 < 0.2*ylen/2.0 || d2 < 0.2*ylen/2.0 || d3 < 0.2*ylen/2.0 || d4 < 0.2*ylen/2.0 || 
					d5 < 0.2*ylen/2.0 || d6 < 0.2*ylen/2.0 || d7 < 0.2*ylen/2.0) m_e[i][j][k] = 1.0;
				else m_e[i][j][k] = -1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("eta phase field 14 corner (eight corner plus six sides) bumps.\n");
}

/*
void fluidbase::init_eta8()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x+y*y +z*z);
				if (r > 2.3 || r < 1.8) m_e[i][j][k] = 1.0;
				else m_e[i][j][k] = -1.0;
			}
		}
	}
	output("Construct eta file; for the best one.");
}
*/

void fluidbase::init_eta8()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x+y*y +z*z);
				m_e[i][j][k] = tanh((1.8 - r)/(sqrt(2.0)*eps)) + tanh((r - 2.3)/(sqrt(2.0)*eps)) + 1.0;
			}
		}
	}
	output("Construct eta file; for the best one.");
}

void fluidbase::init_eta9()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double a, b, c, t;

				a = 1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d1 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = -1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d2 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = -1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d3 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = 1; c = -1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d4 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				double d5 = sqrt(x*x + y*y);
				double d6 = sqrt(x*x + z*z);
				double d7 = sqrt(y*y + z*z);

				if ((d1 < 0.2*ylen/2.0 && d1 > 0.1*ylen/2.0) || (d2 < 0.2*ylen/2.0 && d2 > 0.1*ylen/2.0)  || 
					(d3 < 0.2*ylen/2.0 && d3 > 0.1*ylen/2.0)  || (d4 < 0.2*ylen/2.0 && d4 > 0.1*ylen/2.0)  || 
					(d5 < 0.2*ylen/2.0 && d5 > 0.1*ylen/2.0)  || (d6 < 0.2*ylen/2.0 && d6 > 0.1*ylen/2.0)  || 
					(d7 < 0.2*ylen/2.0 && d7 > 0.1*ylen/2.0) ) m_e[i][j][k] = -1.0;
				else m_e[i][j][k] = 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("eta phase field 14 corner (eight corner plus six sides) colorful bumps.\n");
}

void fluidbase::init_eta10()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double a, b, c, t;

				a = 1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d1 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = -1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d2 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = -1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d3 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = 1; c = -1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d4 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				double d5 = sqrt(x*x + y*y);
				double d6 = sqrt(x*x + z*z);
				double d7 = sqrt(y*y + z*z);

				if ((d1 < 0.2*ylen/2.0 && d1 > 0.13*ylen/2.0) || (d2 < 0.2*ylen/2.0 && d2 > 0.13*ylen/2.0)  || 
					(d3 < 0.2*ylen/2.0 && d3 > 0.13*ylen/2.0)  || (d4 < 0.2*ylen/2.0 && d4 > 0.13*ylen/2.0)  || 
					(d5 < 0.2*ylen/2.0 && d5 > 0.13*ylen/2.0)  || (d6 < 0.2*ylen/2.0 && d6 > 0.13*ylen/2.0)  || 
					(d7 < 0.2*ylen/2.0 && d7 > 0.13*ylen/2.0) ) m_e[i][j][k] = -1.0;
				else m_e[i][j][k] = 1.0;
			}
		}
	}

	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x+y*y +z*z);
				if (r < 1.0) m_e[i][j][k] = 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Construct eta file; for the new2.");
}

void fluidbase::init_eta11()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x+y*y +z*z);
				if (r > 2.175) m_e[i][j][k] = -1.0;
				else m_e[i][j][k] = 1.0;
			}
		}
	}
	output("Construct eta file; for the new3.");
}

void fluidbase::init_eta12()
{
	double R = 2.04;
//	double R = 2.45;
	double r = 0.6283;
//	double r = 0.45;
	double s = R/sqrt(3.0);
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1 = sqrt((x-R)*(x-R)+y*y +z*z);
				double r2 = sqrt((x+R)*(x+R)+y*y +z*z);
				double r3 = sqrt(x*x+(y-R)*(y-R) +z*z);
				double r4 = sqrt(x*x+(y+R)*(y+R) +z*z);
				double r5 = sqrt(x*x+y*y + (z-R)*(z-R));
				double r6 = sqrt(x*x+y*y + (z+R)*(z+R));
				double r7 = sqrt((x-s)*(x-s) + (y-s)*(y-s) + (z-s)*(z-s));
				double r8 = sqrt((x+s)*(x+s) + (y-s)*(y-s) + (z-s)*(z-s));
				double r9 = sqrt((x-s)*(x-s) + (y+s)*(y+s) + (z-s)*(z-s));
				double r10 = sqrt((x+s)*(x+s) + (y+s)*(y+s) + (z-s)*(z-s));
				double r11 = sqrt((x-s)*(x-s) + (y-s)*(y-s) + (z+s)*(z+s));
				double r12 = sqrt((x+s)*(x+s) + (y-s)*(y-s) + (z+s)*(z+s));
				double r13 = sqrt((x-s)*(x-s) + (y+s)*(y+s) + (z+s)*(z+s));
				double r14 = sqrt((x+s)*(x+s) + (y+s)*(y+s) + (z+s)*(z+s));
				if (r1 < r || r2 < r || r3 < r || r4 < r || r5 < r || r6 < r || r7 < r || r8 < r
					 || r9 < r || r10 < r || r11 < r || r12 < r || r13 < r || r14 < r) 
					 m_e[i][j][k] = -1.0;
				else m_e[i][j][k] = 1.0;
			}
		}
	}
	output("Construct eta file; for the new4.");
}

void fluidbase::init_eta13() 
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((z-0.5*ylen/2.0)*(z-0.5*ylen/2.0)+y*y +x*x);
				m_e[i][j][k] = -tanh((0.1*ylen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	output("Construct eta file; a small ball to cut the big surface.");
}

void fluidbase::init_1sphere()
{
	double h = zlen/nz;
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r;
				r = sqrt(x*x + y*y + z*z);
				m_d[i][j][k] = tanh((0.5 * ylen / 2.0 - r) / (sqrt(2.0) * eps));
				//m_d[i][j][k] = sin(y);
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a sphere.\n");
}

void fluidbase::init_1sphere_right()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r;
				r = sqrt((x+0.5*xlen/2.0)*(x+0.5*xlen/2.0) + y*y + z*z);
				m_d[i][j][k] = tanh((0.3*ylen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a small sphere at left.\n");
}

void fluidbase::init_rightsphere()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r;
				r = sqrt((x+0.3*xlen/2.0)*(x+0.3*xlen/2.0) + y*y + z*z);
				m_d[i][j][k] = tanh((0.3*xlen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a small sphere at left.\n");
}

void fluidbase::init_etaleftsphere()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r;
				r = sqrt((x-0.3*xlen/2.0)*(x-0.3*xlen/2.0) + y*y + z*z);
				m_e[i][j][k] = tanh((0.3*xlen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a small sphere at right.\n");
}

void fluidbase::init_1sphere_small_right()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r;
				r = sqrt((x+0.5*xlen/2.0)*(x+0.5*xlen/2.0) + y*y + z*z);
				m_d[i][j][k] = tanh((0.25*ylen/2.0-r)/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a very small sphere at left.\n");
}

void fluidbase::init_2sphere_small()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r, r2;
				r = 0.28*zlen/2.0 - sqrt((z+0.35*xlen/2.0)*(z+0.35*xlen/2.0) + y*y + x*x);
				r2 = 0.28*zlen/2.0 - sqrt((z-0.35*xlen/2.0)*(z-0.35*xlen/2.0) + y*y + x*x);
				m_d[i][j][k] = tanh((r)/(sqrt(2.0)*eps)) + tanh((r2)/(sqrt(2.0)*eps)) + 1;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as two very small sphere at two sides.\n");
}

void fluidbase::init_2spheres()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2;
				r1 = sqrt(x*x + y*y + (z-0.4*zlen/2.0)*(z-0.4*zlen/2.0));
				r2 = sqrt(x*x + y*y + (z+0.5*zlen/2.0)*(z+0.5*zlen/2.0));
				m_d[i][j][k] = tanh((0.3*zlen/2.0-r1)/(sqrt(2.0)*eps)) + tanh((0.2*zlen/2.0-r2)/(sqrt(2.0)*eps)) + 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a larger sphere and a smaller sphere.\n");
}

void fluidbase::init_3spheres()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2, r3;
				r1 = sqrt(x*x + (y-0.5*ylen/2.0)*(y-0.5*ylen/2.0) + z*z);
				r2 = sqrt((x-0.25*sqrt(3.0)*xlen/2.0)*(x-0.25*sqrt(3.0)*xlen/2.0) + (y+0.25*ylen/2.0)*(y+0.25*ylen/2.0) + z*z);
				r3 = sqrt((x+0.25*sqrt(3.0)*xlen/2.0)*(x+0.25*sqrt(3.0)*xlen/2.0) + (y+0.25*ylen/2.0)*(y+0.25*ylen/2.0) + z*z);
				m_d[i][j][k] = tanh((0.3*ylen/2.0-r1)/(sqrt(2.0)*eps)) 
					+ tanh((0.3*ylen/2.0-r2)/(sqrt(2.0)*eps)) + tanh((0.3*ylen/2.0-r3)/(sqrt(2.0)*eps)) + 2.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 3 spheres with same size and on one plan.\n");
}

void fluidbase::init_3spheres_vertical()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2, r3;
				r1 = sqrt(x*x + (z+0.5*ylen/2.0)*(z+0.5*ylen/2.0) + y*y);
				r2 = sqrt((x-0.27*sqrt(3.0)*xlen/2.0)*(x-0.27*sqrt(3.0)*xlen/2.0) + (z-0.25*ylen/2.0)*(z-0.25*ylen/2.0) + y*y);
				r3 = sqrt((x+0.25*sqrt(3.0)*xlen/2.0)*(x+0.25*sqrt(3.0)*xlen/2.0) + (z-0.25*ylen/2.0)*(z-0.25*ylen/2.0) + y*y);
				m_d[i][j][k] = tanh((0.3*ylen/2.0-r1)/(sqrt(2.0)*eps)) 
					+ tanh((0.3*ylen/2.0-r2)/(sqrt(2.0)*eps)) + tanh((0.3*ylen/2.0-r3)/(sqrt(2.0)*eps)) + 2.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 3 spheres with same size and on one vertical plan.\n");
}

void fluidbase::init_ciga_vertical()
{
	input("../ciga", m_d);
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a vertical ciga, need middle size of eps = 2h.\n");
}

void fluidbase::init_ciga_long_thin()
{
	input("../ciga", m_d);
	ARR3D tmp;
	malloc_array(tmp);
	copy_array(tmp, m_d);

	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				m_d[i][j][k] = -1;
			}
		}
	}
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				int a, b;
				a = j - ny/2;
				b = k - nx/2;
				m_d[i][a*2/3+ny/2][b*2/3+nx/2] = tmp[i][j][k];
			}
		}
	}
	free_array(tmp);
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a vertical long and thin ciga, need middle size of eps = 2h.\n");
}

void fluidbase::init_ciga_horizontal()
{
	input("../ciga", m_d);
	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = i; k < nx-1; k++) {
				swap(m_d[i][j][k], m_d[k][j][i]);
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a horizontal ciga, need middle size of eps = 2h.\n");
}

void fluidbase::init_ciga_left()
{
	input("../ciga", m_d);
	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = nx/4; k < nx; k++) {
				m_d[i][j][(k+nx*3/4)%nx] = m_d[i][j][k];
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a left ciga, need middle size of eps = 2h.\n");
}

void fluidbase::init_ciga_two()
{
	input("../ciga", m_d);
	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = nx/4; k < nx; k++) {
				m_d[i][j][(k+nx*3/4)%nx] = m_d[i][j][k];
			}
		}
	}

	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = nx/2; k < nx; k++) {
				m_d[i][j][k] = m_d[i][j][nx-1-k];
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as two cigas, need middle size of eps = 2h.\n");
}

void fluidbase::init_torus_small_ball()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + y*y);
				if (r == 0) {
					r = 1;
				}
				double theta = acos(x/r);
				if (y < 0) theta = 2*PI - theta;
				double d = 0.25*ylen/2.0 - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (y - 0.5*ylen/2.0*sin(theta))*(y - 0.5*ylen/2.0*sin(theta)) + (z-0.25*zlen/2.0)*(z-0.25*zlen/2.0));
				double d2 = 0.15*zlen/2 - sqrt(x*x + y*y + (z+0.2*zlen/2)*(z+0.2*zlen/2));
				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps)) + tanh(d2/(sqrt(2.0)*eps)) + 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a torus with a small ball.\n");
}


void fluidbase::init_torus_big_ball()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + y*y);
				if (r == 0) {
					r = 1;
				}
				double theta = acos(x/r);
				if (y < 0) theta = 2*PI - theta;
				double d = 0.25*ylen/2.0 - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (y - 0.5*ylen/2.0*sin(theta))*(y - 0.5*ylen/2.0*sin(theta)) + (z-0.3*zlen/2.0)*(z-0.3*zlen/2.0));
				double d2 = 0.35*zlen/2 - sqrt(x*x + y*y + (z+0.4*zlen/2)*(z+0.4*zlen/2));
				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps)) + tanh(d2/(sqrt(2.0)*eps)) + 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a torus with a big ball.\n");
}

void fluidbase::init_torus_bing()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + y*y);
				if (r == 0) {
					r = 1;
				}
				double theta = acos(x/r);
				if (y < 0) theta = 2*PI - theta;
				double d = 0.25*ylen/2.0 - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (y - 0.5*ylen/2.0*sin(theta))*(y - 0.5*ylen/2.0*sin(theta)) + (z-0.3*zlen/2.0)*(z-0.3*zlen/2.0));
				double d2 = sqrt((x*x + y*y)/9 + (z+0.4*zlen/2)*(z+0.4*zlen/2))*3;
				d2 = (0.7*ylen/2.0-d2)/0.5;
				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps)) + tanh(d2/(sqrt(2.0)*eps)) + 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a torus with a pancake.\n");
}

void fluidbase::init_two_torus_ring()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k] + 0.25*xlen/2.0;
				double r = sqrt(x*x + y*y);
				if (r == 0) {
					r = 1;
				}
				double theta = acos(x/r);
				if (y < 0) theta = 2*PI - theta;
				double d = 0.125*ylen/2.0 - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (y - 0.5*ylen/2.0*sin(theta))*(y - 0.5*ylen/2.0*sin(theta)) + z*z);

				x = xx[k] - 0.25*xlen/2.0;
				r = sqrt(x*x + z*z);
				if (r == 0) {
					r = 1;
				}
				theta = acos(x/r);
				if (z < 0) theta = 2*PI - theta;
				double d2 = 0.125*ylen/2.0 - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (z - 0.5*ylen/2.0*sin(theta))*(z - 0.5*ylen/2.0*sin(theta)) + y*y);

				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps)) + tanh(d2/(sqrt(2.0)*eps)) + 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as two torus ring, eps should be small.\n");
}

void fluidbase::init_3spheres_touch()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2, r3;
				r1 = sqrt((x-0.3*xlen/2)*(x-0.3*xlen/2) + y*y + (z-0.4*xlen/2)*(z-0.4*xlen/2));
				r2 = sqrt((x+0.3*xlen/2)*(x+0.3*xlen/2) + y*y + z*z);
				r3 = sqrt((x-0.3*xlen/2)*(x-0.3*xlen/2) + y*y + (z+0.4*xlen/2)*(z+0.4*xlen/2));
				m_d[i][j][k] = tanh((0.4*xlen/2.0-r1)/(sqrt(2.0)*eps)) 
					+ tanh((0.4*xlen/2.0-r2)/(sqrt(2.0)*eps)) + tanh((0.4*xlen/2.0-r3)/(sqrt(2.0)*eps)) + 2.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 3 spheres with same size touched together.\n");
}

void fluidbase::init_4spheres_touch()
{
	int i, j, k;
	double a = 0.4 / (1.5*sqrt(2.0)) * xlen/2;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2, r3, r4;
				r1 = sqrt((x-(-2)*a)*(x-(-2)*a) + (y - (-sqrt(2.)/2)*a)*(y - (-sqrt(2.)/2)*a) + z*z);
				r2 = sqrt((x-(1)*a)*(x-(1)*a) + (y - (-sqrt(2.)/2)*a)*(y - (-sqrt(2.)/2)*a) + (z - (sqrt(3.)*a))*(z - (sqrt(3.)*a)));
				r3 = sqrt((x-(1)*a)*(x-(1)*a) + (y - (-sqrt(2.)/2)*a)*(y - (-sqrt(2.)/2)*a) + (z - (-sqrt(3.)*a))*(z - (-sqrt(3.)*a)));
				r4 = sqrt(x*x + (y - (3*sqrt(2.0)/2)*a)*(y - (3*sqrt(2.0)/2)*a) + z*z);
				m_d[i][j][k] = tanh((0.4*xlen/2.0-r1)/(sqrt(2.0)*eps)) 
					+ tanh((0.4*xlen/2.0-r2)/(sqrt(2.0)*eps)) 
					+ tanh((0.4*xlen/2.0-r3)/(sqrt(2.0)*eps))
					+ tanh((0.4*xlen/2.0-r4)/(sqrt(2.0)*eps)) + 3.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 4 spheres with same size mixed together.\n");
}

void fluidbase::init_3spheres_curv_touch()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2, r3;
				r1 = sqrt((x-0.1*xlen/2)*(x-0.1*xlen/2) + y*y + (z-0.5*xlen/2)*(z-0.5*xlen/2));
				r2 = sqrt((x+0.1*xlen/2)*(x+0.1*xlen/2) + y*y + z*z);
				r3 = sqrt((x-0.1*xlen/2)*(x-0.1*xlen/2) + y*y + (z+0.5*xlen/2)*(z+0.5*xlen/2));
				m_d[i][j][k] = tanh((0.35*xlen/2.0-r1)/(sqrt(2.0)*eps)) 
					+ tanh((0.35*xlen/2.0-r2)/(sqrt(2.0)*eps)) + tanh((0.35*xlen/2.0-r3)/(sqrt(2.0)*eps)) + 2.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 3 spheres with same size on a curve and touch.\n");
}

void fluidbase::init_3spheres_curv()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2, r3;
				r1 = sqrt((x-0.1*xlen/2)*(x-0.1*xlen/2) + y*y + (z-0.65*xlen/2)*(z-0.65*xlen/2));
				r2 = sqrt((x+0.1*xlen/2)*(x+0.1*xlen/2) + y*y + z*z);
				r3 = sqrt((x-0.1*xlen/2)*(x-0.1*xlen/2) + y*y + (z+0.65*xlen/2)*(z+0.65*xlen/2));
				m_d[i][j][k] = tanh((0.25*xlen/2.0-r1)/(sqrt(2.0)*eps)) 
					+ tanh((0.25*xlen/2.0-r2)/(sqrt(2.0)*eps)) + tanh((0.25*xlen/2.0-r3)/(sqrt(2.0)*eps)) + 2.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 3 spheres with same size on a curve.\n");
}

void fluidbase::init_8spheres()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double rs, r[8];
				rs = xlen/2/6;
				r[0] = sqrt((x-xlen/6)*(x-xlen/6) + (y+xlen/10)*(y+xlen/10) + z*z);
				r[1] = sqrt(x*x + y*y + (z-xlen/6)*(z-xlen/6));
				r[2] = sqrt((x-xlen/6*2)*(x-xlen/6*2) + y*y + (z+xlen/6)*(z+xlen/6));
				r[3] = sqrt((x-xlen/6)*(x-xlen/6) + (y+xlen/10)*(y+xlen/10) + (z-xlen/6*2)*(z-xlen/6*2));

				r[4] = sqrt((x+xlen/6)*(x+xlen/6) + (y-xlen/10)*(y-xlen/10) + z*z);
				r[5] = sqrt(x*x + y*y + (z+xlen/6)*(z+xlen/6));
				r[6] = sqrt((x+xlen/6*2)*(x+xlen/6*2) + y*y + (z-xlen/6)*(z-xlen/6));
				r[7] = sqrt((x+xlen/6)*(x+xlen/6) +(y-xlen/10)*(y-xlen/10) + (z+xlen/6*2)*(z+xlen/6*2));

				m_d[i][j][k] = 7;
				for (int s = 0; s < 8; s++) {
					m_d[i][j][k] += tanh((rs-r[s])/(sqrt(2.0)*eps));
				}
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 8 spheres on a plane.\n");
}

void fluidbase::init_5spheres()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double rs, r[5];
				rs = xlen*0.12;
				r[0] = sqrt((x+xlen/6)*(x+xlen/6) + (y-xlen*0.241)*(y-xlen*0.241) + z*z);
				r[1] = sqrt((x-xlen/6)*(x-xlen/6) + (y-xlen*0.241)*(y-xlen*0.241) + z*z);
				r[2] = sqrt((x-xlen*0.3)*(x-xlen*0.3) + (y+xlen*0.072)*(y+xlen*0.072) + (z+zlen*0.139)*(z+zlen*0.139));
				r[3] = sqrt(x*x + (y+xlen*0.222)*(y+xlen*0.222) + z*z);

				r[4] = sqrt((x+xlen*0.208)*(x+xlen*0.208) + (y+xlen*0.034)*(y+xlen*0.034) + (z-zlen*0.248)*(z-zlen*0.248));

				m_d[i][j][k] = 4;
				for (int s = 0; s < 5; s++) {
					m_d[i][j][k] += tanh((rs-r[s])/(sqrt(2.0)*eps));
				}
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 5 spheres not on a plane.\n");
}

/*
void fluidbase::init_5spheres()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double rs, r[5];
				rs = xlen/2/4;

				r[0] = sqrt((x+xlen*0.1815)*(x+xlen*0.1815) + (y-xlen*0.241)*(y-xlen*0.241) + z*z);
				r[1] = sqrt((x-xlen*0.1815)*(x-xlen*0.1815) + (y-xlen*0.241)*(y-xlen*0.241) + z*z);
				r[2] = sqrt((x-xlen*0.305)*(x-xlen*0.305) + (y+xlen*0.075)*(y+xlen*0.075) + (z+zlen*0.146)*(z+zlen*0.146));
				r[3] = sqrt((x+xlen*0.00247)*(x+xlen*0.00247) + (y+xlen*0.215)*(y+xlen*0.215) + (z+zlen*0.0)*(z+zlen*0.0));
				r[4] = sqrt((x+xlen*0.218)*(x+xlen*0.218) + (y+xlen*0.036)*(y+xlen*0.036) + (z-zlen*0.243)*(z-zlen*0.243));

				m_d[i][j][k] = 4;
				for (int s = 0; s < 5; s++) {
					m_d[i][j][k] += tanh((rs-r[s])/(sqrt(2.0)*eps));
				}
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 5 spheres not on a plane.\n");
}
*/
/*
void fluidbase::init_5spheres()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double rs, r[5];
				rs = xlen/2/4;
				//0.178, 0.172, 
				r[0] = sqrt((x+xlen*0.184)*(x+xlen*0.184) + (y-xlen/4)*(y-xlen/4) + z*z);
				r[1] = sqrt((x-xlen*0.184)*(x-xlen*0.184) + (y-xlen/4)*(y-xlen/4) + z*z);

				m_d[i][j][k] = 1;
				for (int s = 0; s < 2; s++) {
					m_d[i][j][k] += tanh((rs-r[s])/(sqrt(2.0)*eps));
				}
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 8 spheres on a plane.\n");
}
*/

void fluidbase::init_25spheres()
{
	int i, j, k, s1, s2;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				int nc = 3;
				double rs, r[9];
				rs = xlen/2/nc*0.85;
				for (s1 = 0; s1 < nc; s1++) {
					for (s2 = 0; s2 < nc; s2++) {
						double cx = -xlen/2 + xlen/nc*(s1+0.5);
						double cy = -ylen/2 + ylen/nc*(s2+0.5);
						r[s1*nc+s2] = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + z*z);
					}
				}
				m_d[i][j][k] = nc*nc-1;
				for (int s = 0; s < nc*nc; s++) {
					m_d[i][j][k] += tanh((rs-r[s])/(sqrt(2.0)*eps));
				}
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 25 spheres on a plane.\n");
}

void fluidbase::init_26spheres()
{
	int i, j, k, s1, s2;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				int nc = 3;
				double rs, r[9+4];
				rs = xlen/2/nc*0.85;
				for (s1 = 0; s1 < nc; s1++) {
					for (s2 = 0; s2 < nc; s2++) {
						double cx = -xlen/2 + xlen/nc*(s1+0.5);
						double cy = -ylen/2 + ylen/nc*(s2+0.5);
						r[s1*nc+s2] = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + z*z);
					}
				}
				for (s1 = 0; s1 < (nc-1); s1++) {
					for (s2 = 0; s2 < nc-1; s2++) {
						double cx = -xlen/2 + xlen/nc*(s1+1);
						double cy = -ylen/2 + ylen/nc*(s2+1);
						double cz = zlen/nc*0.707;
						r[nc*nc+s1*(nc-1)+s2] = sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz));
					}
				}

				m_d[i][j][k] = nc*nc-1+(nc-1)*(nc-1);
				for (int s = 0; s < nc*nc+(nc-1)*(nc-1); s++) {
					m_d[i][j][k] += tanh((rs-r[s])/(sqrt(2.0)*eps));
				}

			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as 26 spheres on a plane.\n");
}

//not finish
void fluidbase::init_cube()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double d1, d2, d3;
				d1 = fabs(x) - 0.5*xlen/2.0;
				d2 = fabs(y) - 0.5*ylen/2.0;
				d3 = fabs(z) - 0.5*zlen/2.0;
				m_d[i][j][k] = -tanh(d1/(sqrt(2.0)*eps)) -tanh(d2/(sqrt(2.0)*eps)) -tanh(d3/(sqrt(2.0)*eps)) -2;
				if (m_d[i][j][k] < -1.0) m_d[i][j][k] = -1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a cube.\n");
}


void fluidbase::init_9_cylinder()
{
	int i, j, k;
	double r = zlen/10;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double d1, d2, d3, d4, d5, d6, d7, d8, d9;
				d1 = sqrt(x*x+z*z);
				d2 = sqrt((x-3*r)*(x-3*r) + (z-3*r)*(z-3*r));
				d3 = sqrt((x-3*r)*(x-3*r) + (z+3*r)*(z+3*r));
				d4 = sqrt(x*x + (z-3*r)*(z-3*r));
				d5 = sqrt((x-3*r)*(x-3*r) + z*z);
				d6 = sqrt((x+3*r)*(x+3*r) + z*z);
				d7 = sqrt((x+3*r)*(x+3*r) + (z+3*r)*(z+3*r));
				d8 = sqrt((x+3*r)*(x+3*r) + (z-3*r)*(z-3*r));
				d9 = sqrt(x*x + (z+3*r)*(z+3*r));
				m_d[i][j][k] = tanh((d1-r)/(sqrt(2.0)*eps)) + tanh((d2-r)/(sqrt(2.0)*eps)) + tanh((d3-r)/(sqrt(2.0)*eps)) 
					+ tanh((d4-r)/(sqrt(2.0)*eps)) + tanh((d5-r)/(sqrt(2.0)*eps)) + tanh((d6-r)/(sqrt(2.0)*eps)) 
					+ tanh((d7-r)/(sqrt(2.0)*eps)) + tanh((d8-r)/(sqrt(2.0)*eps)) + tanh((d9-r)/(sqrt(2.0)*eps)) 
					- 8.0;
				m_d[i][j][k] *= -1.0; 

				if (fabs(y) > ylen/2*0.8) m_d[i][j][k] = -1.0;
				if (m_d[i][j][k] < -1.0) m_d[i][j][k] = -1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as 9 cylinder.\n");
}


void fluidbase::init_4_bings()
{
	int i, j, k;
	double r = zlen/16;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];

				if (fabs(x - r*2) < r || fabs(x - r*6) < r || fabs(x + r*2) < r || fabs(x + r*6) < r) m_d[i][j][k] = 1.0;
				else m_d[i][j][k] = -1.0;

				if (sqrt(y*y + z*z) > r*6) m_d[i][j][k] = -1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as 4 bings.\n");
}

void fluidbase::init_tian()
{
	int i, j, k;
	double r = zlen/12;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				m_d[i][j][k] = -1.0;
				if (x*x + y*y < r*r || x*x + z*z < r*r || z*z + y*y < r*r) m_d[i][j][k] = 1.0;
				if (x*x + (y-4*r)*(y-4*r) < r*r || x*x + (z-4*r)*(z-4*r) < r*r) m_d[i][j][k] = 1.0;
				if (y*y + (x-4*r)*(x-4*r) < r*r || y*y + (z-4*r)*(z-4*r) < r*r) m_d[i][j][k] = 1.0;
				if (z*z + (x-4*r)*(x-4*r) < r*r || z*z + (y-4*r)*(y-4*r) < r*r) m_d[i][j][k] = 1.0;
				if (x*x + (y+4*r)*(y+4*r) < r*r || x*x + (z+4*r)*(z+4*r) < r*r) m_d[i][j][k] = 1.0;
				if (y*y + (x+4*r)*(x+4*r) < r*r || y*y + (z+4*r)*(z+4*r) < r*r) m_d[i][j][k] = 1.0;
				if (z*z + (x+4*r)*(x+4*r) < r*r || z*z + (y+4*r)*(y+4*r) < r*r) m_d[i][j][k] = 1.0;

				if ((x-4*r)*(x-4*r) + (y-4*r)*(y-4*r) < r*r || (x-4*r)*(x-4*r) + (z-4*r)*(z-4*r) < r*r || (y-4*r)*(y-4*r) + (z-4*r)*(z-4*r) < r*r) m_d[i][j][k] = 1.0;
				if ((x+4*r)*(x+4*r) + (y+4*r)*(y+4*r) < r*r || (x+4*r)*(x+4*r) + (z+4*r)*(z+4*r) < r*r || (y+4*r)*(y+4*r) + (z+4*r)*(z+4*r) < r*r) m_d[i][j][k] = 1.0;

				if ((x-4*r)*(x-4*r) + (y+4*r)*(y+4*r) < r*r || (x-4*r)*(x-4*r) + (z+4*r)*(z+4*r) < r*r || (y-4*r)*(y-4*r) + (z+4*r)*(z+4*r) < r*r) m_d[i][j][k] = 1.0;
				if ((x+4*r)*(x+4*r) + (y-4*r)*(y-4*r) < r*r || (x+4*r)*(x+4*r) + (z-4*r)*(z-4*r) < r*r || (y+4*r)*(y+4*r) + (z-4*r)*(z-4*r) < r*r) m_d[i][j][k] = 1.0;

				if (fabs(x) > 5*r || fabs(y) > 5*r || fabs(z)>5*r) m_d[i][j][k] = -1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a tian.\n");
}





void fluidbase::init_ellipse()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y) + z*z/0.49);
				m_d[i][j][k] = tanh((0.7*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere with radius 0.7.\n");
}

void fluidbase::init_ellipse2()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y)/4 + z*z)*2;
				m_d[i][j][k] = tanh((0.83*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere the 2nd kind with radius 0.83.\n");
}

void fluidbase::init_ellipse3()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y)/25 + z*z)*5;
				m_d[i][j][k] = tanh((0.8*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere the 3rd kind 5:1 with radius 0.8.\n");
}

void fluidbase::init_ellipse4()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y)/9 + z*z)*3;
				m_d[i][j][k] = tanh((0.8*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere the 4th kind 3:1 with radius 0.8.\n");
}

void fluidbase::init_ellipse5()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y)/2.75/2.75 + z*z)*2.75;
				m_d[i][j][k] = tanh((0.8*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere the 5th kind 2.75:1 with radius 0.8.\n");
}

void fluidbase::init_ellipse_left()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(((x+xlen/4.0)*(x+xlen/4.0) + y*y)/2.75/2.75 + z*z)*2.75;
				m_d[i][j][k] = tanh((0.4*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an small ellipsoid at left 2.75:1 with radius 0.4.\n");
}

void fluidbase::init_ellipse6()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y)/2.75/2.75 + (z-zlen/2.0 + 0.5*ylen/2.0)*(z-zlen/2.0 + 0.5*ylen/2.0))*2.75;
				m_d[i][j][k] = tanh((0.7*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere the 6th kind 2.75:1 with radius 0.7.\n");
}

void fluidbase::init_ellipse7()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y)/1.2/1.2 + (z + 0.2*ylen/2.0)*(z + 0.2*ylen/2.0))*1.2;
				m_d[i][j][k] = tanh((0.5*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere the 7th kind 1.2:1 with radius 0.5.\n");
}

void fluidbase::init_ellipse8()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt((x*x + y*y)/4 + (z-zlen/2.0 + 0.5*ylen/2.0)*(z-zlen/2.0 + 0.5*ylen/2.0))*2;
				m_d[i][j][k] = tanh((0.5*ylen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as an ellipse sphere the 8th kind 2:1 with radius 0.5.\n");
}

void fluidbase::init_2ellipse()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1 = sqrt((x*x + y*y)/9 + (z-1.1)*(z-1.1))*3;
				double r2 = sqrt((x*x + y*y)/9 + (z+1.1)*(z+1.1))*3;
				m_d[i][j][k] = tanh((0.8*ylen/2.0-r1)/(0.5*sqrt(2.0)*eps)) + tanh((0.8*ylen/2.0-r2)/(0.5*sqrt(2.0)*eps)) + 1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as two ellipses.\n");
}

void fluidbase::init_triflower()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + y*y);
				if (r == 0) {
					m_d[i][j][k] = tanh(0.7/3*ylen/2.0/(0.5*sqrt(2.0)*eps));
					continue;
				}
				double theta = acos(y/r);
				double a = r/(cos(3*theta)+2);
				double d = (0.7/3*ylen/2.0 - a) * 3 / (cos(3*theta)+2);
				m_d[i][j][k] = tanh(d/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a triflower.\n");
}

void fluidbase::init_gtorus()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + y*y);
				if (r == 0) {
					m_d[i][j][k] = tanh(-0.5*ylen/2.0/(sqrt(2.0)*eps));
					continue;
				}
				double theta = acos(x/r);
				if (y < 0) theta = 2*PI - theta;
				double d = 0.1*ylen/2.0*(0.75*cos(theta)+2.25) - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (y - 0.5*ylen/2.0*sin(theta))*(y - 0.5*ylen/2.0*sin(theta)) + z*z);
				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a gtorus.\n");
}

void fluidbase::init_gtorus2()
{
//	double s = 2*xlen/2.0/11.0;
	double s = 2*xlen/2.0/13.0;
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k] + s;
				double r = sqrt(x*x + y*y);
				double theta = 0;
				if (r != 0) theta = acos(x/r);
				if (y < 0) theta = 2*PI - theta;
				double d1 = s*(cos(theta)+2.0) - sqrt((x - 2.5*s*cos(theta))*(x - 2.5*s*cos(theta)) + (y - 2.5*s*sin(theta))*(y - 2.5*s*sin(theta)) + z*z);
				theta += PI;
				double d2 = s*(cos(theta)+2.0) - sqrt((x - 2.5*s*cos(theta))*(x - 2.5*s*cos(theta)) + (y - 2.5*s*sin(theta))*(y - 2.5*s*sin(theta)) + z*z);
				double d = d1;
				if (fabs(d2) < fabs(d1)) d = d2;
				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as the second kind of gtorus.\n");
}

void fluidbase::init_vtorus()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + z*z);
				if (r == 0) {
					m_d[i][j][k] = tanh(-0.5*ylen/2.0/(sqrt(2.0)*eps));
					continue;
				}
				double theta = acos(x/r);
				if (z < 0) theta = 2*PI - theta;
				double d = 0.25*ylen/2.0 - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (z - 0.5*zlen/2.0*sin(theta))*(z - 0.5*zlen/2.0*sin(theta)) + y*y);
				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a vtorus.\n");
}

void fluidbase::init_htorus()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + y*y);
				if (r == 0) {
					m_d[i][j][k] = tanh(-0.5*ylen/2.0/(sqrt(2.0)*eps));
					continue;
				}
				double theta = acos(x/r);
				if (y < 0) theta = 2*PI - theta;
				double d = 0.25*ylen/2.0 - sqrt((x - 0.5*xlen/2.0*cos(theta))*(x - 0.5*xlen/2.0*cos(theta)) + (y - 0.5*ylen/2.0*sin(theta))*(y - 0.5*ylen/2.0*sin(theta)) + z*z);
				m_d[i][j][k] = tanh(d/(sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/4;
	output("Starting as a htorus.\n");
}

void fluidbase::init_hulu()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = sqrt(x*x + y*y + z*z/9)*3;
				m_d[i][j][k] = tanh((0.8*zlen/2.0-r)/(0.5*sqrt(2.0)*eps));
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a hulu.\n");
}

void fluidbase::init_art_cherry()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1, r2;
				r1 = sqrt(x*x + y*y + z*z);
				r2 = sqrt(x*x + y*y + (z-0.15*ylen/2.0)*(z-0.15*ylen/2.0));
				m_d[i][j][k] = tanh((0.6*ylen/2.0-r1)/(sqrt(2.0)*eps)) 
							+ tanh((r2-0.4*ylen/2.0)/(sqrt(2.0)*eps)) - 1;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/5;
	output("Starting as a artifical cherry.\n");
}

void fluidbase::init_cherry()
{
	input("../cherry", m_d);
	if (m_angle < 0) m_angle = PI/5;
	output("Starting as a cherry.\n");
}

void fluidbase::init_drop()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				m_d[i][j][k] = -1.0;
			}
		}
	}
	input("../drop", m_d);
	if (m_angle < 0) m_angle = PI/5;
	output("Starting as a drop. eps must be 1.5h\n");
}

void fluidbase::init_disk()
{
//	input("../disksmall", m_d);
	input("../diskflat", m_d);
	if (m_angle < 0) m_angle = PI/5;
	output("Starting as a dimped disk.\n");
}

void fluidbase::init_cherry_inverse()
{
	input("../cherry", m_d);
	int i, j, k;
	for (i = 0; i < nz/2; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				swap(m_d[i][j][k], m_d[nz-i][j][k]);
			}
		}
	}
	if (m_angle < 0) m_angle = PI/5;
	output("Starting as a inverse cherry.\n");
}

void fluidbase::init_Gyroid()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r = cos(x)*sin(y) + cos(y)*sin(z) + cos(z)*sin(x);
				m_d[i][j][k] = r>0?1.0:-1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a Gyroid.\n");
}

void fluidbase::init_3cylinder()
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double r1 = sqrt(x*x+y*y);
				double r2 = sqrt(x*x+z*z);
				double r3 = sqrt(y*y+z*z);
				double r = PI/4;
				m_d[i][j][k] = (r1>r && r2>r && r3>r)?1.0:-1.0;
			}
		}
	}
	if (m_angle < 0) m_angle = PI/6;
	output("Starting as a Gyroid.\n");
}

void fluidbase::init_curvature1()
{
	int i, j, k;
	malloc_array(C);
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				C[i][j][k] = 1.0;
			}
		}
	}
	output("Unified spontaneous curvature 1.0.\n");
}

void fluidbase::init_curvature2()
{
	int i, j, k;
	malloc_array(C);
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double a, b, c, t;

				a = 1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d1 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = -1; b = 1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d2 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = -1; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d3 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				a = 1; b = 1; c = -1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d4 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				if (d1 < 0.15*ylen/2.0 || d2 < 0.15*ylen/2.0 || d3 < 0.15*ylen/2.0 || d4 < 0.15*ylen/2.0) C[i][j][k] = 6.0;
				else C[i][j][k] = 0.0;
			}
		}
	}
	output("Eight corner spontaneous curvature 6.0.\n");
}

void fluidbase::init_curvature3()
{
	int i, j, k;
	malloc_array(C);
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				double a, b, c, t;

				a = 0; b = 0; c = 1;
				t = (a*x+b*y+c*z)/(a*a+b*b+c*c);
				double d1 = sqrt((x-t*a)*(x-t*a) + (y-t*b)*(y-t*b) + (z-t*c)*(z-t*c));

				if (d1 < 0.15*ylen/2.0) C[i][j][k] = 6.0;
				else C[i][j][k] = 0.0;
			}
		}
	}
	output("up-down spontaneous curvature 6.0.\n");
}

void fluidbase::init_fluid_zero()
{
	zero_array(m_u);
	zero_array(m_v);
	zero_array(m_w);
	output("with no flush; ");
}

void fluidbase::init_fluid_zaxis()
{
	double maxspeed = 2*ylen/2.0*10;
	zero_array(m_u);
	zero_array(m_v);
	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				m_w[i][j][k] = maxspeed;
			}
		}
	}
	output("with z-axis fluid; ");
}

void fluidbase::init_fluid_flush()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*60; //old  maxspeed = 2*ylen/2.0*40;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * exp(-(x*x+y*y)*3/r/r);
			}
		}
	}
	output("with fluid flush 60; ");
}

void fluidbase::init_fluid_flush_middle()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*25;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * exp(-(x*x+y*y)*3/r/r);
			}
		}
	}
	output("with fluid flush middle speed; ");
}
/*
void fluidbase::init_fluid_flush_cos()
{
	double maxspeed = 2*ylen/2.0*5;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * (cos(x/2)*cos(y/2));
			}
		}
	}
	output("with fluid flush cos slow; ");
}
*/

void fluidbase::init_fluid_flush_cos()
{
	double maxspeed = 2*ylen/2.0*10;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
//				m_w[i][j][k] = maxspeed *  (cos(x)+1)*(cos(y)+1);
				m_w[i][j][k] = maxspeed * (cos(x/2)*cos(y/2));
			}
		}
	}
	output("with fluid flush cos slow; ");
}

void fluidbase::init_fluid_flush_cos_fast()
{
	double maxspeed = 2*ylen/2.0*40;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * (cos(x/2)*cos(y/2));
			}
		}
	}
	output("with fluid flush cos fast; ");
}

void fluidbase::init_fluid_flush_slow()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*10;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * exp(-(x*x+y*y)*3/r/r);
			}
		}
	}
	output("with fluid flush slow; ");
}

void fluidbase::init_fluid_flushxz()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*40;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * (exp(-(x*x)*3/r/r) - 0.5);
			}
		}
	}
	output("with fluid flush jet at x-z plan; ");
}

void fluidbase::init_fluid_flushxz2()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*10;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = -maxspeed * (exp(-(x*x)*3/r/r));
				m_w[i][j][k] = 0.0;
			}
		}
	}
	output("with fluid flush jet at x-z plan 2; ");
}

void fluidbase::init_fluid_flushp()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*40;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * exp(-((x+0.5*xlen/2.0)*(x+0.5*xlen/2.0)+y*y)*3/r/r);
			}
		}
	}
	output("with fluid flush jet at left point; ");
}

void fluidbase::init_fluid_flushp_slow()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*10;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * exp(-((x+0.5*xlen/2.0)*(x+0.5*xlen/2.0)+y*y)*3/r/r);
			}
		}
	}
	output("with fluid flush jet at left point -- slow; ");
}

void fluidbase::init_fluid_flushp_middle()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*30;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * exp(-((x+0.5*xlen/2.0)*(x+0.5*xlen/2.0)+y*y)*3/r/r);
			}
		}
	}
	output("with fluid flush jet at left point -- middle; ");
}

void fluidbase::init_fluid_curve()
{
	double r = 0.25*ylen/2.0;
	double maxspeed = 2*ylen/2.0*40;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * cos(x/2)*cos(y/2) * cos(1.5*x)*cos(1.5*y);
			}
		}
	}
	output("with fluid flush curve; ");
}

void fluidbase::init_fluid_tanh()
{
	double r = 1;
	double maxspeed = 2*ylen/2.0*40;

	int i, j, k;
	double x, y, z;
	for (i = 0; i < nz; i++) {
		z = zz[i];
		for (j = 0; j < ny; j++) {
			y  = yy[j];
			for (k = 0; k < nx; k++) {
				x = xx[k];
				m_u[i][j][k] = 0.0;
				m_v[i][j][k] = 0.0;
				m_w[i][j][k] = maxspeed * tanh((xlen/2 - fabs(x))/r)*tanh((ylen/2 - fabs(y))/r);
			}
		}
	}
	output("with fluid flush tanh; ");
}

void fluidbase::reload()
{
	if (nflush != -1) {
		input("u", m_u);
		input("v", m_v);
		input("w", m_w);
	}
	input("d", m_d);
	input("e", m_e);
}


void fluidbase::avg_mesh(ARR3D source, ARR3D d)
{
	extern double N[16][16];
	int i, j, k, r, s, t;
	for (i = 0; i < nz; i++) {
		double z = (2.0*i/nz - 1)*fabs(zz[0]);
		if (zz[0] < zz[1]) {
			for (t = 0; t < nz; t++)  
				if (z < zz[t]) break;
		}
		else {
			for (t = 0; t < nz; t++)  
				if (z > zz[t]) break;
		}
		t --;
		for (j = 0; j < ny; j++) {
			double y = (2.0*j/ny - 1)*fabs(yy[0]);
			if (yy[0] < yy[1]) {
				for (s = 0; s < ny; s++)  
					if (y < yy[s]) break;
			}
			else {
				for (s = 0; s < ny; s++)  
					if (y > yy[s]) break;
			}
			s --;
			for (k = 0; k < nx; k++) {
				double x = (2.0*k/nx - 1)*fabs(xx[0]);
				if (xx[0] < xx[1]) {
					for (r = 0; r < nx; r++)  if (x < xx[r]) break;
				}
				else {
					for (r = 0; r < nx; r++)  if (x > xx[r]) break;
				}
				r --;

				double tempx = (xx[r] - x) / (xx[r] - xx[r+1]); 
				double tempy = (yy[s] - y) / (yy[s] - yy[s+1]);
				double tempz = (zz[t] - z) / (zz[t] - zz[t+1]);
				double cor[2][2][2];
				int ii, jj, kk;
				for (ii = 0; ii < 2; ii++)
					for (jj = 0; jj < 2; jj++)
						for (kk = 0; kk < 2; kk++)
							cor[ii][jj][kk] = source[(t+ii)%nz][(s+jj)%ny][(r+kk)%nx];

				double ret = 0;
				double px, py, pz;
				for (ii = 0; ii < 2; ii++) {
					if (ii == 0) pz = 1-tempz;
					else pz = tempz;
					for (jj = 0; jj < 2; jj++) {
						if (jj == 0) py = 1-tempy;
						else py = tempy;
						for (kk = 0; kk < 2; kk++) {
							if (kk == 0) px = 1-tempx;
							else px = tempx;
							ret += cor[ii][jj][kk]*px*py*pz;
						}
					}
				}
				d[i][j][k] = ret;
			}
		}
	}
}

void fluidbase::openmp_info()
{
#if defined (_OPENMP) 

	output("******** OPENMP INFO ********\n");
	sprintf(sout, "CPUs = %d \nMax Threads = %d\n", omp_get_num_procs(), omp_get_max_threads());
#pragma omp parallel
	{
		char sout[100];
		sprintf(sout, "Thread #%d\n", omp_get_thread_num());
		output(sout);
	}

	output(sout);
	output("******** *********** ********\n");

#endif
}

//this function save all the other informations, such as surface fix size,  eps, angle, 
void fluidbase::output_info()
{
	char fname[100];
	sprintf(fname, "%s%dx%dx%d.info", folderpath, nx, ny, nz);
	FILE * f = fopen(fname, "wb");
	if (!f) {
		output("Can not open the information file for write: ");
		output(fname);
		output("\n");
		return;
	}
	fwrite(&time_current, 1, sizeof(double), f);
	fwrite(&time_step, 1, sizeof(double), f);
	fwrite(&alpha, 1, sizeof(double), f);
	fwrite(&beta, 1, sizeof(double), f);
	fwrite(&m_alpha3, 1, sizeof(double), f);
	fwrite(&alpha_1, 1, sizeof(double), f);
	fwrite(&beta_1, 1, sizeof(double), f);
	fclose(f);
}

void fluidbase::load_info()
{
	char fname[100];
	sprintf(fname, "%s%dx%dx%d.info", folderpath, nx, ny, nz);
	FILE * f = fopen(fname, "rb");
	if (!f) {
		output("Starting a new experiment.... start time is 0.\n");
		return;
	}
	else {
		output("Continue the existing experiment....");
	}
	fread(&time_current, 1, sizeof(double), f);
	fread(&time_step, 1, sizeof(double), f);
	fread(&alpha, 1, sizeof(double), f);
	fread(&beta, 1, sizeof(double), f);
	fread(&m_alpha3, 1, sizeof(double), f);
	fread(&alpha_1, 1, sizeof(double), f);
	fread(&beta_1, 1, sizeof(double), f);
	fclose(f);

	if (time_step < 1e-13) time_step = 1e-7;
	if (time_current == 0.0) time_current = 1e-13;

	time_start = time_current;
	char outs[100];
	sprintf(outs, "Set start time to be %12.8f.\n", time_start);
	output(outs);
}

//can not used in parallel routine.
void fluidbase::shift(ARR3D u, int offset)
{
	ARR3D tmp1;
	malloc_array(tmp1);
	int i, j, k;

	if (offset > 0) {
		int rare = nz - offset;
		for (i = rare; i < nz; i++) {
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					tmp1[i-rare][j][k] = u[i][j][k];
				}
			}
		}
		for (i = nz - 1; i >=offset; i--) {
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					u[i][j][k] = u[i-offset][j][k];
				}
			}
		}
		for (i = 0; i < offset; i++) {
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					u[i][j][k] = tmp1[i][j][k];
				}
			}
		}
	}
	else if (offset < 0) {
		int rare = nz + offset;
		for (i = rare; i < nz; i++) {
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					tmp1[i][j][k] = u[i-rare][j][k];
				}
			}
		}
		for (i = 0; i < rare; i++) {
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					u[i][j][k] = u[i-offset][j][k];
				}
			}
		}
		for (i = rare; i < nz; i++) {
			for (j = 0; j < ny; j++) {
				for (k = 0; k < nx; k++) {
					u[i][j][k] = tmp1[i][j][k];
				}
			}
		}
	}
	free_array(tmp1);
}

bool fluidbase::check_folder()
{
	char inifile[100];
	sprintf(inifile, "%sini%d.dat", RESULT_PATH, folderindex);
	FILE * f = fopen(inifile, "r");
	if (!f) {
		return false;
	}
	int ifold;
	fscanf(f, "%d", &ifold);
	fclose(f);

	sprintf(inifile, "%s%d/", RESULT_PATH, ifold);
	if (strcmp(inifile, folderpath)) {
		output("User terminate!\n");
		return true;
	}
	else return false;
}

void fluidbase::changeshape(ARR3D s, ARR3D d)
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		double z = zz[i];
		for (j = 0; j < ny; j++) {
			double y  = yy[j];
			for (k = 0; k < nx; k++) {
				double x = xx[k];
				int kk = nx/2 + (k - nx/2)*(1 + 0.2 * sin(y));
				int jj = ny/2 + (j - ny/2)*(1 + 0.2 * sin(z));
				int ii = nz/2 + (i - ny/2)*(1 + 0.2 * sin(x));
				
				if (ii >= nz || jj >= ny || kk >= nx || ii<0 || jj <0 || kk < 0) d[i][j][k] = -1.0;
				else d[i][j][k] = s[ii][jj][kk];
			}
		}
	}
}

void fluidbase::scale(ARR3D s, ARR3D d, double kx, double ky, double kz)
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				int kk = nx/2 + (k - nx/2)/kx;
				int jj = ny/2 + (j - ny/2)/ky;
				int ii = nz/2 + (i - ny/2)/kz;
				
				if (ii >= nz || jj >= ny || kk >= nx || ii<0 || jj <0 || kk < 0) d[i][j][k] = -1.0;
				else d[i][j][k] = s[ii][jj][kk];
			}
		}
	}
}

void fluidbase::shiftright(ARR3D s, ARR3D d, int offset)
{
	int i, j, k;
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nx; k++) {
				int kk = k - offset;
				int jj = j;
				int ii = i;
				
				if (ii >= nz || jj >= ny || kk >= nx || ii<0 || jj <0 || kk < 0) d[i][j][k] = -1.0;
				else d[i][j][k] = s[ii][jj][kk];
			}
		}
	}
}
