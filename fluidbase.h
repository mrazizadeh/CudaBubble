// fluidbase.h: interface for the fluidbase class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FLUIDBASE_H__B20AB0B6_7E23_43CD_BC79_9DDB6B550BA2__INCLUDED_)
#define AFX_FLUIDBASE_H__B20AB0B6_7E23_43CD_BC79_9DDB6B550BA2__INCLUDED_

#include "tool.h"
#include "Calc.h"

#ifdef WIN32
//#	define RESULT_PATH "F:\\work\\bfluid3d\\results\\"
#	define RESULT_PATH "./results/"
#else
#	define RESULT_PATH "./results/"
#endif

#define TMP_ARRAY_MAX 40  //large number does not matter, dynamic memory allocation  

class fluidbase  
{
	friend class mem_array;
public:
	void shiftright(ARR3D s, ARR3D d, int offset);
	void scale(ARR3D s, ARR3D d, double kx, double ky, double kz);
	void changeshape(ARR3D s, ARR3D x);
	virtual ~fluidbase();
	fluidbase(int folderindex);

protected:
	bool check_folder();
	void shift(ARR3D u, int offset);
	void load_info();
	void avg_mesh(ARR3D s, ARR3D d);
	void release_tmp_array(int i);
	void release_tmp_complex_array(int i);
	DARR3D  get_tmp_array(int &i);
	cufftDoubleComplex* get_tmp_complex_array(int& ind);
	int tmp_array_len, tmp_complex_array_len;

	void output_info();
	void output(const char * s);
	void output(char * name, double x);
	void output(char *name, double x1, double x2);
	void output(char *name, double x1, double x2, double x3, double x4 = 0.0, double x5 = 0.0, double x6 = 0.0, double x7 = 0.0, double x8 = 0.0, double x9 = 0.0);
	void output_graph(char * name, ARR3D u, ARR3D phase, double times = 0.0, double zero = 0.0, double rotateangle = PI/4, double diff = 1.0);
	void output_graph_two_bubble(char * name, ARR3D u, ARR3D e, double times = 0.0, double zero = 0.0, double rotateangle = PI/4, double diff = 1.0);
	void output_vector3d(char * name, ARR3D u, ARR3D v, ARR3D w);
	double time_start;
	double time_step, fluid_time_step, fixed_time_step;
	double time_current;
	double alpha, beta, m_alpha3;
	double alpha_1, beta_1;
	double m_k, m_c, m_linet;  //bending rigity and line tension
	char folderpath[200];
	void input(char * name, ARR3D u);
	void output(char *name, ARR3D u);
	bool malloc_array(ARR3D  &x);
	int nx, ny, nz;
	double xlen, ylen, zlen;   //mesh size
	ARR3D m_u, m_v, m_w;
	ARR3D m_d, m_e;
	double *xx, *yy, *zz;

	ARR3D C;
	double eps, xi;
	double m_angle;
	int b_shape;

	void reload();
	void init();

	int nflush, ninit, einit, ncur;
	double surctrl;
	//the following are init parameters.
	double M1, M2, M3, M4, M5;
	double eta;   //surface tension coefficient
	double re;
	double ga;
	double tmax; 
	int output_interval;
	double gweight;
	char * sout;
	int folderindex;

private:
	void openmp_info();
	void init_eta1();
	void init_eta2();
	void init_eta3();
	void init_eta4();
	void init_eta5();
	void init_eta6();
	void init_eta8();
	void init_eta9();
	void init_eta10();
	void init_eta11();
	void init_eta12();
	void init_eta13();

	void init_1sphere();
	void init_1sphere_right();
	void init_rightsphere();
	void init_etaleftsphere();
	void init_2spheres();
	void init_3spheres();
	void init_3spheres_curv();
	void init_3spheres_curv_touch();
	void init_3spheres_touch();
	void init_4spheres_touch();
	void init_cube();
	void init_ellipse();
	void init_ellipse2();
	void init_ellipse3();
	void init_ellipse4();
	void init_ellipse5();
	void init_ellipse_left();
	void init_ellipse6();
	void init_ellipse7();
	void init_ellipse8();
	void init_triflower();
	void init_gtorus();
	void init_gtorus2();
	void init_vtorus();
	void init_htorus();
	void init_2ellipse();
	void init_hulu();
	void init_cherry();
	void init_drop();
	void init_disk();
	void init_3spheres_vertical();
	void init_ciga_vertical();
	void init_ciga_horizontal();
	void init_torus_small_ball();
	void init_torus_big_ball();
	void init_torus_bing();
	void init_two_torus_ring();
	void init_art_cherry();
	void init_cherry_inverse();
	void init_1sphere_small_right();
	void init_2sphere_small();
	void init_ciga_left();
	void init_ciga_two();
	void init_ciga_long_thin();
	void init_8spheres();
	void init_25spheres();
	void init_26spheres();
	void init_5spheres();
	void init_9_cylinder();
	void init_4_bings();
	void init_tian();
	void init_cross1();
	void init_cross2();
	void init_Gyroid();
	void init_3cylinder();

	void init_curvature3();
	void init_curvature2();
	void init_curvature1();

	void init_fluid_zero();
	void init_fluid_flush();
	void init_fluid_flush_middle();
	void init_fluid_flush_slow();
	void init_fluid_flushp();
	void init_fluid_flushp_slow();
	void init_fluid_flushp_middle();
	void init_fluid_flush_cos();
	void init_fluid_flush_cos_fast();
	void init_fluid_flushxz();
	void init_fluid_flushxz2();
	void init_fluid_zaxis();
	void init_fluid_curve();
	void init_fluid_tanh();
	bool tmp_array_used[TMP_ARRAY_MAX];
	DARR3D base_tmp_array[TMP_ARRAY_MAX];
	bool tmp_complex_array_used[TMP_ARRAY_MAX];
	cufftDoubleComplex * base_tmp_complex_array[TMP_ARRAY_MAX];

};

#endif // !defined(AFX_FLUIDBASE_H__B20AB0B6_7E23_43CD_BC79_9DDB6B550BA2__INCLUDED_)
