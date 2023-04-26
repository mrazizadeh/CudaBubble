#include "stdafx.h"
#include "fluidbase.h"
//#include "axisfluidbase.h"

#ifdef __cplusplus
	extern "C" {
#endif // __cplusplus

//#include "JpegLib/jpeglib.h"
//#include "Jpegfile.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#ifdef __cplusplus
	}
#endif // __cplusplus

struct zpoint 
{
	int x; int y; int z;
	int comp;
	double n[3];
	bool isedge;
	bool iscut;
}; 

void get_around(ARR3D u, double d[2][2][2], int x, int y, int z, int nz, int ny, int nx, double times)
{
	int i, j, k;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			for (k = 0; k < 2; k++)
				d[i][j][k] = u[(z+i+nz)%nz][(y+j+ny)%ny][(x+k+nx)%nx];
}

/*
void read_N()
{
	FILE * f = fopen("f:\\work\\n.txt", "r");
	for(int i = 0; i < 16; i++) {
		for (int j = 0; j < 16; j++) fscanf(f, "%e", &(N[i][j]));
	}
	fclose(f);
}
*/
double get_value(ARR3D u, double x, double y, double z, int nz, int ny, int nx, double times)
{
	int wx, wy, wz;
	wx = int(nx*times)/2;
	wy = int(ny*times)/2;
	wz = int(nz*times)/2;
	int ix = int((x+wx)/times);
	int iy = int((y+wy)/times);
	int iz = int((z+wz)/times);
	double d[2][2][2];
	get_around(u, d, ix, iy, iz, nz, ny, nx , times);

	double tempx = (x+wx)/times-ix;
	double tempy = (y+wy)/times-iy;
	double tempz = (z+wz)/times-iz;

	double ret = 0;
	int i, j, k;
	double px, py, pz;
	for (i = 0; i < 2; i++) {
		if (i == 0) pz = 1-tempz;
		else pz = tempz;
		for (j = 0; j < 2; j++) {
			if (j == 0) py = 1-tempy;
			else py = tempy;
			for (k = 0; k < 2; k++) {
				if (k == 0) px = 1-tempx;
				else px = tempx;
				ret += d[i][j][k]*px*py*pz;
			}
		}
	}

	return ret;
}

void get_n(ARR3D u, double n[3], double x, double y, double z, int nz, int ny, int nx, double times)
{
	int wx, wy, wz;
	wx = int(nx*times)/2;
	wy = int(ny*times)/2;
	wz = int(nz*times)/2;
	int ix = int((x+wx)/times);
	int iy = int((y+wy)/times);
	int iz = int((z+wz)/times);
	double d[2][2][2];
	get_around(u, d, ix, iy, iz, nz, ny, nx , times);

	double tempx = (x+wx)/times-ix;
	double tempy = (y+wy)/times-iy;
	double tempz = (z+wz)/times-iz;

	int i, j, k;
	double px, py, pz;

	n[0] = 0;
	for (i = 0; i < 2; i++) {
		if (i == 0) pz = 1-tempz;
		else pz = tempz;
		for (j = 0; j < 2; j++) {
			if (j == 0) py = 1-tempy;
			else py = tempy;
			n[0] += (d[i][j][1] - d[i][j][0])*py*pz;
		}
	}

	n[1] = 0;
	for (i = 0; i < 2; i++) {
		if (i == 0) pz = 1-tempz;
		else pz = tempz;
		for (k = 0; k < 2; k++) {
			if (k == 0) px = 1-tempx;
			else px = tempx;
			n[1] += (d[i][1][k] - d[i][0][k])*px*pz;
		}
	}

	n[2] = 0;
	for (j = 0; j < 2; j++) {
		if (j == 0) py = 1-tempy;
		else py = tempy;
		for (k = 0; k < 2; k++) {
			if (k == 0) px = 1-tempx;
			else px = tempx;
			n[2] += (d[1][j][k] - d[0][j][k])*px*py;
		}
	}

	if (n[0] == 0.0 && n[1] == 0.0 && n[2] == 0.0) {
		n[0] = n[1] = n[2] = 1.0;
	}

}

/*
double R[3][3] = {    0.5000,         0,    0.8660,
    0.7500,    0.5000,   -0.4330,
	-0.4330,    0.8660,    0.2500};
*/

double RANGLE;

void get_zero_points(ARR3D u, ARR3D phase, int nz, int ny, int nx, double times, zpoint * zp, int & nzp, double zero)
{
	int i, j, k;
	int x1, y1, z1;
	double x, y, z;

	double R[3][3] = {	1.0,         0.0,    0.0,
					0.000,    cos(RANGLE),		-sin(RANGLE),
					0.0,    +sin(RANGLE),    cos(RANGLE)};

	int wx, wy, wz;
	wx = int((nx-1)*times)/2;
	wy = int(((ny-1)*cos(RANGLE) + (nz*2-1)*sin(RANGLE))*times/2);
	wz = int(((ny-1)*sin(RANGLE) + (nz*2-1)*cos(RANGLE))*times/2);

	int rx, ry, rz;
	rx = int(nx*times)/2;
	ry = int(ny*times)/2;
	rz = int(nz*2*times)/2;

	for (y1 = -wy; y1 <= wy; y1++) {
		for (x1 = -wx; x1 <= wx; x1++) {
			for (z1 = -wz; z1 <= wz; z1++)
			{
				bool isedge = false;
				x = R[0][0]*x1 + R[0][1]*y1 + R[0][2]*z1;
				if (x < -rx + 1 || x >= rx) continue;
				y = R[1][0]*x1 + R[1][1]*y1 + R[1][2]*z1;
				if (y < -ry + 1 || y >= ry) continue;
				z = R[2][0]*x1 + R[2][1]*y1 + R[2][2]*z1;
				if (z < -rz + 1 || z >= rz) continue;

				if (y < fabs(x) && y+1>fabs(x)) isedge = true;

				double tempx = x / times;
				double tempy = y / times;
				double tempz = z / times;
				int ix = int((x+rx)/times);
				int iy = int((y+ry)/times);
				int iz = int((z+rz)/times);
				if (ix < 0 || ix > nx || iy < 0 || iy > ny || iz < 0 || iz > 2*nz) continue;
				ix %= nx;  iy %= ny;   iz %= nz;
				if (z >= 0) z -= nz*times/2;
				else z += nz*times/2;
				//if (fabs(u[iz][iy][ix]) > 0.9) continue;

				double d[2][2][2];
				get_around(u, d, ix, iy, iz, nz, ny, nx, times);

				bool l = false;
				for (i = 0; i < 2; i ++) 
					for (j = 0; j < 2; j++) 
						for (k = 0; k < 2; k++) 
							if (!l && (d[i][j][k]-zero)*(d[0][0][0]-zero) < 0) l = true;
				if (!l) continue;

				double v[7];
				v[0] = get_value(u, x, y, z, nz, ny, nx, times) - zero;
				if (fabs(v[0]) > 0.9 /times) continue;

				v[1] = get_value(u, x+1, y, z, nz, ny, nx, times) - zero;
				v[2] = get_value(u, x-1, y, z, nz, ny, nx, times) - zero;
				v[3] = get_value(u, x, y+1, z, nz, ny, nx, times) - zero;
				v[4] = get_value(u, x, y-1, z, nz, ny, nx, times) - zero;
				v[5] = get_value(u, x, y, z+1, nz, ny, nx, times) - zero;
				v[6] = get_value(u, x, y, z-1, nz, ny, nx, times) - zero;
				if (v[1]*v[2] > 0 && v[3]*v[4] > 0 && v[5]*v[6] > 0) continue;

				zp[nzp].x = x1;
				zp[nzp].y = y1;
				zp[nzp].z = wz-z1;
				get_n(u, zp[nzp].n, x, y, z, nz, ny, nx, times);
				zp[nzp].isedge = isedge;
				if (y > 0 && y > fabs(x)) zp[nzp].iscut = true;
				else zp[nzp].iscut = false;
				if (phase == NULL) zp[nzp].comp = -1;
				else {
					double v = get_value(phase, x, y, z, nz, ny, nx, times);
					//if (fabs(v) > 0.9) continue;
					if (v > 0.0) zp[nzp].comp = 1;
//					if (phase[iz][iy][ix] > 0.0) zp[nzp].comp = 1;
					else zp[nzp].comp = -1;
				}
				nzp ++;
			}
		}
	}
}

void get_zero_points_two_bubble(ARR3D u, ARR3D e, ARR3D phase, int nz, int ny, int nx, double times, zpoint * zp, int & nzp, double zero)
{
	int i, j, k;
	int x1, y1, z1;
	double x, y, z;

	double R[3][3] = {	1.0,         0.0,    0.0,
					0.000,    cos(RANGLE),		-sin(RANGLE),
					0.0,    +sin(RANGLE),    cos(RANGLE)};

	int wx, wy, wz;
	wx = int((nx-1)*times)/2;
	wy = int(((ny-1)*cos(RANGLE) + (nz*2-1)*sin(RANGLE))*times/2);
	wz = int(((ny-1)*sin(RANGLE) + (nz*2-1)*cos(RANGLE))*times/2);

	int rx, ry, rz;
	rx = int(nx*times)/2;
	ry = int(ny*times)/2;
	rz = int(nz*2*times)/2;

	for (y1 = -wy; y1 <= wy; y1++) {
		for (x1 = -wx; x1 <= wx; x1++) {
			for (z1 = -wz; z1 <= wz; z1++)
			{
				bool isedge = false;
				x = R[0][0]*x1 + R[0][1]*y1 + R[0][2]*z1;
				if (x < -rx + 1 || x >= rx) continue;
				y = R[1][0]*x1 + R[1][1]*y1 + R[1][2]*z1;
				if (y < -ry + 1 || y >= ry) continue;
				z = R[2][0]*x1 + R[2][1]*y1 + R[2][2]*z1;
				if (z < -rz + 1 || z >= rz) continue;

				if (y < fabs(x) && y+1>fabs(x)) isedge = true;

				double tempx = x / times;
				double tempy = y / times;
				double tempz = z / times;
				int ix = int((x+rx)/times);
				int iy = int((y+ry)/times);
				int iz = int((z+rz)/times);
				if (ix < 0 || ix > nx || iy < 0 || iy > ny || iz < 0 || iz > 2*nz) continue;
				ix %= nx;  iy %= ny;   iz %= nz;
				if (z >= 0) z -= nz*times/2;
				else z += nz*times/2;
				//if (fabs(u[iz][iy][ix]) > 0.9) continue;

				double du[2][2][2], de[2][2][2];
				get_around(u, du, ix, iy, iz, nz, ny, nx, times);
				get_around(e, de, ix, iy, iz, nz, ny, nx, times);

				bool lu = false, le = false;
				for (i = 0; i < 2; i ++) 
					for (j = 0; j < 2; j++) 
						for (k = 0; k < 2; k++) { 
							if (!lu && (du[i][j][k]-zero)*(du[0][0][0]-zero) < 0) lu = true;
							if (!le && (de[i][j][k]-zero)*(de[0][0][0]-zero) < 0) le = true;
						}
				if (!lu && !le) continue;

				double vu[7], ve[7];
				if (lu) vu[0] = get_value(u, x, y, z, nz, ny, nx, times) - zero;
				if (le) ve[0] = get_value(e, x, y, z, nz, ny, nx, times) - zero;
				if (fabs(vu[0]) > 0.9 /times) lu = false;
				if (fabs(ve[0]) > 0.9 /times) le = false;

				if (lu) {
					vu[1] = get_value(u, x+1, y, z, nz, ny, nx, times) - zero;
					vu[2] = get_value(u, x-1, y, z, nz, ny, nx, times) - zero;
					vu[3] = get_value(u, x, y+1, z, nz, ny, nx, times) - zero;
					vu[4] = get_value(u, x, y-1, z, nz, ny, nx, times) - zero;
					vu[5] = get_value(u, x, y, z+1, nz, ny, nx, times) - zero;
					vu[6] = get_value(u, x, y, z-1, nz, ny, nx, times) - zero;
					if (vu[1]*vu[2] > 0 && vu[3]*vu[4] > 0 && vu[5]*vu[6] > 0) lu = false;
				}

				if (le) {
					ve[1] = get_value(e, x+1, y, z, nz, ny, nx, times) - zero;
					ve[2] = get_value(e, x-1, y, z, nz, ny, nx, times) - zero;
					ve[3] = get_value(e, x, y+1, z, nz, ny, nx, times) - zero;
					ve[4] = get_value(e, x, y-1, z, nz, ny, nx, times) - zero;
					ve[5] = get_value(e, x, y, z+1, nz, ny, nx, times) - zero;
					ve[6] = get_value(e, x, y, z-1, nz, ny, nx, times) - zero;
					if (ve[1]*ve[2] > 0 && ve[3]*ve[4] > 0 && ve[5]*ve[6] > 0) le = false;
				}

				if (!lu && !le) continue;

				zp[nzp].x = x1;
				zp[nzp].y = y1;
				zp[nzp].z = wz-z1;
				if (lu) get_n(u, zp[nzp].n, x, y, z, nz, ny, nx, times);
				else get_n(e, zp[nzp].n, x, y, z, nz, ny, nx, times);
				zp[nzp].isedge = isedge;
				if (y > 0 && y > fabs(x)) zp[nzp].iscut = true;
				else zp[nzp].iscut = false;
				if (lu && le) zp[nzp].comp = 0;
				else if (lu) zp[nzp].comp = 1;
				else zp[nzp].comp = -1;

				nzp ++;
			}
		}
	}
}

double getangle(double *n, double x1, double y1, double z1)
{

	double x, y, z;
		x = x1;
		y = y1;
		z = z1;
		/*
	if (ROTATE) {
		x = R[0][0]*x1 + R[0][1]*y1 + R[0][2]*z1;
		y = R[1][0]*x1 + R[1][1]*y1 + R[1][2]*z1;
		z = R[2][0]*x1 + R[2][1]*y1 + R[2][2]*z1;
	}
	*/
	return fabs(3.1415926/2 - fabs(acos((n[0]*x+n[1]*y+n[2]*z)/sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])/sqrt(x*x+y*y+z*z)) - 3.1415926/2)) *2 / 3.1415927;
}


void draw_graph(unsigned char * buf, zpoint * zp, int nzp, int nz, int ny, int nx, double times)
{
	int wx, wy, wz;
	wx = int((nx-1)*times)/2;
	wy = int(((nx-1)*cos(RANGLE) + (ny*2-1)*sin(RANGLE))*times/2);
	wz = int(((nx-1)*sin(RANGLE) + (ny*2-1)*cos(RANGLE))*times/2);
	int width = wx*2;
	int height = wz*2;

	int i;

	for (i = 0; i < nzp; i++)
	{
		int r, g, b;
		r=g=b=0;
		double angle1;
		double angle2;
		double angle3;
		angle1 = getangle(zp[i].n, 0, 1, 0);
		angle2 = getangle(zp[i].n, 2, 4, 8);
		angle3 = getangle(zp[i].n, -5, 1, -0.1);
		
		//r += min(1.0 / (sqrt(double(zp[i].y * zp[i].y + zp[i].x*zp[i].x + (zp[i].z- wz)*(zp[i].z- wz)))/n/times) * 20, 200); 
		//r += int(angle1 * angle1 * 150);

		double tmp = (1+double(zp[i].y)/wy) ;
		
		r += int(angle1 * angle1 * 150 + (10+ pow(1- angle2, 4) * 150) * tmp); 
		g += int((30+ pow(1- angle2, 4) * 150) * tmp);
		b += int((30+ pow(1- angle2, 4) * 150) * tmp + (double(zp[i].x)/wx + 2.0)*60);
		
		r = min(r, 255);
		g = min(g, 255);
		b = min(b, 255);

		int offset;
		if (!zp[i].iscut) {
			offset = zp[i].z * width*2 + zp[i].x + width / 2;
			buf[offset*3] = r;
			buf[offset*3+1] = g;
			buf[offset*3+2] = b;

			if (zp[i].isedge) {
				buf[offset*3] = min(255, r+100);
				buf[offset*3+1] = min(255, g+100);
				buf[offset*3+2] = min(255, b+100);
			}
			/*
			int offset2 = offset + wh*width*2;
			buf[offset2*3] = buf[offset*3];
			buf[offset2*3+1] = buf[offset*3+1];
			buf[offset2*3+2] = buf[offset*3+2];
			*/
		}
		offset = zp[i].z * width*2 + zp[i].x + width / 2 + width;
		buf[offset*3] = r;
		buf[offset*3+1] = g;
		buf[offset*3+2] = b;
		/*
		int offset2 = offset + wh*width*2;
		buf[offset2*3] = buf[offset*3];
		buf[offset2*3+1] = buf[offset*3+1];
		buf[offset2*3+2] = buf[offset*3+2];
		*/
	}

	//draw a line in the middle of two graph.
	for (i = 0; i < height; i++) {
		int offset = i * width*2 + width;
		buf[offset*3] = 255;
		buf[offset*3+1] = 255;
		buf[offset*3+2] = 255;
	}
}

//graph with deep
void draw_graph1(unsigned char * buf, zpoint * zp, int nzp, int nz, int ny, int nx, double times, double redrate = 1.0, double bluerate = 0.0)  
{
	int wx, wy, wz;
	wx = int((nx-1)*times)/2;
	wy = int(((nx-1)*cos(RANGLE) + (ny*2-1)*sin(RANGLE))*times/2);
	wz = int(((nx-1)*sin(RANGLE) + (ny*2-1)*cos(RANGLE))*times/2);
	int width = wx*2;
	int height = wz*2;

	int i;

	for (i = 0; i < nzp; i++)
	{
		int r, g, b;
		r=g=b=0;
		double angle1;
		double angle2;
		double angle3;
		angle1 = getangle(zp[i].n, 0, 1, 0);
		angle2 = getangle(zp[i].n, 2, 4, 8);
		angle3 = getangle(zp[i].n, -5, 1, -0.1);
		
		//r += min(1.0 / (sqrt(double(zp[i].y * zp[i].y + zp[i].x*zp[i].x + (zp[i].z- wz)*(zp[i].z- wz)))/n/times) * 20, 200); 
		//r += int(angle1 * angle1 * 150);

//		double tmp = (1+double(zp[i].y)/wy) ;
		double tmp = 1.0;
		
//		r += int((30+ pow(1- angle2, 4) * 150) * tmp);
		r += int(angle1 * angle1 * 100 + (10+ pow(1- angle2, 4) * 150) * tmp); 
		g += int((30+ pow(1- angle2, 4) * 150) * tmp);
		b += int((30+ pow(1- angle2, 4) * 150) * tmp + 120);
		
		double rate = ((zp[i].y*cos(RANGLE) + zp[i].z*sin(RANGLE))/wx + 1.0)/2.0;
		r = int(r * rate);
		g = int(g * rate);
		b = int(b * rate);

		if (zp[i].comp == 1) {
			int tmp = r;
			r = b;
			b = tmp;
			//balance difference
			int r1 = r, b1 = b;
			b = b1 * redrate + r1*(1-redrate);
			r = r1 * redrate + b1*(1-redrate);
			//project
			tmp = sqrt((r1*r1 + b1*b1+1.0) / (r*r+b*b+1.0));
			r *= tmp;
			b *= tmp;
		}
		else if (zp[i].comp == -1) {
			//balance difference
			int r1 = r, b1 = b;
			r = b1 * bluerate + r1*(1-bluerate);
			b = r1 * bluerate + b1*(1-bluerate);
			//project
			tmp = sqrt((r1*r1 + b1*b1+1.0) / (r*r+b*b+1.0));
			r *= tmp;
			b *= tmp;
		}
		else if (zp[i].comp == 0) {
			r = min(r+50, 255);
			b = max(b-50, 0);
		}

		r = min(r, 255);
		g = min(g, 255);
		b = min(b, 255);

		int tt = r;
		r = b;
		b = tt;

		int offset;
		if (!zp[i].iscut) {
			offset = zp[i].z * width*2 + zp[i].x + width / 2;
			buf[offset*3] = r;
			buf[offset*3+1] = g;
			buf[offset*3+2] = b;

			if (zp[i].isedge) {
				buf[offset*3] = min(255, r+100);
				buf[offset*3+1] = min(255, g+100);
				buf[offset*3+2] = min(255, b+100);
			}
			/*
			int offset2 = offset + wh*width*2;
			buf[offset2*3] = buf[offset*3];
			buf[offset2*3+1] = buf[offset*3+1];
			buf[offset2*3+2] = buf[offset*3+2];
			*/
		}
		offset = zp[i].z * width*2 + zp[i].x + width / 2 + width;
		buf[offset*3] = r;
		buf[offset*3+1] = g;
		buf[offset*3+2] = b;
		/*
		int offset2 = offset + wh*width*2;
		buf[offset2*3] = buf[offset*3];
		buf[offset2*3+1] = buf[offset*3+1];
		buf[offset2*3+2] = buf[offset*3+2];
		*/
	}

	//draw a line in the middle of two graph.
	for (i = 0; i < height; i++) {
		int offset = i * width*2 + width;
		buf[offset*3] = 255;
		buf[offset*3+1] = 255;
		buf[offset*3+2] = 255;
	}
}

//name is the graph's name
//u is the graph data
//times is to the original size.
//phase is the component data
void fluidbase::output_graph(char *name, ARR3D u, ARR3D phase, double times, double zero, double rotateangle, double diff)
{
	RANGLE = rotateangle;
 
	char fname[100];
	sprintf(fname, "%s%s%dx%dx%d_%08d.jpg", folderpath, name, nx, ny, nz, roundx(time_current/TIME_STEP_OUTPUT));
#ifndef _DEBUG
	if (times <= 0.1) times = 400.0/nx;
#else
	if (times <= 0.1) times = 128.0/nx;
#endif
	int wx, wy, wz;
	wx = int((nx-1)*times)/2;
	wy = int(((nx-1)*cos(RANGLE) + (ny*2-1)*sin(RANGLE))*times/2);
	wz = int(((nx-1)*sin(RANGLE) + (ny*2-1)*cos(RANGLE))*times/2);
	int m_width = wx*2;
	int m_height = wz*2;

	unsigned char * m_buf = new unsigned char [6*m_width*m_height];
	memset(m_buf, 0, 6*m_width*m_height);

	zpoint * zp = new zpoint [m_width*m_height * 15];
	int nzp = 0;

	ARR3D v;
	malloc_array(v);
	/*avg_mesh(u, v);*/ copy_array_s(v, u);
//------------------------------------------
#ifdef _TESTDRAWGRAPH_1_
	get_zero_points(v, phase, nz, ny, nx, times, zp, nzp, zero);
	FILE *f = fopen("f:\\temp\\zp.dat", "w+b");
	fwrite(&nzp, 1, sizeof(int), f);
	fwrite(zp, nzp, sizeof(zpoint), f);
	fclose(f);
#else
	#ifdef _TESTDRAWGRAPH_
		FILE *f = fopen("f:\\temp\\zp.dat", "rb");
		fread(&nzp, 1, sizeof(int), f);
		fread(zp, nzp, sizeof(zpoint), f);
		fclose(f);
	#else
		get_zero_points(v, phase, nz, ny, nx, times, zp, nzp, zero);
	#endif
#endif
//------------------------------------------
	free_array(v);

	draw_graph1(m_buf, zp, nzp, nz, ny, nx, times);

	//JpegFile jf;
	//jf.RGBToJpegFile(fname, m_buf, m_width*2, m_height, true, 100);
	stbi_write_jpg(fname, m_width * 2, m_height, 3, m_buf, 100);

	delete[] zp;
	delete[] m_buf;
}

//name is the graph's name
//u is the graph data
//times is to the original size.
//phase is the component data
void fluidbase::output_graph_two_bubble(char *name, ARR3D u, ARR3D e, double times, double zero, double rotateangle, double diff)
{
	RANGLE = rotateangle;
 
	char fname[100];
	sprintf(fname, "%s%s%dx%dx%d_%08d.jpg", folderpath, name, nx, ny, nz, roundx(time_current/TIME_STEP_OUTPUT));
#ifndef _DEBUG
	if (times <= 0.1) times = 400.0/nx;
#else
	if (times <= 0.1) times = 128.0/nx;
#endif
	int wx, wy, wz;
	wx = int((nx-1)*times)/2;
	wy = int(((nx-1)*cos(RANGLE) + (ny*2-1)*sin(RANGLE))*times/2);
	wz = int(((nx-1)*sin(RANGLE) + (ny*2-1)*cos(RANGLE))*times/2);
	int m_width = wx*2;
	int m_height = wz*2;

	unsigned char * m_buf = new unsigned char [6*m_width*m_height];
	memset(m_buf, 0, 6*m_width*m_height);

	zpoint * zp = new zpoint [m_width*m_height * 15];
	int nzp = 0;

	ARR3D v, w;
	malloc_array(v);
	malloc_array(w);
	/*avg_mesh(u, v);*/ copy_array_s(v, u);
	copy_array_s(w, e);
//------------------------------------------
	get_zero_points_two_bubble(v, w, NULL, nz, ny, nx, times, zp, nzp, zero);
//------------------------------------------
	free_array(v);

	draw_graph1(m_buf, zp, nzp, nz, ny, nx, times);

//	JpegFile jf;
//	jf.RGBToJpegFile(fname, m_buf, m_width*2, m_height, true, 100);
	stbi_write_jpg(fname, m_width * 2, m_height, 3, m_buf, 100);

	delete[] zp;
	delete[] m_buf;
}

	
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/
/********************************************************************************************************************/

void get_around(double ** data, double d[4][4], double x, double y, double z, int m, int n, double times)
{
	int i, j;
	double tempr = sqrt(x*x + y*y) / times + 0.5*(n-1);
	double tempz = z / times;
	while (tempz < 0) tempz += m;
	while (tempz >= m) tempz -= m;
	int ir = int(tempr);
	int iz = int(tempz);
				
	for (i = 0; i < 4; i++ )
		for (j = 0; j < 4; j++)
			d[i][j] = data[((iz+ j -1)+m)%m][(ir + i -1 + n)%n];
}


double N[16][16] = {
0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,1.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,
-1.1102230e-016,0.0000000e+000,0.0000000e+000,5.5511151e-017,-3.3333333e-001,-5.0000000e-001,1.0000000e+000,-1.6666667e-001,1.1102230e-016,0.0000000e+000,0.0000000e+000,-5.5511151e-017,-5.5511151e-017,0.0000000e+000,0.0000000e+000,2.7755576e-017,
0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,5.0000000e-001,-1.0000000e+000,5.0000000e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,
0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,-1.6666667e-001,5.0000000e-001,-5.0000000e-001,1.6666667e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,-2.7755576e-017,0.0000000e+000,0.0000000e+000,0.0000000e+000,
-1.1102230e-016,-3.3333333e-001,0.0000000e+000,-2.7755576e-017,0.0000000e+000,-5.0000000e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,1.0000000e+000,0.0000000e+000,0.0000000e+000,5.5511151e-017,-1.6666667e-001,0.0000000e+000,0.0000000e+000,
1.1111111e-001,1.6666667e-001,-3.3333333e-001,5.5555556e-002,1.6666667e-001,2.5000000e-001,-5.0000000e-001,8.3333333e-002,-3.3333333e-001,-5.0000000e-001,1.0000000e+000,-1.6666667e-001,5.5555556e-002,8.3333333e-002,-1.6666667e-001,2.7777778e-002,
-1.6666667e-001,3.3333333e-001,-1.6666667e-001,-5.5511151e-017,-2.5000000e-001,5.0000000e-001,-2.5000000e-001,0.0000000e+000,5.0000000e-001,-1.0000000e+000,5.0000000e-001,0.0000000e+000,-8.3333333e-002,1.6666667e-001,-8.3333333e-002,2.7755576e-017,
5.5555556e-002,-1.6666667e-001,1.6666667e-001,-5.5555556e-002,8.3333333e-002,-2.5000000e-001,2.5000000e-001,-8.3333333e-002,-1.6666667e-001,5.0000000e-001,-5.0000000e-001,1.6666667e-001,2.7777778e-002,-8.3333333e-002,8.3333333e-002,-2.7777778e-002,
0.0000000e+000,5.0000000e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,-1.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,5.0000000e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,
-1.6666667e-001,-2.5000000e-001,5.0000000e-001,-8.3333333e-002,3.3333333e-001,5.0000000e-001,-1.0000000e+000,1.6666667e-001,-1.6666667e-001,-2.5000000e-001,5.0000000e-001,-8.3333333e-002,-2.7755576e-017,0.0000000e+000,0.0000000e+000,-2.7755576e-017,
2.5000000e-001,-5.0000000e-001,2.5000000e-001,0.0000000e+000,-5.0000000e-001,1.0000000e+000,-5.0000000e-001,0.0000000e+000,2.5000000e-001,-5.0000000e-001,2.5000000e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,0.0000000e+000,
-8.3333333e-002,2.5000000e-001,-2.5000000e-001,8.3333333e-002,1.6666667e-001,-5.0000000e-001,5.0000000e-001,-1.6666667e-001,-8.3333333e-002,2.5000000e-001,-2.5000000e-001,8.3333333e-002,-2.7755576e-017,0.0000000e+000,0.0000000e+000,-6.9388939e-018,
0.0000000e+000,-1.6666667e-001,0.0000000e+000,-2.7755576e-017,0.0000000e+000,5.0000000e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,-5.0000000e-001,0.0000000e+000,0.0000000e+000,0.0000000e+000,1.6666667e-001,0.0000000e+000,0.0000000e+000,
5.5555556e-002,8.3333333e-002,-1.6666667e-001,2.7777778e-002,-1.6666667e-001,-2.5000000e-001,5.0000000e-001,-8.3333333e-002,1.6666667e-001,2.5000000e-001,-5.0000000e-001,8.3333333e-002,-5.5555556e-002,-8.3333333e-002,1.6666667e-001,-2.7777778e-002,
-8.3333333e-002,1.6666667e-001,-8.3333333e-002,-2.7755576e-017,2.5000000e-001,-5.0000000e-001,2.5000000e-001,0.0000000e+000,-2.5000000e-001,5.0000000e-001,-2.5000000e-001,0.0000000e+000,8.3333333e-002,-1.6666667e-001,8.3333333e-002,0.0000000e+000,
2.7777778e-002,-8.3333333e-002,8.3333333e-002,-2.7777778e-002,-8.3333333e-002,2.5000000e-001,-2.5000000e-001,8.3333333e-002,8.3333333e-002,-2.5000000e-001,2.5000000e-001,-8.3333333e-002,-2.7777778e-002,8.3333333e-002,-8.3333333e-002,2.7777778e-002
};


double get_value(double ** data, double x, double y, double z, int m, int n, int times)
{
	double d[4][4];
	get_around(data, d, x, y, z, m, n , times);

	double tempr = sqrt(x*x + y*y) / times + 0.5*(n-1);
	double tempz = z / times;
	while (tempz < 0) tempz += m;
	while (tempz >= m) tempz -= m;
	tempr -= int(tempr); 
	tempz -= int(tempz); 

	return d[1][1]*(1-tempr)*(1-tempz) + d[1][2]*(1-tempr)*tempz + d[2][1]*tempr*(1-tempz) + d[2][2]*tempr*tempz;

	double coef[4][4];
	int i, j, k, l;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			coef[i][j] = 0;
			for(k = 0; k < 4; k++)
				for (l = 0; l < 4; l++) 
				{
					coef[i][j] += N[i*4+j][k*4+l] * d[k][l];
				}
		}
			
	double ret = 0;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			double tmp = coef[i][j];
			int ii;
			for (ii = 0; ii < i; ii++) tmp *= tempr;
			for (ii = 0; ii < j; ii++) tmp *= tempz;
			ret += tmp;
		}

	return ret;
}

double get_n(double ** data, double x, double y, double z, int m, int n, int times)
{
	double d[4][4];
	get_around(data, d, x, y, z, m, n , times);

	double tempr = sqrt(x*x + y*y) / times + 0.5*(n-1);
	double tempz = z / times;
	while (tempz < 0) tempz += m;
	while (tempz >= m) tempz -= m;
	tempr -= int(tempr); 
	tempz -= int(tempz); 

	double coef[4][4];
	int i, j, k, l;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			coef[i][j] = 0;
			for(k = 0; k < 4; k++)
				for (l = 0; l < 4; l++) 
				{
					coef[i][j] += N[i*4+j][k*4+l] * d[k][l];
				}
		}
			
	double ret1 = 0;
	for (i = 1; i < 4; i++)
		for (j = 0; j < 4; j++) {
			ret1 += coef[i][j] * i*pow(tempr,(i-1)) * pow(tempz,j);
		}

	double ret2 = 0;
	for (i = 0; i < 4; i++)
		for (j = 1; j < 4; j++) {
			ret2 += coef[i][j] * pow(tempr,i) * j * pow(tempz,(j-1));
		}

	double ret;
	if (ret1 == 0) ret = 1e10;
	else ret = ret2/ret1;
	return ret;
}

/*
double R[3][3] = {    0.5000,         0,    0.8660,
    0.7500,    0.5000,   -0.4330,
	-0.4330,    0.8660,    0.2500};
*/

void get_zero_points(double ** u, double ** phase, int ny, int nx, double times, zpoint * zp, int & nzp, double zero)
{
	int nz = ny;
	ny = nx;

	int i, j;
	int x1, y1, z1;
	double x, y, z;

	double R[3][3] = {	1.0,         0.0,    0.0,
					0.000,    cos(RANGLE),		-sin(RANGLE),
					0.0,    +sin(RANGLE),    cos(RANGLE)};

	int wx, wy, wz;
	wx = int((nx-1)*times)/2;
	wy = int(((ny-1)*cos(RANGLE) + (nz*2-1)*sin(RANGLE))*times/2);
	wz = int(((ny-1)*sin(RANGLE) + (nz*2-1)*cos(RANGLE))*times/2);

	double rx, ry, rz;
	rx = (nx*times)/2;
	ry = (ny*times)/2;
	rz = nz*times;

	for (y1 = -wy; y1 <= wy-1; y1++) {
		for (x1 = -wx; x1 <= wx-1; x1++) {
			for (z1 = -wz; z1 <= wz-1; z1++)
			{
				bool isedge = false;
				x = R[0][0]*x1 + R[0][1]*y1 + R[0][2]*z1;
				if (x < -rx || x >= rx-times) continue;
				y = R[1][0]*x1 + R[1][1]*y1 + R[1][2]*z1;
				if (y < -ry || y >= ry-times) continue;
				if (x*x + y*y >= rx*rx) continue;
				z = R[2][0]*x1 + R[2][1]*y1 + R[2][2]*z1;
				if (z/times < -nz+2 || z/times > nz-1 ) continue;

				if (y < fabs(x) && y+1 > fabs(x)) isedge = true;

				//if (z < 0) z += rz;
				//if (fabs(u[iz][iy][ix]) > 0.9) continue;

				double d[4][4];
				get_around(u, d, x, y, z, nz, nx, times);

				bool l = false;
				for (i = 0; i < 4; i ++) for (j = 0; j < 4; j++) if (!l && (d[i][j]-zero)*(d[1][1]-zero) < 0) l = true;
				if (!l) continue;

				double v[7];
				v[0] = get_value(u, x, y, z, nz, nx, times) - zero;
				if (fabs(v[0] - zero) > 0.5/times) continue;

				v[1] = get_value(u, x+1, y, z, nz, nx, times) - zero;
				v[2] = get_value(u, x-1, y, z, nz, nx, times) - zero;
				v[3] = get_value(u, x, y+1, z, nz, nx, times) - zero;
				v[4] = get_value(u, x, y-1, z, nz, nx, times) - zero;
				v[5] = get_value(u, x, y, z+1, nz, nx, times) - zero;
				v[6] = get_value(u, x, y, z-1, nz, nx, times) - zero;
				if (v[1]*v[2] > 0 && v[3]*v[4] > 0 && v[5]*v[6] > 0) continue;

				zp[nzp].x = x1;
				zp[nzp].y = y1;
				zp[nzp].z = wz-z1-1;
				zp[nzp].n[0] = x / sqrt(double(x*x + y*y));
				zp[nzp].n[1] = y / sqrt(double(x*x + y*y));
				zp[nzp].n[2] = get_n(u, x, y, z, nz, nx, times);
				zp[nzp].isedge = isedge;
				if (y > 0 && y > fabs(x)) zp[nzp].iscut = true;
				else zp[nzp].iscut = false;
				if (phase == NULL) zp[nzp].comp = -1;
				else {
					double v = get_value(phase, x, y, z, nz, nx, times);
					//if (fabs(v) > 0.9) continue;
					if (v > 0.0) zp[nzp].comp = 1;
//					if (phase[iz][iy][ix] > 0.0) zp[nzp].comp = 1;
					else zp[nzp].comp = -1;
				}
				nzp ++;
			}
		}
	}
}

