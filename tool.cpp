#include "tool.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void copyfile(char* filefrom, char* fileto)
{
	FILE* f1 = fopen(filefrom, "rb");
	FILE* f2 = fopen(fileto, "w+b");
	fseek(f1, 0, SEEK_END);
	int len = ftell(f1);
	char* buf = new char[len];
	fseek(f1, 0, SEEK_SET);
	fread(buf, len, 1, f1);
	fclose(f1);
	fwrite(buf, len, 1, f2);
	fclose(f2);
	delete[] buf;
}


int roundx(double x)
{
	return int(floor(x + 0.5));
}

void write_array(ARR3D d, char* fname)
{
	int dnz = ((int*)d)[-1];
	int dny = ((int*)d)[-2];
	int dnx = ((int*)d)[-3];
	FILE* f = fopen(fname, "w+b");

	int i, j;
	for (i = 0; i < dnz; i++) {
		for (j = 0; j < dny; j++) {
			fwrite(d[i][j], sizeof(double), dnx, f);
		}
	}
	fclose(f);
}

void read_array(ARR3D d, char* fname)
{
	int dnz = ((int*)d)[-1];
	int dny = ((int*)d)[-2];
	int dnx = ((int*)d)[-3];
	FILE* f = fopen(fname, "rb");
	if (!f) {
		printf("Can not load file %s! set as zeros\n", fname);
		zero_array(d);
		return;
	}

	int i, j;
	for (i = 0; i < dnz; i++) {
		for (j = 0; j < dny; j++) {
			fread(d[i][j], sizeof(double), dnx, f);
		}
	}
	fclose(f);
}

//every where can use this function except in single openmp thread.
void zero_array(ARR3D  x)
{
	int nz = ((int*)x)[-1];
	int ny = ((int*)x)[-2];
	int nx = ((int*)x)[-3];
	int binpal = 0;
#ifdef _OPENMP
	binpal = omp_in_parallel();
#endif
	int i, j, k;
	if (binpal) {
#pragma omp for
		for (i = 0; i < nz; i++) {
			for (j = 0; j < ny; j++) {
				double* xij = x[i][j];
				for (k = 0; k < nx; k++) {
					xij[k] = 0.0;
				}
			}
		}
	}
	else {
		memset(x[0][0], 0, (x[1][0] - x[0][0]) * nz * sizeof(double));
	}
}

void copy_array(ARR3D d, ARR3D s)
{
	int dnz = ((int*)d)[-1];
	int dny = ((int*)d)[-2];
	int dnx = ((int*)d)[-3];
	int snz = ((int*)s)[-1];
	int sny = ((int*)s)[-2];
	int snx = ((int*)s)[-3];
	if (dnx != snx || dny != sny || dnz != snz) {
		printf("error in copy_array 3d!\n");
		exit(1);
	}
	int binpal = 0;
#if defined (_OPENMP)
	binpal = omp_in_parallel();
#endif
	int i, j, k;
	if (binpal) {
#pragma omp for
		for (i = 0; i < dnz; i++) {
			for (j = 0; j < dny; j++) {
				double* dij = d[i][j];
				double* sij = s[i][j];
				for (k = 0; k < dnx; k++) {
					dij[k] = sij[k];
				}
			}
		}
	}
	else {
		for (i = 0; i < dnz; i++) {
			for (j = 0; j < dny; j++) {
				double* dij = d[i][j];
				double* sij = s[i][j];
				for (k = 0; k < dnx; k++) {
					dij[k] = sij[k];
				}
			}
		}
	}
}

//serial one, special use.
//every where can use this function.
void copy_array_s(ARR3D d, ARR3D s)
{
	int dnz = ((int*)d)[-1];
	int dny = ((int*)d)[-2];
	int dnx = ((int*)d)[-3];
	int snz = ((int*)s)[-1];
	int sny = ((int*)s)[-2];
	int snx = ((int*)s)[-3];
	if (dnx != snx || dny != sny || dnz != snz) {
		printf("error in copy_array 3d!\n");
		exit(1);
	}
	int i, j, k;
	for (i = 0; i < dnz; i++) {
		for (j = 0; j < dny; j++) {
			double* dij = d[i][j];
			double* sij = s[i][j];
			for (k = 0; k < dnx; k++) {
				dij[k] = sij[k];
			}
		}
	}
}
