#include "tool.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

bool malloc_array(ARR3D  &x, int nx, int ny, int nz)
{
	int count = (nx+CACHELINE/4)*ny*nz + (ny+CACHELINE)*nz + nz + ARSUF;
	double * buf = new double [count];
	buf += ARSUF;
	((int*)buf)[-1] = nz;
	((int*)buf)[-2] = ny;
	((int*)buf)[-3] = nx;
	x = (ARR3D)(buf);
	buf = buf + (ny+CACHELINE)*nz + nz;

	int i, j;
	for (i = 0; i < nz; i++) {
		x[i] = (double **)(x + nz + (ny+CACHELINE) * i);
	}
	for (i = 0; i < nz; i++) {
		for (j = 0; j < ny; j++) {
			x[i][j] = buf + (nx+CACHELINE/4)*ny*i + (nx+CACHELINE/4)*j;
		}
	}
	return true;
}

//every where can use this function.
void free_array(ARR3D& x)
{
	if (x) delete[]((double*)x - ARSUF);
	x = NULL;
}


__global__ void config_array(double* buf, int nx, int ny, int nz)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	if (index == 0) {
		((int*)buf)[-1] = nz;
		((int*)buf)[-2] = ny;
		((int*)buf)[-3] = nx;
	}
	ARR3D x = (ARR3D)buf;
	buf = buf + (ny + CACHELINE) * nz + nz;

	if (index < nz) {
		int i = index;
		x[i] = (double**)(x + nz + (ny + CACHELINE) * i);
	}

	if (index < ny * nz) {
		int i = index / ny;
		int j = index % ny;
		double ** p = (double**)(x + nz + (ny + CACHELINE) * i);
		p[j] = buf + (nx + CACHELINE / 4) * ny * i + (nx + CACHELINE / 4) * j;
	}
}

bool cuda_malloc_array(ARR3D& x, int nx, int ny, int nz)
{
	int count = (nx + CACHELINE / 4) * ny * nz + (ny + CACHELINE) * nz + nz + ARSUF;
	double* buf;
	cudaMalloc((void **)&buf, count*sizeof(double));
	buf += ARSUF;
	config_array<<<(ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>> (buf, nx, ny, nz);
	x = (ARR3D)buf;
	return true;
}

bool cuda_copyin_array(ARR3D hx, ARR3D dx, int nx, int ny, int nz)
{
	double* hbuf = (double*)hx;
	double* dbuf = (double*)dx;
	hbuf += (ny + CACHELINE) * nz + nz;
	dbuf += (ny + CACHELINE) * nz + nz;
	cudaMemcpy(dbuf, hbuf, (nx + CACHELINE / 4) * ny * nz * sizeof(double), cudaMemcpyHostToDevice);
	return true;
}

bool cuda_copyout_array(ARR3D hx, ARR3D dx, int nx, int ny, int nz)
{
	double* hbuf = (double*)hx;
	double* dbuf = (double*)dx;
	hbuf += (ny + CACHELINE) * nz + nz;
	dbuf += (ny + CACHELINE) * nz + nz;
	cudaMemcpy(hbuf, dbuf, (nx + CACHELINE / 4) * ny * nz * sizeof(double), cudaMemcpyDeviceToHost);
	return true;
}

