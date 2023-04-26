#include "Calc.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <cufftw.h>

Calc::Calc(int nx1, int nz1, double lenx1)
{
    nx = nx1;
    ny = nx1;
    nz = nz1;
    xlen = lenx1; ylen = lenx1;
    zlen = lenx1;
    //	zlen = lenx1*nz1/nx1; //such that a average mesh.

    if (nx1 % 2 == 1) {
        printf("The grid nx = %d should be a even number!\n", nx1);
        exit(1);
    }
    // CUFFT plan
    cufftPlan3d(&planf, nz, ny, nx, CUFFT_D2Z);
    cufftPlan3d(&planb, nz, ny, nx, CUFFT_Z2D);
    cudaMalloc(&ctmp1, (nx / 2 + 1) * ny * nz * sizeof(cufftDoubleComplex));
    cudaMalloc(&ctmp2, (nx / 2 + 1) * ny * nz * sizeof(cufftDoubleComplex));
    cudaMalloc(&ctmp3, (nx / 2 + 1) * ny * nz * sizeof(cufftDoubleComplex));
    //cudaMalloc(&fftbackcopy, (nx / 2 + 1) * ny * nz * sizeof(cufftDoubleComplex));
    
}

Calc::~Calc()
{
    cufftDestroy(planf);
    cufftDestroy(planb);

    cudaFree(ctmp1);
    cudaFree(ctmp2);
    cudaFree(ctmp3);
    //cudaFree(fftbackcopy);
}

__global__ void normalize_kernel(cufftDoubleComplex *s, double factor, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        s[tid].x *= factor;
        s[tid].y *= factor;
    }
}


void Calc::r2c_3d_f(DARR3D u, cufftDoubleComplex* s)
{
    int N = (nx/2+1) * ny * nz;
    cufftExecD2Z(planf, u, s);
    normalize_kernel << <(N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (s, 1.0/(nx*ny*nz), N);   //normalized
}

//the s data will be changed, I commented the version with backup
void Calc::c2r_3d_b(cufftDoubleComplex* s, DARR3D u)
{
    //cudaMemcpy(fftbackcopy, s, (nx / 2 + 1) * ny * nz * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToDevice);
    //cufftExecZ2D(planb, fftbackcopy, u); //the s data will be changed, so we need to make copy of it.
    cufftExecZ2D(planb, s, u);
}


__global__ void laplace_kernel(cufftDoubleComplex* u, int nx, int ny, int nz, double cx, double cy, double cz, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        int z = tid / (ny * (nx / 2 + 1));
        int y = (tid % (ny * (nx / 2 + 1))) / (nx / 2 + 1);
        int x = (tid % (ny * (nx / 2 + 1))) % (nx / 2 + 1);
        z = z < nz - z ? z : nz - z;
        y = y < ny - y ? y : ny - y;
        x = x < nx - x ? x : nx - x;
        double tmp = -(x * x *cx+ y * y*cy + z * z*cz); 
        u[tid].x *= tmp;
        u[tid].y *= tmp;
    }
}


void Calc::laplace(cufftDoubleComplex* u)
{
    // Multiply the coefficients together and normalize the result
    double cx = 2.0 * PI / (xlen);
    double cy = 2.0 * PI / (ylen);
    double cz = 2.0 * PI / (zlen);
    cx *= cx; cy *= cy; cz *= cz;

    laplace_kernel << <(nz * ny * (nx / 2 + 1) + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (u, nx, ny, nz, cx, cy, cz, nz * ny * (nx / 2 + 1));
    cudaDeviceSynchronize();
}

void Calc::lap(cufftDoubleComplex* fd, DARR3D out)
{
    cudaMemcpy(ctmp1, fd, nz * ny * (nx / 2 + 1) * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToDevice);
    laplace(ctmp1);
    c2r_3d_b(ctmp1, out);
}

void Calc::lap(DARR3D d, DARR3D out)
{
    r2c_3d_f(d, ctmp1);
    laplace(ctmp1);
    c2r_3d_b(ctmp1, out);
}

__global__ void f_kernel(DARR3D f, DARR3D lap, DARR3D d, double eps, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        f[tid] = -lap[tid] * eps + (d[tid] * d[tid] - 1.0) * d[tid] / eps;
    }
}

//calculate f from lap, d, eps.
void Calc::f(DARR3D f, DARR3D lap, DARR3D d, double eps)
{
    f_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (f, lap, d, eps, nz * ny * nx);
    cudaDeviceSynchronize();
}

__global__ void g0_kernel(DARR3D g0, DARR3D d, DARR3D f, double eps, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        g0[tid] = -g0[tid] + (3.0 * d[tid] * d[tid] - 1.0) * f[tid] / (eps * eps);
    }
}

//m_g0 = lap(f);
//m_g0[i][j][k] = -m_g0[i][j][k] + (3.0*m_d[i][j][k]*m_d[i][j][k] - 1.0)*m_f[i][j][k] / (eps*eps);
void Calc::g0(DARR3D g0, DARR3D d, DARR3D f, double eps)
{
    g0_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (g0, d, f, eps, nz * ny * nx);
    cudaDeviceSynchronize();
}


__global__ void sumArray_kernel(double* a, int N, double* result, double* val, double addition, double multiply)
{
    // Allocate shared memory for the block
    __shared__ double sdata[THREADS_PER_BLOCK];

    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        sdata[threadIdx.x] = a[tid];
    }
    else {
        sdata[threadIdx.x] = 0.0;
    }

    // Synchronize threads to ensure all data has been copied
    __syncthreads();

    // Reduce within the block
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            sdata[threadIdx.x] += sdata[threadIdx.x + s];
        }
        __syncthreads();
    }

    // Write the result back to global memory
    if (threadIdx.x == 0) {
        result[blockIdx.x] = sdata[0];
        if (tid == 0) {
            *val = (result[0]+addition) * multiply;
        }
    }
}

__global__ void intf_kernel(double * u, int N, double* ret, double* val, double addition, double multiply)
{
    // Allocate shared memory for the block
    __shared__ double sdata[THREADS_PER_BLOCK];

    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        sdata[threadIdx.x] = u[tid];
    }
    else {
        sdata[threadIdx.x] = 0.0;
    }

    // Synchronize threads to ensure all data has been copied
    __syncthreads();

    // Reduce within the block
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            sdata[threadIdx.x] += sdata[threadIdx.x + s];
        }
        __syncthreads();
    }

    // Write the result back to global memory
    if (threadIdx.x == 0) {
        ret[blockIdx.x] = sdata[0]/N;  //averaged!
        if (tid == 0) {
            *val = (ret[0] + addition) * multiply;
        }
    }
}

//average all and put result in ret1[0], ret2 is temporay use.
void Calc::intf(DARR3D u, DARR3D ret1, DARR3D ret2, double * val, double addition, double multiply)
{
    int n = (nx*ny*nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    intf_kernel << <n, THREADS_PER_BLOCK >> > (u, nx * ny * nz, ret1, val, addition, multiply * xlen * ylen * zlen);
    cudaDeviceSynchronize();

    while (n > 1) {
        DARR3D t;
        t = ret1; ret1 = ret2; ret2 = t;
        int m = n;
        n = (n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        sumArray_kernel << <n, THREADS_PER_BLOCK >> > (ret2, m, ret1, val, addition, multiply * xlen * ylen * zlen);
    }

    cudaDeviceSynchronize();
}

__global__ void nld_kernel1(DARR3D d, DARR3D tmp1, double kappa, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid < N) {
        tmp1[tid] = d[tid] * (d[tid] * d[tid] - 1.0 - kappa);
    }
}

__global__ void nld_kernel2(DARR3D nld, DARR3D tmp1, DARR3D d, DARR3D f, double eps, double D, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    double eps2 = eps * eps;
    double eps3 = eps * eps2;

    if (tid < N) {
        nld[tid] = -nld[tid] / eps + D * tmp1[tid] / eps3 + (3 * d[tid] * d[tid] - D - 1) / eps2 * f[tid];
    }
}
void Calc::nld(DARR3D nldp, DARR3D d, DARR3D f, DARR3D tmp1, double eps, double kappa, double D)
{

    nld_kernel1 << <(nx * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (d, tmp1, kappa, nx*ny*nz);
    cudaDeviceSynchronize();

    lap(tmp1, nldp);

    nld_kernel2 << <(nx * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (nldp, tmp1, d, f, eps, D, nx * ny * nz);
    cudaDeviceSynchronize();
}

//here the factor is actually for x
__global__ void gradient_kernel(cufftDoubleComplex* u, cufftDoubleComplex* u1, cufftDoubleComplex* u2, cufftDoubleComplex* u3, int nx, int ny, int nz, double fac1, double fac2, double fac3)
{
    int N = nz * ny * (nx / 2 + 1);

    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        
        int x = tid % (nx / 2 + 1);
        if (x == nx / 2) u1[tid].x = u1[tid].y = 0;
        else {
            double cx = u[tid].x * fac1;
            double cy = u[tid].y * fac1;
            u1[tid].x = -cy * x;
            u1[tid].y = cx * x;
        }

        int y = (tid % (ny * (nx / 2 + 1))) / (nx / 2 + 1);
        if (y == ny / 2) {
            u2[tid].x = u2[tid].y = 0;
        }
        else {
            y = y < ny - y ? y : y - ny;
            double cx = u[tid].x * fac2;
            double cy = u[tid].y * fac2;
            u2[tid].x = -cy * y;
            u2[tid].y = cx * y;
        }

        int z = tid / (ny * (nx / 2 + 1));
        if (z == nz / 2) u3[tid].x = u3[tid].y = 0;
        else {
            z = z < nz - z ? z : z - nz;
            double cx = u[tid].x * fac3;
            double cy = u[tid].y * fac3;
            u3[tid].x = -cy * z;
            u3[tid].y = cx * z;
        }
    }
}

//we need xlen = ylen = zlen to use here, otherwise, we need change the code...
void Calc::fftgradient(cufftDoubleComplex* uf, DARR3D udx, DARR3D udy, DARR3D udz)
{
    gradient_kernel << <(nz * ny * (nx / 2 + 1) + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (uf, ctmp1, ctmp2, ctmp3, nx, ny, nz, 2.0 * PI / xlen, 2.0 * PI / ylen, 2.0 * PI / zlen); 
    cudaDeviceSynchronize();
    // Transform signal back
    cufftExecZ2D(planb, ctmp1, udx);
    cufftExecZ2D(planb, ctmp2, udy);
    cufftExecZ2D(planb, ctmp3, udz);

    cudaDeviceSynchronize();

//    gradientsq_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (lap, udx, udy, udz, nx, ny, nz);
}

__global__ void s_term_kernel(double* s, double * d, double* dx, double* dy, double* dz, double eps, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        s[tid] = 0.5 * eps * (dx[tid] * dx[tid] + dy[tid] * dy[tid] + dz[tid] * dz[tid])
            + 0.25 / eps * (d[tid] * d[tid] - 1) * (d[tid] * d[tid] - 1);
    }
}

void Calc::s_term(DARR3D s, DARR3D d, DARR3D dx, DARR3D dy, DARR3D dz, double eps)
{
    s_term_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (s, d, dx, dy, dz, eps, nx*ny*nz);
    cudaDeviceSynchronize();
}

__global__ void multi_kernel(double* a, double* b, double* c, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        c[tid] = a[tid]*b[tid];
    }
}

void Calc::multi(DARR3D a, DARR3D b, DARR3D c)
{
    multi_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (a, b, c, nx * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void plus_kernel(DARR3D u, DARR3D nld1, double t, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        u[tid] += t * nld1[tid];
    }
}

void Calc::plus(DARR3D u, DARR3D nld1, double t)
{
    plus_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (u, nld1, t, nx * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void nld1_kernel(DARR3D nld, DARR3D f, DARR3D g, double lamb1, double lamb2, double eta, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        nld[tid] = lamb1 * 0.5 + lamb2 * f[tid] + g[tid] * eta;
    }
}

//nld = lamb1 * 0.5 + lamb2 * f[i][j][k] + g[i][j][k] * eta;
void Calc::nld1(DARR3D nld, DARR3D f, DARR3D g, double lamb1, double lamb2, double eta)
{
    nld1_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (nld, f, g, lamb1, lamb2, eta, nx * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void nld2_kernel(DARR3D nld, DARR3D f, DARR3D g, double lamb1, double lamb2, double eta, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        nld[tid] = -lamb1 * 0.5 - lamb2 * f[tid] * 3.0 / (2.0 * sqrt(2.0)) - g[tid] * eta;
    }
}

//-lamb1 * 0.5 - lamb2 * f * 3.0 / (2.0 * sqrt(2.0)) - g * eta;
void Calc::nld2(DARR3D nld, DARR3D f, DARR3D g, double lamb1, double lamb2, double eta)
{
    nld2_kernel << <(nz * ny * nx + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (nld, f, g, lamb1, lamb2, eta, nx * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void etd1_kernel(cufftDoubleComplex* a, cufftDoubleComplex* fd, cufftDoubleComplex* nld, DARR3D Het, double time_step, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        double c = exp(Het[tid] * time_step);
        a[tid].x = c * fd[tid].x + (c - 1.0) * nld[tid].x / Het[tid];
        a[tid].y = c * fd[tid].y + (c - 1.0) * nld[tid].y / Het[tid];
    }
}

//a[tid] = c * fd[tid] + (c - 1.0) * nld[tid] / Het[tid];
void Calc::etd1(cufftDoubleComplex* d, cufftDoubleComplex* fd, cufftDoubleComplex* nld, DARR3D Het, double time_step)
{
    etd1_kernel << <((nx / 2 + 1) * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (d, fd, nld, Het, time_step, (nx / 2 + 1) * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void etd2_kernel(cufftDoubleComplex* a, cufftDoubleComplex* fd, cufftDoubleComplex* nnld, cufftDoubleComplex* nld, DARR3D Het, double time_step, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        double c = exp(Het[tid] * time_step);
        if (time_step * Het[tid] == 0) c = time_step/2;
        else c = (c - 1.0 - time_step * Het[tid]) / (time_step * Het[tid] * Het[tid]);
        a[tid].x = fd[tid].x + c * (nnld[tid].x - nld[tid].x);
        a[tid].y = fd[tid].y + c * (nnld[tid].y - nld[tid].y);
    }
}

//m_d[tid] = fd[tid] + (c - 1.0 - time_step * Het[tid]) / (time_step * Het[tid] * Het[tid]) * (nnld[tid] - nld[tid]);
void Calc::etd2(cufftDoubleComplex* d, cufftDoubleComplex* fd, cufftDoubleComplex* nnld, cufftDoubleComplex* nld, DARR3D Het, double time_step)
{
    etd2_kernel << <((nx / 2 + 1) * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (d, fd, nnld, nld, Het, time_step, (nx / 2 + 1) * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void etd41_kernel(cufftDoubleComplex* a, cufftDoubleComplex* fd, cufftDoubleComplex* nld, DARR3D Het, double time_step, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        double e1 = exp(Het[tid] * time_step / 2);
        a[tid].x = e1 * fd[tid].x + (e1 - 1.0) * nld[tid].x / Het[tid];
        a[tid].y = e1 * fd[tid].y + (e1 - 1.0) * nld[tid].y / Het[tid];
    }
}

//m_d[tid] = fd[tid] + (c - 1.0 - time_step * Het[tid]) / (time_step * Het[tid] * Het[tid]) * (nnld[tid] - nld[tid]);
void Calc::etd41(DARR3D d, cufftDoubleComplex* fd, cufftDoubleComplex* nld, DARR3D Het, double time_step)
{
    etd41_kernel << <((nx / 2 + 1) * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (ctmp1, fd, nld, Het, time_step, (nx / 2 + 1) * ny * nz);
    cudaDeviceSynchronize();
    c2r_3d_b(ctmp1, d);
}

__global__ void etd42_kernel(cufftDoubleComplex* b, cufftDoubleComplex* fd, cufftDoubleComplex* nlda, DARR3D Het, double time_step, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        double e1 = exp(Het[tid] * time_step / 2);
        b[tid].x = e1 * fd[tid].x + (e1 - 1.0) * nlda[tid].x / Het[tid];
        b[tid].y = e1 * fd[tid].y + (e1 - 1.0) * nlda[tid].y / Het[tid];
    }
}

//m_d[tid] = fd[tid] + (c - 1.0 - time_step * Het[tid]) / (time_step * Het[tid] * Het[tid]) * (nnld[tid] - nld[tid]);
void Calc::etd42(cufftDoubleComplex* b, cufftDoubleComplex* fd, cufftDoubleComplex* nld, DARR3D Het, double time_step)
{
    etd42_kernel << <((nx / 2 + 1) * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (b, fd, nld, Het, time_step, (nx / 2 + 1) * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void etd43_kernel(cufftDoubleComplex* c, cufftDoubleComplex* fd, cufftDoubleComplex* nld, cufftDoubleComplex* nldb, DARR3D Het, double time_step, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        double e1 = exp(Het[tid] * time_step / 2);
        c[tid].x = c[tid].x * e1 + (e1 - 1.0) * (2 * nldb[tid].x - nld[tid].x) / Het[tid];
        c[tid].y = c[tid].y * e1 + (e1 - 1.0) * (2 * nldb[tid].y - nld[tid].y) / Het[tid];
    }
}

//m_d[tid] = fd[tid] + (c - 1.0 - time_step * Het[tid]) / (time_step * Het[tid] * Het[tid]) * (nnld[tid] - nld[tid]);
void Calc::etd43(cufftDoubleComplex* c, cufftDoubleComplex* fd, cufftDoubleComplex* nld, cufftDoubleComplex* nldb, DARR3D Het, double time_step)
{
    etd43_kernel << <((nx / 2 + 1) * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (c, fd, nld, nldb, Het, time_step, (nx / 2 + 1) * ny * nz);
    cudaDeviceSynchronize();
}

__global__ void etd44_kernel(cufftDoubleComplex* d, cufftDoubleComplex* fd, cufftDoubleComplex* nlda, cufftDoubleComplex* nldb, cufftDoubleComplex* nldc, DARR3D Het, double time_step, int N)
{
    // Calculate the thread index
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    // Copy data from global to shared memory
    if (tid < N) {
        double e1 = exp(Het[tid] * time_step / 2);
        double e2 = e1 * e1;
        d[tid].x = e2 * fd[tid].x + 1.0 / (time_step * time_step * Het[tid] * Het[tid] * Het[tid]) * (
            nlda[tid].x * (-4.0 - time_step * Het[tid] + e2 * (4.0 - 3.0 * time_step * Het[tid] + time_step * time_step * Het[tid] * Het[tid])) +
            2 * (nlda[tid].x + nldb[tid].x) * (2.0 + time_step * Het[tid] + e2 * (-2 + time_step * Het[tid])) +
            nldc[tid].x * (-4 - 3 * time_step * Het[tid] - time_step * time_step * Het[tid] * Het[tid] + e2 * (4 - time_step * Het[tid]))
            );
        d[tid].y = e2 * fd[tid].y + 1.0 / (time_step * time_step * Het[tid] * Het[tid] * Het[tid]) * (
            nlda[tid].y * (-4.0 - time_step * Het[tid] + e2 * (4.0 - 3.0 * time_step * Het[tid] + time_step * time_step * Het[tid] * Het[tid])) +
            2 * (nlda[tid].y + nldb[tid].y) * (2.0 + time_step * Het[tid] + e2 * (-2 + time_step * Het[tid])) +
            nldc[tid].y * (-4 - 3 * time_step * Het[tid] - time_step * time_step * Het[tid] * Het[tid] + e2 * (4 - time_step * Het[tid]))
            );
    }
}

//m_d[tid] = fd[tid] + (c - 1.0 - time_step * Het[tid]) / (time_step * Het[tid] * Het[tid]) * (nnld[tid] - nld[tid]);
void Calc::etd44(DARR3D d, cufftDoubleComplex* fd, cufftDoubleComplex* nlda, cufftDoubleComplex* nldb, cufftDoubleComplex* nldc, DARR3D Het, double time_step)
{
    etd44_kernel << <((nx / 2 + 1) * ny * nz + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (ctmp1, fd, nlda, nldb, nldc, Het, time_step, (nx / 2 + 1) * ny * nz);
    cudaDeviceSynchronize();
    c2r_3d_b(ctmp1, d);
}

void Calc::init_fEIF(DARR3D Het, double eps)
{
    double kappa = 2.0;
    double D = 1.5;

    int i, j, k;
    double cx = 2.0 * PI / xlen;
    double cy = 2.0 * PI / ylen;
    double cz = 2.0 * PI / zlen;
    cx *= cx; cy *= cy; cz *= cz;

    ARR3D tmp;
    malloc_array(tmp, nx/2+1, ny, nz);
#pragma omp for 
    for (i = 0; i < nz; i++) {
        int i2 = i < nz - i ? i : nz - i;
        for (j = 0; j < ny; j++) {
            int j2 = j < ny - j ? j : ny - j;
            for (k = 0; k < nx/2+1; k++) {
                int k2 = k < nx - k ? k : nx - k;
                double c = -(i2 * i2 * cz + j2 * j2 * cy + k2 * k2 * cx);
                tmp[i][j][k] = -eps * c * c + kappa * c / eps + D * (c - kappa / (eps * eps)) / eps;// + 0.624*c;
            }
        }
    }
    cudaMemcpy(Het, &tmp[0][0][0], (nx/2+1) * ny * nz * sizeof(double), cudaMemcpyHostToDevice);
    free_array(tmp);
}

void  Calc::fddx(DARR3D u, DARR3D s) {};
void  Calc::fddy(DARR3D u, DARR3D s) {};
void  Calc::fddx_b(DARR3D u, DARR3D s) {};
void  Calc::fddy_b(DARR3D u, DARR3D s) {};
void  Calc::fddxx(DARR3D u, DARR3D s) {};
void  Calc::fddyy(DARR3D u, DARR3D s) {};
void  Calc::fdlap(DARR3D u, DARR3D s) {};
void  Calc::fdgdsq(DARR3D u, DARR3D s) {};
