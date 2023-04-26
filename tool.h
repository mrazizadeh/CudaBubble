#pragma once
#define ARSUF 3  //array surfix, 3 to memory the array size, in order nx, ny, nz.
#define CACHELINE 0   //not use cashline in cuda
#define THREADS_PER_BLOCK 512

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

typedef double*** ARR3D;

bool malloc_array(ARR3D& x, int nx, int ny, int nz);
void free_array(ARR3D& x);
bool cuda_malloc_array(ARR3D& x, int nx, int ny, int nz);
bool cuda_copyin_array(ARR3D hx, ARR3D dx, int nx, int ny, int nz);
bool cuda_copyout_array(ARR3D hx, ARR3D dx, int nx, int ny, int nz);

void copyfile(char* filefrom, char* fileto);
int roundx(double x);

void write_array(ARR3D d, char* fname);
void read_array(ARR3D d, char* fname);
void zero_array(ARR3D  x);
void copy_array(ARR3D d, ARR3D s);
void copy_array_s(ARR3D d, ARR3D s);


inline void swap(double& a, double& b) {
	register double tmp;
	tmp = a; a = b; b = tmp;
}

