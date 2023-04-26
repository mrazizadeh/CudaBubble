# cudaBubble

cudaBubble is to solve the Allen-Cahn phase-field equation for the Elastic Bending Energy model for vesicles.
The code is migrated to CUDA from OpenMP, so that we could take advantage of the fast GPU calculation capacity. 
The speed up in CUDA code is 10 to 100 times faster than serial CPU code. 

The algorithms in this code includes: 
Forward Euler - See the paper "Modeling the Spontaneous Curvature Effects in Static Cell Membrane Deformations by a Phase Field Formulation", or the paper "Simulating the Deformation of Vesicle Membranes under Elastic Bending Energy in Three Dimensions"
Exponential Time Differencing Runge-Kutta Methods - See the paper "Efficient and Stable Exponential Time Differencing Runge-Kutta Methods for Phase Field Elastic Bending Energy Models" 
You can google those paper from google scholar.
All models assume periodic boundary condition.

# Packages Required

CUDA Toolkit - Please google "cuda Toolkit" and install the newest version. Cuda Toolkit includes cuFFT. This code use DFFT for calculation of the derivatives and we use cuFFT for fast FFT calculations. 
For Windows system, cuda Toolkit will be integrated with Visual Studio. 

## vcpkg
You can download vcpkg from Git. And use it to install extra packages you want to use for further development. For example, you can use vcpkg to install FFTW3. Google vcpkg for details.

# System requirements

To use CUDA, I believe that you will need a Nvidia GPU card with your machine.
I do the coding on a Windows 11 system, with Visual Studio, plus CUDA Toolkit.
The whole code should be easy to run on a Linux system, but I did not do that. Maybe I will add it in the future.

# How to run it?

To run the code, you will need a "results" folder. Under the foler, there are some configuration files, named as "ini?.dat", depending on your running parameters. (The default one is ini0.dat).
For example, if you run your code in command line:
>> cudaBubble.exe 3
Then the code will look for "ini3.dat" file and read the first line. (only the first line will be read, so that you can put your comments on the other rows.)

## configuration file
The first line gives the mesh size, parameters, etc. The choice of those parameters can be found in the tempelate ini.dat files.

### Example
First two lines are:
100  128  128  128       1    1    1      1	1e-7       5e-6             1		          2    -1     1        200       0      30      1.0        1   	1     -1		-1	       1.0e4       1.0e4     1.0e4	1.0e4	100.0     30                 0.05	0	30.0	70	10		-1
folder nx  ny  nz   lenx leny lenz   shape	time_step fluid_time_step  fixed_time_step  ninit nflush surctrl  epssize ncur   angle    Re         eta      gamma   alpha    		beta        	M1          M2        M3	M4	tmax    output_interval    gweight	einit	k	c	lineten		alpha_3

### paramters explained
The first one is the folder number, in this case, you will need a subfolder named "100" under the results folder and all the output will be put into this folder.
The next three parameters are for mesh sizes, nx, ny, nz. We suggest you use the same 2^n number for all those mesh sizes. Although nz may be different with nx/ny. But we did not test it thoughoutly.
The next three parameters are the lenth for the intervals on three axis. Usually we fix lenx/leny/lenz to be 1, which is a multiply of 2Pi, i.e., the interval is from \[-PI, PI\]. 
The next parameter is shape: if it is 1, then we will run the shape adjust so that the phase field elastic bending energy lower enough, and thus the phase field function will be more like a tanh function. During the shape adjust, there is no volume/surface area constrain.
if it is 0, no shape adjustment. 
The next parameter is the initial time_step.
The next parameter is the initial time step for solving Navier-Stokes Equation, the original code has been coupled with Navier-Stokes equation to study the dynamics of vesicles in fluid. In this code, that part is removed.
The next parameter is the fixed_time_step, 0 means adjust time step (Trying to double the time step size every 10 iterations. Any step, if energy does not decrease, then cut in half the time step.)
The next parameter is "ninit", which can be choosed from a set of numbers, each one is an initial shape configuration. One can add more into the code. See 'fluidbase.cpp' file.
The next parameter is "nflush", which is for initial fluid velocity type. Keep to be -1 here.
The next parameter is "surctrl", with volume/surface area constrain or not, 1 or 0. 
The next parameter is "epssize", in the code, the eps = "epssize/100" \* h, where 'h' is the grid size. 
The next parameter is "ncur", used for different spontenuous curvature settings. Keep to be 0 here.
The next parameter is "angle", used for the picture output. Every some time interval, we output the data into a set of files, including a picture of the vesicle. In the picture, we draw two period in z-direction. Also two period from left to right. The left side part is cut from the middle with an angle specified here. We cut the vesicle to see the insider.
The next parameter is "Re", Renolds numbe for fluid.
The next parameter is "eta", the elastic bending module coefficient, most keep to be 1.
The next parameter is "gamma", keep to be 1. 
The next two parameters are "alpha", "beta", which is the volume and surface area constrains. If set them to be -1, then automatically use the volume and surface area at time 0 to be the constrain values.
The next four parameters are the penalty coefficients "M1", "M2","M3", "M4", where M1 is the penalty coefficient for volume, M2 is for surface area. M3 and M4 are for other constrains, not used in this code.
The other parameters: 
"tmax" - Maximum time, after reached this time, terminate the code and exit.   
"output_interval" - Every how long time to output a middle results.
"gweight" - pameter to use for a vesicle has weight, not used here.
"einit" - For multi-component, label phase, not used here.
"k,	c, lineten, alpha_3" - for multi-component, with line tention energy, etc. Not used here.

## Resume running
If the code is terminated for whatever reason, if we restart it, it will resume running. 
If you would like to start a new running, you can simply delete all the data files in the output subfolder. 
The code check whether some specific data files exisits or not, if they are existing, then it will read those data files and resume the running.

### Change parameters
You can change some parameters when resume running. Just change the parameters you want in the configuration file, and the code will resume the running using the new parameters.

