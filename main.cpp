#include "belastic.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
	//test for global tools functions
	//	test(); return 1;

	int i = 0;
	if (argc >= 2) {
		i = atoi(argv[1]);
	}

	printf("Elastic bending energy program!\n");

	belastic t(i);
	t.run_static_fEIF();
	//t.run_static();
	//	t.run_d();
	//	t.a_run();
	//	t.tanh_adjust();
	//t.test1();

	return 1;
}

