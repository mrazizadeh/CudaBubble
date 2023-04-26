// belastic.h: interface for the belastic class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include "fluidbase.h"
#include "Calc.h"
#include "phasefield.h"

class belastic : public fluidbase  
{
public:
	belastic(int folderindex);
	virtual ~belastic();

	void run_static();
	void run_static_fEIF();
	void test1();
	double intf(DARR3D s);

private:
	void get_sum_fgqterm(phasefield * pd, DARR3D nld);
	double getstaticenergy(phasefield * pd, bool constrained = true);
	void outputall();
	double outputenergy();
	double gettotalenergy();
	double geterror(phasefield* pd);
	void shape_adjust();

	void get_nld_term(phasefield* pd, cufftDoubleComplex* nld);
	void ETD2RK(DARR3D Het);
	void ETD4RK(DARR3D Het);

	Calc* calc;
	phasefield* m_pd;
	phasefield* m_pa;  //used as a helper for fEIF calculation

	DARR3D tmp;
};

