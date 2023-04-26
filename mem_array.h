// mem_array.h: interface for the mem_array class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MEM_ARRAY_H__A3F94B39_0BFB_443E_9E5D_81D9259D14A3__INCLUDED_)
#define AFX_MEM_ARRAY_H__A3F94B39_0BFB_443E_9E5D_81D9259D14A3__INCLUDED_

#include "fluidbase.h"

class mem_array  
{
public:
	DARR3D get_array();
	cufftDoubleComplex * get_complex_array();
	mem_array(fluidbase * h);
	virtual ~mem_array();

private:
	bool use[TMP_ARRAY_MAX];
	bool use_complex[TMP_ARRAY_MAX];
	fluidbase * m_h;
};

#endif // !defined(AFX_MEM_ARRAY_H__A3F94B39_0BFB_443E_9E5D_81D9259D14A3__INCLUDED_)
