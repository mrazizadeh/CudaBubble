// mem_array.cpp: implementation of the mem_array class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mem_array.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mem_array::mem_array(fluidbase * h)
{
	m_h = h;
	for (int i = 0; i < TMP_ARRAY_MAX; i++) {
		use[i] = false;
		use_complex[i] = false;
	}
}

mem_array::~mem_array()
{
	for (int i = 0; i < TMP_ARRAY_MAX; i++) {
		if (use[i]) m_h->release_tmp_array(i);
		if (use_complex[i]) m_h->release_tmp_complex_array(i);
	}
}

DARR3D mem_array::get_array()
{
	int i;
	DARR3D tmp = m_h->get_tmp_array(i);
	use[i] = true;
	return tmp;
}
cufftDoubleComplex* mem_array::get_complex_array()
{
	int i;
	cufftDoubleComplex* tmp = m_h->get_tmp_complex_array(i);
	use_complex[i] = true;
	return tmp;
}
