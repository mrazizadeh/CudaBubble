// stdafx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently 
//

#if !defined(AFX_STDAFX_H__27C3A855_D524_431A_A00F_3C748BCAF4EA__INCLUDED_)
#define AFX_STDAFX_H__27C3A855_D524_431A_A00F_3C748BCAF4EA__INCLUDED_

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#endif

#if defined (_OPENMP) 
#include <omp.h>
#endif

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#define _CPP_NUMERIC_LIMITS    //used for get off an asm error when compile in IA64.
#include <time.h>
#include <sys/timeb.h>
#include <string.h>
#include "tool.h"
#if defined (_OPENMP) 
#include <omp.h>
#endif


#ifndef WIN32
#	define _finite finite
#	define _snprintf snprintf
#	define _ftime ftime
#	define _timeb timeb
#else
#	include <float.h>
#endif

// TODO: reference additional headers your program requires here

#define TIME_STEP_OUTPUT 1e-6
//#define _OUTPUT_TIME_

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__27C3A855_D524_431A_A00F_3C748BCAF4EA__INCLUDED_)
