//
//  GSSmooth.h
//  COMSOL3DBin
//
//  Routine to Gauss-Seidel smooth a 3D electric field array.
//  Builds a type array from a geometry file.
//
//  Created by Brian Collett on 7/29/15.
//  Copyright (c) 2015 Brian Collett. All rights reserved.
//

#ifndef __COMSOL3DBin__GSSmooth__
#define __COMSOL3DBin__GSSmooth__

#include <stdio.h>
#include "COMSOLData3D.h"

int GSSmooth(const char* fname, CD3Data* dp, int nPass);

#endif /* defined(__COMSOL3DBin__GSSmooth__) */
