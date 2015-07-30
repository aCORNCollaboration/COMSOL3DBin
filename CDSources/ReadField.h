//
//  ReadField.h
//  COMSOLReader2
//
//  This reads a text description of a nested set of fields and creates
//  a single object that describes the complete system.
//
//  The format handled is illustrated below
//  cfield [<name1>]
//    field <name2>
//    cfield <name3>
//      field <name4>
//    end <name3>
//  end [<name1>]
//
//  Note that a field MUST have a file but a cfield need
//  not. If it holds only a disjoint set of fields then the
//  name slot can be empty.
//  In all cases the end tag must match the field tag.
//
//
//  Created by Brian Collett on 3/12/14.
//  Copyright (c) 2014 Brian Collett. All rights reserved.
//

#ifndef _ReadField_h
#define _ReadField_h

#include <stdbool.h>
#include "COMSOLData3D.h"

__BEGIN_DECLS
bool ParseFieldSet(CD3Data* dp, FILE* ifp);
bool ParseCField(CD3Data* dp, FILE* ifp);
bool ParseField(CD3Data* dp, const char* name);
__END_DECLS


#endif
