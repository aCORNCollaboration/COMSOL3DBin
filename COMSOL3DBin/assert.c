/*
* assert.c
*
* Simple assertion system.
* Have a function that does the work, a macro to
* make it easy to use, and a global that controls the
* mode of operation.
* The macro simply adds the file and line number to
* the argument list of the function.
* The function prints an error message and then either
* returns or exits depending on the state of the
* gBC_ASSERT_ABORT.
*/
#include <stdio.h>
#include <stdlib.h>
#include "assert.h"

int gBC_ASSERT_ABORT = 1;

void BC_DO_ASSERT(const char* testStr, const char* fileName, int lineNumber)
{
  fprintf(stderr,"Assertion %s failed at line %d in file %s.\n", testStr, lineNumber, fileName);
  if (gBC_ASSERT_ABORT) {
    exit(-1);
  }
}

