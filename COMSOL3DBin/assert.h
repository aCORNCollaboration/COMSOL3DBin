/*
* assert.h
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
#ifndef __ASSERT_H
#define __ASSERT_H

extern int gBC_ASSERT_ABORT;

void BC_DO_ASSERT(const char* testStr, const char* fileName, int lineNumber);

#define assert(x) if (!(x)) BC_DO_ASSERT("x", __FILE__, __LINE__)
#endif


