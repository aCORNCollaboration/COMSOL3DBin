//
//  ReadField.c
//  COMSOLReader2
//
//  This reads a text description of a nested set of fields and creates
//  a single object that describes the complete system.
//
//  The format handled is illustrated below
//  [fields [<directory>]]
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
//  Modified BCollett 5/8/14 Added support for the fields prefix that sets
//  the working directory for the field files.
//

#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <math.h>
#include "readfield.h"
#ifdef undef //__WX__
#include "FieldViewerApp.h"
#else
#define eprintf printf
#define iprintf printf
#define oprintf printf
#define wprintf printf
#endif
//
//  File globals.
//
static const char* delims = "\r\n\t ,";
static bool FieldInField(const CD3Data* od, const CD3Data* nd);
static bool AddField(CD3Data* od, const CD3Data* nd, const char* linBuff);
static void CheckDir();
static bool SoftPtInBounds(const CD3Data* dp, const double coord[3]);
//
//  ParseFieldSet constucts a complete tree of nested fields from a text
//  description in an open file.
//  On success dp points to the root of the tree. On failure averything is
//  a mess and you should throw up your hands and quit!
//
bool ParseFieldSet(CD3Data* dp, FILE* ifp)
{
  char linBuff[1028];
  char *verb;
  for (;;) {
    //
    //  Read in a line.
    //
    if (fgets(linBuff, 1024, ifp) == NULL) {
/*      eprintf("ParseFieldSet: Unexpected end of input while parsing fields.\n");
      return false;*/
      break;
    }
    //
    //  Extract first word.
    //
    verb = strtok(linBuff, delims);
    //
    // See what we have, if anything.
    //
    if (verb == NULL) return true;
    if (strcmp(verb, "fields") == 0) {
      char* path = strtok(NULL, delims);
      if (strlen(path) > 0) {
        int theErr = chdir(path);
        if (theErr != 0) {
          switch (errno) {
            case ENOTDIR:
              eprintf("Path '%s' is not a directory\n", path);
              break;
              
            case EACCES:
              eprintf("You do not have search access to path '%s'.\n", path);
              break;
              
            case ENAMETOOLONG:
              eprintf("Path '%s' is too long.\n", path);
              break;
              
            case ENOENT:
              eprintf("Component in path '%s' does not exist.\n", path);
              break;
              
            case ELOOP:
              eprintf("Path '%s' has too many symbolic links.\n", path);
              break;
              
            default:
              break;
          }
        }
      }
    } else     if (strcmp(verb, "cfield") == 0) {
      if (!ParseCField(dp, ifp)) {
        eprintf("ParseFieldSet: Failed to find cfield starting at %s\n",
                linBuff);
        return false;
      }
    } else if (strcmp(verb, "field") == 0) {
      const char* iname = strtok(NULL, delims);
      if (!ParseField(dp, iname)) {
        eprintf("ParseFieldSet: Failed to find field starting at %s\n",
                linBuff);
        return false;
      }
    } else {
      eprintf("ParseFieldSet: Expecting 'field' or 'cfield', found %s\n",
              verb);
      return false;
    }
  }
  return true;
}
//
//  This is simple. It reads in a terminal field.
//
bool ParseField(CD3Data* dp, const char* name)
{
  const char* defaultFieldName = "Field name error.";
  char* fieldName = NULL;
  bool success = false;
  size_t nNameChar = -1;
  FILE* nIfp;
  if ((name == NULL) || (strlen(name) == 0)) {
    return false;
  }
  nIfp = fopen(name, "rb");
  if (NULL == nIfp) {
    eprintf("ParseCField: Could not open field file %s, error %d.\n",
            name, errno);
    perror("ParseCField: Main field file--");
    CheckDir();
    return false;
  }
  success = CD3ReadBinary(dp, nIfp);
  fclose(nIfp);
  if (!success) {
    eprintf("ParseCField: ReadBinary could not load field %s.\n",
            name);
    return false;
  }
  oprintf("Loaded field %s.\n",name);
  nNameChar = strlen(name);
  if ((nNameChar < 1) || (nNameChar > 256)) {
    eprintf("ParseCField: field name %s too long or too short.\n",
            name);
    dp->mFieldName = defaultFieldName;
  } else {
    fieldName = (char*) malloc(nNameChar+2);
    if (NULL == fieldName) {
      eprintf("ParseCField: unable to allocate space for name %s.\n",
              name);
      dp->mFieldName = defaultFieldName;
    } else {
      strncpy(fieldName, name, nNameChar);
      fieldName[nNameChar] = 0;
      dp->mFieldName = fieldName;
    }
  }
  return true;
}

//
//  This is more complex.
//  It may or may not read in a field itself but then it looks
//  for child fields and gets them created. The FILE leads to the
//  text description of the field hierarchy.
//
bool ParseCField(CD3Data* dp, FILE* ifp)
{
  char cname[64];
  char linBuff[1028];
  char* verb;
  const char* ename;
  CD3Data* newData;
  //
  //  If we got a name then we start by reading in that field.
  //  In any case we save name info for end matching.
  //
  const char* name = strtok(NULL, delims);
  if (name == NULL) {
    cname[0] = 0;
  } else {
    strncpy(cname, name, 63);
    if (!ParseField(dp, name)) {
      return false;
    }
  }
  //
  //  Now expect a set of daughter fields.
  //
  for (;;) {
    if (fgets(linBuff, 1024, ifp) == NULL) {
      eprintf("ParseCField: Unexpected end of input while parsing cfield.\n");
      return false;
    }
    //
    //  Extract first word.
    //
    verb = strtok(linBuff, delims);
    //
    // See what we have
    //
    if (strcmp(verb, "cfield") == 0) {
      //
      //  Cons up a new field and recurs.
      //
      newData = (CD3Data*) malloc(sizeof(CD3Data));
      if (NULL == newData) {
        eprintf("ParseCField: Failed to get space for new CField.\n");
        return false;
      }
      if (ParseCField(newData, ifp)) {
        if (!AddField(dp, newData, linBuff)) {
          free(newData);
        }
      }
    } else if (strcmp(verb, "field") == 0) {
      //
      //  Get the file name.
      //
      const char* iname = strtok(NULL, delims);
      if (NULL == iname) {
        eprintf("ParseCField: Could not find name of field file in %s.\n",
                linBuff);
        return false;
      }
      //
      //  Cons up a new field. Get it read in, then install here.
      //
      newData = (CD3Data*) malloc(sizeof(CD3Data));
      if (NULL == newData) {
        eprintf("ParseCField: Failed to get space for new CField.\n");
        return false;
      }
      if (ParseField(newData, iname)) {
        if (!AddField(dp, newData, linBuff)) {
          free(newData);
        }
      }
    }  else if (strcmp(verb, "end") == 0) {
      ename = strtok(NULL, delims);
      if (ename == NULL) {
        ename = "";
      }
      if (strcmp(ename, cname) != 0) {
        eprintf("ParseCField: End name %s does not match start name %s\n",
                ename, cname);
        return false;
      }
      break;
    } else {
      eprintf("ParseCField: Expecting 'field' or 'cfield', found %s\n",
              verb);
      return false;
    }
  }
  return true;
}
//
//  This tests that the bounding box of nd lies entirely within the
//  bounds of od. Have to test only extremal corners.
//  This is done with a soft comparison, NOT the exact
//  comparison of PtInBounds.
//
bool FieldInField(const CD3Data* od, const CD3Data* nd)
{
  double coord[3];
  if (od->mField == NULL) {
    return true;    // If no containing field then succeed.
  }
  coord[0] = nd->mMin[0];
  coord[1] = nd->mMin[1];
  coord[2] = nd->mMin[2];
  if (!SoftPtInBounds(od, coord)) {
    SoftPtInBounds(od, coord);
    return false;
  }
  coord[0] = nd->mMax[0];
  coord[1] = nd->mMax[1];
  coord[2] = nd->mMax[2];
  if (!SoftPtInBounds(od, coord)) {
    SoftPtInBounds(od, coord);
    return false;
  }
  return true;
}
//
//  AddField installs a new field (nd) into the next subfield spot in the
//  the old field (od).
//
static bool AddField(CD3Data* od, const CD3Data* nd, const char* linBuff)
{
  //
  //  Check that the new field is totally contained inside this one
  //  and if so install in a free slot.
  //
  if (FieldInField(od, nd)) {
    if (od->mNSubField < kNSub) {
      od->mSubField[od->mNSubField++] = nd;
    } else {
      eprintf("ParseCField: No room for new field at %s.\n",
              linBuff);
      return false;
    }
  } else {
    eprintf("ParseCField: New field %s not contained in old field %s.\n",
            nd->mFieldName, od->mFieldName);
    return false;
  }
  return true;
}


void CheckDir()
{
  DIR *dp;
  struct dirent *ep;
  dp = opendir ("./");
  
  if (dp != NULL)
  {
    iprintf("Contents of current directory.\n");
    while ((ep = readdir (dp)) != NULL)
      iprintf("--%s\n", ep->d_name);
    
    (void) closedir (dp);
  }
  else
    perror ("Couldn't open the directory");
}
//
//  Slightly soft pt in bounds check.
//
bool SoftPtInBounds(const CD3Data* dp, const double coord[3])
{
  int i;
  for (i = 0; i < 3; i++) {
    double eps = (fabs(coord[i]) < 1e-6) ? 1e-6 : 1e-6 * fabs(coord[i]);
    if (((coord[i] - dp->mMin[i]) < -eps) ||
        ((coord[i] - dp->mMax[i]) > eps)) {
      return false;
    }
  }
  return true;
}

