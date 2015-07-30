/*
 *  COMSOLData3D.cpp
 *  COMSOLReader2
 *
 *  This is a specialisation of the COMSOLReader class that only
 *  supports 3D data files containing 3D fields. It checks that
 *  the file not merely has three dimensions but that all three
 *  spatial dimensions are active, that is have more than one value
 *  of that coordinate represented. This means that GetValueAtPoint
 *  can be somewhat simpler.
 *
 *  Created by Brian Collett on 2/22/14.
 *  Copyright (c) 2014 Brian Collett. All rights reserved.
 *  BCollett 3/12/14 Add support for binary storage of 3D files.
 *  They have a header containing information on the structure and
 *  then a single array of doubles as the data organized
 *  data[((((index[2]*dp->mNVal[1])+index[1])*dp->mNVal[1])+index[0])+comp]
 *  where comp is 0 for x, 1 for y, and 2 for z. Restructured mField to match.
 *
 *  BCollett 3/14/14 Merge with 2D version. We keep a flag that tells us
 *  how our internal data are stored and and that allows us to forward
 *  field requests to the right code path.
 *  NOTE that the two types of field now supported are fully 3D fields in
 *  boxes aligned with the Cartesian axes and 2D fields axisymmetric about
 *  the z axis only.
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include "COMSOLData3D.h"

//
//  Forward declarations for file scope helper functions.
//
static CDError Init3D(CD3Data* dp, const CDData* cdp);
static CDError Init2D(CD3Data* dp, const CDData* cdp);
static bool GetAxEAtPoint(const CD3Data* dp, const double coord[3], double* EField);
static bool Get3DEAtPoint(const CD3Data* dp, const double coord[3], double* EField);
static bool Get2DEAtPoint(const CD3Data* dp, const double coord[2], double* EField);

//
//  Define this if you want to bounds check every value.
//
#define CD3BoundsCheck 1

//
//  Global for magic number.
//
uint32_t gCD3Magic = 'CD3B';
//
//  Global for header size.
//
uint32_t gCD3HeadLength = 512;
//
//  Init fills in the data structure using the information in the file.
//  BEWARE: CDInit allocates a lot of storage. We MUST ensure that that
//  storage gets disposed of before we leave. This means that once CDInit
//  has run ALL error must go through an exit code path, hence the rare
//  use of a goto and label.
//
CDError CD3Init(CD3Data* dp, const char* fname)
{
  CDError theErr;
  int dim, nActive = 0;
  //
  //  Start by constructing a CDData from the file.
  //
  CDData cData;
  theErr = CDInit(&cData, fname);
  if (theErr != kCDNoErr) {
    goto ErrorExit;
  }
  //
  //  And filling in default elems in the CD3Data.
  //
  dp->mType = kCD3Error;
  dp->mStride = 0;
  dp->mNSubField = 0;
  dp->mSubField[0] = dp->mSubField[1] = NULL;
  dp->mSubField[2] = dp->mSubField[3] = NULL;
  dp->mField = NULL;
  dp->mFieldName = fname;
  //
  //  Have a real file nicely parsed out. Figure out which kind it was.
  //  Must be either two active dimensions and two expressions or
  //  three active dimensions and three expressions.
  //  I did the tests on separate lines to make debugging a little easier.
  //  Start with needing three actual dimensions.
  //
  if (cData.mNDimension != 3) {
    fprintf(stderr,
            "Expected three dimensions, found %d.\n",
            cData.mNDimension);
    theErr = kCDBadStructure;
    goto ErrorExit;
  }
  //
  //  Then find out how many are active.
  //
  for (dim = 0; dim < 3; dim++) {
    if (cData.mRange[dim].mActive)
      nActive++;
  }
  if (nActive == 2) {
    theErr = Init2D(dp, &cData);
  } else if (nActive == 3) {
    theErr = Init3D(dp, &cData);
  } else {
    fprintf(stderr,
            "Expected two or three active dimensions, found %d.\n",
            nActive);
    theErr = kCDBadStructure;
    goto ErrorExit;
  }
ErrorExit:
  //
  //  This makes sure that we dispose of the remaining storage in CData.
  //  ALL code paths exit through here (hence the label!
  //
  CDFinish(&cData);
  return theErr;
}
//
//  InitFEMM does the same for a file extracted from FEMM
//  and written by Matlab or Octave.
//  This is somewhat simpler because the file must hold 2D data
//  in the form
//  x y Ex Ey
//  stored with the y dimension varying fastest.
//  UNLIKE the COMSOL data we have no information on how many
//  entries there will be in the file. There is no header to
//  help us here. However, the lines are all roughly the same
//  length because they are written with a fixed format. Thus
//  we can get a good estimate of the amount of space needed
//  for the arrays by counting the number of characters in the
//  first line and dividing the file length by that.
//
CDError CD3InitFEMM(CD3Data* dp, const char* fname)
{
  int line = 0;
  char linebuff[1024];
  char* string;
  unsigned int nCharPerLine;
  unsigned int nLine;
  unsigned long fileSize;
  unsigned int nXCopy = 0;
  struct stat st;
  double *xVals, *yVals, *exVals, *eyVals;
  bool xSame = true;
  int row, col;
  //
  //  Make sure that we can read the file.
  //
  FILE* ifp = fopen(fname, "rt");
  if (NULL == ifp) {
    fprintf(stderr, "Failed to open file %s.\n",fname);
    return kCDCantOpenIn;
  }
  string = fgets(linebuff, 1024, ifp);
  if (NULL == ifp) {
    fprintf(stderr, "Could not read from file %s.\n", fname);
    return kCDCantOpenIn;
  }
  //
  //  Use that and the filesize to get the number of lines and thus
  //  the number of entries in each array.
  //
  nCharPerLine = (unsigned int) strlen(string) - 2;
  if (fstat(fileno(ifp), &st) != 0) {
    fprintf(stderr, "Cannot stat file %s.\n", fname);
    return kCDCantOpenIn;
  }
  fileSize = st.st_size;
  nLine = (unsigned int) (fileSize / nCharPerLine);
  printf("File %s has about %d lines.\n", fname, nLine);
  //
  //  Now we can get space for the data.
  //
  xVals = (double *) malloc(nLine * sizeof(double));
  if (xVals == NULL) {
    fprintf(stderr, "Failed to alocate %d slots for x values.\n", nLine);
    return kCDAllocFailed;
  }
  yVals = (double *) malloc(nLine * sizeof(double));
  if (yVals == NULL) {
    fprintf(stderr, "Failed to alocate %d slots for y values.\n", nLine);
    return kCDAllocFailed;
  }
  exVals = (double *) malloc(nLine * sizeof(double));
  if (exVals == NULL) {
    fprintf(stderr, "Failed to alocate %d slots for Ex values.\n", nLine);
    return kCDAllocFailed;
  }
  eyVals = (double *) malloc(nLine * sizeof(double));
  if (eyVals == NULL) {
    fprintf(stderr, "Failed to alocate %d slots for Ey values.\n", nLine);
    return kCDAllocFailed;
  }
  //
  //  Fill in default elems in the CD3Data.
  //
  dp->mType = kCD3Data2;
  dp->mNSubField = 0;
  dp->mSubField[0] = dp->mSubField[1] = NULL;
  dp->mSubField[2] = dp->mSubField[3] = NULL;
  dp->mFieldName = fname;
  //
  //  The ones associated with the structure of the array.
  //  NOTE the input array has entries Ex and Ey but we will
  //  read them in a slice of a 3D array mapping input x to
  //  dimension 1, input y to dimension 2, and setting dimension
  //  0 to 0.0.
  //
  dp->mMin[0] = 0.0;
  dp->mMax[0] = 0.0;
  dp->mMin[1] = dp->mMin[2] = DBL_MAX;
  dp->mMax[1] = dp->mMax[2] = -DBL_MAX;
  //
  //  And read text file in.
  //
  do {
    int nRead = sscanf(linebuff, "%lf %lf %lf %lf",
                       &xVals[line], &yVals[line],
                       &exVals[line], &eyVals[line]);
    if  (nRead < 4) {
      fprintf(ifp, "Only read %d of 4 values on line %d of file %s.\n",
              nRead, line, fname);
      return kCDCantOpenIn;
    }
    //
    //  As we go we count the number of identical x values to find the
    //  first dimension of the data.
    //
    if (xSame) {
      if (xVals[line] == xVals[0]) {
        nXCopy++;
      } else {
        xSame = false;
      }
    }
    //
    //  Track ranges in x an y.
    //
    if (xVals[line] < dp->mMin[1]) {
      dp->mMin[1] = xVals[line];
    }
    if (yVals[line] < dp->mMin[2]) {
      dp->mMin[2] = yVals[line];
    }
    if (xVals[line] > dp->mMax[1]) {
      dp->mMax[1] = xVals[line];
    }
    if (yVals[line] > dp->mMax[2]) {
      dp->mMax[2] = yVals[line];
    }
    line++;
  } while (fgets(linebuff, 1024, ifp) != NULL);
  fclose(ifp);
  //
  //  Check that x min is 0.
  //
  if (dp->mMin[0] != 0) {
    fprintf(stderr,
            "Loading axisymmetric data, x must have min=0.0.");
    return kCDBadStructure;
  }
  //
  //  Make sure that we have a strictly rectangular array.
  //
  if (line % nXCopy != 0) {
    fprintf(stderr, "Error checking rectangular structure. Remainder = %d.\n",
            line % nXCopy);
    return kCDBadStructure;
  }
  //
  //  Now can do nVals and the Deltas.
  //  Note that they have to be set up to reflect the FINAL disposition
  //  of the data, NOT the way they are stored in the text file.
  //
  dp->mNVal[0] = 1;
  dp->mNVal[1] = line / nXCopy;
  dp->mNVal[2] = nXCopy;
  dp->mStride = dp->mNVal[1];
  dp->mDelta[1] = (dp->mMax[1] - dp->mMin[1])/(dp->mNVal[1] - 1);
  dp->mDelta[2] = (dp->mMax[2] - dp->mMin[2])/(dp->mNVal[2] - 1);
  dp->mDelta[0] = dp->mDelta[1];
  //
  //  Redo the limits to be those of the box surrounding
  //  the region that we will interpolate.
  //
  dp->mMin[0] = -dp->mMax[1];
  dp->mMax[0] = dp->mMax[1];
  dp->mMin[1] = -dp->mMax[1];
  //
  //  Now allocate space for the final data.
  //
  dp->mField = (double *) malloc(line * 2 * sizeof(double));
  if (dp->mField == NULL) {
    return kCDAllocFailed;
  }
  //
  //  Copy the data into place.
  //
  for (row = 0; row < dp->mNVal[2]; row++) {
    for (col = 0; col < dp->mNVal[1]; col++) {
      dp->mField[2*(row * dp->mStride + col)+0] = exVals[col * nXCopy + row];
      dp->mField[2*(row * dp->mStride + col)+1] = eyVals[col * nXCopy + row];
    }
  }
  
  return kCDNoErr;
}
//
//  Finish tidies up after us, releasing our storage. Every call of
//  COMSOLDataInit should be balanced by a call to COMSOLDataFinish.
//
void CD3Finish(CD3Data* dp)
{
  if (NULL != dp->mField) {
    free(dp->mField);
  }
}
//
//  Accessor.
//
//  First checks whether a point is inside this field.
//
bool PtInBounds(const CD3Data* dp, const double coord[3])
{
  int i;
  for (i = 0; i < 3; i++) {
    if ((coord[i] < dp->mMin[i]) || (coord[i] > dp->mMax[i])) {
      return false;
    }
  }
  return true;
}
//
//  Then one that hands the work off to helpers who understand the different
//  internal structures.
//  Note that it has to check whether to pass the point off to its daughters
//  first.
//
bool CD3GetEAtPoint(const CD3Data* dp, const double coord[3], double* EField)
{
  int i;
  if (dp->mType > 1) {
    fprintf(stderr, "Invalid field type %d.\n", dp->mType);
    return false;
  }
  //
  //  Search daughters.
  //
  for (i = 0; i < dp->mNSubField; i++) {
    if (PtInBounds(dp->mSubField[i], coord)) {
      return CD3GetEAtPoint(dp->mSubField[i], coord, EField);
    }
  }
  //
  //  Can we handle it ourselves?
  //
  if (PtInBounds(dp, coord)) {
    return (dp->mType == kCD3Data2) ?
    GetAxEAtPoint(dp, coord, EField) :
    Get3DEAtPoint(dp, coord, EField);
  }
  //
  //  Otherwise this was a bust.
  //
  return false;
}
//
//  This is very similar except that intead of returning a field
//  value it returns the name of the file from whose data the
//  the field came.
//
const char* CD3GetNameAtPoint(const CD3Data* dp, const double coord[3])
{
  int i;
  if (dp->mType > 1) {
    fprintf(stderr, "Invalid field type %d.\n", dp->mType);
    return "Invalid field type";
  }
  //
  //  Search daughters.
  //
  for (i = 0; i < dp->mNSubField; i++) {
    if (PtInBounds(dp->mSubField[i], coord)) {
      return dp->mSubField[i]->mFieldName;
    }
  }
  //
  //  Can we handle it ourselves?
  //
  if (PtInBounds(dp, coord)) {
    return dp->mFieldName;
  }
  //
  //  Otherwise this was a bust.
  //
  return "No field found";
}

//
//  File operations.
//  The first writes a complete field to binary with the magic header.
//  It needs a complete CD3Data and a FILE open for writing.
//  The only difference between a 3D and 3D file is the amount of data
//  to write.
//
bool CD3WriteBinary(CD3Data* dp, FILE* ofp)
{
  int i;
  int success = false;
  CD3Header* head = (CD3Header *) malloc(gCD3HeadLength);
  if (head == NULL) {
    fprintf(stderr, "CD3WriteBinary not allocate header.\n");
    return false;
  }
  //
  //  Start by filling in the header fields.
  //
  head->magic = gCD3Magic;
  head->dataOffset = gCD3HeadLength;
  if (gFieldFileName != NULL) {
    strncpy(head->fileName, gFieldFileName, 63);
  }
  if (gModelFileName != NULL) {
    strncpy(head->modelName, gModelFileName, 63);
  }
  head->dp.mType = dp->mType;
  head->dp.mStride = dp->mStride;
  printf("Field type %d, stride %d\n", head->dp.mType, head->dp.mStride);
  head->dp.mNSubField = 0;
  head->dp.mField = 0;
  head->dp.mSubField[0] = head->dp.mSubField[1] = NULL;
  head->dp.mSubField[2] = head->dp.mSubField[3] = NULL;
  for (i = 0; i < 3; i++) {
    head->dp.mNVal[i] = dp->mNVal[i];
    head->dp.mMin[i] = dp->mMin[i];
    head->dp.mMax[i] = dp->mMax[i];
    head->dp.mDelta[i] = dp->mDelta[i];
    printf("Dim %d: %d vals %f to %f by %f\n", i,
           head->dp.mNVal[i],
           head->dp.mMin[i],
           head->dp.mMax[i],
           head->dp.mDelta[i]);
  }
  //
  //  Write header and data to disk.
  //
  if (fwrite(head, 1, gCD3HeadLength, ofp) == gCD3HeadLength) {
    int npoint = dp->mNVal[0] *  dp->mNVal[1] *  dp->mNVal[2];
    switch (dp->mType) {
      case kCD3Data2:
        npoint *= 2;
        break;
        
      case kCD3Data3:
        npoint *= 3;
        break;
        
      default:
        fprintf(stderr, "CD3WriteBinary:Invalid file type %d.\n", dp->mType);
        return false;
    }
    printf("%d = %d * %d * %d\n", npoint, dp->mNVal[0],dp->mNVal[1],dp->mNVal[2]);
    if (fwrite(dp->mField, sizeof(double), npoint, ofp) != npoint) {
      fprintf(stderr, "CD3WriteBinary:Failed to write data.\n");
    } else {
      printf("CD3WriteBinary wrote %d data values.\n", npoint);
      success = true;
    }
  } else {
    fprintf(stderr, "CD3WriteBinary:Failed to write header.\n");
  }
  //
  //  Dispose of header and are done.
  //
  free(head);
  return success;
}
//
//  The second constructs a CD3Data field from a binary file. It is the
//  binary equivalent of CD3Init for text files.
//  Because we allocate storage that must be thrown away even if an
//  error occurs we have to exit with a label!!!
//
bool CD3ReadBinary(CD3Data* dp, FILE* ifp)
{
  int i, npoint, nActive = 0;
  int success = false;
  //
  //  Get space for header, read it in, and make sure it is valid.
  //
  CD3Header* head = (CD3Header *) malloc(gCD3HeadLength);
  if (head == NULL) {
    fprintf(stderr, "CD3ReadBinary: Could not allocate header.\n");
    return false;
  }
  if (fseek(ifp, 0L, SEEK_SET) != 0) {
    fprintf(stderr, "CD3ReadBinary: Could not seek file to zero.\n");
    goto Finish;
  }
  if (fread(head, 1, gCD3HeadLength, ifp) != gCD3HeadLength) {
    fprintf(stderr, "CD3ReadBinary: Could not read header.\n");
    goto Finish;
  }
  if (head->magic != gCD3Magic) {
    fprintf(stderr,
            "CD3ReadBinary: Header magic number %x does not match %x.\n",
            head->magic, gCD3Magic);
    goto Finish;
  }
  //
  //  Copy CD3Data into place doing some checking as we go.
  //
  if ((head->dp.mType == kCD3Data2) || (head->dp.mType == kCD3Data3)) {
    dp->mType = head->dp.mType;
  } else {
    fprintf(stderr,
            "CD3ReadBinary: Invalid field type %d.\n",
            head->dp.mType);
    goto Finish;
  }
  //
  //  Copy header data into dp counting number of inactives.
  //
  for (i = 0; i < 3; i++) {
    dp->mNVal[i] = head->dp.mNVal[i];
    dp->mMin[i] = head->dp.mMin[i];
    dp->mMax[i] = head->dp.mMax[i];
    dp->mDelta[i] = head->dp.mDelta[i];
    if (dp->mNVal[i] > 1) {
      ++nActive;
      if (dp->mMax[i] <= dp->mMin[i]) {
        fprintf(stderr,
                "CD3ReadBinary: Data error. Dimension %d  max <= min.\n",
                i);
        goto Finish;
      }
      if (dp->mDelta[i] <= 0.0) {
        fprintf(stderr,
                "CD3ReadBinary: Data inalid. Dimension %d  delta = %f.\n",
                i, dp->mDelta[i]);
        goto Finish;
      }
    }
  }
  if (((dp->mType == kCD3Data2) && (nActive != 2)) ||
      ((dp->mType == kCD3Data3) && (nActive != 3))) {
    fprintf(stderr,
            "CD3ReadBinary: Data error. "
            "Number of active dims does not match type.\n");
    goto Finish;
  }
  dp->mStride = head->dp.mStride;
  dp->mNSubField = 0;
  dp->mSubField[0] = dp->mSubField[1] = NULL;
  dp->mSubField[2] = dp->mSubField[3] = NULL;
  //
  //  Now figure out how much data we have, get space, and read it in.
  //
  npoint =  dp->mNVal[0] *  dp->mNVal[1] *  dp->mNVal[2];
  switch (dp->mType) {
    case kCD3Data2:
      npoint *= 2;
      break;
      
    case kCD3Data3:
      npoint *= 3;
      break;
      
    default:
      fprintf(stderr, "CD3ReadBinary:Invalid file type %d.\n", dp->mType);
      goto Finish;
  }
  dp->mField = NULL;
  dp->mField = (double *) malloc(npoint * sizeof(double));
  if (dp->mField == NULL) {
    fprintf(stderr,
            "CD3ReadBinary: Failed to allocate %d doubles for data.\n",
            npoint);
    goto Finish;
  }
  if (fread(dp->mField, sizeof(double), npoint, ifp) != npoint) {
    fprintf(stderr,
            "CD3ReadBinary: Failed to read data.\n");
    goto Finish;
  }
  /*
  for (i = 0; i < npoint/3; i++) {
    printf("{%f,%f,%f}\n", dp->mField[3*i], dp->mField[3*i+1], dp->mField[3*i+2]);
  }*/
  success = true;
  //
  //  All exit paths go through here to clean up.
  //
Finish:
  free(head);
  if (dp->mField == NULL) {
    free(dp->mField);
  }
  return success;
}

//
//  Internal file scope helper functions.
//  We use these first two to complete the initialization process once
//  we know how many dimensions the data have.
//  Note that since these do not malloc till the end they don't need a
//  fancy error exit but can just return an error code.
//
static CDError Init3D(CD3Data* dp, const CDData* cdp)
{
  int dim;
  int nVal;                 // Total number of field points
  double* fieldVals[3];     // Point to the individual field component arrays.
  if (cdp->mNExpression != 3) {
    fprintf(stderr,
            "Expected three expressions, found %d.\n",
            cdp->mNExpression);
    return kCDBadStructure;
  }
  //
  //  To be careful we check that the expression names are in the
  //  correct order.
  //
  if (strcmp(cdp->mExprNames[3], "es.Ex") != 0) {
    fprintf(stderr,
            "First expression '%s' should be 'Ex'.\n",
            cdp->mExprNames[3]);
    return kCDBadStructure;
  }
  if (strcmp(cdp->mExprNames[4], "es.Ey") != 0) {
    fprintf(stderr,
            "Second expression '%s' should be 'Ey'.\n",
            cdp->mExprNames[4]);
    return kCDBadStructure;
  }
  //
  //  Now see that all three dimensions are active and fill in the
  //  MVals, ranges, and data.
  //
  for (dim = 0; dim < 3; dim++) {
    if (cdp->mRange[dim].mNVal < 2) {
      return kCDBadStructure;
    }
    dp->mNVal[dim] = cdp->mRange[dim].mNVal;
    dp->mMin[dim] = cdp->mRange[dim].mMin;
    dp->mMax[dim] = cdp->mRange[dim].mMax;
    dp->mDelta[dim] = cdp->mRange[dim].mDelta;
    fieldVals[dim] = cdp->mDStore[dim + 3];
  }
  //
  //  Now we know all the dimensions we know how much storage we need for
  //  the merged data. Get the space and merge the data arrays into the field.
  //  Because we copy the data we let CDFinish throw the originals away.
  //
  nVal = dp->mNVal[0] * dp->mNVal[1] * dp->mNVal[2];
  dp->mField = (double *) malloc(nVal * 3 * sizeof(double));
  if (dp->mField == NULL) {
    return kCDAllocFailed;
  }
  for (dim = 0; dim < nVal; dim++) {
    dp->mField[dim*3] = fieldVals[0][dim];
    dp->mField[dim*3 + 1] = fieldVals[1][dim];
    dp->mField[dim*3 + 2] = fieldVals[2][dim];
  }
  dp->mType = kCD3Data3;
  dp->mFieldName = cdp->mFileName;
  return kCDNoErr;
}
//
//  2D is very similar.
//  It has to make sure that this is a field that can be an axisymmetric
//  slice of a field that is symmetric about the z axis. That means that
//  the min of both x and y coordinates in the data must be 0.
//
static CDError Init2D(CD3Data* dp, const CDData* cdp)
{
  int dim, nVal;
  uint32_t inactiveDim = -4;
  //
  if (cdp->mNExpression != 2) {
    fprintf(stderr, "Expected two expressions, found %d.\n",
            cdp->mNExpression);
    return kCDBadStructure;
  }
  //
  //  Copy range information while we figure out which dimension is inactive
  //  and thus which count to use as the stride.
  //
  for (dim = 0; dim < 3; dim++) {
    if (!cdp->mRange[dim].mActive) {
      inactiveDim = dim;
    }
    dp->mNVal[dim] = cdp->mRange[dim].mNVal;
    dp->mMin[dim] = cdp->mRange[dim].mMin;
    dp->mMax[dim] = cdp->mRange[dim].mMax;
    dp->mDelta[dim] = cdp->mRange[dim].mDelta;
  }
  if (inactiveDim > 1) {
    fprintf(stderr,
            "Invalid inactive dimension: must be x or y (0 or 1) found %d.\n",
            inactiveDim);
    return kCDBadStructure;
  }
  //
  //  Check x and y min and recompute 3D bounding box.
  //
  if ((dp->mMin[0] != 0) || (dp->mMin[1] != 0)) {
    fprintf(stderr,
            "Loading axisymmetric data, x and y must have min=0.0.");
    return kCDBadStructure;
  }
  //
  //  Now can check the field components.
  //
  switch (inactiveDim) {
    case 0:
      dp->mStride = dp->mNVal[1];
      dp->mDelta[0] = dp->mDelta[1];
      dp->mMin[0] = -dp->mMax[1];
      dp->mMax[0] = dp->mMax[1];
      dp->mMin[1] = -dp->mMax[1];
      //
      //  x inActive. Check for Ey and Ez.
      //
      if (strcmp(cdp->mExprNames[3], "Ey") != 0) {
        fprintf(stderr, "First expression '%s' should be 'Ex'.\n",
                cdp->mExprNames[3]);
        return kCDBadStructure;
      }
      if (strcmp(cdp->mExprNames[4], "Ez") != 0) {
        fprintf(stderr, "Second expression '%s' should be 'Ey'.\n",
                cdp->mExprNames[4]);
        return kCDBadStructure;
      }
      break;
      
    case 1:
      dp->mDelta[1] = dp->mDelta[0];
      dp->mStride = dp->mNVal[0];
      dp->mMin[1] = -dp->mMax[0];
      dp->mMax[1] = dp->mMax[0];
      dp->mMin[0] = -dp->mMax[0];
      //
      //  y inActive. Check for Ex and Ez.
      //
      if (strcmp(cdp->mExprNames[3], "Ex") != 0) {
        fprintf(stderr,
                "First expression '%s' should be 'Ex'.\n",
                cdp->mExprNames[3]);
        return kCDBadStructure;
      }
      if (strcmp(cdp->mExprNames[4], "Ez") != 0) {
        fprintf(stderr,
                "Second expression '%s' should be 'Ey'.\n",
                cdp->mExprNames[4]);
        return kCDBadStructure;
      }
      break;
      
  }
  //
  //  Now we know all the dimensions we know how much storage we need for
  //  the merged data. Get the space and merge the data arrays into the field.
  //  Because we copy the data we let CDFinish throw the originals away.
  //
  nVal = dp->mNVal[0] * dp->mNVal[1] * dp->mNVal[2];
  dp->mField = (double *) malloc(nVal * 2 * sizeof(double));
  if (dp->mField == NULL) {
    return kCDAllocFailed;
  }
  for (dim = 0; dim < nVal; dim++) {
    dp->mField[dim*2] = cdp->mDStore[3][dim];
    dp->mField[dim*2 + 1] = cdp->mDStore[4][dim];
  }
  dp->mType = kCD3Data2;
  dp->mFieldName = cdp->mFileName;
  return kCDNoErr;
}

//
//  Internal Accessors know the structure of the data so they can take
//  it to pieces correctly.
//  They assume that the high level routine has done bounds checking.
//
bool Get3DEAtPoint(const CD3Data* dp, const double coord[3], double* EField)
{
  int index[3], i;
  double rc[3], irc[3];                 // Reduced coords and inverses
  double minc[3];                        // Minima of surrounding box
  int idx000, idx001, idx010, idx011, idx100, idx101, idx110, idx111;
  double c00, c01, c10, c11, c0, c1; // Interpolation steps
//  double x = coord[0], y = coord[1], z = coord[2];
  //  printf("[%f,%f,%f]\n",x,y,z);
  for (i = 0; i < 3; i++) {
    index[i] = (coord[i] - dp->mMin[i]) / dp->mDelta[i];
    if (index[i] == dp->mNVal[i]-1) { // Correct if at top edge
      --index[i];
    }
#ifdef CD3BoundsCheck
    if ((coord[i] < dp->mMin[i]) || (coord[i] > dp->mMax[i])) {
      fprintf(stderr, "Get3DEAtPoint: Coordinate %d, %lf, out of range\n",
              i, coord[i]);
      return false;
    }
    if (index[i] >= dp->mNVal[i]-1) {
      fprintf(stderr, "Get3DEAtPoint: Index %d out of range\n", i);
      return false;
    }
#endif
  }
  //
  //  At this point I have found the indices for the three dimensions. These
  //  are normally the indices of the coord BELOW the given coord. Use this
  //  to compute the array indices (idx's) of the eight points that surround
  //  the cell containing the coordinate.
  //  They have form idx<z><y><x>.
  //  NOTE can't do this before we have all three
  //  indices.
  //
  idx000 = (((index[2])*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
            index[0])*3;
  idx001 = (((index[2])*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
            index[0] + 1)*3;
  idx010 = (((index[2])*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
            index[0])*3;
  idx011 = (((index[2])*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
            index[0] + 1)*3;
  idx100 = (((index[2]+1)*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
            index[0])*3;
  idx101 = (((index[2]+1)*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
            index[0] + 1)*3;
  idx110 = (((index[2]+1)*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
            index[0])*3;
  idx111 = (((index[2]+1)*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
            index[0] + 1)*3;
  //
  //  Next have to find where in each dimension of the box the coordinate is.
  //  This produces a set of reduced coords expressed as a fraction of the
  //  way from the lower to upper point.
  //
  for (i = 0; i < 3; i++) {
    minc[i] = dp->mMin[i] + index[i] * dp->mDelta[i];
    rc[i] = (coord[i] - minc[i])/dp->mDelta[i];
    irc[i] = 1.0 - rc[i];
#ifdef CD3BoundsCheck
    //
    //  Put in a tiny amount of slop for rounding error.
    //
    if ((rc[i] < -0.001) || (rc[i] > 1.001)) {
      fprintf(stderr, "Get3DEAtPoint: Reduced coord %d, %lf, out of range.\n",
              i, rc[i]);
      return false;
    }
#endif
  }
  //
  //  Now we can do the interpolation.
  //  This follows the notation of the Wikipedia page quite closely.
  //  Since we are doing the whole field we have to do this three
  //  times.
  //
  //  Ex
  c00 = irc[0]*dp->mField[idx000 + 0] + rc[0]*dp->mField[idx001 + 0];
  c10 = irc[0]*dp->mField[idx010 + 0] + rc[0]*dp->mField[idx011 + 0];
  c01 = irc[0]*dp->mField[idx100 + 0] + rc[0]*dp->mField[idx101 + 0];
  c11 = irc[0]*dp->mField[idx110 + 0] + rc[0]*dp->mField[idx111 + 0];
  c0 = irc[1] * c00 + rc[1] * c10;
  c1 = irc[1] * c01 + rc[1] * c11;
  EField[0] = irc[2] * c0 + rc[2] * c1;
  //
  //  Ey
  //
  c00 = irc[0]*dp->mField[idx000 + 1] + rc[0]*dp->mField[idx001 + 1];
  c10 = irc[0]*dp->mField[idx010 + 1] + rc[0]*dp->mField[idx011 + 1];
  c01 = irc[0]*dp->mField[idx100 + 1] + rc[0]*dp->mField[idx101 + 1];
  c11 = irc[0]*dp->mField[idx110 + 1] + rc[0]*dp->mField[idx111 + 1];
  c0 = irc[1] * c00 + rc[1] * c10;
  c1 = irc[1] * c01 + rc[1] * c11;
  EField[1] = irc[2] * c0 + rc[2] * c1;
  //
  //  Ez
  //
  c00 = irc[0]*dp->mField[idx000 + 2] + rc[0]*dp->mField[idx001 + 2];
  c10 = irc[0]*dp->mField[idx010 + 2] + rc[0]*dp->mField[idx011 + 2];
  c01 = irc[0]*dp->mField[idx100 + 2] + rc[0]*dp->mField[idx101 + 2];
  c11 = irc[0]*dp->mField[idx110 + 2] + rc[0]*dp->mField[idx111 + 2];
  c0 = irc[1] * c00 + rc[1] * c10;
  c1 = irc[1] * c01 + rc[1] * c11;
  EField[2] = irc[2] * c0 + rc[2] * c1;
  //
  return true;
}

//
//  This treats the field as a defining slice for an axi-symmetric
//  field and returns fully 3D values from 3D points.
//
bool GetAxEAtPoint(const CD3Data* dp, const double coord[3], double* EField)
{
  double coord2D[2];
  double field2D[2];
  double sinval = 0.0, cosval = 0.0;
//  double x = coord[0], y = coord[1], z = coord[2];
  //
  //  First check against 3D bounds.
  //
#ifdef CD3BoundsCheck
  if (!PtInBounds(dp, coord)) {
    fprintf(stderr, "GetEAxAtPoint: Coordinate [%lf, %lf, %lf] out of range\n",
            coord[0], coord[1], coord[2]);
    return false;
    }
#endif
  //
  //  Start by mapping from 3d point to 2D point.
  //
  double r = sqrt(coord[0]*coord[0] + coord[1]*coord[1]);
  if (r > 0.0) {
    sinval = coord[1]/r;
    cosval = coord[0]/r;
  }
  coord2D[0] = r;
  coord2D[1] = coord[2];
  //
  //  Get the 2D field.
  //
  if (!Get2DEAtPoint(dp, coord2D, field2D)) {
    return false;
  }
  //
  //  Map back to 3D.
  //
  EField[0] = field2D[0] * cosval;
  EField[1] = field2D[0] * sinval;
  EField[2] = field2D[1];
  return true;
}
//
//  This uses 2D coords to index array as 2D.
//  It does NOT do coordinate checking because it assumes that a higher
//  level routine has already done that.
//
bool Get2DEAtPoint(const CD3Data* dp, const double coord[2], double* EField)
{
  int index[2], i;
  double rc[2], irc[2];                 // Reduced coords and inverses
  double minc[2];                       // Minima of surrounding box
  int idx00, idx01, idx10, idx11;
  double c0 = coord[0], c1 = coord[1];  // Interpolation steps
  //
  //  Compute indices and range check. Limits for r are 0 and xMax or yMax,
  //  those for z are z so I use y and z.
  //
  index[0] = coord[0] / dp->mDelta[1];
  index[1] = (coord[1] - dp->mMin[2]) / dp->mDelta[2];
  for (i = 0; i < 2; i++) {
    if (index[i] == dp->mNVal[i+1]-1) { // Correct if at top edge
      --index[i];
    }
#ifdef CD3BoundsCheck
    if (index[i] >= dp->mNVal[i+1]) {
      fprintf(stderr, "Get2DEAtPoint: Index #%d = %d out of range\n",
              i, index[i]);
      return nan("");
    }
#endif
  }
  //
  //  At this point I have found the indices for the three dimensions. These
  //  are normallly the indices of the coord BELOW the given coord. Use this
  //  to compute the array indices (idx's) of the eight points that surround
  //  the cell containing the coordinate.
  //  They have form idx<z><y><x>.
  //  NOTE can't do this before we have all three
  //  indices.
  //
  idx00 = ((index[1])*dp->mStride + index[0])*2;
  idx01 = ((index[1])*dp->mStride + index[0] + 1)*2;
  idx10 = ((index[1]+1)*dp->mStride + index[0])*2;
  idx11 = ((index[1]+1)*dp->mStride + index[0] + 1)*2;
  //
  //  Next have to find where in each dimension of the box the coordinate is.
  //  This produces a set of reduced coords expressed as a fraction of the
  //  way from the lower to upper point.
  //  Remember that the 2D bounds are NOT the 3D bounds. We use 0->mMax[0|1].
  //
  minc[0] = index[0] * dp->mDelta[1];
  minc[1] = dp->mMin[2] + index[1] * dp->mDelta[2];
  for (i = 0; i < 2; i++) {
    rc[i] = (coord[i] - minc[i])/dp->mDelta[i+1];
    irc[i] = 1.0 - rc[i];
#ifdef CD3BoundsCheck
    //
    //  Put in a tiny amount of slop for rounding error.
    //
    if ((rc[i] < -0.001) || (rc[i] > 1.001)) {
      fprintf(stderr, "Get2DValueAtPoint:Reduced coord %d, %lf, out of range.\n",
              i, rc[i]);
      return false;
    }
#endif
  }
  //
  //  Now we can do the interpolations.
  //
  c0 = irc[0]*dp->mField[idx00] + rc[0]*dp->mField[idx01];
  c1 = irc[0]*dp->mField[idx10] + rc[0]*dp->mField[idx11];
  EField[0] = irc[1] * c0 + rc[1] * c1;
  c0 = irc[0]*dp->mField[idx00+1] + rc[0]*dp->mField[idx01+1];
  c1 = irc[0]*dp->mField[idx10+1] + rc[0]*dp->mField[idx11+1];
  EField[1] = irc[1] * c0 + rc[1] * c1;
  return true;
  
}

/***********************************************************************
 *
 *  Unused routines just in case.
 *
 ***********************************************************************
 double CD3GetValueAtIndex(const CD3Data* dp, unsigned int dim, const unsigned int index[3])
 {
 if (dim > 2) {
 return nan("");
 }
 if ((index[0] > dp->mNVal[0]) || (index[1] > dp->mNVal[1]) ||
 (index[2] > dp->mNVal[2])) {
 return nan("");
 }
 int idx = (((index[2]*dp->mNVal[1])+index[1])*dp->mNVal[1])+index[0];
 return dp->mField[3*idx + dim];
 
 }
 double CD3GetValueAtPoint(const CD3Data* dp, unsigned int dim, const double coord[3])
 {
 int index[3], i;
 double rc[3], irc[3];                 // Reduced coords and inverses
 double minc[3];                        // Minima of surrounding box
 int idx000, idx001, idx010, idx011, idx100, idx101, idx110, idx111;
 double c00, c01, c10, c11, c0, c1, c; // Interpolation steps
 if (dim > 2) {
 return nan("");
 }
 for (i = 0; i < 3; i++) {
 index[i] = (coord[i] - dp->mMin[i]) / dp->mDelta[i];
 #ifdef CD3BoundsCheck
 if ((coord[i] < dp->mMin[i]) || (coord[i] > dp->mMax[i])) {
 fprintf(stderr, "CD3GetValueAtPoint:Coordinate %d out of range\n", i);
 return nan("");
 }
 if (index[i] >= dp->mNVal[i]) {
 fprintf(stderr, "CD3GetValueAtPoint:Index %d out of range\n", i);
 return nan("");
 }
 #endif
 }
 //
 //  At this point I have found the indices for the three dimensions. These
 //  are normallly the indices of the coord BELOW the given coord. Use this
 //  to compute the array indices (idx's) of the eight points that surround
 //  the cell containing the coordinate.
 //  They have form idx<z><y><x>.
 //  NOTE can't do this before we have all three
 //  indices.
 //
 idx000 = (((index[2])*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
 index[0])*3;
 idx001 = (((index[2])*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
 index[0] + 1)*3;
 idx010 = (((index[2])*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
 index[0])*3;
 idx011 = (((index[2])*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
 index[0] + 1)*3;
 idx100 = (((index[2]+1)*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
 index[0])*3;
 idx101 = (((index[2]+1)*dp->mNVal[1] + (index[1]))*dp->mNVal[0] +
 index[0] + 1)*3;
 idx110 = (((index[2]+1)*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
 index[0])*3;
 idx111 = (((index[2]+1)*dp->mNVal[1] + (index[1]+1))*dp->mNVal[0] +
 index[0] + 1)*3;
 //
 //  Next have to find where in each dimension of the box the coordinate is.
 //  This produces a set of reduced coords expressed as a fraction of the
 //  way from the lower to upper point.
 //
 for (i = 0; i < 3; i++) {
 minc[i] = dp->mMin[i] + index[i] * dp->mDelta[i];
 rc[i] = (coord[i] - minc[i])/dp->mDelta[i];
 irc[i] = 1.0 - rc[i];
 #ifdef CDBoundsCheck
 if ((rc < 0.0) || (rc > 1.0)) {
 fprintf(stderr, "CD3GetValueAtPoint:Reduced coord %d out of range.\n", i);
 return nan("");
 }
 #endif
 }
 //
 //  Now we can do the interpolation.
 //  This follows the notation of the Wikipedia page quite closely.
 //
 c00 = irc[0]*dp->mField[idx000 + dim] + rc[0]*dp->mField[idx001 + dim];
 c10 = irc[0]*dp->mField[idx010 + dim] + rc[0]*dp->mField[idx011 + dim];
 c01 = irc[0]*dp->mField[idx100 + dim] + rc[0]*dp->mField[idx101 + dim];
 c11 = irc[0]*dp->mField[idx110 + dim] + rc[0]*dp->mField[idx111 + dim];
 c0 = irc[1] * c00 + rc[1] * c10;
 c1 = irc[1] * c01 + rc[1] * c11;
 c = irc[2] * c0 + rc[2] * c1;
 return c;
 }
 
 double CD3GetExAtPoint(const CD3Data* dp, const double coord[3])
 {
 return CD3GetValueAtPoint(dp, 0, coord);
 }
 
 double CD3GetEyAtPoint(const CD3Data* dp, const double coord[3])
 {
 return CD3GetValueAtPoint(dp, 1, coord);
 }
 
 double CD3GetEzAtPoint(const CD3Data* dp, const double coord[3])
 {
 return CD3GetValueAtPoint(dp, 2, coord);
 }

 */
/*  Helpers added 7/29/15 */

void CD3ClipPt(const CD3Data* dp, const double coord[3], double newCoord[3])
{
  for (int i = 0; i < 3; i++) {
    newCoord[i] = (coord[i] < dp->mMin[i]) ? dp->mMin[i] :
    ((coord[i] > dp->mMax[i]) ? dp->mMax[i] : coord[i]);
  }
}

//
//  Map checks whether a point is in range and if it is fills in the
//  indices of the corresponding entry and returns true, else false.
//
bool CD3Map(const CD3Data* dp, const double coord[3], uint32_t newIndices[3])
{
  if (PtInBounds(dp, coord)) {
    for (int i = 0; i < 3; i++) {
      newIndices[i] = (uint32_t) ((coord[i] - dp->mMin[i]) / dp->mDelta[i] + 0.5);
    }
    return true;
  }
  return false;
}
//
//  This translates an index trio into a single index.
//
uint64_t CD3IndexAt(const CD3Data* dp, uint32_t ix, uint32_t iy, uint32_t iz)
{
  //
  //  Do this step by step to avoid overflow.
  //
  uint64_t index = iz;
  index = index * dp->mNVal[1] + iy;
  index = index * dp->mNVal[0] + ix;
  return index;
}


