#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#include "sginfo.h"


static int str_ibegin(const char *s1, const char *s2) /* string ignore-case */
{                                                     /* begin              */
  char     u1, u2;

  while (*s1 && *s2)
  {
    u1 = toupper(*s1++);
    u2 = toupper(*s2++);
    if      (u1 < u2) return -1;
    else if (u1 > u2) return  1;
  }
  if (*s2) return -1;
  return 0;
}


int BuildSgInfo(T_SgInfo *SgInfo, const char *SgName)
{
  int                VolLetter;
  const T_TabSgName  *tsgn;


  /* look for "VolA", "VolI", or "Hall"
   */

  while (*SgName && isspace(*SgName)) SgName++;

  VolLetter = -1;

  if      (isdigit(*SgName))
    VolLetter = 'A';
  else if (str_ibegin(SgName, "VolA") == 0)
  {
    VolLetter = 'A';
    SgName += 4;
  }
  else if (   str_ibegin(SgName, "VolI") == 0
           || str_ibegin(SgName, "Vol1") == 0)
  {
    VolLetter = 'I';
    SgName += 4;
  }
  else if (str_ibegin(SgName, "Hall") == 0)
  {
    VolLetter = 0;
    SgName += 4;
  }

  while (*SgName && isspace(*SgName)) SgName++;

  /* default is "VolA"
   */

  if (VolLetter == -1)
    VolLetter = 'A';

  /* if we don't have a Hall symbol do a table look-up
   */

  tsgn = NULL;

  if (VolLetter)
  {
    tsgn = FindTabSgNameEntry(SgName, VolLetter);
    if (tsgn == NULL) return -1; /* no matching table entry */
    SgName = tsgn->HallSymbol;
  }

  /* Allocate memory for the list of Seitz matrices and
     a supporting list which holds the characteristics of
     the rotation parts of the Seitz matrices
   */

  SgInfo->MaxList = 192; /* absolute maximum number of symops */

  SgInfo->ListSeitzMx
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));

  if (SgInfo->ListSeitzMx == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  SgInfo->ListRotMxInfo
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));

  if (SgInfo->ListRotMxInfo == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  /* Initialize the SgInfo structure
   */

  InitSgInfo(SgInfo);
  SgInfo->TabSgName = tsgn; /* in case we know the table entry */

  /* Translate the Hall symbol and generate the whole group
   */

  ParseHallSymbol(SgName, SgInfo);
  if (SgError != NULL) return -1;

  /* Do some book-keeping and derive crystal system, point group,
     and - if not already set - find the entry in the internal
     table of space group symbols
   */

  return CompleteSgInfo(SgInfo);
}


int main(int argc, char *argv[])
{
  T_SgInfo  SgInfo;


  if (argc == 2)
  {
    if (BuildSgInfo(&SgInfo, argv[1]) != 0)
      fprintf(stderr, "%s\n", SgError);
    else
    {
      ListSgInfo(&SgInfo, 1, 0, stdout);

      if (SgError)
        fprintf(stderr, "%s\n", SgError);
    }
  }

  return 0;
}
