#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>


#define SGCOREDEF__
#include "sginfo.h"


int str_icmp(char *s, char *t)
{
  char     cs, ct;

  while (*s || *t)
  { cs = toupper(*s++);
    ct = toupper(*t++);
    if (cs < ct) return -1;
    if (cs > ct) return  1;
  }
  return 0;
}


static char *SgInfoCpy(T_SgInfo *Target, T_SgInfo *Source)
{
  T_RTMx       *ListSeitzMx;
  T_RotMxInfo  *ListRotMxInfo;


  ListSeitzMx   = Target->ListSeitzMx;
  ListRotMxInfo = Target->ListRotMxInfo;

  memcpy(Target, Source, sizeof (*Target));

  Target->ListSeitzMx   = ListSeitzMx;
  Target->ListRotMxInfo = ListRotMxInfo;

  memcpy(Target->ListSeitzMx, Source->ListSeitzMx,
         Source->nList * sizeof (*(Target->ListSeitzMx)));

  if (Target->ListRotMxInfo != NULL)
  {
    if (Source->ListRotMxInfo == NULL)
      return "Internal Error: SgInfoCpy()";

    memcpy(Target->ListRotMxInfo, Source->ListRotMxInfo,
           Source->nList * sizeof (*(Target->ListRotMxInfo)));
  }

  return NULL;
}


static const char *progn = "sgtest";


static void progerror(const char *message)
{
  fflush(stdout);
  fprintf(stderr, "%s: %s\n", progn, message);
  exit(1);
}


static void NotEnoughCore(void)
{
  progerror("Not enough core");
}


static void usage(void)
{
  fprintf(stderr, "usage: %s [-v] [-light] [#Start [#End]]\n", progn);
  exit (1);
}


int main(int argc, char *argv[])
{
  int                F_Verbose, F_light, SgNumStart, SgNumEnd;
  int                i, n, pos_hsym, OrSh[3];
  int                nNextBasis, iNextBasis;
  T_SgInfo           SpgrInfo0, SpgrInfoN, SpgrInfo, BC_SgInfo;
  const T_TabSgName  *TSgN, *ConvTSgN;
  T_RTMx             CBMx, InvCBMx;


  F_Verbose = 0;
  F_light   = 0;

  SgNumStart = -1;
  SgNumEnd   = -1;

  for (i = 1; i < argc; i++)
  {
    if      (str_icmp(argv[i], "-v") == 0)
      F_Verbose = 1;
    else if (str_icmp(argv[i], "-light") == 0)
      F_light   = 1;
    else if (SgNumStart == -1)
    {
          n = sscanf(argv[i], "%d", &SgNumStart);
      if (n != 1 || SgNumStart < 1) usage();
    }
    else if (SgNumEnd   == -1)
    {
          n = sscanf(argv[i], "%d", &SgNumEnd);
      if (n != 1 || SgNumEnd < 1) usage();
    }
    else
      usage();
  }

  if (SgNumEnd < 0) {
    if (SgNumStart < 0) {
      SgNumStart =   1;
      SgNumEnd   = 230;
    }
    else
      SgNumEnd   = SgNumStart;
  }

  SpgrInfo0.MaxList = 192;

  SpgrInfo0.ListSeitzMx
    = malloc(SpgrInfo0.MaxList * sizeof (*SpgrInfo0.ListSeitzMx));
  if (SpgrInfo0.ListSeitzMx == NULL) NotEnoughCore();

  SpgrInfo0.ListRotMxInfo
    = malloc(SpgrInfo0.MaxList * sizeof (*SpgrInfo0.ListRotMxInfo));
  if (SpgrInfo0.ListRotMxInfo == NULL) NotEnoughCore();

  SpgrInfoN.MaxList = 192;

  SpgrInfoN.ListSeitzMx
    = malloc(SpgrInfoN.MaxList * sizeof (*SpgrInfoN.ListSeitzMx));
  if (SpgrInfoN.ListSeitzMx == NULL) NotEnoughCore();

  SpgrInfoN.ListRotMxInfo
    = malloc(SpgrInfoN.MaxList * sizeof (*SpgrInfoN.ListRotMxInfo));
  if (SpgrInfoN.ListRotMxInfo == NULL) NotEnoughCore();

  SpgrInfo.MaxList = 192;

  SpgrInfo.ListSeitzMx
    = malloc(SpgrInfo.MaxList * sizeof (*SpgrInfo.ListSeitzMx));
  if (SpgrInfo.ListSeitzMx == NULL) NotEnoughCore();

  SpgrInfo.ListRotMxInfo
    = malloc(SpgrInfo.MaxList * sizeof (*SpgrInfo.ListRotMxInfo));
  if (SpgrInfo.ListRotMxInfo == NULL) NotEnoughCore();

  BC_SgInfo.MaxList = 192;

  BC_SgInfo.ListSeitzMx
    = malloc(BC_SgInfo.MaxList * sizeof (*BC_SgInfo.ListSeitzMx));
  if (BC_SgInfo.ListSeitzMx == NULL) NotEnoughCore();

  BC_SgInfo.ListRotMxInfo
    = malloc(BC_SgInfo.MaxList * sizeof (*BC_SgInfo.ListRotMxInfo));
  if (BC_SgInfo.ListRotMxInfo == NULL) NotEnoughCore();

  for (TSgN = TabSgName; TSgN->HallSymbol; TSgN++)
  {
    if (TSgN->SgNumber < SgNumStart)
      continue;

    if (TSgN->SgNumber > SgNumEnd)
      break;

    InitSgInfo(&SpgrInfo0);

    pos_hsym = ParseHallSymbol(TSgN->HallSymbol, &SpgrInfo0);

    if (SgError != NULL)
    {
      fprintf(stdout, "    %s\n", TSgN->HallSymbol);
      for (i = 0; i < pos_hsym; i++) putc('-', stdout);
      fprintf(stdout, "---^\n");
      fprintf(stdout, "%s\n", SgError);
      SgError = NULL;
      continue;
    }
                     i = VolAPointGroups[TSgN->SgNumber];
    if (   PG_Number(i) >= PG_Number(PG_4)
        && PG_Number(i) <= PG_Number(PG_6_mmm)
        && TSgN->Extension[0] != 'R')
      nNextBasis = 3;
    else
      nNextBasis = 1;

    for (iNextBasis = 0; iNextBasis < nNextBasis; iNextBasis++)
    {
      fprintf(stdout, "TSgN  ");
      PrintTabSgNameEntry(TSgN, 0, 0, stdout);

      i = PG_Index(VolAPointGroups[TSgN->SgNumber]);

      fprintf(stdout, " : %s", PG_Names[i]);
      fprintf(stdout, " : %s", PG_Names[PG_Index(LG_Code_of_PG_Index[i])]);

      if      (iNextBasis == 0)
      {
        if (nNextBasis != 1)
          fprintf(stdout, " : z->z");

        SgError = SgInfoCpy(&SpgrInfoN, &SpgrInfo0);
      }
      else if (iNextBasis == 1)
      {
        fprintf(stdout, " : z->x");

        for (i = 0; i < 9;  i++) {
             CBMx.s.R[i] = CRBF * RMx_3_111[i];
          InvCBMx.s.R[i] = CRBF * RMx_3i111[i];
        }

        for (i = 0; i < 3; i++) {
             CBMx.s.T[i] = 0;
          InvCBMx.s.T[i] = 0;
        }

        InitSgInfo(&SpgrInfoN);
        TransformSgInfo(&SpgrInfo0, &CBMx, &InvCBMx, &SpgrInfoN);
      }
      else if (iNextBasis == 2)
      {
        fprintf(stdout, " : z->y");

        for (i = 0; i < 9;  i++) {
             CBMx.s.R[i] = CRBF * RMx_3i111[i];
          InvCBMx.s.R[i] = CRBF * RMx_3_111[i];
        }

        for (i = 0; i < 3; i++) {
             CBMx.s.T[i] = 0;
          InvCBMx.s.T[i] = 0;
        }

        InitSgInfo(&SpgrInfoN);
        TransformSgInfo(&SpgrInfo0, &CBMx, &InvCBMx, &SpgrInfoN);
      }

      putc('\n', stdout);
      fflush(stdout);

          SgError = SgInfoCpy(&SpgrInfo, &SpgrInfoN);
      if (SgError)
      {
        fprintf(stdout, "@ (%d) %s => %s\n",
          TSgN->SgNumber, TSgN->HallSymbol, SgError);

        SgError = NULL;
        goto Next_TSgN;
      }

      for (OrSh[0] = 0; OrSh[0] < 12; OrSh[0]++)
      for (OrSh[1] = 0; OrSh[1] < 12; OrSh[1]++)
      for (OrSh[2] = 0; OrSh[2] < 12; OrSh[2]++)
      {
            SgError = SgInfoCpy(&SpgrInfo, &SpgrInfoN);
        if (SgError)
        {
          fprintf(stdout, "@ (%d) %s => %s\n",
            TSgN->SgNumber, TSgN->HallSymbol, SgError);

          exit(1);
        }

        if (F_light == 0)
          for (i = 0; i < 3; i++)
            SpgrInfo.OriginShift[i] = OrSh[i];

        if (CompleteSgInfo(&SpgrInfo) != 0)
        {
          if (SgError != NULL)
          {
            fprintf(stdout, "@ (%d) %s (%d %d %d) => %s\n",
              TSgN->SgNumber, TSgN->HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SgError);

            SgError = NULL;
            goto Next_TSgN;
          }
        }

        if (F_Verbose)
          fprintf(stdout, "Hall Symbol  %s\n", SpgrInfo.HallSymbol);

        ConvTSgN = FindReferenceSpaceGroup(&SpgrInfo, &CBMx, &InvCBMx);

        if (SgError)
        {
          fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => %s\n",
            TSgN->SgNumber, TSgN->HallSymbol,
            OrSh[0], OrSh[1], OrSh[2],
            SpgrInfo.HallSymbol, SgError);

          SgError = NULL;
          goto Next_TSgN;
        }
        else
        {
          if (F_Verbose) {
            PrintMapleRTMx(   &CBMx, CRBF, CTBF,    "CBMx", stdout);
            PrintMapleRTMx(&InvCBMx, CRBF, CTBF, "InvCBMx", stdout);
          }

          InitSgInfo(&BC_SgInfo);

              i = TransformSgInfo(&SpgrInfo, &CBMx, &InvCBMx, &BC_SgInfo);
          if (i == 0)
              i = CompleteSgInfo(&BC_SgInfo);
          if (i != 0) {
            fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => (%d) %s => %s\n",
              TSgN->SgNumber, TSgN->HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SpgrInfo.HallSymbol,
              ConvTSgN->SgNumber, ConvTSgN->HallSymbol,
              SgError);

            SgError = NULL;
            goto Next_TSgN;
          }
          else if (BC_SgInfo.TabSgName == NULL) {
            fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => (%d) %s => %s\n",
              TSgN->SgNumber, TSgN->HallSymbol,
              OrSh[0], OrSh[1], OrSh[2],
              SpgrInfo.HallSymbol,
              ConvTSgN->SgNumber, ConvTSgN->HallSymbol,
              BC_SgInfo.HallSymbol);

            goto Next_TSgN;
          }
          else
          {
            if (   BC_SgInfo.TabSgName->SgNumber != ConvTSgN->SgNumber
                || BC_SgInfo.TabSgName->SgNumber != TSgN->SgNumber)
            {
              fprintf(stdout, "@ (%d) %s (%d %d %d) => %s => (%d) %s => ",
                TSgN->SgNumber, TSgN->HallSymbol,
                OrSh[0], OrSh[1], OrSh[2],
                SpgrInfo.HallSymbol,
                ConvTSgN->SgNumber, ConvTSgN->HallSymbol);

              PrintTabSgNameEntry(BC_SgInfo.TabSgName, 0, 0, stdout);
              putc('\n', stdout);

              goto Next_TSgN;
            }

            if (F_Verbose)
            {
              fprintf(stdout, "Change of Basis => ");
              PrintTabSgNameEntry(BC_SgInfo.TabSgName, 0, 0, stdout);
              putc('\n', stdout);
            }
          }
        }

        if (F_light)
          goto Next_iNextBase;
      }

      Next_iNextBase:;
    }

    Next_TSgN:;
  }

  return 0;
}
