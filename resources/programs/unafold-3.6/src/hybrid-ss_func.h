double Eh(int i, int j)
{
  double energy;
  int loopSize = j - i - 1;
  int k;

  if (loopSize < TURN)
    return 0.0;

  if (i <= g_len && g_len < j)
    return 0.0;
  else if (i > g_len)
    {
      i -= g_len;
      j -= g_len;
    }

#if ENABLE_FORCE
  if (!ssOK(i + 1, j - 1))
    return 0.0;
#endif

  if (loopSize <= 30)
    energy = g_hairpinLoop[loopSize - 1];
  else
    energy = g_hairpinLoop[29] * pow(g_misc[12], log((double) loopSize / 30)) / g_scalen[loopSize - 30];

  if (loopSize > 3)
    energy *= g_tstackh[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
  else
    energy *= auPenalty(i, j);

  if (loopSize == 3)
    {
      struct triloop* loop;
      if (numTriloops)
	if ((loop = bsearch(g_seq + i, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
	  energy *= loop->energy;
    }
  else if (loopSize == 4)
    {
      struct tloop* loop;
      if (numTloops)
	if ((loop = bsearch(g_seq + i, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
	  energy *= loop->energy;
    }
  else if (loopSize == 6)
    {
      struct hexaloop* loop;
      if (numHexaloops)
	if ((loop = bsearch(g_seq + i, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
	  energy *= loop->energy;
    }

  /* GGG */
  if (i >= 3 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[i] == 2 && g_seq[j] == 3)
    energy *= g_misc[8];

  /* poly-C */
  if (loopSize == 3 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1 && g_seq[i + 3] == 1)
    energy *= g_misc[11];
  else
    {
      for (k = 1; k <= loopSize; ++k)
	if (g_seq[i + k] != 1)
	  return energy;
      energy *= pow(g_misc[9], loopSize) * g_misc[10];
    }

  return energy;
}

double Es(int i, int j)
{
#ifdef DEBUG
  if (i >= j)
    fputs("Error in Es(): i isn't less than j\n", stderr);
#endif

  if (i == g_len || j == g_len + 1)
    return 0.0;

  if (i > g_len)
    i -= g_len;
  if (j > g_len)
    j -= g_len;

  return g_stack[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
}

#ifdef REDUCED2
double Ebi_old(int, int, int, int);

double Ebi(int i, int j, int ii, int jj)
{
  double value;

  value = Ebi_old(i, j, ii, jj);
  if (ii - i == 2 && j - jj == 2)
    {
      if (value < Es(i, j) * Es(i + 1, j - 1))
	return 0.0;
    }
  else
    {
      if (value < Es(i, j) * Ebi_old(i + 1, j - 1, ii, jj))
	return 0.0;
      if (value < Ebi_old(i, j, ii - 1, jj + 1) * Es(ii - 1, jj + 1))
	return 0.0;
    }
  return value;
}

double Ebi_old(int i, int j, int ii, int jj)
#else
double Ebi(int i, int j, int ii, int jj)
#endif
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
  if (ii >= jj)
    fputs("Error in Ebi(): jj isn't greater than ii\n", stderr);

  if ((i <= g_len && g_len < ii) || (jj <= g_len && g_len < j))
    return 0.0;
#endif

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > g_maxLoop)
    return 0.0;

#ifdef DEBUG
  if (i > g_len)
    i -= g_len;
  if (ii > g_len)
    ii -= g_len;
  if (j > g_len)
    j -= g_len;
  if (jj > g_len)
    jj -= g_len;
#endif

#if ENABLE_FORCE
  if (loopSize1 && !ssOK(i + 1, ii - 1))
    return 0.0;
  if (loopSize2 && !ssOK(jj + 1, j - 1))
    return 0.0;
#endif

#ifdef REDUCED_INTERIOR
  if (loopSize1 && loopSize2)
    {
      if (basePairIndex(g_seq[i + 1], g_seq[j - 1]) < 6)
	return 0.0;
      if (basePairIndex(g_seq[ii - 1], g_seq[jj + 1]) < 6)
	return 0.0;
    }
#endif

#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Ebi() called with nonsense\n", stderr);
      return 1;
    }
  else
#endif
  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return g_bulgeLoop[0] * g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]] * g_scale * g_scale;
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] * auPenalty(i, j) * auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] * pow(g_misc[12], log((double) loopSize2 / 30)) / g_scalen[loopSize2 - 30] * auPenalty(i, j) * auPenalty(ii, jj);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] * g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]] * g_scale * g_scale;
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] * auPenalty(i, j) * auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] * pow(g_misc[12], log((double) loopSize1 / 30)) / g_scalen[loopSize1 - 30] * auPenalty(i, j) * auPenalty(ii, jj);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(g_seq[jj], g_seq[ii])][basePairIndex(g_seq[j], g_seq[i])][g_seq[jj + 1]][g_seq[ii - 1]][g_seq[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[i + 2]][g_seq[j - 2]];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] * pow(g_misc[12], log((double) (loopSize1 + loopSize2) / 30)) / g_scalen[loopSize1 + loopSize2 - 30];
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy *= g_tstacki[g_seq[i]][g_seq[j]][0][0];
	  loopEnergy *= g_tstacki[g_seq[jj]][g_seq[ii]][0][0];
	}
      else
	{
	  loopEnergy *= g_tstacki[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
	  loopEnergy *= g_tstacki[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy *= exp(-asPenalty / RT);

      return loopEnergy;
    }
}

void setDangle5(int i, int* upst, int* dnst)
{
  upst[i - 1] = i - 1;
  dnst[i - 2] = i;
}

void setDangle3(int j, int* upst, int* dnst)
{
  upst[j] = j;
  dnst[j - 1] = j + 1;
}

void setBI(int i, int j, int ii, int jj, int* upst, int* dnst)
{
  int loopSize1, loopSize2;

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;

#ifdef DEBUG
  if (loopSize1 < 0 || loopSize2 < 0 || (loopSize1 == 0 && loopSize2 == 0))
    {
      fputs("Error: setBI() called with nonsense\n", stderr);
      return;
    }
  else
#endif

  if ((loopSize1 == 0 && loopSize2 == 1) || (loopSize2 == 0 && loopSize1 == 1))
    {
      upst[ii - 1] = i;
      dnst[i - 1] = ii;
      upst[j - 1] = jj;
      dnst[jj - 1] = j;
    }
  else if (loopSize1 && loopSize2 && (loopSize1 > 2 || loopSize2 > 2))
    {
      upst[i] = i;
      upst[ii - 1] = ii - 1;
      dnst[i - 1] = i + 1;
      dnst[ii - 2] = ii;
      upst[jj] = jj;
      upst[j - 1] = j - 1;
      dnst[jj - 1] = jj + 1;
      dnst[j - 2] = j;
    }
}

double QBI(int i, int j)
{
  int d, ii, jj;
  double energy = 0.0;

  for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop; --d)
    for (ii = i + 1; ii < j - d && ii <= g_len; ++ii)
      {
	jj = d + ii;
	if (Qprime(ii, jj) != 0.0)
	  energy += Ebi(i, j, ii, jj) * Qprime(ii, jj);
      }

  return energy;
}

double QBI2(int i, int j)
{
  int d, ii, jj;
  double energy = 0.0;

  for (d = j - i - 3; d >= 1 && d >= j - i - 2 - g_maxLoop; --d)
    for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii <= g_len; ++ii)
      {
	jj = d + ii;
	if (Qprime(ii, jj) != 0.0)
	  energy += Ebi(jj - g_len, ii, j - g_len, i) * Qprime(ii, jj);
      }

  return energy;
}

double QBI_noI(int i, int j)
{
  int d, ii, jj;
  double energy = 0.0;

  for (d = j - i - 3; d >= TURN + 3 && d >= j - i - 2 - g_maxLoop; --d)
    for (ii = i + 1; ii < j - d && ii < g_len; ++ii)
      {
	jj = d + ii;
	if (Qprime(ii, jj) != 0.0)
	  energy += Ebi(i, j, ii, jj) * Es(ii, jj) * Qprime(ii + 1, jj - 1);
      }

  return energy;
}

double QBI2_noI(int i, int j)
{
  int d, ii, jj;
  double energy = 0.0;

  for (d = j - i - 3; d >= 3 && d >= j - i - 2 - g_maxLoop; --d)
    for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii < g_len; ++ii)
      {
	jj = d + ii;
	if (Qprime(ii, jj) != 0.0)
	  energy += Ebi(jj - g_len, ii, j - g_len, i) * Es(ii, jj) * Qprime(ii + 1, jj - 1);
      }

  return energy;
}

double* calloc2(int n)
{
  return xcalloc(n * (n - 1) / 2, sizeof(double));
}

double* calloc2_double(int n)
{
  return xcalloc(n * n, sizeof(double));
}
