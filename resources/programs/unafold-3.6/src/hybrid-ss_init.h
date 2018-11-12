/* int helixLength(int i, int j)
{
  int k, length;

  if (Qprime(i, j) == 0.0)
    return 0;

  length = 1;
  for (k = 1; i + k < j - k && Qprime(i + k, j - k) != 0.0; ++k);
  length += k - 1;
  for (k = 1; i > k && j + k <= g_len && Qprime(i - k, j + k) != 0.0; ++k);
  length += k - 1;

  return length;
} */

void prefilter()
{
  char** in;
  int i, j, k, count;

  in = xcalloc(g_len, sizeof(char*));
  for (i = 1; i <= g_len; ++i)
    in[i - 1] = xcalloc(g_len, 1);

  for (i = 1; i <= g_len - g_prefilter2 + 1; ++i)
    for (j = g_len; j >= g_prefilter2 && j >= i; --j)
      {
	count = 0;
	for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	  if (Qprime(i + k, j - k) != 0.0)
	    ++count;
	if (count >= g_prefilter1)
	  for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	    ++in[i + k - 1][j - k - 1];
      }

  for (i = 1; i <= g_len; ++i)
    {
      for (j = g_len; j >= i; --j)
	if (!in[i - 1][j - 1])
	  Qprime(i, j) = 0.0;
      free(in[i - 1]);
    }
  free(in);
}

void initializeMatrices()
{
  int i, j, k;
  struct constraintListNode *top, *newTop;

  /* Q' is initialized to 0 iff base pair is illegal; 1 otherwise
     Q and Q1 are always initialized to 0 */
  for (i = 1; i <= g_len; ++i)
    for (j = i; j <= g_len; ++j)
      if (j - i < TURN + 1 || (basePairIndex(g_seq[i], g_seq[j]) == 6 && !g_allPairs))
	Q(i, j) = Qprime(i, j) = Q1(i, j) = 0.0;
      else if (j - i > g_maxBP)
	Q(i, j) = Qprime(i, j) = Q1(i, j) = 0.0;
      else
	{
	  Q(i, j) = Q1(i, j) = 0.0;
	  Qprime(i, j) = 1.0;
	}

  if (g_bpFile)
    {
      FILE* bp;

      for (i = 1; i <= g_len; ++i)
	for (j = i; j <= g_len; ++j)
	  Qprime(i, j) = 0.0;
    
      if (!(bp = fopen(g_bpFile, "rt")))
	{
	  perror(g_bpFile);
	  exit(EXIT_FAILURE);
	}

      while (fscanf(bp, "%d%d%d", &i, &j, &k) == 3)
	for (--k; k >= 0; --k)
	  Qprime(i + k, j - k) = 1.0;

      fclose(bp);
    }

  top = prohibitList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len &&
	  top->k >= 1 && top->k <= g_len && top->l >= 1 && top->l <= g_len)
	for (i = top->i; i <= top->j; ++i)
	  for (j = top->k; j <= top->l; ++j)
	    {
	      if (i <= j)
		Qprime(i, j) = 0.0;
	      else
		Qprime(j, i) = 0.0;
	    }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len)
	for (k = 0; k < top->k; ++k)
	  {
	    if (top->i + k <= top->j - k)
	      Qprime(top->i + k, top->j - k) = 0.0;
	    else
	      Qprime(top->j - k, top->i + k) = 0.0;
	  }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j == 0)
	for (k = 0; k < top->k; ++k)
	  {
	    for (j = 1; j <= top->i + k; ++j)
	      Qprime(j, top->i + k) = 0.0;
	    for (j = top->i + k; j <= g_len; ++j)
	      Qprime(top->i + k, j) = 0.0;
	  }
      else if (top->l == 0 && top->j >= 1 && top->j <= g_len && top->i == 0)
	for (k = 0; k < top->k; ++k)
	  {
	    for (i = 1; i <= top->j + k; ++i)
	      Qprime(i, top->j + k) = 0.0;
	    for (i = top->j + k; i <= g_len; ++i)
	      Qprime(top->j + k, i) = 0.0;
	  }

      newTop = top->next;
      free(top);
      top = newTop;
    }

#if ENABLE_FORCE
  for (i = 0; i <= g_len + 1; ++i)
    for (j = 0; j <= g_len + 1; ++j)
      ssOK(i, j) = 1;
  top = forceList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len)
	for (i = 0; i <= g_len + 1; ++i)
	  for (j = i; j <= g_len + 1; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (i <= top->i + k && top->i + k <= j)
		ssOK(i, j) = 0;

      if (top->j >= 1 && top->j <= g_len)
	{
	  if (top->i == 0)
	    {
	      for (i = 0; i <= g_len + 1; ++i)
		for (j = i; j <= g_len + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j + k && top->j + k <= j)
		      ssOK(i, j) = 0;
	    }
	  else
	    {
	      for (i = 0; i <= g_len + 1; ++i)
		for (j = i; j <= g_len + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j - k && top->j - k <= j)
		      ssOK(i, j) = 0;
	    }
	}

      if (top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len)
	{
	  for (i = 1; i <= g_len; ++i)
	    for (k = 0; k < top->k; ++k)
	      if (i != top->i + k && i <= top->j - k)
		/* Qprime(i, top->j - k) = Qprime(top->j - k, i) = 0.0; */
		Qprime(i, top->j - k) = 0.0;
	  for (j = 1; j <= g_len; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (j != top->j - k && top->i + k <= j)
		/* Qprime(top->i + k, j) = Qprime(j, top->i + k) = 0.0; */
		Qprime(top->i + k, j) = 0.0;
	}

      newTop = top->next;
      free(top);
      top = newTop;
    }
#endif

  /* if (g_prefilter && !g_allPairs)
    for (i = 1; i <= g_len; ++i) 
      for (j = i; j <= g_len; ++j)
	if (helixLength(i, j) <= g_prefilter)
	Qprime(i, j) = 0.0;*/
 
  prefilter();

  for (i = 1; i <= g_len; ++i)
    for (j = g_len + 1; j < i + g_len; ++j)
      if (Qprime(j - g_len, i) == 0.0)
	Q(i, j) = Q1(i, j) = Qprime(i, j) = 0.0;
      else
	{
	  Q(i, j) = Q1(i, j) = 0.0;
	  Qprime(i, j) = 1.0;
	}
}
