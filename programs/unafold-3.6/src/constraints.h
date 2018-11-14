  if (constraints)
    {
      char c;
      FILE* aux;
      if (constraintsFile)
	{
	  if (!(aux = fopen(constraintsFile, "rt")))
	    {
	      perror(constraintsFile);
	      return EXIT_FAILURE;
	    }
	}
      else
	{
	  buffer = xmalloc(strlen(g_prefix) + 5);
	  sprintf(buffer, "%s.aux", g_prefix);
	  if (!(aux = fopen(buffer, "rt")))
	    {
	      perror(buffer);
	      return EXIT_FAILURE;
	    }
	  free(buffer);
	}

      while (fscanf(aux, " %c", &c) == 1)
	{
	  newTop = xmalloc(sizeof(struct constraintListNode));
	  newTop->l = 0;
	  if ((c == 'R' && fscanf(aux, "%d-%d %d-%d", &newTop->i, &newTop->j, &newTop->k, &newTop->l) != 4) ||
	      (c != 'R' && fscanf(aux, "%d%d%d", &newTop->i, &newTop->j, &newTop->k) != 3))
	    {
	      free(newTop);
	      break;
	    }
	  if (c == 'P' || c == 'R')
	    {
	      newTop->next = prohibitList;
	      prohibitList = newTop;
	    }
	  else if (c == 'F')
	    {
	      newTop->next = forceList;
	      forceList = newTop;
	    }
	  else
	    free(newTop);
	}

      fclose(aux);
    }
