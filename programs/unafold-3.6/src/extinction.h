#ifndef EXTINCTION_H
#define EXTINCTION_H

double g_extinction[20];

void loadExtinctionDat(int NA)
{
  int i;
  FILE* file;
  char* buffer;

  file = fopen("extinction.DAT", "rt");

  if (!file && getenv("UNAFOLDDAT"))
    {
      buffer = xmalloc(strlen(getenv("UNAFOLDDAT")) + 16);
      strcpy(buffer, getenv("UNAFOLDDAT"));
      strcat(buffer, "/");
      strcat(buffer, "extinction.DAT");
      file = fopen(buffer, "rt");
      free(buffer);
    }

  if (!file)
    {
      buffer = xmalloc(strlen(PKGDATADIR) + 15);
      strcpy(buffer, PKGDATADIR);
      strcat(buffer, "extinction.DAT");
      if (!(file = fopen(buffer, "rt")))
	{
	  perror("extinction.DAT");
	  exit(EXIT_FAILURE);
	}
      free(buffer);
    }

  if (NA)
    for (i = 0; i < 20; ++i)
      fscanf(file, "%*g%lg", &g_extinction[i]);
  else
    for (i = 0; i < 20; ++i)
      fscanf(file, "%lg%*g", &g_extinction[i]);

  fclose(file);
}

double xi1(unsigned char base)
{
  if (base == 4)
    return 0.0;
  return g_extinction[base + 16];
}

double xi2(unsigned char base1, unsigned char base2)
{
  /* if one base is unknown, we just take the absorbance of the other base */
  if (base1 == 4)
    return xi1(base2);
  if (base2 == 4)
    return xi1(base1);
  return g_extinction[4 * base1 + base2];
}

#endif /* EXTINCTION_H */
