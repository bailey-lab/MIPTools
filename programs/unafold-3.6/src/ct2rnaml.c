#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#define COMMENT "Output from UNAFold"
#define DTD_DESC "-//University of Montreal//RNAML 1.1//EN"
#define DTD_URL "http://www-lbit.iro.umontreal.ca/rnaml/current/rnaml.dtd"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#if HAVE_STRING_H
# include <string.h>
#endif

#include "xmalloc.h"

struct ctLine
{
  int seqnum;
  char base;
  int previous;
  int next;
  int paired;
  int historic;
  int upst;
  int dnst;
  int strand;
};

struct ctHead
{
  int num_bases;
  int has_dG;
  int has_initial;
  float dG_val;
  float initial_val;
  char name[250];
  int is_rna, is_dna;
  int num_strands;
};

struct ssLine
{
  int seqnum;
  char base;
  double x;
  double y;
  int paired;
};

int readCtFile(FILE*, struct ctHead*, struct ctLine**);
int isValidCt(struct ctHead*, struct ctLine*);
int readSsFile(FILE*, struct ssLine*, int);
int isMatch(int, struct ctLine*, struct ssLine*);
int diff(int, struct ctLine*, struct ctLine*);
void makeValidId(char*);

int main(int argc, char** argv)
{
  char buffer[100];
  struct ctHead head, oldHead;
  int seqnew_bool, seq_num;
  int i, model_num, strand, strandBegin, strandEnd, helixLength;
  struct ctLine *lines, *oldLines;
  struct ssLine* ssLines;
  FILE *ctFile, *ssFile;
  FILE *output;
  
  if (argc == 1)
    {
      printf("Usage: %s file.ct\n", argv[0]);
      return 0;
    }

  if (!strcmp(argv[1] + strlen(argv[1]) - 3, ".ct"))
    argv[1][strlen(argv[1]) - 3] = 0;

  strcpy(buffer, argv[1]);
  strcat(buffer, ".ct");
  if (!(ctFile = fopen(buffer, "rt")))
    {
      perror(buffer);
      return 1;
    }

  strcpy(buffer, argv[1]);
  strcat(buffer, ".rnaml");
  if (!(output = fopen(buffer, "wt")))
    {
      perror(buffer);
      return 1;
    }
  
  fputs("<?xml version=\"1.0\"?>\n", output);
  fputs("<!DOCTYPE rnaml PUBLIC \"" DTD_DESC "\" \"" DTD_URL "\">\n\n", output);
  fputs("<rnaml version=\"1.1\" comment=\"" COMMENT "\">\n\n", output);
  fputs("  <analysis id=\"UNAFold\">\n", output);
  fputs("    <program>\n", output);
  fputs("      <prog-name>UNAFold</prog-name>\n", output);
  fputs("      <prog-version>" PACKAGE_VERSION "</prog-version>\n", output);
  fputs("    </program>\n", output);
  fputs("  </analysis>\n\n", output);
  
  seqnew_bool = 1;
  seq_num = 1;
  model_num = 0;
  lines = NULL;
  ssLines = NULL;
  
  while (1)
    {
      memcpy(&oldHead, &head, sizeof(struct ctHead));
      memset(&head, 0, sizeof(struct ctHead));
      oldLines = lines;

      if (!readCtFile(ctFile, &head, &lines) || !isValidCt(&head, lines))
	break;

      makeValidId(head.name);

      if ((seq_num != 1 || model_num != 0) && (head.num_bases != oldHead.num_bases || strcmp(head.name, oldHead.name) || diff(head.num_bases, lines, oldLines)))
	{
	  seqnew_bool = 1;
	  ++seq_num;
	  model_num = 1;
	}
      else
	++model_num;

      strcpy(buffer, argv[1]);
      strcat(buffer, ".ss");
      ssFile = fopen(buffer, "rt");
      if (!ssFile)
	{	
	  sprintf(buffer, "%s_%d.ss", argv[1], model_num);
	  ssFile = fopen(buffer, "rt");
	}
      if (ssFile)
	{
	  free(ssLines);
	  ssLines = xmalloc(head.num_bases * sizeof(struct ssLine));
	  if (!readSsFile(ssFile, ssLines, head.num_bases))
	    {
	      free(ssLines);
	      ssLines = NULL;
	    }
	  if (!isMatch(head.num_bases, lines, ssLines))
	    {
	      fprintf(stderr, "Warning: %s is inconsistent and will be ignored\n", buffer);
	      free(ssLines);
	      ssLines = NULL;
	    }
	}
      
      
      /*this is the beginning of the output to the file*/
      if (seqnew_bool)
	{  
	  if (seq_num != 1)
	    {
	      fputs("    </structure>\n", output);
	      fputs("  </molecule>\n", output);
	    }
	  
	  if (seq_num % 10 == 1)
	    fprintf(output, "\n  <!-- %dst sequence and its associated foldings ******************* -->\n", seq_num);
	  else if (seq_num % 10 == 2)
	    fprintf(output, "\n  <!-- %dnd sequence and its associated foldings ******************* -->\n", seq_num);
	  else if (seq_num % 10 == 3)
	    fprintf(output, "\n  <!-- %drd sequence and its associated foldings ******************* -->\n", seq_num);
	  else
	    fprintf(output, "\n  <!-- %dth sequence and its associated foldings ******************* -->\n", seq_num);
	  
	  if (head.is_dna && head.is_rna)
	    fputs("  <!-- Warning: type unknown as sequence contains both Ts and Us -->\n", output);
	  if (head.is_dna && !head.is_rna)
	    fprintf(output, "  <molecule id=\"%s\" type=\"dna\">\n", head.name);
	  else if (head.is_rna && !head.is_dna)
	    fprintf(output, "  <molecule id=\"%s\" type=\"rna\">\n", head.name);
	  else
	    fprintf(output, "  <molecule id=\"%s\">\n", head.name);

	  for (strand = 1; strand <= head.num_strands; ++strand)
	    {
	      for (strandBegin = 0; lines[strandBegin].strand < strand; ++strandBegin);
	      for (strandEnd = strandBegin; strandEnd < head.num_bases && lines[strandEnd].strand == strand; ++strandEnd);

	      if (lines[strandBegin].previous == lines[strandEnd - 1].seqnum && 
		  lines[strandEnd - 1].next == lines[strandBegin].seqnum)
		fprintf(output, "    <sequence circular=\"true\" strand=\"%d\">\n", strand);
	      else
		fprintf(output, "    <sequence circular=\"false\" strand=\"%d\">\n", strand);
	      fputs("      <seq-data>\n        ", output);

	      for (i = 0; strandBegin + i < head.num_bases && i < strandEnd; ++i)
		if (lines[strandBegin + i].strand == strand)
		  {
		    fprintf(output, "%c", lines[strandBegin + i].base);
		    if (i % 60 == 59)
		      fputs("\n        ", output);
		    else if (i % 10 == 9)
		      fputs(" ", output);
		  }
	      fputs("\n      </seq-data>\n", output);
	      fputs("    </sequence>\n\n", output);
	    }
	  fputs("    <structure analysis-ids=\"UNAFold\">\n", output);
	  seqnew_bool = 0;
	}

      fprintf(output, "      <!-- Start of folding %d for sequence %d ************************ -->\n", model_num, seq_num);
      fprintf(output, "      <model id=\"_%d.%d\">\n", seq_num, model_num);
      fputs("        <model-info>\n", output);
      
      if (head.has_dG)
	{
	  if (head.has_initial)
	    fprintf(output, "          <free-energy comment=\"refined\">%.2f</free-energy>\n", head.dG_val);
	  else
	    fprintf(output, "          <free-energy>%.2f</free-energy>\n", head.dG_val);
	}      
      if (head.has_initial)
	fprintf(output, "          <free-energy comment=\"initial\">%.2f</free-energy>\n", head.initial_val);
      
      fputs("        </model-info>\n\n", output);
      
      fputs("        <str-annotation>\n", output);
      for (i = 0; i < head.num_bases; i++)
	if (lines[i].seqnum < lines[i].paired)
	  {
	    helixLength = 1;
	    if (lines[i].paired > 1)
	      while (i + helixLength < head.num_bases && lines[i + helixLength].paired == lines[i].paired - helixLength)
		++helixLength;

	    if (helixLength > 1)
	      fputs("          <helix>\n", output);
	    else
	      fputs("          <base-pair>\n", output);
	    fputs("            <base-id-5p>\n", output);
	    fputs("              <base-id>\n", output);
	    if (head.num_strands > 1)
	      fprintf(output, "                <strand>%d</strand>\n", lines[i].strand);
	    fprintf(output, "                <position>%d</position>\n", i + 1) ;
	    fputs("              </base-id>\n", output);
	    fputs("            </base-id-5p>\n", output);
	    fputs("            <base-id-3p>\n", output);
	    fputs("              <base-id>\n", output);
	    if (head.num_strands > 1)
	      fprintf(output, "                <strand>%d</strand>\n", lines[lines[i].paired - 1].strand);
	    fprintf(output, "                <position>%d</position>\n", lines[i].paired);
	    fputs("              </base-id>\n", output);
	    fputs("            </base-id-3p>\n", output);
	    if (helixLength > 1)
	      {
		fprintf(output, "            <length>%d</length>\n", helixLength);
		fputs("          </helix>\n", output);  
		i += (helixLength - 1);
	      }
	    else
	      fputs("          </base-pair>\n", output);
	  }
      for (i = 0; i < head.num_bases; i++)
	{
	  if (lines[i].seqnum < lines[i].dnst)
	    /* Don't bother outputting stacking where it's implied by a <helix>. */
	    if (lines[i].dnst != i + 2 || lines[i].paired - 1 != lines[i + 1].paired)
	      {
		fputs("          <base-stack>\n", output);
		fputs("            <base-id>\n", output);
		if (head.num_strands > 1)
		  fprintf(output, "              <strand>%d</strand>\n", lines[i].strand);
		fprintf(output, "              <position>%d</position>\n", i + 1) ;
		fputs("            </base-id>\n", output);
		fputs("            <base-id>\n", output);
		if (head.num_strands > 1)
		  fprintf(output, "              <strand>%d</strand>\n", lines[lines[i].dnst - 1].strand);
		fprintf(output, "              <position>%d</position>\n", lines[i].dnst) ;
		fputs("            </base-id>\n", output);
		fputs("          </base-stack>\n", output);
	      }
	  if (lines[i].seqnum < lines[i].upst)
	    {
	      fputs("          <base-stack>\n", output);
	      fputs("            <base-id>\n", output);
	      if (head.num_strands > 1)
		fprintf(output, "              <strand>%d</strand>\n", lines[i].strand);
	      fprintf(output, "              <position>%d</position>\n", i + 1) ;
	      fputs("            </base-id>\n", output);
	      fputs("            <base-id>\n", output);
	      if (head.num_strands > 1)
		fprintf(output, "              <strand>%d</strand>\n", lines[lines[i].upst - 1].strand);
	      fprintf(output, "              <position>%d</position>\n", lines[i].upst) ;
	      fputs("            </base-id>\n", output);
	      fputs("          </base-stack>\n", output);
	    }
	}
      fputs("        </str-annotation>\n", output);
      if (ssLines)
	{
	  fputs("\n        <secondary-structure-display>\n", output);
	  for (i = 0; i < head.num_bases; ++i)
	    {
	      fputs("          <ss-base-coord>\n", output);
	      fputs("            <base-id>\n", output);
	      if (head.num_strands > 1)
		fprintf(output, "              <strand>%d</strand>\n", lines[i].strand);
	      fprintf(output, "              <position>%d</position>\n", i + 1);
	      fputs("            </base-id>\n", output);
	      fprintf(output, "            <coordinates>%g %g</coordinates>\n", ssLines[i].x, ssLines[i].y);
	      fputs("          </ss-base-coord>\n", output);
	    }
	  fputs("        </secondary-structure-display>\n", output);
	}
      fputs("      </model>\n", output);
      fprintf(output, "      <!-- End of folding %d for sequence %d ************************ -->\n\n", model_num, seq_num);
      free(oldLines);
    }	
  
  if (model_num)
    {
      fputs("    </structure>\n", output);
      fputs("  </molecule>\n\n", output);
    }
  fputs("</rnaml>\n", output);
  
  return 0;
}

int readCtFile(FILE* file, struct ctHead* head, struct ctLine** lines)
{
  int i;
  int currentStrand = 1;
  char current_line[300];

  if (!fgets(current_line, 300, file))
    return 0;

  if (sscanf(current_line, "%d dG = %g [initially %g] %[^\n]",
	     &head->num_bases, &head->dG_val, &head->initial_val, head->name) == 4)
    {
      head->has_dG = 1;
      head->has_initial = 1;
    }
  else if (sscanf(current_line, "%d dG = %g dH = %*g %[&\n]", &head->num_bases, &head->dG_val, head->name) == 3)
    head->has_dG = 1;
  else if (sscanf(current_line, "%d dG = %g %[^\n]", &head->num_bases, &head->dG_val, head->name) == 3)
    head->has_dG = 1;
  else if (sscanf(current_line, "%d ENERGY = %g [initially %g] %[^\n]",
		  &head->num_bases, &head->dG_val, &head->initial_val, head->name) == 4)
    {
      head->has_dG = 2;
      head->has_initial = 1;
    }
  else if (sscanf(current_line, "%d ENERGY = %g %[^\n]", &head->num_bases, &head->dG_val, head->name) == 3)
    head->has_dG = 2;
  else if (sscanf(current_line, "%d %[^\n]", &head->num_bases, head->name) == 2)
    ;
  else
    return 0;

  if (head->num_bases <= 0)
    return 0;

  while (head->name[strlen(head->name) - 1] == ' ')
    head->name[strlen(head->name) - 1] = '\0';

  *lines = xmalloc(head->num_bases * sizeof(struct ctLine));

  for (i = 0; i < head->num_bases; i++)
    {
      if (!fgets(current_line, 100, file))
	return 0;
      // The last two columns are optional, but we should read at least six fields or the file is corrupt.
      if (sscanf(current_line, "%d %c %d %d %d %d %d %d", &(*lines)[i].seqnum, &(*lines)[i].base, &(*lines)[i].previous, &(*lines)[i].next, &(*lines)[i].paired, &(*lines)[i].historic, &(*lines)[i].upst, &(*lines)[i].dnst) < 6)
	return 0;

      if ((*lines)[i].seqnum != i + 1)
	{
	  puts("The sequence does not go up in sequentially increasing numeric order");
	  return 0;
	}
      
      if (i > 0 && (*lines)[i].previous != i && (*lines)[i - 1].next != i + 1)
	++currentStrand;

      (*lines)[i].strand = currentStrand;

      if (toupper((*lines)[i].base) == 'T')
	head->is_dna = 1;
      else if (toupper((*lines)[i].base) == 'U')
	head->is_rna = 1;
    }

  head->num_strands = currentStrand;

  return 1;
}

int isValidCt(struct ctHead* head, struct ctLine* lines)
{
  int i;

  for (i = 0; i < head->num_bases; ++i) {
    if (lines[i].paired && lines[lines[i].paired - 1].paired != i + 1)
      {
	printf("Pairing information is inconsistent between bases %d and %d", i + 1, lines[i].paired);
	return 0;
      }
    if (lines[i].upst && lines[lines[i].upst - 1].dnst != i + 1)
      {
	printf("Stacking information is inconsistent between bases %d and %d", i + 1, lines[i].upst);
	return 0;
      }
    if (lines[i].dnst && lines[lines[i].dnst - 1].upst != i + 1)
      {
	printf("Stacking information is inconsistent between bases %d and %d", i + 1, lines[i].dnst);
	return 0;
      }
  }

  return 1;
}

int readSsFile(FILE* file, struct ssLine* lines, int n)
{
  int i;
  char current_line[300];

  for (i = 0; i < n; ++i)
    {
      if (!fgets(current_line, 300, file))
	return 0;
      if (sscanf(current_line, "%d %c %lg %lg %*d %d", &lines[i].seqnum, &lines[i].base, &lines[i].x, &lines[i].y, &lines[i].paired) != 5)
	return 0;
      if (lines[i].seqnum != i + 1)
	return 0;
    }

  return 1;
}

int isMatch(int n, struct ctLine* ctLines, struct ssLine* ssLines)
{
  int i;

  for (i = 0; i < n; ++i)
    if (ctLines[i].base != ssLines[i].base || ctLines[i].paired != ssLines[i].paired)
      return 0;

  return 1;
}

int diff(int size, struct ctLine* new, struct ctLine* old)
{
  int i;

  for (i = 0; i < size; ++i)
    if (new[i].base != old[i].base)
      return 1;

  return 0;
}

void makeValidId(char* str)
{
  char buffer[250];

  if (('0' <= *str && *str <= '9') || *str == '-' || *str == '.')
    {
      strcpy(buffer, str);
      strcpy(str, "_");
      strcat(str, buffer);
    }

  for (; *str; ++str)
    if (!isalnum(*str) && *str != '_' && *str != '-' && *str != '.')
      *str = '_';
}
