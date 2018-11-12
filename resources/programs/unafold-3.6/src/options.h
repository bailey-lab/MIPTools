int option_code;

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"NA", required_argument, NULL, 'n'},
  {"tmin", required_argument, NULL, 't'},
  {"tinc", required_argument, NULL, 'i'},
  {"tmax", required_argument, NULL, 'T'},
  {"suffix", required_argument, NULL, 's'},
  {"allpairs", no_argument, &option_code, 11},
  {"maxloop", required_argument, &option_code, 1},
  {"output", required_argument, NULL, 'o'},
  {"debug", no_argument, NULL, 'd'},
  {"nodangle", no_argument, &option_code, 2},
  {"simple", no_argument, &option_code, 3},
  {"sodium", required_argument, NULL, 'N'},
  {"magnesium", required_argument, NULL, 'M'},
  {"polymer", no_argument, NULL, 'p'},
  {"prefilter", required_argument, &option_code, 8},
  {"nopostfilter", no_argument, &option_code, 9},
  {"prohibit", required_argument, NULL, 'r'},
#if ENABLE_FORCE
  {"force", required_argument, NULL, 'f'},
#endif
  {"energyOnly", no_argument, NULL, 'E'},
  {"noisolate", no_argument, NULL, 'I'},
  {"mfold", optional_argument, NULL, 'F'},
  {"quiet", no_argument, NULL, 'q'},
  {"tracebacks", required_argument, NULL, 'k'},
  {"scale", required_argument, &option_code, 10},
  {"zip", no_argument, NULL, 'z'},
  {"maxbp", required_argument, NULL, 'm'},
  {"constraints", optional_argument, NULL, 'c'},
  {"basepairs", required_argument, NULL, 'b'},
  {"circular", no_argument, &option_code, 13},
  {"stream", no_argument, &option_code, 12},
  {NULL, 0, NULL, 0}
};

#define OPTION_TAKES_TWO   1
#define OPTION_DEBUG       2
#define OPTION_SIMPLE      4
#define OPTION_NOISOLATE   8
#define OPTION_MIN        16
#define OPTION_QUIET      32
#define OPTION_TRACEBACK  64
#define OPTION_ZIP       128
#define OPTION_NODANGLE  256
#define OPTION_MAXBP     512
#define OPTION_CIRCULAR 1024

void usage(char* program, int options)
{
  if (options & OPTION_TAKES_TWO)
    printf("Usage: %s [options] file1 file2\n", program);
  else
    printf("Usage: %s [options] file\n", program);
  puts("");
  puts("Options:");
  puts("-V, --version");
  puts("-h, --help");
  puts("-n, --NA=(RNA | DNA) (defaults to RNA)");
  if (options & OPTION_MIN)
    puts("-t, --tmin=<minimum temperature> (defaults to 37)");
  else
    puts("-t, --tmin=<minimum temperature> (defaults to 0)");
  puts("-i, --tinc=<temperature increment> (defaults to 1)");
  if (options & OPTION_MIN)
    puts("-T, --tmax=<maximum temperature> (defaults to 37)");
  else
    puts("-T, --tmax=<maximum temperature> (defaults to 100)");
  puts("-s, --suffix=<free energy suffix>");
  puts("-o, --output=<prefix>");
  if (options & OPTION_DEBUG)
    puts("-d, --debug");
  puts("-N, --sodium=<[Na+] in M> (defaults to 1)");
  puts("-M, --magnesium=<[Mg++] in M> (defaults to 0)");
  puts("-p, --polymer");
  puts("-r, --prohibit=<i,j,k>");
#if ENABLE_FORCE
  puts("-f, --force=<i,j,k>");
#endif
  puts("-E, --energyOnly");
  if (options & OPTION_NOISOLATE)
    puts("-I, --noisolate");
  if (options & OPTION_ZIP)
    puts("-z, --zip");
  if (options & OPTION_MIN)
    puts("-F, --mfold[=<P,W,MAX>] (defaults to 5,*,100; W determined by sequence length)");
  if (options & OPTION_QUIET)
    puts("-q, --quiet");
  if (options & OPTION_TRACEBACK)
    puts("-k, --tracebacks=<number of tracebacks>");
  if (options & OPTION_MAXBP)
    puts("-m, --maxbp=<maximum basepair distance>");
  puts("-c, --constraints[=<name of constraint file>] (defaults to prefix.aux)");
  puts("-b, --basepairs=<name of basepairs file>");
  if (options & OPTION_CIRCULAR)
    puts("    --circular");
  puts("");
  puts("Obscure options:");
  puts("    --allpairs");
  puts("    --maxloop=<maximum bulge/interior loop size> (defaults to 30)");
  if (options & OPTION_NODANGLE)
    puts("    --nodangle");
  if (options & OPTION_SIMPLE)
    puts("    --simple");
  puts("    --prefilter=<value1, value2>");
  if (!(options & OPTION_MIN))
    puts("    --nopostfilter");
  if (options & OPTION_QUIET)
    puts("    --stream");
  puts("");
  puts("Report bugs to " PACKAGE_BUGREPORT);
  exit(EXIT_SUCCESS);
}
