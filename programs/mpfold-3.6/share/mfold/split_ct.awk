# Split a ct file into individual files.

BEGIN	{
	count = 0; i = 0; n = 0;
	if ( ARGC-- != 3 ) {
	  print "usage:";
	  print "nawk -f split.awk file-to-split prefix-for-split-files";
	  exit 1;
	  }
	CurrentFile = ARGV[2];
 	}
{
	i++;
	if (i==(n+1) || n==0) {
	n = $1; i = 0;
	if ( CurrentFile != ARGV[2] ) close(CurrentFile);
	count++;
	CurrentFile = ARGV[2] "_" count ".ct";
	print "CurrentFile=",CurrentFile;
	};
	print $0 > CurrentFile;
}

END {	if ( count ) close(CurrentFile); }


