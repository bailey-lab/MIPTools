#!/usr/bin/perl

use strict 'vars';
#use bytes;


#USEFUL ENVIRONMENTAL VARIABLES#
#print "$^O\n";
#print "$] \n";
#print "$0 \n";
#print "H$ENV{HOME} S\n";
#print "P$ENV{PATH} S$ENV{SHELL}\n";

#LOAD MODULES #########################################################
use Getopt::Long;
use Tk;
use Cwd ;
use Tk::Dialog;
use Tk::Balloon;
use Tk::BrowseEntry;

####THINGS TO ADD/DO #############################################################

#1) A jump and mark function to search multiscreen views
#2) Change to drawing boxes instead of really thick lines for subjects and extras
#    This should allow some cool stuff using Illustator effects such as lighting.

############### PARASIGHT EVOLUTION  (some of the puncuated events) ##########################
#030226 decomtaminate algorithm added to remove high copies over a certain range and all related pairwise.
#030122 added better help descriptions for options and more pod help documentation
#021207 fixed printing finally (how obtuse)
#021203 fixed SscaleC line issues (a little more intutitive now)
#021029 fixed precode issues
#021029 returns are usable now in text_text options
#020821 remove Storable module usage (versions incompatibities were the last straw)
#020411 add FILE option for precode so that it can work with Win32 limited command prompt
#020410 start version 7.2-- clean up for publication
#011026 added support for extracting sequences from files and by fastacmd
#011015 fixed subseq drawing error
#011014 added description array for every variable
#011014 added four execute commands that can  be combined with columns
#011014 modified more variable names changing sub-> seq
#011002 added graph1 and graph2 option
#011002 fixed extra Left Click popup window
#011002 removed verbose output 
#011002 removed defn from extra option--just seq begin end
#011002 added default template locations for unix
#011002 change sub_labelhit_col to sub_labelhit_col
#010822 added menu view of alignments
#010814 fix the option template--it didn't appear to be loading over the defaults
#010807 create my own pop up window--ballon help is slow  and adding overhead
#010529 added extra level order changing and quick color
#010525 added quick color and level order changing for send to background and foreground
#010523 added template option to save menu
#010523 added large amounts of documenation and more popup help windows
#010523 change random subject color so that subject sequence will have consistant color
#010523 added searching through conditional coloring as requested by jules (better command processing)
#010522 added arrows to objects (technological breakthough!)
#010519 added color variation based on conditional pseudo-conditional statements calle hitcondition
#011002 changed save to be *.psa (alignment array) *.pse (extra table) and *.pso (options)
#010517 fixed sequence drawing so that numerical only sequence names can be used#
#011002 changed highlight color to yellow
#011002 added ability to specify arrangement using a file with seq and begin position
#010506 revamped option display to index cards
#010506 fixed seq1 and seq2 so that subject labeling is proper order
#010505 minor cosmetic changes and error fixes
#010505 add new options: -option, -die, and -precode
#010116 add program stats section for program_manager.pl
#001221 fixed -in when to stop alignments reloading on subsequent executions#
#001211 fixed a whole lot of little bugs#
#001011 added vertical scales
#001015 cleaned up interface and added new options

#######################################
#######################################
#####hard-wired defaults###############
######################################


#DECLARE GLOBAL VARIABLES#######################################
use vars qw($n1a $n2a $n1 $n2 $loaded $tmp $b1f $e1f $b2f $e2f);
use vars qw($x_max $x_min $opt_t $opt_m $margin $l1 $l2 %subscaleC @subscaleC);

use vars qw($column_header_display);
use vars qw($first_pass $widest_line);
use vars qw($mw $balloon $ballooni $canvas $frame $fontsize $scale_x $scrolledcanvas $output $scale);
use vars qw($file %deleted_pairwise);
use vars qw(%acc %accsub @acc_order @acc_ordersub @l %lpos %msghash);
use vars qw(@m %mh  @mheader $mstring);
use vars qw(@e %eh @eheader $estring);
use vars qw(@g1 @g2);
use vars qw(%pairwise2delete);
use vars qw(%opt %colheader %newopt $options %iinfo %optdesc);
use vars qw($filepath $optionpath );
use vars qw($canvas_width $bp_per_pixel);

###########################################################
############# PROGRAM DESCRIPTION #########################
###########################################################

use vars qw($program $pversion $pdescription $pgenerate $pusage);
$program = "$0";
$program =~ s/^.*\///;
### program stats ###
$pversion='7.6';
$pdescription = "$program (ver:$pversion) displays pairwise alignments and accompaning annontation in a wide variety of formats";
$pgenerate= 'jeff: labmates: genetics: public:';
$pusage="$program [-in file/-align file/-showseq file] [other data] [other options]";
### program stats end ###

###########################################################
############PRE-DEFAULTS###################################

#template paths#
use vars qw($template_path $export_text_path);
$export_text_path='untitled.txt';
$template_path='~/.PARASIGHT:/people/PARASIGHT';
$template_path=~ s/~/$ENV{HOME}/g;

#default show#
use vars qw($default_show);
$default_show='ALL';



##########################################################
#########################HELP#############################
##########################################################

if (! defined $ARGV[0]) {
print "USAGE 
$pusage
DESCRIPTION\n$pdescription

MAIN DATA INPUT COMMANDS
-in [filepath]  load a saved parasight view (*.pso *.psa *.pse ...)
-align [filepath1:filepath2:etc] load pairwise alignment table(s)
-showseq [ALL | file | seqname(s):]  names of sequences to display
   no colon = load as file of names
   colon(:) = parse as list of colon-delimited seq names
ADDITIONAL DATA INPUT COMMANDS
-extra [filepath1:filepath2:etc] loads extra simple sequence feature(s)
   (e.g. exons, introns, and repeats)
-graph1 [filepath1:filepath2:etc]
   graphs a set of values on a scale above the sequence at given positions 
   (e.g. moving windows such as %GC)
-graph2 [filepath1:filepath2:etc]
   adds another line to the graph scale	
OPTIONAL COMMANDS
-arrangeseq  [oneperline | sameline | file:filepath] arranges sequences
   file:filepath allows specific positions to be assigned to sequences
-arrangesub [oneperline|stagger|subscaleS|subscaleV] (default stagger)
   easier to manipulate from user interface then at command line
-colorsub [NONE | RESET | seqrandom | hitrandom | hitcondition]
-options ['opt1=>value1,opt2=>value2']  change any parasight option
   (for on/off, yes/no or true/false options use 1 and 0 as input)
   (e.g. 'canvas_width=>500,seq_tick_on=>1,-arrangeseq=>oneperline')
-showsub [ALL|file|seqnames:] names of subjects to display
   no colon = load as file of names
   colon(:) = parse as list of colon-delimited seq names
-template [filepath] load a template file containing options to apply
-showseqqueryonly [switch] will only draw sequences in first (blast query) position
-quiet [switch] decreases screen output
ADVANCED COMMANDS
-minload [switch] loads only the pairwise relavant to current -showseq
   (quicker when just certain sequences are needed from large files)
-precode 'perl code commands to execute after first screen draw'
-die  parasight exits after executing precode
FULL DOCUMENTATION
-help
HINT FOR BEGINNERS - GET YOUR DATA LOADED AND THEN MANIPULATE IT THROUGH THE GUI OPTIONS!
";
exit;  
}


######################################################################################################
############################OPTION HANDLING ##########################################################
######################################################################################################
if ( &GetOptions(\%newopt, "in=s", "align=s", "showseq=s","showsub=s","color=s",'colorsub=s', "extra=s","arrangeseq=s","arrangesub=s","options=s", 
				"template=s","graph1=s","graph2=s","showqueryonly","quiet","minload", "precode=s","die","help") ) {
	print "Command line arguments parsed sucessfully\n" if !$newopt{'quiet'};
} else {
	die "\nCommand line arguments were not sucessfully parsed\n";
}
###EXECUTE POD DOCMENTATION/LONG HELP ####
if ($newopt{'help'}) {
	system "perldoc $0\n";
	exit;
}


if (!$newopt{'in'} && !$newopt{'align'} && !$newopt{'extra'} && !$newopt{'graph1'}) {
	print "*******************************************************************************\n" if !$newopt{'quiet'};
	print "**WARNING: No major options (-in, -align,-graph or -extra) were provided! *****\n" if !$newopt{'quiet'};
	print "*******************************************************************************\n" if !$newopt{'quiet'};

}

###########################################################
############ load options into from a template file #######
############ loaded directly into %opt ####################
###########################################################

$options.='';
if ( $newopt{'template'} ) {
	###require file to end in .pst###
	$newopt{'template'}=~s/\.pst//;
	$newopt{'template'}.='.pst';
	my @dir=split /[:;] */,$template_path;
	if (-e $newopt{'template'} ) {
		&load_option_template($newopt{'template'});
		$newopt{'template'}='';
	} else {
		###check for and if found load from default directories#
		foreach(@dir) {
		 	my $p="$_/$newopt{'template'}";
		 	if ( -e $p) {
		 		#load and clear template#
		 		&load_option_template($p);
				$newopt{'template'}='';
				last;
			}
		}
	}
	##template should be blank if found#
	if ($newopt{'template'} ne '') {
		print "ERROR:Can't locate requested option template ($newopt{'template'})\n";
		foreach (@dir) {
			opendir (DIR, "$_" ) || die "Can't read template directory ($_)!\n";
			my @templates=grep { /\.pst$/} readdir DIR;
			print "***VALID templates to choose from in ($_) are: [", join (", ", @templates), "]\n";
		}
		die "\n";
	}
}

#############################################
######## parse command line options #########
#############################################

if ( $newopt{'options'} ) {
	if ($newopt{'options'} !~ /=>/ ) {
		####add file loading with format option=>value###
		open (OPTIONS, $newopt{'options'}) || die "Can't read ($newopt{'options'}) ($!)!\n";
		###this had yet to be implemented###
		$newopt{"options"}='';
		while (<OPTIONS>) {
			next if /^#/;
			chomp;
			s/,$//;
			$newopt{'options'}.= "$_ ,";
		}
			
	}
	$options.=$newopt{'options'};
	delete $newopt{'options'};
}

##########################################
###########parse other major options #####

$newopt{'showseq'}=$default_show if !defined $newopt{'showseq'} &&  !defined $newopt{'in'};
$newopt{'showsub'}='ALL' if !defined $newopt{'showsub'} && !defined $newopt{'in'};
if ($newopt{'arrangeseq'}=~/file:(.+)/ ) {
	$newopt{'arrange_file'}=$1;
	if (open (TEST, $newopt{'arrange_file'} ) ) {
		 close TEST;
		$newopt{'arrangeseq'}='file';
	} else {
		warn "-arrangeseq $newopt{'arrangeseq'} could not be opened\n";
		$newopt{'arrangeseq'}='';
	}
}
$newopt{'arrangeseq'}='oneperline' if !defined $newopt{'showsub'} && !defined $newopt{'in'} &&!defined $newopt{'template'};
$newopt{'arrangesub'}='stagger' if !defined $newopt{'arrangesub'} && !defined $newopt{'in'}&&!defined $newopt{'template'};
$newopt{'color'}='None' if !defined $newopt{'color'} && !defined $newopt{'in'}&&!defined $newopt{'template'};
$newopt{'colorsub'}='None' if !defined $newopt{'colorsub'} && !defined $newopt{'in'}&&!defined $newopt{'template'};

################################################################
#######PARSE NEW OPTIONS OUT INTO PROPER FORMAT#################
if ( $newopt{'arrangesub'} ) { 
	if ( $newopt{'arrangesub'} =~/(subscale[NC]+):([A-Za-z_0-9#]+):([\-0-9.:]+)$/  ) {
		$newopt{'arrangesub'}=$1;
		($newopt{'sub_scale_col'},$newopt{'sub_scale_col2'})=split '#',$2;
		#print "$1|$2|$3|$4|$5\n";
		($newopt{'sub_scale_lines'},$newopt{'sub_scale_min'},$newopt{'sub_scale_max'},$newopt{'sub_scale_step'})= split ':',$3;
	}
	#print "$newopt{'sub_scale_min'},$newopt{'sub_scale_max'},$newopt{'sub_scale_lines'}\n"
}


#############################################################
####################LOAD OLD OPTIONS OR GET DEFAULTS ########
if ($newopt{'in'}) {
	$newopt{'in'}=~s/\.ps[aeo]?$//;
	$newopt{'in'}.='.psa';
	&load_parasight_table( $newopt{'in'} );
	$opt{'align'}=   $newopt{'align'};
	$opt{'extra'} =  $newopt{'extra'};
	$opt{'graph1'} = $newopt{'graph1'};
	$opt{'graph2'} = $newopt{'graph2'};
} 

################################################################
#####general option built-in defaults and descriptions #########
################################################################
# any of these options will be overridden with -option command #
################################################################
#### 1 begin is just a place holder to check for to make sure array is synced
#### 2 is name of the option  for %opt array
#### 3 is the default value of the option 
#### 4 is the description of the option for %optdesc array

my @todefine=(
#'just_pairwise_regions' , 0 ,  'THIS OPTION IS INACTIVATED AND NONFUNCTIONAL'  , 


'begin',	'alignment_col' , 0,  '[integer] column for the first query (first sequence) in a parsed pairwise alignment. Blank/zero hides option from popup menu. The sequence will contain dashes for gaps.'  , 
'begin',	'alignment_col2' , 0,  '[integer] column for subject (second pairwise position) sequence alignment. Blank/zero hides option from popup menu. The sequence will contain dashes for gaps.'  , 
'begin',	'alignment_wrap' , 50,  '[integer] line width in aligned characters (bases/amino acids/dashes) for displaying any alignments'  , 
'begin',	'arrangeseq' ,  'oneperline' ,  '[oneperline|sameline|file] determines the arrangement of sequences that are currently being shown with -showseq. Choices: oneperline = each sequence placed on a separate line; sameline = sequences are place one after the other on the same line; file = load a file with exact positions in terms of line number and base position within the colorsub_hitcond_tests variable'  , 
'begin',	'arrangesub' , 'stagger' ,  '[stagger|oneperline|subscaleC|subscaleN] basic'  , 
'begin',	'arrangesub_stagger_spacing' , 40000 ,  '[integer] bases of spacing between for sequences placed on the same sub line.  Sequences separated by less than this distance from each other will be placed on separate sub lines). This option is useful for providing space for a label.'  , 
'begin',	'canvas_bpwidth' , 250000 ,  '[integer] number of bases that the width of screen represents (not including indentations).  This is the  number of bases per line across'  , 
'begin',	'canvas_indent_left' , 60 ,  '[integer] pixels to indent from the left-side of screen window image (the drawing areas is the canvas in Tk)'  , 
'begin',	'canvas_indent_right' , 30 ,  '[integer] pixels to indent from the right-side of screen window image  (the drawing areas is the canvas in Tk)'  , 
'begin',	'canvas_indent_top' , 40 ,  '[integer] pixels to indent from top of screen before drawing sequence lines (it does not take into account graphs or extras)  (the drawing areas is the canvas in Tk)'  , 
'begin',	'color' , 'None' ,  '(not implemented yet)'  , 
'begin',	'colorsub' , 'None' ,  '[NONE|RESET|hitrandom|seqrandom|hitconditional] Choices for coloring subs: NONE=no coloring routines; RESET=clear all assigned colors to pairwise; hitrandom=randomly color each hit/pairwise a different color; seqrandom=randomly color each defined seequence; hitconditional=color each hit based on pseudo-perl if than statements found in the varaible'  , 
'begin',	'colorsub_hitcond_col' , 34 ,  '[integer] column against which to test conditional statements in pairwise data (does not work on extra items or graphs)'  , 
'begin',	'colorsub_hitcond_tests' , 'red if <2; orange if <0.99; yellow if <0.98; green if <0.97; blue if <0.96; purple if <0.95; brown if <0.94; grey if <0.93; black if <0.92; pink if <0.91' ,  '[fake code] conditional statements to color pairwise hits based on the values in the column colorsub_hitcond_col (format for tests: color [= or < or >] value; )'  , 
'begin',	'execute' , '',  '[external system command] to execute on Control-Shift-Click Left Button'  , 
'begin',	'execute2' , '',  '[external system command] to execute on Control-Shift-Click Middle Button'  , 
'begin',	'execute2_array' , 'm',  '[e|m] extra or pairwise array to use in execute2 command'  , 
'begin',	'execute2_desc', '','[text] description to display in right-click menu for execute2 command',
'begin',	'execute3' , '',  '[external system command] to execute on Control-Shift-Click Right Button'  , 
'begin',	'execute3_array' , 'm',  '[e|m] extra or pairwise array to use in execute3 command'  , 
'begin',	'execute3_desc', 'widget','[text] description to display in right-click menu for execute3 command',
'begin',	'execute4' , '',  '[external system command] to execute from within right-click menu only'  , 
'begin',	'execute4_array' , 'm',  '[e|m] extra or pairwise array to use in execute4 command'  , 
'begin',	'execute4_desc', '','[text] description to display in right-click menu for execute command',
'begin',	'execute_array' , 'e',  '[e|m] extra or pairwise array to use in execute command'  , 
'begin',	'execute_desc', '','[text] description to display in right-click menu for execute command',
'begin',	'extra_arrow_diag' , 5 ,  '[integer] distance from point of arrow to wing/elbow of arrow'  , 
'begin',	'extra_arrow_on' , 1 ,  '[0|1] toggles arrows for extras off and on'  , 
'begin',	'extra_arrow_para' , 5 ,  '[integer] pixel distance from point of arrow along the line'  , 
'begin',	'extra_arrow_perp' , 4 ,  '[integer] pixel distance from base on line to wing of arrow'  , 
'begin',	'extra_color' , 'purple' ,  '[color] default for extra object'  , 
'begin',	'extra_label_col' , 10 ,  '[integer] column to take values to use for extra labels'  , 
'begin',	'extra_label_col_pattern' , '' ,  '[regular expression] pattern to match (and extract via parentheses) replacing current value.  Allows display of only part of the data found in a column.'  , 
'begin',	'extra_label_color' , 'purple' ,  '[color] default for the labels of extra objects'  , 
'begin',	'extra_label_fontsize' , 6 ,  '[integer] font size (in points) the labels of extra objects'  , 
'begin',	'extra_label_offset' , 2 ,  '[integer] horizontal offset for extra labels (left is negative, right is positive)', 
'begin',	'extra_label_on' , 1 ,  '[0|1] toggles the text label for extra objects off and on'  , 
'begin',	'extra_label_test_col' , '' ,  '[integer] column to test for a pattern--if pattern matched then extra not drawn'  , 
'begin',	'extra_label_test_pattern' , '' ,  '[regular expression] pattern to match in order to NOT draw the matching extra object'  , 
'begin',	'extra_offset' , -4 ,  '[integer] default vertical offset of extra object (negative = up; positive = down)'  , 
'begin',	'extra_on' , 1 ,  '[0|1] toggles off and on the display of all extras'  , 
'begin',	'extra_width' , 6 ,  '[integer] default width (horizontal thickness of extra object'  , 
'begin',	'fasta_blastdb' , 'htg:nt',  '[database names] for sequence fastacmd lookups ',
'begin',	'fasta_directory' , '.:fastax',  '[directories] to search for fasta files corresponding to sequence names in order to extract subsequences on command (names of files must be same as names of sequences)'  , 
'begin',	'fasta_fragsize' , 400000,  '[integer] fragment size for sequences in fasta directory.  Useful for quick lookups in long sequences like chromosomes.  If this is non-zero than fragments of files are searched for in the fasta_directory (nomenclature of fragmented files end with _###, e.g. chr1_000, chr1_001, etc.)'  , 
'begin',	'fasta_on',	1, 	'[0|1]  off|on turns fasta extraction on and off', 
'begin',	'fasta_wrap',	50, 	'[integer] line width in characters for fasta files created', 
'begin',	'filename_color' ,  'grey'  ,  '[color] of text label for the filename'  , 
'begin',	'filename_offset' ,  -10 ,  '[integer] vertical offset of text label for filename (up is negative, down is positive)'  , 
'begin',	'filename_offset_h' ,  0  ,  '[integer] horizontal offset of text label for filename (left is negative, right is positive)'  , 
'begin',	'filename_on' ,  1  ,  '[0|1] toggle off and on display of designated filename/parasight name (initially defined by -in if empty)'  , 
'begin',	'filename_pattern' ,  ''  ,  '[regular expression]  pattern to match in the filename.  Useful for removing the path.  (Although if using graphical interface, it is easier to change the filename.) '  ,
'begin',	'filename_size' ,  10  ,  '[integer] point size of text label shown for the filename'  , 
'begin',	'filter1_col' , '', '[integer]  column that contains data with which to filter pairwise' ,  
'begin',	'filter1_max' , '' ,  '[float] limit for value in filter1_col above which pairwise are NOT drawn'  , 
'begin',	'filter1_min' , '' ,  '[float] limit for value in filter1_col below which pairwise are NOT drawn'  , 
'begin',	'filter2_col' , '', '[integer] column that contains data with which to filter pairwise' ,  
'begin',	'filter2_max' , '' ,  '[float] limit for value in filter2_col above which pairwise are NOT drawn'  , 
'begin',	'filter2_min' , '' ,  '[float] limit for value in filter2_col  below which pairwise are NOT drawn'  , 
'begin',	'filterextra1_col' , '', '[integer]  column number that contains data with which to filter extra objects--sequences are removed before arrange functions are executed' ,  
'begin',	'filterextra1_max' , '' ,  '[float] limit for value in filterextra1_col above which extras are NOT drawn'  , 
'begin',	'filterextra1_min' , '' ,  '[float] limit for value in filterextra1_col below which extras are NOT drawn'  , 
'begin',	'filterextra2_col' , '', '[integer] column number that contains data with which to filter extras' ,  
'begin',	'filterextra2_max' , '' ,  '[float] limit for value in filterextra2_col above which extras are NOT drawn'  , 
'begin',	'filterextra2_min' , '' ,  '[float] limit for value in filterextra2_col below which extras are NOT drawn'  , 
'begin',	'filterpre1_col' , '' ,  '[integer] column that contains data with which to prefilter pairwise--prefiltering removes pairwise before any arranging  (normal filtering removes pairwise after filtering)'  , 
'begin',	'filterpre1_max' , '' ,  '[float] limit for value in filterpre1_col above which pairwise are NOT drawn or arranged'  , 
'begin',	'filterpre1_min' , '' ,  '[float] limit for value in filterpre1_col below which pairwise are NOT drawn or arranged'  , 
'begin',	'filterpre2_col' , '' ,  '[integer] column that contains data with which to prefilter pairwise—prefilter removes pairwise before any arranging'  , 
'begin',	'filterpre2_max' , '' ,  '[float] limit for value in filterpre2_col above which pairwise are NOT drawn or arranged'  , 
'begin',	'filterpre2_min' , '' ,  '[float] limit for value in filterpre2_col below which pairwise are NOT drawn or arranged'  , 
'begin',	'gif_anchor' , 'center' ,  '[center|nw|ne|sw|se|e|w|n] positioning of background gif relative to draw point gif_x and gif_y'  , 
'begin',	'gif_on' , 0 ,  'displays a gif image in background (the image will not print out in postscript)'  , 
'begin',	'gif_path' , '' ,  '[file path] of gif image to display in background--image does not make it into the Postscript file file'  , 
'begin',	'gif_x' , int($opt{'window_width'}/2) ,  '[integer] background picture pixel x coordinate position (top of image is zero)'  , 
'begin',	'gif_y' , 0 ,  '[integer] background gif y coordinate position (0 is top of screen)'  , 
'begin',	'graph1_label_color' , 'blue' ,  '[color] for graph1 labels (left side axis)'  , 
'begin',	'graph1_label_decimal' , 2 ,  '[integer] number of decimal points to round graph1 labels (left side axis)'  , 
'begin',	'graph1_label_fontsize' , 10 ,  '[integer] point size of graph1 labels (left side axis)'  , 
'begin',	'graph1_label_multiplier' , 1 ,  '[float] multiplier for graph1 labels (left side axis)'  , 
'begin',	'graph1_label_offset' , 1 ,  '[integer] horizontal offset for graph1 labels (left side axis)'  , 
'begin',	'graph1_label_on' , 1 ,  '[0|1] toggles on labels for graph1 scale (left side axis)'  , 
'begin',	'graph1_line_color' , 'blue' ,  '[color] for graph1 connecting lines'  , 
'begin',	'graph1_line_on' , 1 ,  '[0|1] toggles graph1 connecting line off and on'  , 
'begin',	'graph1_line_smooth' , 0 ,  '[0|1] toggles on and off smoothing function for connecting line',
'begin',	'graph1_line_width' , 1 ,  '[integer] width for graph1 connecting line'  , 
'begin',	'graph1_max' , 100 ,  '[integer] maximum value of graph1 scale'  , 
'begin',	'graph1_min' , -5 ,  '[integer] minimum value of graph1 scale'  , 
'begin',	'graph1_on' , 0 ,  '[0|1] toggles off and on graph1'  , 
'begin',	'graph1_point_fill_color' , 'blue' ,  '[color] to fill points with for graph1'  ,
'begin',	'graph1_point_on' , 1 ,  '[0|1] toggles point drawing on and off for graph1'  ,
'begin',	'graph1_point_outline_color' , 'blue' ,  '[color] to outline point with for graph1'  ,
'begin',	'graph1_point_outline_width' , 1 ,  '[integer] thickness of point outline for graph1'  ,
'begin',	'graph1_point_size' , 2 ,  '[integer] pixel radius size for drawing graph1 points'  ,
'begin',	'graph1_tick_color' , 	'black' ,  '[color] of tick marks for graph1 scale'  , 
'begin',	'graph1_tick_length' , 	6 ,  '[integer] length of tick marks for graph1 scale'  , 
'begin',	'graph1_tick_offset' , 1 ,  '[integer] horizontal offset of tick marks for graph1 scale'  , 
'begin',	'graph1_tick_on' , 1 ,  '[0|1] toggles tick marks for graph1 scale off and on'  , 
'begin',	'graph1_tick_width' , 	3 ,  '[integer] thickness of tick marks for graph1 scale'  , 
'begin',	'graph1_vline_color' , 'black' ,  '[color] of vertical line for graph1 scale on left'  , 
'begin',	'graph1_vline_on' , 1 ,  '[0|1} toggles on and off vertical line for graph1 scale on left'  , 
'begin',	'graph1_vline_width' , 2 ,  '[integer] vertical line width for graph1 scale on left'  , 
'begin',	'graph2_label_color' , 'red' ,  '[color] of graph2 scale labels'  , 
'begin',	'graph2_label_decimal' , 2 ,  '[integer] number of decimal point to round graph2 scale label'  , 
'begin',	'graph2_label_fontsize' , 10 ,  '[integer] point size of graph2 scale labels'  , 
'begin',	'graph2_label_multiplier' , 1 ,  '[float] graph2 scale label multiplier'  , 
'begin',	'graph2_label_offset' , 8 ,  '[integer] horizontal offset of graph2 scale labels'  , 
'begin',	'graph2_label_on' , 1 ,  '[0|1] toggles graph2 scale labels off and n'  , 
'begin',	'graph2_line_color' , 'red' ,  '[color] of graph2 connecting lines'  , 
'begin',	'graph2_line_on' , 1 ,  '[0|1] toggles graph2 connecting lines off and on'  , 
'begin',	'graph2_line_smooth' , 0 ,  '[0|1] toggles graph2 connecting line smoothing off and on',
'begin',	'graph2_line_width' , 1 ,  '[integer] thickness of graph2 connecting lines'  , 
'begin',	'graph2_max' , 1000 ,  '[integer] maximum value for graph2 scale'  , 
'begin',	'graph2_min' , -1000 ,  '[integer] minimum value for graph2 scale'  , 
'begin',	'graph2_on' , 0 ,  '[0|1] toggles graph2_on'  , 
'begin',	'graph2_point_fill_color' , 'red' ,  '[color] of interior of graph2 points'  ,
'begin',	'graph2_point_on' , 1 ,  '[0|1] toggles graph2 point drawing on and off'  ,
'begin',	'graph2_point_outline_color' , 'red' ,  '[color] of graph2 point outline'  ,
'begin',	'graph2_point_outline_width' , 1 ,  '[integer] thickness of graph2 point outline'  ,
'begin',	'graph2_point_size' , 2 ,  '[integer] radius size of graph 2 points'  ,
'begin',	'graph2_tick_color' , 	'black' ,  '[color] of graph2 vertical scale ticks'  , 
'begin',	'graph2_tick_length' , 	6 ,  '[integer] length of graph2 vertical scale ticks'  , 
'begin',	'graph2_tick_offset' , 5 ,  '[integer] horizontal offset of graph2 vertical scale ticks'  , 
'begin',	'graph2_tick_on' , 1 ,  '[0|1] toggles graph2 vertical scale ticks on and off'  , 
'begin',	'graph2_tick_width' , 	3 ,  '[integer] thickness of graph2 vertical scale ticks'  , 
'begin',	'graph2_vline_color' , 'black' ,  '[color] of graph2 vertical scale line'  , 
'begin',	'graph2_vline_on' , 1 ,  '[0|1] toggles graph2 vertical scale line off and on'  , 
'begin',	'graph2_vline_width' , 2 ,  '[integer] thickness of graph2 vertical scale line'  , 
'begin',	'graph_scale_height' , 80 ,  '[integer] pixel height of shared graph scale'  , 
'begin',	'graph_scale_hline_color', 	'black' 	,  '[color] of horizontal shared graph scale lines'  , 
'begin',	'graph_scale_hline_on' ,	 	1 , 		 '[0|1] toggles off and on the shared horizontal interval lines of the graph scales'  , 
'begin',	'graph_scale_hline_width' , 	1 ,  '[integer] width of shared horizontal shared graph scale lines'  , 
'begin',	'graph_scale_indent' , -20 ,  '[integer] indentation for placing gscale above (or even below) the sequence line'  , 
'begin',	'graph_scale_interval' , 4 ,  '[integer] number of intervals'  , 
'begin',	'graph_scale_on' , 0 ,  '[0|1] toggles off and on the graph scales'  , 
'begin',	'help_on',   1, '[0|1] toggles off and on the popup help messages',
'begin',	'help_wrap',   50, '[integer] line width in characters for popup help menus',
'begin',	'mark_advanced' ,  ''  ,  "code for an advanced marking algorithm. Allowing for more complex searches.  Data foreach pair or extra is accessed using an array reference \$c.  Therefore to access column 4 \$\$c[4] would work."  , 
'begin',	'mark_array' , 'm',  '[e|m] default array to search (m is alignment/e is extra)(m is historical)'  , 
'begin',	'mark_col' ,  ''  ,  '[integer] column to search for given pattern in order to mark matches with a color'  , 
'begin',	'mark_col2' ,  ''  ,  '[integer] second column to search for pattern in order to mark matches with a color'  , 
'begin',	'mark_color' ,  'red'  ,  '[color] to  mark objects with'  , 
'begin',	'mark_pairs' ,  0  ,  '[0|1] toggles the coloring/marking of sub(jects) off and on'  , 
'begin',	'mark_pattern' ,  'AC002038'  ,  '[regular expression] pattern to search for with mark/find button'  , 
'begin',	'mark_permanent' ,  0  ,  '[0|1] toggles on and off changing the color of objects permanently (if not permanent then on redraw colors will be erased'  , 
'begin',	'mark_subs' ,  1  ,  '[0|1] toggles the coloring/marking of sub(jects) off and on'  , 
'begin',	'pair_inter_color' , 'red' ,  '[color] default of inter pairwise and connecting lines'  , 
'begin',	'pair_inter_line_on' , 0 ,  '[0|1] toggles off and on the connecting lines between inter pairwise alignments'  , 
'begin',	'pair_inter_offset' , 0 ,  '[integer] default offset from sequence line of inter pairwise (up is negative, down is positive)'  , 
'begin',	'pair_inter_on' , 1 ,  '[0|1] toggles off and on the inter pairwise alignments normally drawn on top of sequence line'  , 
'begin',	'pair_inter_width' , 13 ,  '[integer] width of inter pairwise'  , 
'begin',	'pair_intra_color' , 'blue' ,  '[color] default of intra pairwise and connecting lines'  , 
'begin',	'pair_intra_line_on' , 0 ,  '[0|1] toggles connecting lines between  intra pairwise off and on'  , 
'begin',	'pair_intra_offset' , 0 ,  '[integer] default offset from seuqence'  , 
'begin',	'pair_intra_on' , 1 ,  '[0|1] toggles off and on the intra pairwise'  , 
'begin',	'pair_intra_width' , 9 ,  '[integer] width of intra pairwise'  , 
'begin',	'pair_level' , 'NONE' ,  '[NONE|inter_over_intra|intra_over_inter] determines which pairwise type appears above the other--NONE leaves the appearance to the order of the pairwise in the inputted alignment or parasight.psa table'  , 
'begin',	'pair_type_col' , '' ,  '[integer] column number to determine pairwise type for sequence 1, which is checked against sequence 2.  If match then intra if no match then inter.  (Useful on sequence names that contain chromosome assignment.)'  , 
'begin',	'pair_type_col2' , '' ,  '[integer] column to determine pairwise type for sequence 2 in row '  , 
'begin',	'pair_type_col2_pattern' , '' ,  '[regular expression] to extract pairwise type determing value with parentheses'  , 
'begin',	'pair_type_col_pattern' , '' ,  '[regular expression] to extract pairwise type determining value with parentheses'  , 
'begin',	'popup_format' , 'text' ,  '[text|number] determines whether column numbers or text headers are shown in popup window'  , 
'begin',	'popup_max_len' , 300 ,  '[integer] character length for fields in the popup menu (allows long definitions or sequences be excluded)'  , 
'begin',	'print_command' , 'lpr -P Rainbow {}' ,  '[string] print command with brackets {} representing file name.  This is a system command executed to drive a printer.  I have never been able to get DOS to work.  This is setup for Unix on our system.  Rainbow is our color printer name.  It will fail in MSWin'  , 
'begin',	'print_multipages_high' , 1 ,  '[integer] height in number of pages for the print/postscript all command'  , 
'begin',	'print_multipages_wide' , 1 ,  '[integer] width in number of pages for print/postscript all command'  , 
'begin',	'printer_page_length' , '11i' ,  '[special] physical page length (longest dimension of paper) in inches for printer (requires number followed by units with i=inches or c=cm)'  , 
'begin',	'printer_page_orientation' , 1 ,  '[0|1] toggles printer page orientation (1=landscape 0=portrait)'  , 
'begin',	'printer_page_width' , '8i' ,  '[special] physical page width in inches for printer (requires number followed by units i=inches or c=cm)'  , 
'begin',	'quick_color' , 'purple' ,  '[color] for the quick color function Shift-Button3 and Shift-Double Click Button3'  , 
'begin',	'seq_color' , 'black' ,  '[color] of sequence (All sequences take this color.  There is currently no way to color sequences individually.)'  , 
'begin',	'seq_label_color' , 'black' ,  '[color] of sequence name text'  , 
'begin',	'seq_label_fontsize' , 12 ,  '[integer] font size (in points) for all sequence names'  , 
'begin',	'seq_label_offset' , -4 ,  '[integer] vertical offset of sequence names (up is negative, down is positive)'  , 
'begin',	'seq_label_offset_h' , 0 ,  '[integer] horizontal offset of sequence names'  , 
'begin',	'seq_label_on' , 1 ,  '[0|1] toggles off and on the display of sequence name labels'  , 
'begin',	'seq_label_pattern' , '' ,  '[regular expression] to match in sequence name for display purposes--parentheses must be used to denote the part of match to display'  , 
'begin',	'seq_line_spacing_btwn' , 250 ,  '[integer] pixels to separate  sequence lines from each other (roughly equivalent to spacing between text paragraphs if you consider a wrapping line of sequences to be a paragraph)'  , 
'begin',	'seq_line_spacing_wrap' , 200 ,  '[integer] pixels to space between a wrapping line of sequences (roughly equivaelent to spacing between the lines within a text paragraph)'  , 
'begin',	'seq_spacing_btwn_sequences' ,  10000 ,  '[integer] bases to separate sequences drawn within the same line (roughly equivalent to spacing between words of a text paragraph)'  , 
'begin',	'seq_tick_b_color' , 'black' ,  '[color] for begin tick marks'  , 
'begin',	'seq_tick_b_label_anchor' , 'ne' ,  '[center|n|w|s|e|nw|ne|sw|se] anchor point for begin tick mark labels'  , 
'begin',	'seq_tick_b_label_color' , 'black' ,  '[valid color] of tick mark label at the beginning of sequence'  , 
'begin',	'seq_tick_b_label_fontsize' , 9 ,  '[integer] font size (in points) for label at beginning of sequence'  , 
'begin',	'seq_tick_b_label_multiplier' , 0.001 ,  '[float] scaling factor for begin tick mark labels'  , 
'begin',	'seq_tick_b_label_offset' , 2 ,  '[integer] vertical offset for begin tick mark label'  , 
'begin',	'seq_tick_b_label_offset_h' , 0 ,  '[integer] horizontal offset for begin tick mark labels'  , 
'begin',	'seq_tick_b_label_on' , 1 ,  '[0|1] toggles  off and on the beginning tick mark labels'  , 
'begin',	'seq_tick_b_length' , 10 ,  '[integer] length of begin tick marks'  , 
'begin',	'seq_tick_b_offset' , 0 ,  '[integer] vertical offset for begin tick marks'  , 
'begin',	'seq_tick_b_on' , 1 ,  '[0|1] toggles off and on the begin tick marks'  , 
'begin',	'seq_tick_b_width' , 2 ,  '[integer] width of begin tick marks'  , 
'begin',	'seq_tick_bp' , 20000 ,  '[integer] tick mark interval'  , 
'begin',	'seq_tick_color' , 'black' ,  '[color] of interval tick marks'  , 
'begin',	'seq_tick_e_color' , 'black' ,  '[valid color] for end tick marks'  , 
'begin',	'seq_tick_e_label_anchor' , 'nw' ,  '[center|n|w|s|e|nw|ne|se|sw] anchor point for end tick mark labels'  , 
'begin',	'seq_tick_e_label_color' , 'black' ,  '[valid color] for end tick mark labels'  , 
'begin',	'seq_tick_e_label_fontsize' , 9 ,  '[integer] font size (in points)  for end tick mark labels'  , 
'begin',	'seq_tick_e_label_multiplier' , 0.001 ,  '[float] scaling factor for end tick mark labels'  , 
'begin',	'seq_tick_e_label_offset' , 2 ,  '[integer] vertical offset for end tick mark labels'  , 
'begin',	'seq_tick_e_label_offset_h' , 0 ,  '[integer] horizontal offset for end tick mark labels'  , 
'begin',	'seq_tick_e_label_on' , 1 ,  '[0|1] toggles end tick labels off and on'  , 
'begin',	'seq_tick_e_length' , 10 ,  '[integer] length of end tick marks'  , 
'begin',	'seq_tick_e_offset' , 0 ,  '[integer] vertical offset for ending tick marks'  , 
'begin',	'seq_tick_e_on' , 1 ,  '[0|1] toggles off and on the ending tick marks'  , 
'begin',	'seq_tick_e_width' , 2 ,  '[integer] width of end tick marks'  , 
'begin',	'seq_tick_label_anchor' , 'n' ,  '[center|n|s|w|e|nw|sw|ne|se] anchor of text from tick mark draw point'  , 
'begin',	'seq_tick_label_color' , 'black' ,  '[color] for interval tick mark'  , 
'begin',	'seq_tick_label_fontsize' , 9 ,  '[integer] font size (in points) for interval tick mark label'  , 
'begin',	'seq_tick_label_multiplier' , 0.001 ,  '[float] scaling factor for the interval tick label'  , 
'begin',	'seq_tick_label_offset' , 2 ,  '[integer] vertical offset of sequence interval tick mark labels'  , 
'begin',	'seq_tick_label_on' , 1 ,  '[0|1] toggles off and on the interval tick labels'  , 
'begin',	'seq_tick_length' , 10 ,  '[integer] length of interval tick marks'  , 
'begin',	'seq_tick_offset' , 0 ,  '[integer] vertical offset for interval tick marks'  , 
'begin',	'seq_tick_on' , 1 ,  '[0|1] toggles off and on the interval sequence tick marks'  , 
'begin',	'seq_tick_whole', 0, '[0|1] toggles whether numbering is for each individual sequence (0) or continious across multiple accession on same line (useful when analyzing chromosomes in multiple fragments)' ,
'begin',	'seq_tick_width' , 2 ,  '[integer] width of interval tick marks'  , 
'begin',	'seq_width' , 3 ,  '[integer] width of sequence line'  , 
'begin',	'showqueryonly' , 0 ,  '[0|1] toggles the display of just the first sequence in a pairwise data (i.e.first column in an alignment file).  For most parsing this is equivalent to the Blast query position'  , 
'begin',	'sub_arrow_diag' , 5 ,  '[integer] distance between arrow point to wing/edge of arrow'  , 
'begin',	'sub_arrow_on' , 0 ,  '[0|1] toggles off and on the directional/orientation arrows for subjects'  , 
'begin',	'sub_arrow_paral' , 5 ,  '[integer]  distance between arrow point to base of arrow'  , 
'begin',	'sub_arrow_perp' , 4 ,  '[integer] distance from base end to wing tip of arrow '  , 
'begin',	'sub_color' , 'lightgreen' ,  '[color] default of sub(ject) objects (all other coloring schemes over ride default)'  , 
'begin',	'sub_initoffset' , 30 ,  '[integer] pixel indent from top of subscales to associated sequence line (increasing pushes scales further below associated sequence)'  , 
'begin',	'sub_labelhit_col' , 13 ,  '[integer] column to use for labeling each hit/pairwise (label will be drawn at beginning of each hit sub)'  , 
'begin',	'sub_labelhit_color' , 'black' ,  'color of pairwise hit label text'  , 
'begin',	'sub_labelhit_offset' , 0 ,  '[integer] horizontal offset for hit label'  , 
'begin',	'sub_labelhit_on' , 0 ,  '[0|1] turns on individual labeling of each pairwise hit'  , 
'begin',	'sub_labelhit_pattern' , '0?([0-9.]{4})' ,  '[regular expression] to match in data from column'  , 
'begin',	'sub_labelhit_size' , 9 ,  '[integer] font size (in points) for hit label'  , 
'begin',	'sub_labelseq_col' , 0 ,  '[integer] column to use for the beginning sub label'  , 
'begin',	'sub_labelseq_col2' , 4 ,  '[integer] column for second position sequence in alignment table pairwise row'  , 
'begin',	'sub_labelseq_col2_pattern' , '' ,  '[regular expression] pattern to match in data from sub label sequence column 2'  , 
'begin',	'sub_labelseq_col_pattern' ,'',  '[regular expression] pattern to match in data from sub label sequence column (use parenthesis to denote data within the match to display)' ,  
'begin',	'sub_labelseq_color' , 'black' ,  '[color] of text label for sub objects'  , 
'begin',	'sub_labelseq_offset' , 0 ,  '[integer] horizontal offset label'  , 
'begin',	'sub_labelseq_on' , 1 ,  '[0|1] toggles overall begin sequence label for sub(ject) label off and on'  , 
'begin',	'sub_labelseq_size' , 6 ,  '[integer] font size (in points) for begin label sequence'  , 
'begin',	'sub_labelseqe_col' , 4 ,  '[integer] column to use for the end subject label'  , 
'begin',	'sub_labelseqe_col2' , 0 ,  '[integer] column for second position in alignment table pairwise row'  , 
'begin',	'sub_labelseqe_col2_pattern' , '' ,  '[regular expression] pattern to match in data from column'  , 
'begin',	'sub_labelseqe_col_pattern' , '' ,  '[regular expression] pattern to match in data from column'  , 
'begin',	'sub_labelseqe_color' , 'black' ,  '[valid color] of label text'  , 
'begin',	'sub_labelseqe_offset' , 0 ,  '[integer] horizontal offset for label'  , 
'begin',	'sub_labelseqe_on' , 0 ,  '[0|1] toggles off and on the overall sub(ject) label at end of last hit/pairwise'  , 
'begin',	'sub_labelseqe_size' , 6 ,  '[integer] font size (in points) for end subject label'  , 
'begin',	'sub_line_spacing' , 9 ,  '[integer] pixels per line determining the spacing between subs placed on different lines'  , 
'begin',	'sub_on' , 1 ,  '[0|1] toggles sub(ject) display off and on (these are the pairwise representations drawn below the sequence line) For BLAST searches these traditionally represent the subject sequences found in a database search.'  , 
'begin',	'sub_scale_categoric_string' , '' ,  '[string] list of comma delimited category names'  , 
'begin',	'sub_scale_col' , '' ,  '[integer] column for value to arrange pairwise hit on sub scale (subscale)'  , 
'begin',	'sub_scale_col2' , '' ,  '[integer] column for second position sequence in alignment pairwise (only used if defined)'  , 
'begin',	'sub_scale_col2_pattern' , '' ,  '[regular expression] pattern to match in column 2'  , 
'begin',	'sub_scale_col_pattern' , '' ,  '[regular expression] pattern to match in column'  , 
'begin',	'sub_scale_hline_color' , 'grey' ,  '[valid color] for horizontal sub scale lines'  , 
'begin',	'sub_scale_hline_on' , 1 ,  '[0|1] toggles off and on the horizontal scale lines for sub scale'  , 
'begin',	'sub_scale_hline_width' , 1 ,  '[integer] width of horizontal sub scale lines'  , 
'begin',	'sub_scale_label_color' , 'black' ,  '[color] for sub scale axis label'  , 
'begin',	'sub_scale_label_fontsize' , 12 ,  '[integer] font size (in points) for sub scale axis label'  , 
'begin',	'sub_scale_label_multiplier' , 100 ,  '[integer] multiplication factor for sub scale label'  , 
'begin',	'sub_scale_label_offset' , 1 ,  '[integer] horizontal offset for sub scale axis tick marks'  , 
'begin',	'sub_scale_label_on' , 1 ,  '[0|1] toggles off and on sub scale axis tick mark labels'  , 
'begin',	'sub_scale_label_pattern' , '' ,  '[regular expression] pattern to match in sub scale label'  , 
'begin',	'sub_scale_lines' , 10 ,  '[integer] number of lines (or interval steps) to plot for stagger or cscale (automatically set for subscaleC)'  , 
'begin',	'sub_scale_max' , 1.00 ,  '[float] maximum value to place on the sub scale (automatically set for subscaleC)'  , 
'begin',	'sub_scale_min' , 0.80 ,  '[float] minimum value to place on the sub scale (automatically set for subscaleC)'  , 
'begin',	'sub_scale_on' , 0 ,  '[0|1] toggles sub scale on and off'  , 
'begin',	'sub_scale_step' , 0.01 ,  '[float] value to increment between each step (automatically set to -1 for subscaleC, 1 reverses subscaleC)'  , 
'begin',	'sub_scale_tick_color' , 'black' ,  '[color] for sub scale axis tick marks'  , 
'begin',	'sub_scale_tick_length' , 9 ,  '[integer] length of sub axis tick marks'  , 
'begin',	'sub_scale_tick_offset' , 4 ,  '[integer] offset of sub scale axis tick marks'  , 
'begin',	'sub_scale_tick_on' , 1 ,  '[0|1] toggles off and on the sub scale axis at horizontal tick positions'  , 
'begin',	'sub_scale_tick_width' , 3 ,  '[integer] width of sub scale axis tick marks'  , 
'begin',	'sub_scale_vline_color' , 'black' ,  '[color] for vertical axis line of sub scale'  , 
'begin',	'sub_scale_vline_offset' , -5 ,  '[integer] horizontal offset for subject axis line'  , 
'begin',	'sub_scale_vline_on' , 1 ,  '[0|1] toggles off and on the vertical axis line for sub scale'  , 
'begin',	'sub_scale_vline_width' , 2 ,  '[integer] width of sub scale axis line'  , 
'begin',	'sub_width' , 8 ,  '[integer] default width (thickness) of sub objects'  , 
'begin',	'template_desc_on' , 1 ,  '[0|1] toggles off and on wether  descriptions, such as this one, are saved in a template file with each option variable'  , 
'begin',	'text2_anchor' , 'nw' ,  '[center|n|w|s|e|nw|ne|se|sw] anchor point for end tick mark labels'  , 
'begin',	'text2_color' , 'red' ,  '[color] for end tick mark labels'  , 
'begin',	'text2_offset' , 0 ,  '[integer] vertical offset for end tick mark labels'  , 
'begin',	'text2_offset_h' , 0 ,  '[integer] horizontal offset for end tick mark labels'  , 
'begin',	'text2_on' , 1 ,  '[0|1] toggles end tick labels off and on'  , 
'begin',	'text2_size' , 20 ,  '[integer] font size (in points)  for end tick mark labels'  , 
'begin',	'text2_text' , '' ,  '[text] to display within a parasight view (useful for automation)'  , 
'begin',	'text_anchor' , 'nw' ,  '[center|n|w|s|e|nw|ne|se|sw] anchor point for end tick mark labels'  , 
'begin',	'text_color' , 'red' ,  '[color] for end tick mark labels'  , 
'begin',	'text_fontsize' , 20 ,  '[integer] font size (in points)  for end tick mark labels'  , 
'begin',	'text_offset' , 0 ,  '[integer] vertical offset for end tick mark labels'  , 
'begin',	'text_offset_h' , 0 ,  '[integer] horizontal offset for end tick mark labels'  , 
'begin',	'text_on' , 1 ,  '[0|1] toggles end tick labels off and on'  , 
'begin',	'text_text' , '' ,  '[text] to display within a parasight view (useful for automation)'  , 
'begin',	'window_font_size' ,  9  ,  '[integer] font size for parasight in general (not implemented)'  , 
'begin',	'window_height', 550 ,  '[integer] pixel height of main window on the initial start up'  , 
'begin',	'window_width' , 800 ,  '[integer] pixel width of the main window on the initial start up'  , 
#
);

for(my $i=0; $i<@todefine; $i+=4) {
	my ($b,$o,$v,$d)=@todefine[$i..$i+3];
	die "Error in Option Array ($o,$v,$d,$b) $b ne begin !\n" if $b ne 'begin';
	$opt{$o} = $v if !defined $opt{$o};
	$optdesc{$o}=$d;
}


if ($options) {
	#print "$options\n";
	my @tmparray= split / *[=>,]+ */,$options;
	my %moreoptions=@tmparray;
	#print keys %moreoptions, "\n";
	foreach my $k (keys %moreoptions) {
		my $clean = substr($k,1);
		if (!defined $opt{$clean} ) {
			print "VALID OPTIONS are: ";
			foreach (sort keys %opt) { print "-$_ ";}
			print "\nBAD OPTION: -$clean is an invalid option for parasight.\nValid options are listed above.\n";
			exit;
		}
		#print "OPTIONS:$k\n";
		$opt{$clean}=$moreoptions{$k};
	}


}

%opt = (%opt,%newopt);

warn "Newline character not tolerated in -text_text.\nUse \\\\n if entering to get \\n n the input!\n" if $opt{'text_text'} =~/\n/mg;

###################################################################
###################################################################
############################ GUI CREATION #########################
###################################################################

###########################MAIN WINDOW ############################
###################################################################
	$mw=MainWindow->new;
	my $cwd=cwd();
	$optionpath="$cwd/";
	$filepath='';
	$filepath=$opt{'in'};
	$filepath="$cwd/" if !$filepath;
	if ($opt{'in'} ) {
		$mw->title("PARASIGHT: $opt{'in'}");
	} else {
		$mw->title("PARASIGHT: New");
	}
	$mw->setPalette('lightgrey');
	$mw->configure(-background=>'darkgrey');
	$balloon=$mw->Balloon(-initwait =>600, -background => '#ffff9d', -font=>"Courier 10" );

	#####frame #######
	$frame = $mw->Frame(-relief => 'groove', -bd => 2);
	
	
	$tmp=$frame->Menubutton(-text => "File",
				-menuitems=> [
							#['command' => "Load Parasight",
				#			   -command => sub{
				#						my ($dir,$name);
								#		$name=$filepath;
								#		($dir,$name)= ($1,$2) if $filepath=~ /^(.*)\/(.*)$/;
								#		my @filetypes= (['parasight', '.psa'],['All Files', '*']);
							   	#   my $file = $mw->getOpenFile(-title=> 'LOAD PARASIGHT FILES', -filetypes => \@filetypes,
									#					-initialdir=>$dir, -initialfile=>$name);
									#	if ($file eq "") { return; }
									#	$filepath=$file;
									#	&load_parasight_table($filepath);
									#	$opt{'in'}=$filepath;
									##}
							  #],
							  ['command' => "Saving Parasight",
								-command =>  sub{
										my ($dir,$name);
										$name=$filepath;
										print "$filepath";
										($dir,$name)= ($1,$2) if $filepath=~ /^(.*)\/(.*)$/;
										if ($^O=~/MSWin/) {$dir =~ s/\//\\/mg;}
										#print "POSITION FP:$filepath D($dir)  N($name)\n";
										my @filetypes= (['parasight', '.psa'],['All Files', '*']);
							   	   my $file = $mw->getSaveFile( -title=>'SAVE PARASIGHT FILES',
							   	   			-filetypes => \@filetypes, 
													 -initialdir=> $dir,
												   -initialfile=>$name
												  );
										if ($file eq "") { print "CANCELLED!\n" if $opt{'quiet'}; return; }
										$filepath=$file;
										print "SAVING: $filepath\n" if !$opt{'quiet'};
										&save_parasight_table($filepath);
										$opt{'in'}=$filepath;
										$mw->title("PARASIGHT: $filepath");
									}
								],
							  ['command' => "Load Option Template",
								-command =>  sub{
										my ($dir,$name);
										$name=$optionpath;
										($dir,$name)= ($1,$2) if $optionpath=~ /^(.*)\/(.*)$/;
										if ($^O=~/MSWin/) {$dir =~ s/\//\\/mg;}
										print "POSITION FP:$filepath D($dir)  N($name)\n";
										my @filetypes= (['option template', ['.pst']],['All Files', '*']);
							   	   my $file = $mw->getOpenFile( -title=>'LOAD OPTION TEMPLATE',
							   	   			-filetypes => \@filetypes, 
												    -initialdir=>$dir, -initialfile=>$name);
										if ($file eq "") { return; }
										$optionpath=$file;
										&load_option_template($optionpath);
									}
								],
							  ['command' => "Save Option Template",
								-command =>  sub{
										my ($dir,$name);
										$name=$optionpath;
										($dir,$name)= ($1,$2) if $optionpath=~ /^(.*)\/(.*)$/;
										if ($^O=~/MSWin/) {$dir =~ s/\//\\/mg;}
										my @filetypes= (['option template', ['.pst']],[	'All Files', '*']);
							   	   my $file = $mw->getSaveFile( -title=>'SAVE OPTION TEMPLATE',
							   	   		#	-filetypes => \@filetypes, 
												    -initialdir=>$dir, -initialfile=>$name);
										if ($file eq "") { return; }
										$optionpath=$file;
										&save_option_template($optionpath);
									}
							 	]])->pack(-side => 'left');
							 
	$balloon->attach($tmp,-justify => 'left',-msg=>"Click for menu to save and load\nparasight files\n"
																.  "  or option template files.");

	$tmp=$frame->Menubutton(-text => "Print",
				-menuitems=> [
				['command' => "Print Screen", -command => [\&print_screen,1]],
				['command' => "Postscript Screen", -command => [\&print_screen,0]],
				['command' => "Print All",-command =>  [\&print_all,1]],
				['command' => "Postscript All",-command =>  [\&print_all,0]]

			])->pack(-side => 'left');
	$tmp->separator;
	$tmp->checkbutton(-label=> 'landscape',, -variable => \$opt{'printer_page_orientation'});

	$balloon->attach($tmp,-justify => 'left',-msg=>"Click to see menu of output choices.\nOn Windows print may not work\n");
	$tmp=$frame->Menubutton(-text => "Order",
				-menuitems=> [
				['command' => "pair inter    => raise", -command => sub{$canvas->raise('inter');}],
				['command' => "pair inter    => lower", -command => sub{$canvas->lower('inter');}],
				['command' => "pair intra    => raise", -command => sub{$canvas->raise('intra');}],
				['command' => "pair intra    => lower", -command => sub{$canvas->lower('intra');}],
				['command' => "sub           => raise", -command => sub{$canvas->raise('sub');}],
				['command' => "sub           => lower", -command => sub{$canvas->lower('sub');}],
				['command' => "sub label     => raise", -command => sub{$canvas->raise('subl');}],
				['command' => "sub label     => lower", -command => sub{$canvas->lower('subl');}],
				['command' => "extra         => raise", -command => sub{$canvas->raise('ex');}],
				['command' => "extra         => lower", -command => sub{$canvas->lower('ex');}],
				['command' => "extra label   => raise", -command => sub{$canvas->raise('exl');}],
				['command' => "extra label   => lower", -command => sub{$canvas->lower('exl');}],
				['command' => "graph line    => raise", -command => sub{$canvas->raise('gl');}],
				['command' => "graph line    => lower", -command => sub{$canvas->lower('gl');}],


				['command' => "subscale        => raise", -command => sub{$canvas->raise('ss');}],
				['command' => "subscale        => lower", -command => sub{$canvas->lower('ss');}],
				['command' => "subscale label  => raise", -command => sub{$canvas->raise('ssl');}],
				['command' => "subscale label  => lower", -command => sub{$canvas->lower('ssl');}],
				['command' => "gscale        => raise", -command => sub{$canvas->raise('gs');}],
				['command' => "gscale        => lower", -command => sub{$canvas->lower('gs');}],
				['command' => "gscale label  => raise", -command => sub{$canvas->raise('gsl');}],
				['command' => "gscale label  => lower", -command => sub{$canvas->lower('gsl');}],
				['command' => "sequence      => raise", -command => sub{$canvas->raise('seq');}],
				['command' => "sequence      => lower", -command => sub{$canvas->lower('seq');}],
				['command' => "sequence name => raise", -command => sub{$canvas->raise('seqn');}],
				['command' => "sequence name => lower", -command => sub{$canvas->lower('seqn');}],
				['command' => "tick          => raise", -command => sub{$canvas->raise('tick');}],
				['command' => "tick          => lower", -command => sub{$canvas->lower('tick');}],
				['command' => "tick label    => raise", -command => sub{$canvas->raise('tickl');}],
				['command' => "tick label    => lower", -command => sub{$canvas->lower('tickl');}]
			])->pack(-side => 'left');
	$balloon->attach($tmp,-justify => 'left',-msg=>"Allows the drawn objects to be raised\nand lowered relative to each other.");
	$tmp=$frame->Menubutton(-text => "Misc",
				-menuitems=> [
				['command' => "color transfer sub -> pair", -command => sub{
							my $c=$mh{'color'};
							my $s=$mh{'scolor'};
							for(my $i=0;$i<@m;$i++) { $m[$i][$c]=$m[$i][$s] }
					}],
				['command' => "color transfer pair -> sub", -command => sub{
							my $c=$mh{'color'};
							my $s=$mh{'scolor'};
							for(my $i=0;$i<@m;$i++) { $m[$i][$s]=$m[$i][$c] }
					}],
				['command' => "color transfer pair -> sub", -command => sub{$canvas->raise('inter');}],
			])->pack(-side => 'left');
	$balloon->attach($tmp,-justify => 'left',-msg=>"A Hodge-Podge of Misc functions including\nraising and lowering objects relative to each other.");
	$tmp=$frame->Button(-text => 'Options',
				   -command => [\&indexcard_options])->pack(-side=>'left',-anchor => 'e');
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to Display PopUp Window of All Options");
	
	
	$tmp=$frame->Button(-text => 'PrintScreen',
				   -command => [\&print_screen,1])->pack(-side=>'left',-anchor => 'e');
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to Print Current Visible Screen View\n"
	        . "(Options->Misc to change printing function)");

	$tmp=$frame->Button(-text=>'DeZoom',-command=> sub{$canvas->scale("all",0,0,1/$scale,1/$scale);
											$canvas->configure(-scrollregion=>[$canvas->bbox("all")]);
											#print $canvas->cget(-height),"  ",$canvas->cget(-width),"\n";
											$scale=1;})
				->pack(-side=>'left',-anchor=>'e');
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to reset to normal scale after zooming.\n"
	         . "to Zoom In use Control-Button1 (Left-Click).\n"
	         .  "to Zoom Out use  Control-Button3 (Right-Click).");

	$tmp=$frame->Button(-text=>'FitLongLine', -command => \&fitlongestline  
		)->pack(-side=>'left',-anchor=>'e');
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to set bp width\nto longest sequence line.");


	$tmp=$frame->Button(-text => 'Find', -command=>\&mark_window)
				->pack(-side=>'left',-anchor => 'e');	
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to bring up window\nfor searching and color-marking results.");
	$tmp=$frame->Button(-text => 'Quick Color', -foreground=>$opt{'quick_color'},-background=>'white',
			-activebackground=>$opt{'quick_color'}, -activeforeground=>'white')
				->pack(-side=>'left',-anchor => 'e');	
	$tmp->configure(-command=>
						 [sub{ my $b= $_[0];
							  $mw->grabRelease();
							  my $color = $b->chooseColor(-title=>'Choose New Quick Color',
									-initialcolor=> $opt{'quick_color'});
							  if (defined $color) {
									$b->configure(-foreground=> $color,-activebackground=>$color);
									$opt{'quick_color'}=$color;
							  }
							  $mw->grab();
							  $mw->raise();
					  		}, $tmp] );
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to set a quick color.\n"
	        ."To quickly color a seq, pair, sub or extra object\n use Shift-Right-Click\n"
	        ."To remove color (return to default) ===> Shift-Double-Right-Click \n"
	        ."(NOTE:for pairs and labels color where information is unknown, removing color leaves it black until a redraw!");
#	$tmp=$frame->Button(-text => 'Help')
#				->pack(-side=>'left',-anchor => 'e');	
#	$tmp->configure(-command=>
	#					 sub {my $text="This is the simple text version of the internal #documentation\nBetter foramts can be accessed with parsight -h, perldoc, or pod2html\n";
	#					 	$text.=`pod2text $0`;
	#					 	print "$0\n";
	#					 	&export_text(\$text,"Internal Help Documenation ($0)");
	#				  		} );
	#$#balloon->attach($tmp,-justify => 'left',-msg=>"Press to set  quick color.\n"
	#        ."To quickly color a seq, pair, sub or extra object\n ==>Shift-Right-Click\n"
	#        ."To remove color (return to default) ===> Shift-Double-Right-Click \n"
	#        ."(NOTE:for pairs and labels color were information is unknown, removing color leaves black!");
	$tmp=$frame->Button(-borderwidth => 3,-activebackground=>'black',-activeforeground=>'white',
				-background=>'white',-command=>\&redraw,-text => "Redraw")->pack(-side=>'right');
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to redraw the entire screen.\n"
	        ."This is quicker then the blue button,\nbut any blue options may not be changed."
				);
	$tmp=$frame->Button(-borderwidth => 3,-activebackground=>'blue',-background=>'#bbe8ff',
			-command=>\&reshowNredraw,-text => "R R & R")->pack(-side=>'right');
	$balloon->attach($tmp,-justify => 'left',-msg=>"Press to Reshow, Rearrange and Redraw\n"
	        ."This is slower than the just Redraw\nas many initial calculations are redone."
			);

	$frame->pack(-side => 'top', -expand=> 0,-fill=>'x', -anchor=>'w');	
	##########################################################
	##########################################################
	############## create a canvas ###########################
	$scrolledcanvas = $mw->Scrolled('Canvas',-height=>$opt{'window_height'}, -width=>$opt{'window_width'}, -background => 'white')
			->pack(-side=>'top',-fill => 'both', -expand => 1, -anchor =>'n');
	##could add 
	##$scrolledcanvas->Tk::bind("<Configure>", [sub { print "SIZE CHANGED $_[1], $_[2]\n"; }, Ev{'h'}, Ev{'w'} ]);
	$canvas= $scrolledcanvas->Subwidget('canvas');

############################################################
#############  highlighting on mouse over #############
##########################################################
	my ($oldhighlightcolor,$oldhighlightid);
	my @highlightcolors=('yellow', 'orange','yellow','pink','yellow','lightblue');
	$canvas->Tk::bind('<Motion>',
			[sub{
					my ($canv, $x, $y)=@_;
					my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );

					#my @close=$canv->find("closest",$cx,$cy) ;
					my @tags = $canv->gettags( "current");
					my $id=$tags[0];
					foreach my $tag (@tags) {
						next if $tag !~/^[MSE]/;  #currently sequence gets sS  not S
						$id=$tag;
					}
					#print "$id\n";
					my $color=$canvas->itemcget($id,-fill);
					if ($oldhighlightid ne $id && $oldhighlightcolor ne '') {
						#print "FIX $oldhighlightid ($id) $oldhighlightcolor=>$color\n";
						$canvas->itemconfigure($oldhighlightid ,-fill=>$oldhighlightcolor);
						$oldhighlightid = '';
					}
					if ($id =~ /^[MSE]/ && $color ne '' ) {
							#print   "$oldhighlightid ($id) $color -> oldcolor\n";
							$oldhighlightcolor=$color if $oldhighlightid ne $id;
							$oldhighlightid=$id;
							#print "ENTER:$id:$c\n";
							push @highlightcolors, shift @highlightcolors;
							$canvas->itemconfigure( $id,-fill=>$highlightcolors[0]);					
					}
				}
		,Ev('x'),Ev('y')]);	



	#############################################
	########### popup window ####################

	$canvas->Tk::bind('<Button-1>',
				[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						#my @close=$canv->find("closest",$cx,$cy) ;
						#my @tags = $canv->gettags( "$close[0]");
						my @tags=$canv->gettags("current");
						my $id=$tags[0];
						foreach my $tag (@tags) {
							next if $tag !~/^[MSE]/;
							$id=$tag;
						}
						my ($type, $numb) = $id=~/([A-Za-z]+)(\d+)/;
						my $text.="($id) No Associated Data\n";
						if ($type =~ /^M|S/) {
							my @r=@{ $m[$numb] };  #find the row of the data assocated with this item#
							my $span1=$r[2]-$r[1]+1;
							my $span2=$r[6]-$r[5]+1;
							$span2=-($r[5]-$r[6]+1) if ($r[5]>$r[6]);
							
							for (1..3,5..7) {next if !/[0-9.]+$/; $r[$_]=&commify($r[$_]);}
							$text="$id\n";
							$text.= "Sa: $r[0]  $r[1] - $r[2] ($span1) len:$r[3]\n";
							$text.= "Sb: $r[4]  $r[5] - $r[6] ($span2) len:$r[7]\n";
							for (my $j=8; $j<@r; $j++) { 
								next if $r[$j] eq '';
								my $t=$r[$j];
								if (length($t)> $opt{'popup_max_len'} ) {
									#print "LEN:,",length($t),"\n";
									$t='(too long)';
								}
								if ($opt{'popup_format'} eq 'text') {
									$text.="\[$mheader[$j]\]$t " ;
								} else {
									$text.="\[$j\]$t " ;
								}
							}		
						} elsif ($type =~/^E/ ) {
							my @r=@{ $e[$numb] };
							for (my $i=1;$i<3;$i++) {$r[$i]=&commify($r[$i]);}
							$text="$id\n";
							$text.="$r[0] $r[1]-$r[2]\n"; 
							for (my $j=3; $j<@r; $j++) { 
								next if $r[$j] eq '';
								my $t=$r[$j];
								if (length($t)>$opt{'popup_max_len'}) {
									$t='(too long)';
								}
								if ($opt{'popup_format'} eq 'text') {
									$text.="\[$eheader[$j]\]$t " ;
								} else {
									$text.="\[$j\]$t " ;
								}
							}		
						}
						next if $text=~/No Associated Data/m;
						my $displayedcolor=$canvas->itemcget($id, -fill);
						my $poptext=$canvas->createText($cx, $cy,,-anchor=>'nw',-justify=>'left',-width=>400,-fill=>'black',-text=>$text);
						my ($l,$r,$t,$b)=$canvas->bbox($poptext);
						my $poprect=$canvas->createRectangle($l-2, $r-2, $t+2, $b+2,,-fill=>'#ffff9d');
						$canvas->addtag("POP$poprect",'withtag',$poprect);
						$canvas->addtag("POP$poprect",'withtag',$poptext);
						$canvas->lower($poprect,$poptext);
						$canvas->bind("POP$poprect","<Leave>", sub {$canvas->delete("POP$poprect"); } );
						$canvas->bind("POP$poprect","<Shift-B1-Motion>",
									[sub{
										my ($canv, $x, $y)=@_;
										my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
										$canvas->move( "POP$poprect", $cx - $iinfo{'lastX'}, $cy-$iinfo{'lastY'} );
										$iinfo{'lastX'}=$cx;  
										$iinfo{'lastY'}=$cy;
									}
									,Ev('x'),Ev('y')]);	
						
					}
			,Ev('x'),Ev('y')]);	
	#######################################################33
	################# popup menu ############################
	$canvas->Tk::bind('<Button-3>',
			[sub{
				my ($canv, $x, $y)=@_;
				#my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
				#print "$cx $cy =>";
				#$canv->delete ("2");
				#my @close=$canv->find("closest",$cx,$cy) ;

				#my @close=$canv->find("all") ;
				#print join ("X",@close),"X\n";
				#$canv->itemconfigure( $close[0],-fill=>"yellow");
				
				my @tags = $canv->gettags( "current");
				#print join ("X",@tags),"X\n";
				my $id=$tags[0];
				foreach my $tag (@tags) {
					next if $tag !~/^[MSE]/;
					$id=$tag;
				}
				my ($type, $numb) = $id=~/([A-Za-z]+)(\d+)/;
				my ($apnt,$hhashp, $hpnt);
				if ( $type=~/M|S/) { 
					$apnt = \@m; $hhashp =\%mh; $hpnt=\@mheader;
				} elsif ($type =~/E/) {
					$apnt =\@e; $hhashp =\%eh; ; $hpnt =\@eheader;
				} else {
					return;
				}
				my $menu = $canvas->Menu(-relief => 'groove', -tearoff => 0,
				-menuitems => [ 	['command' => "**$id**"],
										['command' => 'color', -command => [\&color_change, $id, $apnt, $hhashp] ],
										['command' => 'edit', -command => [\&edit, $id,$apnt,$hpnt] ],
										#['command' => 'quick_capture', -command => [\&quick_capture, $id, $apnt, $hhashp] ],
									#['command' => 'quick_edit',-command => [\&quick_edit, $id, $apnt, $hhashp] ]
					# 
						] );
				if ( $type=~/^M|S/ ) { 
					$menu->command ( -label => "decontaminate via seq 1 ($m[$numb][0])",
							-command=> [\&decontaminate_high_copy_repeat,$m[$numb][0],$m[$numb][1],$m[$numb][2] ]);
					$menu->command ( -label => "decontaminate vua seq 2 ($m[$numb][4])",
							-command=> [\&decontaminate_high_copy_repeat,$m[$numb][4],$m[$numb][5],$m[$numb][6] ]);
					if ($opt{'alignment_col'} !=0 && $opt{'alignment_col2'} !=0) {
						$menu->command ( -label => 'show alignment', -command=> [\&show_alignment, $numb]  );
						$menu->command ( -label => '  aln seq 1 w   -',-command=> [\&alignment_internal, $numb,0,'with']  );
						$menu->command ( -label => '            w/o -',-command=> [\&alignment_internal, $numb,0,'without']  );
						$menu->command ( -label => '  aln seq 2  w   -',-command=> [\&alignment_internal, $numb,4,'with']    );
						$menu->command ( -label => '             w/o -' ,-command=> [\&alignment_internal, $numb,4,'without'] );
					}
					if ($opt{'fasta_on'}==1) {
						$menu->command ( -label => "seq 1 external ($m[$numb][0])",
							-command=> [\&extract_sequence,$m[$numb][0],$m[$numb][1],$m[$numb][2] ]);
						$menu->command ( -label => '         modify +/- F/R',
							-command=> [\&extract_sequence,$m[$numb][0],$m[$numb][1],$m[$numb][2],'modify' ]);
						$menu->command ( -label => "seq 2 external ($m[$numb][4])",
							-command=> [\&extract_sequence,$m[$numb][4],$m[$numb][5],$m[$numb][6] ]);
						$menu->command ( -label => '         modify +/- F/R' ,
							-command=> [\&extract_sequence,$m[$numb][4],$m[$numb][5],$m[$numb][6],'modify']);
					}
				}
				
				if ( $type=~/^E/ ) { 
					$menu->command ( -label => "seq 1 external ($m[$numb][0])",
						-command=> [\&extract_sequence,$e[$numb][0],$e[$numb][1],$e[$numb][2] ]);
					$menu->command ( -label => '         modify +/- F/R',
						-command=> [\&extract_sequence,$e[$numb][0],$e[$numb][1],$e[$numb][2],'modify' ]);
				}
				
				my $cnt=1;
				foreach my $n (('',2,3,4)) {
					my $command=$opt{"execute$n"};
					my $desc=$opt{"execute$n"."_desc"};
					my $array=$opt{"execute$n"."_array"};
					my $col=$opt{"execute$n"."_col"};
					#print "x$n  xx$type  xxx$array\n";
					next if $type eq 'E' && $array ne 'e';
					next if $type ne 'E' && $array ne 'm';
					next if $command eq '';
					$menu->command ( -label => "CMD$cnt)$desc" ,
						-command=> sub { &execute_execute($id,$n)} 
					  );
					$cnt++;
				}
				
				
 	$menu-> Popup(-popover => 'cursor', -popanchor => 'nw');

   break;					}
			,Ev('x'),Ev('y')]);


sub alignment_internal {
	my $numb=shift;
	my $column=shift;
	my $type=shift;
	my $seq_col='';
	my $header='';
	my $title='';
	if ($column==0) {
		$seq_col=$opt{'alignment_col'};
		$title.= "First (query) alignment sequence";
		$header= ">$m[$numb][0].$m[$numb][1]-$m[$numb][2]\n";
	} elsif ( $column==4) {
		$title.= "Second (sujbect) alignment sequence";
		$header= ">$m[$numb][4].$m[$numb][5]-$m[$numb][6]\n";
		$seq_col=$opt{'alignment_col2'};
	} else {
		warn "BAD USE OF alignment_internal subroutine\n";
	}
	my $seq=$m[$numb][$seq_col];
	#print "($seq)\n";
	if ($type =~ /without/ ) {
		$seq=~ tr/-//d;
		$title .= ' without dashes';
	} else {
		$title .= ' with dashes';
	}
	#print "$header($seq)\n";
	$title .= " for M$numb";
	$seq=&fasta_format_wrap($seq,50);
	$seq = $header. $seq;
	if ($seq_col==0) {
		$seq= "Prealigned sequences must be included as columns\nin the align file for this option to work!\nThese sequences must contain indel dashes\nas the alignment is  not recalculated."
	}
	
	&export_text( \$seq, $title);         
	
}
	

sub extract_sequence {
	my $seqname=shift;
	my $begin=shift;
	my $end=shift;
	my $modify=shift;# modify=modify
	my $orient='F';
	if ($begin>$end) {
		$orient ='R';
		($end,$begin)=($begin,$end);
	}
	if ($modify) {
		my ($b,$e) = ($begin,$end);
		my $mwx = new MainWindow;
		$mwx->title("MODIFY SEQUENCE EXTRACTION");
		$mwx->setPalette('lightgrey');
		$mwx->configure(-background=>'darkgrey');
		my $f = $mwx->Frame(
				-borderwidth => 4,
				-relief => 'raised',
		  )->pack(-expand=>1,-fill=>'both');
		$f->Label(-text => "$seqname")->pack(-side=> 'left',-anchor => 'e');
		$f->Label(-text => " begin")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$b ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');
		$f->Label(-text => " end")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$e ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');
		$f->Optionmenu(-textvariable=>\$orient, -options => ['F','R'] )->pack(-side => 'left',-anchor=>'e');
		$f->Button(-borderwidth => 3,-activebackground=>'white',-background=>'#bbe8ff',
				-command=>sub { $mwx->destroy;},-text => "Extract")->pack(-side=>'right');
		$f->Button(-borderwidth => 3,-activebackground=>'white',-background=>'#bbe8ff',
				-command=>sub { $b=0; $e=0; $mwx->destroy;},-text => "Cancel")->pack(-side=>'right');
		$mwx->waitWindow;
		$begin = $b if $b;
		$end = $e if $e;

	}
	print "$seqname ($begin-$end) ($orient)\n" if !$opt{'quiet'};
	my $seqpath='';
	my $fragged=0;
	foreach my $fdir ( (split /:/,$opt{'fasta_directory'}) ) {
		#print "DIRECTORY:($fdir)\n";
		if (-d $fdir && opendir (DIR, $fdir) ) {
			my @matches=grep { /^$seqname/ } readdir DIR;
			my @single= grep {/^$seqname$/ } @matches;
			if (@single ==1 ) {
				#print "Perfect file ($single[0]) \n";
				$seqpath="$fdir/$single[0] ";
			} elsif ( $matches[0]=~ /^$seqname.?_\d+$/ ) {
				#print "FRACTIONATED $matches[0]\n";
				$seqpath="$fdir/$matches[0]";
				$fragged=1;
			}
		}
		last if $seqpath;
	
	}
	#search in blastdbs with fastacmd#
	if (! $seqpath ) {
		my $db = $opt{'fasta_blastdb'};
		$db =~ s/:/ /mg;
		#hunt for it in a fastacmd#
		if (open (FASTA, "fastacmd -d '$db' -s $seqname |") ) {
			my $head=<FASTA>;
			
			if ($head =~ />/ && $head =~/$seqname/) {
				#print "EXTRACTING SEQUENCE FROM BLASTDB\n";
				open (OUTFASTA, ">tmpfasta") || die "Can't create tmpfasta\n";
				print OUTFASTA $head;
				while (<FASTA>) {
					print OUTFASTA;
				}
				$seqpath='tmpfasta';
			}
		}
	}
	print "SEARCH RESULT:$seqname ($seqpath) $begin-$end FRAG($fragged)\n" if !$opt{'quiet'};
	my $seq='';
	return if $seqpath eq '';
	if ($fragged == 0 ) {
		print "$seqpath $begin $end\n" if !$opt{'quiet'};;
		$seq=&fasta_getsubseq_whole($seqpath,$begin,$end);
	} else {
		#get fragged sequence#
		$seqpath=~s/_\d+$//;
		$seq=&fasta_getsubseq_frac($seqpath,$begin,$end,$opt{'fasta_fragsize'});
	}
	#print "($seq)\n";
	if ($seq eq '') {
		print "ERROR:No sequence extracted";
		return;
	}
	$seq=&fasta_format_wrap($seq,50);
	$seq = ">$seqname.$begin.$end.$orient\n$seq";
	
	&export_text( \$seq, "sequence extracted", "$seqname.$begin.$end.$orient",);     	#print error message otherwise#

	print "DONE\n";
	############# 
	#export it###


}


{
my $removals_count=0;

my @ocolor =qw(#ffffd4fdd4fd #beb700000000 #ffff00000000 #ffff63d60000 #ffff8f1846e9 #f4bb0000ffff #fffffa5d0000 #c941cb01c10e #764547e80ec8 #12f1ffff0000 #d2cef78ca2e3 #b0a39aa7855c #e2d0ff1effff #778cfbe7ffff #000000000000 #9374d1c1ffff #8f5bab84ffff #00000000ffff #be518cf0beb7 #d2f00000ffff #22d000007fff #d915d915d915 #b374b374b374 #70a370a370a3 #dfffb9488869);
#my @ocolor=qw(red yellow green pink blue purple orange #cd0d32083208 cyan #86e4ce2ceb01 #0068cd0d0000 #cd0daa2f7d15
#						#a207cd0d5a04 #9b1ecd0d9b05 brown #9fb7b6dccdd2 honeydew tan);
sub decontaminate_high_copy_repeat {
	my $seqname=shift;
	my $begin=shift;
	my $end=shift;
	my $flank=500;
	my $max_bpalign=$end-$begin+2*$flank +1;
	my $min_fraction_pair=0.7;
	my $min_fraction_region=0;
	my $join_distance=300;
	
	my $modify=shift;# modify=modify
	
	
	if ($begin>$end) {
		($end,$begin)=($begin,$end);
	}
	#flip the initial search around so that it will be in relation to this sequence region
	$modify=1;
	if ($modify) {
		my ($b,$e, $fl, $ml) = ($begin,$end, $flank, $max_bpalign);
		my $mwx = new MainWindow;
		$mwx->title("DELETE HIGH COPY REPEAT");
		$mwx->setPalette('lightgrey');
		$mwx->configure(-background=>'darkgrey');
		my $f = $mwx->Frame(	-borderwidth => 4,-relief => 'raised'  )->pack(-expand=>1,-fill=>'both');
		$f->Label(-text => "This will irrevocably delete the pairwise from parasight table!\nIt currently errors on the side of removing too much sequence\nand thus pairwise segmental dups may be lost\nThus, this is a tool for initially identifying unknown high copy repeats for better masking!")->pack(-side=> 'left',-anchor => 'e');
		$f = $mwx->Frame(	-borderwidth => 4,	-relief => 'raised',	  )->pack(-expand=>1,-fill=>'both');
		$f->Label(-text => "$seqname")->pack(-side=> 'left',-anchor => 'e');
		$f->Label(-text => " begin")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$b ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');
		$f->Label(-text => " end")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$e ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');
		$f->Label(-text => " flank")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$fl ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');
		$f->Label(-text => " max bp align")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$ml ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');
		$f = $mwx->Frame(	-borderwidth => 4,	-relief => 'raised'	  )->pack(-expand=>1,-fill=>'both');
		$f->Label(-text => " min fraction of region")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$min_fraction_region ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');
		$f->Label(-text => " min fraction of pairwise piece")->pack(-side=> 'left',-anchor => 'e');
		$f->Entry(-textvariable => \$min_fraction_pair ,-width=> 10)->pack(-side=> 'left', -anchor=> 'e');

		$f->Button(-borderwidth => 3,-activebackground=>'white',-background=>'red',
				-command=>sub { $mwx->destroy;},-text => "Remove")->pack(-side=>'right');
		$f->Button(-borderwidth => 3,-activebackground=>'white',-background=>'#bbe8ff',
				-command=>sub { $b=0; $e=0; $fl=0; $mwx->destroy;},-text => "Cancel")->pack(-side=>'right');
		$mwx->waitWindow;
		$begin = $b;
		$end = $e;
		$flank=$fl;
		$max_bpalign=$ml;

	}
	return if $begin ==0;
	$removals_count++;
	my $fremoval_count=substr ("00000".$removals_count, -4);
	print "REMOVAL#$fremoval_count\nCOLLECTING ALL PAIRWISE FOR $seqname ($begin-$end) ($flank) ($max_bpalign)\n" if !$opt{'quiet'};
	#set up regions array with initial region in the + orientation
	my @regions= ( {'seq'=>$seqname,'b'=>$begin,'e'=>$end, 'orient'=>+1,'first'=>1 } );

	#####################################
	##trace out all of the connections###
	##mark all for deletion##############
	my %deleted=();	
	my $one_pass_only=0;
	my @processed_regions=();
	my $first_pass=1;
	REGIONPASS: while (@regions > 0) {
		my @new_regions=();
		my $mlen=@m;
		foreach my $r (@regions) {
			my $begin=$$r{'b'}-$flank;
			$begin=1 if $begin < 1;
			my $end= $$r{'e'}+$flank;
			print "EXAMINING FOR PAIRS MATCHING: $$r{'seq'} $$r{'b'}-$$r{'e'}\n";
			for (my $i=0; $i< $mlen; $i++) { 
				my $refi=$m[$i];
				next if exists $pairwise2delete{$i};
				my @c=();
				next if $$r{'seq'} ne $$refi[0] && $$r{'seq'} ne $$refi[4];
			#	print "SEQUENCE MATCH $$r{'seq'} equals $$refi[0] and/or $$refi[4]\n";
				###adjust c so that seq1 is region and seq2 is match in positive orientation
				###orientation is in integer for so multiplication can take place.
			   if ( $$r{'seq'} eq $$refi[0] ) {
			   	next if $end < $$refi[1] || $begin > $$refi[2]; 
			   	#we have overlap of seq1#
			   	@c=@$refi;
			   	if ($c[5]>$c[6]) {
			   		($c[5],$c[6])=($c[6],$c[5]);
			   		$c[8]=-1;
			   	} else {
			   		$c[8]=+1;
			   	}
			   }
			   if ($$r{'seq'} eq $$refi[4] ) {
			   	next if ($end < $$refi[5] && $end < $$refi[6]) || ($begin > $$refi[5] && $begin > $$refi[6] );
			   	##we have ovlerlapof seq2#
			   	@c=@$refi;
			   	($c[0],$c[1],$c[2],$c[3],$c[4],$c[5],$c[6],$c[7])=($c[4],$c[5],$c[6],$c[7],$c[0],$c[1],$c[2],$c[3]);
			   	if ($c[1]>$c[2]) {
			   		($c[1],$c[2])=($c[2],$c[1]);
			   		$c[8]=-1;
			   	} else {
			   		$c[8]=+1;
			   	}
			   }
			   next if @c==0;
			   ###################
				###is within my region###
			#	print "  region contained $begin-$end contains: $c[0] $c[1]-$c[2]  $c[8]\n";
				my $span=$c[2]-$c[1]+1;
				next if $span > $max_bpalign;
				#print "   not to big\n";
				my $fraction_pair=0;
				my $fraction_region=0;
				if ($c[1] <= $begin && $c[2] >= $end ) {
					#print "   region surround by pairwise\n";
					#print "R             $begin - $end \n";
					#print "P    $c[1]                             $c[2]\n";
					$fraction_pair=($end-$begin+1)/$span;
					$fraction_region=1;
				} elsif ($begin <= $c[1] && $end >= $c[2]) {   
					#print "   pairwise completely contained\n";
					#print "R   $begin                            $end\n";
					#print "P           $c[1]  $c[2]            \n";
					$fraction_pair=1;
					$fraction_region=($c[2]-$c[1]+1)/($end-$begin+1);
				} elsif ( $begin <= $c[2] && $c[2] <= $end ) {
					#print "    pairwise stop contained--start outside\n";
					#print "R             $begin             $end\n";
					#print "P  $c[1]                $c[2]\n";
					my $overlap=$c[2]-$begin+1;
					$fraction_pair=$overlap/$span;
					$fraction_region=$overlap/($end-$begin+1);
				} elsif ( $begin <= $c[1] && $c[1] <= $end ) {
					#print "    pairwise start contained--top outside\n";
					#print "R $begin           $end\n";
					#print "P         $c[1]           $c[2]\n";
					my $overlap=$end-$c[1]+1;
					$fraction_pair=$overlap/$span;
					$fraction_region=$overlap/($end-$begin+1);
					
				} else {
			  		print "($c[0],$c[1],$c[2],$c[3],$c[4],$c[5],$c[6],$c[7])\n";
					&warnNpause("Your algorithm is screwy! I shouldn't be here!\n");
				}
				#print "Rfrac:$fraction_region  Pfrac:$fraction_pair\n";
				next if $fraction_pair < $min_fraction_pair;
				next if $fraction_region < $min_fraction_region; #this should stay low except for specialized purposes
				my $relative_orient=$$r{'orient'} * $c[8];
				#print "Orient relative to first is $relative_orient\n";
				print "Delete pairwise $i\n";
				$deleted{$i}=$refi;
				$pairwise2delete{$i}=$refi;
				push @new_regions, {'seq'=>$c[4],'b'=>$c[5],'e'=>$c[6],  'orient'=>$relative_orient ,"first"=>$first_pass};
				#recolor Sa Sb and M#
				$$refi[$mh{'color'}]='black';
				my $owidth=3;
				$owidth=20 if $one_pass_only==1;
				$canvas->itemconfigure("Sa$i", -fill=>'black',-outline=>$ocolor[0],-width=>$owidth );		
				$canvas->itemconfigure("Sb$i", -fill=>'black' ,-outline=>$ocolor[0],-width=>$owidth );		
				$canvas->itemconfigure("Ma$i",  -fill=>'black',-outline=>$ocolor[0],-width=>$owidth );		
					
			} #@m loop
			push @processed_regions, $r;
			if ($one_pass_only) {
				push @processed_regions, @new_regions;
				last REGIONPASS;
			}
		} #region #loop
		$first_pass=0;
		@regions=@new_regions;#
		#should have same orientations to merge
		#merge any close together regions for subsequent searches
		
		#my $pause=<STDIN>;
	} #while region looop
	print "REMOVING (", scalar keys %deleted, ") PAIRWISE\n";
	my $fname="repeat$fremoval_count.$seqname.b$begin.e$end.f$flank.m$max_bpalign.pf$min_fraction_pair.rf$min_fraction_region";
	#my $fname="repeat$fremoval_count";
	print "$fname\n";
	mkdir "decontamination" if !-d 'decontamination';
	open (OUT, ">decontamination/$fname.pairs") || die "Can't create decontamination/$fname.pairs!\n";
	foreach (reverse sort {$a <=> $b} keys %deleted) {
		print OUT join("\t", @{$m[$_]}), "\n";
	}
	close OUT;
	open (OUT, ">decontamination/$fname.coordinates") || die "Can't create decontamination/$fname.coordinates!\n";
	print OUT "seq\tbegin\tend\tname\torient\n";
	###### catalog ######
	@processed_regions= sort {$$a{'seq'} cmp $$b{'seq'} || $$a{'b'} <=> $$b{'b'} || $$a{'e'} <=> $$b{'e'}  } @processed_regions ;
	for (my $i=0; $i < @processed_regions-1; $i++ ) {
		my $ci=$processed_regions[$i];
		my $cj=$processed_regions[$i+1];
		next if $$ci{'seq'} ne  $$cj{'seq'};
		if ($$ci{'e'} + $join_distance > $$cj{'b'} ) {
			$$ci{'e'} = $$cj{'e'} if $$cj{'e'} > $$ci{'e'};
			splice(@processed_regions,$i+1,1);
			$i--;
			next;	
		}	
	}
	
	foreach (@processed_regions) {
		my $b=$$_{'b'};
		my $e=$$_{'e'};
		($b,$e)=($e,$b) if ($$_{'orient'}<0) ;
		print OUT "$$_{'seq'}\t$b\t$e\t$$_{'seq'}.$b.$e\t$$_{'orient'}\t$$_{'first'}\n";
	}
	close OUT;
	push @ocolor, (shift @ocolor);
}


}



########################################################
######Zoom in #######################################

	$canvas->Tk::bind('<Control-Button-1>',
				[sub{	my ($canv, $x, $y)=@_; $scale*=2;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						$canvas->scale("all",$cx,$cy,2,2);
						my @box=$canvas->bbox("all");
						#print "$box[0] $box[1] $box[2] $box[3]\n";
						$canvas->configure(-scrollregion=>\@box);
						$canvas->configure(-height=>$box[3], -width=>$box[2]);						
						#print $canvas->cget(-height),"  ",$canvas->cget(-width),"\n";
					}
			,Ev('x'),Ev('y')]);
	#######Zoom out ####
	$canvas->Tk::bind('<Control-Button-3>',
				[sub {	my ($canv, $x, $y)=@_; $scale*=0.5;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						$canvas->scale("all",$cx,$cy,0.5,0.5);
						my @box=$canvas->bbox("all");
						#print "$box[0] $box[1] $box[2] $box[3]\n";
						$canvas->configure(-scrollregion=>\@box);
						$canvas->configure(-height=>$box[3], -width=>$box[2]);
						#print $canvas->cget(-height),"  ",$canvas->cget(-width),"\n";
					}
			,Ev('x'),Ev('y')]);


#############################################
####################3 move any object #############
	$canvas->Tk::bind('<Shift-ButtonPress-1>', 
				[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						$iinfo{'lastX'}=$cx;  
						$iinfo{'lastY'}=$cy;
						#print "START move\n";
						
					}
			,Ev('x'),Ev('y')]);	
			
	$canvas->Tk::bind('<Shift-B1-Motion>', 
				[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						$canvas->move( 'current', $cx - $iinfo{'lastX'}, $cy-$iinfo{'lastY'} );
						$iinfo{'lastX'}=$cx;  
						$iinfo{'lastY'}=$cy;
					}
			,Ev('x'),Ev('y')]);	

############################################################
############# any delete a tag ##################################

	$canvas->Tk::bind('<Alt-Double-Button-1>',
			[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						#my @close=$canv->find("closest",$cx,$cy) ;
						my @tags = $canv->gettags( "current");
						my $id=$tags[0];
						foreach my $tag (@tags) {
							next if $tag !~/^[MSE]/;
							$id=$tag;
						}
						my ($type, $numb) = $id=~/([A-Za-z]+)(\d+)/;
						if ($type =~ /^[MS]/) {
								$canvas->delete("M$numb");
								$canvas->delete("Sa$numb");
								$canvas->delete("Sb$numb");
								$m[$numb][$mh{'hide'}]=1;
							
						
						} else {
							$canvas->delete('current');
						}
				}
		,Ev('x'),Ev('y')]);	

	#####################################################	
	########################################3
	########## QUICK COLOR OBJECTS###############
		
	$canvas->Tk::bind('<Shift-Button-3>',
			[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						my @close=$canv->find("closest",$cx,$cy) ;
						my @tags = $canv->gettags( "$close[0]");
						my $id=$close[0];
						foreach my $tag (@tags) {
							next if $tag !~/^[MSE]/;
							$id=$tag;
						}
						my ($type, $numb) = $id=~/([A-Za-z]+)(\d+)/;
						my $text.="Nothing of interest\n";
						my $color=$opt{'quick_color'};
						if ($type =~ /M/) {
							$m[$numb][$mh{'color'}]=$color;
							$canvas->itemconfigure("M$numb", -fill=>$color );
							$oldhighlightcolor=$color if $oldhighlightid eq "M$numb" ;
						} elsif ($type=~/S/) {
							$oldhighlightcolor=$color if $oldhighlightid =~ /^S[ab]$numb/ ;
							$m[$numb][$mh{'scolor'}]=$color;
							$canvas->itemconfigure("Sa$numb", -fill=>$color);
							$canvas->itemconfigure("Sb$numb", -fill=>$color);
						} elsif ($type=~/E/) {
							$oldhighlightcolor=$color if $oldhighlightid eq "E$numb" ;
							$e[$numb][$eh{'color'}]=$color;
							$canvas->itemconfigure("E$numb", -fill=>$color);
							#$oldhighlightcolor=$color if $oldhighlightid eq "E$numb" ;

						} else {
							$canvas->itemconfigure("$id",-fill,$color);
						}
				}
		,Ev('x'),Ev('y')]);	


	$canvas->Tk::bind('<Shift-Double-Button-3>',
			[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						my @close=$canv->find("closest",$cx,$cy) ;
						my @tags = $canv->gettags( "$close[0]");
						my $id=$close[0];
						foreach my $tag (@tags) {
							next if $tag !~/^[MSE]/;
							$id=$tag;
						}
						my ($type, $numb) = $id=~/([A-Za-z]+)(\d+)/;
						if ($type =~ /M/) {
							my $color='black';  #  no way to know if inter or intra
							$oldhighlightcolor=$color if $oldhighlightid eq "M$numb" ;
							$m[$numb][$mh{'color'}]='';
							$canvas->itemconfigure("M$numb", -fill=>$color );
						} elsif ($type=~/S/) {

							my $color=$opt{'sub_color'};
							$color=$accsub{$m[$numb][4]}{'color'} if defined $accsub{$m[$numb][4]}{'color'};
							$m[$numb][$mh{'scolor'}]='';
							$oldhighlightcolor=$color if $oldhighlightid =~ /^S[ab]$numb/ ;
							$canvas->itemconfigure("Sa$numb", -fill=>$color );		
							$color=$accsub{$m[$numb][0]}{'color'} if defined $accsub{$m[$numb][4]}{'color'};
							$canvas->itemconfigure("Sb$numb", -fill=>$color);


						} elsif ($type=~/E/) {
								my $color=$opt{'extra_color'};  #  no way to know if inter or intra
								$oldhighlightcolor=$color if $oldhighlightid eq "E$numb" ;
								$e[$numb][$eh{'color'}]='';
								$canvas->itemconfigure("E$numb", -fill=>$color );

						} else {
							$canvas->itemconfigure ("$id", -fill => black);
						}
				}
		,Ev('x'),Ev('y')]);	

#############################################################
########## raise and lower tags ############################
	$canvas->Tk::bind('<Alt-Button-1>',
						[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						$canv->lower('current');
						}
			,Ev('x'),Ev('y')]);
		
	$canvas->Tk::bind('<Alt-Button-3>',
						[sub{
						my ($canv, $x, $y)=@_;
						my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
						$canv->raise('current');
					}
			,Ev('x'),Ev('y')]);
			

##############################################
#############################################
################ user defined ###############

	$canvas->Tk::bind('<Shift-Control-Button-1>',
		[sub{
						
			my ($canv, $x, $y)=@_;
			my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
			my @tags = $canv->gettags( "current");
			my $id=$tags[0];
			foreach my $tag (@tags) {
				next if $tag !~/^[MSE]/;
				$id=$tag;
			}
			&execute_execute($id,'');
					}
			,Ev('x'),Ev('y')
		]);
	$canvas->Tk::bind('<Shift-Control-Button-3>',
		[sub{
						
			my ($canv, $x, $y)=@_;
			my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
			my @tags = $canv->gettags( "current");
			my $id=$tags[0];
			foreach my $tag (@tags) {
				next if $tag !~/^[MSE]/;
				$id=$tag;
			}
			&execute_execute($id,'2');
					}
			,Ev('x'),Ev('y')
		]);
	$canvas->Tk::bind('<Double-Button-1>',
		[sub{
						
			my ($canv, $x, $y)=@_;
			my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
			my @tags = $canv->gettags( "current");
			my $id=$tags[0];
			foreach my $tag (@tags) {
				next if $tag !~/^[MSE]/;
				$id=$tag;
			}
			&execute_execute($id,'');
					}
			,Ev('x'),Ev('y')
		]);
	$canvas->Tk::bind('<Double-Button-2>',
		[sub{
						
			my ($canv, $x, $y)=@_;
			my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
			my @tags = $canv->gettags( "current");
			my $id=$tags[0];
			foreach my $tag (@tags) {
				next if $tag !~/^[MSE]/;
				$id=$tag;
			}
			&execute_execute($id,'2');
					}
			,Ev('x'),Ev('y')
		]);	
	$canvas->Tk::bind('<Shift-Control-Button-2>',
		[sub{
						
			my ($canv, $x, $y)=@_;
			my ($cx,$cy)=($canv->canvasx($x),$canv->canvasy($y) );
			my @tags = $canv->gettags( "current");
			my $id=$tags[0];
			foreach my $tag (@tags) {
				next if $tag !~/^[MSE]/;
				$id=$tag;
			}
			&execute_execute($id,'3');
					}
			,Ev('x'),Ev('y')
		]);	
	
sub execute_execute {
	#this subroutine executes a predefined command substituting 
	#data from a column via replacing occurence of {}
	my $type=shift;
	my $n=shift;
	my $command=$opt{"execute$n"};
	my $desc=$opt{"execute$n"."_desc"};
	my $array=$opt{"execute$n"."_array"};
	my ($numb) = $type =~ /(\d+)/;
	
	print "row$n  pairORextra$type  arrayname$array\n";
	next if $type =~ /E/ && $array ne 'e';
	next if $type !~ /E/ && $array ne 'm';
	my @c=@{$$array[$numb]};
	print "ARRAY", join ("==>",@c), "\n";
	###evaluate and execute command while trapping runtime errors###
	print "PRE ($command)\n" if !$opt{'quiet'};
	$command = '$command = " ' . $command . ' "';
	print "PREEVAL ($command)\n" if !$opt{'quiet'};
	eval ($command);
	if ($@) {
		warn $@ ;
	} else {
		print "SYSTEM CALL($command)\n" if !$opt{'quiet'};
		system "$command";
	}
}

	################################################
	########### FIRST DRAW #########################
	################################################
	$first_pass=1;
	&reshowNredraw;
	$mw->deiconify();
	$mw->raise();
	$mw->waitVisibility();
	
	#################################################
	######### EXECUTE SUPPLIED PRECODE ##############
	my $to_die=$opt{'die'};
	$mw->after(50, sub {
		if ($opt{'precode'}) {
			print "PRECODE\n" if !$opt{'quiet'};
			if ($opt{'precode'} =~ /^file:(\S+)/ ) {
				my $prefile=$1;
				print "PRECODE IS FILE.  LOADING NOW ...\n" ;
				open (PRECODE, $prefile) || die "Can't read precode file ($prefile)!";
				$opt{'precode'}='';
				while (<PRECODE>) {
					s/\r\n/\n/;
					$opt{'precode'}.=$_;
				}
			}
			my $precodetext=$opt{'precode'};
			###scan precode for valid $opt{'variables'}###
			while ( $precodetext =~ /opt\{["']([a-z_]+)["']\}/mg ) {
				#my $value=$1;
				#print "$value\n";
				warn "\$opt{\"$1\"} is not a valid option!\n" if ! defined $opt{$1};
				#my $pause=<STDIN>;
			
			}
			###run the precode and watch for errors###
			eval $precodetext; 
			if ($@) {
				warn "PRECODE REPORTED A SYNTAX ERROR!\n$@\n" ;
				my @precode_array =split /\n/ , $precodetext;
				if ( $@ =~ /line (\d+)/ ) {
					warn "LINE $1 WAS: $precode_array[$1]\n";
				}
			}
			
		}

		if ($to_die ) {
			print "\n-die has been envoked--quitting now\n" if !$opt{'quiet'};;
			exit;
		}
	} );
	################################################
	########## MAIN LOOP ###########################
 	MainLoop;



#################################################################################33
####################### SUBROUTINES #############################################
#################################################################################

sub find_column_options {
	#######################
	####figures out name of column given a number and vis-versa
	#####m only #######
	my @farray=qw(sub_scale_col sub_scale_col2 sub_labelseq_col sub_labelseq_col2 
	sub_labelseqe_col sub_labelseqe_col2 alignment_col alignment_col2
	sub_labelhit_col sub_labelhit_col2 filter1_col  filter2_col 
	 filterpre1_col 
	filterpre2_col pair_type_col colorsub_hitcond_col pair_type_col2);
   for my $n (@farray) {
		if ($opt{$n} =~/[a-zA-z]/) {
			if ( defined $mh{$opt{$n}} ) {
				$opt{$n}=$mh{$opt{$n}};
			} else {
   			$opt{$n}.='??' if $opt{$n} !~/\?/;
			}
		}
		$colheader{("$n"."_header")}='( )';
		$colheader{("$n"."_header")}="($mheader[$opt{$n}])" if  $opt{$n} =~/^\d+$/;
	}
	##################################
	#### e only #####################
	@farray=qw(filterextra1_col  filterextra2_col);
   for my $n (@farray) {
		if ($opt{$n} =~/[a-zA-z]/) {  
			if ( defined $eh{$opt{$n}} ) {
				$opt{$n}=$eh{$opt{$n}};
			} else {
   			$opt{$n}.='??' if $opt{$n} !~/\?/;
			}
		}
		$colheader{("$n"."_header")}='( )';
		$colheader{("$n"."_header")}="($eheader[$opt{$n}])" if  $opt{$n} =~/^\d+$/;
	}
	#####################################
	#### m and e ########################
	my %fhash = ( mark_col=> mark_array, mark_col2 => mark_array
			) ;
	
   for my $n (keys %fhash) {
   	if ($opt{$fhash{$n}}=~/[psm]/ ) {
   		###change text to number###
   		if ($opt{$n} =~/[a-zA-z]/) {
				if ( defined $mh{ $opt{$n} }) {
					$opt{$n}=$mh{$opt{$n}};
				} else {
					$opt{$n}.='??' if $opt{$n} !~/\?/;
				}
			}
			$colheader{("$n"."_header")}='( )';
			$colheader{("$n"."_header")}="($mheader[$opt{$n}])" if  $opt{$n} =~/^\d+$/;
   	} elsif ($opt{$fhash{$n}}=~/e/ ) {
   		###extras###
   		if ($opt{$n} =~/[a-zA-z]/) {
				if ( defined $eh{$opt{$n}}) {
					$opt{$n}=$eh{$opt{$n}};
				}else {
					$opt{$n}.='??' if $opt{$n} !~/\?/;
				}
			}
			$colheader{("$n"."_header")}='( )';
			$colheader{("$n"."_header")}="($eheader[$opt{$n}])" if  $opt{$n} =~/^\d+$/;
   	} else {
			$colheader{("$n"."_header")}='( )';
     	}
	}
}


sub mark_window {
	#,'#bbe8ff' light blue
	# '#ffff9b' light yellow
	#$ballooni->attach($tmp, -msg => &balloon_format_var('sub_arrow_on') ) if $opt{'help_on'};
#some code from gbarr#
	my $mw = new MainWindow;
	$ballooni = $mw->Balloon(-initwait =>600, -background => '#ffff9d' ,-font=>'Courier 8');
	$mw->title("MARK / FIND OBJECTS");
	my $frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-fill=>'x');
	my ($tl,$te); #tmp varaibles to toss
	($tl,$te) = &fast_lentry($frame,$tmp,"pattern", 'mark_pattern', 50 , \&mark);	
	$te->pack(-fill=>'x',-expand=>1);
	#$te->configure(-background=>'lightgray');
	$frame->Button(-borderwidth => 4,-command=>[\&mark,'markall'],-text => "Mark All")->pack(-side=>'right');

	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	$frame->Button(-borderwidth => 4,-command=>[\&mark, 'findnext'],-text => "Find Next")->pack(-side=>'right');
	&fast_lentry($frame,$tmp,"col", 'mark_col', 5 , \&find_column_options);	
	$frame->Label(-textvariable =>\$colheader{'mark_col_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"col2", 'mark_col2', 5 , \&find_column_options);	
	$frame->Label(-textvariable =>\$colheader{'mark_col2_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"color", 'mark_color', 15 , \&mark);	
	
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-fill=>'x');	
	$tmp=$frame->Checkbutton(-text=>'Permanent Color',-variable => \$opt{'mark_permanent'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('mark_permanent') ) if $opt{'help_on'};
	$frame->Button(-borderwidth => 4,-command=>\&mark_remove,-text => "Clear")->pack(-side=>'right');
	my $text_variable='pairs';

	my $opt_menu=$frame->Optionmenu(-textvariable=>\$text_variable,
			-variable=>\$opt{'mark_array'}, -options => [['pairs','m'], ['extras','e']],
			)->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($opt_menu, -msg => &balloon_format_var('mark_array') ) if $opt{'help_on'};
	my $lab=$frame->Label(-text=>'Mark/Find:')->pack(-side=> 'left',-anchor => 'e');	
	my $mark_pairs=$frame->Checkbutton(-text=>'pairs',-variable => \$opt{'mark_pairs'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($opt_menu, -msg => &balloon_format_var('mark_pairs') ) if $opt{'help_on'};
	my $mark_subs=$frame->Checkbutton(-text=>'subs',-variable => \$opt{'mark_subs'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($opt_menu, -msg => &balloon_format_var('mark_subs') ) if $opt{'help_on'};
	$opt_menu->configure(-command=> 
			sub {
				if ($opt{'mark_array'} eq 'e') {
					#$lab->configure(-background=>'darkgray');
					$mark_pairs->configure(-state=>'disabled');
					$mark_subs->configure(-state=>'disabled');
				} else {
					#$lab->configure(-background=>'lightgray');
					$mark_pairs->configure(-state=>'normal');
					$mark_subs->configure(-state=>'normal');
				}
			}
		);
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-fill=>'x');
	($tl,$te)=&fast_lentry($frame,$tmp,"advanced", 'mark_advanced', 50 , \&mark,);	
	$te->pack(-fill=>'x',-expand=>1);

}


sub mark_windowx {
	#,'#bbe8ff' light blue
	# '#ffff9b' light yellow
	#some code from gbarr#
	my $mw = new MainWindow;
	$mw->title("MARK / FIND OBJECTS");
	my $frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-fill=>'x');
	my ($tl,$te); #tmp varaibles to toss
	($tl,$te) = &fast_lentry($frame,$tmp,"pattern", 'mark_pattern', 50 , \&mark);	
	$te->pack(-fill=>'x',-expand=>1);
	#$te->configure(-background=>'lightgray');
	$frame->Button(-borderwidth => 4,-command=>[\&mark,'markall'],-text => "Mark All")->pack(-side=>'right');

	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	$frame->Button(-borderwidth => 4,-command=>[\&mark, 'findnext'],-text => "Find Next")->pack(-side=>'right');
	&fast_lentry($frame,$tmp,"col", 'mark_col', 5 , \&find_column_options);	
	$frame->Label(-textvariable =>\$colheader{'mark_col_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"col", 'mark_col2', 5 , \&find_column_options);	
	$frame->Label(-textvariable =>\$colheader{'mark_col2_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"color", 'mark_color', 15 , \&mark);	
	
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-fill=>'x');	
	$tmp=$frame->Checkbutton(-text=>'Permanent Color',-variable => \$opt{'mark_permanent'})->pack(-side=>'left',-anchor=>'e');
	#$ballooni->attach($tmp, -msg => &balloon_format_var('mark_permanent') ) if $opt{'help_on'};
	$frame->Button(-borderwidth => 4,-command=>\&mark_remove,-text => "Clear")->pack(-side=>'right');
	my $text_variable='pairs';

	my $opt_menu=$frame->Optionmenu(-textvariable=>\$text_variable,
			-variable=>\$opt{'mark_array'}, -options => [['pairs','m'], ['extras','e']],
			)->pack(-side=> 'left',-anchor => 'e');
	#$ballooni->attach($opt_menu, -msg => &balloon_format_var('mark_array') ) if $opt{'help_on'};
	my $lab=$frame->Label(-text=>'Mark/Find:')->pack(-side=> 'left',-anchor => 'e');	
	my $mark_pairs=$frame->Checkbutton(-text=>'pairs',-variable => \$opt{'mark_pairs'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($mark_pairs, -msg => &balloon_format_var('mark_pairs') ) if $opt{'help_on'};

	my $mark_subs=$frame->Checkbutton(-text=>'subs',-variable => \$opt{'mark_subs'})->pack(-side=>'left',-anchor=>'e');
	$opt_menu->configure(-command=> 
			sub {
				if ($opt{'mark_array'} eq 'e') {
					#$lab->configure(-background=>'darkgray');
					$mark_pairs->configure(-state=>'disabled');
					$mark_subs->configure(-state=>'disabled');
				} else {
					#$lab->configure(-background=>'lightgray');
					$mark_pairs->configure(-state=>'normal');
					$mark_subs->configure(-state=>'normal');
				}
			}
		);
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-fill=>'x');
	($tl,$te)=&fast_lentry($frame,$tmp,"advanced", 'mark_advanced', 50 , \&mark,);	
	$te->pack(-fill=>'x',-expand=>1);

}

{
my %tmp_marks=();
my $array_position_next=0;
my $query='naNa';

sub mark {
	my $command=shift;
	#print "$command\n";
	if ( $command eq 'markall' ) {
		#reset start position
		$array_position_next=0;
		$query='naNa';
	}
	if ( $command eq 'findnext') {
		if ($query ne $opt{'mark_pattern'} ) {
			print "RESETING QUERY\n" if !$opt{'quiet'};;
			$array_position_next=0;
			$query=$opt{'mark_pattern'};
		}
	}
	my $a = eval("\\\@$opt{'mark_array'}");
	#print "$opt{'mark_array'}=>$a=>",scalar @$a," (\@m)\n";
	my $color=$opt{'mark_color'};
	my $total_matches=0;
	for (my $i=$array_position_next; $i < @$a; $i++) {
		my $numb=$i;
		my $match=0;
		$match=1 if $opt{'mark_col'} =~/^\d+$/ &&  $$a[$i][$opt{'mark_col'}]=~/$opt{'mark_pattern'}/;
		#print "$i col:$$a[$i][$opt{'mark_col'}] (",$$a[$i][$opt{'mark_col'}]=~/$opt{'mark_pattern'}/,")match\n";
		$match=1 if   $opt{'mark_col2'} =~/^\d+$/  && $$a[$i][$opt{'mark_col2'}] =~ /$opt{'mark_pattern'}/ ;
		#print "$i col:$$a[$i][$opt{'mark_col2'}] $match\n";
		if ($opt{'mark_advanced'} ne '' ) {
			my $c=$$a[$i];
			my $result=eval($opt{'mark_advanced'});
			$match=1 if $result !=0;
		}
		next if $match==0;
		
		if ($opt{'mark_array'} eq 'm') {
			if ($opt{'mark_pairs'} ==1) {
				if ($opt{'mark_permanent'}) {
					$m[$numb][$mh{'color'}]=$color;
				} else {
					$tmp_marks{"M$numb"} = $canvas->itemcget("M$numb",-fill) if !defined $tmp_marks{"M$numb"};
				}
				$canvas->itemconfigure("M$numb", -fill=>$color );
			}
			if ($opt{'mark_subs'} ==1) {
				if ($opt{'mark_permanent'}) {
					$m[$numb][$mh{'scolor'}]=$color;
				} else {
					$tmp_marks{"Sa$numb"} = $canvas->itemcget("Sa$numb",-fill)if !defined $tmp_marks{"Sa$numb"};
					$tmp_marks{"Sb$numb"} = $canvas->itemcget("Sb$numb",-fill)if !defined $tmp_marks{"Sb$numb"};
				}
				$canvas->itemconfigure("Sa$numb", -fill=>$color);
				$canvas->itemconfigure("Sb$numb", -fill=>$color);
			}
		} elsif ($opt{'mark_array'} eq 'e') {
			if ($opt{'mark_permanent'} ==1) {
				$e[$numb][$eh{'color'}]=$color;
			} else {
				$tmp_marks{"E$numb"} = $canvas->itemcget("E$numb",-fill)if !defined $tmp_marks{"E$numb"};
			}
			$canvas->itemconfigure("E$numb", -fill=>$color);
		} 
		$total_matches++;
		if ($command eq 'findnext') {
			#i found one so I am done#
			#can i center on the object#
			#print "FIND $i\n";
			$array_position_next=$i+1;
			return;
		}
		
	}
	if ($command eq 'markall') {
		print "TOTAL /$opt{'mark_pattern'}/ found was $total_matches\n" if !$opt{'quiet'};;
	}
	if ($command eq 'findnext') {
		print "End of objects without match! Postion reset to beginning!\n"  if !$newopt{'quiet'};
		$array_position_next=0;
	}
} #close sub


sub mark_remove {
	foreach my $o (keys %tmp_marks) {
		$canvas->itemconfigure($o,-fill=> $tmp_marks{$o});
	}
	%tmp_marks=();
}

}  #close private variables


sub indexcard_options {
	#,'#bbe8ff' light blue	# '#ffff9b' light yellow	
	#the index cards were based off of some code from gbarr.  thanks!#
	my $mw = new MainWindow;
	$mw->title("PARASIGHT OPTIONS");
	$ballooni = $mw->Balloon(-initwait =>600, -background => '#ffff9d' ,-font=>'Courier 8');

	my $current;
	$mw->setPalette('lightgrey');

	$mw->configure(-background=>'darkgrey');
	my $bf= $mw->Frame()->pack(-side=>'top',-anchor=>'w');
	$bf->Button(-borderwidth => 4,-background=>'#bbe8ff',-command=>\&reshowNredraw,-text => "Reshow, Rearrange & Redraw")->pack(-side=>'left');
	$bf->Button(-borderwidth => 4,-background=>'white',-command=>\&redraw,-text => "Redraw Only")->pack(-side=>'left');

	my $f = $mw->Frame(
			-borderwidth => 4,
			-relief => 'raised',
	  )->pack(-expand=>1,-fill=>'both');
   my %br;
	my @l=();
	my $depth=0;
	 my $tf = $f->Frame(
		-borderwidth => 0,
		-relief => 'flat'
	 )->pack(-side => 'top',-fill=>'x');

	#print (pop @{$tf->configure(-background)}), "\n";

	foreach ( "MAIN\n","SEQ\nPAIRS","SUB\n","EXTRA\n","GRAPH\n","FILTER","MISC\n") {
		 my $label= $tf->Label(
			-text => $_,
			-borderwidth => 0,
			-relief=>'sunken',
			-background=>'grey',
			-padx => 5, 
			-anchor =>  'w',
		 )->pack(-side => 'left');	 
		 $br{$label}=$depth;
		 $depth+=2;
		 push @l, $label;
	}
	$depth+=2;
	#print "DEPTH:$depth\n";
	foreach (@l) { $_->configure(-pady=>$depth-$br{$_}, -borderwidth=>$br{$_},-padx=>$depth-$br{$_}); }
   
	my $minimize = $tf->Label(
			-text =>    " \n ",
			-borderwidth => $depth,
			-background=>'darkgrey',

			-relief => 'sunken',
			-padx =>2, -pady => 0,
			-anchor =>  'w'
	)->pack(-side => 'right',-fill=> 'x',-expand=>1);	 

	my %c;
	foreach my $i (@l) {
		$i->bind('<1>', [ 	sub {
			#print "CUR$current ==== $i\n";
			my $i=shift;
			return if $current eq $i;
			$br{$i}=-2;
			foreach (@l) {
				$br{$_}+=2 if ($br{$_} < $depth) ;
				$_->configure(-borderwidth=>$br{$_}
							, -pady=>$depth-$br{$_},-padx=>$depth-$br{$_},
					);
				if ( $br{$_}==0 ) {
					$_->configure(-background=>'lightgrey');
				} else {
					$_->configure(-background=>'grey');
				}

			}
			$c{$current}->packForget;
			$current = $i;
			$c{$current}->pack(-side=>'bottom',-fill=>'x',-expand=>1);  
		}, $i]);
		
		$c{$i}= $f->Frame(-borderwidth => 0,-relief => 'raised');

	}

  	$current=$l[0];
  	$current->configure(-background=>'lightgrey');
	#print "COLOR",join (" ",$current->configure(-background)),"\n";
  	#print "CURRENT$current\n";
  	$c{$current}->pack(-side=>'bottom',-fill=>'x',-expand=>1);
	&card_main($c{$l[0]});
	&card_seq($c{$l[1]});
	&card_sub_subscale($c{$l[2]});;
	&card_extra($c{$l[3]});
	&card_graph($c{$l[4]});
	&card_filter($c{$l[5]});
	&card_misc($c{$l[6]});
}

sub card_main {
#	$tmp->configure(-background=>'#ffff9d');
	my $mw=shift;
	my ($frame, $tmp);
	########################################################################
	$frame = $mw->Frame(-borderwidth=>1);
	$frame->Label(-text=> "MAIN DATA INPUT",-borderwidth=>1)->pack(-side=>'left',-anchor=>'w');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');

	$frame = $mw->Frame(-borderwidth=>1);
	&fast_lentry($frame,$tmp,"Filename", 'filename', 25 , \&redraw);		
	$tmp=&doublelabel($frame,$tmp," -in", \$opt{'in'}, 20 , \&reshowNredraw,);
	$ballooni->attach($tmp, -justify => 'left',
	  -msg=>( &help_format( "-in [ parsight path] will display the last loaded or last saved parasight file.  This field can not be directly altered by the user." ))) if $opt{'help_on'};
	
	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');

	$frame = $mw->Frame(-borderwidth=>1);
	$tmp=$frame->Checkbutton(-text=>'  on',-variable => \$opt{'filename_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('filename_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"size", 'filename_size', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"color", 'filename_color', 7 , \&redraw);
	&fast_lentry($frame,$tmp,"offset v:",'filename_offset', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"h:", 'filename_offset', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"pattern",'filename_pattern',8, \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');


	
	$frame = $mw->Frame(-borderwidth=>1);
	$tmp=&fast_lentry($frame,$tmp,"-align", 'align', 65 , \&reshowNredraw,'#bbe8ff');
	$ballooni->attach($tmp, -justify => 'left',
	  -msg=>("-align [filepath1:filepath2:etc]                      \n"
	        ." *more files containing alignments may be added       \n"
			  ." *must be in miropeats format tab-delimited format    \n"
			  ."  (name1 b1 e1 seqlen1 name2 b2 e2 seqlen2 [optcolumns\n" 
			  ." (the text will disappear on sucessful load)          \n" )) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame(-borderwidth=>1);
	$tmp=&fast_lentry($frame,$tmp,"-extra", 'extra', 65 , \&reshowNredraw,'#bbe8ff');
	$ballooni->attach($tmp, -justify => 'left',
	  -msg=>("-extra [filepath1:filepath2:etc]                     \n"
	        ." *more files containing extra sequence annotation     \n"
			  ." *the first 3 columns must be (seqname begin end)     \n"
			  ." (the text will disappear on sucessful load)          \n" )) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame(-borderwidth=>1);
	$tmp=&fast_lentry($frame,$tmp,"-showseq", 'showseq', 65 , \&reshowNredraw,'#bbe8ff');
	$ballooni->attach($tmp, -justify => 'left',
	  -msg=>("Designates sequences to draw:                               \n"
	        ." colon-delimited seq name(s) (e.g. name1:name2:name3)       \n"
			  ." ALL (no-colon) will draw all sequences                     \n"
			  ." if just one sequence name needs ending-colon (e.g. name:)  \n"
			  ." no-colon  assumes that it is a file to open                \n"
	        ."  (seq lengths can be specified in 2nd column)              \n" )) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame(-borderwidth=>1);
	$frame->Label(-text=> "                 show query (1st) sequence only ")->pack(-side=>'left',-anchor=>'w');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'showqueryonly'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('shoqueryonly') ) if $opt{'help_on'};
	$frame->Label(-text=> "(works with ALL only)")->pack(-side=>'left',-anchor=>'w');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	


	$frame = $mw->Frame(-borderwidth=>1);
	$tmp=&fast_lentry($frame,$tmp,"-showsub", 'showsub', 60 , \&reshowNredraw,'#bbe8ff');
	$ballooni->attach($tmp, -justify => 'left',
	  -msg=>("Chooses the subject sequences to be drawn:           \n"
	        ." colon-delimited seq name(s) (e.g. name1:name2:name3)\n"
			  ." ALL (no-colon) will draw all sequences              \n"
			  ." for just one sequence name need colon (e.g. name:)  \n"
			  ." no-colon  assumes that it is a file to open         \n"
	        ."  (seq lengths can be specified in 2nd column)       \n" ))if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame(-borderwidth=>1);
	
	

	$frame = $mw->Frame(-borderwidth=>1);
	&fast_lentry($frame,$tmp,"SCREEN   Indent: left", 'canvas_indent_left', 5 , \&redraw);
	&fast_lentry($frame,$tmp,"right", 'canvas_indent_right', 5 , \&redraw);
	&fast_lentry($frame,$tmp,"top", 'canvas_indent_top', 5 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame = $mw->Frame(-borderwidth=>1);
	&fast_lentry($frame,$tmp,"Window Pixel Width", 'window_width', 6 , \&redraw);
	&fast_lentry($frame,$tmp,"Screen bp Width:", 'canvas_bpwidth', 10 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');


	$frame = $mw->Frame(-relief => 'groove', -bd => 6 );	
	$frame->Label(-text=> "ALIGN ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$mstring,-wraplength=>450,-justify=>'left')->pack(-side=> 'top',-anchor => 'w');
	$frame->pack(-side => 'top', -anchor => 'w');
	$frame->Label(-text=> "EXTRA ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$estring, -wraplength=>450,-justify=>'left' )->pack(-side=> 'top',-anchor => 'w');
	$frame->pack(-side => 'top', -anchor => 'w');
}

sub card_seq { 	
	my $mw=shift;
	my ($frame, $tmp);
	########################
	#############################
	######frame for sequence options#########
	$frame = $mw->Frame();	
	$tmp=$frame->Label(-text => "SEQUENCE", -background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("Sequence options that need explainations:             \n"
	        ."  *spacing is somewhat confusing a line of sequence   \n"
	        ."   can wrap to form the equivalent of a paragraph     \n"
	        ."   spacing between lines (between paragraphs) and     \n"
	        ."   line wrapping spacing can be set to different      \n"
	        ."   values                                             \n"
	        ."                                                      \n"
	        ."                                                      \n"
			  ."                                                      \n")) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame();	
	
	$tmp=$frame->Label(-text => " -arrangeseq",-background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("Arranges sequences to draw:                                              \n"
	  		  ."  onerperline = each sequence is drawn on a separate line (paragraph)    \n"
	  		  ."  sameline    = sequences are drawn on same line with set spacing between\n"
	  		  ."  file        = exact line and begin positions are designated in a file  \n")) if $opt{'help_on'};

	$tmp=$frame->Optionmenu(-background=> '#bbe8ff',-textvariable=>\$opt{'arrangeseq'}, -options => ['oneperline','sameline',
	'file'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('arrangeseq') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"arrange file", 'arrange_file', 30 , \&redraw,'#bbe8ff');
	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
#	$frame = $mw->Frame();	
	
#	$tmp=$frame->Label(-text => " -color",-background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
#	$ballooni->attach($tmp,-justify => 'left',
#	  -msg=>("Choose color schemes for sequence and pairwise:    \n"
#			  ."--not yet implemented!!!!                          \n"));
#	$frame->Optionmenu(-background=> '#bbe8ff',-textvariable=>\$opt{'color'},
#	         -options => ['NONE','???'] )->pack(-side => 'left',-anchor=>'e');
#	$ballooni->attach($tmp, -msg => &balloon_format_var('mark_pairs') ) if $opt{'help_on'};

	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	



	$frame = $mw->Frame();		
	&fast_lentry($frame,$tmp," Color", 'seq_color', 10 , \&redraw,'#bbe8ff');
	&fast_lentry($frame,$tmp,"Width", 'seq_width', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');		
	
	$frame = $mw->Frame();	
	$frame->Label(-text => " Names")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'seq_label_on'})
		->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'seq_label_color', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"size", 'seq_label_fontsize', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset v:", 'seq_label_offset', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"h:", 'seq_label_offset_h', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"pattern", 'seq_label_pattern', 10 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');		

	$frame = $mw->Frame();
	&fast_lentry($frame,$tmp," Spacing: btwn lines (pixels)", 'seq_line_spacing_btwn', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"wrap within line (pixels)", 'seq_line_spacing_wrap', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"btwn sequences (bp)", 'seq_spacing_btwn_sequences', 8 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	####frame ticks #######
	$frame = $mw->Frame();
	$tmp=$frame->Label(-text => "SEQUENCE TICK MARKS")->pack(-side=> 'left',-anchor => 'e');
	$tmp->configure(-background=>'#ffff9d');
	$ballooni->attach($tmp, -justify => 'left',
	  -msg=>("TICK MARKS   (self explainatory)             \n"
	        ." interval tick are equally spaced            \n"
			  ." begin ticks appear at beginning of sequence \n"
			  ." end ticks appear at end of sequence         \n" )) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$tmp=$frame->Checkbutton(-text=>'whole/continuous line numbering',-variable => \$opt{'seq_tick_whole'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_whole') ) if $opt{'help_on'};

	$frame = $mw->Frame();
	$frame->Label(-text => " Interval Tick:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'seq_tick_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"width", 'seq_tick_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"length", 'seq_tick_length', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"color", 'seq_tick_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'seq_tick_offset', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"bp inteval", 'seq_tick_bp', 8 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => "      label:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'seq_tick_label_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"size", 'seq_tick_label_fontsize', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"color", 'seq_tick_label_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"offset",'seq_tick_label_offset', 3 , \&redraw);	
	$frame->Label(-text => "anchor") -> pack(-side=>'left',-anchor=>'e');
	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'seq_tick_label_anchor'}, -options => ['n','e','w','s','ne','nw','se','sw'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_label_anchor') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"scaling",'seq_tick_label_multiplier', 8 , \&redraw);	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => " Begin Tick:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'seq_tick_b_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_b_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"width", 'seq_tick_b_width', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"length", 'seq_tick_b_length', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"color", 'seq_tick_b_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'seq_tick_b_offset', 4 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => "      label:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'seq_tick_b_label_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_b_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"size", 'seq_tick_b_label_fontsize', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"color", 'seq_tick_b_label_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"offset v:",'seq_tick_b_offset', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"h:",'seq_tick_b_label_offset_h', 3 , \&redraw);	
	$frame->Label(-text => "anchor") -> pack(-side=>'left',-anchor=>'e');
	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'seq_tick_b_label_anchor'}, -options => ['n','e','w','s','ne','nw','se','sw'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_b_label_anchor') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"scaling",'seq_tick_b_label_multiplier', 8 , \&redraw);	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => " End Tick:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'seq_tick_e_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_e_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"width", 'seq_tick_e_width', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"length", 'seq_tick_e_length', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'seq_tick_e_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => "      label:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'seq_tick_e_label_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_e_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"size", 'seq_tick_e_label_fontsize', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"color", 'seq_tick_e_label_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"offset v:",'seq_tick_e_label_offset', 3 , \&redraw);	
	&fast_lentry($frame,$tmp,"h:",'seq_tick_e_label_offset_h', 3 , \&redraw);	
	$frame->Label(-text => "anchor") -> pack(-side=>'left',-anchor=>'e');
	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'seq_tick_e_label_anchor'}, -options => ['n','e','w','s','ne','nw','se','sw'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('seq_tick_e_label_anchor') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"scaling",'seq_tick_e_label_multiplier', 8 , \&redraw);	

	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	

	######frame for pair options#########
	$frame = $mw->Frame();
	$tmp=$frame->Label(-text => "PAIRS")->pack(-side=> 'left',-anchor => 'e');
	$tmp->configure(-background=>'#ffff9d');
	$ballooni->attach($tmp, -justify => 'left',
	  -msg=>("PAIR(WISE)                                       \n"  
	        ." *pairwise determination default is intra-screen \n"
			  ."  else define equality using 2 columns here      \n"
			  ." *draw level determines who is on top in picture \n"
			  )) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');		
	$frame = $mw->Frame();	
	&fast_lentry($frame,$tmp,"DEFINING col", 'pair_type_col', 8 ,\&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'pair_type_col_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"pattern", 'pair_type_col_pattern', 10 , \&redraw,);
	&fast_lentry($frame,$tmp,"col2", 'pair_type_col2', 8 , \&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'pair_type_col2_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"pattern", 'pair_type_col2_pattern', 10 , \&redraw,);
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	
	$frame = $mw->Frame();
	$tmp=$frame->Label(-text => "DRAW LEVEL:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'pair_level'},
	         -options => ['NONE','inter_over_intra','intra_over_inter'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('pair_level') ) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');		
	
	$frame = $mw->Frame();
	$frame->Label(-text => "INTRA-PAIRS:  ")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'pair_intra_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('pair_intra_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'pair_intra_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'pair_intra_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'pair_intra_offset', 3 , \&redraw);
	$tmp=$frame->Checkbutton(-text=>'lines on',-variable => \$opt{'pair_intra_line_on'})->pack(-side=>'left',-anchor=>'e');	
	$ballooni->attach($tmp, -msg => &balloon_format_var('pair_intra_line_on') ) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	
	
	#####frame for inter ############# 
	$frame = $mw->Frame();
	$frame->Label(-text => "INTER-PAIRS:   ")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'pair_inter_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('pair_inter_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'pair_inter_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'pair_inter_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'pair_inter_offset', 3 , \&redraw);
	$tmp=$frame->Checkbutton(-text=>'lines on',-variable => \$opt{'pair_inter_line_on'})->pack(-side=>'left',-anchor=>'e');	
	$ballooni->attach($tmp, -msg => &balloon_format_var('pair_inter_line_on') ) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	
	
	$frame = $mw->Frame();
	$frame->Label(-text => "ALIGNMENT TEXT:   ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"file name col", 'pair_align_path_col', 5 , \&redraw);
	&fast_lentry($frame,$tmp,"base path", 'pair_align_path_col', 30 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	
}


sub card_sub_subscale {
	my $mw=shift;
	my ($frame,$tmp);
	#####SUBJECT BELOW #####################
	$frame = $mw->Frame();
	$tmp=$frame->Label(-text => "-arrangesub",-background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("Arranges subject hits below sequence                                      \n"
	  		  ."    oneperline - each subject accession is placed on a separate line      \n"
	  		  ."    stagger - subjects are staggered(placed on same line if nonoverlapping\n"
	  		  ."              (can use arrangesub to choose column to sort the arrange    \n"
	  		  ."    subscaleN - plots subjects on a continuous numerical scale            \n"
	  		  ."    subscaleC - plots subjects on a categorical (noncontinous) scale      \n"
	  		  ."    *subscaleC#CHR_oo21 plots subjects on basis of chromosome assignment  \n"
	  		  ."                 for Jim Kent assembly but now just good example          \n"
	  		  ."    *subscaleN#ident?? graphs each hit on basis of percent similarity     \n"));
	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'arrangesub'},-background=>'#bbe8ff' ,
			-options => ['oneperline','stagger','subscaleN','subscaleC',
			'*subscaleN#ident90', '*subscaleN#ident85', '*subscaleN#ident80',
			'*subscaleC#CHR_oo21']
		)->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('arrangesub') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"Arrangesub Column", 'arrangesub_col', 10 , \&reshowNredraw,'#bbe8ff');
	$tmp=$frame->Checkbutton(-text=>'reverse sort',-variable => \$opt{'arrangesub_rev_on'})->pack(-side=>'left',-anchor=>'e');	
	$ballooni->attach($tmp, -msg => &balloon_format_var('arrangesub_rev_on') ) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	

	
	$frame = $mw->Frame();	
	$tmp = $frame->Label(-text => "-colorsub",-background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("Choose coloring scheme for subjects                                        \n"
	  		  ."    NONE = no color scheme                                                 \n"
	  		  ."    seqrandom = randomly assign a color to all pairwise for a sequence     \n"	  		
	  		  ."    hitrandom = randomly assign a color to each indvidual pairwise         \n"
	  		  ."    hitconditional = assigns color based on conditional tests of a column  \n"
	  		  ."    RESET = removes color of hits (which overrides subject colors)          "));

	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'colorsub'}, -background=>'#bbe8ff' ,
	       -options => ['NONE','seqrandom','hitrandom','hitconditional','RESET'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('colorsub') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"Default Color", 'sub_color', 10 , \&reshowNredraw,'#bbe8ff');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame();	
	&fast_lentry($frame,$tmp,"Condition Column", 'colorsub_hitcond_col', 10 , \&find_column_options,'#bbe8ff');
	$frame->Label(-textvariable =>\$colheader{'colorsub_hitcond_col_header'})->pack(-side=> 'left',-anchor => 'e');	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	
	$frame = $mw->Frame();	
	&fast_lentry($frame,$tmp,"Condition Tests", 'colorsub_hitcond_tests', 70 , \&reshowNredraw,'#bbe8ff');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	


	$frame = $mw->Frame();
	$tmp=$frame->Label(-text => "SUBJECTS:",-background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("Subjects are the matching regions drawn below the sequence                \n"
	  		  ."  Subject comes from blast (query hitting subjects) additionally          \n"
	  		  ."     they are all drawn below the sequence sub-sequence.                  \n"
	  		  ."  col will be used for seq1 and seq2 if col2 is empty.  If col2 is not    \n"
	  		  ."     empty then col2 will be used for seq2 data.                          \n"
	  		  ."  Patterns must be regular expressions with () enclosing part of match    \n"
	  		  ."    to extract from the string.                                              "));


	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'sub_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"Width", 'sub_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"Line: space:", 'sub_line_spacing', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"init indent", 'sub_initoffset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => "   Arrow")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable=> \$opt{'sub_arrow_on'})
		->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_arrow_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"line paral (-)", 'sub_arrow_paral', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"line diag (/)", 'sub_arrow_diag', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"line perp (|)", 'sub_arrow_perp', 7 , \&redraw);	
	$frame->pack(-side => 'top', -anchor => 'w');

	#########SUBJECT LABELS ##################
	$frame = $mw->Frame();
	$frame->Label(-text => "Label Sequence Begin:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_labelseq_on'})->pack(-side=>'left',-anchor=>'e');	
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_labelseq_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"col", 'sub_labelseq_col', 8 , \&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'sub_labelseq_col_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"col2", 'sub_labelseq_col2', 8 , \&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'sub_labelseq_col2_header'})->pack(-side=> 'left',-anchor => 'e');
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	
	###
	$frame = $mw->Frame();
	&fast_lentry($frame,$tmp,"  pattern", 'sub_labelseq_col_pattern', 14 , \&redraw);
	&fast_lentry($frame,$tmp,"pattern2", 'sub_labelseq_col2_pattern', 14 , \&redraw);
	&fast_lentry($frame,$tmp,"color", 'sub_labelseq_color', 10 , \&redraw);
	&fast_lentry($frame,$tmp,"size", 'sub_labelseq_size', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'sub_labelseq_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => "Label Sequence End:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_labelseqe_on'})->pack(-side=>'left',-anchor=>'e');	
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_labelseqe_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"col", 'sub_labelseqe_col', 8 , \&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'sub_labelseqe_col_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"col2", 'sub_labelseqe_col2', 8 , \&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'sub_labelseqe_col2_header'})->pack(-side=> 'left',-anchor => 'e');
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	
	###
	$frame = $mw->Frame();
	&fast_lentry($frame,$tmp,"  pattern", 'sub_labelseqe_col_pattern', 14 , \&redraw);
	&fast_lentry($frame,$tmp,"pattern2", 'sub_labelseqe_col2_pattern', 14 , \&redraw);
	&fast_lentry($frame,$tmp,"color", 'sub_labelseqe_color', 10 , \&redraw);
	&fast_lentry($frame,$tmp,"size", 'sub_labelseqe_size', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'sub_labelseqe_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1, -anchor=>'w');	


	##### HIT LABELS #######################################
	$frame = $mw->Frame();
	$frame->Label(-text => "Label Each Hits: ")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_labelhit_on'})->pack(-side=>'left',-anchor=>'e');	
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_labelhit_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"col", 'sub_labelhit_col', 8 , \&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'sub_labelhit_col_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"col2", 'sub_labelhit_col2', 8 , \&find_column_options);
	$frame->Label(-textvariable =>\$colheader{'sub_labelhit_col2_header'})->pack(-side=> 'left',-anchor => 'e');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame();
	&fast_lentry($frame,$tmp,"  pattern", 'sub_labelhit_col_pattern', 14 , \&redraw);
	&fast_lentry($frame,$tmp,"pattern2", 'sub_labelhit_col2_pattern', 14 , \&redraw);
	&fast_lentry($frame,$tmp,"color", 'sub_labelhit_color', 10 , \&redraw);
	&fast_lentry($frame,$tmp,"size", 'sub_labelhit_size', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'sub_labelhit_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	
	$frame = $mw->Frame();
	$tmp=$frame->Label(-text => "subscale ", -background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("subscale creates a numeric or categorical plot of sub-subjects below sequence:\n"
	  		  ."(N)umeric scale plots numeric values in a continous fashion                 \n"
	  		  ."(C)ategoric scale overlays a name for each step of the numeric scale        \n"
	  		  ."    categorical names are separated by commas (e.g. X,Y,1,2,3,4)            \n"
	  		  ."col and col2 choose columns containing y values to plot                     \n"
	  		  ."   if col2 is empty col will be used for seq1 and seq2 of pairwise          \n"
	  		  ."   if col2 filled then seq2 will use col2 while col will be for seq1        \n"
	  		  ."                                                                             "));
	&fast_lentry($frame,$tmp,"col", 'sub_scale_col', 8 , \&find_column_options,'#bbe8ff');
	$frame->Label(-textvariable=> \$colheader{'sub_scale_col_header'})->pack(-side=>'left',-anchor=>'e');
	&fast_lentry($frame,$tmp,"pattern", 'sub_scale_col_pattern', 12 , \&reshowNredraw,'#bbe8ff');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame();
	&fast_lentry($frame,$tmp,"col2", 'sub_scale_col2', 8 , \&find_column_options,'#bbe8ff');
	$frame->Label(-textvariable=> \$colheader{'sub_scale_col2_header'})->pack(-side=>'left',-anchor=>'e');
	&fast_lentry($frame,$tmp,"pattern2", 'sub_scale_col2_pattern', 12, \&reshowNredraw,'#bbe8ff');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	
	
	$frame = $mw->Frame();
	$frame->Label(-text => " (N)umeric ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"step", 'sub_scale_step', 5 , \&reshowNredraw,'#bbe8ff');
	&fast_lentry($frame,$tmp,"min", 'sub_scale_min', 8 , \&reshowNredraw,'#bbe8ff');
	&fast_lentry($frame,$tmp,"max", 'sub_scale_max', 8 , \&reshowNredraw,'#bbe8ff');
	&fast_lentry($frame,$tmp,"lines", 'sub_scale_lines', 5 , \&reshowNredraw,'#bbe8ff');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame();
	$frame->Label(-text => "(C)ategoric ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"names", 'sub_scale_categoric_string', 60 , \&reshowNredraw,'#bbe8ff');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	
	
	#######################SCALES ##############################3
	$frame = $mw->Frame();	
	$tmp=$frame->Label(-text => "DRAW SCALES", -background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("subscale creates a numeric or categorical plot of sub-subjects below sequence:\n"
	  		  ."(N)umeric scale plots numeric values in a continous fashion                 \n"
	  		  ."(C)ategoric scale overlays a name for each step of the numeric scale        \n"
	  		  ."col and col2 choose columns containing y values to plot                     \n"
	  		  ."   if col2 is empty col will be used for seq1 and seq2 of pairwise          \n"
	  		  ."   if col2 filled then seq2 will use col2 while col will be for seq1        \n"
	  		  ."                                                                             "));

	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_scale_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_scale_on') ) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => "Vertical Line")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_scale_vline_on'}, -command => \&redraw)->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_scale_vline_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'sub_scale_vline_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'sub_scale_vline_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'sub_scale_vline_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	#######################
	$frame = $mw->Frame();
	$frame->Label(-text => "Horizontal Lines: ")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_scale_hline_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_scale_hline_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'sub_scale_hline_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'sub_scale_hline_width', 8 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();
	$frame->Label(-text => "Tick Marks")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_scale_tick_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_scale_tick_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'sub_scale_tick_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'sub_scale_tick_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"length", 'sub_scale_tick_length', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'sub_scale_tick_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();	
	$frame->Label(-text => "Tick Labels")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'sub_scale_label_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('sub_scale_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'sub_scale_label_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"size", 'sub_scale_label_fontsize', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'sub_scale_label_offset', 4 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame();	
	&fast_lentry($frame,$tmp,"     value_scaling", 'sub_scale_label_multiplier', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"pattern", 'sub_scale_label_pattern', 14 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

 #################################
	$frame = $mw->Frame(-relief => 'groove', -bd => 6 );	
	$frame->Label(-text=> "ALIGN ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$mstring,-wraplength=>450,-justify=>'left')->pack(-side=> 'top',-anchor => 'w');
	$frame->pack(-side => 'top', -anchor => 'w');
}


sub card_extra {
	my $mw=shift;
	my ($frame,$tmp);
	$frame = $mw->Frame();
	$frame->Label(-text => "EXTRA")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'extra_on'})
		->pack(-side=>'left',-anchor=>'e');	
	$ballooni->attach($tmp, -msg => &balloon_format_var('extra_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"default color", 'extra_color', 10 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'extra_width', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'extra_offset', 4 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame();
	$frame->Label(-text => "   Label:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable=> \$opt{'extra_label_on'})
		->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('extra_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"col", 'extra_label_col', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"pattern", 'extra_label_col_pattern', 15 , \&redraw);
		
	&fast_lentry($frame,$tmp,"color", 'extra_label_color', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"size", 'extra_label_fontsize', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'extra_label_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame = $mw->Frame();
	&fast_lentry($frame,$tmp,"  test before label:col", 'extra_label_test_col', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"pattern", 'extra_label_test_pattern', 15 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');		
	$frame = $mw->Frame();
	$frame->Label(-text => "   Arrow")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable=> \$opt{'extra_arrow_on'})
		->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('extra_arrow_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"line paral (-)", 'extra_arrow_para', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"line diag (/)", 'extra_arrow_diag', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"line perp (|)", 'extra_arrow_perp', 7 , \&redraw);	
	$frame->pack(-side => 'top', -anchor => 'w');



	$frame = $mw->Frame(-relief => 'groove', -bd => 6 );	
	$frame->Label(-text=> "EXTRA ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$estring, -wraplength=>450,-justify=>'left' )->pack(-side=> 'top',-anchor => 'w');
	$frame->pack(-side => 'top', -anchor => 'w');
	
}

sub card_graph {
	my $mw=shift;
	my ($frame,$tmp);
	#####SUBJECT BELOW #####################
	$frame = $mw->Frame();
	$tmp=$frame->Label(-text => "GRAPH:",-background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("These are general options controlling the graph scale                     \n"
	  		  ."  pixel height of scale determines the breath of the common graph scale   \n"
	  		  ."  # of intervals determines the number of horizontal lines drawn          \n"
	  		  ."  inital indent determines how far above the line bottom of scale is drawn\n"
	  		  ."     empty then col2 will be used for seq2 data.                          \n"
	  		  ."  Horizontal line is common to both graph1 and graph 2                    \n"
        ));


	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph_scale_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph_scale_on') ) if $opt{'help_on'};

	&fast_lentry($frame,$tmp,"pixel height:", 'graph_scale_height', 5 , \&redraw);
	&fast_lentry($frame,$tmp,"# of intervals", 'graph_scale_interval', 5 , \&redraw);
	&fast_lentry($frame,$tmp,"init indent", 'graph_scale_indent', 5 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	###
	$frame = $mw->Frame()->pack(-side => 'top', -anchor => 'w');
	$frame->Label(-text => "   Common Horz Line:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable=> \$opt{'graph_scale_hline_on'})
		->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph_scale_hline_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"width", 'graph_scale_hline_width', 7 , \&redraw);	
	&fast_lentry($frame,$tmp,"color", 'graph_scale_hline_color', 7 , \&redraw);	

	
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$tmp=$frame->Label(-text => "GRAPH1", -background => '#ffff9d')->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("GRAPH1  reates a numeric or categorical plot of sub-subjects below sequence:\n"
	  		  ."(N)umeric scale plots numeric values in a continous fashion                 \n"
	  		  ."(C)ategoric scale overlays a name for each step of the numeric scale        \n"
	  		  ."    categorical names are separated by commas (e.g. X,Y,1,2,3,4)            \n"
	  		  ."col and col2 choose columns containing y values to plot                     \n"
	  		  ."   if col2 is empty col will be used for seq1 and seq2 of pairwise          \n"
	  		  ."   if col2 filled then seq2 will use col2 while col will be for seq1        \n"
	  		  ."                                                                             "));
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph1_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph1_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"Scale Values: min", 'graph1_min', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"Scale Values: min", 'graph1_max', 8 , \&redraw);
		
	#######################SCALES ##############################3
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "Point:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'graph1_point_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph1_point_on') ) if $opt{'help_on'};
	#not a supported choice
	#$frame->Optionmenu(-textvariable=>\$opt{'graph1_point_shape'},-background=>'#bbe8ff' ,
	#		-options => ['Circle']
	#	)->pack(-side => 'left',-anchor=>'e');
	&fast_lentry($frame,$tmp,"point size", 'graph1_point_size', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"fill color", 'graph1_point_fill_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"outline: color", 'graph1_point_outline_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph1_point_outline_width', 3 , \&redraw);

	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "Connecting Line:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'graph1_line_on'},)->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph1_line_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph1_line_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph1_line_width', 3 , \&redraw);
	$frame->Label(-text => "smoothing")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph1_line_smooth'}, )->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph1_line_smooth') ) if $opt{'help_on'};

	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "Vertical Line")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph1_vline_on'},)->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph1_vline_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph1_vline_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph1_vline_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'graph1_vline_offset', 3 , \&redraw);

	#######################
	$frame = $mw->Frame();
	$frame->Label(-text => "Tick Marks")->pack(-side=> 'left',-anchor => 'w');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph1_tick_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph1_tick_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph1_tick_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"length", 'graph1_tick_length', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph1_tick_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'graph1_tick_offset', 3 , \&redraw);
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	

	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame->Label(-text => "Tick Labels")->pack(-side=> 'left',-anchor => 'w');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph1_label_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph1_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph1_label_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"size", 'graph1_label_fontsize', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'graph1_label_offset', 4 , \&redraw);

	$frame = $mw->Frame()->pack(-side=> 'top',-anchor => 'w');
	&fast_lentry($frame,$tmp,"     value_scaling", 'graph1_label_multiplier', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"decimal points", 'graph1_label_decimal', 14 , \&redraw);	
	
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$tmp=$frame->Label(-text => "GRAPH2", -background => '#ffff9d')->pack(-side=> 'left',-anchor => 'w');
	$ballooni->attach($tmp,-justify => 'left',
	  -msg=>("GRAPH2  reates a numeric or categorical plot of sub-subjects below sequence:\n"
	  		  ."(N)umeric scale plots numeric values in a continous fashion                 \n"
	  		  ."(C)ategoric scale overlays a name for each step of the numeric scale        \n"
	  		  ."    categorical names are separated by commas (e.g. X,Y,1,2,3,4)            \n"
	  		  ."col and col2 choose columns containing y values to plot                     \n"
	  		  ."   if col2 is empty col will be used for seq1 and seq2 of pairwise          \n"
	  		  ."   if col2 filled then seq2 will use col2 while col will be for seq1        \n"
	  		  ."                                                                             "));
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph2_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph2_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"Scale Values: min", 'graph2_min', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"Scale Values: min", 'graph2_max', 8 , \&redraw);
		
	#######################SCALES ##############################3
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "Point:")->pack(-side=> 'left',-anchor => 'w');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'graph2_point_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph2_point_on') ) if $opt{'help_on'};
	#not a supported choice
	#$frame->Optionmenu(-textvariable=>\$opt{'graph2_point_shape'},-background=>'#bbe8ff' ,
	#		-options => ['Circle']
	#	)->pack(-side => 'left',-anchor=>'e');
	&fast_lentry($frame,$tmp,"point size", 'graph2_point_size', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"fill color", 'graph2_point_fill_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"outline: color", 'graph2_point_outline_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph2_point_outline_width', 3 , \&redraw);

	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "Connecting Line:")->pack(-side=> 'left',-anchor => 'w');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'graph2_line_on'},)->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph2_line_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph2_line_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph2_line_width', 3 , \&redraw);
	$frame->Label(-text => "smoothing")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph2_line_smooth'}, )->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph2_line_smooth') ) if $opt{'help_on'};

	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "Vertical Line")->pack(-side=> 'left',-anchor => 'w');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph2_vline_on'},)->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph2_vline_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph2_vline_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph2_vline_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'graph2_vline_offset', 3 , \&redraw);

	#######################
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame->Label(-text => "Tick Marks")->pack(-side=> 'left',-anchor => 'w');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph2_tick_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph2_tick_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph2_tick_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"length", 'graph2_tick_length', 4 , \&redraw);
	&fast_lentry($frame,$tmp,"width", 'graph2_tick_width', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'graph2_tick_offset', 3 , \&redraw);

	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	$frame->Label(-text => "Tick Labels")->pack(-side=> 'left',-anchor => 'w');
	$tmp=$frame->Checkbutton(-text=>'ON',-variable => \$opt{'graph2_label_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('graph2_label_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"color", 'graph2_label_color', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"size", 'graph2_label_fontsize', 3 , \&redraw);
	&fast_lentry($frame,$tmp,"offset", 'graph2_label_offset', 4 , \&redraw);

	$frame = $mw->Frame()->pack(-side=> 'top',-anchor => 'w');
	&fast_lentry($frame,$tmp,"     value_scaling", 'graph2_label_multiplier', 8 , \&redraw);
	&fast_lentry($frame,$tmp,"decimal points", 'graph2_label_decimal', 14 , \&redraw);

   #################################
	#$frame = $mw->Frame(-relief => 'groove', -bd => 6 )->pack(-side=> 'top',-anchor => 'e');	
	#$frame->Label(-text=> "GRAPH input files are fixed:sequence position value")->pack(-side=>'top',-anchor=>'w');
}



sub card_filter {
	my $mw=shift;
	my ($frame,$tmp);
	$frame = $mw->Frame();		
	$frame->Label(-text => "ALIGN FILTERS")->pack(-side=> 'left',-anchor => 'e');
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');

	$frame = $mw->Frame();	
	$tmp=$frame->Checkbutton(-text=>'reset to show all',-variable => \$opt{'pfilter_reset'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('pfilter_reset') ) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	
	$frame = $mw->Frame();	
	$frame->Label(-text => "PRE-ARRANGE-1 ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"col", 'filterpre1_col', 5 ,\&find_column_options,'#bbe8ff');
	$frame->Label(-textvariable => \$colheader{'filterpre1_col_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"min", 'filterpre1_min', 6 , \&redraw,'#bbe8ff');
	&fast_lentry($frame,$tmp,"max", 'filterpre1_max', 6 , \&redraw,'#bbe8ff');	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame = $mw->Frame();		
	$frame->Label(-text => "PRE-ARRANGE-2: ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"col", 'filterpre2_col', 5 , \&find_column_options,'#bbe8ff');
	$frame->Label(-textvariable => \$colheader{'filterpre2_col_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"min", 'filterpre2_min', 6 , \&redraw,'#bbe8ff');
	&fast_lentry($frame,$tmp,"max", 'filterpre2_max', 6 , \&redraw,'#bbe8ff');	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	
	$frame = $mw->Frame();	
	$frame->Label(-text => "POST-ARRANGE-1: ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"col", 'filter1_col', 5 , \&find_column_options);
	$frame->Label(-textvariable => \$colheader{'filter1_col_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"min", 'filter1_min', 6 , \&redraw);
	&fast_lentry($frame,$tmp,"max", 'filter1_max', 6 , \&redraw);	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame = $mw->Frame();	
	$frame->Label(-text => "POST-ARRANGE-2: ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"col", 'filter2_col', 5 , \&find_column_options);
	$frame->Label(-textvariable => \$colheader{'filter2_col_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"min", 'filter2_min', 6 , \&redraw);
	&fast_lentry($frame,$tmp,"max", 'filter2_max', 6 , \&redraw);	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');

	$frame = $mw->Frame();	
	$frame->Label(-text => "EXTRA POST-ARRANGE-1: ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"col", 'filterextra1_col', 5 , \&find_column_options);
	$frame->Label(-textvariable => \$colheader{'filterextra1_col_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"min", 'filterextra1_min', 6 , \&redraw);
	&fast_lentry($frame,$tmp,"max", 'filterextra1_max', 6 , \&redraw);	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame = $mw->Frame();	
	$frame->Label(-text => "EXTRA POST-ARRANGE-2: ")->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"col", 'filterextra2_col', 5 , \&find_column_options);
	$frame->Label(-textvariable => \$colheader{'filterextra2_col_header'})->pack(-side=> 'left',-anchor => 'e');
	&fast_lentry($frame,$tmp,"min", 'filterextra2_min', 6 , \&redraw);
	&fast_lentry($frame,$tmp,"max", 'filterextra2_max', 6 , \&redraw);	
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');

	$frame = $mw->Frame(-relief => 'groove', -bd => 6 );	
	$frame->Label(-text=> "ALIGN ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$mstring,-wraplength=>450,-justify=>'left')->pack(-side=> 'top',-anchor => 'w');
	$frame->Label(-text=> "EXTRA ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$estring, -wraplength=>450,-justify=>'left' )->pack(-side=> 'top',-anchor => 'w');
	$frame->pack(-side => 'top', -anchor => 'w');

}

sub card_misc {
	my $mw=shift;
	my ($frame, $tmp);
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');;
	$frame->Label(-text => "PRINTING") -> pack(-side=>'left',-anchor=>'e');
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');;
	&fast_lentry($frame,$tmp,"Printer:", 'print_command', 12, \&focus_ok );
	&fast_lentry($frame,$tmp,"Multipage # wide:", 'print_multipages_wide', 12, \&focus_ok );
	&fast_lentry($frame,$tmp,"# high:", 'print_multipages_high', 12, \&focus_ok );
	$tmp=$frame->Checkbutton(-text=>'landscape',-variable => \$opt{'printer_page_orientation'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('printer_page_orientation') ) if $opt{'help_on'};

	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');	;
	&fast_lentry($frame,$tmp,"Page width:", 'printer_page_width', 12, \&focus_ok );
	&fast_lentry($frame,$tmp,"height:", 'printer_page_length', 12, \&focus_ok );
		
	$frame = $mw->Frame();
	$frame->Label(-text => "BACKGROUND GIF") -> pack(-side=>'left',-anchor=>'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'gif_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('gif_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"path", 'gif_path', 35, \&redraw );
	&fast_lentry($frame,$tmp,"x", 'gif_x', 6, \&redraw );
	&fast_lentry($frame,$tmp,"y", 'gif_y', 6, \&redraw );
	$frame->Label(-text => "anchor") -> pack(-side=>'left',-anchor=>'e');
	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'gif_anchor'}, -options => ['center','n','e','w','s','ne','nw','se','sw'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('gif_anchor') ) if $opt{'help_on'};
	$frame->pack(-side => 'top', -expand=> 1,-anchor=>'w');	
	
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	$frame->Label(-text =>"LEFT-CLICK MENU OPTIONS")->pack(-side=> 'left',-anchor => 'e');	

	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	&fast_lentry($frame,$tmp,"ALIGNMENTS Bases: query col", 'alignment_col', 5 , \&find_column_options);	
	$frame->Label(-textvariable =>\$colheader{'alignment_col_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"subject col", 'alignment_col2', 5 , \&find_column_options);	
	$frame->Label(-textvariable =>\$colheader{'alignment_col2_header'})->pack(-side=> 'left',-anchor => 'e');	
	&fast_lentry($frame,$tmp,"line wrap width", 'alignment_wrap', 5 );	

	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	$frame->Label(-text =>"FASTA EXTRACT:")->pack(-side=> 'left',-anchor => 'e');	
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'fasta_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('fasta_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"Fasta directory: ", 'fasta_directory', 40 , \&find_column_options);	
	&fast_lentry($frame,$tmp,"frag size", 'fasta_fragsize', 7 , \&find_column_options);	
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	&fast_lentry($frame,$tmp,"Fasta blast db(s): ", 'fasta_blastdb', 40 , \&find_column_options);	
	&fast_lentry($frame,$tmp,"fasta wrap width: ", 'fasta_wrap', 4 , \&find_column_options);	

	
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	&fast_lentry($frame,$tmp,"Execute 1 command", 'execute', 35 , \&find_column_options);	
	&fast_lentry($frame,$tmp,"desc", 'execute_desc', 15 , \&find_column_options);	
	my $tv='pairs';
	$tv='extras' if $opt{'execute_array'} eq 'e';
	my $opt_menu=$frame->Optionmenu(-textvariable=>\$tv,
			-variable=>\$opt{'execute_array'}, -options => [['pairs','m'], ['extras','e']],
			-command => \&find_column_options
			)->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($opt_menu, -msg => &balloon_format_var('execute_array') ) if $opt{'help_on'};
	
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	&fast_lentry($frame,$tmp,"Execute 2", 'execute2', 35 , \&find_column_options);	
	&fast_lentry($frame,$tmp,"desc", 'execute2_desc', 15 , \&find_column_options);	
	my $tv2='pairs';
	$tv2='extras' if $opt{'execute2_array'} eq 'e';
	my $opt_menu=$frame->Optionmenu(-textvariable=>\$tv2,
			-variable=>\$opt{'execute2_array'}, -options => [['pairs','m'], ['extras','e']],
			-command => \&find_column_options
			)->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($opt_menu, -msg => &balloon_format_var('execute2_array') ) if $opt{'help_on'};

	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	&fast_lentry($frame,$tmp,"Execute 3", 'execute3', 35 , \&find_column_options);	
	&fast_lentry($frame,$tmp,"desc", 'execute3_desc', 15 , \&find_column_options);	
	my $tv3='pairs';
	$tv3='extras' if $opt{'execute3_array'} eq 'e';
	my $opt_menu=$frame->Optionmenu(-textvariable=>\$tv3,
			-variable=>\$opt{'execute3_array'}, -options => [['pairs','m'], ['extras','e']],
			-command => \&find_column_options
			)->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($opt_menu, -msg => &balloon_format_var('execute3_array') ) if $opt{'help_on'};
	
	$frame= $mw->Frame()->pack(-side=>'top',-anchor=>'w',-expand=>1,-fill=>'x');	
	&fast_lentry($frame,$tmp,"Execute 4", 'execute4', 35 , \&find_column_options);	
	&fast_lentry($frame,$tmp,"desc", 'execute4_desc', 15 , \&find_column_options);	
	my $tv4='pairs';
	$tv4='extras' if $opt{'execute4_array'} eq 'e';
	my $opt_menu=$frame->Optionmenu(-textvariable=>\$tv4,
			-variable=>\$opt{'execute4_array'}, -options => [['pairs','m'], ['extras','e']],
			-command => \&find_column_options
			)->pack(-side=> 'left',-anchor => 'e');
	$ballooni->attach($opt_menu, -msg => &balloon_format_var('execute4_array') ) if $opt{'help_on'};
	
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "POPUP Window Desc: header format") -> pack(-side=>'left',-anchor=>'e');
	$tmp=$frame->Optionmenu(-textvariable=>\$opt{'popup_format'}, -options => ['number','text'] )->pack(-side => 'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('popup_format') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"max entry length", 'popup_max_len', 5 );	
	
	$frame = $mw->Frame()->pack(-side => 'top', -expand=> 1,-anchor=>'w');
	$frame->Label(-text => "HELP POPUP: ") -> pack(-side=>'left',-anchor=>'e');
	$tmp=$frame->Checkbutton(-text=>'on',-variable => \$opt{'help_on'})->pack(-side=>'left',-anchor=>'e');
	$ballooni->attach($tmp, -msg => &balloon_format_var('help_on') ) if $opt{'help_on'};
	&fast_lentry($frame,$tmp,"line wrap width", 'help_wrap', 5 );	

	$frame = $mw->Frame(-relief => 'groove', -bd => 6 );	
	$frame->Label(-text=> "ALIGN ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$mstring,-wraplength=>450,-justify=>'left')->pack(-side=> 'top',-anchor => 'w');
	$frame->Label(-text=> "EXTRA ARRAY ")->pack(-side=>'top',-anchor=>'w');
	$frame->Label(-textvariable => \$estring, -wraplength=>450,-justify=>'left' )->pack(-side=> 'top',-anchor => 'w');
	$frame->pack(-side => 'top', -anchor => 'w');

}
 	




sub global_edit {
	my $aname=$opt{'edit_arrayname'};
	return if !($aname eq 'm' || $aname eq 'e');
	my $ap =  eval ( '\@'.$aname);
	my @a = @{$ap};
	for (my $i=0;$i<@a;$i++) {
		my $do_it=1;
		my @c = @{ $a[$i] };
		foreach my $v ( qw(A B C D E F) ) {
			#print "$v: $opt{'edit'}{$v}{'col'}\n";
			next if $opt{'edit'}{$v}{'col'}  eq '';
			my $variable = '$c['. $opt{'edit'}{$v}{'col'} . "]";
			my $operation = $opt{'edit'}{$v}{'op'};
			#print "===>$variable $operation\n";
			my $e = eval ( ($variable . ' '.$operation) );
			if ($e eq '') {
				$do_it = 0;
				last;
			}
		}
		# if modify @c;
		#if delete remove this bugger
		# if 
		# if del
		# if modify return 
		#
		
		
	}	
}

sub edit_options { 	
	my $mw=MainWindow->new;
	### frame spacing ####
	$frame = $mw->Frame(-relief => 'groove', -bd => 4 );		
	$frame->Label(-text => "ARRAY")->pack(-side=> 'left',-anchor => 'e');
	my $button=$frame->Button(-text => "COLOR:$opt{'edit_color'}",-background=>$opt{'edit_color'})									
					 ->pack(-side=>'right',-anchor => 'e');
	$button->configure(-command=>[sub { my $b= $_[0]; my $mw=$_[1];
						 my $color = $b->chooseColor(-title=>'Choose New Color',
												-initialcolor=> $opt{'edit_color'} );
						$mw->raise();
						if (defined $color) {
								$b->configure(-background=> $color, -text => "COLOR:$color");
								$opt{'edit_color'}=$color;
								
						}
						
			}, $button, $mw] );
						


	$tmp=$frame->Entry(-textvariable => \$opt{'edit_arrayname'} ,-width=>3)->pack(-side=> 'left', -anchor=> 'e');
	$tmp->bind("<Return>", sub{ my $a= $opt{'edit_arrayname'}; 
											#print "XXXX",$a,"YYYYY\n";
										  	my $ah= '\@'."$a"."header";
										  	#print "XXXX",$ah,"YYYYY";
										  	my $ahp= eval( $ah);
											next if ! defined $ahp;
											$column_header_display='';
										  	for (my $i=0; $i< @{$ahp}; $i++) {
										  		$column_header_display.= "$i)$$ahp[$i]  ";
											}
											#print "\n";
							});
	#sub edit_print_array {
		
	$frame->Button(-text=>'Delete',-command=> [\&global_edit, "Delete"])->pack(-side=>'right',-anchor=>'e');	
	
	$frame->Button(-text=>'Hide',-command=> [\&global_edit, "Hide"])->pack(-side=>'right',-anchor=>'e');	
	
	$frame->Button(-text=>'Test and Modify',-command=> [\&global_edit, 'Test&Mod'])->pack(-side=>'right',-anchor=>'e');	
	$frame->pack(-side => 'top', -expand=> 1);
	
	#########TEST AND OPERATION FRAME ##################
	$frame = $mw->Frame(-relief => 'groove', -bd => 4 );	
	$frame->Label(-text => "colA:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Entry(-textvariable => \$opt{'edit_A_col'},-width=>10)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "opA:")->pack(-side=> 'left',-anchor => 'e');
	$tmp = $frame->Entry(-textvariable => \$opt{'edit_A_op'},-width=>20)->pack(-side=> 'left', -anchor=> 'e');
	$frame->pack(-side => 'top', -expand=> 1);	

	$frame->Label(-text => "colB:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Entry(-textvariable => \$opt{'edit'}{'B'}{'col'},-width=>10)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "opB:")->pack(-side=> 'left',-anchor => 'e');
	$tmp = $frame->Entry(-textvariable => \$opt{'edit'}{'B'}{'op'},-width=>20)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "colC:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Entry(-textvariable => \$opt{'edit'}{'C'}{'col'},-width=>10)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "opC:")->pack(-side=> 'left',-anchor => 'e');
	$tmp = $frame->Entry(-textvariable => \$opt{'edit'}{'C'}{'op'},-width=>20)->pack(-side=> 'left', -anchor=> 'e');

	$frame->pack(-side => 'top', -expand=> 1);	

	
	####TEST AND OPERATION FRAME 2 ###################
	$frame = $mw->Frame(-relief => 'groove', -bd => 4 );	
	$frame->Label(-text => "colD:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Entry(-textvariable => \$opt{'edit'}{'D'}{'col'},-width=>10)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "opD:")->pack(-side=> 'left',-anchor => 'e');
	$tmp = $frame->Entry(-textvariable => \$opt{'edit'}{'D'}{'op'},-width=>20)->pack(-side=> 'left', -anchor=> 'e');
	$frame->pack(-side => 'top', -expand=> 1);	

	$frame->Label(-text => "colE:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Entry(-textvariable => \$opt{'edit'}{'E'}{'col'},-width=>10)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "opE:")->pack(-side=> 'left',-anchor => 'e');
	$tmp = $frame->Entry(-textvariable => \$opt{'edit'}{'E'}{'op'},-width=>20)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "colF:")->pack(-side=> 'left',-anchor => 'e');
	$tmp=$frame->Entry(-textvariable => \$opt{'edit'}{'F'}{'col'},-width=>10)->pack(-side=> 'left', -anchor=> 'e');
	
	$frame->Label(-text => "opF:")->pack(-side=> 'left',-anchor => 'e');
	$tmp = $frame->Entry(-textvariable => \$opt{'edit'}{'F'}{'op'},-width=>20)->pack(-side=> 'left', -anchor=> 'e');
	$frame->pack(-side => 'top', -expand=> 1);
	
	$frame = $mw->Frame(-relief => 'groove', -bd => 4 );	
	$frame->Label(-textvariable => \$column_header_display, -wraplength=>800 )->pack(-side=> 'left',-anchor => 'e');
	$frame->pack(-side => 'top', -expand=> 1);
}

############################################################	
##############################DRAWING ######################
sub update {
	$mw->update;
	$canvas->update;
	print "VIEW UPDATE\n" if !$newopt{'quiet'};
}

sub reshowNredraw {
	if ( (keys %pairwise2delete) > 0) {
		foreach (reverse sort {$a <=> $b} keys %pairwise2delete) {
			print "$_ ";
			splice(@m,$_,1);
		}
		%pairwise2delete=();
		
	}
	#print "ALIGN_UPDATE\n" ;
	&align_update;
	$opt{'align'}='';
	#print "show_UPDATE\n";
	&show_update;
	&show_update_subject;	
	#print "extra UPDATE\n";
	&extra_update;
	$opt{'extra'}='';
	#print "graph UPDATE\n";
	&graph_update;
	$opt{'graph'}='';
	&find_column_options;
	print "DISPLAY CALCULATIONS\n"  if !$newopt{'quiet'};
	&show_calculations;

	print "**** REDRAWING ****\n"  if !$newopt{'quiet'};
	foreach my $ac (@acc_order) {
		foreach my $s ( keys % {$acc{$ac}{'sub'}} ) {
		$accsub{$s}{'acc'}{$ac}{'qmin_f'}=99999999999999;
		$accsub{$s}{'acc'}{$ac}{'qmax_f'}=-10;
		}
	}
	&accession_begin_end;
	&arrange_seq_lines;
	&arrange_subjects if $opt{'sub_on'};
	&colorsub;
	&generate_column_header_strings;
	
   &redraw;
   $first_pass=0;
}

sub redraw {
	####RESET CERTAIN VALUES###
	&find_column_options;
	%msghash=();
	$scale=1;

	$canvas->delete("all");
	if ($opt{'filename_on'}) {
		$opt{'filename'}=$filepath if $filepath && !$opt{'filename'};
		my $fn=$opt{'filename'};
		$fn=$filepath if $filepath && !$opt{'filename'};
		$fn=$1 if $fn =~ /$opt{'filename_pattern'}/ && $opt{'filename_pattern'};
		#print "FILENAME $fn\n";
		$canvas->createText(1+$opt{'filename_offset_h'},1+$opt{'filename_offset'}, -text=>$fn, -fill=>$opt{'filename_color'},-font=>"Courier $opt{'filename_size'}", -anchor=>"nw");
	}
	if ($opt{'gif_on'}) {
		my $photo=$canvas->createImage($opt{'gif_x'},$opt{'gif_y'},-image=> $mw->Photo(-file => $opt{'gif_path'}),
		    -anchor => $opt{'gif_anchor'});
	}
	$canvas_width = $opt{'window_width'} -$opt{'canvas_indent_left'} -$opt{'canvas_indent_right'};
	$bp_per_pixel = $opt{'canvas_bpwidth'} / $canvas_width;
	
	&draw_sequences; #draw sequences also draws graph scales
	&draw_scales;
	&draw_pairwise; #draw pairwise also draws subs#
	&draw_extra;
	&draw_subject_labels if$opt{'sub_on'} && (  $opt{'sub_labelseq_on'} || $opt{'sub_labelseqe_on'} ) && $opt{'arrangesub'} !~/statistic/;
	
	&draw_graph;
	my $drawtext=eval( '"'. $opt{'text_text'} . '"' );
	if ($opt{'text_on'} ) {
		$canvas->createText($opt{'text_offset_h'},$opt{'text_offset'},
				-text=> $drawtext, 
				-fill=>$opt{'text_color'},
				-font=>"Courier $opt{'text_fontsize'}", 
				-anchor=>$opt{'text_anchor'}
		 );
	}
	$drawtext=eval( '"'. $opt{'text2_text'} . '"' );
	if ($opt{'text2_on'} ) {
		$canvas->createText($opt{'text2_offset_h'},$opt{'text2_offset'},
				-text=> $drawtext, 
				-fill=>$opt{'text2_color'},
				-font=>"Courier $opt{'text2_fontsize'}", 
				-anchor=>$opt{'text2_anchor'}
		 );
	}

	
	
	
	$canvas->configure(-scrollregion=>[$canvas->bbox("all")]);
	#$canvas->delete ("2");
}

sub accession_begin_end {
#	if ($opt{'just_pairwise'} == 0 ) {
		#print "ALL OF SEQUENCE PLOTTED\n";
		foreach ( keys %acc) {
			$acc{$_}{'b'}=1 if !defined $acc{$_}{'b'};
			$acc{$_}{'e'}=$acc{$_}{'len'} if !defined $acc{$_}{'e'};
			#print "$_  $acc{$_}{'b'} $acc{$_}{'e'}\n";
		}
}


sub arrange_seq_lines {	
	my $i=0;
	@l=();
	%lpos=();
	if ($opt{'arrangeseq'} eq 'oneperline') {
		foreach ( @acc_order ) {
			@{$l[$i]{'acc'}}=($_);
			$i++;
		}
	} elsif ($opt{'arrangeseq'} eq 'sameline') { 
		@{$l[0]{'acc'}}=(@acc_order );
	} elsif  ($opt{'arrangeseq'} eq 'file' ) {
		print "ARRANGING USING FILE ($opt{'arrange_file'}\n";
		open (ARRANGE, $opt{'arrange_file'} ) || die "Can't read ($opt{'arrange_file'})!\n";
		#file that is inputted contains seqname begin#
		#newlines are represented NEWLINE#
		my $head=<ARRANGE>;
		print "HEAD$head";
		while (<ARRANGE>) {
			s/\r\n/\n/;chomp; chomp;
			my @c=split /\t/;
			print "$_\n";
			if ($c[0] =~ /NEWLINE/) {$i++;next;}
			next if $c[0] !~ /\S+/;
			die "$opt{'arrange_file'} contains $_ which has no digit data for position\n" if $c[1] !~ /^\d+$/;
			push @{$l[$i]{'acc'}},$c[0];
			$lpos{$c[0]}=$c[1];
		}
		close ARRANGE;
		#####################################
		my $count=0;
		foreach (keys %acc) {
			if (! defined $lpos{$_}) {
				#print "WARNING\n";
				$count++;
				&warnNpause( "ERROR:$_ does not have  line and base start position in ($opt{'arrange_file'})!\n");
			}
		}
		if ($count > 0) {
			print uc "WARNING missing position data for show seq in ($opt{'arrange_file'})\n";
			print uc "Display Results will be unpredicable!!!!!!!\n";
			print "PRESS RETURN AND CONTINUE AT OWN RISK. (Cntl-C to Quit)  ";
			my $pause=<STDIN>;
		}
	} else {
		die "Can not arrange with -arrangeseq $opt{'arrangeseq'}!\n";
	}
}


sub draw_sequences {
	$widest_line=0;
	my %sub_scale_hash=();
	my %graph_scale_hash=();
	print "DRAWING SEQUENCES...\n"  if !$newopt{'quiet'};
	#variables for changing globally eventually
	my $wrap_on=1;
	my $color=$opt{'seq_color'};
	my $thickness=$opt{'seq_width'};
	###the subroutine itself
	my $liney=$opt{'canvas_indent_top'};
	foreach (my $l=0; $l<@l; $l++) {
		$l[$l]{'liney'} = $liney;
		$l[$l]{'maxbp'}=0;
		my $bp_x=1;
		foreach my $a( @{ $l[$l]{'acc'} }) {
			#####initial calc of line and bp
			$acc{$a}{'l'} = $l;
			if ( defined $lpos{$a} ) {
					#print "LPOS$a => $lpos{$a}\n";
					$bp_x=$lpos{$a};
			}
			$acc{$a}{'xb'}=$bp_x;
			$bp_x+=($acc{$a}{'e'}-$acc{$a}{'b'});
			$acc{$a}{'xe'} = $bp_x;
			$l[$l]{'max'}=$bp_x if $bp_x > $l[$l]{'max'};
			$bp_x+= $opt{'seq_spacing_btwn_sequences'}+1 if  !defined $lpos{$a};
			

			my $xb=$acc{$a}{'xb'};
			my $xe=$acc{$a}{'xe'};
			$widest_line=$xe if $xe>$widest_line;
			my $begin=$acc{$a}{'b'};
			&draw_line_horz_pieces($l,$xb,$l,$xe, "s$a", $color, $thickness);
			my $tagname = "s$a";  #can't go inside a sub!!!
			$canvas->addtag('seq','withtag',$tagname);

			my ($x1,$y1)=&linexbp2xy($l,$xb);
			my ($x2,$y2)=&linexbp2xy($l,$xe);
				
			if ($opt{'seq_tick_whole'} == 0) {
				if ($opt{'seq_tick_b_on'}) {
					my $tag=$canvas->createLine($x1,$y1+$opt{'seq_tick_b_offset'}, $x1,$y1+$opt{'seq_tick_b_length'}+$opt{'seq_tick_b_offset'},
							-width => $opt{'seq_tick_b_width'}, -fill => $opt{'seq_tick_b_color'});
					$canvas->addtag("tick","withtag",$tag);
					if ($opt{'seq_tick_b_label_on'} ) {
						$tag=$canvas->createText($x1+$opt{'seq_tick_b_label_offset_h'}, $y1+$opt{'seq_tick_b_length'} + $opt{'seq_tick_b_label_offset'},
								-text=>int($begin*$opt{'seq_tick_b_label_multiplier'}), 
								-fill=>$opt{'seq_tick_b_label_color'},
								-font=>"Courier $opt{'seq_tick_b_label_fontsize'}", 
								-anchor=>$opt{'seq_tick_b_label_anchor'}
						 );
						$canvas->addtag("tickl","withtag",$tag);
					}
				}
				if ($opt{'seq_tick_e_on'}) {
					my $tag=$canvas->createLine($x2,$y2+$opt{'seq_tick_e_offset'}, $x2,$y2+$opt{'seq_tick_e_length'}+$opt{'seq_tick_e_offset'},
							-width => $opt{'seq_tick_e_width'}, -fill => $opt{'seq_tick_e_color'});
					$canvas->addtag("tick","withtag",$tag);
					if ($opt{'seq_tick_e_label_on'} ) {
						$tag=$canvas->createText($x2+$opt{'seq_tick_e_label_offset_h'}, $y2+$opt{'seq_tick_e_length'} + $opt{'seq_tick_e_label_offset'},
								-text=>int(($begin+$xe-$xb)*$opt{'seq_tick_e_label_multiplier'} ), 
								-fill=>$opt{'seq_tick_e_label_color'},
								-font=>"Courier $opt{'seq_tick_e_label_fontsize'}", 
								-anchor=>$opt{'seq_tick_e_label_anchor'}
						 ) ;
						$canvas->addtag("tickl","withtag",$tag);
					}
				}

						#print "=======>   ($x1,$y1)    ($x2,$y2)\n";
				if ($opt{'seq_tick_on'} ) {
					for (my $t=$xb-1+$opt{'seq_tick_bp'}-(($begin) % $opt{'seq_tick_bp'}); $t<= $xe; $t+=$opt{'seq_tick_bp'}) {
						my ($xt, $yt)=&linexbp2xy($l, $t);
						my $tag=$canvas->createLine($xt,$yt+$opt{'seq_tick_offset'}, $xt,$yt+$opt{'seq_tick_length'}+$opt{'seq_tick_offset'},
								-width => $opt{'seq_tick_width'}, -fill => $opt{'seq_tick_color'});
						$canvas->addtag('tick',"withtag",$tag);
						if ($opt{'seq_tick_label_on'} ) {
							$tag=$canvas->createText($xt, $yt+$opt{'seq_tick_length'} + $opt{'seq_tick_label_offset'},
								-text=>  ($t-$xb+1+$begin)*$opt{'seq_tick_label_multiplier'}, 
								-fill=> $opt{'seq_tick_label_color'},
								-font=>"Courier $opt{'seq_tick_label_fontsize'}", 
								-anchor=>$opt{'seq_tick_label_anchor'}
							);
							$canvas->addtag('tickl','withtag',$tag);
						}
					}
				}
			}
		#	foreach (qw(seq_label_on seq_label_pattern)) {
		#		print "$_=>$opt{$_}####\n";
		#	}
			if ($opt{'seq_label_on'}) {
		#		print "A:$a\n";
				my $name=$a;
				$name=$1 if $opt{'seq_label_pattern'} ne '' && $name=~/$opt{'seq_label_pattern'}/ ;
				my $tag=$canvas->createText($x1+$opt{'seq_label_offset_h'}, $y1-6+$opt{'seq_label_offset'},
						-text=>$name, 
						-font=>"Courier $opt{'seq_label_fontsize'}",
						-fill=> $opt{'seq_label_color'}, -anchor=>"sw");
				$canvas->addtag("seqn","withtag",$tag);
		#		print "SEQN:$name\n";
			}	

		}
		$liney+= int(($bp_x-1)/$opt{'canvas_bpwidth'}) * $opt{'seq_line_spacing_wrap'} + $opt{'seq_line_spacing_btwn'};
	}
		#my $pause=<STDIN>;
}

sub draw_scales {
	print "DRAWING SCALES FOR GRAPHS AND SUBJECTS...\n"  if !$newopt{'quiet'};
	foreach (my $l=0; $l<@l; $l++) {
		my ($x1,$y1)=&linexbp2xy($l,1);
		my ($x2,$y2)=&linexbp2xy($l,$l[$l]{'max'});
		my $xb=1;
		my $xe=$l[$l]{'max'};
		#print "$l $y1 $y2\n";
		if ($opt{'seq_tick_whole'} == 1) {
			if ($opt{'seq_tick_b_on'}) {
				my $tag=$canvas->createLine($x1,$y1+$opt{'seq_tick_b_offset'}, $x1,$y1+$opt{'seq_tick_b_length'}+$opt{'seq_tick_b_offset'},
						-width => $opt{'seq_tick_b_width'}, -fill => $opt{'seq_tick_b_color'});
				$canvas->addtag("tick","withtag",$tag);
				if ($opt{'seq_tick_b_label_on'} ) {
					$tag=$canvas->createText($x1+$opt{'seq_tick_b_label_offset_h'}, $y1+$opt{'seq_tick_b_length'} + $opt{'seq_tick_b_label_offset'},
							-text=>int($xb*$opt{'seq_tick_b_label_multiplier'}), 
							-fill=>$opt{'seq_tick_b_label_color'},
							-font=>"Courier $opt{'seq_tick_b_label_fontsize'}", 
							-anchor=>$opt{'seq_tick_b_label_anchor'}
					 );
					$canvas->addtag("tickl","withtag",$tag);
				}
			}
			if ($opt{'seq_tick_e_on'}) {
				my $tag=$canvas->createLine($x2,$y2+$opt{'seq_tick_e_offset'}, $x2,$y2+$opt{'seq_tick_e_length'}+$opt{'seq_tick_e_offset'},
						-width => $opt{'seq_tick_e_width'}, -fill => $opt{'seq_tick_e_color'});
				$canvas->addtag("tick","withtag",$tag);
				if ($opt{'seq_tick_e_label_on'} ) {
					$tag=$canvas->createText($x2+$opt{'seq_tick_e_label_offset_h'}, $y2+$opt{'seq_tick_e_length'} + $opt{'seq_tick_e_label_offset'},
							-text=>int(($xe +0.000001)*$opt{'seq_tick_e_label_multiplier'}) , 
							-fill=>$opt{'seq_tick_e_label_color'},
							-font=>"Courier $opt{'seq_tick_e_label_fontsize'}", 
							-anchor=>$opt{'seq_tick_e_label_anchor'}
					 ) ;
					$canvas->addtag("tickl","withtag",$tag);
				}
			}

					#print "=======>   ($x1,$y1)    ($x2,$y2)\n";
			if ($opt{'seq_tick_on'} ) {
				for (my $t=$xb-1+$opt{'seq_tick_bp'}; $t<= $xe; $t+=$opt{'seq_tick_bp'}) {
					my ($xt, $yt)=&linexbp2xy($l, $t);
					my $tag=$canvas->createLine($xt,$yt+$opt{'seq_tick_offset'}, $xt,$yt+$opt{'seq_tick_length'}+$opt{'seq_tick_offset'},
							-width => $opt{'seq_tick_width'}, -fill => $opt{'seq_tick_color'});
					$canvas->addtag('tick',"withtag",$tag);
					if ($opt{'seq_tick_label_on'} ) {
						$tag=$canvas->createText($xt, $yt+$opt{'seq_tick_length'} + $opt{'seq_tick_label_offset'},
							-text=>  ($t-$xb+1)*$opt{'seq_tick_label_multiplier'}, 
							-fill=> $opt{'seq_tick_label_color'},
							-font=>"Courier $opt{'seq_tick_label_fontsize'}", 
							-anchor=>$opt{'seq_tick_label_anchor'}
						);
						$canvas->addtag('tickl','withtag',$tag);
					}
				}
			}
	
			#if ($opt{'seq_label_on'}) {
			#	my $name=$_;
			#	$name=$1 if $opt{'seq_label_pattern'} ne '' && $name=~/$opt{'seq_label_pattern'}/ ;
			#	my $tag=$canvas->createText($x1+$opt{'seq_label_offset_h'}, $y1-6+$opt{'seq_label_offset'},
			#			-text=>$name, 
			#			-font=>"Courier $opt{'seq_label_fontsize'}",
			#			-fill=> $opt{'seq_label_color'}, -anchor=>"sw");
			#	$canvas->addtag("seqn","withtag",$tag);
			#}	

		}
		if ( $opt {'graph_scale_on'} ) {
			my $intercept1=$opt{'graph1_min'};
			my $slope1=($opt{'graph1_max'}-$opt{'graph1_min'})/$opt{'graph_scale_height'};
			my $scale1= $opt{'graph_scale_height'} / ( $opt{'graph1_max'}-$opt{'graph1_min'} );
			
			my $intercept2=$opt{'graph2_min'};
			my $slope2=($opt{'graph2_max'}-$opt{'graph2_min'})/$opt{'graph_scale_height'};
			my $scale2= $opt{'graph_scale_height'} / ( $opt{'graph2_max'}-$opt{'graph2_min'} );
			
			my $step = $opt{'graph_scale_height'}/$opt{'graph_scale_interval'};
			
			for(my $i=$y1; $i<=$y2; $i+=$opt{'seq_line_spacing_wrap'}) {
				#next if defined $graph_scale_hash{$i};
				#$graph_scale_hash{$i}='done';
				#print "GSCALE $i\n";
				my $xbegin= $opt{'canvas_indent_left'}+$opt{'sub_scale_vline_offset'};
				my $xend = $opt{'window_width'}-$opt{'canvas_indent_right'};
				if ($i == $y2) { $xend=$x2};
				#print "DRAWING VLINE\n";
				if ($opt{'graph1_on'} && $opt{'graph1_vline_on'}) {
					my $tag=$canvas->createLine ( $xbegin,$i+$opt{'graph_scale_indent'},
						$xbegin,$i-$opt{'graph_scale_height'}+$opt{'graph_scale_indent'},
						-width => $opt{'graph1_vline_width'}, -fill => $opt{'graph1_vline_color'}
					);
					$canvas->addtag('gs','withtag',$tag);
				}
				if ($opt{'graph2_on'} && $opt{'graph2_vline_on'}) {
					my $tag=$canvas->createLine ( $xend,$i+$opt{'graph_scale_indent'},
						$xend,$i-$opt{'graph_scale_height'}+$opt{'graph_scale_indent'},
						-width => $opt{'graph2_vline_width'}, -fill => $opt{'graph2_vline_color'}
					);
					$canvas->addtag('gs','withtag',$tag);
				}
				for (my $j=0; $j<=$opt{'graph_scale_height'};$j+=$step) {
					my $y=$i+$opt{'graph_scale_indent'} -$j;
				#	print "$j  => $l => $y\n";
					if ($opt{'graph_scale_hline_on'}) {
						my $tag=$canvas->createLine($xbegin, $y,$xend,$y,
									-width => $opt{'graph_scale_hline_width'}, -fill => $opt{'graph_scale_hline_color'});
							$canvas->addtag('gs','withtag',$tag);
					}
					if ($opt{'graph1_on'} && $opt{'graph1_tick_on'}) {
						#print "DRAWING GSCALE1 TICK\n";
						my $tag=$canvas->createLine($xbegin-$opt{'graph1_tick_length'}+$opt{'graph1_tick_offset'}, $y
								,$xbegin+$opt{'graph1_tick_offset'},$y,
							-width => $opt{'graph1_tick_width'}, -fill => $opt{'graph1_tick_color'}
						);							
						$canvas->addtag('gs','withtag',$tag);
					}
					if ($opt{'graph2_on'} && $opt{'graph2_tick_on'}) {
						#print "DRAWING GSCALE1 TICK\n";
						my $tag=$canvas->createLine($xend-$opt{'graph2_tick_length'}+$opt{'graph2_tick_offset'}, $y
								,$xend+$opt{'graph2_tick_offset'},$y,
							-width => $opt{'graph2_tick_width'}, -fill => $opt{'graph2_tick_color'}
						);							
						$canvas->addtag('gs','withtag',$tag);
					}
					if ($opt{'graph1_on'} && $opt{'graph1_label_on'}) {
						my $label=$j*$slope1+$intercept1; ;
						$label=int($label* $opt{'graph1_label_multiplier'}*10**$opt{'graph1_label_decimal'}+0.00000000000001 )/10**$opt{'graph1_label_decimal'};
						my $tag=$canvas->createText(
							$xbegin-$opt{'graph1_tick_length'}+$opt{'graph1_tick_offset'}-2+$opt{'graph1_label_offset'}, $y,
							-text=>$label, 
							-font=>"Courier $opt{'graph1_label_fontsize'}", 
							-fill=> $opt{'graph1_label_color'},-anchor=>"e");
						$canvas->addtag('gsl','withtag',$tag);
					}
					if ($opt{'graph2_on'} && $opt{'graph2_label_on'}) {
						my $label=$j*$slope2+$intercept2; ;
						$label=int($label* $opt{'graph2_label_multiplier'}*10**$opt{'graph2_label_decimal'}+0.00000000000001 )/10**$opt{'graph2_label_decimal'};
						my $tag=$canvas->createText(
							$xend-$opt{'graph2_tick_length'}+$opt{'graph2_tick_offset'}-2+$opt{'graph2_label_offset'}, $y,
							-text=>$label, 
							-font=>"Courier $opt{'graph2_label_fontsize'}", 
							-fill=> $opt{'graph2_label_color'},-anchor=>"w");
						$canvas->addtag('gsl','withtag',$tag);
					}
				}
			}
		}
		if ( $opt{'sub_scale_on'} ) {
			for(my $i=$y1; $i<=$y2; $i+=$opt{'seq_line_spacing_wrap'}) {
				######################################################################
				my $scale= $opt{'sub_scale_lines'} / ( $opt{'sub_scale_max'}-$opt{'sub_scale_min'} );
				my $step = $opt{'sub_scale_step'};
				my $xbegin= $opt{'canvas_indent_left'}+$opt{'sub_scale_vline_offset'};
				my $xend = $opt{'window_width'}-$opt{'canvas_indent_right'};
				if ($i == $y2) { $xend=$x2};

				#########draw vertical lines of scale #####################
				#top is $opt{'sub_intoffset'}
				#bottom is $opt{'sub_scale_lines'}	
				if ($opt{'sub_scale_vline_on'}) {
					my $tag=$canvas->createLine ( $xbegin,$i+$opt{'sub_initoffset'},
						$xbegin,$i+$opt{'sub_scale_lines'}*$opt{'sub_line_spacing'}+$opt{'sub_initoffset'},
						-width => $opt{'sub_scale_vline_width'}, -fill => $opt{'sub_scale_vline_color'}
					);
					$canvas->addtag('ss','withtag',$tag);
				}
				#########draw horizontal scale lines  and ticks and labels#######################
				my ($min,$max)=($opt{'sub_scale_min'},$opt{'sub_scale_max'});
				($min,$max)=($max,$min) if $min > $max;
				for (my $j=$opt{'sub_scale_min'}; $j<=$max && $j>=$min;$j+=$step) {
				#need max to min step1#
					my $l = ($opt{'sub_scale_max'} - $j)*$scale;
					my $y= $l*$opt{"sub_line_spacing"}+$opt{'sub_initoffset'} +$i;
					#print "$j  => $l => $y\n";
					if ($opt{'sub_scale_hline_on'}) {
						my $tag=$canvas->createLine($xbegin, $y,$xend,$y,
							-width => $opt{'seq_tick_width'}, -fill => $opt{'sub_scale_hline_color'});
						$canvas->addtag('ss','withtag',$tag);
					}
					if ($opt{'sub_scale_tick_on'}) {
						my $tag=$canvas->createLine($xbegin-$opt{'sub_scale_tick_length'}+$opt{'sub_scale_tick_offset'}, $y
							,$xbegin+$opt{'sub_scale_tick_offset'},$y,
						-width => $opt{'sub_scale_tick_width'}, -fill => $opt{'sub_scale_tick_color'});							
						$canvas->addtag('ss','withtag',$tag);

					}
					if ($opt{'sub_scale_label_on'}) {
						my $label;
						if ($opt{'arrangesub'}=~/subscaleN/) {
							$label=int(($j+0.000000000001)*1000000)/1000000*$opt{'sub_scale_label_multiplier'};
						} else {
							$label=$subscaleC[$j];
						}
						my $tag=$canvas->createText($xbegin-$opt{'sub_scale_tick_length'}+$opt{'sub_scale_tick_offset'}-2+$opt{'sub_scale_label_offset'}, $y,
							-text=>$label, 
						-font=>"Courier $opt{'sub_scale_label_fontsize'}", 
						-fill=> $opt{'sub_scale_label_color'},-anchor=>"e");
						$canvas->addtag('ssl','withtag',$tag);
					}
				}
				#########draw vertical lines of scale #####################
				#top is $opt{'sub_intoffset'}
				#bottom is $opt{'sub_scale_lines'}
				if ($opt{'sub_scale_vline_on'}) {
					my $tag=$canvas->createLine ( $xbegin,$i+$opt{'sub_initoffset'},
							$xbegin,$i+$opt{'sub_scale_lines'}*$opt{'sub_line_spacing'}+$opt{'sub_initoffset'},
							-width => $opt{'sub_scale_vline_width'}, -fill => $opt{'sub_scale_vline_color'}
					) ;
					$canvas->addtag('ss','withtag',$tag);

				}
			}
		}
	}
}


sub draw_pairwise {
	print "DRAWING PAIRWISE AND SUBJECTS...\n"  if !$newopt{'quiet'};
	my $size=scalar(@m);
	#my $m_pointer=\@m;
	#my $m_point_point= \$m_pointer;
	#print "ARRAYPOINTER M:$m_pointer\n";
	my ($n1,$b1,$e1,$l1,$n2,$b2,$e2,$l2,$skip_this ) ;
	for (my $i=0; $i < $size;$i++) {
		($n1,$b1,$e1,$l1,$n2,$b2,$e2,$l2 ) = @{$m[$i]};
	#	print "$i) $n1 $b1 $e1 $l1 $n2 $b2 $e2 $l2\n";
		print "..$i\n" if $i % 1000 ==0 && $i !=0 &&  !$newopt{'quiet'};
		next if $m[$i][$mh{'hide'}];
		my $defined1=(defined $acc{$n1});
		my $defined2=(defined $acc{$n2});
		###calculate whether sequence is present###
		if ($defined1) {
			if ( 		($b1 < $acc{$n1}{'b'}  && $e1 < $acc{$n1}{'b'} ) 
					|| ($b1 > $acc{$n1}{'e'}  && $e1 > $acc{$n1}{'e'} ) ) {
								$defined1=0 ;
			} else {
				$b1=$acc{$n1}{'b'} if $b1 < $acc{$n1}{'b'} ;
				$e1=$acc{$n1}{'e'} if $e1 > $acc{$n1}{'e'} ;
			}
		}
		if ($defined2) {	
			#print "$b2-$e2\n";
			if (   	($b2 < $acc{$n2}{'b'}  && $e2 < $acc{$n2}{'b'}  ) 
					|| ($b2 > $acc{$n2}{'e'}  && $e2 > $acc{$n2}{'e'}  ) ) { 
					 	 $defined2=0 ;
			} else {
				if ($b2 < $e2 ) {
					$b2=$acc{$n2}{'b'} if $b2 < $acc{$n2}{'b'} ;
					$e2=$acc{$n2}{'e'} if $e2 > $acc{$n2}{'e'} ;
				} else {
					#print "else $b2-$e2\n";
					$b2=$acc{$n2}{'e'} if $b2 > $acc{$n2}{'e'} ;
					$e2=$acc{$n2}{'b'} if $e2 < $acc{$n2}{'b'} ;
					
					
				}
			}
				
		}
		next if  !( ($defined1 && defined $accsub{$n2} ) || ($defined2 && defined $accsub{$n1}) );
		next if  ( !$defined1 && !$defined2 );	
		###convert to proper xposition ####
		#print "$n1=> $acc{$n1}{'xb'} $acc{$n1}{'b'}-$acc{$n1}{'e'} ($b1-$e1)\n";
		#print "$n2=> $acc{$n2}{'xb'}  $acc{$n2}{'b'}-$acc{$n2}{'e'} ($b2-$e2)\n";
		if (defined $acc{$n1}) {
			$b1= $acc{$n1}{'xb'}+$b1 -$acc{$n1}{'b'};
			$e1= $acc{$n1}{'xb'}+$e1 -$acc{$n1}{'b'};
		}
		if (defined $acc{$n2}) {
			$b2= $acc{$n2}{'xb'}+$b2 -$acc{$n2}{'b'};
			$e2= $acc{$n2}{'xb'}+$e2 -$acc{$n2}{'b'};	
		}
		#print "ACC ($n1) $acc{$n1}{'b'}-$acc{$n1}{'e'} $b1-$e1\n";
 		#############################################################
		##########FILTERING TO OCCUR#################################
		#############################################################
		$skip_this=0;
		if ($opt{'filter1_col'} =~/^\d+$/) {
			$skip_this=1 if $m[$i][$opt{'filter1_col'}]< $opt{'filter1_min'} && $opt{'filter1_min'} ne '';
			$skip_this=1 if $m[$i][$opt{'filter1_col'}]> $opt{'filter1_max'} && $opt{'filter1_max'} ne '';
		}
		if ($opt{'filter2_col'} =~/^\d+$/) {
			$skip_this=1 if $m[$i][$opt{'filter2_col'}]< $opt{'filter2_min'} && $opt{'filter2_min'} ne '';
			$skip_this=1 if $m[$i][$opt{'filter2_col'}]> $opt{'filter2_max'} && $opt{'filter2_max'} ne '';
		}
		next if $skip_this;
		##############################################################
		#################DRAW SUBJECTS################################
		##############################################################
		#print "$n1 => $acc{$n1}{'xb'} $b1  $acc{$n1}{'b'} ($b1-$e1)\n";
		#print "$n2=> $acc{$n2}{'xb'} $b1  $acc{$n2}{'b'} ($b2-$e2)\n";

		if ($opt{'sub_on'}==1) {
			my $color=$opt{'sub_color'};			
			if ($defined1 && defined $accsub{$n2} ) {
				########N2 is the SUB#####################
				#next if !$defined1;  #don't draw if if $n1 isn't being drawn
				if ($m[$i][$mh{'scolor'}] ) {
					$color=$m[$i][$mh{'scolor'}];
				} else {
					$color=$accsub{$n2}{'color'} if $accsub{$n2}{'color'};
				}
				my $offset= $accsub{$n2}{'acc'}{$n1}{'line'}*$opt{'sub_line_spacing'}+$opt{'sub_initoffset'};
				$offset = $m[$i][$mh{'sline'}] * $opt{'sub_line_spacing'}+$opt{'sub_initoffset'}if $m[$i][$mh{'sline'}] ne '';
				my $width= $opt{'sub_width'};
				my ($start,$stop)=($b1,$e1);
				($stop,$start)=($start,$stop) if $start> $stop;
				#print "$n1:$n2 C$color W$width O$offset\n";
				($stop,$start)=($start,$stop) if $start> $stop;
				$accsub{$n2}{'acc'}{$n1}{'qmin_f'}=$start if $start < $accsub{$n2}{'acc'}{$n1}{'qmin_f'} ;
				$accsub{$n2}{'acc'}{$n1}{'qmax_f'}=$stop if $stop > $accsub{$n2}{'acc'}{$n1}{'qmax_f'};
				#print "$start $stop\n";
				if ( $opt{'sub_arrow_on'} ) {
					my $arrow='first';
					$arrow='last' if $b2<$e2;
					#print "$acc{$n1}{'l'} $start $acc{$n1}{'l'} $stop\n";
					&draw_line_horz_pieces($acc{$n1}{'l'},$start,$acc{$n1}{'l'},$stop, "Sa$i"
						,$color, $width, $offset
						,$arrow,$opt{'sub_arrow_paral'},$opt{'sub_arrow_diag'},$opt{'sub_arrow_perp'}); 
				
				} else {
					&draw_line_horz_pieces($acc{$n1}{'l'},$start,$acc{$n1}{'l'},$stop, "Sa$i",
						$color, $width, $offset);
				}
				if ($opt{'sub_labelhit_on'}==1) {
					#print "WRITING LABEL FOR N2:$n2\n";
					my $label='';
					my $xblabel=$acc{$n1}{'xb'}+$accsub{$n2}{'acc'}{$n1}{'qmin'};
					if ( $opt{'sub_labelhit_col2'} ne '') {
						$label = $m[$i][ $opt{'sub_labelhit_col2'} ];
						$label=$1 if $opt{'sub_labelhit_col2_pattern'} && $label =~ /$opt{'sub_labelhit_col2_pattern'}/;
					} elsif ($opt{'sub_labelhit_col'} ne '') {
						
						$label = $m[$i][ $opt{'sub_labelhit_col'} ];
						$label=$1 if $opt{'sub_labelhit_col_pattern'}  && $label =~ /$opt{'sub_labelhit_col_pattern'}/;
						#print "  DRAW $label\n";
					}
					my ($xl, $yl)=&linexbp2xy($acc{$n1}{'l'},$start);
					my $tag=$canvas->createText($xl, $yl+$offset,
							-text=>$label, -fill => $opt{'sub_labelhit_color'},
							-font=>"Courier $opt{'sub_labelhit_size'}", -anchor=>"e");	
					$canvas->addtag('subl','withtag',$tag);
				}
				my $tagname = "Sa$i";  #can't go inside a sub!!!
				$canvas->addtag("sub","withtag",$tagname);
			}
			#print "ACCN2:(",defined $acc{$n2},")(", defined $accsub{$n1}, ")\n";
			if ($defined2 && defined $accsub{$n1} ) {
			#next if !defined2;
			############N1 is the SUB########################
				$color=$opt{'sub_color'};
				if ($m[$i][$mh{'scolor'}] ){
					$color=$m[$i][$mh{'scolor'}];
				} else {
					$color=$accsub{$n1}{'color'} if $accsub{$n1}{'color'};
				}
				my $offset= $accsub{$n1}{'acc'}{$n2}{'line'}*$opt{'sub_line_spacing'}+$opt{'sub_initoffset'};
				$offset = $m[$i][$mh{'sline'}] * $opt{'sub_line_spacing'}+$opt{'sub_initoffset'}if $m[$i][$mh{'sline'}] ne '';
				my $width= $opt{'sub_width'};
				my ($start,$stop)=($b2,$e2);
				#print "$n1:$n2 C$color W$width O$offset\n";
				($stop,$start)=($start,$stop) if $start> $stop;
				$accsub{$n1}{'acc'}{$n2}{'qmin_f'}=$start if $start < $accsub{$n1}{'acc'}{$n2}{'qmin_f'} ;
				$accsub{$n1}{'acc'}{$n2}{'qmax_f'}=$stop if $stop > $accsub{$n1}{'acc'}{$n2}{'qmax_f'};
				if ($opt{'sub_arrow_on'} ) {
					my $arrow='first';
					$arrow='last' if $b2<$e2;
					&draw_line_horz_pieces($acc{$n2}{'l'},$start,$acc{$n2}{'l'},$stop, "Sb$i",
						$color, $width, $offset, 
						,$arrow,$opt{'sub_arrow_paral'},$opt{'sub_arrow_diag'},$opt{'sub_arrow_perp'}); 
				
				} else {
					&draw_line_horz_pieces($acc{$n2}{'l'},$start,$acc{$n2}{'l'},$stop, "Sb$i",
					"$color", $width, $offset);
				}
				
				#print "Zoom...\n";
				if ( $opt{'sub_labelhit_on'}==1) {
					#print "WRITING LABEL FOR N1:$n1\n";
					my $label='';
					my $xblabel=$acc{$n2}{'xb'}+$accsub{$n1}{'acc'}{$n2}{'qmin'};
					if ( $opt{'sub_labelhit_col'} ne '' ) {
						$label = $m[$i][ $opt{'sub_labelhit_col'} ];
						$label=$1 if $opt{'sub_labelhit_col_pattern'} && $label=~/$opt{'sub_labelhit_col_pattern'}/;
					}
					my ($xl, $yl)=&linexbp2xy($acc{$n2}{'l'},$start);
					my $tag=$canvas->createText($xl, $yl+$offset,
							-text=>$label, 
							-font=>"Courier $opt{'sub_labelhit_size'}", -anchor=>"e");			
					$canvas->addtag('subl','withtag',$tag);
					$acc{$n2}{'labeldrawn'}=1;
				}
				my $tagname = "Sb$i";  #can't go inside a sub!!!
				$canvas->addtag("sub","withtag",$tagname);
			}

		}
		
		##############################################
		####DETERMINE IF INTER OR INTRA PICTURE ######
		my $pairtype = 'inter';
		#print "$n1:($acc{$n1}{'e'})$n2:($acc{$n2}{'e'})";
		$pairtype='intra' if (defined $acc{$n1} && defined $acc{$n2});	
		if ($opt{'pair_type_col'} ne '') {
			if ($opt{'pair_type_col2'} ne '') {
				#print "PAIR TYPE CHECK\n";
				my $text1=$m[$i][$opt{'pair_type_col'}];
				$text1=$1 if $opt{'pair_type_col_pattern'} && $text1=~/$opt{'pair_type_col_pattern'}/;
				my $text2=$m[$i][$opt{'pair_type_col2'}];
				$text2=$1 if $opt{'pair_type_col2_pattern'} && $text2=~/$opt{'pair_type_col2_pattern'}/;
				if ($text1 eq $text2) { $pairtype='intra'} else {$pairtype='inter'}
			
			} else {
				my $text1=$m[$i][$opt{'pair_type_col'}];
				$text1=$1 if $opt{'pair_type_col_pattern'} && $text1=~/$opt{'pair_type_col_pattern'}/;
				$pairtype='intra' if $text1=~/intra/i;
				$pairtype='inter' if $text1=~/inter/i;
			}
		
		}
		######FILTER THE PAIRWISE AND HIDE IF NECESSARY ######
		#####################################################
		#########determine other characteristics ###########
		my $color = $opt{"pair_$pairtype".'_color'};
		$color = $m[$i][$mh{'color'}] if $m[$i][$mh{'color'}];
		my $width = $opt{"pair_$pairtype".'_width'};
		$width = $m[$i][$mh{'width'}] if $m[$i][$mh{'width'}];
		my $offset = $opt{"pair_$pairtype".'_offset'};
		$offset =$m[$i][$mh{'offset'}] if $m[$i][$mh{'offset'}];		
		
		#### SKIP IF PAIRTYPE NOT TO BE DISPLAYED ######	
		next if !$opt{'pair_intra_on'} && $pairtype eq 'intra';
		next if !$opt{'pair_inter_on'} && $pairtype eq 'inter';
		####draw the lines #######
		my ($start,$stop)=($b1,$e1);
		($stop,$start)=($start,$stop) if $start> $stop;
		&draw_line_horz_pieces($acc{$n1}{'l'},$start,$acc{$n1}{'l'},$stop, "M$i",
				$color, $width, $offset)
					if $defined1;	
		($start,$stop)=($b2,$e2);
		($stop,$start)=($start,$stop) if $start> $stop;
		&draw_line_horz_pieces($acc{$n2}{'l'},$start,$acc{$n2}{'l'},$stop, "M$i",
				$color, $width, $offset)
					if $defined2;
		if ( 	($pairtype eq 'intra' && $opt{'pair_intra_line_on'}) 
				|| ($pairtype eq 'inter' && $opt{'pair_inter_line_on'} && $defined1 && $defined2 )
			) {			
			my $acc_e1=$acc{$n1}{'e'}+$acc{$n1}{'xb'}-1;
			my $acc_e2=$acc{$n2}{'e'}+$acc{$n2}{'xb'}-1;
			
				
			#print " $acc{$n1}{'xb'} <= $b1   && $b1 <= $acc_e1 && $acc{$n2}{'xb'} <= $b2   && $b2 <= $acc_e2\n";
			#print " $acc{$n1}{'xb'} <= $e1   && $e1 <= $acc_e1 && $acc{$n2}{'xb'} <= $e2   && $e2 <= $acc_e2\n";
			if ($acc{$n1}{'xb'} <= $b1   && $b1 <= $acc_e1 && $acc{$n2}{'xb'} <= $b2   && $b2 <= $acc_e2) {
				my ($b1x, $b1y)=&linexbp2xy($acc{$n1}{'l'},$b1);
				my ($b2x, $b2y)=&linexbp2xy($acc{$n2}{'l'},$b2);
				my $line;
				if ($b1y==$b2y) {
					#print "draw begin on same\n";
					$line= $canvas->createLine($b1x,$b1y, ($b1x+$b2x)/2, $b1y-0.66*$opt{'seq_line_spacing_wrap'}, 
								$b2x,$b2y,-width => 1, -fill => $color);
				} else {
					#print " begin drawing different\n";
					$line= $canvas->createLine($b1x,$b1y,$b2x,$b2y,-width => 1, -fill => $color);
				}
				$canvas->addtag("M$i", 'withtag',$line);	
			}
			if ( $acc{$n1}{'xb'} <= $e1   && $e1 <= $acc_e1 && $acc{$n2}{'xb'} <= $e2   && $e2 <= $acc_e2) {
				my ($e1x, $e1y)=&linexbp2xy($acc{$n1}{'l'},$e1);
				my ($e2x, $e2y)=&linexbp2xy($acc{$n2}{'l'},$e2);
				my $line;
				if ($e1y==$e2y) {
					$line= $canvas->createLine($e1x,$e1y, ($e1x+$e2x)/2, $e1y-0.66*$opt{'seq_line_spacing_wrap'}, 
							$e2x,$e2y,-width => 1, -fill => $color);
				} else {
					$line= $canvas->createLine($e1x,$e1y,$e2x,$e2y,-width => 1, -fill => $color);

				}		
				$canvas->addtag("M$i", 'withtag',$line);	
			}
		} elsif ($pairtype eq 'inter' && $opt{'pair_inter_line_on'} ) {
			my ($b1x,$b1y,$b2x,$b2y);
			if ($defined1 ) {
				#draw n1#
				($b1x, $b1y)=&linexbp2xy($acc{$n1}{'l'},$b1);
				($b2x, $b2y)=&linexbp2xy($acc{$n1}{'l'},$e1);
			} else {
				#draw n2#
				($b1x, $b1y)=&linexbp2xy($acc{$n2}{'l'},$b2);
				($b2x, $b2y)=&linexbp2xy($acc{$n2}{'l'},$e2);			
			}
			my $line= $canvas->createLine($b1x,$b1y, $b1x, $b1y-0.90*$opt{'seq_line_spacing_wrap'}, 
					$b2x,$b2y,-width => 1, -fill => $color);
			$canvas->addtag("M$i", 'withtag',$line);				
		}	
		
		my $tagname = "M$i";  #can't go inside a sub!!!
		if ($opt{'pair_level'} eq 'inter_over_intra' ) {
			if ($pairtype eq 'inter') {$canvas->raise($tagname);} else { $canvas->lower($tagname); }
		} elsif ($opt{'pair_level'} eq 'intra_over_inter' ) {
			if ($pairtype eq 'inter') {$canvas->lower($tagname);} else { $canvas->raise($tagname); }
		
		}
		$canvas->addtag("$pairtype","withtag",$tagname);
	}
	print "..$size\n"  if !$newopt{'quiet'};
}

sub draw_subject_labels { 
	foreach my $ac (@acc_order) {
		foreach my $s ( keys % {$acc{$ac}{'sub'}} ) {
			#print "$ac $s\n";
			my $ref=$m[$accsub{$s}{'acc'}{$ac}{'eghit'}];
			my $label='';
			#print " $s $a $accsub{$s}{'acc'}{$ac}{'qmin'} $accsub{$s}{'acc'}{$ac}{'qmax'}\n ";
			next if $accsub{$s}{'acc'}{$ac}{'qmin'} > $acc{$ac}{'e'} 
					|| $accsub{$s}{'acc'}{$ac}{'qmax'} < $acc{$ac}{'b'};
			if($opt{'sub_labelseq_on'} && $opt{'sub_labelseq_col'} ne '') {
				#print "TRYING TO LABEL SUBJECTS $$ref[4]\n";
				if ( $$ref[4] eq $s && $opt{'sub_labelseq_col2'} ne '' ) {
					$label = $$ref[ $opt{'sub_labelseq_col2'} ];
					$label=$1 if $opt{'sub_labelseq_col2_pattern'} && $label =~ /$opt{'sub_labelseq_col2_pattern'}/;
				} elsif ($opt{'sub_labelseq_col'} ne '') {

					$label = $$ref[ $opt{'sub_labelseq_col'} ];
					$label=$1 if $opt{'sub_labelseq_col_pattern'}  && $label =~ /$opt{'sub_labelseq_col_pattern'}/;
					#print "  DRAW $label\n";
				}
				my $xblabel= $accsub{$s}{'acc'}{$ac}{'qmin'}+ $acc{$ac}{'xb'}-$acc{$ac}{'b'};
				my $offset= $accsub{$s}{'acc'}{$ac}{'line'}*$opt{'sub_line_spacing'}+$opt{'sub_initoffset'};
				my ($xl, $yl)=&linexbp2xy($acc{$ac}{'l'},$xblabel);
				my $tag=$canvas->createText($xl- $opt{'sub_labelseq_offset'}, $yl+$offset,
						-text=>$label, -fill => $opt{'sub_labelseq_color'},
						-font=>"Courier $opt{'sub_labelseq_size'}", -anchor=>"e");			
				$canvas->addtag('subl','withtag',$tag);

			}
			if($opt{'sub_labelseqe_on'} && $opt{'sub_labelseqe_col'} ne '') {
				#print "TRYING TO LABEL SUBJECTS $$ref[4]\n";
				if ( $$ref[4] eq $s && $opt{'sub_labelseqe_col2'} ne '' ) {
					$label = $$ref[ $opt{'sub_labelseqe_col2'} ];
					$label=$1 if $opt{'sub_labelseqe_col2_pattern'} && $label =~ /$opt{'sub_labelseqe_col2_pattern'}/;
				} elsif ($opt{'sub_labelseqe_col'} ne '') {

					$label = $$ref[ $opt{'sub_labelseqe_col'} ];
					$label=$1 if $opt{'sub_labelseqe_col_pattern'}  && $label =~ /$opt{'sub_labelseqe_col_pattern'}/;
					#print "  DRAW $label\n";
				}
				my $xblabel= $accsub{$s}{'acc'}{$ac}{'qmax'}+ $acc{$ac}{'xb'}-$acc{$ac}{'b'};
				my $offset= $accsub{$s}{'acc'}{$ac}{'line'}*$opt{'sub_line_spacing'}+$opt{'sub_initoffset'};
				my ($xl, $yl)=&linexbp2xy($acc{$ac}{'l'},$xblabel);
				my $tag=$canvas->createText($xl- $opt{'sub_labelseqe_offset'}, $yl+$offset,
						-text=>$label, -fill => $opt{'sub_labelseqe_color'},
						-font=>"Courier $opt{'sub_labelseqe_size'}", -anchor=>"w");			
				$canvas->addtag('subl','withtag',$tag);

			}

		}	
	}
}

sub draw_extra {
	print "DRAWING EXTRA...\n"  if !$newopt{'quiet'};
	my $size=scalar(@e);
	my $skip_this=0;
	for (my $i=0; $i < $size;$i++) {
		print "..$i\n" if $i % 1000 ==0 && $i !=0 &&  !$newopt{'quiet'};
		my $ep=$e[$i];
		my ($n1,$b1,$e1,$color,$offset,$width,$arrow)=@{$ep}[0..6];
		next if !defined($acc{$n1});
		$skip_this=0;
		if ($opt{'filter1_col'} =~/^\d+$/) {
			$skip_this=1 if $$ep[$opt{'filterextra1_col'}]< $opt{'filterextra1_min'} && $opt{'filterextra1_min'} ne '';
			$skip_this=1 if $$ep[$opt{'filterextra1_col'}]> $opt{'filterextra1_max'} && $opt{'filterextra1_max'} ne '';
		}
		if ($opt{'filter2_col'} =~/^\d+$/) {
			$skip_this=1 if $$ep[$opt{'filterextra2_col'}]< $opt{'filterextra2_min'} && $opt{'filterextra2_min'} ne '';
			$skip_this=1 if $$ep[$opt{'filterextra2_col'}]> $opt{'filterextra2_max'} && $opt{'filterextra2_max'} ne '';
		}
		next if $skip_this;
		if ( 		($b1 < $acc{$n1}{'b'}  && $e1 < $acc{$n1}{'b'} ) 
				|| ($b1 > $acc{$n1}{'e'}  && $e1 > $acc{$n1}{'e'} ) ) {
						next;
		} else {
			$b1=$acc{$n1}{'b'} if $b1 < $acc{$n1}{'b'} ;
			$e1=$acc{$n1}{'e'} if $e1 > $acc{$n1}{'e'} ;
		}
		#print "$n1  $b1  $e1 C$color   W$width   O$offset ..../n";	
		$color=$opt{'extra_color'} if ! $color;
		$width=$opt{'extra_width'} if !$width;
		$offset=$opt{'extra_offset'} if $offset eq '';

		$b1= $acc{$n1}{'xb'}+$b1 -$acc{$n1}{'b'};
		$e1= $acc{$n1}{'xb'}+$e1 -$acc{$n1}{'b'};		
		my $tag='';
		if (defined $acc{$n1} ) {
			my ($start,$stop)=($b1,$e1);
			($stop,$start)=($start,$stop) if $start> $stop;
			if ($opt{'extra_arrow_on'} ) {
				if ($arrow eq 'F') { $arrow = 'last'
				} elsif ($arrow eq 'R') { $arrow = 'first'
				} else { $arrow = 'none'}
#				#print "DRAWING ARROW:$arrow\n";
				&draw_line_horz_pieces($acc{$n1}{'l'},$start,$acc{$n1}{'l'},$stop, "E$i", $color,$width,$offset
						,$arrow,$opt{'extra_arrow_para'},$opt{'extra_arrow_diag'},$opt{'extra_arrow_perp'}); #$thick *$opt{'pair_intra_width'});	
			} else {
				&draw_line_horz_pieces($acc{$n1}{'l'},$start,$acc{$n1}{'l'},$stop, "E$i", $color,$width,$offset); #$thick *$opt{'pair_intra_width'});	

			}						
			
			if ( $opt{'extra_label_on'}==1) {
				#print "DRAW LABEL\n";
				if ( $opt{'extra_label_test_col'} && $opt{'extra_label_test_pattern'} ) {
					#print "SKIPPING\n";
					next if $e[$i][ $opt{'extra_label_test_col'}] !~ /$opt{'extra_label_test_pattern'}/;
				}
				
				my $label='';
				my $xblabel=$start;
				if ( $opt{'extra_label_col'} ne '' ) {
					#print "FIGURE OUT LABEL $opt{'extra_label_col'}\n";
					$label = $e[$i][ $opt{'extra_label_col'} ];
					$label=$1 if $opt{'extra_label_col_pattern'} && $label=~/$opt{'extra_label_col_pattern'}/;
				}
				#print "DRAWING LABEL $label\n";
				my ($xl, $yl)=&linexbp2xy($acc{$n1}{'l'},$start);
				$tag=$canvas->createText($xl+$opt{'extra_label_offset'}, $yl+$offset,
						-text=>$label, 
						-font=>"Courier $opt{'extra_label_fontsize'}", -anchor=>"e",-fill=>$opt{'extra_label_color'});
				$canvas->addtag('exl','withtag',$tag);
				
			}			
			my $tagname = "E$i";  #can't go inside a sub!!!
			$canvas->addtag("ex","withtag",$tagname);
		}
	
	}
	print "..$size\n"  if !$newopt{'quiet'};
	
}

sub draw_graph {
	foreach my $numb (1,2) {
		my $array = "g$numb";
		my $size= scalar(@$array);
		next if $size==0;
		print "GRAPHING $numb..\n"  if !$newopt{'quiet'};
		my $p = \@$array;
		my $bpos;
		my $offset_slope=  -$opt{'graph_scale_height'}/($opt{"graph$numb"."_max"} - $opt{"graph$numb"."_min"});
		my $offset_intercept=$opt{'graph_scale_indent'} - $offset_slope*$opt{"graph$numb"."_min"};
		my @line;
		my $line_on=$opt{"graph$numb"."_line_on"};
		my $line_color=$opt{"graph$numb"."_line_color"};
		my $line_width=$opt{"graph$numb"."_line_width"};
		my $line_smooth=$opt{"graph$numb"."_line_smooth"};
		my $point_on=$opt{"graph$numb"."_point_on"};		
		my $point_shape=$opt{"graph$numb"."_point_shape"};		
		my $point_size=$opt{"graph$numb"."_point_size"};
		my $point_fill_color=$opt{"graph$numb"."_point_fill_color"};
		my $point_outline_color=$opt{"graph$numb"."_point_outline_color"};
		my $point_outline_width=$opt{"graph$numb"."_point_outline_width"};
		#print "$point_on  $point_shape $point_size $point_fill_color\n";
		my ($seq,$position,$val);
		my ($x1,$y1,$offset);
		for (my $i=0; $i < $size;$i++ ) {
			print "..$i\n" if $i % 1000 ==0 && $i !=0  && !$newopt{'quiet'};
			$seq=$$p[$i][0];
			$position=$$p[$i][1];
			$val=$$p[$i][2];
			next if !defined($acc{$seq});
			next if $position < $acc{$seq}{'b'};
			next if $position > $acc{$seq}{'e'};
			### calculate the graph positions###
			$bpos = $acc{$seq}{'xb'} + $position - $acc{$seq}{'b'};
			($x1,$y1) = &linexbp2xy($acc{$seq}{'l'},$bpos);
			$offset=$val*$offset_slope + $offset_intercept;
			$offset+=$y1;
			if ($line_on) {
				push @line, $x1,$offset if $val ne '' ;
				my ($x2,$y2);
				if ($i< $size-1) {
					($x2,$y2)=&linexbp2xy($acc{$seq}{'l'},$acc{$seq}{'xb'} + $$p[$i+1][1] - $acc{$seq}{'b'})
				}
				#print "$y1 ne $y2   $seq ne $$p[$i+1][0] $$p[$i+1][1] > $acc{$seq}{'e'}\n";
				if ( $i == $size-1 || $y1 ne $y2 || $seq ne $$p[$i+1][0] 
					|| $$p[$i+1][1] > $acc{$seq}{'e'} || $val eq '') {
				#draw lines
					#print "DRAWING LINES\n";
					 if (@line!=0) {
						my $tag=$canvas->createLine( @line,-fill=>$line_color,-width=>$line_width,-smooth=>$line_smooth,);
						$canvas->addtag('gl','withtag',$tag);
					}
					@line=();
				}
			}
			#print "$x1-$point_size,$offset-$point_size,$x1+$point_size, $offset+$point_size\n";
			#print " -fill => $point_fill_color,-outline=>$point_outline_color,-width=>$point_outline_width \n";
			if ($point_on && $val ne '') {
				$canvas->createOval($x1-$point_size, $offset-$point_size,
				$x1+$point_size, $offset+$point_size,
					 -fill => $point_fill_color,
					-outline=>$point_outline_color,
					-width=>$point_outline_width );
			}	
		}
		print "..$size\n"  if !$newopt{'quiet'};

	}
}


sub linexbp2xy {		
	my $line=shift;
	my $xbp= (shift) -1;
	my $x = ( $xbp % $opt{'canvas_bpwidth'})/$bp_per_pixel +$opt{'canvas_indent_left'};
	my $y = $l[$line]{'liney'} + int($xbp/$opt{'canvas_bpwidth'}) * $opt{'seq_line_spacing_wrap'};
	return $x, $y;
	
}

#the main line/rectangle drawing routine for seqs, pairs, subs, and extras
sub draw_line_horz_pieces {
	my ($l1,$xbp1,$l2,$xbp2, $tagname, $color, $width,$offset,$arrow,$a1,$a2,$a3) = @_;
	my ($x1,$y1)=&linexbp2xy($l1,$xbp1);
	my ($x2,$y2)=&linexbp2xy($l2,$xbp2);
	#print "=======>   ($x1,$y1)    ($x2,$y2)\n";
	my $line;
	for(my $i=$y1; $i<=$y2; $i+=$opt{'seq_line_spacing_wrap'}) {
		my ($xb,$xe)=($x1,$x2);
		$xb=$opt{'canvas_indent_left'} if $i > $y1;
		$xe=$opt{'canvas_indent_left'}+$canvas_width if $i< $y2;
		if (defined $arrow && $arrow ne 'none' ) {
			$line= $canvas->createLine($xb,$i+$offset,$xe,$i+$offset,-width => $width, 
					-arrow=>$arrow, -arrowshape=>[$a1,$a2,$a3],
					-fill => $color);
		} else {
			#$line= $canvas->createLine($xb,$i+$offset,$xe,$i+$offset,-width => $width, -fill => $color);
			$line= $canvas->createRectangle($xb,$i+$offset-$width/2, $xe, $i+$offset+$width/2 ,-fill=>$color, -outline => undef);
		}
		$canvas->addtag($tagname, 'withtag',$line);	
	}
}

sub arrange_subjects {
	###########################################################################
	############### sub categories ############################################
	if ($opt{'arrangesub'} =~/subscaleC#CHR_oo21/ ) {
		$opt{'arrangesub'}='subscaleC';
		$opt{'sub_scale_categoric_string'}='UK,NA,UL,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y';
		$opt{'sub_scale_max'}=0;
		$opt{'sub_scale_min'}=26;
		$opt{'sub_scale_lines'}=26;
		$opt{'sub_scale_step'}=-1;
		$opt{'sub_scale_col'}='0';
		$opt{'sub_scale_col_pattern'}='^([A-Z0-9]+)';
		$opt{'sub_scale_col2'}='4';
		$opt{'sub_scale_col2_pattern'}='^([A-Z0-9]+)';
		$opt{'sub_on'}=1;
		&find_column_options;
	} elsif ($opt{'arrangesub'} =~/subscaleN#ident(\d+)/ ) {
		$opt{'arrangesub'}='subscaleN';
		$opt{'sub_scale_categoric_string'}='';
		$opt{'sub_scale_max'}=1;
		$opt{'sub_scale_min'}=$1/100;
		$opt{'sub_scale_lines'}=20;
		$opt{'sub_scale_step'}=0.01;
		$opt{'sub_scale_col'}='FRACBPMATCH';
		$opt{'sub_scale_col_pattern'}='';
		$opt{'sub_scale_col2'}='';
		$opt{'sub_scale_col2_pattern'}='';
		&find_column_options;
	}
   ###########################################################################
   ############### main categories ###########################################
	if ($opt{'arrangesub'}=~ /oneperline/ ) {
		print "ARRANGING SUBJECTS ONE PER LINE\n" if !$newopt{'quiet'};
		for (my $i=0; $i< @m; $i++) {	$m[$i][$mh{'sline'}] ='';	}
		my $line=0;
		foreach my $ac (@acc_order) {
			#print "$ac\n";
			my @subnames = keys %{ $acc{$ac}{'sub'} } ; 
			if ($opt{'arrangesub_col'}) {
				print "SORTING ONEPERLINE ARRANGE USING COLUMN $opt{'arrangesub_col'}\n"  if !$newopt{'quiet'};
				foreach (@subnames) {
					#print "($accsub{$_}{'acc'}{$ac}{'eghit'})   $m[$acc{$ac}{'sub'}{$_}{'eghit'}][$opt{'arrangesub_col'}],$_] \n";
				
				}
				
				@subnames= map { $_->[1] }
							  sort {$a->[0] <=> $b->[0] }
							  map {    [  $m[$accsub{$_}{'acc'}{$ac}{'eghit'}][$opt{'arrangesub_col'}],$_  ]     }
							  @subnames;
				
			
			}
			@subnames=reverse @subnames if ($opt{'arrangesub_rev_on'}==1) ;
			$acc{$ac}{'subnum'}=0;
			foreach my $s (@subnames) {
				#print
				#foreach my $a (keys %{$accsub{$s}{'acc'}}) {
					$accsub{$s}{'acc'}{$ac}{'line'}=$acc{$ac}{'subnum'}++;
					$accsub{$s}{'acc'}{$ac}{'color'}=$opt{'sub_color'};
				#}
			}
		}
	} elsif ($opt{'arrangesub'} =~/stagger/ ) {
		print "ARRANGING SUBJECTS IN A STAGGER \n"  if !$newopt{'quiet'};
		my $spacer=$opt{'arrangesub_stagger_spacing'};
		for (my $i=0; $i< @m; $i++) {	$m[$i][$mh{'sline'}] ='';	}
		
		foreach my $ac (@acc_order) {
			#print "*staggering:$ac\n";
			my @tmpsub;
			my @subnames= keys % {$acc{$ac}{'sub'}};
			foreach my $s (@subnames) {
				my %tmp= ( 'name'=>$s, 'qmin'=>$accsub{$s}{'acc'}{$ac}{'qmin'}-$spacer,
							'qmax'=>$accsub{$s}{'acc'}{$ac}{'qmax'}+$spacer );
				push @tmpsub , \%tmp;
			}
			@tmpsub = sort { $$a{'qmin'} <=> $$b{'qmin'} } @tmpsub;
			my @endp=();
			#print "ACCS=>",scalar(@tmpsub),"\n";
			for (my $i=0; $i<@tmpsub; $i++) {
				$accsub{$tmpsub[$i]{'name'}}{'acc'}{$ac}{'color'}=$opt{'sub_color'};
				my $placed=0;
				for (my $j=0; $j<=@endp;$j++) {
					#print "$i:SEARCHING FOR EMPTY SLOT ",scalar(@endp),"\n";
					if ($endp[$j] < $tmpsub[$i]{'qmin'} ) {
						$accsub{$tmpsub[$i]{'name'}}{'acc'}{$ac}{'line'}=$j+1;
						$endp[$j]=$tmpsub[$i]{'qmax'};
						$placed=1;
						last;
					}
				}
				if ($placed==0) {
					#print "PUSHING $i\n";
					push @endp,$tmpsub[$i]{'qmax'};
					$accsub{$tmpsub[$i]{'name'}}{'acc'}{$ac}{'line'}=scalar(@endp);
				}
			}
			#my $pause=<STDIN>;
		}	
	}	
	
	if ($opt{'arrangesub'} =~/subscaleC/ ) {
		
		%subscaleC=();
		@subscaleC=split /[:, ]+/, $opt{'sub_scale_categoric_string'};
		for(my $i=0;$i<@subscaleC;$i++) { $subscaleC{$subscaleC[$i]}=$i; }
		$opt{'sub_scale_lines'}=scalar @subscaleC;
		if ($opt{'sub_scale_step'} ==1 ) {
			$opt{'sub_scale_max'}=$opt{'sub_scale_lines'};
			$opt{'sub_scale_min'}=0;
		} else {
			$opt{'sub_scale_step'}=-1 ;
			$opt{'sub_scale_min'}=$opt{'sub_scale_lines'};
			$opt{'sub_scale_max'}=0;	
		}
		
		my $upper_bound=$opt{'sub_scale_max'};
		my $lower_bound=$opt{'sub_scale_min'};
		my $scale= $opt{'sub_scale_lines'} / ($upper_bound-$lower_bound);
		my $col = $opt{'sub_scale_col'};
		my $col2 = $opt{'sub_scale_col2'};
		#print "U:$upper_bound  L:$lower_bound LN:$opt{'sub_scale_lines'} S:$scale  F:$field\n";
		for (my $i=0; $i< @m; $i++) {	$m[$i][$mh{'sline'}] ='';	}
		foreach my $s (@acc_ordersub) {
			my $line=0;
			my $cval='';
			foreach my $a (keys %{$accsub{$s}{'acc'}}) {
				my $eg=$accsub{$s}{'acc'}{$a}{'eghit'};
				if ( $col2 eq '' || $m[$eg][0] eq $s ) {
					$cval=$m[$eg][$col]	;
					##add if minside###
					$cval=$1 if $cval=~/$opt{'sub_scale_col_pattern'}/;
				} else {
					$cval=$m[$eg][$col2]	;
					$cval=$1 if $cval=~/$opt{'sub_scale_col2_pattern'}/;						
				}
				$line=$subscaleC{$cval};
				$accsub{$s}{'acc'}{$a}{'line'}= ($upper_bound - $line)*$scale;
				if ($line) {
					$accsub{$s}{'acc'}{$a}{'color'}='purple'; #$opt{'sub_color'};
				} else {
					$accsub{$s}{'acc'}{$a}{'color'}='darkgrey';
				}
			}
		}
	} elsif ($opt{'arrangesub'} =~/subscaleN/ ) {
		my $upper_bound=$opt{'sub_scale_max'};
		my $lower_bound=$opt{'sub_scale_min'};
		my $scale= $opt{'sub_scale_lines'} / ($upper_bound-$lower_bound);
		my $col = $opt{'sub_scale_col'};
		#print "U:$upper_bound  L:$lower_bound LN:$opt{'sub_scale_lines'} S:$scale  F:$field\n";
		for (my $i=0; $i< @m; $i++) {
			#print "hello";
			$m[$i][$mh{'sline'}] = ($upper_bound - $m[$i][$col])*$scale;
		}
	}
	
}

sub colorsub {
	#$opt{'colorsub'}='hitrandom';
	if ($opt{'colorsub'}=~/NO CHANGE/) {
	
	} elsif ($opt{'colorsub'}=~/RESET/) {
		my $col=$mh{'scolor'};
		foreach (@m) { $$_[$col]='';}
	} elsif ($opt{'colorsub'}=~/hitrandom/ ) {
		my @colors = qw(#cd0d32083208 #86e4ce2ceb01 orange #0068cd0d0000  purple #cd0daa2f7d15 #a207cd0d5a04 #9b1ecd0d9b05 brown);
		my $col=$mh{'scolor'};
		foreach (@m) { 
			$$_[$col]=$colors[0];
			push(@colors,shift(@colors)); #first shall be last
		}
	} elsif ($opt{'colorsub'}=~/hitconditional/ ) {
		 print "COLORING SUBJECTS CONDITIONALLY...\n"  if !$newopt{'quiet'};
		 my $col=$mh{'scolor'};
		 my @colors=split / *[,;] */, $opt{'colorsub_hitcond_tests'};
		 my @tests=();
		 foreach (@colors) {
		 	my @c=split / *if */;
		 	$c[0] =~s/ +//mg;
		 	if (@c !=2) {
		 		warn "$c[0] $c[1] bad color extract from ($_)!\n";
		 		return;
		 	}
		 	next if $c[0] !~/^[A-Za-z0-9]+$/;
		 	push @tests , \@c;
		 }
		foreach my $row (@m) { 
			foreach my $c (@tests) {
				if ( eval( '$$row[$opt{"colorsub_hitcond_col"}]' ." $$c[1]" ) ) {
		 			$$row[$col]=$$c[0];
				}
			}
		}
	} elsif ( $opt{'colorsub'}=~/seqrandom/ ) {
		print "COLORING...\n"  if !$newopt{'quiet'};
		my @colors = qw(#cd0d32083208 orange #0068cd0d0000 #86e4ce2ceb01 purple #cd0daa2f7d15 #a207cd0d5a04 #9b1ecd0d9b05 brown);
		foreach my $ac (@acc_order) {
			my @subnames= keys % { $acc{$ac}{'sub'} };
			foreach my $s (@subnames) {
				$accsub{$s}{'color'}=$colors[0];
				push(@colors,shift(@colors)); #first shall be last
			}
		}
	} 
}




sub show_calculations {
	#print "showKEYS: ",join(" ",keys %acc),"\n";

	for (my $i=0; $i<@m; $i++) {
		##########FILTERING TO OCCUR#################################
		#############################################################
		$m[$i][$mh{'hide'}]=0 if $opt{'pfilter_reset'} == 1; #checkbox
		if ($opt{'filterpre1_col'} =~/^\d+$/) {
			$m[$i][$mh{'hide'}]=1 if $m[$i][$opt{'filterpre1_col'}]< $opt{'filterpre1_min'} && $opt{'filterpre1_min'} ne '';
			$m[$i][$mh{'hide'}]=1 if $m[$i][$opt{'filterpre1_col'}]> $opt{'filterpre1_max'} && $opt{'filterpre1_max'} ne '';
		}
		if ($opt{'filterpre2_col'} =~/^\d+$/) {
			$m[$i][$mh{'hide'}]=1 if $m[$i][$opt{'filterpre2_col'}]< $opt{'filterpre2_min'} && $opt{'filterpre2_min'} ne '';
			$m[$i][$mh{'hide'}]=1 if $m[$i][$opt{'filterpre2_col'}]> $opt{'filterpre2_max'} && $opt{'filterpre2_max'} ne '';
		}
		next if $m[$i][$mh{'hide'}]==1;
		############################################################
		my $n1= $m[$i][0];
		my $n2= $m[$i][4];
		##############QUERY########################
		#print "$m[$i][0]  $m[$i][4]\n";
		if ( defined $acc{ $n1 }  ) { 
			$acc{$n1}{'len'}=$m[$i][3];
			$acc{$n1}{'sub'}{$n2}++ if defined $accsub{ $n2};
			
			#print "show:seq1 $m[$i][0] => $acc{$m[$i][0]}{'len'}\n"; 
		} 
		if ( defined $acc{$n2}  ) { 
			#print "show:seq2 $m[$i][4] => $acc{$m[$i][4]}{'len'}\n"; 
			$acc{$n2}{'len'}=$m[$i][7]; 
			$acc{$n2}{'sub'}{$n1}++ if defined $accsub{ $n1};
		}
		###########SUBJECT ###########################
		if ( defined $accsub{ $n1 }  && defined $acc{$n2} ) { 
			$accsub{$n1}{'len'}=$m[$i][3];
			my ($b,$e)=($m[$i][5],$m[$i][6] );
			($b,$e)=($e,$b) if $b> $e;
			if (! defined $accsub{$n1}{'acc'}{$n2}{'qmin'}) {
				$accsub{$n1}{'acc'}{$n2}{'qmin_f'} = 999999999;
				$accsub{$n1}{'acc'}{$n2}{'qmax_f'} = -10;
				$accsub{$n1}{'acc'}{$n2}{'qmin'} = 999999999 ;
				$accsub{$n1}{'acc'}{$n2}{'qmax'} =-10 ;
				$accsub{$n1}{'acc'}{$n2}{'eghit'} = $i;
			}
			if ( $b<$accsub{$n1}{'acc'}{$n2}{'qmin'} ) { $accsub{$n1}{'acc'}{$n2}{'qmin'}= $b ;}
			if ($e>$accsub{$n1}{'acc'}{$n2}{'qmax'} ) { $accsub{$n1}{'acc'}{$n2}{'qmax'}= $e ; }
			
			#print "show:seq1 $m[$i][0] => $acc{$m[$i][0]}{'len'}\n"; 
		} 
		
		if ( defined $accsub{ $n2 }  && defined $acc{$n1} ) { 
			my ($b,$e) = ( $m[$i][1],$m[$i][2] );
			$accsub{$n2}{'len'}=$m[$i][7]; 
			$accsub{$n2}{'desc'}=$m[$i][$mh{'DEFN2'}]; 
			if (! defined $accsub{$n2}{'acc'}{$n1}{'qmin'} ) {
				$accsub{$n2}{'acc'}{$n1}{'qmin_f'}=99999999999999;
				$accsub{$n2}{'acc'}{$n1}{'qmax_f'}=-10;
				$accsub{$n2}{'acc'}{$n1}{'qmin'}= 99999999999999 ;
				$accsub{$n2}{'acc'}{$n1}{'qmax'}=-10;
				$accsub{$n2}{'acc'}{$n1}{'eghit'} = $i;
			}
			if ( $b<$accsub{$n2}{'acc'}{$n1}{'qmin'} ) { $accsub{$n2}{'acc'}{$n1}{'qmin'}= $b ; }
			if ( $e>$accsub{$n2}{'acc'}{$n1}{'qmax'} ) { $accsub{$n2}{'acc'}{$n1}{'qmax'}= $e ; }	
		}
		
	}
	#print "showKEYS: ",join(" ",keys %acc),"\n";
	foreach my $s (keys %accsub) {
		foreach my $a (keys %{	$accsub{$s}{'acc'} } ) {
			#print "$s => $a \n";
			#print "$s:$a $accsub{$s}{'acc'}{$a}{'qmin'} = $acc{$a}{'b'}\n";
			if ( defined $acc{$a}{'b'} && $accsub{$s}{'acc'}{$a}{'qmin'} < $acc{$a}{'b'}
						&& $accsub{$s}{'acc'}{$a}{'qmax'} >= $acc{$a}{'b'}  ) {
				  $accsub{$s}{'acc'}{$a}{'qmin'} = $acc{$a}{'b'};
					#print "$s:$a $accsub{$s}{'acc'}{$a}{'qmin'} = $acc{$a}{'b'}\n";
			}
			#print "$s:$a $accsub{$s}{'acc'}{$a}{'qmax'} = $acc{$a}{'e'}\n";
			if (defined $acc{$a}{'e'} && $accsub{$s}{'acc'}{$a}{'qmax'} > $acc{$a}{'e'} 
					&& $accsub{$s}{'acc'}{$a}{'qmin'} <= $acc{$a}{'e'} ) {
				  $accsub{$s}{'acc'}{$a}{'qmax'} = $acc{$a}{'e'};
			#print "$s:$a $accsub{$s}{'acc'}{$a}{'qmax'} = $acc{$a}{'e'}\n";
			}
		}
	}
	#my $pause=<STDIN>;

}


sub fitlongestline { 
				$opt{'canvas_bpwidth'}=$widest_line+100;
				&redraw;
}  
#################################################################################
################ FILE HANDLING AND DATA PARSING #################################
#################################################################################


sub align_update {
	if (!defined $mh{'name1'}) {
		###INITALIZE IF $MH NOT DEFINED
		%mh=(name1=>0,begin1=>1,end1=>2,len1=>3,name2=>4,begin2=>5,end2=>6,len2=>7) ;
		@mheader=qw(name1 begin1 end1 len1 name2 begin2 end2 len2);
	}
	my $col=scalar(@mheader);
	######THIS IS WHERE TO ADD ADDITIONAL REQUIRED FIELDS#
	$mh{'color'}=$col++   if !defined $mh{'color'};
	$mheader[$mh{'color'}]='color';
	$mh{'offset'}=$col++  if !defined $mh{'offset'};
	$mheader[$mh{'offset'}]='offset';
	$mh{'width'}=$col++   if !defined $mh{'width'};
	$mheader[$mh{'width'}]='width';
	$mh{'display'}=$col++ if !defined $mh{'display'};
	$mheader[$mh{'display'}]='display';	
	$mh{'sline'}=$col++ if !defined $mh{'sline'};
	$mheader[$mh{'sline'}] = 'sline';
	$mh{'scolor'}=$col++ if !defined $mh{'scolor'};
	$mheader[$mh{'scolor'}] = 'scolor';	
	$mh{'hide'}=$col++ if !defined $mh{'hide'};
	$mheader[$mh{'hide'}] = 'hide';	
	

	return if $opt{'align'} eq '';
	my @files;
	(@files=split":",$opt{'align'}) || ( $files[0]=$opt{'align'} );	
	foreach my $f (@files) {
		print "LOADING ALIGN $f...\n"  if !$newopt{'quiet'};
		if (!open (ALIGN, $f) ) {
			 $opt{'align'}.='(bad)';
			 die "Can't open $f!\n";
		}
		my $line=<ALIGN>;
		$line=~s/\r\n/\n/;
		close ALIGN;
		my ($tmp, @head_file);
		if ($line=~ /\t/) {
			#some form of tab delimited#
			chomp $line;
			@head_file=split "\t", $line;
			$tmp = &import_blastout($f);
			
		} else {
			#traditional mirrorpeat output#
			@head_file=();  #no extra columns;
			$tmp = &import_miropeat($f);
		}
		for (my $i=8; $i<@head_file; $i++) {
			if (!defined $mh{$head_file[$i]} ) {
				my $l=scalar(@mheader);
				$mheader[$l]=$head_file[$i];
				$mh{$head_file[$i]}=$l;
			}
		}
		#################################
		###this is really slow I think###
		foreach (@{$tmp}) {
			my @c=@{$_};
			next if $c[0] eq '';
			my @cc;
			$cc[0]=$c[0]; $cc[1]=$c[1]; $cc[2]=$c[2]; $cc[3]=$c[3];
			$cc[4]=$c[4]; $cc[5]=$c[5]; $cc[6]=$c[6]; $cc[7]=$c[7];
			for (my $i=8; $i<@head_file; $i++) {
				$cc[$mh{$head_file[$i]}]=$c[$i];
			}
			
			push @m , \@cc;
		}
	} #files loop
	print "   ===> ", scalar(@m)," total pairwise comparisons to display\n"  if !$newopt{'quiet'};

}

sub show_update {
	#####show (NONE, blank current, file
	#print "FORCING...($opt{'showseq'})\n";
	#assume a sequence name if file path doesn't exist or if input contains a colon or comma#
	if ( $opt{'showseq'} eq 'ALL') {
		print "BUILDING show data for ALL\n"  if !$newopt{'quiet'};
		%acc=();
		@acc_order=();
		foreach ( @m ) {
			$acc{$$_[0]}{'len'}="0";
			$acc{$$_[4]}{'len'}="0" if !$opt{'showqueryonly'};
		}
		@acc_order= sort keys %acc;
	} elsif ($opt{'showseq'} =~ /[:,]/ ||  ! -e $opt{'showseq'} )  {
		warn "ASSUMMING: $opt{'showseq'} is a sequence name rather than file\n" if $opt{'showseq'} !~ /[:,]/;
		my $show = ":$opt{'showseq'}";
		$show=~s/\s$//;
		$show.=':' if $show !~ /:$/;
		%acc=();
		@acc_order=();
		foreach (split /:/	,$opt{'showseq'}) {
			print "ACC_ORDER ($_)\n"  if !$newopt{'quiet'};
			if (/,/) {
				warn "ERROR: Bad format for -showseq file commas ($_)!\n" if !/^\S+,\d*,\d+,\d+$/ && !/^\S+,\d+$/;
				my @c=split /,/;
				if ($c[1] eq '') { $acc{$c[0]}{'len'}=0 } else { $acc{$c[0]}{'len'}= $c[1] }
				$acc{$c[0]}{'b'}=$c[2];
				$acc{$c[0]}{'e'}=$c[3];
				push @acc_order,$c[0];
				
			} else {
				$acc{$_}{'len'}='0';
				push @acc_order,$_;
			}
		}
		
	} elsif ($opt{'showseq'} ne '') {
		%acc=();
		@acc_order=();
		open(IN,$opt{'showseq'} )|| die "Can't open  assummed showseq file --showseq ($opt{'showseq'})\n";
		print "LOADING show $opt{'showseq'}...\n"  if !$newopt{'quiet'};
		my $header=<IN>;
		while (<IN>) {
			s/\r\n/\n/;
			chomp; chomp;
			my @c = split "\t";
			next if $c[0] eq '';
			$acc{$c[0]}{'len'}="0";
			if ($c[1] =~/^\d+$/ ) {
				$acc{$c[0]}{'len'}=$c[1] ;
			} elsif ($c[2] =~ /\w/ ) {
				die "SHOWSEQ FILE ($opt{'showseq'}) ERROR!: sequence ($c[0]) subseq begin ($c[2]) is not an integer position\n";
			}
			if ($c[2] =~/^\d+$/ ) {
				$acc{$c[0]}{'b'}=$c[2] ;
				&warnNpause("SHOWSEQ FILE ($opt{'showseq'}) WARNING!:  sequence ($c[0]) subseq begin ($c[2]) more than the length ($c[1])!") 	if $c[2] > $c[1];
			} elsif ($c[2] =~ /\w/ ) {
				die "SHOWSEQ FILE ($opt{'showseq'}) ERROR!:  sequence ($c[0]) subseq begin ($c[2]) is not an integer position\n";
			}
			if ($c[2] =~/^\d+$/ ) {
				$acc{$c[0]}{'e'}=$c[3] ;
				&warnNpause("SHOWSEQ FILE ($opt{'showseq'}) WARNING!:  sequence ($c[0]) subseq end($c[3]) more than the length ($c[1])!\n") if $c[3] > $c[1];
			} elsif ($c[2] =~ /\w/ ) {
				die "SHOWSEQ FILE ($opt{'showseq'}) ERROR!:  sequence ($c[0]) subseq end  ($c[3]) is not an integer position!\n";
			}
		
			push @acc_order, $c[0];
		}
		close IN;
	} else {
		#return if  keys (%acc);
		###This is the first time!!!
		die "IMPROPER show $opt{'showseq'}...\n";
	}	
	print "   ===> ",scalar (@acc_order),  " total sequences to display\n"  if !$newopt{'quiet'};
}

sub show_update_subject {
	#####show (NONE, blank current, file
	print "FORCINGSUB...$opt{'showsub'} \n"  if !$newopt{'quiet'};
	if ($opt{'showsub'} =~ /\:/ ) {
		my $show = $opt{'showsub'};
		%accsub=();
		@acc_ordersub=();
		for (my $i=0; $i < @m; $i++) {
			if ( defined $acc{$m[$i][4]} && $show =~/$m[$i][0]/ ) {
				$accsub{$m[$i][0]}{'len'}="0";
			}
			if ( defined $acc{$m[$i][0]} && $show =~ /$m[$i][4]/ ) {
				$accsub{$m[$i][4]}{'len'}="0";
			}
			@acc_ordersub=keys %accsub;
			if ( $show =~/^SORT/ ) {
				@acc_ordersub=sort @acc_order;
			}
		}
	} elsif ( $opt{'showsub'} eq 'ALL') {
		print "BUILDING show FILE BY ALL\n";
		%accsub=();
		@acc_ordersub=();
		for (my $i=0; $i < @m; $i++) {
			#print "$m[$i][0]:::$m[$i][3]\n";
			$accsub{$m[$i][0]}{'len'}="0" if defined $acc{$m[$i][4]};
			$accsub{$m[$i][4]}{'len'}="0" if defined $acc{$m[$i][0]};
		}
		@acc_ordersub= sort keys %accsub;
		
	} elsif ($opt{'showsub'} ne '') {
		%accsub=();
		@acc_order=();
		open(IN,$opt{'showsub'} )|| die "Can't open accession show display --showseq ($opt{'showseq'}\n";
		print "    ===>loading file ($opt{'showsub'})...\n"  if !$newopt{'quiet'};
		my $header=<IN>;
		while (<IN>) {
			s/\r\n/\n/;
			chomp; chomp;
			my @c = split "\t";
			next if $c[0] eq '';
			#print "$c[0]  => $c[1]\n";
			$accsub{$c[0]}{'len'}="0";
			push @acc_ordersub, $c[0];
		}
		close IN;
	} else {
		return if  keys (%accsub);
		die "Unable to choose subjects to choose ALL use -showsub ALL\n";
	}	
	print "   ===> ",scalar(@acc_ordersub)," total sub sequences to display\n"  if !$newopt{'quiet'};
}
sub extra_update {
	if (!keys %eh) {
		#must have these seq, begin end#
		%eh=(seq=>0,begin=>1,end=>2,color=>3,offset=>4,width=>5, orient=>6);	
		@eheader=qw(seq begin end color offset width orient);		
	}
	#ADD NEW CALCULATED VARIABLES AS BELOW WITH ORIENT#
	my $col = scalar(@eheader);
	$eh{'orient'}=$col++   if !defined $eh{'orient'};
	$eheader[$eh{'orient'}]='orient';		
	return if $opt{'extra'} eq '';
	my @files;
	@files=split ":",$opt{'extra'} || ($files[0]=$opt{'extra'});
	foreach my $f (@files) {
		print "LOADING EXTRA FILE ($f)\n"  if !$newopt{'quiet'};
		if (! open (IN,"$f")  ) {
			$opt{'extra'}.="(bad name)";
			die "Can't open extra file ($f) [$!]\n";
		}
		my $header=<IN>;
		$header=~s/\r\n/\n/;
		chomp $header;
		my @head = split "\t",$header;
		for (my $i=3; $i<@head; $i++) {
			if (!defined $eh{$head[$i]}) {
				my $l= scalar(@eheader);
				$eheader[$l]=$head[$i];
				$eh{$head[$i]}=$l;
			}
			#print "$head[$i] ===> "
			
		}
		#Sprint "PAUSE\n"; my $pause=<STDIN>;
		my $test_show='';
		$test_show=":$opt{'showseq'}" if  $opt{'minload'} && $opt{'showseq'}=~/:$/;
		print "$test_show\n";
		while (<IN>) {
			s/\r\n/\n/;
			chomp; chomp;
			my @c = split "\t";
			next if $c[0] eq '' || $c[1] eq '';
			if ($test_show) {
				next if $test_show !~ /\:$c[0]\:/;
			
			}
			####
			my @cc=();
			$cc[0]=$c[0];$cc[1]=$c[1];$cc[2]=$c[2];
			for (my $i=3; $i<@head; $i++) {
				if ($head[$i] eq 'orient' && $c[$i] ne '') {
					$c[$i]=uc$c[$i];
					if ($c[$i] eq 'R' || $c[$i] eq 'F') {
					} elsif ($c[$i]=~/^PLUS|^POSITIVE|^FORWARD|^1|^\+/ ) { #gets 1 or +1
						$c[$i]='F';
					} elsif ($c[$i]=~/^MINUS|NEGATIVE|REVERSE|^\-/ ) {  #this gets -1 too
						$c[$i]='R';
					} else {
						&warnNpause("Unknown designation for orientation ($c[$i])!\n" );
					}				
				}
				$cc[$eh{$head[$i]}]=$c[$i];
			}
			push @e, \@cc;		
		}
		close IN;
	}
	print "   ===> ", scalar(@e)," total extra sequence features to display\n"  if !$newopt{'quiet'};

}	
sub graph_update {
	#ADD NEW DATA TO GRAPH#
	#	@g1 and @g2  #
	print "UPDATING GRAPH DATA\n"  if !$newopt{'quiet'};
	foreach my $g ( 1, 2 ) {
		next if $opt{"graph$g"} eq '';
		my @files;
		@files=split /:/, $opt{"graph$g"};
		my $array = "g$g";
		my $p = \@$array;
		foreach my $f (@files) {
			print "LOADING GRAPH $f\n"  if !$newopt{'quiet'};
			if (! open (IN,"$f")  ) {
				$opt{"graph$g"}.="(bad)";
				die "Can't read graph file ($f) [$!]\n";
			}
			my $header=<IN>;
			my $last_pos=0;
			my $last_seq='';
			while (<IN>) {
				s/\r\n/\n/;
				chomp;
				my @c = split /\t/;
				next if $c[0] eq '' || $c[1] eq '';
				my @cc=($c[0],$c[1],$c[2]);
				warn "Bad non-numerical format in ($f)  for $cc[1]\n" if $cc[1] !~/^[-0-9.]+$/;
				warn "Bad non-numerical format in ($f)  for $cc[2]\n" if $cc[2] !~/^[-0-9.]+$/ && $cc[2] ne '';
				if ($cc[1] < $last_pos && $last_seq eq $cc[0]  ) {
					warn "graph data not in sequential order\n";
				}
				push @$p, \@cc;
				$last_pos=$cc[1];
				$last_seq=$cc[0];
			}
		}
		print "   ===> ",scalar(@$p)," total points to display for graph$g\n"  if !$newopt{'quiet'};
		$opt{"graph$g"}='';
	}

}	


sub import_blastout {
	my $file = shift;
	my $showr=":$opt{'showseq'}";
	if ($opt{'showseq'} !~/\:/) {
		if ($opt{'showseq'} ne 'ALL') {
			$showr=':';
			print " PREPROCESSING show\n"  if !$newopt{'quiet'};
			open (show, "$opt{'showseq'}") || die "Can't open show file ($opt{'showseq'}\n";
			my $header=<show>;
			while (<show>) {
				s/\r\n/\n/;
				chomp;
				my @c=split "\t";
				#print "$c[0]\n";
				$showr.="$c[0]:";
			}
		}
	}
	print "$showr\n";
	open (IN, "$file")  || die "Can't open blast $file \n";
	my $header=<IN>;
	my @array;
	my $i;
	while (<IN>) {
		s/\r\n/\n/;
		chomp;
		my @c=split "\t";
		if ($opt{'showseq'} ne 'ALL') {
			next if $c[1]!~/^\d+$/;  #check for proper format
			next if $opt{'minload'} && $showr!~/$c[0]:/ && $showr !~ /$c[4]:/;
		}
		@{$array[$i] }=  @c;
		#print $array[$i][0], " ";
		$i++;
	}
	close IN;
	print "DONE\n";
	return \@array;
}


sub import_miropeat { 
	my $file = $_[0];
	open (IN, "$file")  || die "Can't open mirorepeat output $file \n";
	my $line ="";
	$line =<IN> until ( $line=~/^\./ );;
	$line=~s/^.//;
	return undef if $line=~/not find/;
	my $i =0;
	my @array=();
	while ($line !~ /^Graphic/ ) {
		$line=~s/\r\n/\n/;
		chomp $line;
		@{$array[$i] }= split " ",$line;
		$line=<IN>;
		$i++;
	}
	close IN;
	return \@array;
}


sub show_alignment {
	my $i=shift;
	my $orient = 'F';
	if ($m[$i][5]>$m[$i][6]) { $orient ='R'}
	#print "$i)))$m[$i][$opt{'alignment_col'}]\n";
	my $text;
	$$text="Prealigned sequences must be included as columns\nin the align file for this option to work!\nThese sequences must contain indel dashes\nas the alignment is  not recalculated.";
	if ( $opt{'alignment_col'}!=0 && $opt{'alignment_col2'} !=0 ) {
		$text=&alignment_format( -seq1=>$m[$i][$opt{'alignment_col'}], -seq2=>$m[$i][$opt{'alignment_col2'}],
					-name1=>$m[$i][0], -begin1=>$m[$i][1], -end1=>$m[$i][2],
					 -name2=>$m[$i][4], -begin2=>$m[$i][5],-end2=>$m[$i][6],
					 -orient2=>'F', -width=>$opt{'alignment_wrap'}); 
	}
	#print "$$text\n";
	&export_text( $text, "Formatted Alignment M$i" );         
}


sub save_parasight_table {
 	my $name = $_[0];
	$name=~ s/\.ps[a-z]?$//;  
	print "SAVING BASENAME ($name)\n";

	#####SAVE .psa ####
	if ( !open(OUT, ">$name.psa")   ) {
		print "WARNING: Can't Save file ($name.psa)\n";
		 return;
	}
	print OUT join("\t",@mheader),"\n";
	for (my $i=0; $i< @m; $i++) {
		print OUT join ( "\t",@{$m[$i]} ),"\n";
	}
	close OUT;
	#####SAVE .pse ####
	if ( !open(OUT, ">$name.pse")   ) {
		print "WARNING: Can't Save file ($name.pse)\n";
		 return;
	}
	print OUT join("\t",@eheader),"\n";
	for (my $i=0; $i< @e; $i++) {
		print OUT join ( "\t",@{$e[$i]} ),"\n";
	}
	close OUT;
	####SAVE .psg if it exists####
	my $max=@g1;
	$max=@g2 if @g2 > $max;
	if ($max > 0 ) {
		if (  !open (OUT, ">$name.psg")  )  {
			print "WARNING: Can't Save file ($name.psg)\n";
			return;
		}
		print OUT "g1seq\tg1point\tg1value\tg2seq\tg2point\tg2value\n";
		for (my $i=0; $i < $max; $i++) {
			if (defined $g1[$i]) {
				print OUT "$g1[$i][0]\t$g1[$i][1]\t$g1[$i][2]";
			} else {
				print OUT "\t\t\t";
			}
			if (defined $g2[$i]) {
				print OUT "\t$g2[$i][0]\t$g2[$i][1]\t$g2[$i][2]\n";
			} else {
				print OUT "\t\t\t\n";
			}
		}
		close OUT;
	}
	##### SAVE.pso ######
	###removing nstore###
	if (  !open (OUT, ">$name.pso")  )  {
		print "WARNING: Can't save parsight option file ($name.pso)\n";
		return;
	}
	foreach (keys %acc) {
		print OUT "#ACC||||$_||||$acc{$_}\n";
	}
	foreach (@acc_order) {
		print OUT "#ORDER||||$_\n";
 	}
	foreach (sort keys %opt) {
		print OUT "#OPT||||$_||||$opt{$_}\n";
	}
	close OUT;
}





sub load_option_template {
	my $f=shift;
	if (!open (OPTION, "$f")) {
		print "WARNING: Can't read option file ($f)\n";
		return;
	}
	my $line=<OPTION>;
	while (<OPTION>) {
		next if !/=>/;
		s/\r\n/\n/;chomp; chomp;
		my @c=split / *=> */;
		$c[1]='' if /=>$/;
		if (@c !=2) {
			print "WARNING: IMPORPER OPTION FORMAT:($_)\n";
			return;
		}
		next if $_=~/^#/;
		next if $_ eq 'showseq' || $_ eq 'in' || $_ eq 'align' || $_ eq 'extra' || $_ eq 'showsub' || $_ eq 'graph1' || $_ eq 'graph2'; 
		$opt{$c[0]}=$c[1];
	}
	close OPTION;
	print "TEMPLATE ($f) LOADED.\n  (Use RESHOW AND REDRAW to effect changes)!\n" if !$opt{'quiet'};

}

sub save_option_template {
	my $f=shift;
	$f.='.pst' if ($f!~/\.pst$/);
	if (!open (OPTION, ">$f")) {
		print "WARNING: Can't write option file ($f)\n";
		return;
	}
	print OPTION "OPTION FILE   option=>value\n";
	foreach  (sort keys %opt) {
		next if $_ eq 'showseq' || $_ eq 'in' || $_ eq 'align' || $_ eq 'extra' || $_ eq 'showsub' || $_ eq 'graph1' || $_ eq 'graph2'; 
		print OPTION "$_=>$opt{$_}\n";
		print OPTION "# $optdesc{$_}\n" if defined $optdesc{$_} && $opt{'template_desc_on'};
	}
	close OPTION;
	print "OPTION TEMPLATE ($f) CREATED AND SAVED!\n" if !$opt{'quiet'};;

}

sub load_parasight_table {
	my $name=$_[0];
	print "LOADING PARASIGHT FILE:$name\n";
	$name=~ s/\.ps[aeog]?$//;
	#system "ls";
	foreach ('.psa','.pso','.pse') {
		die "\nERROR: Could not locate Parasight save file ($name$_)!\n" if  ! -e "$name$_";
	}
	##### LOAD.psa ####
	if ( !open(IN, "$name.psa")   ) {
		print "WARNING: Can't Load file ($name.psa)\n";
		 return;
	}
	my $line=<IN>; $line=~ s/\r\n/\n/; chomp $line;
	@m=();%mh=();@mheader=();
	@mheader = split "\t", ($line);
	for(my $i=0;$i<@mheader;$i++) {$mh{$mheader[$i]}=$i;}
	my $count=0;
	while (<IN>) {
		s/\r\n/\n/; chomp;
		next if $_ eq '';
		my @c = split "\t";
		#print "loading $c[0]\n";
		$m[$count++]=\@c;
	}
	close IN;
	
	##### LOAD .pse ####
	if ( !open(IN, "$name.pse")   ) {
		print "WARNING: Cannot load parsight extra file ($name.pse)\n";
		return;
	}
	my $line=<IN>;
	chomp $line;
	@e=(); %eh=(); @eheader=();
	@eheader = split "\t", ($line);
	for(my $i=0;$i<@eheader;$i++) {$eh{$eheader[$i]}=$i;}
	my $count=0;
	while (<IN>) {
		chomp;
		next if $_ eq '';
		my @c = split /\t/;
		$e[$count++]=\@c;
	}
	
	close IN;	
	
	####### LOAD .psg ######
	@g1=();
	@g2=();
	if (open (IN, "$name.psg")   )  {
		my $header=<IN>;
		while (<IN>) {
			s/\r\n/\n/;
			chomp;
			#print $_, "\n";
			my @c = split /\t/;
			#print @c, "\n";
			push @g1, [ $c[0], $c[1], $c[2] ] if $c[0] ne '';
			push @g2, [ $c[3], $c[4], $c[5] ] if $c[3] ne '';
		}
		close IN;
	} else {
	}
	
	%acc=();
	@acc_order=();
	%opt=();
	if ( !open(IN, "$name.pso")   ) {
		print "WARNING: cannot load parasight option file ($name.pso)!" if !$opt{'quiet'};
		return;
	}
	while (<IN>) {
		s/\r\n//;
		chomp; chomp;
		my @c= split /\|\|\|\|/;
		if ($c[0] eq '#ACC') {
			$acc{$c[1]}=$c[2];
		} elsif ($c[0] eq '#ORDER') {
			push @acc_order, $c[1];
		} elsif ($c[0] eq '#OPT') {
			$opt{$c[1]}=$c[2];
		} else {
			print "WARNING: unidentified line in parasight option file ($_)!\n";
		}
		
	}
	close IN;
}


############################################
############################################
############################################
############################################

sub balloon_format_var {
	my $name=shift;
	my $desc=$optdesc{$name};
	my $text=&help_format ( $desc );
	$text="$text (variable:$name)";
	
	return $text;
}

sub fast_lentry {
	my ($frame,$tmp,$text,$var, $width, $bindsub,$color,$desc)=@_;
	##currently $tmp is an extraneous variable which could be replaced by the next nessary input##
	$color ||='white';
	my $varname='';
	if ($var !~/^SCALAR\(/ ) {
		$desc=$optdesc{$var};
		$varname=$var;
		$var=\$opt{$var};
	}
	my $label=$frame->Label(-text => $text)->pack(-side=> 'left',-anchor => 'e');
	my $tmp=$frame->Entry(-textvariable => $var ,-width=> $width, -background=>$color)->pack(-side=> 'left', -anchor=> 'e');
   $tmp->bind("<Return>", $bindsub) if $bindsub; 
	if ($opt{'help_on'}) {
		my $text=&help_format($desc);
		$ballooni->attach( $tmp, -justify=>'left',-msg=> "$text(variable:$varname)" );
	}

	return $label,$tmp;
}

sub help_format {
	my $desc=shift;
	my @c=split /\n/,$desc;
	my $max=0;
	foreach (@c) { $max = length ($_) if $max < length($_); }
	$max=$opt{'help_wrap'} if $max > $opt{'help_wrap'};
	my $text='';
	#simple wordwrap algorithm
	foreach my $line (@c) { 
		while (length ($line) ) {
			my $len=length($line);
			$len=$max if $len > $max;
			my $lens=$len;
			if ($lens == $max ) {
				$lens-- until substr($line,$lens-1,1) eq ' ' || $lens ==0;
				$lens=$len if $lens==0;
			}				
			my $piece = substr($line,0, $lens);
			substr($line,0,$lens)='';
			$text.=$piece;
			$text.= ' ' x ($max - $lens);
			$text.="\n";
		}
		
	}	
	return $text;
}	


sub generate_column_header_strings {
	$mstring='';
	for (my $i=0; $i< @mheader; $i++) {
		$mstring.= "$i)$mheader[$i]  ";
	}
	$estring='';
	for (my $i=0; $i< @eheader; $i++) {
		$estring.= "$i)$eheader[$i]  ";
	}	
}


sub color_change {
	my ($id,$ap,$headhashp) = @_;
	my ($type, $numb) = $id=~ /([A-Z])(\d+)/;
	my $oldcolor = $canvas->itemcget($id,-fill);
	my $newcolor= $canvas-> chooseColor(-title=> 'Choose New Color',
			-initialcolor=> $oldcolor);
	if (defined $newcolor) {
		$$ap[$numb][$$headhashp{'color'}]=$newcolor;
		$canvas->itemconfigure($id,-fill=>$newcolor);
	} 
}


sub fasta_format_wrap {
	my ($whole_seq, $width, $numb_on, $plus_num) =@_;
	my $m=0;
	my $broken_seq="";
	for ($m=0; $m+$width<length ($whole_seq); $m+=$width) {
		$broken_seq .= substr($whole_seq, $m, $width);
		$broken_seq .= "  ".($m+$width+$plus_num) if ($numb_on);
		$broken_seq .= "\n";
	}
	if ($m <=length($whole_seq)-1) {
		$broken_seq .=  substr($whole_seq, $m) ;
		$broken_seq .= "  ".length($whole_seq+$plus_num) if $numb_on;
		$broken_seq .= "\n";
	}
	return $broken_seq;
}


sub edit {
	my $id=$_[0];
	print "EDIT $id\n";
	my $ap=$_[1];
	my $ahp=$_[2];
	my @names=();
	my @data=();
	my ($type, $numb) = $id=~/([A-Z])(\d+)/;
	for (my $i=0; $i <@$ahp; $i++) {
		push @names, $$ahp[$i];
		push @data,$$ap[$numb][$i];
		#print "$$ahp[$i]  => $$ap[$numb][$i]\n";
	}
	my ($button, $pdata) = &data_dialog ("Edit Sequence Feature $type$numb", \@names, \@data);
	return if $button==0;
	@data=@{$pdata};
	for(my $i=0; $i<@$ahp ;$i++ ) {
		$$ap[$numb][$i]=$data[$i];

	}
}

sub doublelabel {
	my ($frame,$tmp,$text,$var)=@_;
	my $label=$frame->Label(-text=>$text)->pack(-side=>'left',-anchor=>'e');
	$frame->Label(-textvariable=>$var, -relief=>'sunken',-background=>'grey')->pack(-side=>'left',-anchor=>'e');
	return $label;

}

sub data_dialog {
	my $save;
	my ($title, $pn,$pd) = @_;
	#print "PD $pd\n";
	my $tl= $mw->Toplevel();
	$tl->title("$title");
	my $frame=$tl->Frame->pack(-side => 'bottom');
	$frame->Button(-text =>"Cancel",
				-command => sub {$save=0; $tl->destroy();}
		) -> pack(-side => 'left');
	$frame->Button(-text =>"Save",
				-command => sub { $save=1; $tl->destroy();}
		) -> pack(-side => 'bottom');
	
	my $text = $tl->Scrolled("Text", -width => 60, -height => scalar( @{$pn})*1.5+6,-wrap => 'none')
			->pack(-expand => 1, -fill => 'both');

	for(my $i=0; $i< @$pn; $i++) {
		my $w= $text->Label(-text=> $$pn[$i], -width => 20);
		$text->windowCreate('end', -window => $w);
		$w= $text->Entry(-width=>40, -textvariable => \$$pd[$i]);
		#print "$$pd[$i]\n";
		$text->windowCreate('end', -window => $w);
		$text->insert('end', "\n");
	}
	$text->configure(-state => 'disabled');
	$tl->deiconify();
	$tl->raise();
	$tl->grab();
	$tl->waitWindow();
	return ($save, $pd);


}

sub warnNpause {
	my $message=shift;
	chomp $message;
	warn "$message\n";;
	print "(PRESS ENTER TO CONTINUE!)"; 
	my $pause=<STDIN>;

}

sub focus_ok {
		my $dialog=$mw->Dialog(-title=>'Sorry',-text=>"This is not implemented!",
			-bitmap=>'question',-buttons => ['OK'],-default_button=>'NOK', );
}

sub  print_screen {
	my $printme=shift;
	my $newname=shift;
	$newname ||=0;
	my $filebase='screen';
	
	###############################################
	my @box=$canvas->bbox("all");
	print "$box[0] $box[1] $box[2] $box[3]\n";
	$canvas->configure(-scrollregion=>\@box);
	###use scroll region to figure out current visible position###
	my(@xview) = $canvas->xview;
	my(@yview) = $canvas->yview;
	my(@scrollregion) = @{$canvas->cget(-scrollregion)};
	my $x1 = $xview[0] * ($scrollregion[2]-$scrollregion[0]) + $scrollregion[0];
	my $y1 = $yview[0] * ($scrollregion[3]-$scrollregion[1]) + $scrollregion[1];
	my $x2 = $xview[1] * ($scrollregion[2]-$scrollregion[0]) + $scrollregion[0];
	my $y2 = $yview[1] * ($scrollregion[3]-$scrollregion[1]) + $scrollregion[1];
	
	print "($x1,$y1) ($x2,$y2)\n";
	my $canvas_height=$y2-$y1+1; ###vertival height
	my $canvas_width=$x2-$x1+1;  ###horizontal width
	
	print "H(Y) $canvas_height   W(X) $canvas_width\n";
	
	my ($dimension,$scalevalue)=print_determine_proper_scale_dimension($canvas_width,$canvas_height);
	
	if ($newname eq '1') {
		$filebase=&get_filename($filebase);
		return if $filebase eq '';
	} elsif ($newname !=0 ) {
		$filebase=$newname;
	}
	$filebase =~s/\.ps$//;
	$filebase.=".ps";
	print "$dimension\n";
	if ($dimension eq 'height') {
		$canvas->postscript(-file=>"$filebase",
			-colormode => 'color',
			-rotate => $opt{'printer_page_orientation'},
			-pageheight=> $scalevalue,
			-x=>$x1,
			-y=>$y1,
			-width=>$canvas_width,
			-height=>$canvas_height
		);
	} else {
		$canvas->postscript(-file=>"$filebase",
			-colormode => 'color',
			-rotate => $opt{'printer_page_orientation'},
			-pagewidth=>$scalevalue,
			-x=>$x1,
			-y=>$y1,
			-width=>$canvas_width,
			-height=>$canvas_height		
		);
	}
	print "Postscript file ($filebase) generated...\n";
	if ($printme) {
		&print_ps_file( "$filebase");
	}
}

sub print_determine_proper_scale_dimension {
	my $pixel_w=shift;
	my $pixel_h=shift;
	my $page_w=$opt{'printer_page_width'};
	my $page_h=$opt{'printer_page_length'};
	my $page_h_val=$page_h;
	$page_h_val =~ s/i//;
	my $page_w_val=$page_w;
	$page_w_val =~ s/i//;
	print "CW:$pixel_w, CH:$pixel_h  \n";
	if ($opt{'printer_page_orientation'}==1 ) {
		#landscape
		$page_w=$opt{'printer_page_length'};
		$page_h=$opt{'printer_page_width'};
	}
	print "PW:$page_w  PH:$page_h\n";
	if ( $pixel_w/$pixel_h < $page_w_val/$page_h_val ) {
		print "return 'height',$page_h\n";
		return 'height',$page_h;
	} else {
		print "return 'width', $page_w;\n";
		return 'width', $page_w;
	}
	
}


sub print_all {
	my $printme=shift;
	my $basename=shift;
	$basename ||= 'allout';
	my @box=$canvas->bbox("all");
	my $canvas_width=($box[2]-$box[0])+1;
	my $canvas_height=($box[3]-$box[1])+1;
	my $pages_x=$opt{'print_multipages_wide'};
	my $pages_y=$opt{'print_multipages_high'};
	my $page_w_pixels=$canvas_width/$pages_x;
	my $page_h_pixels=$canvas_height/$pages_y;
	my ($dimension,$value) = &print_determine_proper_scale_dimension($page_w_pixels,$page_h_pixels);
	print "$box[0] $box[1] $box[2] $box[3]\n";
	my $xpage=0;
	my $ypage=0;
	for (my $x= $box[0]; $x<$box[2]-1; $x+=$page_w_pixels) {
		$xpage++;
		$xpage=substr("000$xpage",-2);
		for (my $y=$box[1]; $y<$box[3]-1; $y+=$page_h_pixels) {
			$ypage++;
			$ypage=substr("000$ypage",-2);
			#print "$x ($xpage) $y ($ypage)    ", ($box[2]-$box[0])/$pages_x," ",($box[3]-$box[1])/$pages_y,"\n";
			if ($dimension eq 'height') {
				$canvas->postscript(-file=>"$basename.$xpage.$ypage.ps",
					-colormode => 'color',
					-rotate => $opt{'printer_page_orientation'},
					-x => $x, 
					-y => $y, 
					-width => $page_w_pixels,
					-height => $page_h_pixels,
					-pageheight=> $value,
				); 
			} else {
				$canvas->postscript(-file=>"$basename.$xpage.$ypage.ps",
					-colormode => 'color',
					-rotate => $opt{'printer_page_orientation'},
					-x => $x, 
					-y => $y, 
					-width => $page_w_pixels,
					-height => $page_h_pixels,
					-pagewidth=> $value,
				);  
			}
			print "Postscript output generated $basename.$xpage.$ypage.ps...\n";
			if ($printme) {
				&print_ps_file( "$basename.$xpage.$ypage.ps") ;	
				unlink "$basename.$xpage.$ypage.ps";
			}
		}
	}
}

sub get_filename (   ) {
	my $intialname=shift;
	my ($dir,$name);
	$name=$filepath;
	($dir,$name)= ($1,$2) if $filepath=~ /^(.*)\/(.*)$/;
	$name=$intialname;
	#print "POSITION FP:$filepath D($dir)  N($name)\n";
	my @filetypes= (['postscript', ['.ps']],['All Files', '*']);
	my $file = $mw->getSaveFile( -title=>'SAVE PARASIGHT FILES',-filetypes => \@filetypes, 
				 -initialdir=>$dir, -initialfile=>$name);
	return $file;
}

sub print_ps_file {
	my $file = shift;
	my $os=$^O;
	#print "$os\n";
	if ( $os =~ /MSWin/ ) {
		print "OS is $os.\n";
		print "Sorry, I haven't figure MSWin printing yet.\n";
		print "Grab the postscript ($file) and do it yourself!\n";
	} else {
		my $command=$opt{'print_command'};
		if ($command =~ s/[{][}]/$file/ ) {
			#skip ahead
		} else {
			$command.= " $file";
		}
		print "Printing $file\n";
		system $command;
		print "  ...done\n";
	}

}

sub commify {
	my $text=reverse $_[0];
	$text =~s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}
###################################################################################
#########################EXPORT AND EXTRACT TO TEXT ###############################
{

sub export_text {
	my $top = $mw->Toplevel();
	my $text=$_[0];
	my $name=$_[1];
	my $save_name=$_[2];
	print "SAVE_NAME:($save_name)\n";
	$top->title($name);
	my $textbox = $top->Scrolled("Text")->pack( -side=>'top' , -fill =>'both', -expand=> '1');
	my $f = $top->Frame->pack(-side => 'bottom', -fill => 'x');
	
	$f->Button(-text=>'Save As', -command =>  sub { 
				#$info= "Getting File Name For Save...";
				my $name=$export_text_path;
				my ($dir)='';
				$dir=$1 if $name=~ s/(^.*\/)//;
				$name=$save_name if $save_name;
				my @filetypes= (['All Files', '*'],['Text', '.txt']);
				my $file = $mw->getSaveFile( -filetypes => \@filetypes, -initialdir=>$dir,
											-initialfile=>$name, -defaultextension => '*',
											-title => 'Save As Unix Text');
				if ($file eq "") {print "Save Canceled!\n"; return; }
				open(FILE, ">$file") || die "Can't write text file ($file) [$!]\n";
				binmode(FILE);
				print FILE $textbox->get("0.0", 'end');
				close FILE;
				$export_text_path=$1		if $file=~ s/^(.*\/)//;
				$export_text_path.='untitled.txt';
			})->pack(-side=>'left');
	$f->Button(-text=>'Append To', -command =>  sub { 
				#$info= "Getting File Name To Append To...";
				my $name=$export_text_path;
				#print "ETP:$export_text_path\n";
				my ($dir)='';
				$dir=$1 if $name=~ s/(^.*\/)//;
		
				#print "NAME=$dir, $name\n";
				my @filetypes= (['Text', '.txt'],['All Files', '*']);
				my $file = $mw->getSaveFile( -filetypes => \@filetypes, -initialdir=>$dir,
										-initialfile=>$name, -defaultextension => '*',
										-title => 'Append Unix Text To');
				if ($file eq "") {print "Save Canceled!\n"; return; }
				open(FILE, ">>$file") || die "Can't write text file ($file) [$!]\n";
				binmode(FILE);
				print FILE  $textbox->get("0.0",'end');
				close FILE;
				$export_text_path=$1		if $file=~ s/^(.*\/)//;
				$export_text_path.='untitled.txt';
				#print "FINAL $export_text_path\n";
			})->pack(-side=>'left');
	##put text in
	$textbox->insert("end", $$text);

} #close sub

} #close private variable

sub alignment_format {
	my %args=( -seq1=>'seq1', -seq2=>'seq2', -begin1=>1, -end1=>'',-begin2=>1,-end2=>'',
				-name1=>'', -name2=>'', -orient2=>'F', -width=>60, @_);
	my $filehandle='';
	my ($seq1,$b1,$name1,$seq2,$b2,$name2,$orient2,$e1,$e2);
	($seq1,$b1,$name1,$e1)=($args{'-seq1'},$args{'-begin1'},$args{'-name1'}, $args{'-end1'});
	($seq2,$b2,$name2,$orient2,$e2)
			=($args{'-seq2'},$args{'-begin2'},$args{'-name2'},$args{'-orient2'},$args{'-end2'});
	my $width=$args{'-width'};
	my ($pos1,$pos2)=($b1,$b2);
	my ($pos1_del, $pos2_del)=(0,0);
	#####HEADER ##########
	$filehandle.="$name1 ($b1 to $e1) vs $name2 ($b2 to $e2)\n";
	$filehandle .="General alignment with ".length($seq1)." spaces\n";
	#$filehandle.= " ";
	####ALIGNMENT ########
	my $m=0;
	for ($m=0; $m+$width<length ($seq1); $m+=$width) {
		
		my $ss1=substr($seq1, $m, $width);
		my $ss2=substr($seq2, $m, $width);
		my $filler='';
		my ($num1,$num2)=('','');
		for (my $i=0; $i<length($ss1); $i++) {
			my ($b1,$b2)=(substr($ss1,$i,1),substr($ss2,$i,1) );
			$pos1++ if $b1 ne '-';
			$pos1_del=0 if $b1 ne '-';
			if ($b2 ne '-') {
				$pos2_del=0;
				if ($orient2 eq 'F') { $pos2++;} else {$pos2--;}
			}
			if ($pos1 % 10 ==0 && $pos1_del==0) {
				$pos1_del=1;
				my $add=' 'x($width+20) . $pos1;
				my $add_len= length ($num1) -$i-11;
				$num1.=substr($add,$add_len);
			}
			if ($pos2 % 10 ==0 && $pos2_del==0) {
				$pos2_del=1;
				my $add=' 'x($width+20) . $pos2;
				my $add_len= length ($num2) -$i-11;
				$num2.=substr($add,$add_len);
			}
			if ( $b1 eq  $b2 ) {
				if ($b1 ne 'N' && $b2 ne 'N') {$filler.='|'; } else { $filler.=' ';}
			} else {
				if ( $b1 eq '-' || $b2 eq '-' ) 	{$filler.=' ';} else {	$filler.='*';}
			}
		}
		$filehandle.= "\n";
		$filehandle.= " $num1\n";
		$filehandle.= substr($name1.' 'x10,0,9)." $ss1\n";
		$filehandle.= substr( $m+1 .' 'x10,0 ,10). "$filler\n";
		$filehandle.= substr($name2.' 'x10,0,9)." $ss2\n";
		$filehandle.= " $num2\n";

	}
	if ($m <=length($seq1)) {
		my $ss1=substr($seq1, $m);
		my $ss2=substr($seq2, $m);
		my $filler='';
		my ($num1,$num2)=('','');
		for (my $i=0; $i<length($ss1); $i++) {
			my ($b1,$b2)=(substr($ss1,$i,1),substr($ss2,$i,1) );
			$pos1++ if $b1 ne '-';
			$pos1_del=0 if $b1 ne '-';
			if ($b2 ne '-') {
				$pos2_del=0;
				if ($orient2 eq 'F') { $pos2++;} else {$pos2--;}
			}
			if ($pos1 % 10 == 0 && $pos1_del == 0) {
				$pos1_del=1;
				my $add=' 'x ($width+20) . $pos1;
				my $add_len= length ($num1) -$i-11;
				$num1.=substr($add,$add_len);
			}
			if ($pos2 % 10 == 0 && $pos2_del == 0) {
				$pos2_del=1;
				my $add=' 'x ($width+20) . $pos2;
				my $add_len= length ($num2) -$i-11;
				$num2.=substr($add,$add_len);
			}
			if ( $b1 eq  $b2 ) {
				if ($b1 ne 'N' && $b2 ne 'N') {$filler.='|'; } else { $filler.=' ';}
			} else {
				if ( $b1 eq '-' || $b2 eq '-' ) 	{$filler.=' ';} else {	$filler.='*';}
			}
		}		
		$filehandle.= "\n";
		$filehandle.=" $num1\n";
		$filehandle.= substr($name1.' 'x10,0,9). " $ss1\n";
		$filehandle.= substr( $m+1 .' 'x10,0 ,10). "$filler\n";
		$filehandle.= substr($name2.' 'x10,0,9)." $ss2\n";
		$filehandle.= " $num2\n";
	}
	return \$filehandle;
}


#################################################
#################################################
####### FASTASUBSEQ SUBROUTINES #################
#################################################

sub fasta_getsubseq_whole {
	my $file=shift;
	my $begin=shift;
	my $end=shift;
	my $position=0;
	my $subseq='';
	open (FASTA,"$file") || die "Can not open fasta file($file)!\n";
	my $header = <FASTA>;
	$header =~ s/\r\n/\n/;
	chomp $header;
	while ( <FASTA> ) {
		s/\r\n/\n/;
		chomp; 
		s/\s+//;
		$position += length;
		if ($position >= $begin) {
			my $start= $begin - ($position- length ) -1;
			$start=0 if $start<0;
			#print "START:$start\n";
			$subseq.= substr($_,$start);

		}
		last if $position >= $end;
	}
	close FASTA;
	if ($end > $position) {
		warn "Trim ($begin to $end []) too long! Length is $position\n" if $end > $position;
	 	return '';
	 }
	$subseq = substr($subseq, 0, $end -$begin+1 );
	return $subseq;
}

#####create header ######
sub find_file {
	my ($paths, $names) = @_;
	my @paths=split ":",$paths;
	my @names=split ":",$names;
	for my $path (@paths) {
	 	for my $name (@names) {
	 		my $p = $path;
	 		$p .= "/$name";
	 		#print "TESTING $p\n";
	 		return $p if (-e $p);
	 	}
	}
	return "";
}

###this program uses a global variable called $subseq
sub fasta_getsubseq_frac {
	print "GETTING A FRACTIONATED FILE\n";
	my $opt_f=shift;
	my $opt_b=shift;
	my $opt_e=shift;
	my $fasta_frag_size=shift;
	my $bfasta=int(($opt_b-1)/$fasta_frag_size);
	my $bpos=$opt_b-$bfasta*$fasta_frag_size;
	my $efasta=int(($opt_e-1)/$fasta_frag_size);
	my $epos=$opt_e-$efasta*$fasta_frag_size;
	print "$bfasta:$bpos:$efasta:$epos\n";
	my $subseq='';
	my $header='';
	for (my $i=$bfasta; $i<=$efasta; $i++) {
		my $path= $opt_f."_" . substr("000$i",-3);
		my $begin=1;
		my $end=$fasta_frag_size;
		$begin=$bpos if $i==$bfasta;
		$end=$epos if $i==$efasta;
		if ($begin != 1 || $end != $fasta_frag_size) {
			print "EXTRACTING SUBSEQ:$path:$begin:$end(I:$i)\n";
			$subseq.= &getsubseq_frac($path,$begin,$end,$fasta_frag_size);
		} else {
			print "EXTRACTING WHOLESEQ:$path:$begin:$end(I:$i)\n";
			$subseq.= &getwholeseq_frac ($path, $fasta_frag_size);
		}
		print "SUBSEQ LEN:",length($subseq);
	}
	return $subseq;
}


sub getsubseq_frac {
	my $fasta=shift;
	my $begin=shift;
	my $end=shift;
	my $fasta_frag_size=shift;
	my $subseq='';
	my $position=0;
	#this sub returns nothing because it is built for speed#
	open (FASTA,$fasta) || die "Can not open fasta file $fasta\n";
	my $header=<FASTA>;
	die "Improper fasta $fasta\n" if $header !~/^>/;
	while ( <FASTA> ) {
		s/\r\n/\n/;
		chomp; chomp;
		s/\s+//;
		$position += length;
		if ($position >= $begin) {
			my $start= $begin - ($position- length ) -1;
			$start=0 if $start<0;
			#print "START:$start\n";
			$subseq.= substr($_,$start);

		}
		last if $position >= $end;
	}
	die "Extracting $fasta ($begin to $end)failed! Request to long!(end length:$position)!\n" if $end > $position;
	$subseq = substr($subseq, 0, $end -$begin+1 );
	close FASTA;
	return $subseq
}

sub getwholeseq_frac {
	my $position=0;
	my $fasta=shift;
	my $fasta_frag_size=shift;
	my $subseq='';
	open (FASTA,$fasta) || die "Can not open fasta file $fasta\n";
	my $header=<FASTA>;
	die "Improper fasta $fasta\n" if $header !~ /^>/;
	while ( <FASTA> ) {
		s/\r\n/\n/;
		chomp;
		s/\s+//;
		$position += length;
		$subseq.=$_;
	}
	die "End $position not equal to $fasta_frag_size\n!" if $fasta_frag_size !=$position
}

sub breaksequence {
	my ($whole_seq, $width, $numb_on, $plus_num) =@_;
	my $m=0;
	my $broken_seq="";
	for ($m=0; $m+$width<length ($whole_seq); $m+=$width) {
		$broken_seq .= substr($whole_seq, $m, $width);
		$broken_seq .= "  ".($m+$width+$plus_num) if ($numb_on);
		$broken_seq .= "\n";
	}
	if ($m <=length($whole_seq)-1) {
		$broken_seq .=  substr($whole_seq, $m) ;
		$broken_seq .= "  ".length($whole_seq+$plus_num) if $numb_on;
		$broken_seq .= "\n";
	}
	return $broken_seq;
}




############################################################
############### DOCUMENTATION ##############################
############################################################
############################################################


#pod commands  I<italic> B<bold>  C<code>  S<no break space> E<lt>  E<34> E<gt>
#  L<name>  man page
# F<filename>
# X<index entry>
# Z<zero width text>

=pod

=head1 NAME

B<parasight> (version 7.4)

=head1 SYNOPSIS

 parasight -align alignment.table

=over 5

This simply loads the file F<alignment.table> containing either a table of tab-delimited alignments or miropeats standard output (see below).

=back

 parasight -align AC002038.blast.parse -showseqqueryonly

=over 5

This simply loads the file F<AC002038.blast.parse> containing parse blastdata and displays the hits relative to the query sequence.  This is the number one use in our lab for parasight.

=back

 parasight -align AC002038.blast.parse -showseq AC002038.1: 
     -extra repeats -template bacblastview.pst 
     -options 'seq_color=>red, canvas_width=>1000'

=over 5

This draws the blast output from a search with AC002038.1 formatted with the options contained in the template file F<bacblastview.pst>.  It uses B<-options> to modify the screen width and sequence color.

=back

 parasight -in saved.parasight -showseq AC002304:AC002035: 
     -arrangeseq sameline -template template.file 
     -options 'seq_color=>red,extra_arrow_on=>1'

=over 5

This loads a previously saved parasight view (the files F<saved.parsight.psa>, F<saved.parasight.pse> and F<saved.parasight.pso>, shows only 2 of the sequences (AC002034 and AC002035), arranges or places these two sequences on the same line, loads a template file of options to reformat the view, and modifies two options directly (sequence color and turns on arrows for annotation extra).

=back


=head1 DESCRIPTION

Parasight is a generalized pairwise alignment viewer originally developed for analyzing segmental duplications (or paralogy) within the human genome. It is designed to display the positions and relationships of pairwise alignments within sequeunce(s).  It provides both interactive analysis as well as publication quality postscript output.  Parasight can arrange and color alignments on the basis of any other included data such as size, percent similarity or even species designation. It can also display the position of any type of simple sequence annotation such as repeats and exons. Finally, it can graph numerical data in relation to the seqeunce such as windows of percentage GC content.  Parasight has been used to analyze output from programs such as B<BLAST>, B<miropeats>, and B<pip maker> from the scale of whole-genomes to (such as segmental duplications in the human genome) to the analysis of a single protein searched against a database of interest.  If it is pairwise data, parasight can display the output. 

Parasight functions on both Unix and Windows platforms. The program is written in Perl using the graphical Perl/Tk module. It was designed to be extremely flexible and thus a price is paid in terms of speed as well as the complexity, i.e. the large number of options). However, the numerous options makes it more likely that parasight can do what you want it to do.  Although not necessary for basic interactive use, an understanding of regular expressions and a familiarity Perl is helpful in order to utilize the program more fully.  Parasight and its options are fully accessible through the GUI interface, the command line, or via loadable templates making anlyses flexible and automatable. Most users of parasight load their data into the program and then format the view interactively using the extensive options menu (and now templates).  Programmers will be the most likely to use command line manipulation of internal options. Parasight has been used and tested extensively on both Linux and MS Windows.  The Unix version is the most extensively tested and all options should be available.  Windows lacks some of the more advanced options due to incompatiblities/inflexabilities of Bill Gates' operating system.


To extol parasight's strengths:

=over 5

=item 1 Flexible

The RAM is the limit when it comes to loading data.  Other than the most basic description of a pairwise alignment or extra sequence feature Parasight makes absolutely no assumptions about your data allowing the user to analyze what they are interested in analyzing. Technically parasight can't tell a bp apart from an inch or DNA from Protein.

=item 2 Formatible

Parasight has a pletora of options, if I have ever needed it parasight has it.  Every parasight option is available is availble from the command line, the GUI interface or from a saved template file of options.  Thus, the basic user and the programmer can completely tailor their parasight views to their exact needs.

=item 3 Interactive

Parasight can interact with the user and with other programs. The user can interaction format the parasight image via the  GUI option menu as well as edit the data. Parasight allows the user to print screen shots or dump a postscript of the entire image.  Popup windows over alignments and extra sequence features display the objects data.  In addition to scaling with the options menu, users can zoom in and out to gain an appreciation of the detail. Parasight has the ability to link to (or execute) other programs allowing the viewing of web pages or associated sequence alignments at the bp level.

=item 4 Programmable

Parasight has the ability to accept additional Perl code from the command line or via a file, which allows for more complex formatting or for the execution of commands such as searching or printing.  Combined with the -die option this allows for powerful batch processes such as generating PostScript images of 30,000 BACs (if you are so inclined).

=back


=head1 UNDERSTANDING DATA CATEGORIZATION/CLASSIFICATION

A basic understanding of the logic behind parasight is useful in understanding data input and manipulation.  Data falls into three categories, pairwise alignments, extra sequence annotation and graph data. The graphical option menu is organized on the basis of which data is being manipulated. Alignments, the core of parasight, have two forms of display: pairs and subs.  Pairs are representations of the pairwise alignments that are normally drawn atop the seqeunce.  For each alignment the pairs representing it can be connected by lines to show their relationship.  Thus, for pairs relationships are only visable if the both sequences containing the pairwise are drawn.  Subs are representations of pairwise in relation to only one of the sequences mimicing blast type output.  Sub is for sub-sequence (as they are drawn below the seqeunce) or subjects (if you are examining BLAST results). Extras are simple sequence annotations that have one beginning and one end such as introns, LINEs, SINES, motifs, etc.  In the case of a gene the intron exon structure can not be drawn as one object, but only as individual exons and indvidual introns. (A -gene data structure is planned for the distant future.)  The last data type is graph data Graph data is a plot of sequence positions (x-axis) versus a numerical value (y-axis).

Below is a crude schematic (the best I could do in POD) of a typical display with Sequence(-), Pairs(P), Subs(S), Extras(E), and Graph(G) data.

                                                                 G
                   G                   G         G                   G 
          G               G     G                               G           
               G                                          G                 
               
           EEE   EEEE          EEEE         EEEE      EEE     EEEEEEE
    S0001--------PPPPPPPPPPPPP-----------PPPPPPPPPPPPPP--PPPPPPPP-----
    SEQ04        SSSSSSSSSSSSS                           SSSSSSSS
    SEQ02        SSSSSS   SSSS           SSSSSSSSS
    SEQ03           SSSSSSSSSSS              SSSSSSSSSSS  SSSSS    
                             
                            
 
 

=head1 COMMAND LINE OPTIONS

The main command line arguments available when parsight is executed can be divided into three main headings:  L<data input|/"DATA INPUT">, L<reloading a saved parasight view|/"RELOADING A SAVED PARASIGHT VIEW">, and L<changing view options|/"CHANGING DISPLAY OPTIONS">.


=head2 DATA INPUT

While many types of data can be loaded and displayed, the absolute minimum input is simply the length of a sequence to be drawn.  The length of a sequence can be provided in B<-showseq> file (see below).  Usually the lengths are supplied as part of the pairwise alignment file.  If no alignments are being drawn then the user must supply the lengths with the B<-showseq> option.  Of course examining a line representing a sequence is pretty boring--even if it is decorated with tick marks--so parasight has a few other data input options:

=head3 Option: -align

B<-align [I<filepath1:filepath2:filepath3:etc>]> loads  files containing pairwise alignments.  The files must be either the saved standard output for B<miropeats> (Jeremy Parsons) or a tab-delimited format akin to miropeats standard output. The tab-delimited format is simply a file where the first 8 columns contain the pairwise coordinates and lengths of the two similar sequences. The align file is assumed to have a descriptive header in the first row.  Hence, the first alignment will be lost (and loaded as the header) if no header actually is present.

=head4 Miropeats standard output

An example of Jeremy Parson's Miropeats standard output is:

 ## Minimum repeat length set to 300.
 
 
         ICAass  Version 2.1
         =======
 
 
 Indexing all the sequences now. This may take a few minutes.
 
 Total of 1 sequences indexed
 The sorted index is being saved to the file cluster.index.7507 
 .AC002038 118 731 161973 AC002038 44681 44068 161973
 AC002038 1299 1788 161973 AC002038 47175 47664 161973
 AC002038 22870 23591 161973 AC002038 39920 40641 161973
 AC002038 46067 46524 161973 AC002038 26363 26820 161973
 AC002038 46067 47435 161973 AC002038 26363 27731 161973
 AC002038 46699 47435 161973 AC002038 26995 27731 161973
 AC002038 47175 47664 161973 AC002038 1299 1788 161973
 Graphic ready for printing - type the command shown below to print:
 lp threshold300

=head4 Example tab-delimited alignment file

An example of a tab-delimited align file consisting of parsed BLAST output with additional data columns is shown below:

 name1 begin1 end1  len1  name2 begin2 end2  len2  similarity  transversions
 S001  1322   20001 20001 S002  1      18064 18064 0.945632    125
 S001  1322   20001 20001 S003  1      21010 21010 0.980581    143
 S002  1      18064 18064 S003  100    21010 21010 0.999587    7
 S002  1      18064 18064 S004  1      19041 19041 0.989587    43
 S002  1      12191 18064 S005  1      12141 18073 0.997548    17
 S002  12799  18064 18064 S005  12809  18073 18073 0.998548    3

The format consists of 2 sequence names, the coordinates of the pairwise similarity, and the overall lengths of the named sequences (C<name1 begin1 end1 len1 name2 begin2 end2 len2>), which must always be in this given order.  The first row of the alignment file contains column header names. For the first 8 rows the header rows are ignored and thus it is necessary to place these eiqht columns in the exact order given above.  The only data from these first 8 rows that may be omitted is the overall lengths of the sequence.  However, these columns must still be present and contain no value (empty cell in Excel).  Also, the lengths of the sequence must then be provided via B<-showseq>.  Any additional columns (such as the C<similarity> and C<transversions> columns above) are kept within the internal alignment data table. This additional data can be used to format and filter the parasight views generated. Additional columns that are created if not present in the alignment file are C<color>, C<width>, C<offset>, C<sline>, C<scolor>, and C<hide>. These are case sensitive and all lower case. ( C<color> contains the color of the pairwise.  C<width> is the width or thickness for the bar representing the pairwise.  C<offset> is the offset of the subject object.  C<scolor> is the color of a subject object. C<hide> does not display a pairwise if it is equal to 1. If the values for these columns are not inputed or blank then teh default values for the options are used.  (NOTE: It is usually simpler to modify these formatting columns in saved parasight tables using programs such as Excel.)  

=head3 Option: -showseq

B<-showseq [I<filepath | seqname1[,length,begin,end]:seqname2[,length,begin,end]:etc>]> displays only the designated sequences.   With no colon a filename is assumed and the program attempts to load it. If a colon is found in the option then it is assumed that the input is a colon separated group of sequence names.  Optional length and begin and end positions can be designated as well. This information may be given on the command line using commas after the sequence name or be contained within a file (tab-delimited).  For analysis such as BLAST searches where you just want to display the query sequence it is easier to use the short cut option B<-showseqqueryonly> in combination with B<-showseq> C<ALL>.  

=head4 Example showseq file

The format of the tab-delimited show file is shown below:

 seqname    length   begin   end
 S001       10000    50      1000
 S002       15432    1000    15432
 
An example of the data above as a command line entry is:

 -showseq S001,10000,50,1000:S002,15432,1,15432:

Or with just the lengths to display the entire sequence:

 -showseq S001,10000:S002,15432:
 
Or if the lengths are in the alignment file, just the begin and end positions can be designated (note the double comma to skip the sequence length):

 -showseq S001,,50,1000:S002,,1000,15432

The only required column is the sequence name. The names must be exactly the same as in the alignment file and extra file. The lengths, begins, and ends, are optional for the sequences to be drawn.  However, if length, begin or end are used, they must appear in the proper columns.  Begin and end must always be found in the 3rd and 4th columns and thus a blank sequence length must be provided for column 2 when sequence length is not designated but begin and end are. If lengths are not supplied in the alignment files or the show file then errors will occur. Lengths designated by -showseq always supercede lengths found in alignments.

B<-extra [I<filepath1:filepath2:etc>]> loads any 'extra'  sequence annotation/feature that can be expressed as a continuous block of the displayed sequence (i.e. a simple begin and end position). This can include features such as high copy repeats, introns, exons, and genes (if you don’t care about introns/exon structure).  The simplest extra file contains 3 columns in the given order C<seqname begin end>.

=head4 Example of an tab-delimited extra file

An example of a tab-delimited extra file is shown below:

 seqname   begin   end    name   color  offset
 S001      50      1000   exon1  blue   -10
 S001      5000    5500   exon2  blue   -10
 S002      5000    9000   LINE1  red    -20


Columns added if not present are C<color>, C<offset>, C<width>, and C<orient>.
(Again it is usually simpler to modify formatting data after the fact unless it is generated beforehand (e.g. C<orient>).  In terms of formatting simple color names should work: black, red, green.  Orientation should be either 'F' or 'R' (capitalized). Plus and Minus will I<NOT> work.  However, during the initial loading of the extras they as well as other common designations for orientation are automatically changed to 'F' and 'R'.  Other columns may be added that give additional information such as names and descriptions. 


=head3 Option: -graph1 or -graph2

B<-graph1 I<file(:s)> and -graph2 I<file(:s)>> loads simple graphing data in the form of C<seqname position value>.  The C<position> within the sequence in bases and the C<value> must be numerical (floating point). It plots the points and/or a connecting line.  The graph is numerical values of sequence positon on the x-axis versus y-axis numerical values.  -graph1 is used to generate one plot.  -graph2 is used to plot another line or set of data points on the same scales.  The left and right axis can be scaled for different ranges.  Thus, Alu content and GC content can be graphed at the same time. The left axis shows the scale for -graph1 and the right axis the scale for -graph2. No header is required for the input file. An unknown value can be designated with an empty value position. An empty value position causes a discontinuous line to be drawn. Only the first 3 columns are loaded all additional columns are ignored.  Graphing is built for speed-not flexiblity. The only flexibility is in scaling and formating the axes.

=head4 Example of a graph file

 seqname  position  value
 chr1      5000      0.43
 chr1     10000      0.65
 chr1     15000      0.73
 chr1     20000      0.65

=head2 RELOADING A SAVED PARASIGHT VIEW

Data can be saved directly as parasight formatted files.

=head3 Option: -in

B<-in [I<base filepath>]> loads a previously saved parasight dataset. Data is saved in 4 separate files (F<basefilename.psa>, F<basefilename.pse>, F<basefilename.pso>) and F<basefilename.psg>. Each file is editable text. The F<.psa>, F<.pse>, and F<.pso> are required even if there are no alignments and/or extras. These extensions to the basefilename are automatically searched for in the given path.

The F<.psa> and F<.pse> are tab-delimited tables containing the alignment and extra data, respectively. These tables are easily edited with any text editor; spreadsheets such as Excel are particularly useful in modifying these tables since the data is separated into columns and calculations can easily be done to modify the data as necessary.

The F<.pso> file contains all of the current option information. It is saved as text which can be modified by the end user.  It has a similar format to the template files.  

The F<.psg> file contains all of the current data to graph.  The -graph1 data is stored in the first 3 columns and the -graph2 data is stored in the next 3 columns (4 to 6).  The graph file, unlike the alignment and extra files is only created if graph data has been loaded.  Thus, a missing *.psg will will not generate an error.  For each set of 3 columns, column one is the seqeunce, column 2 is the position on the sequence, and column 3 is the value of to plot on the y-axis.

IMPORTANT:  All data necessary for a parasight view is contained in these files I<EXCEPT> for any B<-showseq> files or B<-arrangeseq> files.  These files must still be accessible in the same relative path positions in order for the saved file to be loaded properly.  In other words, only the file names to a show and arrange files are saved and that data must be reloaded. If the files get moved then the link will be broken and their paths will need to be altered.

=head2 CHANGING DISPLAY OPTIONS

Option arguments modify and format the parasight view. All of these options may either be change from the command line, a template file, or interactively within the program OPTION menu. The interactive menu is the easiest way to learn and template files the easiest way to apply a set of options again and again.  Changing options at the command line follows a set order of precedence--whereby old options loaded from a previous parasight view (B<-in>) are overridden by an option template file (B<-template>), both of which are overridden by any options specified in (B<-option>) command. All of these are overridden by direct command line options such as B<-arrangeseq>, B<-colorsub>, and B<-showsub>.

PRECEDENCE: B<internal default ---E<gt> -in ---E<gt> -template ---E<gt> -option ---E<gt> commandline> 

=head3 Option: -template

B<-template [I<filepath>]> loads an option template file.  This allows a user to quickly format future parasight views so that they are just like the saved one.  It is created using the save option template in the file menu.  When loading a template, if the file is not found in the current or specified path, hard-wired default template directories are searched.  For our lab one directory contains templates shared among multiple users.  And a user specific directory for an individuals PARASIGHT files.  The $template_path variable contains the paths.  To modify them you must modify the code. The current setting is '~/.PARASIGHT:/people/PARASIGHT'.  The search is left to right and first one found is first one used.    Template directory as it is currently set does not work for WINDOWS.  The tilde must be removed as it only works on Unix where the HOME directory (~) is designated by the environmental variables..

The template is an standard text file so a user can modify the values easily. It is created using the save option template in the file menu.  0 and 1 are used for on and off as well as yes and no values.  An empty string is simply a line return right after the (=E<gt>)  A line beginning with ### is ignored and is used to give descriptions of the values.  Be careful about adding blank space, it is a good idea to edit with normally unseen characters such as spaces and line breaks visualized.

=head3 Option: -options

B<-options [I<'opt1=E<gt>value1,opt2=E<gt>value2'>]> is a list of options to modify. All of the underlying options are available; however, there are probably many that you will never have reason to modify, but they are all listed in appendix A for completeness.

=head3 Option: -showsub

B<-showsub [I<filepath | seqname1:seqname2:etc>]> This option shows only the designated subs to be drawn under sequences.  Multiple sequence names can be directly entered with colon delimitation.  If no colon is present then the input will be treated as a file containing a list of subs and will be loaded.  Default is C<ALL>, which displays all possible subs.

=head3 Option: -arrangeseq

B<-arrangeseq [I<oneperline | sameline |file:filename>]> This option arranges the sequences in a specified manner.  

=over 5

I<oneperline> draws each sequence on a separate line that may wrap if needed.  

I<sameline> draws all of the sequences on the same line with a given amount of spacing between them.  

I<file:filename> uses the data in the file to arrange the sequences in user defined pattern.  The file consists of two columns C<seqname> and C< position > in current line.  To start a new line C<NEWLINE> is typed alone.  The example below places the chromosomes on 3 lines. 

=head4 Example Arrange File


	acc	start
   chr1	400000000	
   chr6	1668388704	
   chr7	1870803946	
   chr8	2057427852	
   chr9	2230204273	
   NEWLINE		
   chr22	1	
   NEWLINE		
   chr10	400000000	
   chr11	565589288	
   chr12	736372841	
   chr13	900655330	
   chr14	1040400228	
   chr15	1167353549	

=back

B<-arrangesub [I<oneperline|stagger|subscale|cscale>]> This option arranges subs below the drawn sequences.  The name came from blast subjects, but you can also think of them in terms of sub (beneath) the sequence.

=over 5

I<oneperline> = each sub sequence is placed on its own line underneath the drawn sequence.  The ordering of sequences can be altered by choosing a column to sort on (arrangesub_col) 

I<stagger> = multiple subjects are placed on same line only when non-overlapping. The spacing required between the beginning and end of two subs can be varied.  This spacing gives room for labels.  The ordering starts in terms of other sequences with hits closest to the beginning of the sequence of interest under which the subs are being drawn

I<subscaleN> = subjects are places on a numerical scale based on given column values. Tricky so avoid setting up from command line--use the GUI and then save a template from that.

I<subscaleC> = subjects are placed on categorical scale based on column values. Tricky so avoid setting up from command line.  Use a template or the GUI.

Note: the best way to figure the scales out is to experiment with them interactively in the options menu .There are specific modifications of subscaleN and subscaleC that are included as choices.  They are denoted by a preceding asterisk and were developed to display breakdowns of percent similarity and chromosome position (for mostly oudated draft versions of the genome) However,they may be instructive to the new user. New views are now simply done via a template rather than adding adding even more choices.
=back

B<-color >  I<***not implemented***> When implemented it will color the pairwise sequences and connecting lines. Currently, coloring is only based inter and intrachromosomal designation. (As of yet the need hasn't really arisen.)  For consistency this should be called colorseq.

B<-colorsub [I<NONE|RESET|seqrandom|hitrandom|hitconditional>]> This option provides color schemes for the subs drawn below the sequence.  

=over 5


I<NONE> does not change the color and leaves hit colors intact. Hit colors are stored within each pairwise in the table.  Subject colors are stored transiently.  Hit colors over-ride subject colors.  To remove hit colors use RESET. 


I<RESET> removes hit (individual pairwise) colors, which override any assigned subject colors.  For example, if you use hitrandom and then try to switch to seqrandom, nothing will change.  This is because hitrandom colors are still stored in the internal alignment table and they take precedence over the subject color scheme.  Thus, this intermediate C<RESET> step is required to clear the hit colors. CAUTION: if you use RESET all of your manual coloring will be wiped out. (NOTE: This is because hit colors reside in the same column C<scolor>as manually modified sub colors.  The column C<color> defines the pairwise color--overriding inter and intra colors.)  Sorry, this is part of the program that could be simplified if I ever have a chance to gut it.

I<seqrandom> randomly assigns colors to the various sequences that are displayed as subs. (There is a random set of 20 odd colors that are cycled through.) 

I<hitrandom> randomly assigns colors to each individual hit or pairwise alignment. (There is a random set of 20 odd colors that are cycled through.)

I<hitconditional> allows for each pairwise to be assigned a color based on pseudo-Perl code by using a series of conditional statements that test a single alignment column.  Basic syntax is C<[color] [test] [value];>, where color= color to set, test is =, E<gt>, or E<lt>, and value is some numerical value.

=back

B<-minload> is a switch to load only the alignments and extras for the sequences that will be drawn as designated by -showseq. It is very useful for increasing the speed of the program when there are a large number of alignments that will not be drawn in the current view.  Why load the genome if you only want to look at chromosome 22?

B<-precode I<['Perl code']>> This code is executed after the initial drawing of objects. It allows automation for batch processes when combined with die option. (See Advanced option section below for details.)

B<-die> parasight quits after executing the precode option (See Advanced option section below for details.)

=head1 INTERACTIVE MENUS

=head2 RESHOW, REARRANGE, REDRAW

This part is to answer why there is a blue and white button for updating the drawing.  For beginners, I simply suggest using the blue R,R&R (Reshow, Rearrange, and Redraw) button.  For extremely large data sets; however, the Reshow, and Rearrange calculations can take a significant amount of time.  Thus, if you are just changing the spacing of tick marks it is handy to skip the sequence and arrangement calculations.  However, for simple views of BAC BLAST output stick with the blue button.  

=head2 OPTION MENU

The option menu has popup help (over yellow text) and most options are self-explanatory. If in doubt try changing an option and see what happens.  I have tried to adhere to a semi-logical naming convention when ever possible.  Blue color coding is to show whether a variable will require reshow and rearrangement before taking effect. The menu is subdivided into 6 main parts: MAIN, SEQ/PAIRS, SUBS, EXTRA, GRAPH, FILTER, and MISC.  The organization trys to follow the organization of the data in parasight.

The MAIN menu allows access to important command line options like C<-showseq> and C<-showsub>.  Also, basic screen properties such as size of the window and the number of bases for the width of the screen.

The SEQ/PAIRS portion allows manipulation of the sequence and assocaiated tick marks.  Pairs and their designation as inter and intrachromomal as well as connecting lines are controled from this part of the menu as well.

The SUB portion of course is all about the manipulation of subs.  This is some of the more complex data manipulation.

The EXTRA portion is about the options relating to the extra data.

The GRAPH portion is for the graph data.  Try turning everything on when you first test out this feature.

The FILTER portion allows for the filtering/removal of pairwise and extras based on data in a given column of numerical data. 

The MISC portion allows for the setting of options controling printing, the display of alignments, the extraction of sequence, and the execution of other programs. 


=head2 FILE DROP-DOWN MENU

This is the only place where the save parsight command is found.  All data and options are saved. A few files are not saved--see information about B<-in>.  Loading must be done at the command line. Additionally, template files (F<*.pst>) may be saved and loaded through this menu.  After loading a template file the screen must be R,R,& R.

=head2 PRINT DROP-DOWN MENU

The print menu allows for the generation of a postscript file and its subsequent transmission to a printer if the option print_command is properly set.  The postscript file can consist of the visible screen (screen) or the entire parasight drawing (all).  If the all option is chosen then the number of pages (vertically and horizontally) across which to print the image is set with the option print_multipages_wide and print_multipages_high.  The postcript files are encapsulated and can be easily turned into PDF files with software such as Adobe Distiller or imported into Adobe Illustrator.  Also, word has a special eps import option which was handy when writing my dissertation.

=head2 ORDER DROP-DOWN MENU

The order menu on the main drop down menu bar allows the order or level of objects in the display to be changed.  You can either send objects all the way to the background or the foreground.


=head2 MISC DROP-DOWN MENU

Currently it contains the ability to transfer colors for alignments in order to allow syncing of colors between the pairs and the subs.  It requires a redraw to see the effect after choosing one of these options.  This is really the only way to currently go outside of the inter vs intra coloring schemes for pairs.

=head1 SCREEN MANIPULATION

In addition to gazing lovingly at the pretty images after formating them using the option menu, direct manipulation of the display once drawn can be accomplished with various commands.

=head2 MOUSE BUTTON FUNCTIONS 

(see APPENDIX B: for table of mouse functions)

First when you mouse over an object it will shimmer with a number of bright colors. the shimmering object represents the object you will select if you click on it.Most of the mouse commands work on sequence, pairwise, extra, and subjects.  Tick Marks and Labels are immune except for the ALT buttons. The middle mouse button is not used since some systems like my home PC lack them (and I don’t have the dexterity to precisely click both Left and Right at the exact same time which is the usual substitute).

B<DATA POPUP WINDOW (Left-Click)> This pops up a simple window displaying all data for an alignment or an extra object.  Use Shift-Drag to move the popup window if it is obscured or obscuring data. Formatting options for this popup window are found under MISC tab of the OPTIONS menu.

B<OPTIONS POPUP (Right-Click)> Brings up a popup menu of options, which includes a variety of commands such as choosing colors and editing the underlying data.  If the actual alignments are present in the alignment table, the alignments can be viewed.  If the underlying sequence files are available, subsequences representing objects can be extracted.

B<ZOOM IN AND OUT (Control-Left-Click and Control-Right-Click)>
Zooming can be accomplished with Control held down at the same time as a mouse click. The left mouse clicked in conjunction with the control key will zoom in two fold centered at the point of the click.  The right mouse has the opposite effect and zooms out.  The B<DeZoom> button on the main window returns the scaling to normal.

B<MOVE OBJECT TO FOREGROUND OR BACKGROUND (Alt-Left-Click and Alt-Right-Click)> This causes the object clicked on to move all the way to the foreground or the background. The left mouse button moves it to the foreground. The right mouse button moves the object to the background.

B<MOVE (nonpermanent) ANY OBJECT (Shift-Left Drag)> 
Allows for the movement of object in the drawing--even tick marks and sequence lines.  It is non-permanent but it is useful for removing tick marks or names before you print or create a PostScript file.

B<QUICK COLOR (Shift-Right-Click to color and Shift-Right-Double-Click to uncolor)>
Allow for rapid coloring of objects.  Shift-Right Click causes the object's color to change to that of the Quick Color Button on the Main Window.  Shift-Double-Click-Button attempts to remove the color and leave the default color.  In the case of Pairs, black is assigned to the object as inter and intra chromosomal colors can not be reassigned until a Redraw.  Coloring of all other objects (i.e. not extras and not alignments) are not saved or stored and consequently revert to normal as soon as the image is redrawn.


B<HIDE SEQUENCE OR EXTRA (Alt-Right-Double-Click)> This will hide sequences from view (i.e it will disappear from view). To unhide sequences you must use the pre-filter in the filter options.  (For which I should add a command line!).




=head1 APPENDIX A: LIST OF VALID OPTIONS WITH INTERNAL DEFAULTS

=over 3

=item * C<alignment_col>=E<gt>C<0>S<     >[integer] column for the first query (first sequence) in a parsed pairwise alignment. Blank/zero hides option from popup menu. The sequence will contain dashes for gaps.  

=item * C<alignment_col2>=E<gt>C<0>S<     >[integer] column for subject (second pairwise position) sequence alignment. Blank/zero hides option from popup menu. The sequence will contain dashes for gaps.  

=item * C<alignment_wrap>=E<gt>C<50>S<     >[integer] line width in aligned characters (bases/amino acids/dashes) for displaying any alignments  

=item * C<arrangeseq>=E<gt>C<oneperline>S<     >[oneperline|sameline|file] determines the arrangement of sequences that are currently being shown with -showseq. Choices: oneperline = each sequence placed on a separate line; sameline = sequences are place one after the other on the same line; file = load a file with exact positions in terms of line number and base position within the colorsub_hitcond_tests variable  

=item * C<arrangesub>=E<gt>C<stagger>S<     >[stagger|oneperline|subscaleC|subscaleN] basic  

=item * C<arrangesub_stagger_spacing>=E<gt>C<40000>S<     >[integer] bases of spacing between for sequences placed on the same sub line.  Sequences separated by less than this distance from each other will be placed on separate sub lines). This option is useful for providing space for a label.  

=item * C<canvas_bpwidth>=E<gt>C<250000>S<     >[integer] number of bases that the width of screen represents (not including indentations).  This is the  number of bases per line across  

=item * C<canvas_indent_left>=E<gt>C<60>S<     >[integer] pixels to indent from the left-side of screen window image (the drawing areas is the canvas in Tk)  

=item * C<canvas_indent_right>=E<gt>C<30>S<     >[integer] pixels to indent from the right-side of screen window image  (the drawing areas is the canvas in Tk)  

=item * C<canvas_indent_top>=E<gt>C<40>S<     >[integer] pixels to indent from top of screen before drawing sequence lines (it does not take into account graphs or extras)  (the drawing areas is the canvas in Tk)  

=item * C<color>=E<gt>C<None>S<     >(not implemented yet)  

=item * C< colorsub>=E<gt>C<None>S<     >[NONE|RESET|hitrandom|seqrandom|hitconditional] Choices for coloring subs: NONE=no coloring routines; RESET=clear all assigned colors to pairwise; hitrandom=randomly color each hit/pairwise a different color; seqrandom=randomly color each defined seequence; hitconditional=color each hit based on pseudo-perl if than statements found in the variable  

=item * C< colorsub_hitcond_col >=E<gt>C<34>S<     >[integer] column against which to test conditional statements in pairwise data (does not work on extra items or graphs)  

=item * C< colorsub_hitcond_tests>=E<gt>C<red if E<lt>2; orange if E<lt>0.99; yellow if E<lt>0.98; green if E<lt>0.97; blue if E<lt>0.96; purple if E<lt>0.95; brown if E<lt>0.94; grey if E<lt>0.93; black if E<lt>0.92; pink if E<lt>0.91>S<     >[fake code] conditional statements to color pairwise hits based on the values in the column colorsub_hitcond_col (format for tests: color [= or E<lt> or E<gt>] value; )  

=item * C<execute>=E<gt>C<>S<     >[external system command] to execute on Control-Shift-Click Left Button  

=item * C<execute2>=E<gt>C<>S<     >[external system command] to execute on Control-Shift-Click Middle Button  

=item * C<execute2_array>=E<gt>C<m>S<     >[e|m] extra or pairwise array to use in execute2 command  

=item * C<execute2_desc>=E<gt>C<>S<     >[text] description to display in right-click menu for execute2 command

=item * C<execute3>=E<gt>C<>S<     >[external system command] to execute on Control-Shift-Click Right Button  

=item * C<execute3_array>=E<gt>C<m>S<     >[e|m] extra or pairwise array to use in execute3 command  

=item * C<execute3_desc>=E<gt>C<widget>S<     >[text] description to display in right-click menu for execute3 command

=item * C<execute4>=E<gt>C<>S<     >[external system command] to execute from within right-click menu only  

=item * C<execute4_array>=E<gt>C<m>S<     >[e|m] extra or pairwise array to use in execute4 command  

=item * C<execute4_desc>=E<gt>C<>S<     >[text] description to display in right-click menu for execute command

=item * C<execute_array>=E<gt>C<e>S<     >[e|m] extra or pairwise array to use in execute command  

=item * C<execute_desc>=E<gt>C<>S<     >[text] description to display in right-click menu for execute command

=item * C<extra_arrow_diag>=E<gt>C<5>S<     >[integer] distance from point of arrow to wing/elbow of arrow  

=item * C<extra_arrow_on>=E<gt>C<1>S<     >[0|1] toggles arrows for extras off and on  

=item * C<extra_arrow_para>=E<gt>C<5>S<     >[integer] pixel distance from point of arrow along the line  

=item * C<extra_arrow_perp>=E<gt>C<4>S<     >[integer] pixel distance from base on line to wing of arrow  

=item * C<extra_color>=E<gt>C<purple>S<     >[color] default for extra object  

=item * C<extra_label_col>=E<gt>C<10>S<     >[integer] column to take values to use for extra labels  

=item * C<extra_label_col_pattern>=E<gt>C<>S<     >[regular expression] pattern to match (and extract via parentheses) replacing current value.  Allows display of only part of the data found in a column.  

=item * C<extra_label_color>=E<gt>C<purple>S<     >[color] default for the labels of extra objects  

=item * C<extra_label_fontsize>=E<gt>C<6>S<     >[integer] font size (in points) the labels of extra objects  

=item * C<extra_label_offset>=E<gt>C<2>S<     >[integer] horizontal offset for extra labels (left is negative, right is positive)

=item * C<extra_label_on>=E<gt>C<1>S<     >[0|1] toggles the text label for extra objects off and on  

=item * C<extra_label_test_col>=E<gt>C<>S<     >[integer] column to test for a pattern--if pattern matched then extra not drawn  

=item * C<extra_label_test_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in order to NOT draw the matching extra object  

=item * C<extra_offset>=E<gt>C<-4>S<     >[integer] default vertical offset of extra object (negative = up; positive = down)  

=item * C<extra_on>=E<gt>C<1>S<     >[0|1] toggles off and on the display of all extras  

=item * C<extra_width>=E<gt>C<6>S<     >[integer] default width (horizontal thickness of extra object  

=item * C<fasta_blastdb>=E<gt>C<htg:nt>S<     >[database names] for sequence fastacmd lookups 

=item * C< fasta_directory>=E<gt>C<.:fastax>S<     >[directories] to search for fasta files corresponding to sequence names in order to extract subsequences on command (names of files must be same as names of sequences)  

=item * C<fasta_fragsize>=E<gt>C<400000>S<     >[integer] fragment size for sequences in fasta directory.  Useful for quick lookups in long sequences like chromosomes.  If this is non-zero than fragments of files are searched for in the fasta_directory (nomenclature of fragmented files end with _###, e.g. chr1_000, chr1_001, etc.)  

=item * C<fasta_on>=E<gt>C<1>S<     >[0|1]  off|on turns fasta extraction on and off

=item * C<fasta_wrap>=E<gt>C<50>S<     >[integer] line width in characters for fasta files created

=item * C<filename_color>=E<gt>C<grey>S<     >[color] of text label for the filename  

=item * C<filename_offset>=E<gt>C<-10>S<     >[integer] vertical offset of text label for filename (up is negative, down is positive)  

=item * C<filename_offset_h>=E<gt>C<0>S<     >[integer] horizontal offset of text label for filename (left is negative, right is positive)  

=item * C<filename_on>=E<gt>C<1>S<     >[0|1] toggle off and on display of designated filename/parasight name (initially defined by -in if empty)  

=item * C<filename_pattern>=E<gt>C<>S<     >[regular expression]  pattern to match in the filename.  Useful for removing the path.  (Although if using graphical interface, it is easier to change the filename.)   

=item * C<filename_size>=E<gt>C<10>S<     >[integer] point size of text label shown for the filename  

=item * C<filter1_col>=E<gt>C<>S<     >[integer]  column that contains data with which to filter pairwise 

=item * C<filter1_max>=E<gt>C<>S<     >[float] limit for value in filter1_col above which pairwise are NOT drawn  

=item * C<filter1_min>=E<gt>C<>S<     >[float] limit for value in filter1_col below which pairwise are NOT drawn  

=item * C<filter2_col>=E<gt>C<>S<     >[integer] column that contains data with which to filter pairwise 

=item * C<filter2_max>=E<gt>C<>S<     >[float] limit for value in filter2_col above which pairwise are NOT drawn  

=item * C<filter2_min>=E<gt>C<>S<     >[float] limit for value in filter2_col  below which pairwise are NOT drawn  

=item * C<filterextra1_col>=E<gt>C<>S<     >[integer]  column number that contains data with which to filter extra objects--sequences are removed before arrange functions are executed 

=item * C<filterextra1_max>=E<gt>C<>S<     >[float] limit for value in filterextra1_col above which extras are NOT drawn  

=item * C<filterextra1_min>=E<gt>C<>S<     >[float] limit for value in filterextra1_col below which extras are NOT drawn  

=item * C<filterextra2_col>=E<gt>C<>S<     >[integer] column number that contains data with which to filter extras 

=item * C<filterextra2_max>=E<gt>C<>S<     >[float] limit for value in filterextra2_col above which extras are NOT drawn  

=item * C<filterextra2_min>=E<gt>C<>S<     >[float] limit for value in filterextra2_col below which extras are NOT drawn  

=item * C<filterpre1_col>=E<gt>C<>S<     >[integer] column that contains data with which to prefilter pairwise--prefiltering removes pairwise before any arranging  (normal filtering removes pairwise after filtering)  

=item * C<filterpre1_max>=E<gt>C<>S<     >[float] limit for value in filterpre1_col above which pairwise are NOT drawn or arranged  

=item * C<filterpre1_min>=E<gt>C<>S<     >[float] limit for value in filterpre1_col below which pairwise are NOT drawn or arranged  

=item * C<filterpre2_col>=E<gt>C<>S<     >[integer] column that contains data with which to prefilter pairwise—prefilter removes pairwise before any arranging  

=item * C<filterpre2_max>=E<gt>C<>S<     >[float] limit for value in filterpre2_col above which pairwise are NOT drawn or arranged  

=item * C<filterpre2_min>=E<gt>C<>S<     >[float] limit for value in filterpre2_col below which pairwise are NOT drawn or arranged  

=item * C<gif_anchor>=E<gt>C<center>S<     >[center|nw|ne|sw|se|e|w|n] positioning of background gif relative to draw point gif_x and gif_y  

=item * C<gif_on>=E<gt>C<0>S<     >displays a gif image in background (the image will not print out in postscript)  

=item * C<gif_path>=E<gt>C<>S<     >[file path] of gif image to display in background--image does not make it into the Postscript file file  

=item * C<gif_x>=E<gt>C<int($opt{window_width}/2)>S<     >[integer] background picture pixel x coordinate position (top of image is zero)  

=item * C<gif_y>=E<gt>C<0>S<     >[integer] background gif y coordinate position (0 is top of screen)  

=item * C<graph1_label_color>=E<gt>C<blue>S<     >[color] for graph1 labels (left side axis)  

=item * C<graph1_label_decimal>=E<gt>C<2>S<     >[integer] number of decimal points to round graph1 labels (left side axis)  

=item * C<graph1_label_fontsize>=E<gt>C<10>S<     >[integer] point size of graph1 labels (left side axis)  

=item * C<graph1_label_multiplier>=E<gt>C<1>S<     >[float] multiplier for graph1 labels (left side axis)  

=item * C<graph1_label_offset>=E<gt>C<1>S<     >[integer] horizontal offset for graph1 labels (left side axis)  

=item * C<graph1_label_on>=E<gt>C<1>S<     >[0|1] toggles on labels for graph1 scale (left side axis)  

=item * C<graph1_line_color>=E<gt>C<blue>S<     >[color] for graph1 connecting lines  

=item * C<graph1_line_on>=E<gt>C<1>S<     >[0|1] toggles graph1 connecting line off and on  

=item * C<graph1_line_smooth>=E<gt>C<0>S<     >[0|1] toggles on and off smoothing function for connecting line

=item * C<graph1_line_width>=E<gt>C<1>S<     >[integer] width for graph1 connecting line  

=item * C<graph1_max>=E<gt>C<100>S<     >[integer] maximum value of graph1 scale  

=item * C<graph1_min>=E<gt>C<-5>S<     >[integer] minimum value of graph1 scale  

=item * C<graph1_on>=E<gt>C<0>S<     >[0|1] toggles off and on graph1  

=item * C<graph1_point_fill_color>=E<gt>C<blue>S<     >[color] to fill points with for graph1  

=item * C<graph1_point_on>=E<gt>C<1>S<     >[0|1] toggles point drawing on and off for graph1  

=item * C<graph1_point_outline_color>=E<gt>C<blue>S<     >[color] to outline point with for graph1  

=item * C<graph1_point_outline_width>=E<gt>C<1>S<     >[integer] thickness of point outline for graph1  

=item * C<graph1_point_size>=E<gt>C<2>S<     >[integer] pixel radius size for drawing graph1 points  

=item * C<graph1_tick_color>=E<gt>C<black>S<     >[color] of tick marks for graph1 scale  

=item * C<graph1_tick_length>=E<gt>C<6>S<     >[integer] length of tick marks for graph1 scale  

=item * C<graph1_tick_offset>=E<gt>C<1>S<     >[integer] horizontal offset of tick marks for graph1 scale  

=item * C<graph1_tick_on>=E<gt>C<1>S<     >[0|1] toggles tick marks for graph1 scale off and on  

=item * C<graph1_tick_width>=E<gt>C<3>S<     >[integer] thickness of tick marks for graph1 scale  

=item * C<graph1_vline_color>=E<gt>C<black>S<     >[color] of vertical line for graph1 scale on left  

=item * C<graph1_vline_on>=E<gt>C<1>S<     >[0|1} toggles on and off vertical line for graph1 scale on left  

=item * C<graph1_vline_width>=E<gt>C<2>S<     >[integer] vertical line width for graph1 scale on left  

=item * C<graph2_label_color>=E<gt>C<red>S<     >[color] of graph2 scale labels  

=item * C<graph2_label_decimal>=E<gt>C<2>S<     >[integer] number of decimal point to round graph2 scale label  

=item * C<graph2_label_fontsize>=E<gt>C<10>S<     >[integer] point size of graph2 scale labels  

=item * C<graph2_label_multiplier>=E<gt>C<1>S<     >[float] graph2 scale label multiplier  

=item * C<graph2_label_offset>=E<gt>C<8>S<     >[integer] horizontal offset of graph2 scale labels  

=item * C<graph2_label_on>=E<gt>C<1>S<     >[0|1] toggles graph2 scale labels off and n  

=item * C<graph2_line_color>=E<gt>C<red>S<     >[color] of graph2 connecting lines  

=item * C<graph2_line_on>=E<gt>C<1>S<     >[0|1] toggles graph2 connecting lines off and on  

=item * C<graph2_line_smooth>=E<gt>C<0>S<     >[0|1] toggles graph2 connecting line smoothing off and on

=item * C<graph2_line_width>=E<gt>C<1>S<     >[integer] thickness of graph2 connecting lines  

=item * C<graph2_max>=E<gt>C<1000>S<     >[integer] maximum value for graph2 scale  

=item * C<graph2_min>=E<gt>C<-1000>S<     >[integer] minimum value for graph2 scale  

=item * C<graph2_on>=E<gt>C<0>S<     >[0|1] toggles graph2_on  

=item * C<graph2_point_fill_color>=E<gt>C<red>S<     >[color] of interior of graph2 points  

=item * C<graph2_point_on>=E<gt>C<1>S<     >[0|1] toggles graph2 point drawing on and off  

=item * C<graph2_point_outline_color>=E<gt>C<red>S<     >[color] of graph2 point outline  

=item * C<graph2_point_outline_width>=E<gt>C<1>S<     >[integer] thickness of graph2 point outline  

=item * C<graph2_point_size>=E<gt>C<2>S<     >[integer] radius size of graph 2 points  

=item * C<graph2_tick_color>=E<gt>C<black>S<     >[color] of graph2 vertical scale ticks  

=item * C<graph2_tick_length>=E<gt>C<6>S<     >[integer] length of graph2 vertical scale ticks  

=item * C<graph2_tick_offset>=E<gt>C<5>S<     >[integer] horizontal offset of graph2 vertical scale ticks  

=item * C<graph2_tick_on>=E<gt>C<1>S<     >[0|1] toggles graph2 vertical scale ticks on and off  

=item * C<graph2_tick_width>=E<gt>C<3>S<     >[integer] thickness of graph2 vertical scale ticks  

=item * C<graph2_vline_color>=E<gt>C<black>S<     >[color] of graph2 vertical scale line  

=item * C<graph2_vline_on>=E<gt>C<1>S<     >[0|1] toggles graph2 vertical scale line off and on  

=item * C<graph2_vline_width>=E<gt>C<2>S<     >[integer] thickness of graph2 vertical scale line  

=item * C<graph_scale_height>=E<gt>C<80>S<     >[integer] pixel height of shared graph scale  

=item * C<graph_scale_hline_color>=E<gt>C<black>S<     >[color] of horizontal shared graph scale lines  

=item * C<graph_scale_hline_on>=E<gt>C<1>S<     >[0|1] toggles off and on the shared horizontal interval lines of the graph scales  

=item * C<graph_scale_hline_width>=E<gt>C<1>S<     >[integer] width of shared horizontal shared graph scale lines  

=item * C<graph_scale_indent>=E<gt>C<-20>S<     >[integer] indentation for placing gscale above (or even below) the sequence line  

=item * C<graph_scale_interval>=E<gt>C<4>S<     >[integer] number of intervals  

=item * C<graph_scale_on>=E<gt>C<0>S<     >[0|1] toggles off and on the graph scales  

=item * C<help_on>=E<gt>C<1>S<     >[0|1] toggles off and on the popup help messages

=item * C<help_wrap>=E<gt>C<50>S<     >[integer] line width in characters for popup help menus

=item * C<mark_advanced>=E<gt>C<     >code for an advanced marking algorithm. Allowing for more complex searches.  Data foreach pair or extra is accessed using an array reference \$c.  Therefore to access column 4 \$\$c[4] would work. 

=item * C<mark_array>=E<gt>C<m>S<     >[e|m] default array to search (m is alignment/e is extra)(m is historical)  

=item * C<mark_col>=E<gt>C<>S<     >[integer] column to search for given pattern in order to mark matches with a color  

=item * C<mark_col2>=E<gt>C<>S<     >[integer] second column to search for pattern in order to mark matches with a color  

=item * C<mark_color>=E<gt>C<red>S<     >[color] to  mark objects with  

=item * C<mark_pairs>=E<gt>C<0>S<     >[0|1] toggles the coloring/marking of sub(jects) off and on  

=item * C<mark_pattern>=E<gt>C<AC002038>S<     >[regular expression] pattern to search for with mark/find button  

=item * C<mark_permanent>=E<gt>C<0>S<     >[0|1] toggles on and off changing the color of objects permanently (if not permanent then on redraw colors will be erased  

=item * C<mark_subs>=E<gt>C<1>S<     >[0|1] toggles the coloring/marking of sub(jects) off and on  

=item * C<pair_inter_color>=E<gt>C<red>S<     >[color] default of inter pairwise and connecting lines  

=item * C<pair_inter_line_on>=E<gt>C<0>S<     >[0|1] toggles off and on the connecting lines between inter pairwise alignments  

=item * C<pair_inter_offset>=E<gt>C<0>S<     >[integer] default offset from sequence line of inter pairwise (up is negative, down is positive)  

=item * C<pair_inter_on>=E<gt>C<1>S<     >[0|1] toggles off and on the inter pairwise alignments normally drawn on top of sequence line  

=item * C<pair_inter_width>=E<gt>C<13>S<     >[integer] width of inter pairwise  

=item * C<pair_intra_color>=E<gt>C<blue>S<     >[color] default of intra pairwise and connecting lines  

=item * C<pair_intra_line_on>=E<gt>C<0>S<     >[0|1] toggles connecting lines between  intra pairwise off and on  

=item * C<pair_intra_offset>=E<gt>C<0>S<     >[integer] default offset from seuqence  

=item * C<pair_intra_on>=E<gt>C<1>S<     >[0|1] toggles off and on the intra pairwise  

=item * C<pair_intra_width>=E<gt>C<9>S<     >[integer] width of intra pairwise  

=item * C<pair_level>=E<gt>C<NONE>S<     >[NONE|inter_over_intra|intra_over_inter] determines which pairwise type appears above the other--NONE leaves the appearance to the order of the pairwise in the inputted alignment or parasight.psa table  

=item * C<pair_type_col>=E<gt>C<>S<     >[integer] column number to determine pairwise type for sequence 1, which is checked against sequence 2.  If match then intra if no match then inter.  (Useful on sequence names that contain chromosome assignment.)  

=item * C<pair_type_col2>=E<gt>C<>S<     >[integer] column to determine pairwise type for sequence 2 in row   

=item * C<pair_type_col2_pattern>=E<gt>C<>S<     >[regular expression] to extract pairwise type determing value with parentheses  

=item * C<pair_type_col_pattern>=E<gt>C<>S<     >[regular expression] to extract pairwise type determining value with parentheses  

=item * C<popup_format>=E<gt>C<text>S<     >[text|number] determines whether column numbers or text headers are shown in popup window  

=item * C<popup_max_len>=E<gt>C<300>S<     >[integer] character length for fields in the popup menu (allows long definitions or sequences be excluded)  

=item * C<print_command>=E<gt>C<lpr -P Rainbow {}>S<     >[string] print command with brackets {} representing file name.  This is a system command executed to drive a printer.  I have never been able to get DOS to work.  This is setup for Unix on our system.  Rainbow is our color printer name.  It will fail in MSWin  

=item * C<print_multipages_high>=E<gt>C<1>S<     >[integer] height in number of pages for the print/postscript all command  

=item * C<print_multipages_wide>=E<gt>C<1>S<     >[integer] width in number of pages for print/postscript all command  

=item * C<printer_page_length>=E<gt>C<11i>S<     >[special] physical page length (longest dimension of paper) in inches for printer (requires number followed by units with i=inches or c=cm)  

=item * C<printer_page_orientation>=E<gt>C<1>S<     >[0|1] toggles printer page orientation (1=landscape 0=portrait)  

=item * C<printer_page_width>=E<gt>C<8i>S<     >[special] physical page width in inches for printer (requires number followed by units i=inches or c=cm)  

=item * C<quick_color>=E<gt>C<purple>S<     >[color] for the quick color function Shift-Button3 and Shift-Double Click Button3  

=item * C<seq_color>=E<gt>C<black>S<     >[color] of sequence (All sequences take this color.  There is currently no way to color sequences individually.)  

=item * C<seq_label_color>=E<gt>C<black>S<     >[color] of sequence name text  

=item * C<seq_label_fontsize>=E<gt>C<12>S<     >[integer] font size (in points) for all sequence names  

=item * C<seq_label_offset>=E<gt>C<-4>S<     >[integer] vertical offset of sequence names (up is negative, down is positive)  

=item * C<seq_label_offset_h>=E<gt>C<0>S<     >[integer] horizontal offset of sequence names  

=item * C<seq_label_on>=E<gt>C<1>S<     >[0|1] toggles off and on the display of sequence name labels  

=item * C<seq_label_pattern>=E<gt>C<>S<     >[regular expression] to match in sequence name for display purposes--parentheses must be used to denote the part of match to display  

=item * C<seq_line_spacing_btwn>=E<gt>C<250>S<     >[integer] pixels to separate  sequence lines from each other (roughly equivalent to spacing between text paragraphs if you consider a wrapping line of sequences to be a paragraph)  

=item * C<seq_line_spacing_wrap>=E<gt>C<200>S<     >[integer] pixels to space between a wrapping line of sequences (roughly equivaelent to spacing between the lines within a text paragraph)  

=item * C<seq_spacing_btwn_sequences>=E<gt>C<10000>S<     >[integer] bases to separate sequences drawn within the same line (roughly equivalent to spacing between words of a text paragraph)  

=item * C<seq_tick_b_color>=E<gt>C<black>S<     >[color] for begin tick marks  

=item * C<seq_tick_b_label_anchor>=E<gt>C<ne>S<     >[center|n|w|s|e|nw|ne|sw|se] anchor point for begin tick mark labels  

=item * C<seq_tick_b_label_color>=E<gt>C<black>S<     >[valid color] of tick mark label at the beginning of sequence  

=item * C<seq_tick_b_label_fontsize>=E<gt>C<9>S<     >[integer] font size (in points) for label at beginning of sequence  

=item * C<seq_tick_b_label_multiplier>=E<gt>C<0.001>S<     >[float] scaling factor for begin tick mark labels  

=item * C<seq_tick_b_label_offset>=E<gt>C<2>S<     >[integer] vertical offset for begin tick mark label  

=item * C<seq_tick_b_label_offset_h>=E<gt>C<0>S<     >[integer] horizontal offset for begin tick mark labels  

=item * C<seq_tick_b_label_on>=E<gt>C<1>S<     >[0|1] toggles  off and on the beginning tick mark labels  

=item * C<seq_tick_b_length>=E<gt>C<10>S<     >[integer] length of begin tick marks  

=item * C<seq_tick_b_offset>=E<gt>C<0>S<     >[integer] vertical offset for begin tick marks  

=item * C<seq_tick_b_on>=E<gt>C<1>S<     >[0|1] toggles off and on the begin tick marks  

=item * C<seq_tick_b_width>=E<gt>C<2>S<     >[integer] width of begin tick marks  

=item * C<seq_tick_bp>=E<gt>C<20000>S<     >[integer] tick mark interval  

=item * C<seq_tick_color>=E<gt>C<black>S<     >[color] of interval tick marks  

=item * C<seq_tick_e_color>=E<gt>C<black>S<     >[valid color] for end tick marks  

=item * C<seq_tick_e_label_anchor>=E<gt>C<nw>S<     >[center|n|w|s|e|nw|ne|se|sw] anchor point for end tick mark labels  

=item * C<seq_tick_e_label_color>=E<gt>C<black>S<     >[valid color] for end tick mark labels  

=item * C<seq_tick_e_label_fontsize>=E<gt>C<9>S<     >[integer] font size (in points)  for end tick mark labels  

=item * C<seq_tick_e_label_multiplier>=E<gt>C<0.001>S<     >[float] scaling factor for end tick mark labels  

=item * C<seq_tick_e_label_offset>=E<gt>C<2>S<     >[integer] vertical offset for end tick mark labels  

=item * C<seq_tick_e_label_offset_h>=E<gt>C<0>S<     >[integer] horizontal offset for end tick mark labels  

=item * C<seq_tick_e_label_on>=E<gt>C<1>S<     >[0|1] toggles end tick labels off and on  

=item * C<seq_tick_e_length>=E<gt>C<10>S<     >[integer] length of end tick marks  

=item * C<seq_tick_e_offset>=E<gt>C<0>S<     >[integer] vertical offset for ending tick marks  

=item * C<seq_tick_e_on>=E<gt>C<1>S<     >[0|1] toggles off and on the ending tick marks  

=item * C<seq_tick_e_width>=E<gt>C<2>S<     >[integer] width of end tick marks  

=item * C<seq_tick_label_anchor>=E<gt>C<n>S<     >[center|n|s|w|e|nw|sw|ne|se] anchor of text from tick mark draw point  

=item * C<seq_tick_label_color>=E<gt>C<black>S<     >[color] for interval tick mark  

=item * C<seq_tick_label_fontsize>=E<gt>C<9>S<     >[integer] font size (in points) for interval tick mark label  

=item * C<seq_tick_label_multiplier>=E<gt>C<0.001>S<     >[float] scaling factor for the interval tick label  

=item * C<seq_tick_label_offset>=E<gt>C<2>S<     >[integer] vertical offset of sequence interval tick mark labels  

=item * C<seq_tick_label_on>=E<gt>C<1>S<     >[0|1] toggles off and on the interval tick labels  

=item * C<seq_tick_length>=E<gt>C<10>S<     >[integer] length of interval tick marks  

=item * C<seq_tick_offset>=E<gt>C<0>S<     >[integer] vertical offset for interval tick marks  

=item * C<seq_tick_on>=E<gt>C<1>S<     >[0|1] toggles off and on the interval sequence tick marks  

=item * C<seq_tick_whole>=E<gt>C<0>S<     >[0|1] toggles whether numbering is for each individual sequence (0) or continious across multiple accession on same line (useful when analyzing chromosomes in multiple fragments) 

=item * C<seq_tick_width>=E<gt>C<2>S<     >[integer] width of interval tick marks  

=item * C<seq_width>=E<gt>C<3>S<     >[integer] width of sequence line  

=item * C<showqueryonly>=E<gt>C<0>S<     >[0|1] toggles the display of just the first sequence in a pairwise data (i.e.first column in an alignment file).  For most parsing this is equivalent to the Blast query position  

=item * C<sub_arrow_diag>=E<gt>C<5>S<     >[integer] distance between arrow point to wing/edge of arrow  

=item * C<sub_arrow_on>=E<gt>C<0>S<     >[0|1] toggles off and on the directional/orientation arrows for subjects  

=item * C<sub_arrow_paral>=E<gt>C<5>S<     >[integer]  distance between arrow point to base of arrow  

=item * C<sub_arrow_perp>=E<gt>C<4>S<     >[integer] distance from base end to wing tip of arrow   

=item * C<sub_color>=E<gt>C<lightgreen>S<     >[color] default of sub(ject) objects (all other coloring schemes over ride default)  

=item * C<sub_initoffset>=E<gt>C<30>S<     >[integer] pixel indent from top of subscales to associated sequence line (increasing pushes scales further below associated sequence)  

=item * C<sub_labelhit_col>=E<gt>C<13>S<     >[integer] column to use for labeling each hit/pairwise (label will be drawn at beginning of each hit sub)  

=item * C<sub_labelhit_color>=E<gt>C<black>S<     >color of pairwise hit label text  

=item * C<sub_labelhit_offset>=E<gt>C<0>S<     >[integer] horizontal offset for hit label  

=item * C<sub_labelhit_on>=E<gt>C<0>S<     >[0|1] turns on individual labeling of each pairwise hit  

=item * C<sub_labelhit_pattern>=E<gt>C<0?([0-9.]{4})>S<     >[regular expression] to match in data from column  

=item * C<sub_labelhit_size>=E<gt>C<9>S<     >[integer] font size (in points) for hit label  

=item * C<sub_labelseq_col>=E<gt>C<0>S<     >[integer] column to use for the beginning sub label  

=item * C<sub_labelseq_col2>=E<gt>C<4>S<     >[integer] column for second position sequence in alignment table pairwise row  

=item * C<sub_labelseq_col2_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in data from sub label sequence column 2  

=item * C<sub_labelseq_col_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in data from sub label sequence column (use parenthesis to denote data within the match to display) 

=item * C<sub_labelseq_color>=E<gt>C<black>S<     >[color] of text label for sub objects  

=item * C<sub_labelseq_offset>=E<gt>C<0>S<     >[integer] horizontal offset label  

=item * C<sub_labelseq_on>=E<gt>C<1>S<     >[0|1] toggles overall begin sequence label for sub(ject) label off and on  

=item * C<sub_labelseq_size>=E<gt>C<6>S<     >[integer] font size (in points) for begin label sequence  

=item * C<sub_labelseqe_col>=E<gt>C<4>S<     >[integer] column to use for the end subject label  

=item * C<sub_labelseqe_col2>=E<gt>C<0>S<     >[integer] column for second position in alignment table pairwise row  

=item * C<sub_labelseqe_col2_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in data from column  

=item * C<sub_labelseqe_col_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in data from column  

=item * C<sub_labelseqe_color>=E<gt>C<black>S<     >[valid color] of label text  

=item * C<sub_labelseqe_offset>=E<gt>C<0>S<     >[integer] horizontal offset for label  

=item * C<sub_labelseqe_on>=E<gt>C<0>S<     >[0|1] toggles off and on the overall sub(ject) label at end of last hit/pairwise  

=item * C<sub_labelseqe_size>=E<gt>C<6>S<     >[integer] font size (in points) for end subject label  

=item * C<sub_line_spacing>=E<gt>C<9>S<     >[integer] pixels per line determining the spacing between subs placed on different lines  

=item * C<sub_on>=E<gt>C<1>S<     >[0|1] toggles sub(ject) display off and on (these are the pairwise representations drawn below the sequence line) For BLAST searches these traditionally represent the subject sequences found in a database search.  

=item * C<sub_scale_categoric_string>=E<gt>C<>S<     >[string] list of comma delimited category names  

=item * C<sub_scale_col>=E<gt>C<>S<     >[integer] column for value to arrange pairwise hit on sub scale (subscale)  

=item * C<sub_scale_col2>=E<gt>C<>S<     >[integer] column for second position sequence in alignment pairwise (only used if defined)  

=item * C<sub_scale_col2_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in column 2  

=item * C<sub_scale_col_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in column  

=item * C<sub_scale_hline_color>=E<gt>C<grey>S<     >[valid color] for horizontal sub scale lines  

=item * C<sub_scale_hline_on>=E<gt>C<1>S<     >[0|1] toggles off and on the horizontal scale lines for sub scale  

=item * C<sub_scale_hline_width>=E<gt>C<1>S<     >[integer] width of horizontal sub scale lines  

=item * C<sub_scale_label_color>=E<gt>C<black>S<     >[color] for sub scale axis label  

=item * C<sub_scale_label_fontsize>=E<gt>C<12>S<     >[integer] font size (in points) for sub scale axis label  

=item * C<sub_scale_label_multiplier>=E<gt>C<100>S<     >[integer] multiplication factor for sub scale label  

=item * C<sub_scale_label_offset>=E<gt>C<1>S<     >[integer] horizontal offset for sub scale axis tick marks  

=item * C<sub_scale_label_on>=E<gt>C<1>S<     >[0|1] toggles off and on sub scale axis tick mark labels  

=item * C<sub_scale_label_pattern>=E<gt>C<>S<     >[regular expression] pattern to match in sub scale label  

=item * C<sub_scale_lines>=E<gt>C<10>S<     >[integer] number of lines (or interval steps) to plot for stagger or cscale (automatically set for subscaleC)  

=item * C<sub_scale_max>=E<gt>C<1.00>S<     >[float] maximum value to place on the sub scale (automatically set for subscaleC)  

=item * C<sub_scale_min>=E<gt>C<0.80>S<     >[float] minimum value to place on the sub scale (automatically set for subscaleC)  

=item * C<sub_scale_on>=E<gt>C<0>S<     >[0|1] toggles sub scale on and off  

=item * C<sub_scale_step>=E<gt>C<0.01>S<     >[float] value to increment between each step (automatically set to -1 for subscaleC, 1 reverses subscaleC)  

=item * C<sub_scale_tick_color>=E<gt>C<black>S<     >[color] for sub scale axis tick marks  

=item * C<sub_scale_tick_length>=E<gt>C<9>S<     >[integer] length of sub axis tick marks  

=item * C<sub_scale_tick_offset>=E<gt>C<4>S<     >[integer] offset of sub scale axis tick marks  

=item * C<sub_scale_tick_on>=E<gt>C<1>S<     >[0|1] toggles off and on the sub scale axis at horizontal tick positions  

=item * C<sub_scale_tick_width>=E<gt>C<3>S<     >[integer] width of sub scale axis tick marks  

=item * C<sub_scale_vline_color>=E<gt>C<black>S<     >[color] for vertical axis line of sub scale  

=item * C<sub_scale_vline_offset>=E<gt>C<-5>S<     >[integer] horizontal offset for subject axis line  

=item * C<sub_scale_vline_on>=E<gt>C<1>S<     >[0|1] toggles off and on the vertical axis line for sub scale  

=item * C<sub_scale_vline_width>=E<gt>C<2>S<     >[integer] width of sub scale axis line  

=item * C<sub_width>=E<gt>C<8>S<     >[integer] default width (thickness) of sub objects  

=item * C<template_desc_on>=E<gt>C<1>S<     >[0|1] toggles off and on wether  descriptions, such as this one, are saved in a template file with each option variable  

=item * C<text2_anchor>=E<gt>C<nw>S<     >[center|n|w|s|e|nw|ne|se|sw] anchor point for end tick mark labels  

=item * C<text2_color>=E<gt>C<red>S<     >[color] for end tick mark labels  

=item * C<text2_offset>=E<gt>C<0>S<     >[integer] vertical offset for end tick mark labels  

=item * C<text2_offset_h>=E<gt>C<0>S<     >[integer] horizontal offset for end tick mark labels  

=item * C<text2_on>=E<gt>C<1>S<     >[0|1] toggles end tick labels off and on  

=item * C<text2_size>=E<gt>C<20>S<     >[integer] font size (in points)  for end tick mark labels  

=item * C<text2_text>=E<gt>C<>S<     >[text] to display within a parasight view (useful for automation)  

=item * C<text_anchor>=E<gt>C<nw>S<     >[center|n|w|s|e|nw|ne|se|sw] anchor point for end tick mark labels  

=item * C<text_color>=E<gt>C<red>S<     >[color] for end tick mark labels  

=item * C<text_fontsize>=E<gt>C<20>S<     >[integer] font size (in points)  for end tick mark labels  

=item * C<text_offset>=E<gt>C<0>S<     >[integer] vertical offset for end tick mark labels  

=item * C<text_offset_h>=E<gt>C<0>S<     >[integer] horizontal offset for end tick mark labels  

=item * C<text_on>=E<gt>C<1>S<     >[0|1] toggles end tick labels off and on  

=item * C<text_text>=E<gt>C<>S<     >[text] to display within a parasight view (useful for automation)  

=item * C<window_font_size>=E<gt>C<9>S<     >[integer] font size for parasight in general (not implemented)  

=item * C<window_height>=E<gt>C<550>S<     >[integer] pixel height of main window on the initial start up  

=item * C<window_width>=E<gt>C<800>S<     >[integer] pixel width of the main window on the initial start up  , 

=back

=head1 APPENDIX B: QUICK REFERENCE

=head2 COMMAND LINE SUMMARY

B<-align [I<filepath1:filepath2:etc>]> load pairwise alignment table(s) (table must be miropeats format) 

B<-arrangeseq  [I<oneperline/sameline/file>]>  (default is oneperline)

 *oneperline = each sequence is placed on a separate wrapping line
 *sameline = the sequences are placed in alphabetical 
    order on the same line
 *file:filepath = arrange file that allows specification
    of line/paragraph and position 

B<-arrangesub [I<oneperline/stagger/subscale/cscale>]> (default stagger)
Arrange subs below the sequence.

 *oneperline = each sequence is placed on its own line
    underneath sequence
 *stagger = multiple subjects are placed on same line 
    only when non-overlapping
 *subscaleN =  pairwise hits are placed on a numerical scale
    based on values in chosen column(s)
 *subscaleC =  pairwise hits are placed on categorical 
    scale based on hash(s)

-color [I<scheme>] B<***not implemented yet, no demand yet***> Use other options for determining inter vs intrachromosal***

B<-colorsub [I<NONE/RESET/seqrandom/hitrandom/hitconditional>]>

 *NONE = does not add a colorsub and does not remove colors 
    for pairwise hits
 *RESET = removes colors for pairwise hits 
    colors for pairwise hits override colors for sequence hits
 *seqrandom = color all pairwise comparisons for a subject the same
 *hitrandom = randomly independently color each pairwise comparison
 *hitconditional = allows coloring based on a conditional statement

B<-extra [I<filepath1:filepath2:etc>]> loads extra sequence feature table(s) Sequence features are annotation that have single begin and end points (e.g. exons, introns, and repeats). The rows must consist of seqname[tab]begin[tab]end.  Further columns may contain optional data.  Columns named C<offset>, C<width>, and C<color> provide extra formatting information.

B<-graph1 [I<filepath1:filepath2:etc>]> Graphs a data set of values above the sequence line. such as %GC. The data scale is found on the left.  The data row format is simply seqname[TAB]begin[TAB]value.  No more, no less.  For regions with out a value a blank will cause the graph line to be disrupted.

B<-graph2 [I<filepath1:filepath2:etc>]>  Creates another graph using the scale on the right axis.  Same parameters as -graph1

B<-in [I<filepath>]>  load a previously saved parasight view. Three files required are *.psa,  *.pse and *.psm  (*.psg needed only if a graph has been used)

B<-options [I<'opt1=E<gt>value1,opt2=E<gt>value2'>]>  *** Allows all of the parasight options to be changed directly ***. One and zero are used for on/off, yes/no and true/false.  Complete access for the programmer using parasight as a displayer   (e.g. 'canvas_width=E<gt>500,seq_tick_on=E<gt>1,graph_scale_on=E<gt>1')

B<-showseq [I<a file or seqname(s):>]>  names of sequences to display

   *ALL = show all files (default) 
   *no colon = load as file of names 
     format each line ( seqname[TAB]length[TAB]begin[TAB]end )
     only sequence name is required other info optional
   *colon(:) = parse as list of colon-delimited seq names
     format: (seqname,length,begin,end:seqname2,length2,begin2,end2)

B<-showseqqueryonly> This toggles the display of only the first sequence in a given row.  This is the usually position for a blast query (hence the name of the option).

B<-showsub [I<file | seqnames: | ALL>]> names of subjects to display

   *ALL: displays all subject sequences (default)
   *no colon = load file containing names (one seqname per line)
   *colon(:) = parse input as list of colon-delimited sequence names

B<-template [filepath]> loads a saved option template file.  Template files can be stored in default directories for easy loading.

B<ADVANCED OPTIONS>

B<-minload>     
   *loads only the relevant pairwise that will be displayed 
   (quicker when just certain sequences are needed from large files)

B<-precode I<'perl code commands to execute after first screen draw'>>
   *an advanced option useful for automating initial tasks

B<-die>  parasight ends after executing precode
   *an advanced option useful in automating tasks


=head2 OPTION PRECEDENCE

B<internal default ---E<gt> -in ---E<gt> -template ---E<gt> -option ---E<gt> commandline> 


=head2 MOUSE FUNCTIONS

 [DBL]=double click  [DRAG]=button hold down and move mouse 
 EXECUTE # = Execute Command (User Defined under MISC options)
 
 KEY            LEFT-BUTTON       MIDDLE BUTTON  RIGHT-CLICK                             
 ---------      -----------       -------------  --------------------
 NONE           Popup Desc                       Menu
 CONTROL        Zoom in                          Zoom out
 SHIFT          Move Object[DRAG]                Quick color; Uncolor[DBL]
 ALTERNATE      Del  Object [DBL]                Lower Object; Raise Object[DBL]
 CONTROL-SHIFT  Execute 1         Execute 2      Execute 3


=head2 COMPACT ALPHABETICAL LIST OF -OPTIONS WITH DEFAULTS

B<alignment_col>=E<gt>C<0>S<  >
B<alignment_col2>=E<gt>C<0>S<  >
B<alignment_wrap>=E<gt>C<50>S<  >
B<arrangeseq>=E<gt>C<oneperline>S<  >
B<arrangesub>=E<gt>C<stagger>S<  >
B<arrangesub_stagger_spacing>=E<gt>C<40000>S<  >
B<canvas_bpwidth>=E<gt>C<250000>S<  >
B<canvas_indent_left>=E<gt>C<60>S<  >
B<canvas_indent_right>=E<gt>C<30>S<  >
B<canvas_indent_top>=E<gt>C<40>S<  >
B<color>=E<gt>C< None>S<  >
B<colorsub>=E<gt>C< None>S<  >
B<colorsub_hitcond_col>=E<gt>C<34>S<  >
B<colorsub_hitcond_tests>=E<gt>C<red if E<lt>2; orange if E<lt>0.99; yellow if E<lt>0.98; green if E<lt>0.97; blue if E<lt>0.96; purple if E<lt>0.95; brown if E<lt>0.94; grey if E<lt>0.93; black if E<lt>0.92; pink if E<lt>0.91>S<  >
B<execute>=E<gt>C<>S<  >
B<execute2>=E<gt>C<>S<  >
B<execute2_array>=E<gt>C<m>S<  >
B<execute2_desc>=E<gt>C<>S<  >
B<execute3>=E<gt>C<>S<  >
B<execute3_array>=E<gt>C<m>S<  >
B<execute3_desc>=E<gt>C<widget>S<  >
B<execute4>=E<gt>C<>S<  >
B<execute4_array>=E<gt>C<m>S<  >
B<execute4_desc>=E<gt>C<>S<  >
B<execute_array>=E<gt>C<e>S<  >
B<execute_desc>=E<gt>C<>S<  >
B<extra_arrow_diag>=E<gt>C<5>S<  >
B<extra_arrow_on>=E<gt>C<1>S<  >
B<extra_arrow_para>=E<gt>C<5>S<  >
B<extra_arrow_perp>=E<gt>C<4>S<  >
B<extra_color>=E<gt>C<purple>S<  >
B<extra_label_col>=E<gt>C<10>S<  >
B<extra_label_col_pattern>=E<gt>C<>S<  >
B<extra_label_color>=E<gt>C<purple>S<  >
B<extra_label_fontsize>=E<gt>C<6>S<  >
B<extra_label_offset>=E<gt>C<2>S<  >
B<extra_label_on>=E<gt>C<1>S<  >
B<extra_label_test_col>=E<gt>C<>S<  >
B<extra_label_test_pattern>=E<gt>C<>S<  >
B<extra_offset>=E<gt>C<-4>S<  >
B<extra_on>=E<gt>C<1>S<  >
B<extra_width>=E<gt>C<6>S<  >
B<fasta_blastdb>=E<gt>C<htg:nt>S<  >
B<fasta_directory>=E<gt>C<.:fastax>S<  >
B<fasta_fragsize>=E<gt>C<400000>S<  >
B<fasta_on>=E<gt>C<1>S<  >
B<fasta_wrap>=E<gt>C<50>S<  >
B<filename_color>=E<gt>C<grey>S<  >
B<filename_offset>=E<gt>C<-10>S<  >
B<filename_offset_h>=E<gt>C<0>S<  >
B<filename_on>=E<gt>C<1>S<  >
B<filename_pattern>=E<gt>C<>S<  >
B<filename_size>=E<gt>C<10>S<  >
B<filter1_col>=E<gt>C<>S<  >
B<filter1_max>=E<gt>C<>S<  >
B<filter1_min>=E<gt>C<>S<  >
B<filter2_col>=E<gt>C<>S<  >
B<filter2_max>=E<gt>C<>S<  >
B<filter2_min>=E<gt>C<>S<  >
B<filterextra1_col>=E<gt>C<>S<  >
B<filterextra1_max>=E<gt>C<>S<  >
B<filterextra1_min>=E<gt>C<>S<  >
B<filterextra2_col>=E<gt>C<>S<  >
B<filterextra2_max>=E<gt>C<>S<  >
B<filterextra2_min>=E<gt>C<>S<  >
B<filterpre1_col>=E<gt>C<>S<  >
B<filterpre1_max>=E<gt>C<>S<  >
B<filterpre1_min>=E<gt>C<>S<  >
B<filterpre2_col>=E<gt>C<>S<  >
B<filterpre2_max>=E<gt>C<>S<  >
B<filterpre2_min>=E<gt>C<>S<  >
B<gif_anchor>=E<gt>C<center>S<  >
B<gif_on>=E<gt>C<0>S<  >
B<gif_path>=E<gt>C<>S<  >
B<gif_x>=E<gt>C< int($opt{window_width}/2)>S<  >
B<gif_y>=E<gt>C<0>S<  >
B<graph1_label_color>=E<gt>C<blue>S<  >
B<graph1_label_decimal>=E<gt>C<2>S<  >
B<graph1_label_fontsize>=E<gt>C<10>S<  >
B<graph1_label_multiplier>=E<gt>C<1>S<  >
B<graph1_label_offset>=E<gt>C<1>S<  >
B<graph1_label_on>=E<gt>C<1>S<  >
B<graph1_line_color>=E<gt>C<blue>S<  >
B<graph1_line_on>=E<gt>C<1>S<  >
B<graph1_line_smooth>=E<gt>C<0>S<  >
B<graph1_line_width>=E<gt>C<1>S<  >
B<graph1_max>=E<gt>C<100>S<  >
B<graph1_min>=E<gt>C<-5>S<  >
B<graph1_on>=E<gt>C<0>S<  >
B<graph1_point_fill_color>=E<gt>C<blue>S<  >
B<graph1_point_on>=E<gt>C<1>S<  >
B<graph1_point_outline_color>=E<gt>C<blue>S<  >
B<graph1_point_outline_width>=E<gt>C<1>S<  >
B<graph1_point_size>=E<gt>C<2>S<  >
B<graph1_tick_color>=E<gt>C<black>S<  >
B<graph1_tick_length>=E<gt>C<6>S<  >
B<graph1_tick_offset>=E<gt>C<1>S<  >
B<graph1_tick_on>=E<gt>C<1>S<  >
B<graph1_tick_width>=E<gt>C<3>S<  >
B<graph1_vline_color>=E<gt>C<black>S<  >
B<graph1_vline_on>=E<gt>C<1>S<  >
B<graph1_vline_width>=E<gt>C<2>S<  >
B<graph2_label_color>=E<gt>C<red>S<  >
B<graph2_label_decimal>=E<gt>C<2>S<  >
B<graph2_label_fontsize>=E<gt>C<10>S<  >
B<graph2_label_multiplier>=E<gt>C<1>S<  >
B<graph2_label_offset>=E<gt>C<8>S<  >
B<graph2_label_on>=E<gt>C<1>S<  >
B<graph2_line_color>=E<gt>C<red>S<  >
B<graph2_line_on>=E<gt>C<1>S<  >
B<graph2_line_smooth>=E<gt>C<0>S<  >
B<graph2_line_width>=E<gt>C<1>S<  >
B<graph2_max>=E<gt>C<1000>S<  >
B<graph2_min>=E<gt>C<-1000>S<  >
B<graph2_on>=E<gt>C<0>S<  >
B<graph2_point_fill_color>=E<gt>C<red>S<  >
B<graph2_point_on>=E<gt>C<1>S<  >
B<graph2_point_outline_color>=E<gt>C<red>S<  >
B<graph2_point_outline_width>=E<gt>C<1>S<  >
B<graph2_point_size>=E<gt>C<2>S<  >
B<graph2_tick_color>=E<gt>C<black>S<  >
B<graph2_tick_length>=E<gt>C<6>S<  >
B<graph2_tick_offset>=E<gt>C<5>S<  >
B<graph2_tick_on>=E<gt>C<1>S<  >
B<graph2_tick_width>=E<gt>C<3>S<  >
B<graph2_vline_color>=E<gt>C<black>S<  >
B<graph2_vline_on>=E<gt>C<1>S<  >
B<graph2_vline_width>=E<gt>C<2>S<  >
B<graph_scale_height>=E<gt>C<80>S<  >
B<graph_scale_hline_color>=E<gt>C<black>S<  >
B<graph_scale_hline_on>=E<gt>C<1>S<  >
B<graph_scale_hline_width>=E<gt>C<1>S<  >
B<graph_scale_indent>=E<gt>C<-20>S<  >
B<graph_scale_interval>=E<gt>C<4>S<  >
B<graph_scale_on>=E<gt>C<0>S<  >
B<help_on>=E<gt>C<1>S<  >
B<help_wrap>=E<gt>C<50>S<  >
B<mark_advanced>=E<gt>C<>S<  >
B<mark_array>=E<gt>C<m>S<  >
B<mark_col>=E<gt>C<>S<  >
B<mark_col2>=E<gt>C<>S<  >
B<mark_color>=E<gt>C<red>S<  >
B<mark_pairs>=E<gt>C<0>S<  >
B<mark_pattern>=E<gt>C<AC002038>S<  >
B<mark_permanent>=E<gt>C<0>S<  >
B<mark_subs>=E<gt>C<1>S<  >
B<pair_inter_color>=E<gt>C<red>S<  >
B<pair_inter_line_on>=E<gt>C<0>S<  >
B<pair_inter_offset>=E<gt>C<0>S<  >
B<pair_inter_on>=E<gt>C<1>S<  >
B<pair_inter_width>=E<gt>C<13>S<  >
B<pair_intra_color>=E<gt>C<blue>S<  >
B<pair_intra_line_on>=E<gt>C<0>S<  >
B<pair_intra_offset>=E<gt>C<0>S<  >
B<pair_intra_on>=E<gt>C<1>S<  >
B<pair_intra_width>=E<gt>C<9>S<  >
B<pair_level>=E<gt>C<NONE>S<  >
B<pair_type_col>=E<gt>C<>S<  >
B<pair_type_col2>=E<gt>C<>S<  >
B<pair_type_col2_pattern>=E<gt>C<>S<  >
B<pair_type_col_pattern>=E<gt>C<>S<  >
B<popup_format>=E<gt>C<text>S<  >
B<popup_max_len>=E<gt>C<300>S<  >
B<print_command>=E<gt>C<lpr -P Rainbow {}>S<  >
B<print_multipages_high>=E<gt>C<1>S<  >
B<print_multipages_wide>=E<gt>C<1>S<  >
B<printer_page_length>=E<gt>C<11i>S<  >
B<printer_page_orientation>=E<gt>C<1>S<  >
B<printer_page_width>=E<gt>C<8i>S<  >
B<quick_color>=E<gt>C<purple>S<  >
B<seq_color>=E<gt>C<black>S<  >
B<seq_label_color>=E<gt>C<black>S<  >
B<seq_label_fontsize>=E<gt>C<12>S<  >
B<seq_label_offset>=E<gt>C<-4>S<  >
B<seq_label_offset_h>=E<gt>C<0>S<  >
B<seq_label_on>=E<gt>C<1>S<  >
B<seq_label_pattern>=E<gt>C<>S<  >
B<seq_line_spacing_btwn>=E<gt>C<250>S<  >
B<seq_line_spacing_wrap>=E<gt>C<200>S<  >
B<seq_spacing_btwn_sequences>=E<gt>C<10000>S<  >
B<seq_tick_b_color>=E<gt>C<black>S<  >
B<seq_tick_b_label_anchor>=E<gt>C<ne>S<  >
B<seq_tick_b_label_color>=E<gt>C<black>S<  >
B<seq_tick_b_label_fontsize>=E<gt>C<9>S<  >
B<seq_tick_b_label_multiplier>=E<gt>C<0.001>S<  >
B<seq_tick_b_label_offset>=E<gt>C<2>S<  >
B<seq_tick_b_label_offset_h>=E<gt>C<0>S<  >
B<seq_tick_b_label_on>=E<gt>C<1>S<  >
B<seq_tick_b_length>=E<gt>C<10>S<  >
B<seq_tick_b_offset>=E<gt>C<0>S<  >
B<seq_tick_b_on>=E<gt>C<1>S<  >
B<seq_tick_b_width>=E<gt>C<2>S<  >
B<seq_tick_bp>=E<gt>C<20000>S<  >
B<seq_tick_color>=E<gt>C<black>S<  >
B<seq_tick_e_color>=E<gt>C<black>S<  >
B<seq_tick_e_label_anchor>=E<gt>C<nw>S<  >
B<seq_tick_e_label_color>=E<gt>C<black>S<  >
B<seq_tick_e_label_fontsize>=E<gt>C<9>S<  >
B<seq_tick_e_label_multiplier>=E<gt>C<0.001>S<  >
B<seq_tick_e_label_offset>=E<gt>C<2>S<  >
B<seq_tick_e_label_offset_h>=E<gt>C<0>S<  >
B<seq_tick_e_label_on>=E<gt>C<1>S<  >
B<seq_tick_e_length>=E<gt>C<10>S<  >
B<seq_tick_e_offset>=E<gt>C<0>S<  >
B<seq_tick_e_on>=E<gt>C<1>S<  >
B<seq_tick_e_width>=E<gt>C<2>S<  >
B<seq_tick_label_anchor>=E<gt>C<n>S<  >
B<seq_tick_label_color>=E<gt>C<black>S<  >
B<seq_tick_label_fontsize>=E<gt>C<9>S<  >
B<seq_tick_label_multiplier>=E<gt>C<0.001>S<  >
B<seq_tick_label_offset>=E<gt>C<2>S<  >
B<seq_tick_label_on>=E<gt>C<1>S<  >
B<seq_tick_length>=E<gt>C<10>S<  >
B<seq_tick_offset>=E<gt>C<0>S<  >
B<seq_tick_on>=E<gt>C<1>S<  >
B<seq_tick_whole>=E<gt>C<0>S<  >
B<seq_tick_width>=E<gt>C<2>S<  >
B<seq_width>=E<gt>C<3>S<  >
B<showqueryonly>=E<gt>C<0>S<  >
B<sub_arrow_diag>=E<gt>C<5>S<  >
B<sub_arrow_on>=E<gt>C<0>S<  >
B<sub_arrow_paral>=E<gt>C<5>S<  >
B<sub_arrow_perp>=E<gt>C<4>S<  >
B<sub_color>=E<gt>C<lightgreen>S<  >
B<sub_initoffset>=E<gt>C<30>S<  >
B<sub_labelhit_col>=E<gt>C<13>S<  >
B<sub_labelhit_color>=E<gt>C<black>S<  >
B<sub_labelhit_offset>=E<gt>C<0>S<  >
B<sub_labelhit_on>=E<gt>C<0>S<  >
B<sub_labelhit_pattern>=E<gt>C<0?([0-9.]{4})>S<  >
B<sub_labelhit_size>=E<gt>C<9>S<  >
B<sub_labelseq_col>=E<gt>C<0>S<  >
B<sub_labelseq_col2>=E<gt>C<4>S<  >
B<sub_labelseq_col2_pattern>=E<gt>C<>S<  >
B<sub_labelseq_col_pattern>=E<gt>C<>S<  >
B<sub_labelseq_color>=E<gt>C<black>S<  >
B<sub_labelseq_offset>=E<gt>C<0>S<  >
B<sub_labelseq_on>=E<gt>C<1>S<  >
B<sub_labelseq_size>=E<gt>C<6>S<  >
B<sub_labelseqe_col>=E<gt>C<4>S<  >
B<sub_labelseqe_col2>=E<gt>C<0>S<  >
B<sub_labelseqe_col2_pattern>=E<gt>C<>S<  >
B<sub_labelseqe_col_pattern>=E<gt>C<>S<  >
B<sub_labelseqe_color>=E<gt>C<black>S<  >
B<sub_labelseqe_offset>=E<gt>C<0>S<  >
B<sub_labelseqe_on>=E<gt>C<0>S<  >
B<sub_labelseqe_size>=E<gt>C<6>S<  >
B<sub_line_spacing>=E<gt>C<9>S<  >
B<sub_on>=E<gt>C<1>S<  >
B<sub_scale_categoric_string>=E<gt>C<>S<  >
B<sub_scale_col>=E<gt>C<>S<  >
B<sub_scale_col2>=E<gt>C<>S<  >
B<sub_scale_col2_pattern>=E<gt>C<>S<  >
B<sub_scale_col_pattern>=E<gt>C<>S<  >
B<sub_scale_hline_color>=E<gt>C<grey>S<  >
B<sub_scale_hline_on>=E<gt>C<1>S<  >
B<sub_scale_hline_width>=E<gt>C<1>S<  >
B<sub_scale_label_color>=E<gt>C<black>S<  >
B<sub_scale_label_fontsize>=E<gt>C<12>S<  >
B<sub_scale_label_multiplier>=E<gt>C<100>S<  >
B<sub_scale_label_offset>=E<gt>C<1>S<  >
B<sub_scale_label_on>=E<gt>C<1>S<  >
B<sub_scale_label_pattern>=E<gt>C<>S<  >
B<sub_scale_lines>=E<gt>C<10>S<  >
B<sub_scale_max>=E<gt>C<1.00>S<  >
B<sub_scale_min>=E<gt>C<0.80>S<  >
B<sub_scale_on>=E<gt>C<0>S<  >
B<sub_scale_step>=E<gt>C<0.01>S<  >
B<sub_scale_tick_color>=E<gt>C<black>S<  >
B<sub_scale_tick_length>=E<gt>C<9>S<  >
B<sub_scale_tick_offset>=E<gt>C<4>S<  >
B<sub_scale_tick_on>=E<gt>C<1>S<  >
B<sub_scale_tick_width>=E<gt>C<3>S<  >
B<sub_scale_vline_color>=E<gt>C<black>S<  >
B<sub_scale_vline_offset>=E<gt>C<-5>S<  >
B<sub_scale_vline_on>=E<gt>C<1>S<  >
B<sub_scale_vline_width>=E<gt>C<2>S<  >
B<sub_width>=E<gt>C<8>S<  >
B<template_desc_on>=E<gt>C<1>S<  >
B<text2_anchor>=E<gt>C<nw>S<  >
B<text2_color>=E<gt>C<red>S<  >
B<text2_offset>=E<gt>C<0>S<  >
B<text2_offset_h>=E<gt>C<0>S<  >
B<text2_on>=E<gt>C<1>S<  >
B<text2_size>=E<gt>C<20>S<  >
B<text2_text>=E<gt>C<>S<  >
B<text_anchor>=E<gt>C<nw>S<  >
B<text_color>=E<gt>C<red>S<  >
B<text_fontsize>=E<gt>C<20>S<  >
B<text_offset>=E<gt>C<0>S<  >
B<text_offset_h>=E<gt>C<0>S<  >
B<text_on>=E<gt>C<1>S<  >
B<text_text>=E<gt>C<>S<  >
B<window_font_size>=E<gt>C<9>S<  >
B<window_height>=E<gt>C<550>S<  >
B<window_width>=E<gt>C<800>


=head1 APPENDIX C: INSTALLATION (WINDOWS OR UNIX)

Parasight has been tested extensively on Solaris, Linux, and MsWindows. Perl is available from www.perl.org. ActiveState (www.activestate.com) has binary versions available for many platforms--particularly useful for Windows installs.  Follow instructions on the choosen sites for installing Perl. Unix installs should be easier simply because you probably have more experience with Perl or you have a network administrator.  Windows installs are quite easy--just like installing any other program.  Once the install is done put parasight program in the Perl bin directory (usually C:\Perl\bin).  If you need to install any Perl modules such as Tk consult the individual OS.  For Windows Active State binary the B<PPM> provides easy searches and installations of modules.  UNIX environments can utilize the B<CPAN> module..  

If there is a strong need a standalone versionsof the program that are package together with all need Perl functions could be generated using ActiveState's PerlApp program.  All needed components are contained within the "packed up" executable for both Linux, Solaris, and Windows.  No installation of Perl is needed.  Note this is not a compiled version, so the run speed will be the same as the non-PerlApp-packaged program.  It is actually just an executable that has collected all of the Perl components required for Parasight to run.

=head1 APPENDIX D: PRECODE HINTS

Precode affords the ability to add additional code to further manipulate parasight.  Extensive use of precode is found in the F<parasight.examples> file. The best way to figure out how to manipulate parasight is to study all of the parasight code. Of course even I am trying to forget most of the code so the following are useful subroutines to abuse:

First the hash variable containing all of the command line options is %opt.  So, if you want to chance arrangesub you have to use the code $opt{'arrangesub'};

=over 5

Useful commands to use when scripting:

 $opt{'x'}

Any normal option can be accessed within the hash %opt.

 &reshowNredraw; &update;
 
These two subroutines will cause the any changes in options to be redrawn and updated on the screen.  While update is not normally used in the internal code (as it is called automatically whenever control is returned to the GUI), it is necessary when a script has control of parsight.

 &print_screen(0, "fileoutpath");

This will print a postscript of the visble screen to the designated file.  If 1 is used for the initial print varaible then the postscript will be sent to the printer. If zero is used only the file is created.

 &print_all (1, "fileoutpath");

This will print a poscript of the entire parasight area to the designated file. If 1 is used for the intial print variable then the postscript will be sent to the printer. If zero is used only the file is created.  Depending upon the multipage options, multiple files may be created.

 &save_parasight_table("basefileoutpath");

To save as parasight formated files which can be reload with the -in "basefileoutpath" name.

 &fitlongestline;

This will force the length of the screen in bases to the length of the longest sequence.  This is most useful for BLAST views.

 $opt{'die'}=0;

This is useful to turn off the die option if you are subsequently saving the parasight files.  Otherwise when you load the saved parasight it will "die" before you get to see it.

 &reshowNredraw; &update; print "PAUSED\n"; my $pause=<STDIN>; 

A useful sequence of commands if you want to pause for the user. 

 $opt{"text_text"}="This is displayed text."; $opt{"text_fontsize"}=16; 	$opt{"text_offset_h"}=10;

Allows for a line of text to be printed within the image.  text2_text allows for a second line.

=back


=head1 APPENDIX E: ADDITIONAL EXAMPLES

 parasight -showseq show.file -extra repeat.file:exon.file

=over 5

This draws the sequences specified in C<show.file> decorated with the repeats and exons specified in C<repeat.file> and C<exon.file>. Note: this example does not contain any alignments so C<show.file> is required in order to specify the lengths of the sequencesto be displayed.

=back

 parasight -in saved  -extra exons:introns 
     -arrangeseq oneperline

=over 5

This loads a saved parasight, adds extra annotation from the files C<exons> and C<introns> annotation.  It arranges subjects one per line below the sequence

=back


=head1 AUTHOR

Jeff Bailey (jab@cwru.edu)

=head1 ACKNOWLEDGEMENTS

This software was developed in the laboratory of Evan Eichler, Department of Genetics,Case Western Reserve University and University Hosiptals, Cleveland.

=head1 COPYRIGHT

Copyright (C) 2001-3 Jeff Bailey. Distribute and modify freely as defined by the GNU General Public License.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut


