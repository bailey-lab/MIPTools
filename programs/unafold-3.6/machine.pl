#! /usr/bin/perl -w

use 5.006_001;
use strict;
use warnings;

sub nullify () {
    s/[+-]?[0-9]+\.[0-9]*/0/g;
    s/[+-]?[0-9]*\.[0-9]+/0/g;
}

my $NO_WC_MISMATCHES = 0;
my $NO_GU_MISMATCHES = 0;
my $NO_GU_GU_STACKS = 0;
my $NO_GU_CLOSE = 0;
my $NO_DOUBLE_GU_CLOSE = 0;
my $NO_GU_BASEPAIRS = 0;
my $NULL_ENERGIES = 0;

for (my $i = 0; $i < @ARGV; ++$i) {
    if ($ARGV[$i] eq '--no-wc-mismatches') {
	++$NO_WC_MISMATCHES;
	splice @ARGV, $i, 1;
	--$i;
    } elsif ($ARGV[$i] eq '--no-gu-mismatches') {
	++$NO_GU_MISMATCHES;
	splice @ARGV, $i, 1;
	--$i;
    } elsif ($ARGV[$i] eq '--no-gu-gu-stacks') {
	++$NO_GU_GU_STACKS;
	splice @ARGV, $i, 1;
	--$i;
    } elsif ($ARGV[$i] eq '--no-gu-close') {
	++$NO_GU_CLOSE;
	splice @ARGV, $i, 1;
	--$i;
    } elsif ($ARGV[$i] eq '--no-double-gu-close') {
	++$NO_DOUBLE_GU_CLOSE;
	splice @ARGV, $i, 1;
	--$i;
    } elsif ($ARGV[$i] eq '--no-gu-basepairs') {
	++$NO_GU_BASEPAIRS;
	splice @ARGV, $i, 1;
	--$i;
    } elsif ($ARGV[$i] eq '--null-energies') {
	++$NULL_ENERGIES;
	splice @ARGV, $i, 1;
	--$i;
    }
}

my ($INFINITY, $inFile, $outFile) = @ARGV;
my ($prefix) = $inFile =~ /^(.+)\./;

open IN, '<', $inFile or die $!;
open OUT, '>', $outFile or die $!;

if ($prefix =~ /triloop/ or $prefix =~ /tloop/ or $prefix =~ /hexaloop/) {
    my $junk = <IN>;
    $junk = <IN>;
    my %loops;
    while (<IN>) {
	nullify if $NULL_ENERGIES;
	my ($loop, $energy) = split or last;
	$loops{$loop} = $energy;
    }
    foreach my $loop (sort keys %loops) {
	print OUT "$loop\t$loops{$loop}\n";
    }
} elsif ($prefix =~ /tstack/) {
    my $i = 0;
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/ or not /\d/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES;
	s/([^\d])\.([^\d])/$1$INFINITY$2/g;

	my @data = split;
	if ($NO_WC_MISMATCHES) {
	    @data[3,7,11,15] = ($INFINITY) x 4 if $i % 4 == 0;
	    @data[2,6,10,14] = ($INFINITY) x 4 if $i % 4 == 1;
	    @data[1,5, 9,13] = ($INFINITY) x 4 if $i % 4 == 2;
	    @data[0,4, 8,12] = ($INFINITY) x 4 if $i % 4 == 3;
	}
	if ($NO_GU_MISMATCHES) {
	    @data[3,7,11,15] = ($INFINITY) x 4 if $i % 4 == 2;
	    @data[2,6,10,14] = ($INFINITY) x 4 if $i % 4 == 3;
	}
	if ($NO_GU_CLOSE or $NO_GU_BASEPAIRS) {
	    if (int($i / 4) == 2) {
		@data[12..15] = ($INFINITY) x 4;
	    } elsif (int($i / 4) == 3) {
		@data[8..11] = ($INFINITY) x 4;
	    }
	}

	print OUT join("\t", @data), "\n";

	++$i;
    }
} elsif ($prefix =~ /stack/ or $prefix =~ /coaxial/) {
    my $i = 0;
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/ or not /\d/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES;
	s/([^\d])\.([^\d])/$1$INFINITY$2/g;

	my @data = split;
	if ($NO_GU_GU_STACKS) {
	    $data[15] = $INFINITY if $i == 10;
	    $data[14] = $INFINITY if $i == 11;
	    $data[11] = $INFINITY if $i == 14;
	    $data[10] = $INFINITY if $i == 15;
	}
	if ($NO_GU_BASEPAIRS) {
	    if ($i % 4 == 2) {
		@data[3,7,11,15] = ($INFINITY) x 4;
	    } elsif ($i % 4 == 3) {
		@data[2,6,10,14] = ($INFINITY) x 4;
	    }
	    if (int($i / 4) == 2) {
		@data[12..15] = ($INFINITY) x 4;
	    } elsif (int($i / 4) == 3) {
		@data[8..11] = ($INFINITY) x 4;
	    }
	}

	print OUT join("\t", @data), "\n";

	++$i;
    }
} elsif ($prefix =~ /dangle/) {
    my $i = 0;
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/ or not /\d/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES;
	s/([^\d])\.([^\d])/$1$INFINITY$2/g;

	my @data = split;
	if ($NO_GU_BASEPAIRS) {
	    if ($i % 4 == 2) {
		@data[12..15] = ($INFINITY) x 4;
	    } elsif ($i % 4 == 3) {
		@data[8..11] = ($INFINITY) x 4;
	    }
	}

	print OUT join("\t", @data), "\n";

	++$i;
    }
} elsif ($prefix =~ /sint2/) {
    my $i = 0;
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/ or not /\d/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES;

	my @data = split;
	if ($i % 4 == 0) {
	    @data[3,7,11,15,19,23] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	} elsif ($i % 4 == 1) {
	    @data[2,6,10,14,18,22] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	} elsif ($i % 4 == 2) {
	    @data[1,5,9,13,17,21] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	    @data[3,7,11,15,19,23] = ($INFINITY) x 6 if $NO_GU_MISMATCHES;
	} else {
	    @data[0,4,8,12,16,20] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	    @data[2,6,10,14,18,22] = ($INFINITY) x 6 if $NO_GU_MISMATCHES;
	}
	@data[16..23] = ($INFINITY) x 8 if ($i >= 16 and $NO_DOUBLE_GU_CLOSE) or $NO_GU_CLOSE or $NO_GU_BASEPAIRS;
	@data = ($INFINITY) x 24 if $i >= 16 and ($NO_GU_CLOSE or $NO_GU_BASEPAIRS);

	print OUT join("\t", @data), "\n";

	++$i;
    }
} elsif ($prefix =~ /asint1x2/) {
    my $i = 0;
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/ or not /\d/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES;

	my @data = split;
	if ($i % 4 == 0) {
	    @data[3,7,11,15,19,23] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	} elsif ($i % 4 == 1) {
	    @data[2,6,10,14,18,22] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	} elsif ($i % 4 == 2) {
	    @data[1,5,9,13,17,21] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	    @data[3,7,11,15,19,23] = ($INFINITY) x 6 if $NO_GU_MISMATCHES;
	} else {
	    @data[0,4,8,12,16,20] = ($INFINITY) x 6 if $NO_WC_MISMATCHES;
	    @data[2,6,10,14,18,22] = ($INFINITY) x 6 if $NO_GU_MISMATCHES;
	}
	if ($i % 4 + int($i / 4) % 4 == 3) {
	    @data = ($INFINITY) x 24 if $NO_WC_MISMATCHES;
	} elsif ($i % 4 + int($i / 4) % 4 == 5) {
	    @data = ($INFINITY) x 24 if $NO_GU_MISMATCHES;
	}
	@data[16..23] = ($INFINITY) x 8 if $NO_GU_CLOSE or $NO_GU_BASEPAIRS;
	@data = ($INFINITY) x 24 if $i >= 64 and ($NO_GU_CLOSE or $NO_GU_BASEPAIRS);
	@data[16..23] = ($INFINITY) x 8 if $i >= 64 and $NO_DOUBLE_GU_CLOSE;

	print OUT join("\t", @data), "\n";

	++$i;
    }
} elsif ($prefix =~ /sint4/) {
    my $i = 0;
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/ or not /\d/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES;

	my @data = split;
	if ($i % 16 == 3 or $i % 16 == 6 or $i % 16 == 9 or $i % 16 == 12) {
	    @data = ($INFINITY) x 16 if $NO_WC_MISMATCHES;
	} elsif($i % 16 == 11 or $i % 16 == 14) {
	    @data = ($INFINITY) x 16 if $NO_GU_MISMATCHES;
	}
	@data[3,6,9,12] = ($INFINITY) x 4 if $NO_WC_MISMATCHES;
	@data[11,14] = ($INFINITY) x 2 if $NO_GU_MISMATCHES;

	@data = ($INFINITY) x 16 if (int($i / 16) == 4 or int($i / 16) == 5 or int($i / 16) == 10 or int($i / 16) == 11 or int($i / 16) == 16 or int($i / 16) == 17 or int($i / 16) >= 22) and ($NO_GU_CLOSE or $NO_GU_BASEPAIRS);
	@data = ($INFINITY) x 16 if (int($i / 16) == 28 or int($i / 16) == 29 or int($i / 16) == 34 or int($i / 16) == 35) and $NO_DOUBLE_GU_CLOSE;

	print OUT join("\t", @data), "\n";

	++$i;
    }
} elsif ($prefix =~ /miscloop/) {
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES and not /1.07/;
	s/([^\d])\.([^\d])/$1$INFINITY$2/g;
	s/\s+/\t/g;
	s/^\s//;
	s/\s$/\n/;
	print OUT;
    }    
} elsif ($prefix =~ /extinction/) {
    while (<IN>) {
	s/[A-z]//g;
	print OUT;
    }
} else {
    while (<IN>) {
	next if /--/ or /[A-z]/ or /=/;
	s/\+/ +/g;
	s/-/ -/g;
	nullify if $NULL_ENERGIES;
	s/([^\d])\.([^\d])/$1$INFINITY$2/g;
	s/\s+/\t/g;
	s/^\s//;
	s/\s$/\n/;
	print OUT;
    }
}

close IN or die $!;
close OUT or die $!;
