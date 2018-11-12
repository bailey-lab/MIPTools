#! /usr/bin/perl -w

use strict;
use warnings;
use IPC::Open3;
use XML::Parser;

sub printXML ($) {
    my $type = ref $_[0];
    $type =~ s/^main:://;

    if ($type eq 'Characters') {
	print $_[0]->{Text};
	return;
    }

    print "<$type";
    foreach my $attribute (sort keys %{$_[0]}) {
	next if $attribute eq 'Kids';
	print " $attribute=\"$_[0]->{$attribute}\"";
    }
    print '>';
    foreach my $kid (@{$_[0]->{Kids}}) {
	::printXML $kid if ref($kid) =~ /^main::/;
    }
    print "</$type>";
}

unless (@ARGV) {
    print "Usage: $0 <testfile>...\n";
    exit 0;
}

my $xp = new XML::Parser(Style => 'Objects');
foreach (@ARGV) {
    my $tree = $xp->parsefile($_); 
    ref $tree->[0] eq 'main::tests' or die 'Bad root element';
    foreach my $test (@{$tree->[0]->{Kids}}) {
	next if ref $test eq 'main::Characters';
	ref $test eq 'main::test' or die 'Bad test element';

	my $command = $test->{command} or die 'No command';
	
	my $stdin = join '', map {
	    $_->{Kids}->[0]->{Text}
	} grep {
	    ref $_ eq 'main::stdin' and ref $_->{Kids}->[0] eq 'main::Characters'
	} @{$test->{Kids}};

	local (*WTR, *RDR, *ERR);
	my $pid = open3(\*WTR, \*RDR, \*ERR, $command);
	print WTR $stdin;
	close WTR or die $!;
	my $stdout = join '', <RDR>;
	close RDR or die $!;
	my $stderr = join '', <ERR>;
	close ERR or die $!;
	waitpid $pid, 0;

	$test->{return} = $? >> 8;

	my ($stdoutNode) = grep {
	    ref $_ eq 'main::stdout'
	} @{$test->{Kids}};
	$stdoutNode->{Kids}->[0]->{Text} = $stdout;
	bless $stdoutNode->{Kids}->[0], 'main::Characters';

	my ($stderrNode) = grep {
	    ref $_ eq 'main::stderr'
	} @{$test->{Kids}};
	$stderrNode->{Kids}->[0]->{Text} = $stderr;
	bless $stderrNode->{Kids}->[0], 'main::Characters';

	for ($stdoutNode->{Kids}->[0]->{Text}, $stderrNode->{Kids}->[0]->{Text}) {
	    s'&'&amp;'g;    #'
	    s'<'&lt;'g;     #'
	    s'>'&gt;'g;     #'
	}
    }

    print qq/<?xml version="1.0" encoding="utf-8"?>\n/;
    print qq/<!DOCTYPE tests SYSTEM "test.dtd">\n/;
    printXML $tree->[0];
    print "\n";
}
