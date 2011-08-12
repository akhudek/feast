#!/usr/bin/perl
use strict;

my $cvs = $ARGV[0];
my $bn_file = $ARGV[1];
my $cxx_file = $ARGV[2];
my $cmd = $ARGV[3];

open(BN, $bn_file) or die("Failed to open $bn_file");
my $build_number = <BN>;
close(BN);
chomp($build_number);

$build_number =~ /^(\d+)/;
my $build_tag = "BUILD$1";

sub write_build_number {
	my( $bn, $bnf ) = @_;
	open( BN, ">$bnf" ) or die( "Failed ot open $bnf for writing." );
	print BN $bn;
	close( BN );
}

if( $cmd =~ /^generate$/i ){
	open( CX, ">$cxx_file" ) or die ("Failed to open $cxx_file for writing.");
	print CX "const char* build_string = \"$build_number\";\n\n";
	close(CX);
	exit 0;
	
} elsif( $cmd =~ /^verify$/i ) {
	if( $build_number =~ /^\d+$/ && -x $cvs ) {
		# check if trunk has changed
		my $cvs_cmd = "$cvs -Q diff -r $build_tag --brief";
		my @cvs_output = grep( /^Index/, `$cvs_cmd 2>&1` );
      @cvs_output = grep( !/Index: CMakeLists.txt/, @cvs_output );
		if( @cvs_output > 0 ) { $build_number .= "+"; }
		&write_build_number( $build_number, $bn_file );
	} elsif( $build_number =~ /^\d+$/ ) {
		warn "Missing CVS, cannot verify build files.";
	}
	
} elsif( $cmd =~ /^freeze$/i ) {
	die( "Missing CVS." ) if ! -x $cvs;

	# check for uncommited changes
	my $cvs_cmd = "$cvs -Q diff --brief";
	my @cvs_output = `$cvs_cmd 2>&1`;
	die( "Cannot freeze build, uncommited changes in source tree." ) if @cvs_output > 0;

	# otherwise tag and write new build number
	$build_number =~ /^(\d+)/;
	$build_number = $1 + 1;
	$build_tag = "BUILD$build_number";
	&write_build_number( $build_number, $bn_file );
	$cvs_cmd = "$cvs commit -m \"Increment build number.\";";
	die( "CVS commit failed") if system($cvs_cmd) != 0;
	$cvs_cmd = "$cvs tag $build_tag";
	die( "CVS tag failed") if system($cvs_cmd) != 0;
	&write_build_number( $build_number, $bn_file );
	exit 0;
}
