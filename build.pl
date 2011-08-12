#!/usr/bin/perl
# 
# Build script for feast.
#

my $cmake = `which cmake`;
if( $cmake eq '' ) {
	print "Please make sure that cmake 2.6 or greater is in your PATH.\n";
	exit(1);
}

# Make build directory and move into it.
`mkdir -p build`;
chdir("build") or die("Failed to create build directory.\n");

# Configure the build system.
my $code = system( "cmake ../src" );
if( $code != 0 ) {
	print("Configuration failed. Ensure all required libraries are installed.\n");
	exit($code);
}

# Do the build.
$code = system( "make" );
if( $code != 0 ) {
	print( "Build failed.\n" );
	exit($code);
}

# Move back and copy binary to bin directory.
chdir("..");
`mkdir -p bin`;
`cp build/feast/feast bin`;

print "\nBuild complete. Binary is in bin/feast.\n\n";


