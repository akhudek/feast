FEAST local aligner
===================

Requirements
------------

* A recent version of g++.
* Cmake 2.6+
* Boost 1.34+ 
* Blitz++

Running the build.pl script will produce a static binary. The build is
tested on Linux 64/32 and Mac 10.6. 

Basic Usage
-----------

To align two sequences you must supply a parameter file such as 2-R8. 
Input is expected to be in fasta format. Repeats should be lowercase
masked. There can be any number of query sequences in the query file,
but only the first sequence from the target file will be used. The
output format is MAF. If no output file is specified, standard output
is used.

`$ feast -t target.fa -q query.mfa -o output.maf -p 2-R8`

Training
--------

To train parameters, you may either use an existing parameter file, or
you can supply an init file that initializes default parameters. See
the supplied sample init file for more information on the format. 
The default training method is Baum-Welch. However, this method does
not work for more than one sub-model at this time. Thus we recommend
switching to Viterbi.  In training mode the output is a parameter set.
Status information is sent to standard error. If you break the 
procedure before it completes, the parameters for the last training
iteration will be in in the output file. We do not support the
standard output stream in this case, you must specify a file.

`$ feast -t target.fa -q query.fa -o trained_parameters \
 -i init_file --train --cut-threshold 50`

IMPORTANT: feast assumes lowercase masked input files. Thus,
input that is completely lowercase will return no results. Use
the --ignore-mask flag to disable masking of repeats and treat
all files as uppercase.

Please contact (akhudek@cs.uwaterloo.ca) if you have further questions.


