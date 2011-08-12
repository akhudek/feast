#ifndef ANNOTATION_FUNCTIONS_H__
#define ANNOTATION_FUNCTIONS_H__
#include <iostream>
#include <string>
#include "annotation.h"
#include "pairwise_dna_alignment.h"

void write_seqa_gff_annotation( std::ostream &out, annotation const &a, pairwise_dna_alignment &align  );

#endif
