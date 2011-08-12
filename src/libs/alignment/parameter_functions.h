#ifndef MAIN_FUNCTIONS_H__
#define MAIN_FUNCTIONS_H__
#include <string>
#include "dna_sequence.h"
#include "nw_model_parameters.h"

nw_model_parameters_ptr load_parameters_from_file( std::string filename );
nw_model_parameters_ptr estimate_parameters( dna_sequence_ptr_vector &sequences, double d, double a, double L);
#endif
