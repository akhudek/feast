# Alignment sub-library.
set(alignment_SRCS 
	simple_nw_model_parameters.cc
	basic_nw_model_parameters.cc 
	nw_model_parameters.cc 
	alignment_functions.cc
	annotation_interval.cc
	annotation_io.cc
	annotation.cc
	mm_alignment.cc
	needleman.cc
	multi_alignment_model.cc
	sm_extension.cc
	parameter_functions.cc
	mm_trainer.cc
)

add_library(alignment ${alignment_SRCS})

