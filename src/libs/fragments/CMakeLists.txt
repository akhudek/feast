set(fragment_SRCS
	fragment_algorithm.cc
	fragment.cc
	scored_fragment.cc
	transformed_fragment.cc
)

add_library(fragments ${fragment_SRCS})