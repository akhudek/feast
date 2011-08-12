set(cape_SRCS
	cape_main.cc
	cape_pairwise_aligner.cc
	custom_fragment_score_worker.cc
	dynamic_fragment_score_worker.cc
)

set(cape_LIBS
	${Boost_LIBRARIES}
   ${COMMON_LIBRARIES}
	fragments
	io
	alignment
	random
	math
	seeding
	sequence
	buildnumber
	btl
)	
add_executable(cape ${cape_SRCS})

target_link_libraries(cape ${cape_LIBS})
