set(rda_SRCS
	rda_main.cc
  krdalign.cc
  rda_functions.cc
)

set(rda_LIBS
   imp_boost_thread
   imp_boost_program_options
   ${COMMON_LIBRARIES}
	alignment
	math
	fragments
	io
	random
	sequence
	seeding
	chaining
	phmm
	buildnumber
	btl
)	
add_executable(rda ${rda_SRCS})

target_link_libraries(rda ${rda_LIBS})
