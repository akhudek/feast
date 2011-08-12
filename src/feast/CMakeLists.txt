set(feast_SRCS
	feast_main.cc
)

set(feast_LIBS
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
add_executable(feast ${feast_SRCS})

target_link_libraries(feast ${feast_LIBS})
