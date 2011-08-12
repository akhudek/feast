set(explore_SRCS
    explore_main.cc
    explore_io.cc
    dot_plot.cc
)

set(explore_LIBS
	btl
	${Boost_LIBRARIES}
	math
	buildnumber
)	
add_executable(explore ${explore_SRCS})

target_link_libraries(explore ${explore_LIBS})