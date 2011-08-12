find_path(Blitz_INCLUDE_DIRS blitz/blitz.h 
	/usr/include 
	/usr/local/include
	/sw/include
)

find_library(Blitz_LIBRARIES blitz
	/usr/lib
	/usr/local/lib
	/sw/lib
)

set(Blitz_FOUND TRUE)

if(NOT Blitz_INCLUDE_DIRS)
	set(Blitz_FOUND FALSE)
endif(NOT Blitz_INCLUDE_DIRS)

if(NOT Blitz_LIBRARIES)
	set(Blitz_FOUND FALSE)
endif(NOT Blitz_LIBRARIES)
