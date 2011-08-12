set(extension_store_test_LIBS
	phmm
	math
)
	
add_executable(extension_store_test main.cc)

target_link_libraries(extension_store_test ${extension_store_test_LIBS})