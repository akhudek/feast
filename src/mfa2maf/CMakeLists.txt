set(mfa2maf_SRCS mfa2maf_main.cc)

set(mfa2maf_LIBS
   imp_boost_thread
   imp_boost_program_options
   ${COMMON_LIBRARIES}
   buildnumber
	 btl
)	
add_executable(mfa2maf ${mfa2maf_SRCS})

target_link_libraries(mfa2maf ${mfa2maf_LIBS})
