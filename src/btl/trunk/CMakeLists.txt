cmake_minimum_required(VERSION 2.6)
project(libbtl)
set(btl_SRCS
	dna_gap.cc			dna_gap_symbol.cc
	dna_iupac.cc		dna_iupac_symbol.cc
	dna.cc				dna_symbol.cc
	mdna.cc				mdna_symbol.cc
	mdna_iupac.cc		mdna_iupac_symbol.cc
	mdna_iupac_gap.cc	mdna_iupac_gap_symbol.cc)

add_library(btl ${btl_SRCS})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS} ./)
