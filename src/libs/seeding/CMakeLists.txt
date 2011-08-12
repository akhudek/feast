set(seeding_SRCS
	seed.cc
	hash_hit_data.cc
	random_seed_factory.cc
)

add_library(seeding ${seeding_SRCS})