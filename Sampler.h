#pragma once

#include <cstdint>
#include <algorithm>

#include <random>

#define ONE_MINUS_EPSILON float(0x1.fffffep-1)

// Better implementation here https://github.com/wjakob/pcg32
class Pcg {
public:
	uint32_t rng_state;

	Pcg(uint32_t seed) {
		rng_state = seed;
	}


	Pcg() {
		rng_state = 0u;
	}

	uint32_t rand_pcg()
	{
		uint32_t state = rng_state;
		rng_state = rng_state * 747796405u + 2891336453u;
		uint32_t word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
		return (word >> 22u) ^ word;
	}

	uint32_t pcg_hash(uint32_t input)
	{
		uint32_t state = input * 747796405u + 2891336453u;
		uint32_t word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
		return (word >> 22u) ^ word;
	}

	float Uniform() {
		return std::min<float>(ONE_MINUS_EPSILON, rand_pcg() * 0x1p-32f);
	}
};