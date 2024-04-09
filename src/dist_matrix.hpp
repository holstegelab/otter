#ifndef DIST_MATRIX_HPP
#define DIST_MATRIX_HPP

#include <cstdint>
#include <vector>


class DistMatrix{
	public:
		uint32_t n;
		std::vector<double> values;
		DistMatrix(uint32_t);
		DistMatrix(uint32_t&, std::vector<double>&);
		void set_dist(uint32_t, uint32_t, double);
		double get_dist(uint32_t, uint32_t) const;
		uint32_t get_min_dist_i(std::vector<uint32_t>&) const;
		uint32_t size()const;
		void print() const;
};
#endif
