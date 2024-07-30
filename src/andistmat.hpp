#ifndef ANDISTMAT_HPP
#define ANDISTMAT_HPP

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
		void get_dists(const std::vector<uint32_t>&, std::vector<double>& dists) const;
		uint32_t get_medoid(std::vector<uint32_t>&) const;
		uint32_t size()const;
};
#endif