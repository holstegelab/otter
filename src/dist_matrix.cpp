#include <vector>
#include <exception>
#include <iostream>
#include <cstddef>
#include <cmath>
#include "dist_matrix.hpp"

DistMatrix::DistMatrix(uint32_t _n): n(_n), d_n((n * (n - 1)) / 2)
{
	values.resize(d_n, 1.0);
}

void DistMatrix::set_dist(uint32_t i, uint32_t j, double d)
{
	if(i == j) throw std::exception();
	int a = i < j ? i : j;
	int b = i > j ? i : j;
	values[(static_cast<std::ptrdiff_t>(2*n-3-(a))*(a)>>1)+(b)-1] = d;
}

double DistMatrix::get_dist(uint32_t i, uint32_t j) const 
{
	if(i == j) throw std::exception();
	int a = i < j ? i : j;
	int b = i > j ? i : j;
	return values[(static_cast<std::ptrdiff_t>(2*n-3-(a))*(a)>>1)+(b)-1];
}

uint32_t DistMatrix::get_min_dist_i(std::vector<uint32_t>& indeces) const
{
	uint32_t min_i = indeces.front();
	//set infinity
	double min_dist_sum = 100000000.0;
	for(const auto& i : indeces){
		double dist_sum = 0.0;
		for(const auto& j : indeces) if(i != j) dist_sum += get_dist(i,j);
		if(dist_sum < min_dist_sum) {
			min_i = i;
			min_dist_sum = dist_sum;
		}
	}
	return min_i;
}

uint32_t DistMatrix::size()const{return values.size();}