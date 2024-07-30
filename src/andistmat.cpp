#include "andistmat.hpp"
#include <vector>
#include <exception>
#include <iostream>
#include <cstddef>
#include <cmath>

DistMatrix::DistMatrix(uint32_t _n): n(_n)
{
	values.resize((n * (n - 1)) / 2, 1.0);
}

DistMatrix::DistMatrix(uint32_t& _n, std::vector<double>& _values):n(_n),values(_values){}

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

void DistMatrix::get_dists(const std::vector<uint32_t>& indeces, std::vector<double>& dists) const
{
	for(uint32_t i = 0; i < indeces.size(); ++i) for(uint32_t j = i+1; j < indeces.size(); ++j) dists.emplace_back(get_dist(i, j));
}

uint32_t DistMatrix::get_medoid(std::vector<uint32_t>& indeces) const
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