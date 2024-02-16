#include "kde.hpp"
#include <utility>
#include <iostream>
#include <cmath>

KDE::KDE(double _h): h(_h){};

double KDE::k(double x) const
{
	return (1/std::sqrt(2*pi))*(std::exp(-(x*x/2)));
}

double KDE::k_h(double x) const
{
	return (1/h)*(k(x/h));
}

double KDE::f(double x) const
{
	double total = 0.0;
	for(const auto& v : values) total += k_h(x - v);
	return total / values.size();
}

void KDE::maximas(const std::vector<double>& densities, std::vector<std::pair<int,double>>& maxs, std::vector<std::pair<int,double>>& mins) const
{
	//int width = 1;
	bool find_maxima = true;
	double last_sum = 0.0;
	int last_sum_i = 1;
	for(int i = 1; i < (int)densities.size() - 1; ++i){
		double sum = 0.0;
		sum += densities[i];
		sum += densities[i - 1];
		sum += densities[i + 1];
		//for(int j = 0; j < width; ++j) sum += densities[i - 1];//(width * j)];
		//for(int j = 0; j < width; ++j) sum += densities[i + 1];//(width * j)];
		if(find_maxima){
			if(sum < last_sum) {
				find_maxima = false;
				maxs.emplace_back(std::make_pair(last_sum_i, last_sum));
			}
		}
		else {
			if(sum > last_sum) {
				find_maxima = true;
				mins.emplace_back(std::make_pair(last_sum_i, last_sum));
			}
		}
		last_sum = sum;
		last_sum_i = i;
	}
	if(find_maxima) maxs.emplace_back(std::make_pair(last_sum_i, last_sum));
	//for(const auto& m : maxs) std::cout << m.first << '\t' << m.second << '\n';
}