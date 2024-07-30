#ifndef ANKDE_HPP
#define ANKDE_HPP

#include <vector>
#include <utility>

class KDE {
	public:
		double h;
		std::vector<double> values;

		KDE(double);
		double f(double) const;
		void maximas(const int&, const std::vector<double>&, std::vector<std::pair<int,double>>&, std::vector<std::pair<int,double>>&) const;
	
	private:
		double pi = 3.14159265358979323846;
		double k(double) const;
		double k_h(double) const;
};

#endif