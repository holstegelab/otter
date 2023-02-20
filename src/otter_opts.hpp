#ifndef OTTEROPTS_HPP
#define OTTEROPTS_HPP

#include <cstdint>

class OtterOpts {
	public:
		int offset;
		int max_alleles;
		int mapq;
		int threads;
		double max_error;
		double bandwidth;
		int flank;
		double min_sim;

		void init_offset(int);
		void init_max_alleles(int);
		void init_mapq(int);
		void init_threads(int);
		void init_max_error(double);
		void init_bandwidth(double);
		void init_flank(int);
		void init_min_sim(double);

};

#endif