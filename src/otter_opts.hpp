#ifndef OTTEROPTS_HPP
#define OTTEROPTS_HPP

#include <cstdint>
#include <string>

class OtterOpts {
	public:
		int offset;
		int max_alleles;
		int mapq;
		int max_cov;
		double min_cov_fraction;
		int min_cov_fraction2_l;
		double min_cov_fraction2_f;
		int threads;
		double max_error;
		double bandwidth;
		int flank;
		double min_sim;
		bool nonprimary;
		bool is_sam;
		bool omitnonspanning;
		bool is_wga;
		std::string read_group;

		void init_offset(int);
		void init_max_alleles(int);
		void init_mapq(int);
		void init_max_cov(int);
		void init_min_cov_fraction(double);
		void init_threads(int);
		void init_max_error(double);
		void init_bandwidth(double);
		void init_flank(int);
		void init_min_sim(double);
		void init_min_cov_fraction2(std::string);

};

#endif