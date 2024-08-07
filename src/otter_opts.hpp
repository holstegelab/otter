#ifndef OTTEROPTS_HPP
#define OTTEROPTS_HPP

#include <cstdint>
#include <string>

class OtterOpts {
	public:
		int offset_l;
		int offset_r;
		int max_alleles;
		int mapq;
		double read_quality;
		int max_cov;
		double min_cov_fraction;
		int min_cov_fraction2_l;
		double min_cov_fraction2_f;
		int threads;
		double max_error;
		double bandwidth_short;
		double bandwidth_long;
		int bandwidth_length;
		int flank;
		double min_sim;
		bool nonprimary;
		bool is_fa;
		bool omitnonspanning;
		bool is_wga;
		bool ignore_haps;
		bool is_debug;
		std::string read_group;
		int dist_length_large;
		double dist_length_frac;
		double dist_length_se;
		double max_cosdis;

		void init_offset(std::string);
		void init_max_alleles(int);
		void init_mapq(int);
		void init_read_quality(double);
		void init_max_cov(int);
		void init_min_cov_fraction(double);
		void init_threads(int);
		void init_max_error(double);
		void init_max_cosdis(double);
		void init_bandwidth(std::string);
		void init_flank(int);
		void init_min_sim(double);
		void init_min_cov_fraction2(std::string);
		void init_dist_length(std::string);

};

#endif