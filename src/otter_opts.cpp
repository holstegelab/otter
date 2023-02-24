#include "otter_opts.hpp"
#include "antimestamp.hpp"
#include <iostream>
#include <ostream>
#include <cstdint>

void OtterOpts::init_offset(int _offset)
{
	if(_offset >= 0) offset = _offset;
	else{
		std::cout << '(' << antimestamp() << "): Invalid offset value: " << _offset <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_max_alleles(int _max_alleles)
{
	if(_max_alleles > 0) max_alleles = _max_alleles;
	else{
		std::cout << '(' << antimestamp() << "): Invalid maximum-alleles value: " << _max_alleles <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_mapq(int _mapq)
{
	if(_mapq >= 0 && _mapq <= 60) mapq = _mapq;
	else{
		std::cout << '(' << antimestamp() << "): Invalid mapq value: " << _mapq <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_max_cov(int _max_cov)
{
	if(_max_cov >= 0) max_cov = _max_cov;
	else{
		std::cout << '(' << antimestamp() << "): Invalid max-coverage value: " << _max_cov <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_threads(int _threads)
{
	if(_threads > 0 && _threads <= 32) threads = _threads;
	else{
		std::cout << '(' << antimestamp() << "): Invalid threads value: " << _threads <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_max_error(double _max_error)
{
	if(_max_error >= 0.0 && _max_error <= 1.0) max_error = _max_error;
	else{
		std::cout << '(' << antimestamp() << "): Invalid max-error value: " << _max_error <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_bandwidth(double _bandwidth)
{
	//process mapq accordingly
	if(_bandwidth >= 0.0 && _bandwidth <= 1.0) bandwidth = _bandwidth;
	else{
		std::cout << '(' << antimestamp() << "): Invalid bandwidth value: " << _bandwidth <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_flank(int _flank)
{
	if(_flank >= 21 && _flank < 10000) flank = _flank;
	else{
		std::cout << '(' << antimestamp() << "): Invalid flanking-sequence size for realignment: " << _flank <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_min_sim(double _min_sim)
{
	if(_min_sim >= 0.0 && _min_sim < 1.0) min_sim = _min_sim;
	else{
		std::cout << '(' << antimestamp() << "): Invalid min-similarity for realignment: " << _min_sim <<  std::endl;
		exit(0);
	}
}