#include "otter_opts.hpp"
#include "antimestamp.hpp"
#include <sstream>
#include <vector>
#include <iostream>
#include <ostream>
#include <cstdint>

void parse_line(const char delim, std::string& line, std::vector<std::string>& values)
{
	auto it = line.begin();
	while(it != line.end()){
		if(std::isspace(*it)) it = line.erase(it);
		else ++it;
	}
	std::string value;
	std::istringstream stringstream(line);
	while(std::getline(stringstream, value, delim)) values.emplace_back(value);
}

bool is_zero_one_range(const double& input)
{
	return input >= 0.0 && input <= 1.0;
}

void OtterOpts::init_offset(std::string tmp)
{
	std::vector<std::string> tmp_inputs;
	parse_line(',', tmp, tmp_inputs);

	if(tmp_inputs.size() == 1){
		offset_l = std::stoi(tmp_inputs[0]);
		offset_r = std::stod(tmp_inputs[0]);
	}
	else if(tmp_inputs.size() == 2) {
		offset_l = std::stoi(tmp_inputs[0]);
		offset_r = std::stod(tmp_inputs[1]);
	}
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid offset value: " << tmp <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_max_alleles(int _max_alleles)
{
	if(_max_alleles >= 0) max_alleles = _max_alleles;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid maximum-alleles value: " << _max_alleles <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_mapq(int _mapq)
{
	if(_mapq >= 0 && _mapq <= 60) mapq = _mapq;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid mapq value: " << _mapq <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_read_quality(double _rq)
{
	if(is_zero_one_range(_rq)) read_quality = _rq;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid read-quality value: " << _rq <<  std::endl;
		exit(0);
	}
}


void OtterOpts::init_max_cov(int _max_cov)
{
	if(_max_cov >= 0) max_cov = _max_cov;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid max-coverage value: " << _max_cov <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_min_cov_fraction(double _min_cov_fraction)
{
	if(is_zero_one_range(_min_cov_fraction)) min_cov_fraction = _min_cov_fraction;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid _min_cov_fraction value: " << _min_cov_fraction <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_threads(int _threads)
{
	if(_threads > 0 && _threads <= 32) threads = _threads;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid threads value: " << _threads <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_max_error(double _max_error)
{
	if(is_zero_one_range(_max_error)) max_error = _max_error;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid max-error value: " << _max_error <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_max_cosdis(double _max_cosdis)
{
	if(is_zero_one_range(_max_cosdis)) max_cosdis = _max_cosdis;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid max cosine-dissimilarity value: " << _max_cosdis <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_bandwidth(double _bandwidth)
{
	if(is_zero_one_range(_bandwidth)) bandwidth = _bandwidth;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid bandwidth value: " << _bandwidth <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_flank(int _flank)
{
	if(_flank >= 21 && _flank < 10000) flank = _flank;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid flanking-sequence size for realignment: " << _flank <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_min_sim(double _min_sim)
{
	if(is_zero_one_range(_min_sim)) min_sim = _min_sim;
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] Invalid min-similarity for realignment: " << _min_sim <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_min_cov_fraction2(std::string tmp)
{
	std::vector<std::string> tmp_inputs;
	parse_line(',', tmp, tmp_inputs);
	if(tmp_inputs.size() == 2) {
		min_cov_fraction2_l = std::stoi(tmp_inputs[0]);
		min_cov_fraction2_f = std::stod(tmp_inputs[1]);
	}
	else{
		std::cerr << '(' << antimestamp() << "): [ERROR] expected two comma-separated values: " << tmp <<  std::endl;
		exit(0);
	}
}

void OtterOpts::init_dist_length(std::string tmp)
{
	std::vector<std::string> tmp_inputs;
	std::string value;
	std::istringstream stringstream(tmp);
	while(std::getline(stringstream, value, ',')) tmp_inputs.emplace_back(value);
	if(tmp_inputs.size() != 3) {
		std::cerr << '(' << antimestamp() << "): [ERROR] expected three comma-separated values: " << tmp <<  std::endl;
		exit(0);
	}
	else {
		dist_length_large = std::stoi(tmp_inputs[0]);
		dist_length_frac = std::stod(tmp_inputs[1]);
		dist_length_se = std::stod(tmp_inputs[2]);
	}
}