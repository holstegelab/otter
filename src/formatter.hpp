#include <string>
#include "otter_opts.hpp"

void output_fa2sam(
	const std::string& read,
	const std::string& chr, 
	const int& start, 
	const int& end, 
	const std::string& seq, 
	const std::string& rg, 
	const int& acov, 
	const int& tc, 
	const int& spanning_l, 
	const int& spanning_r,
	const double& se
	);
