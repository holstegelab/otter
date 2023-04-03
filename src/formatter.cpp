#include "formatter.hpp"
#include <string>
#include <iostream>
#include "otter_opts.hpp"

void output_fa2sam(const std::string& read, const std::string& chr, const int& start, const int& end, const std::string& seq, const std::string& rg, const int& acov, const int& tc, const int& spanning_l, const int& spanning_r)
{
	std::string pseudo_qual(seq.size(), '!');
	if(!rg.empty()) pseudo_qual += "\tRG:Z:" + rg;
	std::cout << read << "\t0\t" << chr << '\t' << start << "\t0\t" << seq.size() << "M\t*\t0\t0\t" << seq << '\t' << pseudo_qual << '\t' << "ta:Z:" << chr << ':' << start << '-' << end;
	if(acov > 0) std::cout << "\tac:i:" << acov;
	std::cout << "\ttc:i:" << tc;
	if(spanning_l >= 0) {
		char sp;
		if(spanning_l && spanning_r) sp = 'b';
		else if(spanning_l) sp = 'l';
		else if(spanning_r) sp = 'r';
		else sp = 'n';
		std::cout << "\tsp:A:" << sp;
	}
	std::cout << '\n';
}