#ifndef OPINTERVAL_HPP
#define OPINTERVAL_HPP

#include "sam.h"
#include "interval_tree.h"
#include <string>
#include <vector>
#include <utility>

class OpInterval{
	public:
		int start;
		int end;
		int op;
		OpInterval();
		OpInterval(const int& s, const int& e, const int& o);
};


void get_op_intervals(
	bam1_t*&, 
	std::vector<std::pair<int,int>>&,
	std::vector<OpInterval>&
);

bool is_contained(
	const int&, 
	const int&, 
	const int&
);

void intersect_bed_op_interval(
	const int&, 
	const int&, 
	const std::vector<Interval<int, OpInterval>>&, 
	int&, 
	int&
);


#endif