#include "opinterval.hpp"
#include "interval_tree.h"
#include "sam.h"
#include <string>
#include <vector>
#include <utility>

OpInterval::OpInterval(){};

OpInterval::OpInterval(const int& s, const int& e, const int& o): start(s), end(e), op(o){};

void get_op_intervals(bam1_t*& alignment, std::vector<std::pair<int,int>>& ref_intervals, std::vector<OpInterval>& op_intervals)
{
	uint32_t *cigar = bam_get_cigar(alignment);
	//current positions in ref and query
	int rpos = alignment->core.pos, rpos_acc = rpos;
	int qpos = 0, qpos_acc = qpos;
	for (uint32_t i = 0; i < alignment->core.n_cigar; ++i) { 
        const int op = bam_cigar_op(cigar[i]);
		const int ol = bam_cigar_oplen(cigar[i]);
		if(op == BAM_CSOFT_CLIP) qpos_acc += ol;
		else if(op == BAM_CMATCH || op ==  BAM_CEQUAL || op == BAM_CDIFF){
			rpos_acc += ol;
			qpos_acc += ol;
		}
		else if(op == BAM_CINS) qpos_acc += ol;
		else if(op == BAM_CDEL) rpos_acc += ol;
		ref_intervals.emplace_back(rpos, rpos_acc);
		op_intervals.emplace_back(qpos, qpos_acc, op);
		rpos = rpos_acc;
		qpos = qpos_acc;

	}
}

bool is_contained(const int& x_1, const int& x_2, const int& y){return y >= x_1 && y <= x_2;}

void intersect_bed_op_interval(const int& rstart, const int& rend, const std::vector<Interval<int, OpInterval>>& intervals, int& qstart, int& qend)
{
	bool is_qstart_set = false;
	bool is_qend_set = false;
	for(const auto& interval : intervals){
		if(is_qstart_set && is_qend_set) break;
		else {
			if(!is_qstart_set || !is_qend_set){
				bool is_start_contained = is_contained(interval.start, interval.stop, rstart);
				bool is_end_contained = is_contained(interval.start, interval.stop, rend);
				if(interval.value.op == BAM_CMATCH || interval.value.op == BAM_CEQUAL || interval.value.op == BAM_CDIFF) {
					if(is_start_contained && !is_qstart_set) {
						qstart = interval.value.start + (rstart - interval.start);
						is_qstart_set = true;
					}
					if(is_end_contained && !is_qend_set){
						qend = interval.value.start + (rend - interval.start);
						is_qend_set = true;
					}
				}
				else if(interval.value.op == BAM_CDEL){
					if(is_start_contained && !is_qstart_set){
						qstart = interval.value.end;
						is_qstart_set = true;
					}
					if(is_end_contained && !is_qend_set){
						qend = interval.value.start;
						is_qend_set = true;
					}
				}
			}
		}
	}
}