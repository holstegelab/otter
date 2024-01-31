#include "parse_bam_alignments.hpp"
#include "anbamfilehelper.hpp"
#include "antimestamp.hpp"
#include "anbed.hpp"
#include <htslib/sam.h>
#include <memory>
#include <string>
#include <vector>
#include <iostream>

ParsingStatus::ParsingStatus(): successful(true),spanning_l(true), spanning_r(true)
{
	alignment_coords = std::make_pair(-1,-1);
};
bool ParsingStatus::is_spanning() const {return spanning_l && spanning_r;};

 bool AlignmentBlock::n_spannable(int n) const
 {
 	int m = 0;
 	for(const auto& s : statuses) if(s.is_spanning()) ++m;
 	return m >= n;
 }

void project_positions(bam1_t* alignment, bool& clipped_l, bool& clipped_r, std::vector<int>& refcoords)
{
	refcoords.resize(alignment->core.l_qseq, -1);
	uint32_t *cigar = bam_get_cigar(alignment);
	//current positions in ref and query
	int rpos = alignment->core.pos;
	int qpos = 0;
	for (uint32_t i = 0; i < alignment->core.n_cigar; ++i) { 
        const int op = bam_cigar_op(cigar[i]);
		const int ol = bam_cigar_oplen(cigar[i]);
		if(op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP){
			if(i == 0) clipped_l = true;
			if(i == alignment->core.n_cigar - 1) clipped_r = true;
			if(op == BAM_CSOFT_CLIP) qpos += ol;
		}
		else if(op == BAM_CMATCH || op ==  7 || op == 8){
			for(int j = 0; j < ol; ++j){
				refcoords[qpos] = rpos;
				++rpos;
				++qpos;
			}
		}
		else if(op == BAM_CINS) qpos += ol;
		else if(op == BAM_CDEL) rpos += ol;
	}
}

void get_breakpoints(const int& start,  const int& end, bam1_t* alignment, ParsingStatus& msg, std::unique_ptr<std::pair<int, int>>& ptr)
{
	//init variables
	bool clipped_l = false, clipped_r = false;
	int qstart_dist = -1, qend_dist = -1;
	int leftmost_q = -1, rightmost_q = -1, leftmost_r = -1, rightmost_r = -1;
	int qstart_q = -1, qend_q = -1, qstart_r = -1, qend_r = -1;
	uint32_t qstart_cigar_i = 0, qend_cigar_i = 0;
	uint32_t *cigar = bam_get_cigar(alignment);
	//current positions in ref and query
	int rpos = alignment->core.pos;
	int qpos = 0;
	//iterate through each cigar operation
	for (uint32_t i = 0; i < alignment->core.n_cigar; ++i) { 
        const int op = bam_cigar_op(cigar[i]);
		const int ol = bam_cigar_oplen(cigar[i]);
		if(op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP){
			if(i == 0) clipped_l = true;
			if(i == alignment->core.n_cigar - 1) clipped_r = true;
			if(op == BAM_CSOFT_CLIP) qpos += ol;
		}
		else if(op == BAM_CMATCH || op ==  7 || op == 8){
			for(int j = 0; j < ol; ++j){
				//update the left-most index
				if(leftmost_q == -1) {
					leftmost_q = qpos;
					leftmost_r = rpos;
				}
				//update the right most index
				if(rightmost_q == -1 || rpos > rightmost_r) {
					rightmost_q = qpos;
					rightmost_r = rpos;
				}
				//distance to start coordinate 
				int cstart_dist = rpos - start;
				//distance to end coordinate
				int cend_dist = end - rpos;
				//new closest start coordinate found
				if(cstart_dist >= 0 && (qstart_dist < 0 || cstart_dist < qstart_dist)){
					qstart_dist = cstart_dist;
					qstart_q = qpos;
					qstart_cigar_i = i;
					qstart_r = rpos;
				}
				//new closest end coordinate found
				if(cend_dist >= 0 && (qend_dist < 0 || cend_dist < qend_dist)){
					qend_dist = cend_dist;
					qend_q = qpos;
					qend_cigar_i = i;
					qend_r = rpos;
				}
				++rpos;
				++qpos;
			}
		}
		else if(op == BAM_CINS) qpos += ol;
		else if(op == BAM_CDEL) rpos += ol;
	}

	//std::cerr << (char*)alignment->data << '\n';
	//std::cerr << start << '-' << end << '\n';
	//std::cerr << "leftmost_r: " << leftmost_r << "; rightmost_r: " << rightmost_r << "; leftmost_q: " << leftmost_q << "; rightmost_q: " << rightmost_q << '\n';
	//std::cerr << "qstart_q: " << qstart_q << ';' << "qend_q: " << qend_q << '\n';
	//std::cerr << "qstart_r: " << qstart_r << ';' << "qend_r: " << qend_r << '\n';
	//std::cerr << qstart_dist << ',' << qend_dist << '\n';

	//alignment does not span either start/end coord
	if(rightmost_r < start || leftmost_r > end){
		qstart_q = -1;
		qend_q = -1;
		msg.successful = false;
		msg.spanning_l = false;
		msg.spanning_r = false;
	}
	//region deleted
	else if(qstart_q > -1 && qend_q > -1 && qstart_q > qend_q){
		qstart_q = -1;
		qend_q = -1;
		msg.successful = true;
		msg.spanning_l = true;
		msg.spanning_r = true;
	} 
	else{
		msg.alignment_coords = std::make_pair(qstart_q, qend_q);
		//readjust if alignment is clipped
		if(leftmost_r > start && clipped_l && qstart_cigar_i == 1) {
			//continue expanding until end of query or end of cigar
			while(qstart_q > 0 && qstart_cigar_i > 0) {
				const int op = bam_cigar_op(cigar[qstart_cigar_i - 1]);
				const int ol = bam_cigar_oplen(cigar[qstart_cigar_i - 1]);
				if(op == BAM_CDEL) --qstart_cigar_i;
				else if(op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP || op == BAM_CINS){
					qstart_q -= ol;
					--qstart_cigar_i;
				}
				else break;
			}
		}
		//readjustment if alignment is clipped
		if(rightmost_r < end && clipped_r && qend_cigar_i == (alignment->core.n_cigar - 1)){
			//continue expanding until end of query or end of cigar
			while(qend_q < ((int)alignment->core.l_qseq - 1) && qend_cigar_i < alignment->core.n_cigar){
				const int op = bam_cigar_op(cigar[qend_cigar_i - 1]);
				const int ol = bam_cigar_oplen(cigar[qend_cigar_i - 1]);
				if(op == BAM_CDEL) ++qend_cigar_i;
				else if(op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP || op == BAM_CINS){
					qend_q += ol;
					++qend_cigar_i;
				}
				else break;
			}
		} 

		if(leftmost_q >= 0 && leftmost_r <= start) msg.spanning_l = true; else msg.spanning_l = false;
		if(rightmost_q >= 0 && rightmost_r >= end) msg.spanning_r = true; else msg.spanning_r = false;
		msg.successful = true;
	}

	if(msg.successful){
		//spanning both sides
		if(msg.spanning_l && msg.spanning_r) ptr.reset(new std::pair<int,int>(qstart_q, qend_q));
		//only left spanning
		else if(msg.spanning_l) ptr.reset(new std::pair<int,int>(qstart_q, alignment->core.l_qseq));
		//only right spanning
		else if(msg.spanning_r) ptr.reset(new std::pair<int,int>(0, qend_q));
		//neither spanning, take whole read
		else ptr.reset(new std::pair<int,int>(0, alignment->core.l_qseq));
	}

}

void parse_alignment(const int& rstart, const int& rend, bam1_t* read, ParsingStatus& msg, std::string& seq)
{
	std::unique_ptr<std::pair<int,int>> query_ptr;
	std::vector<int> projectedcoords;
	get_breakpoints(rstart, rend, read, msg, query_ptr);
	if(msg.successful){
		if((query_ptr->first == -1 && query_ptr->second != -1) || (query_ptr->first != -1 && query_ptr->second == -1)){
			std::cerr << '(' << antimestamp() << "): ERROR: unexpected querty start/end coords found for read " << (char*)read->data << '\n';
			exit(1);
		}
		if(query_ptr->first == -1 || (read->core.l_qseq < (query_ptr->second - query_ptr->first))) seq = "N"; 
		else {
			int l_qsubseq = query_ptr->second - query_ptr->first;
			int l_qsubseq_og = msg.alignment_coords.second - msg.alignment_coords.first;
			msg.alignment_coords.first = msg.alignment_coords.first - query_ptr->first;
			msg.alignment_coords.second = msg.alignment_coords.first +  l_qsubseq_og;
			//quality string
			uint8_t *q = bam_get_seq(read);
			//update array for read seqeunce
			for(int i = 0; i  < l_qsubseq; i++) seq += seq_nt16_str[bam_seqi(q, i + query_ptr->first)];
		}
	}	
}

void parse_alignments(const OtterOpts& params, const BED& bed, const BamInstance& bam_inst, AlignmentBlock& alignment_block)
{
	hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, bed.toBEDstring().c_str());
	if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " <<  bed.toBEDstring() << std::endl;
	else{
		while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
			if(bam_inst.read->core.qual >= params.mapq && (params.nonprimary || !(bam_inst.read->core.flag & BAM_FSECONDARY || bam_inst.read->core.flag & BAM_FSUPPLEMENTARY))){
				std::string name = (char*)bam_inst.read->data;
				std::string seq;
				ParsingStatus msg;
				parse_alignment(bed.start, bed.end, bam_inst.read, msg, seq);
				if(msg.successful){
					alignment_block.names.emplace_back(name);
					alignment_block.seqs.emplace_back(seq);
					alignment_block.statuses.emplace_back(msg);
				}
			}
		}
	}
	hts_itr_destroy(iter);
}

