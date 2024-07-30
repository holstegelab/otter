#include "anseqs.hpp"
#include "antimestamp.hpp"
#include "sam.h"
#include <iostream>
#include <string>
#include <cmath>
#include <array>

const std::string ps_tag = "PS";
const std::string hp_tag = "HP";
const std::string rq_tag = "rq";
const std::string rg_tag = "RG";
const std::string ta_tag = "ta";
const std::string tc_tag = "tc";
const std::string ac_tag = "ac";
const std::string sc_tag = "sc";
const std::string se_tag = "se";
const std::string sp_tag = "sp";
const std::string ic_tag = "ic";

char _get_spanning_tag_value(const bool& is_spanning_l, const bool& is_spanning_r)
{
	if(is_spanning_l && is_spanning_r) return('b');
	else if(is_spanning_l) return('l');
	else if(is_spanning_r) return('r');
	else return('n');
}

/************************
 *** ANSEQ defintions ***
 ************************/
ANSEQ::ANSEQ(){};
ANSEQ::ANSEQ(const std::string& _s): seq(_s){};

/************************
 *** ANALLELE defintions ***
 ************************/
ANALLELE::ANALLELE(){};
ANALLELE::ANALLELE(const std::string& _seq): ANSEQ(_seq), scov(1), acov(1), tcov(1), se(0), ic(1), hpt(-1,-1){};
ANALLELE::ANALLELE(const std::string& _seq, int s, int a, int t, float _se, int i, int h, int p): ANSEQ(_seq), scov(s), acov(a), tcov(t), se(_se), ic(i), hpt(h,p){};

void ANALLELE::stdout_sam(const std::string& name, const std::string& chr, const int& start, const int& end, const std::string& rg, const bool& is_read, const bool& is_spanning_l, const bool& is_spanning_r) const
{
	std::string pseudo_qual(seq.size(), '!');
	std::cout << name << "\t0\t" << chr << '\t' << start << "\t0\t" << seq.size() << "M\t*\t0\t0\t" << seq << '\t' << pseudo_qual;
	if(!rg.empty()) std::cout << '\t' << rg_tag << ":Z:" << rg;
	std::cout  << '\t' << ta_tag << ":Z:" << chr << ':' << start << '-' << end << '\t' << tc_tag << ":i:" << tcov << '\t' << ac_tag << ":i:" << acov << '\t' << sc_tag << ":i:" << scov;
	if(is_read) std::cout << '\t' << sp_tag << ":A:" << _get_spanning_tag_value(is_spanning_l, is_spanning_r);
	std::cout << '\t' << ic_tag << ":i:" << ic;
	std::cout << '\t' << se_tag << ":f:" << se;
	if(hpt.ps >= 0) std::cout << '\t' << ps_tag << ":i:" << hpt.ps;
	if(hpt.hp >= 0) std::cout << '\t' << hp_tag << ":i:" << hpt.hp;
	std::cout << '\n';
}

void ANALLELE::stdout_fa(const std::string& name, const std::string& region, const bool& is_read, const bool& is_spanning_l, const bool& is_spanning_r) const
{
	std::cout << '>' << name << '#' << region << '#' << tc_tag << ":i:" << tcov << '#' << ac_tag << ":i:" << acov << '#' << sc_tag << ":i:" << scov;
	if(is_read) std::cout << '#' << sp_tag << ":A:" << _get_spanning_tag_value(is_spanning_l, is_spanning_r);
	if(hpt.ps >= 0) std::cout << '#' << ps_tag << ":i:" << hpt.ps;
	if(hpt.hp >= 0) std::cout << '#' << hp_tag << ":i:" << hpt.hp;
	std::cout << '\n' << seq << '\n';
}

/************************
 *** HAPLOTAG defintions ***
 ************************/
HAPLOTAG::HAPLOTAG(): ps(-1), hp(-1){};
HAPLOTAG::HAPLOTAG(int _ps, int _hp): ps(_ps),hp(_hp){};
bool HAPLOTAG::operator==(const HAPLOTAG &x) const { return ps == x.ps && hp == x.hp;}
bool HAPLOTAG::operator!=(const HAPLOTAG &x) const { return ps != x.ps || hp != x.hp;}
bool HAPLOTAG::is_defined() const { return ps >= 0 && hp >= 0; }

/************************
 *** ANREAD defintions ***
 ************************/
ANREAD::ANREAD(): rq(0), is_spanning_l(false), is_spanning_r(false){};
ANREAD::ANREAD(const std::string& _seq, const std::string& _name): ANSEQ(_seq), name(_name), rq(0), is_spanning_l(false), is_spanning_r(false){};
ANREAD::ANREAD(const std::string& _seq, const std::string& _name, double _rq, bool _sl, bool _sr, int _ps, int _hp): ANSEQ(_seq), name(_name), rq(_rq), is_spanning_l(_sl), is_spanning_r(_sr), hpt(_ps, _hp){};
bool ANREAD::is_spanning() const { return is_spanning_l && is_spanning_r;}
void ANREAD::set_is_spanning() {is_spanning_l = true; is_spanning_r = true;}

void ANREAD::stdout_sam(const std::string& chr, const int& start, const int& end, const std::string& rg) const
{
	std::string pseudo_qual(seq.size(), '!');
	std::cout << name << "\t0\t" << chr << '\t' << start << "\t0\t" << seq.size() << "M\t*\t0\t0\t" << seq << '\t' << pseudo_qual;
	if(!rg.empty()) std::cout << '\t' << rg_tag << ":Z:" << rg;
	std::cout  << '\t' << ta_tag << ":Z:" << chr << ':' << start << '-' << end << '\t' << sp_tag << ":A:";
	if(is_spanning_l && is_spanning_r) std::cout << 'b';
	else if(is_spanning_l) std::cout << 'l';
	else if(is_spanning_r) std::cout << 'r';
	else std::cout << 'n';
	if(hpt.ps >= 0) std::cout << '\t' << ps_tag << ":i:" << hpt.ps;
	if(hpt.hp >= 0) std::cout << '\t' << hp_tag << ":i:" << hpt.hp;
	std::cout << '\t' << rq_tag << ":f:" << rq;
	std::cout << '\n';
}

void ANREAD::stdout_fa(const std::string& region) const
{
	std::cout << '>' << name << '#' << region;
	std::cout << '#' << sp_tag << ":A:" << _get_spanning_tag_value(is_spanning_l, is_spanning_r);
	if(hpt.ps >= 0) std::cout << '#' << ps_tag << ":i:" << hpt.ps;
	if(hpt.hp >= 0) std::cout << '#' << hp_tag << ":i:" << hpt.hp;
	std::cout << '\n' << seq << '\n';
}

/**************************
 *** KMERVEC defintions ***
 **************************/
KUSAGE::KUSAGE(const std::vector<double>& kcounts): vec(kcounts.size(), 0), vnorm(0)
{
	int total_counts = 0;
	for(const auto& k : kcounts) total_counts += k;
	for(uint32_t i = 0; i < vec.size(); ++i) {
		double value = kcounts[i]/total_counts;
		vec[i] = value;
		vnorm += value*value;
	}
	vnorm = std::sqrt(vnorm);
}

double KUSAGE::cosine_sim(const KUSAGE& that) const
{
	if(vec.size() != that.vec.size()){
		std::cerr << '(' << antimestamp() << "): ERROR: unexpected lengths for two KUSAGE vectors: " << vec.size() << " vs " << that.vec.size() << std::endl;
			exit(1);
	}
	double x_dot_y = 0;
	for(uint32_t i = 0; i < vec.size(); ++i) x_dot_y += vec[i]*that.vec[i];
	return x_dot_y/(vnorm * that.vnorm);
}

double KUSAGE::hsdiv() const
{
	/**
	double acc = 0;
	for(const auto& ku : vec) acc += (ku*ku);
	return acc;
	*/

	double acc = 0;
	for(const auto& ku : vec) if(ku > 0) acc += (ku*std::log(ku));
	acc = -1*acc;
	return std::pow(M_E, acc);

}

void seq2kcounts(const uint32_t& k, const KmerEncoding& encoding, const std::string& seq, std::vector<double>& kmercounts)
{
	uint32_t max_index = (int)std::pow(4, k);
	kmercounts.resize(max_index + 1, 0);
	if(seq.size() >= k){
		for(uint32_t j = 0; j < seq.size() - k + 1; ++j){
			bool is_valid = true;
			for(uint32_t h = 0; h < k; ++h) {
				if(encoding.nt2encoding[*(seq.c_str() + j + h)] == 4) {
					is_valid = false;
					break;
				}
			}
			uint64_t k_j = is_valid ? encoding.kmer2index(seq.c_str() + j, k, max_index) : max_index;
			++kmercounts[k_j];
		}
	}
}

/*******************************
 *** KmerEncoding defintions ***
 *******************************/
KmerEncoding::KmerEncoding()
{
	for(int i = 0; i < 256; ++i) nt2encoding[i] = 4;
	nt2encoding['A'] = nt2encoding['a'] = 0; 
	nt2encoding['C'] = nt2encoding['c'] = 1;
	nt2encoding['G'] = nt2encoding['g'] = 2; 
	nt2encoding['T'] = nt2encoding['t'] = 3;

	for(int i = 0; i < 4; ++i) encoding2nt[i] = 'N';
    encoding2nt[0] = 'A';
    encoding2nt[1] = 'C';
    encoding2nt[2] = 'G';
    encoding2nt[3] = 'T';
}

uint64_t KmerEncoding::kmer2index(const char* nt, uint32_t k, const uint32_t& max_index) const {return _kmer2index(nt+k-1, k, max_index);}

void KmerEncoding::index2kmer(std::string& seq, uint64_t index, uint32_t k) const
{
	if(k == 1) {
	    uint64_t nt_index = encoding2nt[(uint8_t)index];
	    seq[k-1] = nt_index;
	}
	else {
		uint64_t prefixindex = index / 4;
		uint8_t r = (uint8_t)(index % 4);
		uint64_t nt_index = encoding2nt[(uint8_t)r];
		seq[k-1] = nt_index;
		index2kmer(seq, prefixindex, k - 1);
	}
}

uint64_t KmerEncoding::_kmer2index(const char* nt, uint32_t k, const uint32_t& max_index) const
{
	if(k == 0) return 0;
	uint64_t index = (uint64_t)nt2encoding[(uint8_t)(*nt)];
	return 4 * _kmer2index(nt - 1, k - 1, max_index) + index;
}

/***************************
 *** PARSEMSG defintions ***
 ***************************/
/**
 * PARSEMSG class and definitions 
 * 		Used to store the status when parsing the corresponding subsequence of an alignment.
 * 		Locally defined within this source file.
 */
class PARSEMSG {
    public:
    	bool successful;
    	bool spanning_l;
    	bool spanning_r;
        std::pair<int, int> alignment_coords;
    	PARSEMSG();
    	bool is_spanning() const;
    	void transfer_status(ANREAD&) const;
};

PARSEMSG::PARSEMSG(): successful(true),spanning_l(true), spanning_r(true){alignment_coords = std::make_pair(-1,-1);};

bool PARSEMSG::is_spanning() const {return spanning_l && spanning_r;};

void PARSEMSG::transfer_status(ANREAD& anread) const
{
	if(is_spanning()) anread.set_is_spanning();
	else if(spanning_l) anread.is_spanning_l = true;
	else if(spanning_r)anread.is_spanning_r = true;
	anread.ccoords = alignment_coords;
}

/*****************************************
 *** Parse ALIGNMENT BAM AUX functions ***
 *****************************************/
void parse_standard_auxs(bam1_t* alignment, ANREAD& anread)
{
	auto aux_ptr = bam_aux_get(alignment, hp_tag.c_str());
	if(aux_ptr != NULL) anread.hpt.hp = bam_aux2i(aux_ptr);
	aux_ptr = bam_aux_get(alignment, ps_tag.c_str());
	if(aux_ptr != NULL) anread.hpt.ps = bam_aux2i(aux_ptr);
	aux_ptr = bam_aux_get(alignment, rq_tag.c_str());
	if(aux_ptr != NULL) anread.rq = bam_aux2f(aux_ptr);
}

/**
 * 
 */
void project_positions(bam1_t* alignment, bool& clipped_l, bool& clipped_r, std::vector<int>& query_coords)
{
	query_coords.resize(alignment->core.l_qseq, -1);
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
				query_coords[qpos] = rpos;
				++rpos;
				++qpos;
			}
		}
		else if(op == BAM_CINS) qpos += ol;
		else if(op == BAM_CDEL) rpos += ol;
	}
}
/**
 * 
 */
void get_breakpoints(const int& start,  const int& end, bam1_t* alignment, PARSEMSG& msg, std::unique_ptr<std::pair<int, int>>& ptr)
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
/**
 * 
 */
void parse_alignment(const int& rstart, const int& rend, bam1_t* read, PARSEMSG& msg, std::string& seq)
{
	std::unique_ptr<std::pair<int,int>> query_ptr;
	std::vector<int> projectedcoords;
	get_breakpoints(rstart, rend, read, msg, query_ptr);
	if(msg.successful){
		if((query_ptr->first == -1 && query_ptr->second != -1) || (query_ptr->first != -1 && query_ptr->second == -1)){
			std::cerr << '(' << antimestamp() << "): ERROR: unexpected querty start/end coords found for read " << (char*)read->data << std::endl;
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
			if(seq.size() == 0) seq = "N";
		}
	}	
}
/**
 * 
 */
void parse_anreads(const OtterOpts& params, const BED& bed, const BamInstance& bam_inst, std::vector<ANREAD>& anread_block)
{
	hts_itr_t *iter = bam_itr_querys(bam_inst.idx, bam_inst.header, bed.toScString().c_str());
	if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " <<  bed.toScString() << std::endl;
	else{
		while(bam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
			if(bam_inst.read->core.qual >= params.mapq && (params.nonprimary || !(bam_inst.read->core.flag & BAM_FSECONDARY || bam_inst.read->core.flag & BAM_FSUPP))){
				ANREAD anread;
				anread.name = (char*)bam_inst.read->data;
				PARSEMSG msg;
				
				parse_alignment(bed.start, bed.end, bam_inst.read, msg, anread.seq);
				if(msg.successful && (!params.omitnonspanning || msg.is_spanning())){
					msg.transfer_status(anread);
					parse_standard_auxs(bam_inst.read, anread);
					if(anread.rq >= params.read_quality) anread_block.emplace_back(anread);
				}
			}
		}
	}
	bam_itr_destroy(iter);
}

void parse_anallele(const std::string& target_region, const std::map<std::string, int>& sample2index, bam1_t* read, std::vector<ANALLELE>& anallele_block, std::vector<int>& allele_sample_indeces)
{
	//parse region
	std::string parsed_region;
	auto aux_ptr = bam_aux_get(read, ta_tag.c_str());
	if(aux_ptr != NULL) parsed_region = (char*)(aux_ptr+1);
	if(target_region == parsed_region){
		//parse sample-name
		aux_ptr = bam_aux_get(read, rg_tag.c_str());
		std::string sample;
		if(aux_ptr != NULL) sample = (char*)(aux_ptr+1);
		auto it_sample_index = sample2index.find(sample);
		if(it_sample_index == sample2index.end()){
			std::cerr << "(" << antimestamp() << "): ERROR unrecognized sample-name (read-group): " << sample << std::endl;
			exit(1);
		}
		/**parse coverage fields**/
		int tc = 1, ac = 1, sc = 1;
		//total coverage
		aux_ptr = bam_aux_get(read, tc_tag.c_str());
		if(aux_ptr != NULL) tc = bam_aux2i(aux_ptr);
		//allele coverage
		aux_ptr = bam_aux_get(read, ac_tag.c_str());
		if(aux_ptr != NULL) ac = bam_aux2i(aux_ptr);
		//spanning coverage
		aux_ptr = bam_aux_get(read, sc_tag.c_str());
		if(aux_ptr != NULL) sc = bam_aux2i(aux_ptr);
		/**parse haplotag fields**/
		HAPLOTAG hpt;
		aux_ptr = bam_aux_get(read, ps_tag.c_str());
		if(aux_ptr != NULL) hpt.ps = bam_aux2i(aux_ptr);
		aux_ptr = bam_aux_get(read, hp_tag.c_str());
		if(aux_ptr != NULL) hpt.hp = bam_aux2i(aux_ptr);
		/**se field**/
		float se = 0.0;
		aux_ptr = bam_aux_get(read, se_tag.c_str());
		if(aux_ptr != NULL) se = bam_aux2f(aux_ptr);
		/**ic field**/
		int ic = 1;
		aux_ptr = bam_aux_get(read, ic_tag.c_str());
		if(aux_ptr != NULL) ic = bam_aux2i(aux_ptr);
		/**parse sequence**/
		uint8_t *q = bam_get_seq(read);
	    uint32_t l_qseq = read->core.l_qseq > 0 ? read->core.l_qseq : 1;
	    std::string seq(l_qseq, 'N');
	    for(int i = 0; i  < read->core.l_qseq; i++) seq[i] = seq_nt16_str[bam_seqi(q, i)];
	    allele_sample_indeces.emplace_back(it_sample_index->second);
	    anallele_block.emplace_back(seq, sc, ac, tc, se, ic, hpt.ps, hpt.hp);
	}
}

void parse_analleles(const OtterOpts& params, const BamInstance& bam_inst, const BED& bed, const std::map<std::string, int>& sample2index, std::vector<ANALLELE>& anallele_block, std::vector<int>& allele_sample_indeces)
{
	hts_itr_t *iter = bam_itr_querys(bam_inst.idx, bam_inst.header, bed.toScString().c_str());
	if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " <<  bed.toScString() << std::endl;
	else{
		while(bam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
			ANALLELE anallele;
			parse_anallele(bed.toScString(), sample2index, bam_inst.read, anallele_block, allele_sample_indeces);
		}
	}
	bam_itr_destroy(iter);
}

