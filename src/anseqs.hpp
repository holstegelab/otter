#ifndef ANSEQS_HPP
#define ANSEQS_HPP

#include "sam.h"
#include "anbamfilehelper.hpp"
#include "otter_opts.hpp"
#include "anbed.hpp"
#include <string>
#include <array>

extern const std::string rq_tag;
extern const std::string hp_tag;
extern const std::string ps_tag;
extern const std::string rg_tag;
extern const std::string pg_tag;
extern const std::string ta_tag;
extern const std::string ac_tag;
extern const std::string tc_tag;
extern const std::string sp_tag;


class ANSEQ{
	public:
		std::string seq;
		ANSEQ();
		ANSEQ(const std::string&);
};

class HAPLOTAG{
	public:
		int ps;
		int hp;
        HAPLOTAG();
        HAPLOTAG(int, int);
        bool operator==(const HAPLOTAG &x) const;
        bool operator!=(const HAPLOTAG &x) const;
        bool is_defined() const;
};

class ANALLELE : public ANSEQ{
	public:
		int scov;
		int acov;
		int tcov;
		float se;
		int ic;
		HAPLOTAG hpt;

		ANALLELE();
		ANALLELE(const std::string&);
		ANALLELE(const std::string&, int, int, int, float, int, int, int);
		void stdout_sam(const std::string&, const std::string&, const int&, const int&, const std::string&, const bool& is_read = false, const bool& is_spanning_l = false, const bool& is_spanning_r = false) const;
		void stdout_fa(const std::string&, const std::string&, const bool& is_read = false, const bool& is_spanning_l = false, const bool& is_spanning_r = false) const;
};

class ANREAD : public ANSEQ{
	public:
		std::string name;
		double rq;
		//read is spanning: 0 (none); 1 (left); 2 (right)
		bool is_spanning_l;
		bool is_spanning_r;
		HAPLOTAG hpt;
		//left/right most coords before clipping
		std::pair<int,int> ccoords;

		ANREAD();
		ANREAD(const std::string&, const std::string&);
		ANREAD(const std::string&, const std::string&, double, bool, bool, int, int);
		bool is_spanning() const;
		void set_is_spanning();
		void stdout_sam(const std::string&, const int&, const int&, const std::string&) const;
		void stdout_fa(const std::string&) const;


};

class KUSAGE{
	public:
		std::vector<double> vec;
		double vnorm;
		KUSAGE(const std::vector<double>&);
		double cosine_sim(const KUSAGE&) const;
		//hill-shannon diversity
		double hsdiv() const;
};

class KmerEncoding{
	public:
		std::array<uint8_t, 256> nt2encoding;
		std::array<char, 4> encoding2nt;
		KmerEncoding();
		uint64_t kmer2index(const char*, uint32_t, const uint32_t&) const;
		void index2kmer(std::string&, uint64_t, uint32_t) const;
	private:
		uint64_t _kmer2index(const char*, uint32_t, const uint32_t&) const;
};

void parse_anreads(
	const OtterOpts& params, 
	const BED& bed, 
	const BamInstance& bam_inst, 
	std::vector<ANREAD>& anread_block
);

void parse_anallele(
	const std::string& target_region, 
	const std::map<std::string, int>& sample2index, 
	bam1_t* read, 
	std::vector<ANALLELE>& anallele_block, 
	std::vector<int>& allele_sample_indeces
);

void parse_analleles(
	const OtterOpts& params, 
	const BamInstance& bam_inst, 
	const BED& bed, 
	const std::map<std::string, int>& sample2index, 
	std::vector<ANALLELE>& anallele_block, 
	std::vector<int>& allele_sample_indeces
);

void seq2kcounts(
	const uint32_t& k,
	const KmerEncoding& encoding, 
	const std::string& seq, 
	std::vector<double>& kmercounts
);


#endif