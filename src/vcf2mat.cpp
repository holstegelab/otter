#include "antimestamp.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "anbamdb.hpp"
#include "BS_thread_pool.hpp"
#include "anbamfilehelper.hpp"
#include "anseqs.hpp"
#include "angzipiter.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>

void line2columns(const char delim, const std::string& line, std::vector<std::string>& columns)
{
  std::string value;
  std::istringstream stream(line);
  while(std::getline(stream, value, delim)) columns.emplace_back(value);
}

void parse_alleles(const std::string& line, std::string& region, std::vector<std::string>& alleles)
{
	std::string column;
  std::istringstream stream(line);
  int index = 0;
  while(std::getline(stream, column, '\t')) {
  	if (index == 2) region = column;
  	else if(index == 3) alleles.emplace_back(column);
  	else if(index == 4 && column != ".") {
  		if(column == "<DEL>") alleles.emplace_back("N"); else line2columns(',', column, alleles);
  	}
  	++index;
  }
}

double get_gc_content(const KmerEncoding& encoding, const std::string& seq)
{
	double gc = 0;
	for(const auto& nt : seq) {
		uint8_t e = encoding.nt2encoding[nt];
		if(e == 1 || e == 2) gc += 1;
	}
	return gc/seq.size();
}

void vcf2mat(const OtterOpts& params, const std::string& bed, const std::string& vcf, const int& k_l)
{
	std::vector<BED> regions;
	parse_bed_file(bed, regions);
	
	std::mutex stdout_mtx;
	KmerEncoding encoding;
	GZIPiter vcfiter(vcf);
	while(vcfiter.hasNext()){
		std::string line;
		vcfiter.next(line);
		//std::cout << line << '\n';
		if(line.front() != '#'){
			std::vector<std::string> alleles;
			std::string region;
			parse_alleles(line, region, alleles);
			//std::vector<double> means;
			//std::vector<double> sds;
			for(uint32_t i = 0; i < alleles.size(); ++i){
				std::vector<double> kcounts;
				seq2kcounts(k_l, encoding, alleles[i], kcounts);
				KUSAGE kusage(kcounts);
				std::cout << region << '\t' << i << '\t' << get_gc_content(encoding, alleles[i]) << '\t' << alleles[i].size() << '\t' << kusage.hsdiv();
				for(const auto& ku : kusage.vec) std::cout << '\t' << ku;
				std::cout << '\n';
			}
		}
	}
	vcfiter.close();
}