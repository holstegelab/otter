#include "bindings/cpp/WFAligner.hpp"
#include "antimestamp.hpp"
#include "anbed.hpp"
#include "otter_opts.hpp"
#include "anbamdb.hpp"
#include "BS_thread_pool.hpp"
#include "anbamfilehelper.hpp"
#include "anfahelper.hpp"
#include "anseqs.hpp"
#include "andistmat.hpp"
#include "otterclust.hpp"
#include <iostream>
#include <vector>
#include <array>

void output_vcf_header(const std::string& bam, const std::vector<std::string>& sample_index, const std::string& ref_name)
{
	BamInstance bam_inst;
	bam_inst.init(bam, true);
	std::cout << "##fileformat=VCFv4.2\n";
	for(int i = 0; i < bam_inst.header->n_targets; ++i) std::cout << "##contig=<ID=" << bam_inst.header->target_name[i] << ",length=" << bam_inst.header->target_len[i] << ">\n";
	bam_inst.destroy();
	std::cout << 
	//"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Coverage of Region\">\n" <<
	"##INFO=<ID=HSD,Number=R,Type=Float,Description=\"Hill-Shannon Diversity Metric\">\n" <<
	"##ALT=<ID=DEL,Description=\"Deletion\">\n" <<
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" <<
	"##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">\n" <<
	"##FORMAT=<ID=HP,Number=1,Type=Integer,Description=\"Haplotype Identifier\">\n" <<
	"##FORMAT=<ID=TC,Number=1,Type=Integer,Description=\"Total Coverage of Region\">\n" <<
	"##FORMAT=<ID=AC,Number=2,Type=Integer,Description=\"Total Coverage For Each Allele\">\n" <<
	"##FORMAT=<ID=SC,Number=2,Type=Integer,Description=\"Total Coverage of Spanning Reads For Each Allele\">\n" <<
	"##FORMAT=<ID=SE,Number=2,Type=Float,Description=\"Standard Mean Error of Spanning Reads For Each Allele\">\n";
	std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for(const auto& sample : sample_index) if(sample != ref_name) std::cout << '\t' << sample;
	std::cout << '\n';
}

void output_vcf_line(const int& offset_l, const int& offset_r, const BED& region, const SampleIndex& si, const int& ref_allele_index, const std::vector<ANALLELE>& alleles, const std::vector<Genotype>& genotypes, const std::vector<int>& reps, const std::vector<std::optional<std::pair<int,int>>>& sample2localindeces)
{
	std::cout << region.chr << '\t' << (1 + region.start - offset_l) << '\t' << region.toScString() << '\t' << alleles[ref_allele_index].seq << '\t';
	if(reps.size() == 1) std::cout << '.';
	else{
		for(uint32_t i = 1; i < reps.size(); ++i) {
			if(i > 1) std::cout << ',';
			if(alleles[reps[i]].seq == "N") std::cout << "<DEL>"; else std::cout << alleles[reps[i]].seq;
		}
	}
	std::cout << "\t.\t.\tHSD=";
	for(uint32_t i = 0; i < reps.size(); ++i) {
		if(i > 0) std::cout << ',';
		std::cout << genotypes[reps[i]].hsd;
	}
	std::cout << "\tGT:PS:HP:TC:AC:SC:SE";
	for(uint32_t i = 0; i < sample2localindeces.size() - 1; ++i) {
		const auto& sample_indeces = sample2localindeces[i];
		if(!sample_indeces.has_value()) std::cout << "\t./.:.:.:.:.:.:.";
		else {
			const auto& allele_1 = alleles[sample_indeces.value().first];
			const auto& allele_2 = alleles[sample_indeces.value().second];
			if(allele_1.hpt != allele_2.hpt) std::cerr << "(" << antimestamp() << "): [WARNING] mismatching phased information for " << si.index2sample[i] << ": allele1=PS:" << allele_1.hpt.ps << ":HP:" << allele_1.hpt.hp << " allele2=PS:" << allele_1.hpt.ps << ":HP:" << allele_1.hpt.hp << '\n';
			std::cout << '\t' << genotypes[sample_indeces.value().first].gt << '/' << genotypes[sample_indeces.value().second].gt << ':' << allele_1.hpt.ps << ':' << allele_1.hpt.hp << ':' << allele_1.tcov << ':' << allele_1.acov << ',' << allele_2.acov << ':' << allele_1.scov << ',' << allele_2.scov << ':' << allele_1.se << ',' << allele_2.se;
		}
	}
	std::cout << '\n';

}

void genotype_process(const OtterOpts& params, const std::string& bam, const std::vector<BED>& regions, const std::string& reference, const SampleIndex& si, const int& refindex)
{
	BS::thread_pool pool(params.threads);
	std::mutex stdout_mtx;
	pool.parallelize_loop(0, regions.size(),
		[&pool, &stdout_mtx, &params, &bam, &regions, &reference, &si, &refindex](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);
			FaidxInstance faidx_inst;
			if(!reference.empty()) faidx_inst.init(reference);
			for(int region_i = a; region_i < b; ++region_i) {
				const BED& region = regions[region_i];
				//std::cerr << region.toScString() << std::endl;
				std::vector<ANALLELE> anallele_block;
				std::vector<int> allele_sample_indeces;
				parse_analleles(params, bam_inst, region, si.sample2index, anallele_block, allele_sample_indeces);
				if(anallele_block.size() != allele_sample_indeces.size()){
					std::cerr << "(" << antimestamp() << "): [ERROR] expected matching total number of alleles and samples: " << anallele_block.size() << " vs " << allele_sample_indeces.size() << std::endl;
					exit(1);
				}
				if(anallele_block.empty()) std::cerr << "(" << antimestamp() << "): [WARNING] no alleles found for " << region.toScString() << std::endl;
				else{
					int ref_allele_index = -1;
					{
						std::string refseq;
						faidx_inst.fetch(region.chr, region.start - si.offset_l, region.end + si.offset_r - 1, refseq);
						ref_allele_index = allele_sample_indeces.size();
						allele_sample_indeces.emplace_back(refindex);
						anallele_block.emplace_back(refseq);

					}
					std::vector<std::optional<std::pair<int,int>>> sample2localindeces(si.sample2index.size());
					for(int i = 0; i < (int)allele_sample_indeces.size(); ++i){
						auto& pair_indeces = sample2localindeces[allele_sample_indeces[i]];
						if(!pair_indeces.has_value()) pair_indeces.emplace(std::make_pair(i, i));
						else {
							if(i < pair_indeces.value().first) pair_indeces.value().first = i;
							else if(i > pair_indeces.value().second) pair_indeces.value().second = i;
						}
					}
					/**
					for(uint32_t i = 0; i < si.sample2index.size(); ++i) {
						const auto& local = sample2localindeces[i];
						if(local.has_value()) std::cout << i << '\t' << local.value().first << '\t' << local.value().second << '\n';
						else std::cout << "NO LOCAL: " << si.index2sample[i] << std::endl;
					}
					
					std::cout << anallele_block.size() << '\t' << allele_sample_indeces.size() << '\n';
					for(int i = 0; i < (int)anallele_block.size(); ++i){
						std::cout << si.index2sample[allele_sample_indeces[i]] << '\t' << anallele_block[i].se << '\t' << anallele_block[i].seq << '\n';
					}
					*/
					std::vector<Genotype> genotypes(anallele_block.size());
					std::vector<int> gt_reps;
					int acc_gt = anallele_cluster(params.max_error, params.max_cosdis, anallele_block, genotypes, gt_reps);
					
					if(acc_gt != (int)gt_reps.size()){
						std::cerr << "(" << antimestamp() << "): ERROR unexpected representative alleles (" << gt_reps.size() << ") for " << acc_gt << " total alleles" << std::endl;
						exit(1);
					}
					//std::cout << "total gts: " << acc_gt << " total reps: " << gt_reps.size() << std::endl;
					//for(uint32_t i = 0; i < anallele_block.size(); ++i) std::cout << i << '\t' << si.index2sample[allele_sample_indeces[i]] << '\t' << genotypes[i] << '\t' << anallele_block[i].seq << '\n';
					int ref_gt = genotypes[ref_allele_index].gt;
					//std::cout << "ref_gt: "  << ref_gt << '\n';
					
					std::vector<int> gt_reps_centered = gt_reps;
					for(int i = 0; i < (int)gt_reps_centered.size(); ++i){
						if(i == 0) gt_reps_centered[0] = ref_allele_index;
						else if(i <= ref_gt) gt_reps_centered[i] = gt_reps[i-1];
					}
					
					for(uint32_t i = 0; i < anallele_block.size(); ++i){
						if(genotypes[i].gt == ref_gt) genotypes[i].gt = 0;
						else if(genotypes[i].gt < ref_gt) ++genotypes[i].gt;
					}
					
					//for(uint32_t i = 0; i < anallele_block.size(); ++i) std::cout << i << '\t' << si.index2sample[allele_sample_indeces[i]] << '\t' << genotypes[i].gt << '\t' << genotypes[i].gt_l << '\t' << genotypes[i].gt_k << '\t' << anallele_block[i].seq << std::endl;
					stdout_mtx.lock();
					output_vcf_line(si.offset_l, si.offset_r, region, si, ref_allele_index, anallele_block, genotypes, gt_reps_centered, sample2localindeces);
					stdout_mtx.unlock();
				}
			}
			faidx_inst.destroy();
			bam_inst.destroy();
	}).wait();

}

void genotype(OtterOpts& params, const std::string& bam, const std::string& bed, const std::string& reference)
{
	std::string refname = "OTTER_INTREF";
	std::vector<BED> regions;
	parse_bed_file(bed, regions);
	
	SampleIndex si;
	si.init(bam);
	std::cerr << '(' << antimestamp() << "): Found " << si.index2sample.size() << " samples (read-group tags)\n";
	std::cerr << '(' << antimestamp() << "): Using offset of " << si.offset_l << ',' << si.offset_r << '\n' ;
	int refindex = -1;
	{	
		refindex = si.index2sample.size();
		si.index2sample.emplace_back(refname);
		si.sample2index[refname] = refindex;
	}
	if(!params.is_fa) output_vcf_header(bam, si.index2sample, refname);

	genotype_process(params, bam, regions, reference, si, refindex);
}