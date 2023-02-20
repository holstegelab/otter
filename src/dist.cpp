#include "dist_matrix.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "kseq.h"
#include "otter_opts.hpp"
#include <zlib.h>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include <iostream>


KSEQ_INIT(gzFile, gzread)

class Region {
	public:
		std::vector<std::string> seqs;
		bool is_empty() const;
		double bipartite_dist(Region&, wfa::WFAligner&);
};

bool Region::is_empty() const {return seqs.empty();}

double Region::bipartite_dist(Region& that, wfa::WFAligner& aligner)
{
	DistMatrix distmatrix(seqs.size() + that.seqs.size());
	for(int i = 0; i < (int)seqs.size(); ++i){
		for(int j = 0; j < (int)that.seqs.size(); ++j){
			if(seqs[i].size() < that.seqs[j].size()) aligner.alignEnd2End(seqs[i], that.seqs[j]);
			else aligner.alignEnd2End(that.seqs[j],seqs[i]);
			distmatrix.set_dist(i, (int)seqs.size() + j, aligner.getAlignmentScore());
		}
	}

	std::vector<std::pair<int,int>> edges;

	for(int i = 0; i < (int)seqs.size(); ++i){
		int min_j = -1;
		double min_dist = 100000000.0;
		for(int j = 0; j < (int)that.seqs.size(); ++j){
			bool nodes_available = true;
			for(const auto& edge : edges) if(edge.first == i || edge.second == j) nodes_available = false;
			if(nodes_available){
				double dist = distmatrix.get_dist(i, j + (int)seqs.size());
				if(dist < min_dist) {
					min_j = j;
					min_dist = dist;
				}
			}
		}
		if(min_j > -1) edges.emplace_back(std::make_pair(i, min_j));
	}

	double total_sum = 0.0;
	int total_length = 0;

	for(const auto& e : edges) {
		int x_l = seqs[e.first].size();
		int y_l = that.seqs[e.second].size();
		total_sum += distmatrix.get_dist(e.first, e.second + (int)seqs.size());
		total_length += x_l > y_l ? x_l : y_l;
	}

	return total_sum / total_length;

}

void dist(const std::vector<std::string>& fastas, const OtterOpts& params)
{
	KS_FULL_COMMENT = true;
	std::map<std::string, int> regions;
	int regions_i = 0;
	for(uint32_t fasta_i = 0; fasta_i < fastas.size(); ++fasta_i){
		const std::string& fasta = fastas[fasta_i];
		//LOAD FASTA
		gzFile fp = gzopen(fasta.c_str(), "r");
		//sequence pointer
		kseq_t *seq = kseq_init(fp);
		int seq_l;
		while ((seq_l = kseq_read(seq)) >= 0 ) {
			std::string region = seq->name.s;
			if(regions.find(region) == regions.end()) regions.insert({region, regions_i++});
		}
		kseq_destroy(seq);
		gzclose(fp);
	}

	std::vector<std::string> regions_indexed(regions.size());
	for(const auto& r : regions) regions_indexed[r.second] = r.first;

	std::vector<std::vector<Region>> sample_regions(fastas.size());
	for(uint32_t fasta_i = 0; fasta_i < fastas.size(); ++fasta_i){
		std::vector<Region>& local_regions = sample_regions[fasta_i];
		local_regions.resize(regions.size());
		const std::string& fasta = fastas[fasta_i];
		//LOAD FASTA
		gzFile fp = gzopen(fasta.c_str(), "r");
		//sequence pointer
		kseq_t *seq = kseq_init(fp);
		int seq_l;
		while ((seq_l = kseq_read(seq)) >= 0 ) {
			std::string region = seq->name.s;
			std::string local_seq = seq->seq.s;
			local_regions[regions.find(region)->second].seqs.emplace_back(local_seq);
		}
		kseq_destroy(seq);
		gzclose(fp);
		for(auto& region : local_regions) {
			if(!region.seqs.empty() && (int)region.seqs.size() < params.max_alleles) {
				int i = params.max_alleles - region.seqs.size();
				while(i-- > 0) region.seqs.emplace_back(region.seqs.front());
			}
		}
	}

	/**
	for(int i = 0; i < (int)sample_regions.size();++i){
		std::cout << i << ":\n";
		for(int j = 0; j < (int)regions.size(); ++j){
			std::cout << '\t' << j << ":\n";
			const auto& local_regions = sample_regions[i][j];
			for(const auto& s : local_regions.seqs) std::cout << "\t\t" << s << '\n';

		}
	}
	*/

	wfa::WFAlignerEdit aligner(wfa::WFAligner::Score, wfa::WFAligner::MemoryMed);

	for(int i = 0; i < (int)regions.size(); ++i){
		for(int j = 0; j < (int)sample_regions.size(); ++j){
			for(int k = j + 1; k < (int)sample_regions.size(); ++k){
				if(!sample_regions[j][i].is_empty() || sample_regions[k][i].is_empty()){
					double dist = sample_regions[j][i].bipartite_dist(sample_regions[k][i], aligner);
					std::cout << fastas[j] << '\t' << fastas[k] << '\t' << regions_indexed[i] << '\t' << dist << '\n';
					std::cout << fastas[k] << '\t' << fastas[j] << '\t' << regions_indexed[i] << '\t' << dist << '\n';
				}
			}
		}
	}
}