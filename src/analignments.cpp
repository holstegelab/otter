#include "bindings/cpp/WFAligner.hpp"
#include "analignments.hpp"
#include "anseqs.hpp"
#include "anfahelper.hpp"
#include "andistmat.hpp"
#include "anppoa.hpp"
#include <cmath>
#include <vector>
#include <string>

void local_realignment(const std::string& chr, const int& start, const int& end, const int& flank, const double& min_sim, const FaidxInstance& faidx_inst, std::vector<ANREAD>& reads, wfa::WFAligner& aligner)
{
	std::string ref_left;
	std::string ref_right;
	for(int i = 0; i < (int)reads.size(); ++i){
		auto& local_read = reads[i];
		if(!local_read.is_spanning() && (local_read.is_spanning_l || local_read.is_spanning_r)){
			bool left_realignment = local_read.is_spanning_r && local_read.ccoords.first >= flank;
			bool right_realignment = local_read.is_spanning_l && (int)local_read.seq.size() - local_read.ccoords.second >= flank;
			std::string subseq;
			if(left_realignment) {
				if(ref_left.empty()) faidx_inst.fetch(chr, start - flank, start, ref_left);
				subseq = local_read.seq.substr(0, local_read.ccoords.first);
				//aligner.alignEnd2End(ref_left, subseq);
				aligner.alignEnd2End(subseq, ref_left);
			}
			else if(right_realignment){
				if(ref_right.empty()) faidx_inst.fetch(chr, end, end + flank, ref_right);
				subseq = local_read.seq.substr(local_read.ccoords.second);
				//aligner.alignEnd2End(ref_right, subseq);
				aligner.alignEnd2End(subseq, ref_right);
			}

			if(!subseq.empty()){
				std::vector<int> scores(subseq.size(), 0);
				int j = 0;
				for(const auto& op : aligner.getAlignmentCigar()){
					if(op != 'I'){
						int penalty = op == 'M' ? 1 : -1;
						if(penalty > 0){
							if(j == 0) scores[j] = penalty;
							else scores[j] = scores[j-1] + penalty;
						}
						else if(j > 0 && scores[j-1] > 0) scores[j] = scores[j-1] + penalty;
						++j;
					}
				}
				int max_sum_i = 0;
				for(j = 0; j < (int)scores.size(); ++j) if(scores[j] > scores[max_sum_i]) max_sum_i = j;
				int start_i = max_sum_i;
				while(start_i > 0 && scores[start_i] > 0) --start_i;
				if((scores[max_sum_i] / (double)flank) >= min_sim){
					if(left_realignment) local_read.seq = local_read.seq.substr(max_sum_i);
					else if(right_realignment) local_read.seq = local_read.seq.substr(0, local_read.ccoords.second + start_i);
					local_read.set_is_spanning();
				}
			}
		}
	}
}

double align_anreads(wfa::WFAligner& aligner, ANREAD& read_x, ANREAD& read_y)
{
	if(read_x.seq == read_y.seq) return 0.0;
	else if((read_x.is_spanning() && read_y.is_spanning()) || (read_y.is_spanning() && read_x.seq.size() >= read_y.seq.size())){
		//std::cout << "now here again\n" << std::endl;
		bool x_is_smallest = read_x.seq.size() < read_y.seq.size();
		double largest = x_is_smallest ? (double)read_y.seq.size(): (double)read_x.seq.size();
		//if(x_is_smallest) aligner.alignEnd2End(read_x.seq, read_y.seq); else aligner.alignEnd2End(read_y.seq, read_x.seq);
		if(x_is_smallest) aligner.alignEnd2End(read_y.seq, read_x.seq); else aligner.alignEnd2End(read_x.seq, read_y.seq);
		int dist = aligner.getAlignmentScore();
		return dist/largest;
		/**
		if(read_x.is_spanning() && read_y.is_spanning()) return dist/largest;
		//correct for non-spanning situations
		else {
			int size_diff = x_is_smallest ? (int)read_y.seq.size() - (int)read_x.seq.size() : (int)read_x.seq.size() - (int)read_y.seq.size();
			if(dist >= size_diff) return dist(dist - size_diff)/largest;
			else 
		}
		*/
	}
	else if(read_y.is_spanning()){
		//std::cout << "now here\n" << std::endl;
		int length_diff = read_y.seq.size() - read_x.seq.size();
		if(length_diff < 0){
			length_diff = -1 * length_diff;
			if(read_x.is_spanning_l) aligner.alignEndsFree(read_x.seq, 0, 0, read_y.seq, 0, length_diff);
			else if(read_x.is_spanning_r) aligner.alignEndsFree(read_x.seq, 0, 0, read_y.seq, length_diff, 0);
			else aligner.alignEndsFree(read_x.seq, 0, 0, read_y.seq, length_diff/2, length_diff/2);
			return (aligner.getAlignmentScore() / (double)read_x.seq.size());
		}
		else{
			if(read_x.is_spanning_l) aligner.alignEndsFree(read_y.seq, 0, length_diff, read_x.seq, 0, 0);
			else if(read_x.is_spanning_r) aligner.alignEndsFree(read_y.seq, length_diff, 0,read_x.seq, 0, 0);
			else aligner.alignEndsFree(read_y.seq, length_diff/2, length_diff/2, read_x.seq, 0, 0);
			return (aligner.getAlignmentScore() / (double)read_x.seq.size());
		}
	}
	else return -1.0;
}

double get_dist_anreads(const bool& ignore_haps, wfa::WFAligner& aligner, ANREAD& read_x, ANREAD& read_y)
{
	//when haplotags should be ignored
	if(ignore_haps) return align_anreads(aligner, read_x, read_y);
	else{
		//both reads are haplotagged
		if(read_x.hpt.is_defined() && read_y.hpt.is_defined()){
			if(read_x.hpt == read_y.hpt) return 0; else return 1.0;
		}
		//ambiguous, set them apart
		else return 1.0;
	}
}

void fill_dist_matrix(const bool& ignore_haps, wfa::WFAligner& aligner, std::vector<ANREAD>& reads, std::vector<int>& indeces, DistMatrix& distmatrix)
{
	for(uint32_t i = 0; i < indeces.size(); ++i){
		for(uint32_t j = i + 1; j < indeces.size(); ++j){
			distmatrix.set_dist(i, j, get_dist_anreads(ignore_haps, aligner, reads[indeces[i]], reads[indeces[j]]));
		}
	}
}

void invalid_reassignment(const bool& ignore_haps, const double& min_sim, const double& max_error, const int& total_alleles, std::vector<ANREAD>& reads, wfa::WFAligner& aligner, std::vector<int>& labels)
{
	//iterate through each read
	for(int i = 0; i < (int)labels.size(); ++i){
		//process unassigned read; note: assumes unlabeled reads are non-spanning or not haplotagged
		if(labels[i] < 0){
			//std::cout << "startin" << std::endl;
			//vector of maximum sim observed per allele albel
			std::vector<double> max_sim(total_alleles, 0.0);
			auto& read_i = reads[i];
			//iterate through each target read
			for(int j = 0; j < (int)labels.size(); ++j){
				auto& read_j = reads[j];
				//target read is assigned; note: assumes target is spanning
				if(i != j && labels[j] >= 0 && read_j.is_spanning()) {
					//std::cout << "here" << std::endl;
					//std::cout << aligner.getAlignmentStatus() << std::endl;
					//std::cout << read_i.name << '\t' << read_j.name << '\n';
					double dist = get_dist_anreads(true, aligner, read_i, read_j);
					//std::cout << dist << '\n';
					if(dist < 0) {
						std::cerr << "ERROR: unexpected distance for the following alignment:\n";
						std::cerr << read_i.name << '\t' << read_i.is_spanning() << '\n';
						std::cerr << read_i.seq << '\n';
						std::cerr << read_j.name << '\t' << read_j.is_spanning() << '\n';
						std::cerr << read_j.seq << std::endl;
						exit(1);
					}
					double sim = 1 - dist;
					if(sim > max_sim[labels[j]]) max_sim[labels[j]] = sim;
				}
			}
			int max_sim_label = 0;
			for(int j = 1; j < total_alleles; ++j) if(max_sim[j] > max_sim[max_sim_label]) max_sim_label = j;

			int same_max_sim = 0;
			for(const auto& s : max_sim) if(s == max_sim[max_sim_label]) ++same_max_sim;
			if(same_max_sim == 1){
				if(max_sim[max_sim_label] >= min_sim){
					double min_diff = 1.0;
					for(int j = 0; j < total_alleles; ++j) {
						if(max_sim_label != j) {
							double diff = max_sim[max_sim_label] - max_sim[j];
							if(diff < min_diff) min_diff =  diff;
						}
					}
					if(min_diff >= max_error) labels[i] = max_sim_label;
				}
			}
		}
	}
}

double compute_se(const std::vector<double>& values)
{
	if(values.empty()) return -1.0;
	else{
		double u = 0.0, n = 0.0;
		for(const auto& v : values) u += v;
		u /= values.size();
		for(const auto& v : values) n += std::pow(v - u, 2.0);
		return std::sqrt(n/(values.size() - 1)) / std::sqrt(values.size());
	}
	
}

void rapid_consensus(const bool& ignore_haps, std::vector<ANREAD>& reads, std::vector<int>& labels, std::vector<int>& valid_indeces, int& total_alleles, DistMatrix& valid_distmatrix, wfa::WFAligner& aligner, std::vector<ANALLELE>& alleles)
{
	if(valid_indeces.empty()){
		std::cerr << "ERROR: empty vector of valid read-indeces" << std::endl;
		exit(1);
	}
	for(int label = 0; label < total_alleles; ++label){
		/** get indeces of valid reads corresponding to allele label **/
		//indexed to reads
		std::vector<uint32_t> label_indeces_valid_reads;
		//indexed to *valid_indeces*
		std::vector<uint32_t> label_indeces_valid_indeces;
		for(uint32_t i = 0; i < valid_indeces.size(); ++i) {
			if(label == labels[valid_indeces[i]]) {
				label_indeces_valid_reads.emplace_back(valid_indeces[i]);
				label_indeces_valid_indeces.emplace_back(i);
			}
		}
		if(label_indeces_valid_reads.empty()){
			std::cerr << "ERROR: empty vector of valid read-indeces for allele cluster " << label << std::endl;
			exit(1);
		}
		//find representative
		int rep_index_valid_indeces = valid_distmatrix.get_medoid(label_indeces_valid_indeces);
		//std::cout <<"rep_i: " << label << '\t' << rep_index_valid_indeces << std::endl;
		int rep = valid_indeces[rep_index_valid_indeces];
		//std::cout <<"rep: " << label << '\t' << rep << '\t' << rep_index_valid_indeces << std::endl;
		//indexed to reads; excludes rep
		std::vector<int> label_indeces_all_reads;
		for(int i = 0; i < (int)reads.size(); ++i) if(i != rep && labels[i] == label) label_indeces_all_reads.emplace_back(i);
		
		/** construct allele object **/
		auto& local_allele = alleles[label];
		local_allele.tcov = reads.size();
		local_allele.acov = label_indeces_all_reads.size() + 1;
		local_allele.scov = label_indeces_valid_reads.size();
		/** calcualte SE **/
		if(label_indeces_valid_indeces.size() == 1) local_allele.se = 0;
		else if(label_indeces_valid_indeces.size() == 2) local_allele.se = valid_distmatrix.get_dist(label_indeces_valid_indeces[0], label_indeces_valid_indeces[1]);
		else{
			std::vector<double> valid_dists;
			for(const auto& i : label_indeces_valid_indeces) if((int)i != rep_index_valid_indeces) valid_dists.emplace_back(valid_distmatrix.get_dist(i, rep_index_valid_indeces));
			local_allele.se = compute_se(valid_dists);
		}

		int ps = -1, hp = -1;
		bool conflicting = false;
		if(!ignore_haps){
			for(const auto& i : label_indeces_valid_reads) {
				//std::cout << reads[i].name << '\t' << reads[i].seq << '\n';
				if(ps < 0) ps = reads[i].hpt.ps;
				else if(ps != reads[i].hpt.ps) conflicting = true;
				if(hp < 0) hp = reads[i].hpt.hp;
				else if(hp != reads[i].hpt.hp) conflicting = true;
			}
		}

		if(conflicting) {
			std::cerr << "ERROR: conflicting haplotag information:\n";
			for(const auto& i : label_indeces_valid_reads) std::cerr << reads[i].name << '\t' << reads[i].hpt.ps << '\t' << reads[i].hpt.hp << '\n';
			std::cerr << std::endl;
			exit(1);
		}

		auto& rep_read = reads[rep];
		if(!ignore_haps) local_allele.hpt = rep_read.hpt;

		if(label_indeces_all_reads.size() + 1 <= 2) local_allele.seq = reads[label_indeces_valid_reads.front()].seq;
		else{
			PPOA poa;
			poa.init(rep_read.seq);
			for(const auto& i : label_indeces_all_reads){
				auto& read = reads[i];
				//std::cout << read.name << ' ' << read.is_spanning_l << ' ' << read.is_spanning_r << std::endl;
				int length_diff = rep_read.seq.size() - read.seq.size();
				if(read.is_spanning() || length_diff < 0) {
					if(length_diff >= 0) aligner.alignEnd2End(rep_read.seq, read.seq);
					else{
						if(read.is_spanning_l) aligner.alignEndsFree(rep_read.seq, 0, 0, read.seq, 0, -length_diff);
						else if(read.is_spanning_r) aligner.alignEndsFree(rep_read.seq, 0, 0, read.seq, -length_diff, 0);
					}
				}
				else{
					//std::cout << read.seq << std::endl;
					if(read.is_spanning_l) aligner.alignEndsFree(rep_read.seq, 0, length_diff, read.seq, 0, 0);
					else if(read.is_spanning_r) aligner.alignEndsFree(rep_read.seq, length_diff, 0, read.seq, 0, 0);
					else aligner.alignEndsFree(rep_read.seq, length_diff/2, length_diff/2, read.seq, 0, 0);
				}
				std::string cigar = aligner.getAlignmentCigar();
				//std::cout << cigar << std::endl;
				poa.insert_alignment(read.seq, cigar, read.is_spanning_l, read.is_spanning_r);
			}
			//std::cout << "adjusting" << std::endl;
			float c =(label_indeces_all_reads.size()+1)*0.4;
			float t = 0.3;
			if(label_indeces_all_reads.size() + 1 < 4) c = 1.0;
			poa.adjust_weights(c,t);
			//poa.printDOT();
			//std::cout << std::endl;
			poa.consensus(local_allele.seq);
			if(local_allele.seq.empty()) local_allele.seq = "N";
			//std::cout << "consensus: " << local_allele.seq << std::endl;
		}

	}

}

