#include "parse_bam_alignments.hpp"
#include "dist_matrix.hpp"
#include "fasta_helper.hpp"
#include <string>
#include "bindings/cpp/WFAligner.hpp"

double otter_seq_dist(
	wfa::WFAligner& aligner, 
	std::string& x, 
	ParsingStatus& x_status, 
	std::string& y, 
	ParsingStatus& y_status
);

void otter_pairwise_dist(
	const std::vector<int>&,
	AlignmentBlock&,
	wfa::WFAligner&,
	DistMatrix&
);

std::pair<double,double> otter_find_clustering_dist(
	const double&,
	const std::vector<int>&,
	AlignmentBlock&, 
	wfa::WFAligner&,
	DistMatrix&
);

void otter_hclust(
	const int& max_alleles,
	const double& bandwidth,
	const double&, 
	const std::vector<int>&, 
	DistMatrix& distmatrix,
	wfa::WFAligner& aligner, 
	AlignmentBlock& sequences,
	std::vector<int>& cluster_labels
);

void otter_realignment2(
	const std::string& chr, 
	const int& flank, 
	const double& min_sim, 
	std::string& ref_left, 
	std::string&ref_right, 
	AlignmentBlock& sequences, 
	wfa::WFAligner& aligner
);

void otter_realignment(
	const std::string& chr, 
	const int& start, 
	const int& end, 
	const int& flank, 
	const double& min_sim,
	const FaidxInstance&,
	AlignmentBlock& sequences, 
	wfa::WFAligner& aligner
);

void otter_nonspanning_assigment(
	const double&,
	const double&,
	AlignmentBlock& sequences, 
	wfa::WFAligner& aligner, 
	std::vector<int>& labels
);

void otter_rapid_consensus(
	const std::vector<int>& spannable_indeces, 
	const std::vector<int>& labels, const DistMatrix& distmarix, 
	wfa::WFAligner& aligner,
	AlignmentBlock& sequences, 
	std::vector<std::string>& consensus,
	std::vector<std::vector<double>>& ses
);

