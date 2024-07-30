#include "bindings/cpp/WFAligner.hpp"
#include "andistmat.hpp"
#include "anseqs.hpp"
#include <vector>
#include <string>

class ClusteringStatus{
	public:
		int ic;
		int fc;
		std::vector<int> labels;
		ClusteringStatus();
		void set_global_label(int);
};

class Genotype{
	public:
		int gt;
		int gt_l;
		int gt_k;
		double hsd;
		Genotype();
};

class DecisionBound{
	public:
		double dist0;
		double dist1;
		double cut0;
		DecisionBound(double, double, double);
};

DecisionBound otter_find_clustering_dist(
	const int&,
	const double&,
	const double&,
	DistMatrix&
);

void otter_hclust(
	const bool& ignore_haps, 
	const int& max_alleles, 
	const double& bandwidth, 
	const double& max_tolerable_diff, 
	const double& min_cov_fraction, 
	const int& min_cov_fraction2_l, 
	const double& min_cov_fraction2_f,
	const std::vector<int>& indeces, 
	DistMatrix& distmatrix, wfa::WFAligner& aligner, 
	std::vector<ANREAD>& reads, 
	ClusteringStatus& clustering
);

double length_dist(
	const uint32_t& x, 
	const uint32_t& y
);

void generate_kusage(
	const uint32_t& k, 
	const KmerEncoding& encoding, 
	const std::vector<ANALLELE>& alleles, 
	const std::vector<int>& indeces, 
	std::vector<KUSAGE>& kusages
);

int anallele_cluster(
	const double& max_error_l, 
	const double& max_error_c,
	const std::vector<ANALLELE>& alleles, 
	std::vector<Genotype>& genotypes,
	 std::vector<int>& gt_reps
);
