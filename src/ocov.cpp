#include "antimestamp.hpp"
#include "anbamfilehelper.hpp"
#include "anbed.hpp"
#include "BS_thread_pool.hpp"
#include "otter_opts.hpp"
#include <mutex>
#include <iostream>
#include <string>
#include <vector>

void get_cov(const OtterOpts& params, const BED& bed, const BamInstance& bam_inst, int& local_coverage)
{
	hts_itr_t *iter = sam_itr_querys(bam_inst.idx, bam_inst.header, bed.toBEDstring().c_str());
	if(iter == nullptr) std::cerr << "(" << antimestamp() << "): WARNING: query failed at region " <<  bed.toBEDstring() << std::endl;
	else{
		while(sam_itr_next(bam_inst.fp, iter, bam_inst.read) > 0){
			if(bam_inst.read->core.qual >= params.mapq && (params.nonprimary || !(bam_inst.read->core.flag & BAM_FSECONDARY || bam_inst.read->core.flag & BAM_FSUPPLEMENTARY))){
				++local_coverage;
			}
		}
	}
	hts_itr_destroy(iter);
}

void ocov(const OtterOpts& params, const std::string& bed, const std::string& bam)
{
	BS::thread_pool pool(params.threads);
 	
 	std::vector<BED> bed_regions;
 	parse_bed_file(bed, bed_regions);

 	std::cerr << '(' << antimestamp() << "): Processing " << bam << '\n';
	std::mutex std_out_mtx;


	bool is_read_group = params.read_group.size() > 1;
	if(is_read_group) std::cout << "#Sample\tRegion\tCoverage\n";
	else std::cout << "#Region\tCoverage\n";

	pool.parallelize_loop(0, bed_regions.size(), 
		[&pool, &std_out_mtx, &params, &bam, &bed_regions, &is_read_group](const int a, const int b){
			BamInstance bam_inst;
			bam_inst.init(bam, true);
			for(int i = a; i < b; ++i) {
				const BED& local_bed = bed_regions[i];
				BED mod_bed = local_bed;
				mod_bed.start -= (uint32_t)params.offset;
				mod_bed.end += (uint32_t)params.offset;
				int local_coverage = 0;
				get_cov(params, mod_bed, bam_inst, local_coverage);
				std_out_mtx.lock();
				if(is_read_group) std::cout << params.read_group << '\t';
				std::cout << local_bed.toBEDstring() << '\t' << local_coverage << '\n';
				std_out_mtx.unlock();
			}
			bam_inst.destroy();
		}
	).wait();

}