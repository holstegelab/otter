#include "anbamfilehelper.hpp"
#include "antimestamp.hpp"
#include "BS_thread_pool.hpp"
#include <string>
#include <map>
#include <vector>
#include <thread>
#include <iostream>
#include <mutex>

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2))

//BamInstance::ansbam(){};

void BamInstance::init(const std::string& file, bool _index)
{
	index = _index;
	fp = hts_open(file.c_str(), "r");
	header = sam_hdr_read(fp);
	//load index
	if(index) idx = sam_index_load(fp, file.c_str());
	//init read
	read = bam_init1();
};

void BamInstance::destroy()
{
	sam_close(fp);
	fp = nullptr;
	if(index) hts_idx_destroy(idx);
	idx = nullptr;
	bam_hdr_destroy(header);
	header = nullptr;
	bam_destroy1(read);
	read = nullptr;
};

void allocate_bam_threads(const int threads, const std::string& bam, BS::thread_pool& pool, std::map<std::thread::id, int>& thread2index, std::vector<BamInstance>& bam_insts)
{
	//initiate bam file openers for each thread
 	bam_insts.resize(threads, BamInstance());
 	//thread index lock
 	std::mutex thread_int_mutex;
 	int thread_int = 0;
 	//run all threads and log thread id -> index
 	pool.parallelize_loop(0, threads,[&thread_int_mutex, &thread_int, &thread2index](const int a, const int b){
 		//sleep or else not all specified threads will be logged
 		std::this_thread::sleep_for(std::chrono::milliseconds(500));
 		thread_int_mutex.lock();
 		thread2index.insert({ std::this_thread::get_id(), thread_int });
 		++thread_int;
 		thread_int_mutex.unlock();
 	}).wait();
 	
 	for(uint32_t i = 0; i < thread2index.size(); ++i) bam_insts[i].init(bam, true);

 	std::cerr << "(" << antimestamp() << "): Allocated " << thread2index.size() << " bam-file instances" << std::endl;
}

void destroy_bam_threads(std::map<std::thread::id, int>& thread2index, std::vector<BamInstance>& bam_insts)
{
	for(uint32_t i = 0; i < thread2index.size(); ++i) bam_insts[i].destroy();
	std::cerr << "(" << antimestamp() << "): Destroyed " << thread2index.size() << " bam-file instances" << std::endl;
}