#ifndef FASTA_HELPER_HPP
#define FASTA_HELPER_HPP

#include <string>
#include <htslib/faidx.h>

class FaidxInstance {
	public:
		faidx_t* faidx_ptr;
		void init(const std::string&);
		int fetch(const std::string&, int, int, std::string&) const;
		void destroy();
};

#endif