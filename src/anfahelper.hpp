#ifndef ANFAHELPER_HPP
#define ANFAHELPER_HPP

#include "faidx.h"
#include <string>

class FaidxInstance {
	public:
		faidx_t* faidx_ptr;
		void init(const std::string&);
		int fetch(const std::string&, int, int, std::string&) const;
		void destroy();
};

#endif