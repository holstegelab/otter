#ifndef ANGZIPITER_HPP
#define ANGZIPITER_HPP

#include <zlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static const unsigned GZIPITER_BUFLEN = 1000000;

class GZIPiter {
	public:
		GZIPiter(const std::string&);
		bool hasNext();
		void next(std::string&);
		void close();

	private:
		void error(const char* const);
		void scanNext();
		gzFile input_file;
		std::string line;
		int err;
		//set buffer
    	char buf[GZIPITER_BUFLEN];
    	//set pointer to buffer offset location
    	char* offset = buf;
    	//cur pointer to buf
		char* cur = buf;
    	//ending pointer to buf
		char* end = buf;
		char* eol;
};

inline GZIPiter::GZIPiter(const std::string& file): input_file(gzopen(file.c_str(), "rb")){}

inline void GZIPiter::scanNext()
{
	//double check that there is enough rooom
	int bytes_read = sizeof(buf)-(offset-buf);
	if (bytes_read == 0) error("Buffer to small for input line lengths");
	bytes_read = gzread(input_file, offset, bytes_read);
	//problem reading bytes
	if(bytes_read < 0) error(gzerror(input_file, &err));
	//std::cout << "bytes reads: " << bytes_read << '\n';
	//std::cout << "offset: " << (offset - buf) << '\n';
	//starting pointer to buf
	cur = buf;
	//ending pointer to buf
	end = offset+bytes_read;
	offset = buf;
}

inline bool GZIPiter::hasNext()
{
	//std::cout <<"pointer offset: " << (end - cur) << '\n';
	if(end > cur) return true;
	else{
		scanNext();
		//std::cout <<"pointer offset: " << (end - cur) << '\n';
		return end > cur;
	}
}

inline void GZIPiter::next(std::string& output)
{
	//std::cout << "eol pointer: " << (std::find(cur, end, '\n') < end) << "\n";
	if(cur < end && (eol = std::find(cur, end, '\n')) < end){
		//std::cout << "setting\n"
		output = std::string(cur,eol);
		cur = eol + 1;
	}
	else {
		//std::cout << "resetting\n";
		//copy trailing/remaining data to buffer
    	offset = std::copy(cur, end, buf);
    	cur = buf;
    	end = buf;
	}
}

inline void GZIPiter::error(const char* const msg)
{
    std::cerr << msg << "\n";
    exit(255);
}

inline void GZIPiter::close()
{
	if (gzclose(input_file) != Z_OK) error("failed gzclose");
}

inline void parse_line(const std::string& line, const char& delim, std::vector<std::string>& columns)
{
	std::string value;
  	std::istringstream stringstream(line);
  	while(std::getline(stringstream, value, delim)) columns.emplace_back(value);
}

#endif