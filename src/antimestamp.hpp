#ifndef ANTIMESTAMP_HPP
#define ANTIMESTAMP_HPP
#include <string>
#include <ctime>

/**
 * Function to generate current time-stamp as a string
 * 
 * @return std::string
 */
inline std::string antimestamp()
{
    time_t rawtime;
    time (&rawtime);
    std::string t;
    t = (ctime (&rawtime));
    t.erase(t.end() - 1);
    return t;
}

#endif