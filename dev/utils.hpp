#pragma once

namespace VUtils
{
	/** @brief
	 * \attribution Inspired by https://stackoverflow.com/a/16575025 .*/
	bool is_num(std::string line, int base = 10)
	{
	    char* p;
	    strtol(line.c_str(), &p, base);
	    return *p == 0;
	}
};