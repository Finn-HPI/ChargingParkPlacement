#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <memory>

bool contains(std::string& delimiters, char c) {
	bool result = false;
	for (char ch : delimiters) result |= ch == c;
	return result;
}

std::vector<std::string> split(std::string& line, std::string delimiters, bool empty_allowed) {
	std::vector<std::string> tokens;
	size_t start = 0, end = 0;
	while (start < line.size()) {
		while (end < line.size() && !contains(delimiters, line[end])) ++end;
		if (start < end || empty_allowed)
			tokens.push_back(line.substr(start, end - start));
		start = ++end;
	}
	return tokens;
}