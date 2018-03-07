#pragma once

#include <vector>
#include <sstream>

std::vector<double> zeros(const size_t rowCount);
std::vector<std::vector<double>> zeros(const size_t rowCount, const size_t columnCount);
std::string toString(const std::vector<std::vector<double>>& matrix);