#pragma once

#include <vector>
#include <sstream>

std::vector<double> zeros(const size_t rowCount)
{
    return std::vector<double>(rowCount, 0);
}
std::vector<std::vector<double>> zeros(const size_t rowCount, const size_t columnCount)
{
    return std::vector<std::vector<double>>(rowCount, std::vector<double>(columnCount, 0));
}
std::string toString(const std::vector<std::vector<double>>& matrix)
{
    std::ostringstream sstream;

    for(size_t i = 0; i < matrix.size(); i++)
    {
        if(i > 0)
        {
            sstream << '\n';
        }

        for(size_t j = 0; j < matrix[0].size(); j++)
        {
            if(j > 0)
            {
                sstream << ',';
            }

            sstream << matrix[i][j];
        }
    }

    return sstream.str();
}