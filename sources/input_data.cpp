#include "input_data.h"

#include <regex>
#include <stdexcept>

Input_data::Input_data(std::ifstream &input_file) {
    if (!input_file.is_open())
        throw std::runtime_error("Input file is not open!");

    std::string line;
    while (std::getline(input_file, line)) {
        if (line.empty())
            continue;

        std::regex reg("\\s+");

        std::sregex_token_iterator beg(line.begin(), line.end(), reg, -1);
        std::sregex_token_iterator end;

        std::string key = *beg;
        std::vector<std::string> vec(++beg, end);
        keys.insert(std::make_pair(key, vec));
    }
}

std::string Input_data::first(const std::string &key) const {
    return keys.at(key).at(0);
}

std::string Input_data::second(const std::string &key) const {
    return keys.at(key).at(1);
}

std::string Input_data::operator()(const std::string &key, const int &index) const {
    return keys.at(key).at(index);
}

int Input_data::size(const std::string &key) const {
    return keys.at(key).size();
}