#ifndef SOURCES_INPUT_DATA_H
#define SOURCES_INPUT_DATA_H

#include <fstream>
#include <map>
#include <string>
#include <vector>

class Input_data {
   public:
    Input_data(std::ifstream &input_file);

    std::map<std::string, std::vector<std::string>> keys;
    
    std::string first(const std::string &key) const;
    std::string second(const std::string &key) const;
    std::string operator()(const std::string &key, const int &index) const;
    int size(const std::string &key) const;
};

#endif /* SOURCES_INPUT_DATA_H_ */
