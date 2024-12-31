//
// Created by nloyfer on 5/27/19.
//

#include "homog.h"


/**
 *  Main
 *
 */

/*
 * Input Arguments Parsering class
 */
class InputParser {
public:
    InputParser(int &argc, char **argv) {
        for (int i = 1; i < argc; ++i)
            this->tokens.emplace_back(std::string(argv[i]));
    }

    const std::string &getCmdOption(const std::string &option) const {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
            return *itr;
        }
        static const std::string empty_string;
        return empty_string;
    }

    bool cmdOptionExists(const std::string &option) const {
        return std::find(this->tokens.begin(), this->tokens.end(), option)
               != this->tokens.end();
    }

private:
    std::vector <std::string> tokens;
};


std::string get_param_str(InputParser &input, std::string name, std::string defval) {
    std::string param_str = input.getCmdOption(name);
    if (!param_str.empty()) {
        return param_str;
    }
    return defval;
}


void print_help() {
    std::cerr << "" << std::endl;
}

std::vector<float> parse_range(std::string &range_str) {
    std::vector<float> vect;

    std::stringstream ss(range_str);

    for (float i; ss >> i;) {
        vect.push_back(i);
        if (ss.peek() == ',')
            ss.ignore();
    }
    float tmp = -1;
    for (std::size_t i = 0; i < vect.size(); i++) {
        if (vect[i] <= tmp) {
            std::cout << "Invalid range - non monotonic: " << range_str << std::endl;
            throw std::invalid_argument("Invalid range");
        }
        if ((vect[i] < 0) || (vect[i] > 1)) {
            std::cout << "Invalid range - must be in [0,1]: " << vect[i] << std::endl;
            throw std::invalid_argument("Invalid range");
        }
        tmp = vect[i];
//        std::cout << vect[i] << std::endl;
    }
    if ((vect[0] > 0) || (vect[vect.size() - 1] < 1)) {
        std::cout << "Invalid range - must start with 0 and end with 1."  << std::endl;
        throw std::invalid_argument("Invalid range");
    }
    return vect;
}

int main(int argc, char *argv[]) {

    if ((argc < 3)) {
        std::cerr << "Usage: EXEC PAT_PATH OUTPUT_PREFIX [-r RANGE] [-d] [-b BLOCKS_PATH] [-l MIN_LEN] [--gzip]\n"
                  << "-l       Minimal sites in read to consider. Default is l=5.\n"
                  << "-r       Ranges of methylation average, in [0,1]. For example: 0,0.2001,0.8,1.\n"
                  << "-b       Blocks file. No header. Columns are:\n"
                  << "         <chrall_blocks.sorted.tsv, start, end, startCpG, endCpG>. Default is s150.\n"
                  << "         File may be gzipped. Note sometimes bgzip causes problems.\n"
                  << "-d       Debug mode. Only use the first 2,500 blocks in the blocks file."
                  << std::endl;
        return -1;
    }

    InputParser input(argc, argv);

    bool help = input.cmdOptionExists("-h");
    if (help) {
        print_help();
        return 1;
    }
//    std::string pat_path = "/cs/cbio/jon/segmentation_files/cortex_small.TF.pat";//
//    std::string output_path = "";//
//    std::string blocks_path = "/cs/cbio/jon/segmentation_files/cortex_small.TF.homog";//
//    std::string range_str = "0,.201,.8,1";
    std::string pat_path = std::string(argv[1]);//"/cs/cbio/jon/segmentation_files/cortex_small.TF.pat";//
    std::string output_path = std::string(argv[2]);//"";//
    std::string blocks_path = get_param_str(input, "-b", DEFAULT_BLOCKS);//"/cs/cbio/jon/segmentation_files/cortex_small.TF.homog";//
    std::string range_str = get_param_str(input, "-r", "0,.201,.8,1");//"0,.201,.8,1";
    if (!range_str.size()) {
        std::cout << "You must provide a range" << std::endl;
        return 1;
    }
    int min_cpgs = std::stoi(get_param_str(input, "-l", DEFAULT_LEN));//3;//

    bool debug = input.cmdOptionExists("-d");//false;//

    try {
        std::vector<float> range = parse_range(range_str);
        Homog(output_path, blocks_path, range, min_cpgs, debug).parse(pat_path);
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}

