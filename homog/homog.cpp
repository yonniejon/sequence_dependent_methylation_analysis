

#include <queue>
#include "homog.h"




std::vector <std::string> line2tokens(std::string &line) {
    /**
     * Break string line to tokens, return it as a vector of strings
     */
    std::vector <std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while (getline(lineStream, cell, '\t'))
        result.push_back(cell);
    return result;
}

void print_vec(std::vector <std::string> &vec) {
    /** print a vector to stderr, tab separated */
    for (auto &j: vec)
        std::cerr << j << "\t";
    std::cerr << std::endl;
}

auto end_compare = [](Block a, Block b) { return a.end > b.end; };

//Block pop_heap(std::vector<Block> a_vec){
//    std::pop_heap (a_vec.begin(),a_vec.end(), end_compare);
//    Block a = a_vec.back();
//    a_vec.pop_back();
//    return a;
//}
//void push_heap(std::vector<Block> a_vec, Block a_block){
//    a_vec.push_back(a_block);
//    std::push_heap(a_vec.begin(), a_vec.end(), end_compare);
//}

bool intersects_block(int read_start, int read_end, int block_start, int block_end){
    bool is_inside = (read_start >= block_start && read_start < block_end) || (read_end >= block_start && read_end < block_end);
    bool read_contains = read_start <= block_start && read_end >= block_end;
    return is_inside || read_contains;
}


int32_t *Homog::init_array(int len) {
    int *arr = new int32_t[nr_blocks * len];
    std::fill_n(arr, nr_blocks * len, 0);
    return arr;
}

// Constructor
Homog::Homog(std::string in_output_prefix, std::string in_blocks_path, std::vector<float> in_range,
             int in_min_cpgs, bool deb) {
    min_cpgs = in_min_cpgs;
    output_prefix = in_output_prefix;
    blocks_path = in_blocks_path;
    range = in_range;
    debug = deb;

    nr_bins = range.size() - 1;

    // load blocks file
    int r = read_blocks();
    nr_blocks = borders.size();

    // Init arrays to zeros
    counts = init_array(nr_bins);
}

Homog::~Homog() {
    //delete[] counts;
}


int Homog::blocks_helper(std::istream &instream) {
    //Iterate lines
    std::vector <std::string> tokens;
    std::string line;
    int cur_start = 0, cur_end = 0;
    int bi = 0;
    while (std::getline(instream, line)) {
//        std::cerr << ++abc << std::endl;
        

        // skip empty lines and comments
        if (line.empty() || (!(line.rfind("#", 0)))) { continue; }

        tokens = line2tokens(line);
        if (tokens.size() < 5) {
            std::cerr << "Invalid blocks file format. ";
            std::cerr << "Should be: chr start end startCpG endCpG\n";
            abort();
        }

        // skip header, if exists
        if ((bi == 0) && (tokens[0] == "chr")) {
            continue;
        }

        cur_start = std::stoi(tokens[3]);
        cur_end = std::stoi(tokens[4]);

        // If block is invalid, abort:
        if ((cur_end <= cur_start)) {
            std::cerr << "Invalid block: " << cur_start << "\t" << cur_end << std::endl;
            throw std::invalid_argument("Invalid block: endCpG <= startCpG");
        } else if (cur_start < 1) {
            throw std::invalid_argument("Invalid block: startCpG < 1");
        }

        // skip borders with <min_cpgs CpGs
//        if ((cur_end - cur_start < min_cpgs)) {
//            continue;
//        }

        // if block is duplicated, continue
        if ((!borders.empty()) &&    // Can't be dup if it's first
                (borders.back().start == cur_start) &&
            (borders.back().end == cur_end)) {
            // only update the count and continue:
            borders_counts[bi - 1]++;
            continue;
        }
        
        // block isn't dup:
        borders_counts.push_back(1);
        borders.push_back(Block(cur_start, cur_end, bi));
        coords.push_back(tokens[0] + "\t" + tokens[1] + "\t" + tokens[2]);
        

        // make sure there's no overlap or non-monotonic blocks
        if ( (bi > 0) && (borders[bi].start < borders[bi - 1].start) ) {
            std::cerr << "Invalid block start: " << cur_start << std::endl;
            std::cerr << "Make sure blocks are sorted by CpG-Index and monotonic (blockEnd > blockStart).\n";
            throw std::invalid_argument("Invalid blocks");
        }

        if (debug && (bi >= 2500)) {
            break;
        }
        bi++;
    }
    return 0;
}

int Homog::read_blocks() {
    /**
     * Load blocks gzipped file into vector<int> borders_starts, borders_ends.
     */


    std::cout << "loading blocks..." << std::endl;

    if (hasEnding(blocks_path, ".gz")) {
        // Open the gzipped file:
        std::ifstream file(blocks_path, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf <boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(file);
        std::istream instream(&inbuf);
        blocks_helper(instream);
        //Cleanup
        file.close();


    } else {
        std::ifstream instream(blocks_path);
        blocks_helper(instream);
        instream.close();
    }

    if (borders.size() == 0) {
        std::cerr << "Error while loading blocks. 0 blocks found\n";
        throw std::invalid_argument("");
     }
    std::cerr << "nr_blocks: " << borders.size() << ", last block: " ;
    std::cerr << borders[borders.size() - 1].start << "-" << borders[borders.size() - 1].end << std::endl;
    return 0;
}


void Homog::dump(int32_t *data, int width, std::string out_path) {
    /**
     */
    std::ofstream bofs;
    bofs.open(out_path);
    if (!(bofs.is_open())) {
        std::cerr << "could not open output file " << out_path << std::endl;
        return;
    }

    std::cerr << "dumping to " << out_path << "... " << std::endl;
    for (int i = 0; i < nr_blocks; i++) {
        for (int d = 0; d < borders_counts[i]; d++) {
            bofs << coords[i] << SEP << borders[i].start << SEP << borders[i].end << SEP;
            for (int j = 0; j < width; j++) {
                bofs << data[i * width + j];
                if (j < width - 1)
                    bofs << SEP;
            }
            bofs << std::endl;
        }
    }
    bofs.close();

    // gzip:
    system(("bgzip -f " + out_path + " && tabix -p bed " + out_path + ".gz").c_str());
}


void Homog::update_m2(int block_ind, std::string pat, int count) {

    int nrC = 0;
    int nrT = 0;
    for (int i = 0; i < pat.length(); i++) {
        if (pat[i] == 'C')
            nrC++;
        else if (pat[i] == 'T')
            nrT++;
    }

    if (nrC + nrT < min_cpgs) {
//	    printf("nrC + nrT <  min_cpgs. %d + %d < %d\n", nrC, nrT, min_cpgs);
        return;
    }

    float meth = (float) nrC / (float) (nrC + nrT);
    if (meth < range[0])
        return;

    int bin_ind = 0;
    for (bin_ind = 0; bin_ind < range.size() - 1; bin_ind++) {
        if ((meth >= range[bin_ind]) && (meth < range[bin_ind + 1])) {
//	    printf("smaller than %f. returning bin_ind %d\n", range[bin_ind + 1], bin_ind);
            break;
        }
    }
    if (bin_ind == range.size() - 1) {
        bin_ind--;
    }
//    std::cout << "block:" << block_ind << ", pat: " <<  pat << ", meth: " << meth << ", bin_ind:" <<  bin_ind << std::endl;
    counts[block_ind * nr_bins + bin_ind] += count;
}

void Homog::update_m2_blocks( std::string pat, int count, int block_id) {

    int nrC = 0;
    int nrT = 0;
    for (int i = 0; i < pat.length(); i++) {
        if (pat[i] == 'C')
            nrC++;
        else if (pat[i] == 'T')
            nrT++;
    }

    if (nrC + nrT < min_cpgs) {
//	    printf("nrC + nrT <  min_cpgs. %d + %d < %d\n", nrC, nrT, min_cpgs);
        return;
    }

    float meth = (float) nrC / (float) (nrC + nrT);
    if (meth < range[0])
        return;

    int bin_ind = 0;
    for (bin_ind = 0; bin_ind < range.size() - 1; bin_ind++) {
        if ((meth >= range[bin_ind]) && (meth < range[bin_ind + 1])) {
//	    printf("smaller than %f. returning bin_ind %d\n", range[bin_ind + 1], bin_ind);
            break;
        }
    }
    if (bin_ind == range.size() - 1) {
        bin_ind--;
    }
    counts[block_id * nr_bins + bin_ind] += count;
//    for (Block block : cur_block_queue){
//
//    }
}

void Homog::update_block(int block_ind, std::string pat, int32_t count) {

//    std::cerr << block_ind << ") updating: " << pat << std::endl;

    int len = pat.length();
    // skip reads with less then min_cpgs sites:
    if (len < min_cpgs) {
//	printf("%d < %d\n", len, min_cpgs);
        return;
    }

    update_m2(block_ind, pat, count);
}
void Homog::update_blocks(std::string pat, int32_t count, int read_start, int read_end) {

//    std::cerr << block_ind << ") updating: " << pat << std::endl;

    int len = pat.length();
    // skip reads with less then min_cpgs sites:
    if (len < min_cpgs) {
//	printf("%d < %d\n", len, min_cpgs);
        return;
    }
    for (Block block : cur_block_queue){
        if (intersects_block(read_start, read_end, block.start, block.end)){
            if ((read_start >= block.start) && (read_start < block.end)) {
                int head_size = std::min(block.end - read_start, (int) pat.length());
//            pat.substr(head_size);
                std::string temp_pat = pat.substr(0, head_size);
                update_m2_blocks(temp_pat, (int32_t) count, block.id);
            } else if (read_end < block.start) {
                break;
            } else { // read_start < block.start
                int pat_len = 0;
                if (block.end <= read_end){
                    pat_len = block.end - block.start;
                } else {
                    pat_len = read_end - block.start + 1;
                }
//                int pat_len = std::min(block.end, read_end) - block.start + 1;
                std::string temp_pat = pat.substr(block.start - read_start, pat_len);
                update_m2_blocks(temp_pat, (int32_t) count, block.id);
            }
        } else {
            continue;
        }
    }
}

int Homog::proc_line(std::vector <std::string> tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", update counts array -
     * Increase in the corresponding cells.
     * Overall strategy: keep a queue of blocks in sorted by end CpG index.
     * At any given line of the pat file the queue contains all blocks which
     * can intersect the read.
     *
     */
    if (tokens.size() < 4) {
        throw std::invalid_argument("Invalid site in input file. too few columns");
    }

    int read_start = std::stoi(tokens[1]);
    std::string pattern = tokens[2];
    int count = std::stoi(tokens[3]);
    int read_end = read_start + (int) pattern.length() - 1;



//    print_vec(tokens);


    // read starts later than the last border - finish
    //                                          CTCTCT
    // |---|   |-----| |--|  ...  |-----| |-|
    if (read_start > borders[borders.size() - 1].end) {
        std::cerr << "Reads reached the end of blocks: " << read_start << " > "
                  << borders[borders.size() - 1].end << ". Breaking.\n";
        return 1;
    }

    // todo: dump block line on the fly
    // read ends before current block starts: skip read.
    //  CTCTCT
    //         |-----|
    while(cur_block_ind < borders.size() && read_end >= borders[cur_block_ind].start) {
        if (intersects_block(read_start, read_end, borders[cur_block_ind].start, borders[cur_block_ind].end)){
            Block block = Block(borders[cur_block_ind].start, borders[cur_block_ind].end, cur_block_ind);
            cur_block_queue.push_back(block);
            std::push_heap(cur_block_queue.begin(), cur_block_queue.end(), end_compare);
        }
        cur_block_ind++;
    }
    while (cur_block_queue.size() > 0 && cur_block_queue.front().end <= read_start){
        std::pop_heap (cur_block_queue.begin(),cur_block_queue.end(), end_compare);
        cur_block_queue.pop_back();
    }
    if (cur_block_queue.size() > 0){
        update_blocks(pattern, (int32_t) count, read_start, read_end);
        return 0;
    } else {
        return 0;
    }
}

void Homog::debug_print(int ind, std::vector <std::string> tokens) {
    std::cerr << ind << ") " << borders[ind].start << "-" << borders[ind].end << std::endl;
    print_vec(tokens);
}


bool hasEnding(std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void Homog::parse(std::string pat_path) {

    try {

        // Open the gzipped pat file:
//        std::ifstream file(pat_path, std::ios_base::in | std::ios_base::binary);
//        boost::iostreams::filtering_streambuf <boost::iostreams::input> inbuf;
//        inbuf.push(boost::iostreams::gzip_decompressor());
//        inbuf.push(file);
//        std::istream instream(&inbuf);
//        if (!(instream)) {
//            printf("not opened right\n");
//            return;
//        }


        int line_i = 0;
        std::ifstream infile(pat_path);

        for (std::string line_str; std::getline(std::cin, line_str);) {
//        for (std::string line_str; std::getline(infile, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines

            std::vector <std::string> tokens = line2tokens(line_str);
//            std::cerr << line_i << std::endl;
//	    print_vec(tokens); // TODO: del
            if (tokens.size() < 4) {
                throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
            } else if (!(tokens.empty())) {
                int r = proc_line(tokens);
                if (r == 1) {
                    std::cerr << "breaking" << std::endl;
                    break;
                }
            } else {
                std::cerr << "something went wrong... tokens is empty" << std::endl;
            }
            line_i++;
//            std::cerr << line_i << std::endl;
            if (line_i % 10000000 == 0) {
                std::cerr << line_i / 1000000 << "M" << std::endl;
            }
        }
        //file.close();

        // dump reads lengths file:
        dump(counts, nr_bins, output_prefix + ".homog");
        printf("finished %d reads\n", line_i);

    }
    catch (std::exception &e) {
        std::cerr << "failed calculating homog" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}


