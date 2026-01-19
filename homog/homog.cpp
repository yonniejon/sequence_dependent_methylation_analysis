#include <queue>
#include "homog.h"

/**
 * @brief Splits a tab-separated string into a vector of strings.
 * @param line The input string to tokenize.
 * @return A vector containing the tab-separated tokens.
 */
std::vector <std::string> line2tokens(std::string &line) {
    std::vector <std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while (getline(lineStream, cell, '\t'))
        result.push_back(cell);
    return result;
}

/**
 * @brief Utility to print a vector of strings to stderr for debugging.
 * @param vec The vector to print.
 */
void print_vec(std::vector <std::string> &vec) {
    for (auto &j: vec)
        std::cerr << j << "\t";
    std::cerr << std::endl;
}

/**
 * @brief Comparator for a min-heap based on the end position of a Block.
 */
auto end_compare = [](Block a, Block b) { return a.end > b.end; };

/**
 * @brief Checks if a genomic read overlaps with a defined block.
 * @return True if there is any intersection or containment.
 */
bool intersects_block(int read_start, int read_end, int block_start, int block_end){
    bool is_inside = (read_start >= block_start && read_start < block_end) || (read_end >= block_start && read_end < block_end);
    bool read_contains = read_start <= block_start && read_end >= block_end;
    return is_inside || read_contains;
}

/**
 * @brief Allocates and zeroes out an array for storing bin counts.
 * @param len The number of bins (width) per block.
 * @return Pointer to the allocated memory.
 */
int32_t *Homog::init_array(int len) {
    int *arr = new int32_t[nr_blocks * len];
    std::fill_n(arr, nr_blocks * len, 0);
    return arr;
}

/**
 * @brief Constructor for the Homog class.
 * Initializes methylation ranges, loads genomic blocks, and prepares the count matrix.
 */
Homog::Homog(std::string in_output_prefix, std::string in_blocks_path, std::vector<float> in_range,
             int in_min_cpgs, bool deb) {
    min_cpgs = in_min_cpgs;
    output_prefix = in_output_prefix;
    blocks_path = in_blocks_path;
    range = in_range;
    debug = deb;

    nr_bins = range.size() - 1;

    // Load coordinates and CpG indices from the blocks file
    int r = read_blocks();
    nr_blocks = borders.size();

    // Allocate memory for the results matrix (Blocks x Bins)
    counts = init_array(nr_bins);
}

Homog::~Homog() {
    // Memory cleanup handled by the system or specific delete calls if uncommented.
}

/**
 * @brief Parses the input stream of blocks and populates the internal structures.
 * Validates format, CpG indices, and ensures monotonicity.
 */
int Homog::blocks_helper(std::istream &instream) {
    std::vector <std::string> tokens;
    std::string line;
    int cur_start = 0, cur_end = 0;
    int bi = 0;
    while (std::getline(instream, line)) {
        // Skip whitespace and comment lines
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

/**
 * @brief Loads the blocks file, handling both raw text and gzipped inputs.
 */
int Homog::read_blocks() {
    std::cout << "loading blocks..." << std::endl;

    if (hasEnding(blocks_path, ".gz")) {
        std::ifstream file(blocks_path, std::ios_base::in | std::ios_base::binary);
        boost::iostreams::filtering_streambuf <boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(file);
        std::istream instream(&inbuf);
        blocks_helper(instream);
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

/**
 * @brief Writes the accumulated counts to a TSV file and indexes it with tabix.
 */
void Homog::dump(int32_t *data, int width, std::string out_path) {
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

    // Finalize with tabix for genomic accessibility
    system(("bgzip -f " + out_path + " && tabix -p bed " + out_path + ".gz").c_str());
}

/**
 * @brief Internal logic to calculate methylation fraction and update bin counts.
 * Filters by min_cpgs and specified methylation range.
 */
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
        return;
    }

    float meth = (float) nrC / (float) (nrC + nrT);
    if (meth < range[0])
        return;

    int bin_ind = 0;
    for (bin_ind = 0; bin_ind < range.size() - 1; bin_ind++) {
        if ((meth >= range[bin_ind]) && (meth < range[bin_ind + 1])) {
            break;
        }
    }
    if (bin_ind == range.size() - 1) {
        bin_ind--;
    }
    counts[block_ind * nr_bins + bin_ind] += count;
}

/**
 * @brief Helper to update counts for a specific block ID using a pattern string.
 */
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
        return;
    }

    float meth = (float) nrC / (float) (nrC + nrT);
    if (meth < range[0])
        return;

    int bin_ind = 0;
    for (bin_ind = 0; bin_ind < range.size() - 1; bin_ind++) {
        if ((meth >= range[bin_ind]) && (meth < range[bin_ind + 1])) {
            break;
        }
    }
    if (bin_ind == range.size() - 1) {
        bin_ind--;
    }
    counts[block_id * nr_bins + bin_ind] += count;
}

/**
 * @brief Updates a single block index.
 */
void Homog::update_block(int block_ind, std::string pat, int32_t count) {
    int len = pat.length();
    if (len < min_cpgs) {
        return;
    }
    update_m2(block_ind, pat, count);
}

/**
 * @brief Processes a read pattern against all blocks currently in the active queue.
 * Handles sub-string extraction for reads partially overlapping block boundaries.
 */
void Homog::update_blocks(std::string pat, int32_t count, int read_start, int read_end) {
    int len = pat.length();
    // skip reads with less then min_cpgs sites:
    if (len < min_cpgs) {
        return;
    }

    // Check intersection with all blocks in the sliding window queue
    for (Block block : cur_block_queue){
        if (intersects_block(read_start, read_end, block.start, block.end)){
            if ((read_start >= block.start) && (read_start < block.end)) {
                int head_size = std::min(block.end - read_start, (int) pat.length());
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

/**
 * @brief Entry point for parsing a PAT file.
 * Iterates through standard input, processes each line, and dumps results.
 */
void Homog::parse(std::string pat_path) {
    try {
        int line_i = 0;
        std::ifstream infile(pat_path);

        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines

            std::vector <std::string> tokens = line2tokens(line_str);
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

            // Progress reporting
            if (line_i % 10000000 == 0) {
                std::cerr << line_i / 1000000 << "M" << std::endl;
            }
        }

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