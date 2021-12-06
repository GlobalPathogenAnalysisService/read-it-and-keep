#include "CLI11.hpp"
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_map>
#include <zlib.h>
#include "gzstream.h"
#include "minimap.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


const std::string READ_IT_AND_KEEP_VERSION = "0.1.0";

const std::map<unsigned int, std::string> ERRORS = {
    {64, "Error parsing command line options or input file not found"},
    {65, "Error opening reads input file"},
    {66, "Error opening file for writing"},
    {67, "Reads file 1 has more reads than reads file 2"},
    {68, "Reads file 2 has more reads than reads file 1"},
};


void exitWithError(unsigned int errorCode);

class CommandLineOptions {
public:
    CommandLineOptions(int argc, char *argv[]);
    unsigned int minMapLength;
    float minMapLengthPercent;
    std::string refFasta;
    std::string readsIn1;
    std::string readsIn2;
    std::string readsOutprefix1;
    std::string readsOutprefix2;
    std::string outprefix;
    std::string tech;
    bool debug;
private:
    int parseCommandLineOpts(int argc, char *argv[]);
};


class QueryReads {
public:
    QueryReads(std::string& filenameIn, std::string& filenameOutprefix, CommandLineOptions& options);
    ~QueryReads();
    bool readNext();
    void resetForNextRefSeq();
    void map(mm_idx_t* mi, mm_mapopt_t& mopt);
    void clearMapping();
    void write();
    void writeDebug();
    bool isGoodMapping();
    kseq_t* readPtr;
    bool isBeingUsed;

private:
    gzFile filehandleIn_;
    ogzstream filehandleOut_;
    std::ofstream filehandleDebugOut_;
    mm_reg1_t* reg_;
    int n_reg_;
    mm_tbuf_t* tbuf_;
    unsigned int minMapLength_;
    float minMapLengthPercent_;
    bool debug_;
    bool hasQualScores_;
};


class Stats {
public:
    Stats();
    void toStdout();
    unsigned long long int readsIn1;
    unsigned long long int readsIn2;
    unsigned long long int readsKept1;
    unsigned long long int readsKept2;
};



int main(int argc, char *argv[])
{
    CommandLineOptions options(argc, argv);
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    int n_threads = 1;

    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    if (options.tech == "illumina") {
        mm_set_opt("sr", &iopt, &mopt);
    }
    else if (options.tech == "ont") {
        mm_set_opt("map-ont", &iopt, &mopt);
    }

    //mopt.flag |= MM_F_CIGAR; // perform alignment

    QueryReads queryReads1(options.readsIn1, options.readsOutprefix1, options);
    QueryReads queryReads2(options.readsIn2, options.readsOutprefix2, options);
    Stats stats;

    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(options.refFasta.c_str(), &iopt, 0);
    mm_idx_t *mi;
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        queryReads1.resetForNextRefSeq();
        queryReads2.resetForNextRefSeq();
        while (queryReads1.readNext()) {
            stats.readsIn1++;
            bool keptReads = false;
            if (queryReads2.isBeingUsed) {
                if (not queryReads2.readNext()) {
                    std::cerr << "Error reading mate of read '"
                              << queryReads1.readPtr->name.s << "'\n";
                    exitWithError(67);
                }
                stats.readsIn2++;
            }
            queryReads1.map(mi, mopt);
            if (queryReads1.isGoodMapping()) {
                queryReads1.write();
                keptReads = true;
                stats.readsKept1++;
                if (queryReads2.isBeingUsed) {
                    queryReads2.write();
                    stats.readsKept2++;
                }
            }
            else if (queryReads2.isBeingUsed) {
                queryReads2.map(mi, mopt);
                if (queryReads2.isGoodMapping()) {
                    queryReads1.write();
                    queryReads2.write();
                    keptReads = true;
                    stats.readsKept1++;
                    stats.readsKept2++;
                }
                queryReads2.clearMapping();
            }

            queryReads1.clearMapping();
            if (not keptReads) {
                   queryReads1.writeDebug();
                   queryReads2.writeDebug();
            }

            if (stats.readsIn1 % 100000 == 0) {
                std::cerr << "Processed " << stats.readsIn1 << " reads (or read pairs)" << std::endl;
            }
        }
        mm_idx_destroy(mi);
    }
    if (queryReads2.isBeingUsed and queryReads2.readNext()) {
        std::cerr << "Used all reads from first file, but got another read '"
                  << queryReads2.readPtr->name.s
                  << "' from second reads file\n";
        exitWithError(68);
    }
    mm_idx_reader_close(r); // close the index reader
    stats.toStdout();
    return 0;
}

CommandLineOptions::CommandLineOptions(int argc, char *argv[]) {
    minMapLength = 50;
    minMapLengthPercent = 50.0;
    refFasta = "";
    debug = false;
    tech = "illumina";
    auto code = this->parseCommandLineOpts(argc, argv);
    if (code != 0) { // error parsing command line options
        exitWithError(64);
    }
    else if (refFasta == "") { // user must have used --help
        exit(0);
    }
}

void exitWithError(unsigned int errorCode) {
    std::cerr << ERRORS.at(errorCode) << ". Cannot continue." << std::endl;
    exit(errorCode);
}


int CommandLineOptions::parseCommandLineOpts(int argc, char *argv[]) {
    CLI::App app{};

    app.add_option("--tech", tech, "Sequencing technology, must be 'illumina' or 'ont' [illumina]");

    app.add_option("--ref_fasta", refFasta, "Reference genome FASTA filename")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("--reads1", readsIn1, "Name of first reads file")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("--reads2", readsIn2, "Name of second reads file, ie mates file for paired reads")
        ->check(CLI::ExistingFile);

    app.add_option("-o,--outprefix", outprefix, "Prefix of output files")
        ->required();

    app.add_flag("--debug", debug, "Debug mode. More verbose and writes debugging files");

    app.add_option("--min_map_length", minMapLength, "Minimum length of match required to keep a read in bp [50]");

    app.add_option("--min_map_length_pc", minMapLengthPercent, "Minimum length of match required to keep a read, as a percent of the read length [50.0]");

    app.add_flag_callback("-V,--version", []() {
        std::cout << "readItAndKeep version " << READ_IT_AND_KEEP_VERSION << std::endl;
        throw(CLI::Success {});
        },
        "Show version and exit");

    CLI11_PARSE(app, argc, argv);

    if (tech != "illumina" and tech != "ont") {
        std::cerr << "--tech option must be 'illumina' or 'ont', I got: '"
                  << tech << "'\n";
        exitWithError(64);
    }

    if (readsIn2 == "" ) {
        readsOutprefix1 = outprefix + ".reads";
        readsOutprefix2 = "";
    }
    else {
        readsOutprefix1 = outprefix + ".reads_1";
        readsOutprefix2 = outprefix + ".reads_2";
    }
    return 0;
}


QueryReads::QueryReads(std::string& filenameIn, std::string& filenameOutprefix, CommandLineOptions& options) {
    if (filenameIn != "") {
        filehandleIn_ = gzopen(filenameIn.c_str(), "r");
        if (!filehandleIn_) {
            std::cerr << "Error opening input reads file '" << filenameIn << "'\n";
            exitWithError(65);
        }

        // We need to know if the reads file is FASTA or FASTQ format.
        // Decide by: if char of the file being '@' then it's FASTQ, otherwise
        // it's FASTA
        std::string filenameOut = filenameOutprefix;
        if (gzgetc(filehandleIn_) == '@') {
            hasQualScores_ = true;
            filenameOut += ".fastq.gz";
        }
        else {
            hasQualScores_ = false;
            filenameOut += ".fasta.gz";
        }
        gzrewind(filehandleIn_);

        readPtr = kseq_init(filehandleIn_);
        isBeingUsed = true;
        tbuf_ = mm_tbuf_init();
        filehandleOut_.open(filenameOut.c_str());
        if (not filehandleOut_) {
            std::cerr << "Error opening output reads file '" << filenameOut << "'\n";
            exitWithError(66);
        }
    }
    else {
        isBeingUsed = false;
    }
    minMapLength_ = options.minMapLength;
    minMapLengthPercent_ = options.minMapLengthPercent;
    debug_ = options.debug;
    if (isBeingUsed and debug_) {
        std::string filename = filenameOutprefix + ".debug";
        filehandleDebugOut_.open(filename);
        if (not filehandleDebugOut_.is_open()) {
            std::cerr << "Error opening debug output file '" << filename << "'\n";
            exitWithError(66);
        }
    }
}

void QueryReads::resetForNextRefSeq() {
    if (isBeingUsed) {
        mm_tbuf_destroy(tbuf_);
        tbuf_ = mm_tbuf_init();
        gzrewind(filehandleIn_);
        kseq_rewind(readPtr);
    }
}

bool QueryReads::readNext() {
    return (isBeingUsed and kseq_read(readPtr) >= 0);
}

void QueryReads::map(mm_idx_t* mi, mm_mapopt_t& mopt) {
    if (isBeingUsed) {
        reg_ = mm_map(mi, readPtr->seq.l, readPtr->seq.s, &n_reg_, tbuf_, &mopt, 0);
    }
}

void QueryReads::clearMapping() {
    free(reg_);
}

void QueryReads::write() {
    if (hasQualScores_) {
        filehandleOut_ << '@' << readPtr->name.s;
        if (readPtr->comment.l > 0) {
            filehandleOut_ << ' ' << readPtr->comment.s;
        }
        filehandleOut_ << '\n'
                       << readPtr->seq.s << '\n'
                       << "+\n"
                       << readPtr->qual.s << '\n';
    }
    else {
        filehandleOut_ << '>' << readPtr->name.s;
        if (readPtr->comment.l > 0) {
            filehandleOut_ << ' ' << readPtr->comment.s;
        }
        filehandleOut_ << '\n' << readPtr->seq.s << '\n';
    }
}


void QueryReads::writeDebug() {
    if (not isBeingUsed) {
        return;
    }
    filehandleDebugOut_ << "REJECTED_READ\t" << readPtr->name.s
                        << '\t' << readPtr->comment.s
                        << '\t' << readPtr->seq.s
                        << '\t' << readPtr->qual.s << '\n';
}


bool QueryReads::isGoodMapping() {
    int j;
    if (debug_) {
        filehandleDebugOut_ << "MAPPING_COUNT\t" << readPtr->name.s
                            << '\t' << readPtr->comment.s
                            << '\t' << n_reg_ << '\n';
    }
    for (j = 0; j < n_reg_; ++j) {
        mm_reg1_t *r = &reg_[j];
        //assert(r->p); // with MM_F_CIGAR, this should not be NULL
        // note: regardless of the strand, we always have query start < end
        bool ok = (minMapLength_ <= (unsigned int) (r->qe - r->qs) or minMapLengthPercent_ <= 100.0 * (r->qe - r->qs) / readPtr->seq.l);
        if (debug_) {
            filehandleDebugOut_ << "MAPPING\t" << readPtr->name.s
                       << '\t' << readPtr->comment.s
                       << '\t' << r->qs << "-" << r->qe
                       << '\t' << "+-"[r->rev]
                       << "\tpass:" << ok << '\n';
        }
        free(r->p);
        if (ok) {
            return true;
        }
    }
    return false;
}

QueryReads::~QueryReads() {
    if (isBeingUsed) {
        gzclose(filehandleIn_);
        filehandleOut_.close();
        kseq_destroy(readPtr);
        mm_tbuf_destroy(tbuf_);
        if (debug_) {
            filehandleDebugOut_.close();
        }
    }
}


Stats::Stats() {
    readsIn1 = 0;
    readsIn2 = 0;
    readsKept1 = 0;
    readsKept2 = 0;
}


void Stats::toStdout() {
    std::cout << "Input reads file 1\t" << readsIn1 << '\n'
              << "Input reads file 2\t" << readsIn2 << '\n'
              << "Kept reads 1\t" << readsKept1 << '\n'
              << "Kept reads 2\t" << readsKept2 << '\n';
}
