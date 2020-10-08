#include <include/io.hpp>

#include <fstream>
#include <htslib/sam.h>

#define MAX_BUFFER_SIZE 1048576

using namespace std;

string loadFromFile(const string &fileName) {
    ifstream ifs(fileName.c_str(), ifstream::in);
    char buffer[MAX_BUFFER_SIZE];
    string sequence = "";
    while (!ifs.eof()) {
        ifs.getline(buffer, MAX_BUFFER_SIZE);
        if (buffer[0] == '>') {
            //  header += string(buffer) + '\n';
            continue;
        }
        sequence += string(buffer);
    }
    return sequence;
}

void
loadAlignFromSAM(const std::string &filePath, std::list<AlignPair> &aligns) {

    htsFile *inFile = sam_open(filePath.c_str(), "r");
    bam_hdr_t *head = sam_hdr_read(inFile);
    bam1_t *content = bam_init1();

    aligns.clear();
    while (sam_read1(inFile, head, content) >= 0) {
        AlignPair a;
        // TODO: bam_get_seq returns a 2 bit codification
        //       using bam_seqi(s,i) we can get the 4 bit
        //       encoding for  base 'i' of sequence 's'
        //       we then need to convert into a sequence
        //a.first = bam_get_seq(content);
        a.second = content->core.pos;
        aligns.push_back(a);
    }

    bam_destroy1(content);
    bam_hdr_destroy(head);
    sam_close(inFile);
}
