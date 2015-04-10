/**
 * Author: Andreas Wilm
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 *
 */
#ifndef STATS_H
#define STATS_H


namespace clustalw
{
    
using namespace std;


/** Log simple statistics about input and output
 *
 * @section LICENSE
 * FIXME: please add license text
 *
 * @section DESCRIPTION
 * The hashing functions have to be enabled at compile/configure time
 * otherwise they will not be used.
 *
 * 
 */
class Stats {
    
public:
    /* FUNCTIONS */
    
    /** constructor */
    Stats();
    /** destructor */
    ~Stats();

    /**  set the name of the file where stats are logged to
     * @param f filename
     */
    void setStatsFile(string f) {logfilename=f;};
    
    /**
     * return name of used stats file
     * @return file name
     */
    string getStatsFile() {return logfilename;};
    
    /** enable stats logging
     */
    void setEnabled(bool b) {enabled=b;};
    
    /** check if stats logging is enabled
     * @return true if enabled, false otherwise
     */    
    bool isEnabled() {return enabled;};
    
    /** log the command line arguments
     * @param argc number arguments
     * @param argv argument array
     */
    void logCmdLine(int argc, char **argv);
    
    /** log statistics about sequence input
     * @param alnObj: input "alignment" (confusing isn't it?)
     */
    void logInputSeqStats(Alignment *alnObj);
    
    /** log statistics about sequence input
     * @param alnObj: input "alignment" (confusing isn't it?)
     */
    void logAlignedSeqStats(Alignment *alnObj);

    
    /* ATTRIBUTES */


    
private:
    
    /* FUNCTIONS */

    /** Compute pairwise identity of two sequences
     *
     * Adopted from Sean Eddy'ssquid:aligneval.c
     * He defines pairwise identity as fraction:
     * idents / min(len1, len2)
     *
     * Function is case sensitive!
     *
     * @param alnObj alignment object
     * @param s1 index of first sequence
     * @param s2 index of first sequence
     * @return pairwise identity as float
     */
    float pairwiseIdentity(Alignment *alnObj, int s1, int s2);
    
#ifdef HAVE_MHASH_H

    /** compute md5 hash from string
     * @param thread string to hash
     * @return hash string, or NULL on failure
     */
    char * Md5Hash(const char *thread);

    /** compute md5 hash for a sequence
     * @param alnObj alignment object
     * @param s index of sequence
     * @return md5 hash string
     */
    string Md5ForSeq(Alignment *alnObj, int s);
    
    /** compute md5 hash for all sequences
     * sorts and concatenates all sequences beforehand
     * thus sequence order doesn't matter
     * 
     * @param alnObj alignment object
     * @param s index of sequence
     * @return md5 hash string
     */
    string ConcatInputHash(Alignment *alnObj);
#endif

    /* ATTRIBUTES */

    /** the filename used for logging statistics */
    string logfilename;
    /** logging of statistics enabled */
    bool enabled;
};
}
#endif
