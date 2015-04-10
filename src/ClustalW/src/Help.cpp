/**
 * Author: Andreas Wilm
 *
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.
 */
/**
 * Changes:
 * 2007-12-12: Andreas Wilm (UCD): initial implementation
 *
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif

#include <iostream>
#include <vector>
#include "clustalw_version.h"
#include "Help.h"

using namespace std;


Help::Help()
{
    section s;
    string version = CLUSTALW_VERSION;

    // all sections were added in order of appearance in original
    // clustalw_help
    
    // To add new entries use the template below.
    // Paste content as text,
    // escape all " with \",
    // add " at the beginning of the line,
    // and \n" at the end of the line
    // template:
    //  s.marker = "";
    //  s.title = "";
    //  s.content = "";
    //  sections.push_back(s);
    
    s.marker = "NEW";
    s.title = "NEW FEATURES/OPTIONS";
    s.content = "==UPGMA=="        
" \n"
" The UPGMA algorithm has been added to allow faster tree construction. The user now\n"
" has the choice of using Neighbour Joining or UPGMA. The default is still NJ, but the\n"
" user can change this by setting the clustering parameter.\n"
" \n"
" -CLUSTERING=   :NJ or UPGMA\n"
" \n"
"==ITERATION==\n"
"\n"
" A remove first iteration scheme has been added. This can be used to improve the final\n"
" alignment or improve the alignment at each stage of the progressive alignment. During the \n"
" iteration step each sequence is removed in turn and realigned. If the resulting alignment \n"
" is better than the  previous alignment it is kept. This process is repeated until the score\n"
" converges (the  score is not improved) or until the maximum number of iterations is \n"
" reached. The user can  iterate at each step of the progressive alignment by setting the \n"
" iteration parameter to  TREE or just on the final alignment by seting the iteration \n"
" parameter to ALIGNMENT. The default is no iteration. The maximum number of  iterations can \n"
" be set using the numiter parameter. The default number of iterations is 3.\n"
"  \n"
" -ITERATION=    :NONE or TREE or ALIGNMENT\n"
" \n"
" -NUMITER=n     :Maximum number of iterations to perform\n"
" \n"
"==HELP==\n"
" \n"        
" -FULLHELP      :Print out the complete help content\n"
" \n"
"==MISC==\n"
"\n"
" -MAXSEQLEN=n   :Maximum allowed sequence length\n"
" \n"        
" -QUIET         :Reduce console output to minimum\n"
" \n"
" -STATS=file    :Log some alignents statistics to file\n"
"";
    sections.push_back(s);

    
    s.marker = "1";
    s.title = "General help for CLUSTAL W (" + version + ")";
    s.content = "Clustal W is a general purpose multiple alignment program for DNA or proteins.\n"
"\n"
"SEQUENCE INPUT:  all sequences must be in 1 file, one after another.  \n"
"7 formats are automatically recognised: NBRF-PIR, EMBL-SWISSPROT, \n"
"Pearson (Fasta), Clustal (*.aln), GCG-MSF (Pileup), GCG9-RSF and GDE flat file.\n"
"All non-alphabetic characters (spaces, digits, punctuation marks) are ignored\n"
"except \"-\" which is used to indicate a GAP (\".\" in MSF-RSF).  \n"
"\n"
"To do a MULTIPLE ALIGNMENT on a set of sequences, use item 1 from this menu to \n"
"INPUT them; go to menu item 2 to do the multiple alignment.\n"
"\n"
"PROFILE ALIGNMENTS (menu item 3) are used to align 2 alignments.  Use this to\n"
"add a new sequence to an old alignment, or to use secondary structure to guide \n"
"the alignment process.  GAPS in the old alignments are indicated using the \"-\" \n"
"character.   PROFILES can be input in ANY of the allowed formats; just \n"
"use \"-\" (or \".\" for MSF-RSF) for each gap position.\n"
"\n"
"PHYLOGENETIC TREES (menu item 4) can be calculated from old alignments (read in\n"
"with \"-\" characters to indicate gaps) OR after a multiple alignment while the \n"
"alignment is still in memory.\n"
"\n"
"\n"
"The program tries to automatically recognise the different file formats used\n"
"and to guess whether the sequences are amino acid or nucleotide.  This is not\n"
"always foolproof.\n"
"\n"
"FASTA and NBRF-PIR formats are recognised by having a \">\" as the first \n"
"character in the file.  \n"
"\n"
"EMBL-Swiss Prot formats are recognised by the letters\n"
"ID at the start of the file (the token for the entry name field).  \n"
"\n"
"CLUSTAL format is recognised by the word CLUSTAL at the beginning of the file.\n"
"\n"
"GCG-MSF format is recognised by one of the following:\n"
"       - the word PileUp at the start of the file. \n"
"       - the word !!AA_MULTIPLE_ALIGNMENT or !!NA_MULTIPLE_ALIGNMENT\n"
"         at the start of the file.\n"
"       - the word MSF on the first line of the line, and the characters ..\n"
"         at the end of this line.\n"
"\n"
"GCG-RSF format is recognised by the word !!RICH_SEQUENCE at the beginning of\n"
"the file.\n"
"\n"
"\n"
"If 85% or more of the characters in the sequence are from A,C,G,T,U or N, the\n"
"sequence will be assumed to be nucleotide.  This works in 97.3% of cases\n"
"but watch out!\n"
"";
    sections.push_back(s);

    
    s.marker = "2";
    s.title = "Help for multiple alignments";
    s.content = "If you have already loaded sequences, use menu item 1 to do the complete\n"
"multiple alignment.  You will be prompted for 2 output files: 1 for the \n"
"alignment itself; another to store a dendrogram that describes the similarity\n"
"of the sequences to each other.\n"
"\n"
"Multiple alignments are carried out in 3 stages (automatically done from menu\n"
"item 1 ...Do complete multiple alignments now):\n"
"\n"
"1) all sequences are compared to each other (pairwise alignments);\n"
"\n"
"2) a dendrogram (like a phylogenetic tree) is constructed, describing the\n"
"approximate groupings of the sequences by similarity (stored in a file).\n"
"\n"
"3) the final multiple alignment is carried out, using the dendrogram as a guide.\n"
"\n"
"\n"
"PAIRWISE ALIGNMENT parameters control the speed-sensitivity of the initial\n"
"alignments.\n"
"\n"
"MULTIPLE ALIGNMENT parameters control the gaps in the final multiple alignments.\n"
"\n"
"\n"
"RESET GAPS (menu item 7) will remove any new gaps introduced into the sequences\n"
"during multiple alignment if you wish to change the parameters and try again.\n"
"This only takes effect just before you do a second multiple alignment.  You\n"
"can make phylogenetic trees after alignment whether or not this is ON.\n"
"If you turn this OFF, the new gaps are kept even if you do a second multiple\n"
"alignment. This allows you to iterate the alignment gradually.  Sometimes, the \n"
"alignment is improved by a second or third pass.\n"
"\n"
"SCREEN DISPLAY (menu item 8) can be used to send the output alignments to the \n"
"screen as well as to the output file.\n"
"\n"
"You can skip the first stages (pairwise alignments; dendrogram) by using an\n"
"old dendrogram file (menu item 3); or you can just produce the dendrogram\n"
"with no final multiple alignment (menu item 2).\n"
"\n"
"\n"
"OUTPUT FORMAT: Menu item 9 (format options) allows you to choose from 6 \n"
"different alignment formats (CLUSTAL, GCG, NBRF-PIR, PHYLIP, GDE, NEXUS, and FASTA).  \n"
"\n"
"";
    sections.push_back(s);


    s.marker = "3";
    s.title = "Help for pairwise alignment parameters";
    s.content = "A distance is calculated between every pair of sequences and these are used to\n"
"construct the dendrogram which guides the final multiple alignment. The scores\n"
"are calculated from separate pairwise alignments. These can be calculated using\n"
"2 methods: dynamic programming (slow but accurate) or by the method of Wilbur\n"
"and Lipman (extremely fast but approximate). \n"
"\n"
"You can choose between the 2 alignment methods using menu option 8.  The\n"
"slow-accurate method is fine for short sequences but will be VERY SLOW for \n"
"many (e.g. >100) long (e.g. >1000 residue) sequences.   \n"
"\n"
"SLOW-ACCURATE alignment parameters:\n"
"	These parameters do not have any affect on the speed of the alignments. \n"
"They are used to give initial alignments which are then rescored to give percent\n"
"identity scores.  These % scores are the ones which are displayed on the \n"
"screen.  The scores are converted to distances for the trees.\n"
"\n"
"1) Gap Open Penalty:      the penalty for opening a gap in the alignment.\n"
"2) Gap extension penalty: the penalty for extending a gap by 1 residue.\n"
"3) Protein weight matrix: the scoring table which describes the similarity\n"
"                          of each amino acid to each other.\n"
"4) DNA weight matrix:     the scores assigned to matches and mismatches \n"
"                          (including IUB ambiguity codes).\n"
"\n"
"\n"
"FAST-APPROXIMATE alignment parameters:\n"
"\n"
"These similarity scores are calculated from fast, approximate, global align-\n"
"ments, which are controlled by 4 parameters.   2 techniques are used to make\n"
"these alignments very fast: 1) only exactly matching fragments (k-tuples) are\n"
"considered; 2) only the 'best' diagonals (the ones with most k-tuple matches)\n"
"are used.\n"
"\n"
"K-TUPLE SIZE:  This is the size of exactly matching fragment that is used. \n"
"INCREASE for speed (max= 2 for proteins; 4 for DNA), DECREASE for sensitivity.\n"
"For longer sequences (e.g. >1000 residues) you may need to increase the default.\n"
"\n"
"GAP PENALTY:   This is a penalty for each gap in the fast alignments.  It has\n"
"little affect on the speed or sensitivity except for extreme values.\n"
"\n"
"TOP DIAGONALS: The number of k-tuple matches on each diagonal (in an imaginary\n"
"dot-matrix plot) is calculated.  Only the best ones (with most matches) are\n"
"used in the alignment.  This parameter specifies how many.  Decrease for speed;\n"
"increase for sensitivity.\n"
"\n"
"WINDOW SIZE:  This is the number of diagonals around each of the 'best' \n"
"diagonals that will be used.  Decrease for speed; increase for sensitivity.\n"
"";
    sections.push_back(s);


    s.marker = "4";
    s.title = "Help for multiple alignment parameters";
    s.content = "These parameters control the final multiple alignment. This is the core of the\n"
"program and the details are complicated. To fully understand the use of the\n"
"parameters and the scoring system, you will have to refer to the documentation.\n"
"\n"
"Each step in the final multiple alignment consists of aligning two alignments \n"
"or sequences.  This is done progressively, following the branching order in \n"
"the GUIDE TREE.  The basic parameters to control this are two gap penalties and\n"
"the scores for various identical-non-indentical residues.  \n"
"\n"
"1) and 2) The GAP PENALTIES are set by menu items 1 and 2. These control the \n"
"cost of opening up every new gap and the cost of every item in a gap. \n"
"Increasing the gap opening penalty will make gaps less frequent. Increasing \n"
"the gap extension penalty will make gaps shorter. Terminal gaps are not \n"
"penalised.\n"
"\n"
"3) The DELAY DIVERGENT SEQUENCES switch delays the alignment of the most\n"
"distantly related sequences until after the most closely related sequences have \n"
"been aligned.   The setting shows the percent identity level required to delay\n"
"the addition of a sequence; sequences that are less identical than this level\n"
"to any other sequences will be aligned later.\n"
"\n"
"\n"
"\n"
"4) The TRANSITION WEIGHT gives transitions (A <--> G or C <--> T \n"
"i.e. purine-purine or pyrimidine-pyrimidine substitutions) a weight between 0\n"
"and 1; a weight of zero means that the transitions are scored as mismatches,\n"
"while a weight of 1 gives the transitions the match score. For distantly related\n"
"DNA sequences, the weight should be near to zero; for closely related sequences\n"
"it can be useful to assign a higher score.\n"
"\n"
"\n"
"5) PROTEIN WEIGHT MATRIX leads to a new menu where you are offered a choice of\n"
"weight matrices. The default for proteins in version 1.8 is the PAM series \n"
"derived by Gonnet and colleagues. Note, a series is used! The actual matrix\n"
"that is used depends on how similar the sequences to be aligned at this \n"
"alignment step are. Different matrices work differently at each evolutionary\n"
"distance. \n"
"\n"
"6) DNA WEIGHT MATRIX leads to a new menu where a single matrix (not a series)\n"
"can be selected. The default is the matrix used by BESTFIT for comparison of\n"
"nucleic acid sequences.\n"
"\n"
"Further help is offered in the weight matrix menu.\n"
"\n"
"\n"
"7)  In the weight matrices, you can use negative as well as positive values if\n"
"you wish, although the matrix will be automatically adjusted to all positive\n"
"scores, unless the NEGATIVE MATRIX option is selected.\n"
"\n"
"8) PROTEIN GAP PARAMETERS displays a menu allowing you to set some Gap Penalty\n"
"options which are only used in protein alignments.\n"
"";
    sections.push_back(s);


    s.marker = "A";
    s.title = "Help for protein gap parameters.";
    s.content = "1) RESIDUE SPECIFIC PENALTIES are amino acid specific gap penalties that reduce\n"
"or increase the gap opening penalties at each position in the alignment or\n"
"sequence.  See the documentation for details.  As an example, positions that \n"
"are rich in glycine are more likely to have an adjacent gap than positions that\n"
"are rich in valine.\n"
"\n"
"2) 3) HYDROPHILIC GAP PENALTIES are used to increase the chances of a gap within\n"
"a run (5 or more residues) of hydrophilic amino acids; these are likely to\n"
"be loop or random coil regions where gaps are more common.  The residues that \n"
"are \"considered\" to be hydrophilic are set by menu item 3.\n"
"\n"
"4) GAP SEPARATION DISTANCE tries to decrease the chances of gaps being too\n"
"close to each other. Gaps that are less than this distance apart are penalised\n"
"more than other gaps. This does not prevent close gaps; it makes them less\n"
"frequent, promoting a block-like appearance of the alignment.\n"
"\n"
"5) END GAP SEPARATION treats end gaps just like internal gaps for the purposes\n"
"of avoiding gaps that are too close (set by GAP SEPARATION DISTANCE above).\n"
"If you turn this off, end gaps will be ignored for this purpose.  This is\n"
"useful when you wish to align fragments where the end gaps are not biologically\n"
"meaningful.\n"
"";
    sections.push_back(s);

    
    s.marker = "5";
    s.title = "Help for output format options.";
    s.content = "Several output formats are offered. You can choose any (or all 6 if you wish).  \n"
"\n"
"CLUSTAL format output is a self explanatory alignment format.  It shows the\n"
"sequences aligned in blocks.  It can be read in again at a later date to\n"
"(for example) calculate a phylogenetic tree or add a new sequence with a \n"
"profile alignment.\n"
"\n"
"GCG output can be used by any of the GCG programs that can work on multiple\n"
"alignments (e.g. PRETTY, PROFILEMAKE, PLOTALIGN).  It is the same as the GCG\n"
".msf format files (multiple sequence file); new in version 7 of GCG.\n"
"\n"
"Fasta output cis widely used because of it's simplicity. Each sequence name is\n"
"preceeded by a '>'-sign. The sequence itself is printed out in the following lines\n"
"\n"
"PHYLIP format output can be used for input to the PHYLIP package of Joe \n"
"Felsenstein.  This is an extremely widely used package for doing every \n"
"imaginable form of phylogenetic analysis (MUCH more than the the modest intro-\n"
"duction offered by this program).\n"
"\n"
"NBRF-PIR:  this is the same as the standard PIR format with ONE ADDITION.  Gap\n"
"characters \"-\" are used to indicate the positions of gaps in the multiple \n"
"alignment.  These files can be re-used as input in any part of clustal that\n"
"allows sequences (or alignments or profiles) to be read in.  \n"
"\n"
"GDE:  this is the flat file format used by the GDE package of Steven Smith.\n"
"\n"
"NEXUS: the format used by several phylogeny programs, including PAUP and\n"
"MacClade.\n"
"\n"
"GDE OUTPUT CASE: sequences in GDE format may be written in either upper or\n"
"lower case.\n"
"\n"
"CLUSTALW SEQUENCE NUMBERS: residue numbers may be added to the end of the\n"
"alignment lines in clustalw format.\n"
"\n"
"OUTPUT ORDER is used to control the order of the sequences in the output\n"
"alignments.  By default, the order corresponds to the order in which the\n"
"sequences were aligned (from the guide tree-dendrogram), thus automatically\n"
"grouping closely related sequences. This switch can be used to set the order\n"
"to the same as the input file.\n"
"\n"
"PARAMETER OUTPUT: This option allows you to save all your parameter settings\n"
"in a parameter file. This file can be used subsequently to rerun Clustal W\n"
"using the same parameters.\n"
"";
    sections.push_back(s);

    
    s.marker = "6";
    s.title = "Help for profile and structure alignments";
    s.content = "By PROFILE ALIGNMENT, we mean alignment using existing alignments. Profile \n"
"alignments allow you to store alignments of your favourite sequences and add\n"
"new sequences to them in small bunches at a time. A profile is simply an\n"
"alignment of one or more sequences (e.g. an alignment output file from CLUSTAL\n"
"W). Each input can be a single sequence. One or both sets of input sequences\n"
"may include secondary structure assignments or gap penalty masks to guide the\n"
"alignment. \n"
"\n"
"The profiles can be in any of the allowed input formats with \"-\" characters\n"
"used to specify gaps (except for MSF-RSF where \".\" is used).\n"
"\n"
"You have to specify the 2 profiles by choosing menu items 1 and 2 and giving\n"
"2 file names.  Then Menu item 3 will align the 2 profiles to each other. \n"
"Secondary structure masks in either profile can be used to guide the alignment.\n"
"\n"
"Menu item 4 will take the sequences in the second profile and align them to\n"
"the first profile, 1 at a time.  This is useful to add some new sequences to\n"
"an existing alignment, or to align a set of sequences to a known structure.  \n"
"In this case, the second profile would not be pre-aligned.\n"
"\n"
"\n"
"The alignment parameters can be set using menu items 5, 6 and 7. These are\n"
"EXACTLY the same parameters as used by the general, automatic multiple\n"
"alignment procedure. The general multiple alignment procedure is simply a\n"
"series of profile alignments. Carrying out a series of profile alignments on\n"
"larger and larger groups of sequences, allows you to manually build up a\n"
"complete alignment, if necessary editing intermediate alignments.\n"
"\n"
"SECONDARY STRUCTURE OPTIONS. Menu Option 0 allows you to set 2D structure\n"
"parameters. If a solved structure is available, it can be used to guide the \n"
"alignment by raising gap penalties within secondary structure elements, so \n"
"that gaps will preferentially be inserted into unstructured surface loops.\n"
"Alternatively, a user-specified gap penalty mask can be supplied directly.\n"
"\n"
"A gap penalty mask is a series of numbers between 1 and 9, one per position in \n"
"the alignment. Each number specifies how much the gap opening penalty is to be \n"
"raised at that position (raised by multiplying the basic gap opening penalty\n"
"by the number) i.e. a mask figure of 1 at a position means no change\n"
"in gap opening penalty; a figure of 4 means that the gap opening penalty is\n"
"four times greater at that position, making gaps 4 times harder to open.\n"
"\n"
"The format for gap penalty masks and secondary structure masks is explained\n"
"in the help under option 0 (secondary structure options).\n"
"";
    sections.push_back(s);

    
    s.marker = "B";
    s.title = "Help for secondary structure - gap penalty masks";
    s.content = "The use of secondary structure-based penalties has been shown to improve the\n"
"accuracy of multiple alignment. Therefore CLUSTAL W now allows gap penalty \n"
"masks to be supplied with the input sequences. The masks work by raising gap \n"
"penalties in specified regions (typically secondary structure elements) so that\n"
"gaps are preferentially opened in the less well conserved regions (typically \n"
"surface loops).\n"
"\n"
"Options 1 and 2 control whether the input secondary structure information or\n"
"gap penalty masks will be used.\n"
"\n"
"Option 3 controls whether the secondary structure and gap penalty masks should\n"
"be included in the output alignment.\n"
"\n"
"Options 4 and 5 provide the value for raising the gap penalty at core Alpha \n"
"Helical (A) and Beta Strand (B) residues. In CLUSTAL format, capital residues \n"
"denote the A and B core structure notation. The basic gap penalties are\n"
"multiplied by the amount specified.\n"
"\n"
"Option 6 provides the value for the gap penalty in Loops. By default this \n"
"penalty is not raised. In CLUSTAL format, loops are specified by \".\" in the \n"
"secondary structure notation.\n"
"\n"
"Option 7 provides the value for setting the gap penalty at the ends of \n"
"secondary structures. Ends of secondary structures are observed to grow \n"
"and-or shrink in related structures. Therefore by default these are given \n"
"intermediate values, lower than the core penalties. All secondary structure \n"
"read in as lower case in CLUSTAL format gets the reduced terminal penalty.\n"
"\n"
"Options 8 and 9 specify the range of structure termini for the intermediate \n"
"penalties. In the alignment output, these are indicated as lower case. \n"
"For Alpha Helices, by default, the range spans the end helical turn. For \n"
"Beta Strands, the default range spans the end residue and the adjacent loop \n"
"residue, since sequence conservation often extends beyond the actual H-bonded\n"
"Beta Strand.\n"
"\n"
"CLUSTAL W can read the masks from SWISS-PROT, CLUSTAL or GDE format input\n"
"files. For many 3-D protein structures, secondary structure information is\n"
"recorded in the feature tables of SWISS-PROT database entries. You should\n"
"always check that the assignments are correct - some are quite inaccurate.\n"
"CLUSTAL W looks for SWISS-PROT HELIX and STRAND assignments e.g.\n"
"\n"
"FT   HELIX       100    115\n"
"FT   STRAND      118    119\n"
"\n"
"The structure and penalty masks can also be read from CLUSTAL alignment format \n"
"as comment lines beginning \"!SS_\" or \"!GM_\" e.g.\n"
"\n"
"!SS_HBA_HUMA    ..aaaAAAAAAAAAAaaa.aaaAAAAAAAAAAaaaaaaAaaa.........aaaAAAAAA\n"
"!GM_HBA_HUMA    112224444444444222122244444444442222224222111111111222444444\n"
"HBA_HUMA        VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK\n"
"\n"
"Note that the mask itself is a set of numbers between 1 and 9 each of which is \n"
"assigned to the residue(s) in the same column below. \n"
"\n"
"In GDE flat file format, the masks are specified as text and the names must\n"
"begin with \"SS_ or \"GM_.\n"
"\n"
"Either a structure or penalty mask or both may be used. If both are included in\n"
"an alignment, the user will be asked which is to be used.\n"
"";
    sections.push_back(s);

    
    s.marker = "C";
    s.title = "Help for secondary structure - gap penalty mask output options";
    s.content = "The options in this menu let you choose whether or not to include the masks\n"
"in the CLUSTAL W output alignments. Showing both is useful for understanding\n"
"how the masks work. The secondary structure information is itself very useful\n"
"in judging the alignment quality and in seeing how residue conservation\n"
"patterns vary with secondary structure.\n"
"";
    sections.push_back(s);


    s.marker = "7";
    s.title = "Help for phylogenetic trees";
    s.content = "1) Before calculating a tree, you must have an ALIGNMENT in memory. This can be\n"
"input in any format or you should have just carried out a full multiple\n"
"alignment and the alignment is still in memory. \n"
"\n"
"\n"
"*************** Remember YOU MUST ALIGN THE SEQUENCES FIRST!!!! ***************\n"
"\n"
"\n"
"The methods used are NJ (Neighbour Joining) and UPGMA. First\n"
"you calculate distances (percent divergence) between all pairs of sequence from\n"
"a multiple alignment; second you apply the NJ or UPGMA method to the distance matrix.\n"
"\n"
"2) EXCLUDE POSITIONS WITH GAPS? With this option, any alignment positions where\n"
"ANY of the sequences have a gap will be ignored. This means that 'like' will be\n"
"compared to 'like' in all distances, which is highly desirable. It also\n"
"automatically throws away the most ambiguous parts of the alignment, which are\n"
"concentrated around gaps (usually). The disadvantage is that you may throw away\n"
"much of the data if there are many gaps (which is why it is difficult for us to\n"
"make it the default).  \n"
"\n"
"\n"
"\n"
"3) CORRECT FOR MULTIPLE SUBSTITUTIONS? For small divergence (say <10%) this\n"
"option makes no difference. For greater divergence, it corrects for the fact\n"
"that observed distances underestimate actual evolutionary distances. This is\n"
"because, as sequences diverge, more than one substitution will happen at many\n"
"sites. However, you only see one difference when you look at the present day\n"
"sequences. Therefore, this option has the effect of stretching branch lengths\n"
"in trees (especially long branches). The corrections used here (for DNA or\n"
"proteins) are both due to Motoo Kimura. See the documentation for details.  \n"
"\n"
"Where possible, this option should be used. However, for VERY divergent\n"
"sequences, the distances cannot be reliably corrected. You will be warned if\n"
"this happens. Even if none of the distances in a data set exceed the reliable\n"
"threshold, if you bootstrap the data, some of the bootstrap distances may\n"
"randomly exceed the safe limit.  \n"
"\n"
"4) To calculate a tree, use option 4 (DRAW TREE NOW). This gives an UNROOTED\n"
"tree and all branch lengths. The root of the tree can only be inferred by\n"
"using an outgroup (a sequence that you are certain branches at the outside\n"
"of the tree .... certain on biological grounds) OR if you assume a degree\n"
"of constancy in the 'molecular clock', you can place the root in the 'middle'\n"
"of the tree (roughly equidistant from all tips).\n"
"\n"
"5) TOGGLE PHYLIP BOOTSTRAP POSITIONS\n"
"By default, the bootstrap values are correctly placed on the tree branches of\n"
"the phylip format output tree. The toggle allows them to be placed on the\n"
"nodes, which is incorrect, but some display packages (e.g. TreeTool, TreeView\n"
"and Phylowin) only support node labelling but not branch labelling. Care\n"
"should be taken to note which branches and labels go together.\n"
"\n"
"6) OUTPUT FORMATS: four different formats are allowed. None of these displays\n"
"the tree visually. Useful display programs accepting PHYLIP format include\n"
"NJplot (from Manolo Gouy and supplied with Clustal W), TreeView (Mac-PC), and\n"
"PHYLIP itself - OR get the PHYLIP package and use the tree drawing facilities\n"
"there. (Get the PHYLIP package anyway if you are interested in trees). The\n"
"NEXUS format can be read into PAUP or MacClade.\n"
"";
    sections.push_back(s);


    s.marker = "8";
    s.title = "Help for choosing a weight matrix";
    s.content = "For protein alignments, you use a weight matrix to determine the similarity of\n"
"non-identical amino acids.  For example, Tyr aligned with Phe is usually judged \n"
"to be 'better' than Tyr aligned with Pro.\n"
"\n"
"There are three 'in-built' series of weight matrices offered. Each consists of\n"
"several matrices which work differently at different evolutionary distances. To\n"
"see the exact details, read the documentation. Crudely, we store several\n"
"matrices in memory, spanning the full range of amino acid distance (from almost\n"
"identical sequences to highly divergent ones). For very similar sequences, it\n"
"is best to use a strict weight matrix which only gives a high score to\n"
"identities and the most favoured conservative substitutions. For more divergent\n"
"sequences, it is appropriate to use \"softer\" matrices which give a high score\n"
"to many other frequent substitutions.\n"
"\n"
"1) BLOSUM (Henikoff). These matrices appear to be the best available for \n"
"carrying out database similarity (homology searches). The matrices used are:\n"
"Blosum 80, 62, 45 and 30. (BLOSUM was the default in earlier Clustal W\n"
"versions)\n"
"\n"
"2) PAM (Dayhoff). These have been extremely widely used since the late '70s.\n"
"We use the PAM 20, 60, 120 and 350 matrices.\n"
"\n"
"3) GONNET. These matrices were derived using almost the same procedure as the\n"
"Dayhoff one (above) but are much more up to date and are based on a far larger\n"
"data set. They appear to be more sensitive than the Dayhoff series. We use the\n"
"GONNET 80, 120, 160, 250 and 350 matrices. This series is the default for\n"
"Clustal W version 1.8.\n"
"\n"
"We also supply an identity matrix which gives a score of 1.0 to two identical \n"
"amino acids and a score of zero otherwise. This matrix is not very useful.\n"
"Alternatively, you can read in your own (just one matrix, not a series).\n"
"\n"
"A new matrix can be read from a file on disk, if the filename consists only\n"
"of lower case characters. The values in the new weight matrix must be integers\n"
"and the scores should be similarities. You can use negative as well as positive\n"
"values if you wish, although the matrix will be automatically adjusted to all\n"
"positive scores.\n"
"\n"
"\n"
"\n"
"For DNA, a single matrix (not a series) is used. Two hard-coded matrices are \n"
"available:\n"
"\n"
"\n"
"1) IUB. This is the default scoring matrix used by BESTFIT for the comparison\n"
"of nucleic acid sequences. X's and N's are treated as matches to any IUB\n"
"ambiguity symbol. All matches score 1.9; all mismatches for IUB symbols score 0.\n"
" \n"
" \n"
"2) CLUSTALW(1.6). The previous system used by Clustal W, in which matches score\n"
"1.0 and mismatches score 0. All matches for IUB symbols also score 0.\n"
"\n"
"INPUT FORMAT  The format used for a new matrix is the same as the BLAST program.\n"
"Any lines beginning with a # character are assumed to be comments. The first\n"
"non-comment line should contain a list of amino acids in any order, using the\n"
"1 letter code, followed by a * character. This should be followed by a square\n"
"matrix of integer scores, with one row and one column for each amino acid. The\n"
"last row and column of the matrix (corresponding to the * character) contain\n"
"the minimum score over the whole matrix.\n"
"";
    sections.push_back(s);


    s.marker = "9";
    s.title = "Help for command line parameters";
    s.content = "                DATA (sequences)\n"
"\n"
"-INFILE=file.ext                             :input sequences.\n"
"-PROFILE1=file.ext  and  -PROFILE2=file.ext  :profiles (old alignment).\n"
"\n"
"\n"
"                VERBS (do things)\n"
"\n"
"-OPTIONS            :list the command line parameters\n"
"-HELP  or -CHECK    :outline the command line params.\n"
"-FULLHELP           :output full help content.\n"
"-ALIGN              :do full multiple alignment.\n"
"-TREE               :calculate NJ tree.\n"
"-PIM                :output percent identity matrix (while calculating the tree)\n"
"-BOOTSTRAP(=n)      :bootstrap a NJ tree (n= number of bootstraps; def. = 1000).\n"
"-CONVERT            :output the input sequences in a different file format.\n"
"\n"
"\n"
"                PARAMETERS (set things)\n"
"\n"
"***General settings:****\n"
"-INTERACTIVE :read command line, then enter normal interactive menus\n"
"-QUICKTREE   :use FAST algorithm for the alignment guide tree\n"
"-TYPE=       :PROTEIN or DNA sequences\n"
"-NEGATIVE    :protein alignment with negative values in matrix\n"
"-OUTFILE=    :sequence alignment file name\n"
"-OUTPUT=     :CLUSTAL(default), GCG, GDE, PHYLIP, PIR, NEXUS and FASTA\n"
"-OUTORDER=   :INPUT or ALIGNED\n"
"-CASE        :LOWER or UPPER (for GDE output only)\n"
"-SEQNOS=     :OFF or ON (for Clustal output only)\n"
"-SEQNO_RANGE=:OFF or ON (NEW: for all output formats)\n"
"-RANGE=m,n   :sequence range to write starting m to m+n\n"
"-MAXSEQLEN=n :maximum allowed input sequence length\n"
"-QUIET       :Reduce console output to minimum\n"        
"-STATS=      :Log some alignents statistics to file\n"
"\n"
"***Fast Pairwise Alignments:***\n"
"-KTUPLE=n    :word size\n"
"-TOPDIAGS=n  :number of best diags.\n"
"-WINDOW=n    :window around best diags.\n"
"-PAIRGAP=n   :gap penalty\n"
"-SCORE       :PERCENT or ABSOLUTE\n"
"\n"
"\n"
"***Slow Pairwise Alignments:***\n"
"-PWMATRIX=    :Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename\n"
"-PWDNAMATRIX= :DNA weight matrix=IUB, CLUSTALW or filename\n"
"-PWGAPOPEN=f  :gap opening penalty        \n"
"-PWGAPEXT=f   :gap opening penalty\n"
"\n"
"\n"
"***Multiple Alignments:***\n"
"-NEWTREE=      :file for new guide tree\n"
"-USETREE=      :file for old guide tree\n"
"-MATRIX=       :Protein weight matrix=BLOSUM, PAM, GONNET, ID or filename\n"
"-DNAMATRIX=    :DNA weight matrix=IUB, CLUSTALW or filename\n"
"-GAPOPEN=f     :gap opening penalty        \n"
"-GAPEXT=f      :gap extension penalty\n"
"-ENDGAPS       :no end gap separation pen. \n"
"-GAPDIST=n     :gap separation pen. range\n"
"-NOPGAP        :residue-specific gaps off  \n"
"-NOHGAP        :hydrophilic gaps off\n"
"-HGAPRESIDUES= :list hydrophilic res.    \n"
"-MAXDIV=n      :% ident. for delay\n"
"-TYPE=         :PROTEIN or DNA\n"
"-TRANSWEIGHT=f :transitions weighting\n"
"-ITERATION=    :NONE or TREE or ALIGNMENT\n"
"-NUMITER=n     :maximum number of iterations to perform\n"
"-NOWEIGHTS     :disable sequence weighting\n"
"\n"
"\n"
"***Profile Alignments:***\n"
"-PROFILE      :Merge two alignments by profile alignment\n"
"-NEWTREE1=    :file for new guide tree for profile1\n"
"-NEWTREE2=    :file for new guide tree for profile2\n"
"-USETREE1=    :file for old guide tree for profile1\n"
"-USETREE2=    :file for old guide tree for profile2\n"
"\n"
"\n"
"***Sequence to Profile Alignments:***\n"
"-SEQUENCES   :Sequentially add profile2 sequences to profile1 alignment\n"
"-NEWTREE=    :file for new guide tree\n"
"-USETREE=    :file for old guide tree\n"
"\n"
"\n"
"***Structure Alignments:***\n"
"-NOSECSTR1     :do not use secondary structure-gap penalty mask for profile 1 \n"
"-NOSECSTR2     :do not use secondary structure-gap penalty mask for profile 2\n"
"-SECSTROUT=STRUCTURE or MASK or BOTH or NONE   :output in alignment file\n"
"-HELIXGAP=n    :gap penalty for helix core residues \n"
"-STRANDGAP=n   :gap penalty for strand core residues\n"
"-LOOPGAP=n     :gap penalty for loop regions\n"
"-TERMINALGAP=n :gap penalty for structure termini\n"
"-HELIXENDIN=n  :number of residues inside helix to be treated as terminal\n"
"-HELIXENDOUT=n :number of residues outside helix to be treated as terminal\n"
"-STRANDENDIN=n :number of residues inside strand to be treated as terminal\n"
"-STRANDENDOUT=n:number of residues outside strand to be treated as terminal \n"
"\n"
"\n"
"***Trees:***\n"
"-OUTPUTTREE=nj OR phylip OR dist OR nexus\n"
"-SEED=n        :seed number for bootstraps.\n"
"-KIMURA        :use Kimura's correction.   \n"
"-TOSSGAPS      :ignore positions with gaps.\n"
"-BOOTLABELS=node OR branch :position of bootstrap values in tree display\n"
"-CLUSTERING=   :NJ or UPGMA\n"
"";
    sections.push_back(s);


    s.marker = "0";
    s.title = "Help for tree output format options";
    s.content = "Four output formats are offered: 1) Clustal, 2) Phylip, 3) Just the distances\n"
"4) Nexus\n"
"\n"
"None of these formats displays the results graphically. Many packages can\n"
"display trees in the the PHYLIP format 2) below. It can also be imported into\n"
"the PHYLIP programs RETREE, DRAWTREE and DRAWGRAM for graphical display. \n"
"NEXUS format trees can be read by PAUP and MacClade.\n"
"\n"
"1) Clustal format output. \n"
"This format is verbose and lists all of the distances between the sequences and\n"
"the number of alignment positions used for each. The tree is described at the\n"
"end of the file. It lists the sequences that are joined at each alignment step\n"
"and the branch lengths. After two sequences are joined, it is referred to later\n"
"as a NODE. The number of a NODE is the number of the lowest sequence in that\n"
"NODE.   \n"
"\n"
"2) Phylip format output.\n"
"This format is the New Hampshire format, used by many phylogenetic analysis\n"
"packages. It consists of a series of nested parentheses, describing the\n"
"branching order, with the sequence names and branch lengths. It can be used by\n"
"the RETREE, DRAWGRAM and DRAWTREE programs of the PHYLIP package to see the\n"
"trees graphically. This is the same format used during multiple alignment for\n"
"the guide trees. \n"
"\n"
"Use this format with NJplot (Manolo Gouy), supplied with Clustal W. Some other\n"
"packages that can read and display New Hampshire format are TreeView (Mac/PC),\n"
"TreeTool (UNIX), and Phylowin.\n"
"\n"
"3) The distances only.\n"
"This format just outputs a matrix of all the pairwise distances in a format\n"
"that can be used by the Phylip package. It used to be useful when one could not\n"
"produce distances from protein sequences in the Phylip package but is now\n"
"redundant (Protdist of Phylip 3.5 now does this).\n"
"\n"
"4) NEXUS FORMAT TREE. This format is used by several popular phylogeny programs,\n"
"including PAUP and MacClade. The format is described fully in:\n"
"Maddison, D. R., D. L. Swofford and W. P. Maddison.  1997.\n"
"NEXUS: an extensible file format for systematic information.\n"
"Systematic Biology 46:590-621.\n"
"\n"
"5) TOGGLE PHYLIP BOOTSTRAP POSITIONS\n"
"By default, the bootstrap values are placed on the nodes of the phylip format\n"
"output tree. This is inaccurate as the bootstrap values should be associated\n"
"with the tree branches and not the nodes. However, this format can be read and\n"
"displayed by TreeTool, TreeView and Phylowin. An option is available to\n"
"correctly place the bootstrap values on the branches with which they are\n"
"associated.\n"
"";
    sections.push_back(s);


    // mark
    // replace-string " \"
    // replace-regexp ^ "
    // replace-regexp $ \\n"
    
    // std::cout << "exiting Help::Help" << "\n";
}
Help::~Help()
{
    sections.clear();
}


vector<string> Help::ListSectionMarkers()
{
    vector<string> markers;
    for (unsigned int i=0; i<this->sections.size(); i++) {
        markers.push_back(this->sections[i].marker);        
    }
    
    return markers;
}


string Help::GetSection(char marker)
{
    string s(1, marker);
    return GetSection(s);
}

string Help::GetSection(string marker)
{
    for (unsigned int i=0; i<this->sections.size(); i++) {
        if (this->sections[i].marker == marker) {
            return this->sections[i].content;
        }
    }   
    return "";
}

string Help::GetSectionTitle(char marker)
{
    string s(1, marker);
    return GetSectionTitle(s);
}
string Help::GetSectionTitle(string marker)
{
    for (unsigned int i=0; i<this->sections.size(); i++) {
        if (this->sections[i].marker == marker) {
            return this->sections[i].title;
        }
    }   
    return "";
}
