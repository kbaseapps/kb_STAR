/*
   Name of module: STAR

   This KBase module wraps the free open source software STAR: ultrafast universal RNA-seq aligner.
   STAR-2.5.3a

   References:
   https://github.com/alexdobin/STAR/
   https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
*/

module STAR {
    /* 
        A 'typedef' allows you to provide a more specific name for
        a type.  Built-in primitive types include 'string', 'int',
        'float'.  Here we define a type named assembly_ref to indicate
        a string that should be set to a KBase ID reference to an
        Assembly data object.
    */
    typedef string assembly_ref;

    /*
	Arguments for star_generate_indexes

	string runMode: default: alignReads
		type of the run:
		alignReads => map reads
		genomeGenerate => generate genome files
		inputAlignmentsFromBAM => input alignments from BAM. Presently only works with -outWigType
			and -bamRemoveDuplicates.
		liftOver => lift-over of GTF files (-sjdbGTFfile) between genome assemblies using
			chain file(s) from -genomeChainFiles.
	int runThreadN: default: 1
		number of threads to run STAR
	list<string> genomeFastaFiles: path(s) to the fasta files with genomic sequences for genome generation. 
		Only used if runMode==genomeGenerate.These files should be plain text FASTA files, they *cannot* be zipped.
	string sjdbGTFfile: default: -; path to the GTF file with annotations
	int sjdbOverhang: default: 100; int>0: length of the donor/acceptor sequence on each side of the junctions,
		ideally = (mate length - 1)
	list<string> readFilesIn: default: Read1 Read2
		paths to files that contain input read1 (and, if needed, read2)
    */
    typedef structure {
        string reads_ref;
        string assembly_ref;
        string genome_ref;

        string workspace_name;
	string runMode;
	int runThreadN;
	list<string> genomeFastaFiles;
	string sjdbGTFfile;
	int sjdbOverhang;
	
	list<string> readFilesIn;
        string outFileNamePrefix;
    } AlignReadsParams;

    /*
        Here is the definition of the output of the function.  The output
        can be used by other SDK modules which call your code, or the output
        visualizations in the Narrative.  'report_name' and 'report_ref' are
        special output fields- if defined, the Narrative can automatically
        render your Report.
    */
    typedef structure {
	string reads_alignment_ref;

        string report_name;
        string report_ref;
    } AlignReadsResults;
    
    /*
        The actual function is declared using 'funcdef' to specify the name
        and input/return arguments to the function.  For all typical KBase
        Apps that run in the Narrative, your function should have the 
        'authentication required' modifier.
    */
    funcdef run_star(AlignReadsParams params)
        returns (AlignReadsResults output) authentication required;

};
