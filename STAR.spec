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

    /* A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int boolean;

    /*
	Arguments for star_generate_indexes

	string reads_ref, assembly_ref and genome_ref: KBase style variable references
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
	list<string> readFilesIn: default: Read1 Read2
		paths to files that contain input read1 (and, if needed, read2)

	string sjdbGTFfile: default: -; path to the file with annotated transcripts in the standard GTF format
	int sjdbOverhang: default: 100; int>0: length of the donor/acceptor sequence on each side of the junctions,
		ideally = (ReadLength - 1)
	string outFileNamePrefix: you can change the file prefixes using --outFileNamePrefix /path/to/output/dir/prefix.
		By default, this parameter is ./, i.e. all output files are written in the current directory

	@optional sjdbGTFfile
	@optional sjdbOverhang
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
    } STARParams;

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
    } STARResults;
    
    /*
        The actual function is declared using 'funcdef' to specify the name
        and input/return arguments to the function.  For all typical KBase
        Apps that run in the Narrative, your function should have the 
        'authentication required' modifier.
    */
    funcdef run_star(STARParams params)
        returns (STARResults output) authentication required;

};
