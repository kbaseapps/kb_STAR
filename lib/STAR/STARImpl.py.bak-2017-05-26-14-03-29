# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
from Bio import SeqIO
from pprint import pprint, pformat
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport
#END_HEADER


class STAR:
    '''
    Module Name:
    STAR

    Module Description:
    Name of module: STAR

This KBase module wraps the free open source software STAR: ultrafast universal RNA-seq aligner.
STAR-2.5.3a

References:
https://github.com/alexdobin/STAR/
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_STAR.git"
    GIT_COMMIT_HASH = "3867ac434fd223aecda32543fbe1ec3475d484c9"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        
        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']

        #END_CONSTRUCTOR
        pass


    def star_generate_indexes(self, ctx, params):
        """
        The actual function is declared using 'funcdef' to specify the name
        and input/return arguments to the function.  For all typical KBase
        Apps that run in the Narrative, your function should have the 
        'authentication required' modifier.
        :param params: instance of type "GenerateIndexesParams" (Arguments
           for star_generate_indexes string runMode: default: alignReads type
           of the run: alignReads => map reads genomeGenerate => generate
           genome files inputAlignmentsFromBAM => input alignments from BAM.
           Presently only works with -outWigType and -bamRemoveDuplicates.
           liftOver => lift-over of GTF files (-sjdbGTFfile) between genome
           assemblies using chain file(s) from -genomeChainFiles. int
           runThreadN: default: 1 number of threads to run STAR list<string>
           genomeFastaFiles: path(s) to the fasta files with genomic
           sequences for genome generation. Only used if
           runMode==genomeGenerate.These files should be plain text FASTA
           files, they *cannot* be zipped. string sjdbGTFfile: default: -;
           path to the GTF file with annotations int sjdbOverhang: default:
           100; int>0: length of the donor/acceptor sequence on each side of
           the junctions, ideally = (mate length - 1)) -> structure:
           parameter "workspace_name" of String, parameter "runMode" of
           String, parameter "runThreadN" of Long, parameter
           "genomeFastaFiles" of list of String, parameter "sjdbGTFfile" of
           String, parameter "sjdbOverhang" of Long
        :returns: instance of type "STARResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN star_generate_indexes
        #END star_generate_indexes

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method star_generate_indexes return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def star_mapping(self, ctx, params):
        """
        The actual function is declared using 'funcdef' to specify the name
        and input/return arguments to the function.  For all typical KBase
        Apps that run in the Narrative, your function should have the 
        'authentication required' modifier.
        :param params: instance of type "MappingParams" (Arguments for
           star_mapping int runThreadN: default: 1 number of threads to run
           STAR list<string> readFilesIn: default: Read1 Read2 paths to files
           that contain input read1 (and, if needed, read2)) -> structure:
           parameter "workspace_name" of String, parameter "runThreadN" of
           Long, parameter "readFilesIn" of list of String
        :returns: instance of type "STARResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN star_mapping
        #END star_mapping

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method star_mapping return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
