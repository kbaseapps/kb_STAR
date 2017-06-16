# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
import json
import time
from pprint import pprint, pformat

from STAR.Utils.STARUtils import STARUtil
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
    GIT_COMMIT_HASH = "d9f109adadd4d58e98c605b6b99055a41837034e"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block

    def log(self, message, prefix_newline=False):
            print(('\n' if prefix_newline else '') +
                str(time.time()) + ': ' + str(message))

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        
        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.config = config
        self.workspaceURL = config['workspace-url']
        self.scratch = os.path.abspath(config['scratch'])
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        self.shared_folder = config['scratch']

        #END_CONSTRUCTOR
        pass


    def run_star(self, ctx, params):
        """
        The actual function is declared using 'funcdef' to specify the name
        and input/return arguments to the function.  For all typical KBase
        Apps that run in the Narrative, your function should have the 
        'authentication required' modifier.
        :param params: instance of type "STARParams" (Arguments for
           star_generate_indexes obj_ref reads_ref: KBase workspace reference
           for reads to align obj_ref genome_ref: KBase workspace reference
           genome to align reads against string workspace_name: the workspace
           name provided by the narrative for housing output in KBase string
           runMode: default: alignReads type of the run: alignReads => map
           reads genomeGenerate => generate genome files
           inputAlignmentsFromBAM => input alignments from BAM. Presently
           only works with -outWigType and -bamRemoveDuplicates. liftOver =>
           lift-over of GTF files (-sjdbGTFfile) between genome assemblies
           using chain file(s) from -genomeChainFiles int runThreadN:
           default: 1 number of threads to run STAR list<obj_ref>
           genomeFastaFiles: path(s) to the fasta files with genomic
           sequences for genome generation Only used if
           runMode==genomeGenerate. These files should be plain text FASTA
           files, they *cannot* be zipped. list<obj_ref> readFilesIn:
           default: Read1 Read2 paths to files that contain input read1 (and,
           if needed, read2) string sjdbGTFfile: default: -; path to the file
           with annotated transcripts in the standard GTF format int
           sjdbOverhang: default: 100; int>0: length of the donor/acceptor
           sequence on each side of the junctions, ideally = (ReadLength - 1)
           string outFileNamePrefix: you can change the file prefixes using
           --outFileNamePrefix /path/to/output/dir/prefix By default, this
           parameter is ./, i.e. all output files are written in the current
           directory @optional sjdbGTFfile @optional sjdbOverhang) ->
           structure: parameter "reads_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "genome_ref" of type "obj_ref" (An X/Y/Z
           style reference), parameter "workspace_name" of String, parameter
           "runMode" of String, parameter "runThreadN" of Long, parameter
           "genomeFastaFiles" of list of type "obj_ref" (An X/Y/Z style
           reference), parameter "sjdbGTFfile" of String, parameter
           "sjdbOverhang" of Long, parameter "readFilesIn" of list of type
           "obj_ref" (An X/Y/Z style reference), parameter
           "outFileNamePrefix" of String
        :returns: instance of type "STARResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report. output_folder: folder path that holds all files generated
           by STAT report_name: report name generated by KBaseReport
           report_ref: report reference generated by KBaseReport) ->
           structure: parameter "output_folder" of String, parameter
           "alignment_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_star
        self.log('Running run_star with params:\n' + pformat(params))

        token = ctx['token']

        for key, value in params.iteritems():
            if isinstance(value, basestring):
                params[key] = value.strip()

        star_runner = STARUtil(self.config)
        output = star_runner.run_star(params)
        #END run_star

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_star return value ' +
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
