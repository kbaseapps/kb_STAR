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
    GIT_COMMIT_HASH = "00cef08012550b3fc2301a560b058630d94de38d"

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
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']

        if 'workspace-url' in config:
            self.workspaceURL = config['workspace-url']
        if 'shock-url' in config:
            self.__SHOCK_URL = config['shock-url']
        if 'handle-service-url' in config:
            self.__HS_URL = config['handle-service-url']
        self.__CALLBACK_URL = os.environ['SDK_CALLBACK_URL']

        self.__SERVICES = {'workspace_service_url': self.workspaceURL,
                           'shock_service_url': self.__SHOCK_URL,
                           'handle_service_url': self.__HS_URL,
                           'callback_url': self.__CALLBACK_URL}

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
        :param params: instance of type "STARParams" (Will align the input
           reads (or set of reads specified in a SampleSet) to the specified
           assembly or assembly for the specified Genome (accepts Assembly,
           ContigSet, or Genome types) and produces a ReadsAlignment object,
           or in the case of a SampleSet, a ReadsAlignmentSet object obj_ref
           genome_ref: KBase workspace reference genome to align reads
           against obj_ref sampleset_ref: the workspace reference for the
           sampleset of reads to align, referring to either a
           SingleEnd/PairedEnd reads, or a SampleSet, or a ReadsSet input
           string workspace_name: the workspace name provided by the
           narrative for housing output in KBase string alignmentset_name -
           name of the output ReadsAlignment or ReadsAlignmentSet object int
           runThreadN - the number of threads for STAR to use (default to 2)
           string outFileNamePrefix: you can change the file prefixes using
           --outFileNamePrefix /path/to/output/dir/prefix By default, this
           parameter is ./, i.e. all output files are written in current
           directory without a prefix string quantMode: types of
           quantification requested--none/TranscriptomeSAM/GeneCounts int
           outFilterMultimapNmax: max number of multiple alignments allowed
           for a read: if exceeded, the read is considered unmapped, default
           to 20 int alignSJoverhangMin: minimum overhang for unannotated
           junctions, default to 8 int alignSJDBoverhangMin: minimum overhang
           for annotated junctions, default to 1 int outFilterMismatchNmax:
           maximum number of mismatches per pair, large number switches off
           this filter, default to 999 int alignIntronMin: minimum intron
           length, default to 20 int alignIntronMax: maximum intron length,
           default to 1000000 int alignMatesGapMax: maximum genomic distance
           between mates, default to 1000000 @optional outFilterType
           @optional outFilterMultimapNmax @optional outSAMtype @optional
           outSAMattrIHstart @optional outSAMstrandField @optional quantMode
           @optional alignSJoverhangMin @optional alignSJDBoverhangMin
           @optional outFilterMismatchNmax @optional alignIntronMin @optional
           alignIntronMax @optional alignMatesGapMax @optional
           outFileNamePrefix) -> structure: parameter "sampleset_ref" of type
           "obj_ref" (An X/Y/Z style reference), parameter "genome_ref" of
           type "obj_ref" (An X/Y/Z style reference), parameter
           "workspace_name" of String, parameter "runThreadN" of Long,
           parameter "alignmentset_name" of String, parameter "output_name"
           of String, parameter "outFilterType" of String, parameter
           "outSAMtype" of String, parameter "outSAMattrIHstart" of Long,
           parameter "outSAMstrandField" of String, parameter "quantMode" of
           String, parameter "outFilterMultimapNmax" of Long, parameter
           "alignSJoverhangMin" of Long, parameter "alignSJDBoverhangMin" of
           Long, parameter "outFilterMismatchNmax" of Long, parameter
           "alignIntronMin" of Long, parameter "alignIntronMax" of Long,
           parameter "alignMatesGapMax" of Long, parameter
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
