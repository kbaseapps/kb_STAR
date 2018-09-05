# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
import time
from pprint import pformat

from STAR.Utils.STAR_Aligner import STAR_Aligner
#END_HEADER


class STAR:
    '''
    Module Name:
    STAR

    Module Description:
    Name of module: STAR

This KBase module wraps the free open source software STAR: ultrafast universal RNA-seq aligner.
STAR-2.6.1a

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
    VERSION = "0.1.0"
    GIT_URL = "https://github.com/kbaseapps/kb_STAR.git"
    GIT_COMMIT_HASH = "7c1524911d95ef57d241dc223e096949af892d06"

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
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        #END_CONSTRUCTOR
        pass


    def run_star(self, ctx, params):
        """
        :param params: instance of type "AlignReadsParams" (Will align the
           input reads (or set of reads specified in a SampleSet) to the
           specified assembly or assembly for the specified Genome (accepts
           Assembly, ContigSet, or Genome types) and produces a
           ReadsAlignment object, or in the case of a SampleSet, a
           ReadsAlignmentSet object obj_ref genome_ref: KBase workspace
           reference Genome obj_ref readsset_ref: the workspace reference for
           the set of reads to align, referring to either a
           SingleEnd/PairedEnd reads, or a ReadsSet input string
           output_workspace - name or id of the WS to save the results to,
           provided by the narrative for housing output in KBase string
           output_name - name of the output ReadsAlignment or
           ReadsAlignmentSet object int runThreadN - the number of threads
           for STAR to use (default to 2) string outFileNamePrefix: you can
           change the file prefixes using --outFileNamePrefix
           /path/to/output/dir/prefix By default, this parameter is ./, i.e.
           all output files are written in current directory without a prefix
           string quantMode: types of quantification
           requested--none/TranscriptomeSAM/GeneCounts int
           outFilterMultimapNmax: max number of multiple alignments allowed
           for a read: if exceeded, the read is considered unmapped, default
           to 20 int alignSJoverhangMin: minimum overhang for unannotated
           junctions, default to 8 int alignSJDBoverhangMin: minimum overhang
           for annotated junctions, default to 1 int outFilterMismatchNmax:
           maximum number of mismatches per pair, large number switches off
           this filter, default to 999 int alignIntronMin: minimum intron
           length, default to 20 int alignIntronMax: maximum intron length,
           default to 1000000 int alignMatesGapMax: maximum genomic distance
           between mates, default to 1000000 int create_report: = 1 if we
           build a report, 0 otherwise. (default 1) (shouldn not be user set
           - mainly used for subtasks) @optional alignmentset_suffix
           @optional alignIntronMin @optional alignIntronMax @optional
           alignMatesGapMax @optional alignSJoverhangMin @optional
           alignSJDBoverhangMin @optional quantMode @optional outFilterType
           @optional outFilterMultimapNmax @optional outSAMtype @optional
           outSAMattrIHstart @optional outSAMstrandField @optional
           outFilterMismatchNmax @optional outFileNamePrefix) -> structure:
           parameter "readsset_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "genome_ref" of type "obj_ref" (An X/Y/Z
           style reference), parameter "output_workspace" of String,
           parameter "output_name" of String, parameter "alignment_suffix" of
           String, parameter "condition" of String, parameter
           "concurrent_njsw_tasks" of Long, parameter
           "concurrent_local_tasks" of Long, parameter "outSAMunmapped" of
           String, parameter "create_report" of type "bool" (A boolean - 0
           for false, 1 for true. @range (0, 1)), parameter
           "alignmentset_suffix" of String, parameter "alignIntronMin" of
           Long, parameter "alignIntronMax" of Long, parameter
           "alignMatesGapMax" of Long, parameter "alignSJoverhangMin" of
           Long, parameter "alignSJDBoverhangMin" of Long, parameter
           "quantMode" of String, parameter "outFilterType" of String,
           parameter "outFilterMultimapNmax" of Long, parameter "outSAMtype"
           of String, parameter "outSAMattrIHstart" of Long, parameter
           "outSAMstrandField" of String, parameter "outFilterMismatchNmax"
           of Long, parameter "outFileNamePrefix" of String
        :returns: instance of type "AlignReadsResult" (Here is the definition
           of the output of the function.  The output can be used by other
           SDK modules which call your code, or the output visualizations in
           the Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report. output_directory: folder path that holds all output files
           generated by run_star alignmentset_ref: if an alignment set is
           created alignment_objs: for each individual alignment created. The
           keys are the references to the reads object being aligned.
           report_name: report name generated by KBaseReport report_ref:
           report reference generated by KBaseReport) -> structure: parameter
           "output_directory" of String, parameter "report_name" of String,
           parameter "report_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "alignmentset_ref" of type "obj_ref" (An
           X/Y/Z style reference), parameter "alignment_objs" of mapping from
           type "obj_ref" (An X/Y/Z style reference) to type "AlignmentObj"
           (Created alignment object returned. ref = the workspace reference
           of the new alignment object name = the name of the new object, for
           convenience.) -> structure: parameter "ref" of type "obj_ref" (An
           X/Y/Z style reference), parameter "name" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_star
        self.log('Running run_star with params:\n' + pformat(params))
        for key, value in params.iteritems():
            if isinstance(value, basestring):
                params[key] = value.strip()

        star_aligner = STAR_Aligner(self.config, ctx.provenance())

        returnVal = star_aligner.run_align(params)
        #END run_star

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_star return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
