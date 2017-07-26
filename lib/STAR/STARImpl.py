# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
import json
import time
from STAR.Utils.file_util import (
    fetch_reads_from_reference,
    fetch_reads_refs_from_sampleset
)
from KBParallel.KBParallelClient import KBParallel
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
    GIT_COMMIT_HASH = "c77a101595cf73a5a9ecd3f0dd0a7765682dd07e"

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
            self.workspace_url = config['workspace-url']
        if 'shock-url' in config:
            self.__SHOCK_URL = config['shock-url']
        if 'handle-service-url' in config:
            self.__HS_URL = config['handle-service-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.srv_wiz_url = config['srv-wiz-url']
        self.__CALLBACK_URL = os.environ['SDK_CALLBACK_URL']

        self.__SERVICES = {'workspace_service_url': self.workspace_url,
                           'shock_service_url': self.__SHOCK_URL,
                           'handle_service_url': self.__HS_URL,
                           'callback_url': self.callback_url}

        self.scratch = os.path.abspath(config['scratch'])
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        self.shared_folder = config['scratch']

        #END_CONSTRUCTOR
        pass


    def star_align_reads_to_assembly(self, ctx, params):
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
           between mates, default to 1000000 create_report = 1 if we build a
           report, 0 otherwise. (default 1) (shouldn not be user set - mainly
           used for subtasks) @optional output_alignment_filename_extension
           @optional outFilterType @optional outFilterMultimapNmax @optional
           outSAMtype @optional outSAMattrIHstart @optional outSAMstrandField
           @optional quantMode @optional alignSJoverhangMin @optional
           alignSJDBoverhangMin @optional outFilterMismatchNmax @optional
           alignIntronMin @optional alignIntronMax @optional alignMatesGapMax
           @optional outFileNamePrefix) -> structure: parameter
           "readsset_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "genome_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "output_workspace" of String, parameter
           "output_name" of String, parameter
           "output_alignment_filename_extension" of String, parameter
           "runThreadN" of Long, parameter "condition" of String, parameter
           "outFilterType" of String, parameter "outSAMtype" of String,
           parameter "outSAMattrIHstart" of Long, parameter
           "outSAMstrandField" of String, parameter "quantMode" of String,
           parameter "outFilterMultimapNmax" of Long, parameter
           "alignSJoverhangMin" of Long, parameter "alignSJDBoverhangMin" of
           Long, parameter "outFilterMismatchNmax" of Long, parameter
           "alignIntronMin" of Long, parameter "alignIntronMax" of Long,
           parameter "alignMatesGapMax" of Long, parameter
           "outFileNamePrefix" of String, parameter "concurrent_njsw_tasks"
           of Long, parameter "concurrent_local_tasks" of Long, parameter
           "create_report" of type "bool" (A boolean - 0 for false, 1 for
           true. @range (0, 1))
        :returns: instance of type "AlignReadsResult" (Here is the definition
           of the output of the function.  The output can be used by other
           SDK modules which call your code, or the output visualizations in
           the Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report. alignment_ref: can be either an Alignment or AlignmentSet,
           depending on inputs. report_name: report name generated by
           KBaseReport report_ref: report reference generated by KBaseReport)
           -> structure: parameter "alignment_ref" of type "obj_ref" (An
           X/Y/Z style reference), parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN star_align_reads_to_assembly
        star_runner = STARUtil(self.config)

        result = {}
        # pre-process the input parameters
        validated_params = star_runner.process_params(params)
        result = star_runner.star_run_single(validated_params)

        # indexing if not yet existing
        if not hasattr(self.__class__, 'STARGenomeIndex'):
            # convert the input parameters (from refs to file paths, especially)
            params_ret = star_runner.convert_params(validated_params)
            input_params = params_ret.get('input_parameters', None)
            # generate the indices
            idx_ret = star_runner.run_star_indexing(input_params)
            if idx_ret == 0:
                self.__class__.STARGenomeIndex = 'indexed'
        #END star_align_reads_to_assembly

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method star_align_reads_to_assembly return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]

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
           between mates, default to 1000000 create_report = 1 if we build a
           report, 0 otherwise. (default 1) (shouldn not be user set - mainly
           used for subtasks) @optional output_alignment_filename_extension
           @optional outFilterType @optional outFilterMultimapNmax @optional
           outSAMtype @optional outSAMattrIHstart @optional outSAMstrandField
           @optional quantMode @optional alignSJoverhangMin @optional
           alignSJDBoverhangMin @optional outFilterMismatchNmax @optional
           alignIntronMin @optional alignIntronMax @optional alignMatesGapMax
           @optional outFileNamePrefix) -> structure: parameter
           "readsset_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "genome_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "output_workspace" of String, parameter
           "output_name" of String, parameter
           "output_alignment_filename_extension" of String, parameter
           "runThreadN" of Long, parameter "condition" of String, parameter
           "outFilterType" of String, parameter "outSAMtype" of String,
           parameter "outSAMattrIHstart" of Long, parameter
           "outSAMstrandField" of String, parameter "quantMode" of String,
           parameter "outFilterMultimapNmax" of Long, parameter
           "alignSJoverhangMin" of Long, parameter "alignSJDBoverhangMin" of
           Long, parameter "outFilterMismatchNmax" of Long, parameter
           "alignIntronMin" of Long, parameter "alignIntronMax" of Long,
           parameter "alignMatesGapMax" of Long, parameter
           "outFileNamePrefix" of String, parameter "concurrent_njsw_tasks"
           of Long, parameter "concurrent_local_tasks" of Long, parameter
           "create_report" of type "bool" (A boolean - 0 for false, 1 for
           true. @range (0, 1))
        :returns: instance of type "AlignReadsResult" (Here is the definition
           of the output of the function.  The output can be used by other
           SDK modules which call your code, or the output visualizations in
           the Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report. alignment_ref: can be either an Alignment or AlignmentSet,
           depending on inputs. report_name: report name generated by
           KBaseReport report_ref: report reference generated by KBaseReport)
           -> structure: parameter "alignment_ref" of type "obj_ref" (An
           X/Y/Z style reference), parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_star
        self.log('Running run_star with params:\n' + pformat(params))
        for key, value in params.iteritems():
            if isinstance(value, basestring):
                params[key] = value.strip()

        returnVal = {}

        star_runner = STARUtil(self.config)

        # 1. process the input parameters
        validated_params = star_runner.process_params(params)
        input_obj_info = star_runner.determine_input_info(validated_params)

        # indexing if not yet existing
        if not hasattr(self.__class__, 'STARGenomeIndex'):
            # convert the input parameters (from refs to file paths, especially)
            params_ret = star_runner.convert_params(validated_params)
            input_params = params_ret.get('input_parameters', None)
            # generate the indices
            idx_ret = star_runner.run_star_indexing(input_params)
            if idx_ret == 0:
                self.__class__.STARGenomeIndex = 'indexed'

        # 2. Run STAR with index and reads.
        # If there's only one, run it locally right now.
        # If there's more than one:
        #  1. make a list of tasks to send to KBParallel.
        #  2. add a flag to not make a report for each subtask.
        #  3. make the report when it's all done.
        if input_obj_info['run_mode'] == 'single_library':
            star_ret = star_runner.star_run_single(validated_params)
        elif input_obj_info['run_mode'] == 'sample_set':
            star_ret = star_runner.star_run_batch(validated_params)

        if star_ret.get('alignment_ref', None) is not None:
            returnVal['alignment_ref'] = star_ret['alignment_ref']
        if star_ret['report_info'].get('name', None) is not None:
            returnVal['report_name'] = star_ret['report_info']['name']
            returnVal['report_ref'] = star_ret['report_info']['ref']
        #END run_star

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_star return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]

    def get_star_index(self, ctx, params):
        """
        :param params: instance of type "GetSTARIndex" (Provide a reference
           to either an Assembly or Genome to get a STAR index. output_dir is
           optional, if provided the index files will be saved in that
           directory.  If not, a directory will be generated for you and
           returned by this function.  If specifying the output_dir, the
           directory must not exist yet (to ensure only the index files are
           added there). Currently, STAR indexes are cached to a WS object. 
           If that object does not exist, then calling this function can
           create a new object.  To create the cache, you must specify the ws
           name or ID in 'ws_for_cache' in which to create the cached index. 
           If this field is not set, the result will not be cached.  This
           parameter will eventually be deprecated once the big file cache
           service is implemented.) -> structure: parameter "ref" of String,
           parameter "output_dir" of String, parameter "ws_for_cache" of
           String
        :returns: instance of type "GetSTARIndexResult" (output_dir - the
           folder containing the index files from_cache - 0 if the index was
           built fresh, 1 if it was found in the cache pushed_to_cache - if
           the index was rebuilt and successfully added to the cache, this
           will be set to 1; otherwise set to 0) -> structure: parameter
           "output_dir" of String, parameter "from_cache" of type "bool" (A
           boolean - 0 for false, 1 for true. @range (0, 1)), parameter
           "pushed_to_cache" of type "bool" (A boolean - 0 for false, 1 for
           true. @range (0, 1))
        """
        # ctx is the context object
        # return variables are: result
        #BEGIN get_star_index
        #END get_star_index

        # At some point might do deeper type checking...
        if not isinstance(result, dict):
            raise ValueError('Method get_star_index return value ' +
                             'result is not type dict as required.')
        # return the results
        return [result]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
