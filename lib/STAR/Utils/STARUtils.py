import time
import json
import os
import re
import copy
import uuid
import subprocess
import shutil
import sys
import zipfile
from pprint import pprint, pformat
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing

from STAR.Utils.Program_Runner import Program_Runner
from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from KBParallel.KBParallelClient import KBParallel
from kb_QualiMap.kb_QualiMapClient import kb_QualiMap
from SetAPI.SetAPIServiceClient import SetAPI
from ExpressionUtils. ExpressionUtilsClient import ExpressionUtils

from file_util import (
    valid_string,
    check_reference,
    get_unique_names,
    fetch_fasta_from_object,
    fetch_reads_refs_from_sampleset,
    fetch_reads_from_reference,
    extract_geneCount_matrix
)

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class STARUtils:
    STAR_VERSION = 'STAR 2.5.3a'
    STAR_BIN = '/kb/deployment/bin/STAR'
    STAR_IDX_DIR = 'STAR_Genome_index'
    STAR_OUT_DIR = 'STAR_Output'
    PARAM_IN_WS = 'output_workspace'
    PARAM_IN_FASTA_FILES = 'genomeFastaFiles'
    PARAM_IN_OUTFILE_PREFIX = 'outFileNamePrefix'
    PARAM_IN_STARMODE = 'runMode'
    PARAM_IN_THREADN = 'runThreadN'
    PARAM_IN_READS_FILES = 'readFilesIn'
    PARAM_IN_READS = 'readsset_ref'
    PARAM_IN_GENOME = 'genome_ref'
    SET_READS = 'set_reads_refs'

    def __init__(self, scratch_dir, workspace_url, callback_url, srv_wiz_url, provenance):
        self.workspace_url = workspace_url
        self.callback_url = callback_url
        self.srv_wiz_url = srv_wiz_url
        self.au = AssemblyUtil(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.scratch = scratch_dir
        self.working_dir = scratch_dir
        self.prog_runner = Program_Runner(self.STAR_BIN, self.scratch)
        self.provenance = provenance
        self.ws_client = Workspace(self.workspace_url)

        self.parallel_runner = KBParallel(self.callback_url)
        self.qualimap = kb_QualiMap(self.callback_url, service_ver='dev')
        self.set_api_client = SetAPI(self.srv_wiz_url, service_ver='dev')
        self.eu = ExpressionUtils(self.callback_url, service_ver='dev')

    def process_params(self, params):
        """
        process_params: checks params passed to run_star method and set default values
        """
        log('Start validating run_star parameters')
        # check for required parameters
        if params.get(self.PARAM_IN_WS, None) is None:
            raise ValueError(self.PARAM_IN_WS + ' parameter is required')
        if (params.get(self.PARAM_IN_READS, None) is None or
                not valid_string(params[self.PARAM_IN_READS], is_ref=True)):
            raise ValueError("Parameter readsset_ref must be a valid Workspace object reference, "
                      "not {}".format(params.get(self.PARAM_IN_READS, None)))
        if "alignment_suffix" not in params or not valid_string(params["alignment_suffix"]):
            raise ValueError("Parameter alignment_suffix must be a valid Workspace object string, "
                      "not {}".format(params.get("alignment_suffix", None)))

        if params.get(self.PARAM_IN_STARMODE, None) is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'
	if params.get(self.PARAM_IN_GENOME, None) is None:
            raise ValueError(self.PARAM_IN_GENOME +
				' parameter is required for generating genome index')

        if (params.get(self.PARAM_IN_STARMODE, None) is not None and
		params[self.PARAM_IN_STARMODE] != "genomeGenerate"):
            if params.get(self.PARAM_IN_READS, None) is None:
		raise ValueError(self.PARAM_IN_READS +
				' parameter is required for reads mapping')

        if params.get(self.PARAM_IN_THREADN, None) is not None:
            if not isinstance(params[self.PARAM_IN_THREADN], int):
                raise ValueError(self.PARAM_IN_HASH_THREADN + ' must be of type int')
	else:
             params[self.PARAM_IN_THREADN] = 2

        if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
            if params[self.PARAM_IN_OUTFILE_PREFIX].find('/') != -1:
                raise ValueError(self.PARAM_IN_OUTFILE_PREFIX + ' cannot contain subfolder(s).')
        else:
            params[self.PARAM_IN_OUTFILE_PREFIX] = 'star_'

        if params.get('create_report', None) is None:
            params['create_report'] = 0

        return self._setDefaultParameters(params)


    def _setDefaultParameters(self, params_in):
        """set default for this group of parameters
        """
        params = copy.deepcopy(params_in)
        if params.get('outFilterType', None) is None:
            params['outFilterType'] = "\"BySJout\""
        if params.get('outFilterMultimapNmax', None) is None:
            params['outFilterMultimapNmax'] = 20
        if params.get('outSAMtype', None) is None:
            params['outSAMtype'] = 'BAM'
        if params.get('outSAMattrIHstart', None) is None:
            params['outSAMattrIHstart'] = 0
        if params.get('outSAMstrandField', None) is None:
            params['outSAMstrandField'] = 'intronMotif'
        if params.get('outFilterIntronMotifs', None) is None:
            params['outFilterIntronMotifs'] = 'RemoveNoncanonical'
        if params.get(self.SET_READS, None) is None:
            params[self.SET_READS] = self._get_reads_refs_from_setref(params)

        return params

    def _get_genome_gtf_file(self, gnm_ref, gtf_file_dir):
        """
        Get data from genome object ref and return the GTF filename (with path)
        for STAR indexing and mapping.
        STAR uses the reference annotation to guide assembly and for creating alignment
        """
        log("Converting genome {0} to GFF file in folder {1}".format(gnm_ref, gtf_file_dir))
        gfu = GenomeFileUtil(self.callback_url)
        try:
            gfu_ret = gfu.genome_to_gff({self.PARAM_IN_GENOME: gnm_ref,
                                         'is_gtf': 1,
                                         'target_dir': gtf_file_dir
                                      })
        except ValueError as egfu:
            log('GFU getting GTF file raised error:\n')
            pprint(egfu)
            return None
        else:#no exception raised
            return gfu_ret.get('file_path')


    def _construct_indexing_cmd(self, params):
	# STEP 1: construct the command for running `STAR --runMode genomeGenerate...`
        idx_cmd = [self.STAR_BIN]
	idx_cmd.append('--genomeDir')
	idx_cmd.append(params[self.STAR_IDX_DIR])
	idx_cmd.append('--' + self.PARAM_IN_STARMODE)
	idx_cmd.append('genomeGenerate')
	idx_cmd.append('--' + self.PARAM_IN_THREADN)
	idx_cmd.append(str(params[self.PARAM_IN_THREADN]))

	if params.get(self.PARAM_IN_FASTA_FILES, None) is not None:
            idx_cmd.append('--' + self.PARAM_IN_FASTA_FILES)
            for fasta_file in params[self.PARAM_IN_FASTA_FILES]:
                idx_cmd.append(fasta_file)

	# STEP 2: append the standard optional inputs
        if params.get('sjdbGTFfile', None) is not None:
            idx_cmd.append('--sjdbGTFfile')
            idx_cmd.append(params['sjdbGTFfile'])
        if (params.get('sjdbOverhang', None) is not None
		and params['sjdbOverhang'] > 0):
            idx_cmd.append('--sjdbOverhang')
            idx_cmd.append(str(params['sjdbOverhang']))

        #print ('STAR indexing CMD:')
        #print ' '.join(idx_cmd)
        return idx_cmd

    def _construct_mapping_cmd(self, params):
	if params.get(self.PARAM_IN_STARMODE, None) is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'

        # STEP 1: set the working folder housing the STAR output results as well as the reads info
        star_out_dir = ''
	if params.get('align_output', None) is None:
            star_out_dir = self.scratch
	else:
            star_out_dir = params['align_output']

        # STEP 2: construct the command for running STAR mapping
        mp_cmd = [self.STAR_BIN]
	mp_cmd.append('--genomeDir')
	mp_cmd.append(params[self.STAR_IDX_DIR])
	mp_cmd.append('--' + self.PARAM_IN_STARMODE)
	mp_cmd.append(params[self.PARAM_IN_STARMODE])
	mp_cmd.append('--' + self.PARAM_IN_THREADN)
	mp_cmd.append(str(params[self.PARAM_IN_THREADN]))

	if params.get(self.PARAM_IN_READS_FILES, None) is not None:
            #print('Input reads files:\n' + pformat(params[self.PARAM_IN_READS_FILES]))
            mp_cmd.append('--' + self.PARAM_IN_READS_FILES)
            for reads_file in params[self.PARAM_IN_READS_FILES]:
                mp_cmd.append(reads_file)
		readName, readsExtension = os.path.splitext(reads_file)
                #print ('Reads file name-- {}/extension-- {}:'.format(readName, readsExtension))
		if readsExtension == '.gz':
			mp_cmd.append('--readFilesCommand')
			mp_cmd.append('gunzip')
			mp_cmd.append('-c')

		if readsExtension == '.bz2':
			mp_cmd.append('--readFilesCommand')
			mp_cmd.append('bunzip2')
			mp_cmd.append('-c')

        # STEP 3: appending the advanced optional inputs
        mp_cmd.append('--' + self.PARAM_IN_OUTFILE_PREFIX)
        mp_cmd.append(os.path.join(star_out_dir, params[self.PARAM_IN_OUTFILE_PREFIX]))

        if params.get('sjdbGTFfile', None) is not None:
            mp_cmd.append('--sjdbGTFfile')
            mp_cmd.append(params['sjdbGTFfile'])
        if (params.get('sjdbOverhang', None) is not None
		and params['sjdbOverhang'] > 0):
            mp_cmd.append('--sjdbOverhang')
            mp_cmd.append(str(params['sjdbOverhang']))

        if (params.get('outFilterType', None) is not None
                and isinstance(params['outFilterType'], str)):
            mp_cmd.append('--outFilterType')
            mp_cmd.append(params['outFilterType'])
        if (params.get('outFilterMultimapNmax', None) is not None
                and isinstance(params['outFilterMultimapNmax'], int)
                and params['outFilterMultimapNmax'] >= 0):
            mp_cmd.append('--outFilterMultimapNmax')
            mp_cmd.append(str(params['outFilterMultimapNmax']))

        #output sorted file:Aligned.sortedByCoord.out.bam
        #allowed values of --outSAMtype are BAM Unsorted or SortedByCoordinate or both
        if params.get('outSAMtype', None) is not None:
            mp_cmd.append('--outSAMtype')
            mp_cmd.append(params['outSAMtype'])
            if params.get('outSAMtype', None) == 'BAM':
                mp_cmd.append('SortedByCoordinate')

        # 'It is recommended to remove the non-canonical junctions for Cnks runs using
        # --outFilterIntronMotifs RemoveNoncanonical'
        if params.get('outFilterIntronMotifs', None) is not None:
            mp_cmd.append('--outFilterIntronMotifs')
            mp_cmd.append('RemoveNoncanonical')

        if (params.get('outSAMattrIHstart', None) is not None
                and isinstance(params['outSAMattrIHstart'], int)
                and params['outSAMattrIHstart'] >= 0):
            mp_cmd.append('--outSAMattrIHstart')
            mp_cmd.append(str(params['outSAMattrIHstart']))
        if (params.get('outSAMstrandField', None) is not None
                and isinstance(params['outSAMstrandField'], str)):
            mp_cmd.append('--outSAMstrandField')
            mp_cmd.append(params['outSAMstrandField'])

        quant_modes = ["TranscriptomeSAM", "GeneCounts", "Both"]
        if (params.get('quantMode', None) is not None
                and params.get('quantMode', None) in quant_modes):
            mp_cmd.append('--quantMode')
            if params['quantMode'] == "Both":
                mp_cmd.append("TranscriptomeSAM")
                mp_cmd.append("GeneCounts")
            else:
                mp_cmd.append(params['quantMode'])
        if (params.get('alignSJoverhangMin', None) is not None
		and isinstance(params['alignSJoverhangMin'], int)
                and params['alignSJoverhangMin'] > 0):
            mp_cmd.append('--alignSJoverhangMin')
            mp_cmd.append(str(params['alignSJoverhangMin']))
        if (params.get('alignSJDBoverhangMin', None) is not None
                and isinstance(params['alignSJDBoverhangMin'], int)
                and params['alignSJDBoverhangMin'] > 0):
            mp_cmd.append('--alignSJDBoverhangMin')
            mp_cmd.append(str(params['alignSJDBoverhangMin']))
        if (params.get('outFilterMismatchNmax', None) is not None
		and isinstance(params['outFilterMismatchNmax'], int)
                and params['outFilterMismatchNmax'] > 0):
            mp_cmd.append('--outFilterMismatchNmax')
            mp_cmd.append(str(params['outFilterMismatchNmax']))
        if (params.get('alignIntronMin', None) is not None
		and isinstance(params['alignIntronMin'], int)
                and params['alignIntronMin'] > 0):
            mp_cmd.append('--alignIntronMin')
            mp_cmd.append(str(params['alignIntronMin']))
        if (params.get('alignIntronMax', None) is not None
		and isinstance(params['alignIntronMax'], int)
                and params['alignIntronMax'] >= 0):
            mp_cmd.append('--alignIntronMax')
            mp_cmd.append(str(params['alignIntronMax']))
        if (params.get('alignMatesGapMax', None) is not None
		and isinstance(params['alignMatesGapMax'], int)
                and params['alignMatesGapMax'] >= 0):
            mp_cmd.append('--alignMatesGapMax')
            mp_cmd.append(str(params['alignMatesGapMax']))

        #print ('STAR mapping CMD:')
        #print ' '.join(mp_cmd)
        return mp_cmd

    def _exec_indexing(self, params):
        log('Running STAR index generating with params:\n' + pformat(params))

        idx_cmd = self._construct_indexing_cmd(params)

        exitCode = self.prog_runner.run(idx_cmd, self.scratch)

        return exitCode

    def _exec_mapping(self, params):
        log('Running STAR mapping with params:\n' + pformat(params))

        mp_cmd = self._construct_mapping_cmd(params)

        exitCode = self.prog_runner.run(mp_cmd, self.scratch)

        return exitCode

    def _exec_star_pipeline(self, params, rds_files, rds_name, idx_dir, out_dir):
        # build the parameters
        params_idx = self._get_indexing_params(params, idx_dir)
        params_mp = self._get_mapping_params(params, rds_files, rds_name, idx_dir, out_dir)

        # execute indexing and then mapping
        retVal = {}
        try:
            if params[self.PARAM_IN_STARMODE]=='genomeGenerate':
                ret = self._exec_indexing(params_idx)
            else:
		ret = 0
            while( ret != 0 ):
                time.sleep(1)
        except ValueError as eidx:
            log('STAR genome indexing raised error:\n')
            pprint(eidx)
        else:#no exception raised by genome indexing and STAR returns 0, then run mapping
            params_mp[self.PARAM_IN_STARMODE] = 'alignReads'
            try:
                ret = self._exec_mapping(params_mp)
                while( ret != 0 ):
                    time.sleep(1)
            except ValueError as emp:
                log('STAR mapping raised error:\n')
                pprint(emp)
            else:#no exception raised by STAR mapping and STAR returns 0, then move to saving and reporting  
                ret = {'star_idx': star_idx, 'star_output': params_mp.get('align_output')}

        return ret


    def upload_STARalignment(self, input_params, reads_ref, reads_info, output_bam_file):
        """
        Uploads the alignment file + metadata.
        Returns the STAR alignment reference.
        """

        aligner_opts = dict()
        for k in input_params:
            aligner_opts[k] = str(input_params[k])
        pprint(reads_info)

        alignment_name = reads_ref['alignment_output_name']
        align_upload_params = {
            "destination_ref": "{}/{}".format(input_params[self.PARAM_IN_WS], alignment_name),
            "file_path": output_bam_file,
            "assembly_or_genome_ref": input_params[self.PARAM_IN_GENOME],
            "read_library_ref": reads_info['object_ref'],
            "library_type": reads_info['style'],
            "condition": reads_info['condition'],
            "aligned_using": 'STAR',
            "aligner_version":self.STAR_VERSION,
            "aligner_opts": aligner_opts
        }

        pprint(align_upload_params)

        ra_util = ReadsAlignmentUtils(self.callback_url, service_ver='dev')
        rau_upload_ret = ra_util.upload_alignment(align_upload_params)
        alignment_ref = rau_upload_ret["obj_ref"]
        print("STAR alignment uploaded as object {}".format(alignment_ref))
        return rau_upload_ret


    def generate_report_for_single_run(self, run_output_info, params):
        input_ref = run_output_info['upload_results']['obj_ref']
        index_dir = run_output_info['index_dir']
        output_dir = run_output_info['output_dir']
        output_files = self._generate_output_file_list(index_dir, output_dir)

        # first run qualimap
        qualimap_report = self.qualimap.run_bamqc({'input_ref': input_ref})
        qc_result_zip_info = qualimap_report['qc_result_zip_info']

        # create report
        report_text = 'Ran on a single reads library.\n\n'
        alignment_info = self.get_obj_infos(input_ref)[0]
        report_text = 'Created ReadsAlignment: ' + str(alignment_info[1]) + '\n'
        report_text += '                        ' + input_ref + '\n'
        kbr = KBaseReport(self.callback_url)
        report_info = kbr.create_extended_report({'message': report_text,
                                                  'file_links': output_files,
                                                  'objects_created': [{'ref': input_ref,
                                                                       'description': 'ReadsAlignment'}],
                                                  'report_object_name': 'kb_STAR_report_' + str(uuid.uuid4()),
                                                  'direct_html_link_index': 0,
                                                  'html_links': [{'shock_id': qc_result_zip_info['shock_id'],
                                                                  'name': qc_result_zip_info['index_html_file_name'],
                                                                  'label': qc_result_zip_info['name']}],
                                                  'html_window_height': 366,
                                                  'workspace_name': params['output_workspace']
                                                  })
        return report_info #{'report_name': report_info['name'], 'report_ref': report_info['ref']}

    def _get_reads_info(self, reads, readsSet_ref):
        '''
        _get_reads_info:fetches the detailed info for each reads with ref in list reads_refs
        return an object of the following structure:
        {
            "style": "paired", "single", or "interleaved",
            "file_fwd": path_to_file,
            "name": name of the reads,
            "file_rev": path_to_file, only if paired end,
            "object_ref": reads reference for downstream convenience,
            "condition": the condition for the reads.
        }
        '''
        try:
            print("Fetching FASTA file from reads reference {}".format(reads['ref']))
            ret_reads_info = fetch_reads_from_reference(reads['ref'], self.callback_url)
        except ValueError:
            print("Incorrect object type for fetching a FASTA file!")
            raise

        if ret_reads_info.get("file_fwd", None) is None:
            raise RuntimeError("FASTA file fetched from reads {} doesn't seem to exist!".format(reads['ref']))
        else:
            if reads.get('condition', None) is not None:
                ret_reads_info['condition'] = reads['condition']
            else:
                ret_reads_info['condition'] = 'unspecified'
            if reads.get('object_ref', None) != readsSet_ref:
                ret_reads_info[self.PARAM_IN_READS] = readsSet_ref

        return ret_reads_info


    def _get_genome_fasta(self, gnm_ref):
        genome_fasta_files = list()
	if gnm_ref is not None:
            try:
		print("Fetching FASTA file from object {}".format(gnm_ref))
		genome_fasta_file = fetch_fasta_from_object(gnm_ref, self.workspace_url, self.callback_url)
		print("Done fetching FASTA file! Path = {}".format(genome_fasta_file.get("path", None)))
            except ValueError:
		print("Incorrect object type for fetching a FASTA file!")
		raise

            if genome_fasta_file.get("path", None) is None:
		raise RuntimeError("FASTA file fetched from object {} doesn't seem exist!".format(gnm_ref))
            else:
		genome_fasta_files.append(genome_fasta_file["path"])
        return genome_fasta_files


    def convert_params(self, validated_params):
        """
        Convert input parameters with KBase ref format into STAR parameters,
        and add the advanced options.
        """
        params = copy.deepcopy(validated_params)
        params['runMode'] = 'genomeGenerate'

        if validated_params.get('create_report', None) is not None:
                params['create_report'] = validated_params['create_report']
        if validated_params.get('concurrent_local_tasks', None) is not None:
                params['concurrent_local_tasks'] = validated_params['concurrent_local_tasks']
        if validated_params.get('concurrent_njsw_tasks', None) is not None:
                params['concurrent_njsw_tasks'] = validated_params['concurrent_njsw_tasks']
        if validated_params.get('alignmentset_suffix', None) is not None:
                params['alignmentset_suffix'] = validated_params['alignmentset_suffix']

        # Add advanced options from validated_params to params
        sjdbGTFfile = validated_params.get("sjdbGTFfile", None)
	if sjdbGTFfile is not None:
            params['sjdbGTFfile'] = sjdbGTFfile
        else:
            params['sjdbGTFfile'] = self._get_genome_gtf_file(
                                        params[self.PARAM_IN_GENOME],
                                        os.path.join(self.scratch, self.STAR_IDX_DIR))
        if validated_params.get('sjdbOverhang', None) is not None :
            params['sjdbOverhang'] = validated_params['sjdbOverhang']
        else:
            params['sjdbOverhang'] = 100

        quant_modes = ["TranscriptomeSAM", "GeneCounts", "Both"]
        if (validated_params.get('quantMode', None) is not None
                and validated_params.get('quantMode', None) in quant_modes):
            params['quantMode'] = validated_params['quantMode']
        else:
            params['quantMode'] = 'Both'

        return params


    def _get_indexing_params(self, params, star_idx_dir):
        params_idx = {
                'runMode': 'genomeGenerate',
		'runThreadN': params[self.PARAM_IN_THREADN],
		self.STAR_IDX_DIR: star_idx_dir,
                'genomeFastaFiles': params[self.PARAM_IN_FASTA_FILES]
        }
        if params.get('sjdbGTFfile', None) is not None:
            params_idx['sjdbGTFfile'] = params['sjdbGTFfile']
        if params.get('sjdbOverhang', None) is not None :
            params_idx['sjdbOverhang'] = params['sjdbOverhang']

        return params_idx


    def _get_mapping_params(self, params, rds_files, rds_name, idx_dir, out_dir):
        ''' build the mapping parameters'''
        aligndir = out_dir
        if rds_name:
            aligndir = os.path.join(out_dir, rds_name)
            self._mkdir_p(aligndir)
            #print '**********STAR output directory created:{}'.format(aligndir)

        params_mp = copy.deepcopy(params)
        params_mp['runMode'] = 'alignReads'
        params_mp['readFilesIn'] = rds_files
	params_mp[self.STAR_IDX_DIR] = idx_dir
        params_mp['align_output'] = aligndir

        return params_mp


    def determine_input_info(self, validated_params):
        ''' get info on the readsset_ref object and determine if we run once or run on a set
        input info provides information on the input and tells us if we should
        run as a single_library or as a set:
             input_info = {'run_mode': '', 'info': [..], 'ref': '55/1/2'}
        '''
        info = self.get_obj_infos(validated_params[self.PARAM_IN_READS])[0]
        obj_type = self.get_type_from_obj_info(info)
        if obj_type in ['KBaseAssembly.PairedEndLibrary', 'KBaseAssembly.SingleEndLibrary',
                        'KBaseFile.PairedEndLibrary', 'KBaseFile.SingleEndLibrary']:
            return {'run_mode': 'single_library', 'info': info, 'ref': validated_params[self.PARAM_IN_READS]}
        if obj_type == 'KBaseRNASeq.RNASeqSampleSet':
            return {'run_mode': 'sample_set', 'info': info, 'ref': validated_params[self.PARAM_IN_READS]}
        if obj_type == 'KBaseSets.ReadsSet':
            return {'run_mode': 'sample_set', 'info': info, 'ref': validated_params[self.PARAM_IN_READS]}

        raise ValueError('Object type of readsset_ref is not valid, was: ' + str(obj_type))

    def determine_unique_reads_names(self, validated_params):
        infos = self.get_obj_infos(validated_params[self.PARAM_IN_READS])
        return get_unique_names(infos)

    def get_type_from_obj_info(self, info):
        return info[2].split('-')[0]

    def get_name_from_obj_info(self, info):
        return info[1]

    def get_obj_infos(self, ref):
        return self.ws_client.get_object_info3({'objects': [{'ref': ref}]})['infos']

    def get_object_names(self, ref_list):
        """
        From a list of workspace references, returns a mapping from ref -> name of the object.
        """
        obj_ids = list()
        for ref in ref_list:
            obj_ids.append({"ref": ref})
        info = self.ws_client.get_object_info3({"objects": obj_ids})
        name_map = dict()
        # we already have the refs as passed previously, so use those for mapping, as they're in
        # the same order as what's returned.
        for i in range(len(info["infos"])):
            name_map[ref_list[i]] = info["infos"][i][1]
        return name_map


    def _mkdir_p(self, dir):
        """
        _mkdir_p: make directory for given path
        """
        log('Creating a new dir: ' + dir)
        if not dir:
            return
        if not os.path.exists(dir):
            os.makedirs(dir)
        else:
            log('{} has existed, so skip creating.'.format(dir))


    def create_star_dirs(self, star_home):
        '''creating the directories for STAR'''
        # the index directory
        idxdir = os.path.join(star_home, self.STAR_IDX_DIR)
        self._mkdir_p(idxdir)
        # the output directory
        outdir = os.path.join(star_home, self.STAR_OUT_DIR)
        self._mkdir_p(outdir)

        return (idxdir, outdir)


    def _get_reads_refs_from_setref(self, params):
        readsSet_ref = params[self.PARAM_IN_READS]
        reads_refs = list()
        try:
            #print("Fetching reads ref(s) from sample/reads set ref {}".format(readsSet_ref))
            reads_refs = fetch_reads_refs_from_sampleset(
                                    readsSet_ref,
                                    self.workspace_url,
                                    self.callback_url,
                                    params)
            #print("\nDone fetching reads ref(s) from readsSet {}--\nDetails:\n".format(readsSet_ref))
            #pprint(reads_refs)
        except ValueError:
            print("Incorrect object type for fetching reads ref(s)!")
            raise

        return reads_refs

    def _generate_output_file_list(self, idx_dir, out_dir):
        """
        _generate_output_file_list: zip result files and generate file_links for report
        """

        log('start packing result files')

        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        star_index = os.path.join(output_directory, 'star_index.zip')
        star_output = os.path.join(output_directory, 'star_output.zip')
        self.zip_folder(idx_dir, star_index)
        self.zip_folder(out_dir, star_output)

        #star_index = self.zip_folder_withDFU(idx_dir, 'star_index')
        #star_output = self.zip_folder_withDFU(out_dir, 'star_output')

        output_files.append({'path': star_index,
                             'name': os.path.basename(star_index),
                             'label': os.path.basename(star_index),
                             'description': 'Index file(s) generated by STAR'})

        output_files.append({'path': star_output,
                             'name': os.path.basename(star_output),
                             'label': os.path.basename(star_output),
                             'description': 'Output file(s) generated by STAR'})

        return output_files


    def zip_folder_withDFU(self, folder_path, output_name):
        """Zip the contents of an entire folder (with that folder included
        in the archive). Empty subfolders will be included in the archive
        as well.
        """
        output_path = self.dfu.pack_file(
                {'file_path': folder_path + '/' + output_name,
                 'pack': 'zip'})['file_path']

        print "{} created successfully.".format(output_path)

        with zipfile.ZipFile(output_path, "r") as f:
            print 'Checking the zipped file......\n'
            #for info in f.infolist():
                #    print info.filename, info.date_time, info.file_size, info.compress_size
            for fn in f.namelist():
                print fn

        return output_path


    def zip_folder(self, folder_path, output_path):
        """Zip the contents of an entire folder (with that folder included in the archive). 
        Empty subfolders could be included in the archive as well if the commented portion is used.
        """
        with zipfile.ZipFile(output_path, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as ziph:
            for root, folders, files in os.walk(folder_path):
                # Include all subfolders, including empty ones.
                #for folder_name in folders:
                #    absolute_path = os.path.join(root, folder_name)
                #    relative_path = os.path.join(os.path.basename(root), folder_name)
                #    print "Adding {} to archive.".format(absolute_path)
                #    ziph.write(absolute_path, relative_path)
                for f in files:
                    absolute_path = os.path.join(root, f)
                    relative_path = os.path.join(os.path.basename(root), f)
                    print "Adding {} to archive.".format(absolute_path)
                    ziph.write(absolute_path, relative_path)

        print "{} created successfully.".format(output_path)

        with zipfile.ZipFile(output_path, "r") as f:
            print 'Checking the zipped file......\n'
            for info in f.infolist():
                print info.filename, info.date_time, info.file_size, info.compress_size


    def _generate_html_report(self, out_dir, obj_ref):
        """
        _generate_html_report: generate html summary report
        """

        log('start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        star_obj = self.ws_client.get_objects2({'objects':
                                                 [{'ref': obj_ref}]})['data'][0]
        star_obj_info = star_obj['info']
        star_obj_data = star_obj['data']
        star_obj_type = star_obj_info[2]

        Overview_Content = ''
        if re.match('KBaseRNASeq.RNASeqAlignment-\d.\d', star_obj_type):
            Overview_Content += '<br/><table><tr><th>Generated Alignment Object</th>'
            Overview_Content += '<th></th></tr>'
            Overview_Content += '<tr><th>Alignment Name</th><th>Condition</th></tr>'
            Overview_Content += '<tr><td>{} ({})</td>'.format(star_obj_info[1],obj_ref)
            Overview_Content += '<td>{}</td></tr>'.format(star_obj_data['condition'])
            Overview_Content += '</table>'
        elif (re.match('KBaseRNASeq.RNASeqAlignmentSet-\d.\d', star_obj_type)
                or re.match('KBaseSets.ReadsAlignmentSet-\d.\d', star_obj_type)
                or re.match('KBaseSet.RNASeqAlignmentSet-\d.\d', star_obj_type)):
            Overview_Content += '<br/><table><tr><th>Generated AlignmentSet Object</th></tr>'
            Overview_Content += '<tr><td>{} ({})'.format(star_obj_info[1],obj_ref)
            Overview_Content += '</td></tr></table>'
            Overview_Content += '<p><br/></p>'
            Overview_Content += '<table><tr><th>Generated Alignment Objects</th>'
            Overview_Content += '<th></th></tr>'
            Overview_Content += self._fill_html_trs('Alignment Name', star_obj_data)
            Overview_Content += '</table>'
        elif re.match('KBaseRNASeq.RNASeqExpression-\d.\d', star_obj_type):
            Overview_Content += '<br/><table><tr><th>Generated Expression Object</th>'
            Overview_Content += '<th></th></tr>'
            Overview_Content += '<tr><th>Expression Name</th><th>Condition</th></tr>'
            Overview_Content += '<tr><td>{} ({})</td>'.format(star_obj_info[1], obj_ref)
            Overview_Content += '<td>{}</td></tr>'.format(star_obj_data['condition'])
            Overview_Content += '</table>'
        elif re.match('KBaseSets.ExpressionSet-\d.\d', star_obj_type):
            Overview_Content += '<br/><table><tr><th>Generated ExpressionSet Object</th></tr>'
            Overview_Content += '<tr><td>{} ({})'.format(star_obj_info[1], obj_ref)
            Overview_Content += '</td></tr></table>'
            Overview_Content += '<p><br/></p>'
            Overview_Content += '<table><tr><th>Generated Expression Objects</th>'
            Overview_Content += '<th></th></tr>'
            Overview_Content += self._fill_html_trs('Expression Name', star_obj_data)
            Overview_Content += '</table>'

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for STAR'})
        return html_report

    def _fill_html_trs(self, col_caption, obj_data):
        '''
        _fill_html_trs: simple creates an html string that has rows (tr) of td for a table
        '''
        tr_html_str = '<tr><th>{}</th><th>Condition</th></tr>'.format(col_caption)

        for item in obj_data['items']:
            item_obj = self.ws_client.get_objects2({'objects':[{'ref': item['ref']}]})['data'][0]
            item_obj_info = item_obj['info']
            item_obj_data = item_obj['data']
            obj_name = item_obj_info[1]

            tr_html_str += '<tr><td>{} ({})</td>'.format(obj_name, item['ref'])
            tr_html_str += '<td>{}</td></tr>'.format(item_obj_data['condition'])

        return tr_html_str

    def _generate_star_report(self, obj_ref, report_text, html_links, workspace_name, index_dir, output_dir):
        """
        _generate_star_report: generate summary report
        """
        log('creating STAR report')

        output_files = self._generate_output_file_list(index_dir, output_dir)
        output_html_files = self._generate_html_report(output_dir, obj_ref)
        output_html_files += html_links

        star_obj = self.ws_client.get_objects2({'objects':[{'ref': obj_ref}]})['data'][0]
        star_obj_info = star_obj['info']
        star_obj_data = star_obj['data']

        star_obj_type = star_obj_info[2]
        if re.match('KBaseRNASeq.RNASeqAlignment-\d+.\d+', star_obj_type):
            objects_created = [{'ref': obj_ref,
                                'description': 'RNASeqAlignment generated by STAR'}]
        elif (re.match('KBaseRNASeq.RNASeqAlignmentSet-\d+.\d+', star_obj_type)
                or re.match('KBaseSets.ReadsAlignmentSet-\d+.\d+', star_obj_type)
                or re.match('KBaseSet.RNASeqAlignmentSet-\d+.\d+', star_obj_type)):
            objects_created = [{'ref': obj_ref,
                'description': '{} generated by STAR'.format(re.sub(r"-\d+.\d+", "",star_obj_type))}]
            items = star_obj_data['items']
            for item in items:
                objects_created.append({'ref': item['ref'],
                                        'description': 'Alignment generated by STAR'})
        elif re.match('KBaseRNASeq.RNASeqExpression-\d+.\d+', star_obj_type):
            objects_created = [{'ref': obj_ref,
                                'description': 'Expression generated by STAR'}]
        elif re.match('KBaseSets.ExpressionSet-\d+.\d+', star_obj_type):
            objects_created = [{'ref': obj_ref,
                                'description': 'ExpressionSet generated by STAR'}]
            items = star_obj_data['items']
            for item in items:
                objects_created.append({'ref': item['ref'],
                                        'description': 'Expression generated by STAR'})

        report_params = {'message': report_text,
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'objects_created': objects_created,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 366,
                         'report_object_name': 'kb_STAR_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        report_output = kbase_report_client.create_extended_report(report_params)

        return report_output

    def upload_alignment_set(self, alignment_items, alignmentset_name, ws_name):
        """
        Compiles and saves a set of alignment references (+ other stuff) into a
        KBaseRNASeq.RNASeqAlignmentSet.
        Returns the reference to the new alignment set.
        alignment_items: [{
            "ref": alignment_ref,
            "label": condition label.
        }]
        """
        print("Uploading completed alignment set")
        alignment_set = {
            "description": "Alignments using STAR, v.{}".format(self.STAR_VERSION),
            "items": alignment_items
        }
        set_info = self.set_api_client.save_reads_alignment_set_v1({
            "workspace": ws_name,
            "output_object_name": alignmentset_name,
            "data": alignment_set
        })
        return set_info

