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
    fetch_reads_from_reference
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

    def __init__(self, config, provenance):
        self.config = config
        self.workspace_url = config['workspace-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.srv_wiz_url = config['srv-wiz-url']
        self.au = AssemblyUtil(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.scratch = config['scratch']
        self.working_dir = self.scratch
        self.prog_runner = Program_Runner(self.STAR_BIN, self.scratch)
        self.provenance = provenance
        self.ws_client = Workspace(self.workspace_url, token=self.token)

        # from the provenance, extract out the version to run by exact hash if possible
        self.my_version = 'release'
        if len(provenance) > 0:
            if 'subactions' in provenance[0]:
                self.my_version = self.get_version_from_subactions('kb_STAR', provenance[0]['subactions'])
        print('Running kb_STAR version = ' + self.my_version)

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
        if "expression_suffix" not in params or not valid_string(params["expression_suffix"]):
            raise ValueError("Parameter expression_suffix must be a valid Workspace object string, "
                      "not {}".format(params.get("expression_suffix", None)))

        if params.get(self.PARAM_IN_STARMODE, None) is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'
	else:
            if params[self.PARAM_IN_STARMODE] == "genomeGenerate":
		if params.get(self.PARAM_IN_GENOME, None) is None:
                    raise ValueError(self.PARAM_IN_GENOME +
				' parameter is required for generating genome index')
                params['sjdbGTFfile'] = self._get_genome_gtf_file(
                                        params[self.PARAM_IN_GENOME], os.path.join(self.scratch, self.STAR_IDX_DIR))

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

        if params.get('create_report', None) is None:
            params['create_report'] = 0

        return self._setDefaultParameters(params)


    def _setDefaultParameters(self, params):
        """set default for this group of parameters
        """
        if params.get('outFilterType', None) is None:
            params['outFilterType'] = "\"BySJout\""
        if params.get('outFilterMultimapNmax', None) is None:
            params['outFilterMultimapNmax'] = 20
        if params.get('outSAMtype', None) is not None:
            params['outSAMtype'] = "BAM" #SortedByCoordinate
        if params.get('outSAMattrIHstart', None) is None:
            params['outSAMattrIHstart'] = 0
        if params.get('outSAMstrandField', None) is None:
            params['outSAMstrandField'] = "intronMotif"

        return params


    def _get_genome_gtf_file(self, gnm_ref, gtf_file_dir):
        """
        Get data from genome object ref and return the GTF filename (with path)
        for STAR indexing and mapping.
        STAR uses the reference annotation to guide assembly and for creating alignment
        """
        log("Converting genome {0} to GFF file {1}".format(gnm_ref, gtf_file_dir))
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
            print('Input fasta reads files:' + pformat(params[self.PARAM_IN_FASTA_FILES]))
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

        # STEP 3: return idx_cmd
        print ('STAR indexing CMD:')
        print ' '.join(idx_cmd)
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
            print('Input reads files:' + pformat(params[self.PARAM_IN_READS_FILES]))
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
        if (params.get('outSAMtype', None) is not None
                and isinstance(params['outSAMtype'], str)):
            mp_cmd.append('--outSAMType')
            mp_cmd.append(params['outSAMType'])
            mp_cmd.append('SortedByCoordinate')#BAM SortedByCoordinate
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

        # STEP 4: return mp_cmd
        print ('STAR mapping CMD:')
        print ' '.join(mp_cmd)
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


    def upload_STARalignment(self, input_params, reads, output_sam_file):
        """
        Uploads the alignment file + metadata.
        Returns the STAR alignment reference.
        """
        reads_ref = reads.get('readsRefs', None)[0]
        reads_info = reads.get('readsInfo', None)[0]

        aligner_opts = dict()
        for k in input_params:
            aligner_opts[k] = str(input_params[k])
        pprint(reads_info)

        alignment_name = reads_ref['alignment_output_name']
        align_upload_params = {
            "destination_ref": "{}/{}".format(input_params[self.PARAM_IN_WS], alignment_name),
            "file_path": output_sam_file,
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
                                                  'objects_created': [{'ref': input_ref,
                                                                       'description': 'ReadsAlignment'}],
                                                  'report_object_name': 'kb_STAR_' + str(uuid.uuid4()),
                                                  'direct_html_link_index': 0,
                                                  'html_links': [{'shock_id': qc_result_zip_info['shock_id'],
                                                                  'name': qc_result_zip_info['index_html_file_name'],
                                                                  'label': qc_result_zip_info['name']}],
                                                  'workspace_name': params['output_workspace']
                                                  })
        return report_info #{'report_name': report_info['name'], 'report_ref': report_info['ref']}



    def _get_reads(self, params):
        '''fetch the refs and info from the input given by params[self.PARAM_IN_READS], which is a ref
        '''
        reads_refs = list()
        radsset_ref = params[self.PARAM_IN_READS]
	if radsset_ref is not None:
            try:
		print("Fetching reads ref(s) from sample/reads set ref {}".format(radsset_ref))
		reads_refs = fetch_reads_refs_from_sampleset(radsset_ref, self.workspace_url, self.callback_url, params)
		print("Done fetching reads ref(s)!")
            except ValueError:
		print("Incorrect object type for fetching reads ref(s)!")
		raise

        reads_info = list()
	for source_reads in reads_refs:
            try:
                print("Fetching FASTA file from reads reference {}".format(source_reads['ref']))
                ret_reads = fetch_reads_from_reference(source_reads['ref'], self.callback_url)
            except ValueError:
                print("Incorrect object type for fetching a FASTA file!")
                raise

            if ret_reads.get("file_fwd", None) is None:
                raise RuntimeError("FASTA file fetched from reads ref {} doesn't seem to exist!".format(source_reads['ref']))
            else:
                if source_reads.get('condition', None) is not None:
                    ret_reads['condition'] = source_reads['condition']
                else:
                    ret_reads['condition'] = params['condition']
                if ret_reads.get('object_ref', None) != radsset_ref:
                    ret_reads[self.PARAM_IN_READS] = radsset_ref
                reads_info.append(ret_reads)

        return {'readsRefs': reads_refs, 'readsInfo': reads_info}


    def _get_genome_fasta(self, genome_ref):
        genome_fasta_files = list()
	if genome_ref is not None:
            try:
		print("Fetching FASTA file from object {}".format(genome_ref))
		genome_fasta_file = fetch_fasta_from_object(genome_ref, self.workspace_url, self.callback_url)
		print("Done fetching FASTA file! Path = {}".format(genome_fasta_file.get("path", None)))
            except ValueError:
		print("Incorrect object type for fetching a FASTA file!")
		raise

            if genome_fasta_file.get("path", None) is None:
		raise RuntimeError("FASTA file fetched from object {} doesn't seem exist!".format(genome_ref))
            else:
		genome_fasta_files.append(genome_fasta_file["path"])
        return genome_fasta_files


    def convert_params(self, input_params):
        """
        Convert input parameters with KBase ref format into STAR parameters,
        and add the advanced options.
        """
	params = {
            'output_workspace': input_params[self.PARAM_IN_WS],
            'runMode': 'genomeGenerate',
            'alignment_suffix': input_params['alignment_suffix'],
            'expression_suffix': input_params['expression_suffix'],
            self.PARAM_IN_GENOME: input_params[self.PARAM_IN_GENOME],
            'runThreadN': input_params[self.PARAM_IN_THREADN]
	}

        if input_params.get('create_report', None) is not None:
                params['create_report'] = input_params['create_report']
        if input_params.get('concurrent_local_tasks', None) is not None:
                params['concurrent_local_tasks'] = input_params['concurrent_local_tasks']
        if input_params.get('concurrent_njsw_tasks', None) is not None:
                params['concurrent_njsw_tasks'] = input_params['concurrent_njsw_tasks']
        if input_params.get('alignmentset_suffix', None) is not None:
                params['alignmentset_suffix'] = input_params['alignmentset_suffix']
        if input_params.get('expression_set_suffix', None) is not None:
                params['expression_set_suffix'] = input_params['expression_set_suffix']
	# STEP 1: Converting refs to file locations in the scratch area
        reads = self._get_reads(input_params)

        params[self.PARAM_IN_FASTA_FILES] = self._get_genome_fasta(input_params[self.PARAM_IN_GENOME])

        # STEP 2: Add advanced options from input_params to params
        sjdbGTFfile = input_params.get("sjdbGTFfile", None)
	if sjdbGTFfile is not None:
            params['sjdbGTFfile'] = sjdbGTFfile
            if input_params.get('sjdbOverhang', None) is not None :
                params['sjdbOverhang'] = input_params['sjdbOverhang']
            else:
                params['sjdbOverhang'] = 100

        if input_params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
            params[self.PARAM_IN_OUTFILE_PREFIX] = input_params[self.PARAM_IN_OUTFILE_PREFIX]
        else:
            params[self.PARAM_IN_OUTFILE_PREFIX] = 'star_'

        if input_params.get('outFilterType', None) is not None:
            params['outFilterType'] = input_params['outFilterType']
        if input_params.get('outFilterMultimapNmax', None) is not None:
            params['outFilterMultimapNmax'] = input_params['outFilterMultimapNmax']
        if input_params.get('outSAMtype', None) is not None:
            params['outSAMType'] = input_params['outSAMType']
        if input_params.get('outSAMattrIHstart', None) is not None:
            params['outSAMattrIHstart'] = input_params['outSAMattrIHstart']
        if input_params.get('outSAMstrandField', None) is not None:
            params['outSAMstrandField'] = input_params['outSAMstrandField']

        quant_modes = ["TranscriptomeSAM", "GeneCounts", "Both"]
        if (input_params.get('quantMode', None) is not None
                and input_params.get('quantMode', None) in quant_modes):
            params['quantMode'] = input_params['quantMode']
        else:
            params['quantMode'] = 'Both'
        if input_params.get('alignSJoverhangMin', None) is not None:
            params['alignSJoverhangMin'] = input_params['alignSJoverhangMin']
        if (input_params.get('alignSJDBoverhangMin', None) is not None
                and isinstance(input_params['alignSJDBoverhangMin'], int)
                and input_params['alignSJDBoverhangMin'] > 0):
            params['alignSJDBoverhangMin'] = input_params['alignSJDBoverhangMin']
        if input_params.get('outFilterMismatchNmax', None) is not None:
            params['outFilterMismatchNmax'] = input_params['outFilterMismatchNmax']
        if input_params.get('alignIntronMin', None) is not None:
            params['alignIntronMin'] = input_params['alignIntronMin']
        if input_params.get('alignIntronMax', None) is not None:
            params['alignIntronMax'] = input_params['alignIntronMax']
        if input_params.get('alignMatesGapMax', None) is not None:
            params['alignMatesGapMax'] = input_params['alignMatesGapMax']

        return {'input_parameters': params, 'reads': reads}


    def _get_indexing_params(self, params, star_idx_dir):
        params_idx = {
                'runMode': params[self.PARAM_IN_STARMODE],
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
            print '\n**********STAR output directory created: ' + aligndir

        params_mp = {
                'runMode': 'alignReads',#params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'readFilesIn': rds_files,
		self.STAR_IDX_DIR: idx_dir,
		'align_output': aligndir
        }

        if params.get('sjdbGTFfile', None) is not None:
            params_mp['sjdbGTFfile'] = params['sjdbGTFfile']
        if params.get('sjdbOverhang', None) is not None :
            params_mp['sjdbOverhang'] = params['sjdbOverhang']

        if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
            params_mp[self.PARAM_IN_OUTFILE_PREFIX] = params[self.PARAM_IN_OUTFILE_PREFIX]

        if (params.get('outFilterType', None) is not None
                and isinstance(params['outFilterType'], str)):
            params_mp['outFilterType'] = params['outFilterType']
        if (params.get('outFilterMultimapNmax', None) is not None
                and isinstance(params['outFilterMultimapNmax'], int)):
            params_mp['outFilterMultimapNmax'] = params['outFilterMultimapNmax']
        if (params.get('outSAMtype', None) is not None
                and isinstance(params['outSAMtype'], str)):
            params_mp['outSAMtype'] = params['outSAMtype']
        if (params.get('outSAMattrIHstart', None) is not None
                and isinstance(params['outSAMattrIHstart'], int)):
            params_mp['outSAMattrIHstart'] = params['outSAMattrIHstart']
        if (params.get('outSAMstrandField', None) is not None
                and isinstance(params['outSAMstrandField'], str)):
            params_mp['outSAMstrandField'] = params['outSAMstrandField']

        quant_modes = ["TranscriptomeSAM", "GeneCounts", "Both"]
        if (params.get('quantMode', None) is not None
                and params.get('quantMode', None) in quant_modes):
            params_mp['quantMode'] = params['quantMode']
        if (params.get('alignSJoverhangMin', None) is not None
		and isinstance(params['alignSJoverhangMin'], int)
                and params['alignSJoverhangMin'] > 0):
            params_mp['alignSJoverhangMin'] = params['alignSJoverhangMin']
        if (params.get('alignSJDBoverhangMin', None) is not None
                and isinstance(params['alignSJDBoverhangMin'], int)
                and params['alignSJDBoverhangMin'] > 0):
            params_mp['alignSJDBoverhangMin'] = params['alignSJDBoverhangMin']
        if (params.get('outFilterMismatchNmax', None) is not None
		and isinstance(params['outFilterMismatchNmax'], int)
                and params['outFilterMismatchNmax'] > 0):
            params_mp['outFilterMismatchNmax'] = params['outFilterMismatchNmax']
        if (params.get('alignIntronMin', None) is not None
		and isinstance(params['alignIntronMin'], int)
                and params['alignIntronMin'] > 0):
            params_mp['alignIntronMin'] = params['alignIntronMin']
        if (params.get('alignIntronMax', None) is not None
		and isinstance(params['alignIntronMax'], int)
                and params['alignIntronMax'] >= 0):
            params_mp['alignIntronMax'] = params['alignIntronMax']
        if (params.get('alignMatesGapMax', None) is not None
		and isinstance(params['alignMatesGapMax'], int)
                and params['alignMatesGapMax'] >= 0):
            params_mp['alignMatesGapMax'] = params['alignMatesGapMax']

        return params_mp


    def build_single_execution_task(self, rds_ref, params):
        task_params = copy.deepcopy(params)

        task_params[self.PARAM_IN_READS] = rds_ref
        task_params['create_report'] = 0

        if 'condition' in rds_ref:
                task_param['condition'] = rds_ref['condition']
        else:
                task_params['condition'] = 'unspecified'

        return {'module_name': 'STAR',
                'function_name': 'run_star',
                'version': self.my_version,
                #'version': 'dev',
                'parameters': task_params}


    def process_batch_result(self, batch_result, params, reads_refs):
        n_jobs = len(batch_result['results'])
        n_success = 0
        n_error = 0
        ran_locally = 0
        ran_njsw = 0

        set_name = self.get_object_names([params[self.PARAM_IN_READS]])[params[self.PARAM_IN_READS]]
        # reads alignment set items
        alignment_items = []
        alignment_objs = []

        for k in range(0, len(batch_result['results'])):
            reads_ref = reads_refs[k]
            job = batch_result['results'][k]
            result_package = job['result_package']
            if job['is_error']:
                n_error += 1
            else:
                n_success += 1
                output_info = result_package['result'][0]['output_info']
                ra_ref = output_info['upload_results']['obj_ref']
                alignment_items.append({
                        'ref': ra_ref,
                        'label': reads_ref.get(
                                'condition',
                                params.get('condition','unspecified'))
                })
                alignment_objs.append({'ref': ra_ref})

            if result_package['run_context']['location'] == 'local':
                ran_locally += 1
            if result_package['run_context']['location'] == 'njsw':
                ran_njsw += 1

        # Save the alignment set
        output_alignmentset_name = set_name + params['alignmentset_suffix']
        save_result = self.upload_alignment_set(
                        alignment_items,
                        output_alignmentset_name,
                        params['output_workspace'])

        # Report
        report_info = {'name': None, 'ref': None}
        input_ref = save_result['set_ref']

        #run qualimap
        qualimap_report = self.qualimap.run_bamqc({'input_ref': input_ref})
        qc_result_zip_info = qualimap_report['qc_result_zip_info']

        # create the report
        report_text = 'Ran on SampleSet or ReadsSet.\n\n'
        report_text += 'Created ReadsAlignmentSet: ' + str(output_alignmentset_name) + '\n\n'
        report_text += 'Total ReadsLibraries = ' + str(n_jobs) + '\n'
        report_text += '        Successful runs = ' + str(n_success) + '\n'
        report_text += '            Failed runs = ' + str(n_error) + '\n'
        report_text += '       Ran on main node = ' + str(ran_locally) + '\n'
        report_text += '   Ran on remote worker = ' + str(ran_njsw) + '\n\n'

        print('Report text=')
        print(report_text)

        kbr = KBaseReport(self.callback_url)
        report_info = kbr.create_extended_report({'message': report_text,
                                                  'objects_created': alignment_objs,
                                                  'report_object_name': 'kb_STAR_' + str(uuid.uuid4()),
                                                  'direct_html_link_index': 0,
                                                  'html_links': [{'shock_id': qc_result_zip_info['shock_id'],
                                                                  'name': qc_result_zip_info['index_html_file_name'],
                                                                  'label': qc_result_zip_info['name']}],
                                                  'workspace_name': params['output_workspace']
                                                  })

        result = {'alignmentset_ref': save_result['set_ref'],
                'output_info': batch_result,
                'alignment_objs': alignment_objs,
                'report_name': report_info['name'],
                'report_ref': report_info['ref']
        }

        return result


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
        print('Uploading completed alignment set...')
        alignment_set_data = {
            'description': 'Alignments using STAR, v.{}'.format(self.STAR_VERSION),
            'items': alignment_items
        }
        set_info = self.set_api_client.save_reads_alignment_set_v1({
            'workspace': ws_name,
            'output_object_name': alignmentset_name,
            'data': alignment_set_data
        })
        pprint(set_info)

        return set_info


    def determine_input_info(self, validated_params):
        ''' get info on the readsset_ref object and determine if we run once or run on a set '''
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

    def get_version_from_subactions(self, module_name, subactions):
        # go through each sub action looking for
        if not subactions:
            return 'release'  # default to release if we can't find anything
        for sa in subactions:
            if 'name' in sa:
                if sa['name'] == module_name:
                    # local-docker-image implies that we are running in kb-test, so return 'dev'
                    if sa['commit'] == 'local-docker-image':
                        return 'dev'
                    # to check that it is a valid hash, make sure it is the right
                    # length and made up of valid hash characters
                    if re.match('[a-fA-F0-9]{40}$', sa['commit']):
                        return sa['commit']
        # again, default to setting this to release
        return 'dev' #'release'


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


    def get_reads_refs(self, validated_params):
        return fetch_reads_refs_from_sampleset(validated_params[self.PARAM_IN_READS], self.workspace_url, self.callback_url, validated_params)


    # borrowed from kb_stringtie
    def _save_expression(self, output_dir, alignment_ref, workspace_name, gtf_file,
                         expression_suffix):
        """
        _save_expression: save Expression object to workspace
        """

        log('start saving Expression object')

        alignment_data_object = self.ws_client.get_objects2({'objects':
                                                     [{'ref': alignment_ref}]})['data'][0]

        alignment_name = alignment_data_object['info'][1]
        if re.match('.*_*[Aa]lignment', alignment_name):
            expression_obj_name = re.sub('_*[Aa]lignment', expression_suffix, alignment_name)
        else:
            expression_obj_name = alignment_name + expression_suffix
        destination_ref = workspace_name + '/' + expression_obj_name
        upload_expression_params = {'destination_ref': destination_ref,
                                    'source_dir': output_dir,
                                    'alignment_ref': alignment_ref,
                                    'tool_used': 'STAR',
                                    'tool_version': self.STAR_VERSION}

        expression_ref = self.eu.upload_expression(upload_expression_params)['obj_ref']

        return expression_ref


    def _save_expression_set(self, alignment_expression_map, alignment_set_ref, workspace_name,
                             expression_set_suffix):
        """
        _save_expression_set: save ExpressionSet object to workspace
        """

        log('start saving ExpressionSet object')

        items = []
        for alignment_expression in alignment_expression_map:
            items.append({'ref': alignment_expression.get('expression_obj_ref')})

        expression_set_data = {'description': 'ExpressionSet using StringTie', 
                               'items': items}

        alignment_set_data_object = self.ws_client.get_objects2({'objects':
                                                         [{'ref': alignment_set_ref}]})['data'][0]

        alignment_set_name = alignment_set_data_object['info'][1]
        if re.match('.*_*[Aa]lignment_*[Ss]et', alignment_set_name):
            expression_set_name = re.sub('_*[Aa]lignment_*[Ss]et',
                                         expression_set_suffix,
                                         alignment_set_name)
        else:
            expression_set_name = alignment_set_name + expression_set_suffix

        expression_set_save_params = {'data': expression_set_data,
                                      'workspace': workspace_name,
                                      'output_object_name': expression_set_name}

        save_result = self.set_api_client.save_expression_set_v1(expression_set_save_params)
        expression_set_ref = save_result['set_ref']

        return expression_set_ref


    def _process_alignment_object(self, params):
        """
        _process_alignment_object: process KBaseRNASeq.RNASeqAlignment type input object
        """
        log('start processing RNASeqAlignment object\nparams:\n{}'.format(json.dumps(params, 
                                                                                     indent=1)))
        alignment_ref = params.get('alignment_ref')

        alignment_object_info = self.ws_client.get_object_info3({"objects": 
                                                         [{"ref": alignment_ref}]}
                                                         )['infos'][0]
        alignment_name = alignment_object_info[1]

        output_directory = os.path.join(self.scratch, 
                                        alignment_name + '_' + str(int(time.time() * 100)))
        self._mkdir_p(output_directory)

        # input files
        params['input_file'] = self._get_input_file(alignment_ref)
        if not params.get('gtf_file'):
            params['gtf_file'] = self._get_gtf_file(alignment_ref, output_directory)
        else:
            shutil.copy(params.get('gtf_file'), output_directory)
        log('using {} as reference annotation file.'.format(params.get('gtf_file')))

        # output files
        self.output_transcripts = 'transcripts.gtf'
        params['output_transcripts'] = os.path.join(output_directory, self.output_transcripts)

        self.gene_abundances_file = 'genes.fpkm_tracking'
        params['gene_abundances_file'] = os.path.join(output_directory, self.gene_abundances_file)

        command = self._generate_command(params)
        self._run_command(command)

        if not params.get('merge'):
            expression_obj_ref = self._save_expression(output_directory,
                                                       alignment_ref,
                                                       params.get('workspace_name'),
                                                       params['gtf_file'],
                                                       params['expression_suffix'])
        else:
            log('skip generating expression object')
            expression_obj_ref = ''

        returnVal = {'output_directory': output_directory,
                     'expression_obj_ref': expression_obj_ref,
                     'alignment_ref': alignment_ref,
                     'annotation_file': params['gtf_file']}

        return returnVal


    def _process_alignment_set_object(self, params):
        """
        _process_alignment_set_object: process KBaseRNASeq.RNASeqAlignmentSet type input object
        """

        log('start processing AlignmentSet object\nparams:\n{}'.format(json.dumps(params, 
                                                                                  indent=1)))

        alignment_set_ref = params.get('alignment_set_ref')
        alignment_set_object = self.ws_client.get_objects2({'objects':
                                                    [{'ref': alignment_set_ref}]}
                                                    )['data'][0]

        alignment_set_info = alignment_set_object['info']
        alignment_set_data = alignment_set_object['data']

        alignment_set_type = alignment_set_info[2]

        mul_processor_params = []
        if re.match('KBaseRNASeq.RNASeqAlignmentSet-\d.\d', alignment_set_type):
            mapped_alignment_ids = alignment_set_data['mapped_alignments_ids']
            for i in mapped_alignment_ids:
                for sample_name, alignment_id in i.items():
                    alignment_upload_params = params.copy()
                    alignment_upload_params['alignment_ref'] = alignment_id
                    mul_processor_params.append(alignment_upload_params)
        elif re.match('KBaseSets.ReadsAlignmentSet-\d.\d', alignment_set_type):
            items = alignment_set_data['items']
            for item in items:
                alignment_ref = item['ref']
                alignment_upload_params = params.copy()
                alignment_upload_params['alignment_ref'] = alignment_ref
                mul_processor_params.append(alignment_upload_params)

        cpus = min(params.get('num_threads'), multiprocessing.cpu_count())
        pool = Pool(ncpus=cpus)
        log('running _process_alignment_object with {} cpus'.format(cpus))
        alignment_expression_map = pool.map(self._process_alignment_object, mul_processor_params)

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)

        for proc_alignment_return in alignment_expression_map:
            alignment_ref = proc_alignment_return.get('alignment_ref')
            alignment_name = self.ws_client.get_object_info([{"ref": alignment_ref}],
                                                     includeMetadata=None)[0][1]
            self._run_command('cp -R {} {}'.format(proc_alignment_return.get('output_directory'),
                                                   os.path.join(output_directory, 
                                                                alignment_name)))
        if not params.get('merge'):
            expression_obj_ref = self._save_expression_set(alignment_expression_map,
                                                           alignment_set_ref,
                                                           params.get('workspace_name'),
                                                           params['expression_set_suffix'])
        else:
            log('skip generating expression set object')
            expression_obj_ref = ''

        annotation_file_name = os.path.basename(alignment_expression_map[0]['annotation_file'])
        annotation_file_path = os.path.join(output_directory, 
                                            os.listdir(output_directory)[0], 
                                            annotation_file_name)

        returnVal = {'output_directory': output_directory,
                     'expression_obj_ref': expression_obj_ref,
                     'annotation_file': annotation_file_path}

        return returnVal


