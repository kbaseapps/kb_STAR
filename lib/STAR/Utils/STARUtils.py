import time
import json
import os
import re
import copy
import uuid
import errno
import subprocess
import shutil
import sys
import zipfile
from pprint import pprint, pformat

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


class STARUtil:
    STAR_VERSION = 'STAR 2.5.3a'
    STAR_BIN = '/kb/deployment/bin/STAR'
    STAR_IDX_DIR = 'STAR_Genome_index'
    STAR_OUT_DIR = 'STAR_Output'
    GENOME_ANN_GTF = 'genome_annotation.gtf'
    #STAR_DATA = '/kb/module/testReads'
    PARAM_IN_WS = 'output_workspace'
    PARAM_IN_OUTPUT_NAME = 'output_name'
    PARAM_IN_FASTA_REFS = 'genomeFastaFile_refs'
    PARAM_IN_FASTA_FILES = 'genomeFastaFiles'
    PARAM_IN_OUTFILE_PREFIX = 'outFileNamePrefix'
    PARAM_IN_STARMODE = 'runMode'
    PARAM_IN_THREADN = 'runThreadN'
    PARAM_IN_READS_INFO = 'readsInfo'
    PARAM_IN_READS_FILES = 'readFilesIn'

    INVALID_WS_OBJ_NAME_RE = re.compile('[^\\w\\|._-]')
    INVALID_WS_NAME_RE = re.compile('[^\\w:._-]')

    PARAM_IN_READS = 'readsset_ref'
    PARAM_IN_GENOME = 'genome_ref'

    def __init__(self, config, provenance):
        self.config = config
        self.workspace_url = config['workspace-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.srv_wiz_url = config['srv-wiz-url']
        self.ws_client = Workspace(self.workspace_url)
        self.au = AssemblyUtil(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.scratch = config['scratch']
        self.working_dir = self.scratch
        self.STAR_output = ''
        self.STAR_idx = ''
        self.prog_runner = Program_Runner(self.STAR_BIN, self.scratch)
        self.provenance = provenance

        # from the provenance, extract out the version to run by exact hash if possible
        self.my_version = 'release'
        if len(provenance) > 0:
            if 'subactions' in provenance[0]:
                self.my_version = self.get_version_from_subactions('kb_STAR', provenance[0]['subactions'])
        print('Running kb_STAR version = ' + self.my_version)

        self.parallel_runner = KBParallel(self.callback_url)
        self.qualimap = kb_QualiMap(self.callback_url, service_ver="dev")


    def _mkdir_p(self, dir):
        """
        _mkdir_p: make directory for given path
        """
        log('Creating a new dir: ' + dir)
        if not dir:
            return
        if not os.path.exists(dir):
            os.makedirs(dir)

    def process_params(self, params):
        """
        process_params:
                checks params passed to run_star method and set default values
        """
        log('Start validating run_star parameters')

        # STEP 0: creating the directories for STAR
        # the index directory
        idxdir = os.path.join(self.scratch, self.STAR_IDX_DIR)
        self._mkdir_p(idxdir)
        self.STAR_idx = idxdir
        # the output directory
        outdir = os.path.join(self.scratch, self.STAR_OUT_DIR)
        self._mkdir_p(outdir)
        self.STAR_output = outdir

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
	else:
            if params[self.PARAM_IN_STARMODE] == "genomeGenerate":
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

    def _exec_star_pipeline(self, params, rds_files, rds_name):
        # build the parameters
        params_idx = self._get_indexing_params(params)
        params_mp = self._get_mapping_params(params, rds_files, rds_name)

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
                ret = {'star_idx': self.STAR_idx, 'star_output': params_mp.get('align_output')}

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

        ra_util = ReadsAlignmentUtils(self.callback_url, service_ver="dev")
        rau_upload_ret = ra_util.upload_alignment(align_upload_params)
        alignment_ref = rau_upload_ret["obj_ref"]
        print("STAR alignment uploaded as object {}".format(alignment_ref))
        return rau_upload_ret

    # borrowed from kb_stringtie and modified for STAR, for cases when user wants the files
    def _get_output_file_list(self, out_filename, output_dir):
        """
        _get_output_file_list: zip result files and generate file_links for report
        """
	if os.path.splitext(out_filename) != '.zip':
            out_filename += '.zip'
        log('Start packing result files to ' + out_filename)

        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file = os.path.join(output_directory, out_filename)

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(output_dir):
                for file in files:
                    zip_file.write(os.path.join(root, file),os.path.join(os.path.basename(root), file))

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by STAR App in ' + output_dir})

        return output_files


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


    def _generate_extended_report(self, alignObjs, params):
        """
        generate_extended_report: generate a summary STAR report, including index files and alignment output files
        """
        log('Generating summary report...')

        genomeName = self.get_name_from_obj_info(self.get_obj_infos(params[self.PARAM_IN_GENOME])[0])
        created_objects = list()
        index_files = list()
        output_files = list()
        for key, value in alignObjs.iteritems():
            created_objects.append({
                'ref': key,
                'reads': value['readsName'],
                'alignment': value['alignment_name'],
                'description': 'Reads {} aligned to Genome {}'.format(value['readsName'], genomeName)
            })

        index_files = self._get_output_file_list(self.STAR_IDX_DIR, self.STAR_IDX_DIR)
        output_files = self._get_output_file_list(self.STAR_OUT_DIR, self.STAR_OUT_DIR)

        report_params = {
              'message': 'Created a set of {} alignment(s) from the given sample set.'.format(len(alignObjs)),
              'workspace_name': params[self.PARAM_IN_WS],
              'objects_created': created_objects,
              'file_links': index_files + output_files,
              'direct_html_link_index': 0,
              'html_window_height': 366,
              'summary_window_height': 0,#366,
              'report_object_name': 'kb_star_report_' + str(uuid.uuid4())
	}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        report_info = kbase_report_client.create_extended_report(report_params)

        return report_info #{'report_name': report_info['name'], 'report_ref': report_info['ref']}


    def _generate_report(self, alignmentSet, params):
        """
        _generate_report: Creates a brief STAR report.
        """
        print("Creating STAR output report...in workspace " + params[self.PARAM_IN_WS])

        alignmentName = ''
        genomeName = self.get_name_from_obj_info(self.get_obj_infos(params[self.PARAM_IN_GENOME])[0])
        created_objects = list()
        for key, value in alignmentSet.iteritems():
            created_objects.append({
                'ref': key,
                'reads': value['readsName'],
                'alignment': value['alignment_name'],
                'description': 'Reads {} aligned to Genome {}'.format(value['readsName'], genomeName)
            })
            alignmentName = value['alignment_name'] #for now only, will implement the alignmentSet name later

        report_text = 'Created ReadsAlignment: ' + alignmentName + '\n'

        report_client = KBaseReport(self.callback_url, token=self.token)
        report_info = report_client.create({
            'workspace_name': params[self.PARAM_IN_WS],
            "report": {
                "objects_created": created_objects,
                "text_message": report_text
            }
        })
        return report_info #{'report_name': report_info['name'], 'report_ref': report_info['ref']}


    def get_reads(self, params):
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
	# STEP 1: Converting refs to file locations in the scratch area
        reads = self.get_reads(input_params)

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


    def run_star_indexing(self, input_params):
        """
        Runs STAR in genomeGenerate mode to build the index files and directory for STAR mapping.
        It creates a directory as defined by self.STAR_IDX_DIR in the scratch area that houses the index files.
        """
        genome_params = copy.deepcopy(input_params)

        # GTF file -create only once as the index is being generated
        genome_params['sjdbGTFfile'] = self._get_genome_gtf_file(input_params[self.PARAM_IN_GENOME], self.STAR_idx)

        idx_dir = self.STAR_idx
        try:
            os.mkdir(idx_dir)
        except OSError:
            print("Ignoring error for already existing {} directory".format(idx_dir))

        # build the indexing parameters
        params_idx = self._get_indexing_params(genome_params)

        ret = 1
        try:
            if genome_params[self.PARAM_IN_STARMODE]=='genomeGenerate':
                ret = self._exec_indexing(params_idx)
            else:
		ret = 0
            while( ret != 0 ):
                time.sleep(1)
        except ValueError as eidx:
            log('STAR genome indexing raised error:\n')
            pprint(eidx)
        else:
            ret = 0

        return ret


    def _get_indexing_params(self, params):
        params_idx = {
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
		self.STAR_IDX_DIR: self.STAR_idx,
                'genomeFastaFiles': params[self.PARAM_IN_FASTA_FILES]
        }
        if params.get('sjdbGTFfile', None) is not None:
            params_idx['sjdbGTFfile'] = params['sjdbGTFfile']
        if params.get('sjdbOverhang', None) is not None :
            params_idx['sjdbOverhang'] = params['sjdbOverhang']

        return params_idx


    def _get_mapping_params(self, params, rds_files, rds_name):
        ''' build the mapping parameters'''
        aligndir = self.STAR_output
        if rds_name:
            aligndir = os.path.join(self.STAR_output, rds_name)
            self._mkdir_p(aligndir)
            print '\n**********STAR output directory created: ' + aligndir

        params_mp = {
                'runMode': 'alignReads',#params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'readFilesIn': rds_files,
		self.STAR_IDX_DIR: self.STAR_idx,
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


    def run_star_mapping(self, params, rds_files, rds_name):
        """
        Runs STAR in alignReads mode for STAR mapping.
        It creates a directory as defined by self.STAR_OUT_DIR with a subfolder named after the reads
        """
        params_mp = self._get_mapping_params(params, rds_files, rds_name)

        retVal = {}
        params_mp[self.PARAM_IN_STARMODE] = 'alignReads'
        try:
            ret = self._exec_mapping(params_mp)
            while( ret != 0 ):
                time.sleep(1)
        except ValueError as emp:
            log('STAR mapping raised error:\n')
            pprint(emp)
            retVal = {'star_idx': self.STAR_idx, 'star_output': None}
        else:#no exception raised by STAR mapping and STAR returns 0, then move to saving and reporting  
            retVal = {'star_idx': self.STAR_idx, 'star_output': params_mp.get('align_output')}

        return retVal


    def star_run_single(self, validated_params):
        """
        Performs a single run of STAR against a single reads reference. The rest of the info
        is taken from the params dict - see the spec for details.
        """
        log('--->\nrunning STARUtil.run_single\n' +
                'params:\n{}'.format(json.dumps(validated_params, indent=1)))

	# convert the input parameters (from refs to file paths, especially)
        params_ret = self.convert_params(validated_params)
        input_params = params_ret.get('input_parameters', None)
        reads = params_ret.get('reads', None)
        reads_ref = reads.get('readsRefs', None)[0]
        reads_info = reads.get('readsInfo', None)[0]
        rds_name = reads_ref['alignment_output_name'].replace(input_params['alignment_suffix'], '')

        alignment_objs = list()
        alignment_ref = None
        singlerun_output_info = {}
        report_info = {'name': None, 'ref': None}

        if not 'condition' in reads_info:
            reads_info['condition'] = input_params['condition']

        rds_files = list()
        ret_fwd = reads_info["file_fwd"]
        if ret_fwd is not None:
            rds_files.append(ret_fwd)
            if reads_info.get('file_rev', None) is not None:
                rds_files.append(reads_info['file_rev'])

        # After all is set, do the alignment and upload the output.
        star_mp_ret = self.run_star_mapping(input_params, rds_files, rds_name)

        if star_mp_ret.get('star_output', None) is not None:
            if input_params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
                prefix = format(input_params[self.PARAM_IN_OUTFILE_PREFIX])
                output_sam_file = '{}Aligned.out.sam'.format(prefix)
            else:
                output_sam_file = 'Aligned.out.sam'
            output_sam_file = os.path.join(star_mp_ret['star_output'], output_sam_file)

            # Upload the alignment
            #print("Uploading STAR output object...")
            upload_results = self.upload_STARalignment(input_params, reads, output_sam_file)
            alignment_ref = upload_results['obj_ref']
            alignment_obj = {
                "ref": alignment_ref,
                "name": reads_ref['alignment_output_name']
            }
            alignment_objs.append({
                'reads_ref': reads_ref['ref'],
                'AlignmentObj': alignment_obj
            })

            singlerun_output_info['output_dir'] = star_mp_ret['star_output']
            singlerun_output_info['output_sam_file'] = output_sam_file

            singlerun_output_info['upload_results'] = upload_results

            if input_params.get("create_report", 0) == 1:
                report_info = self.generate_report_for_single_run(singlerun_output_info, input_params)

        if ret_fwd is not None:
            os.remove(ret_fwd)
            if reads_info.get('file_rev', None) is not None:
                os.remove(reads_info["file_rev"])

        return {'alignmentset_ref': None,
                'output_info': singlerun_output_info,
                'alignment_objs': alignment_objs,
                'report_name': report_info['name'],
                'report_ref': report_info['ref']
        }


    def star_run_batch(self, validated_params):
        ''' convert the input parameters (from refs to file paths, especially)'''
        reads_refs = fetch_reads_refs_from_sampleset(validated_params[self.PARAM_IN_READS], self.workspace_url, self.callback_url, validated_params)

        # build task list and send it to KBParallel
        tasks = []
        for r in reads_refs:
            tasks.append(self.build_single_execution_task(r['ref'], validated_params))

        batch_run_params = {'tasks': tasks,
                            'runner': 'parallel',
                            'max_retries': 2}

        if validated_params.get('concurrent_local_tasks', None) is not None:
                batch_run_params['concurrent_local_tasks'] = validated_params['concurrent_local_tasks']
        if validated_params.get('concurrent_njsw_tasks', None) is not None:
                batch_run_params['concurrent_njsw_tasks'] = validated_params['concurrent_njsw_tasks']

        results = self.parallel_runner.run_batch(batch_run_params)
        print('Batch run results=')
        pprint(results)

        batch_result = self.process_batch_result(results, validated_params, reads_refs)

        return batch_result


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
        alignment_items = list()
        alignment_objs = list()
        alignments = dict()

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
                                "condition",
                                params.get("condition","unspecified"))
                })
                alignments[reads_ref] = result_package["result"][0]["alignment_objs"][reads_ref]
                alignment_objs += result_package["result"][0]["alignment_objs"]

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
                                                  'objects_created': objects_created,
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
        print("Uploading completed alignment set")
        alignment_set_data = {
            "description": "Alignments using STAR, v.{}".format(self.STAR_VERSION),
            "items": alignment_items
        }
        set_api = SetAPI(self.srv_wiz_url)
        set_info = set_api.save_reads_alignment_set_v1({
            "workspace": ws_name,
            "output_object_name": alignmentset_name,
            "data": alignment_set_data
        })
        pprint(set_info)

        return set_info


    def star_run_sequential(self, reads, input_params, input_obj_info):
        """
        run_star_sequential: run the STAR app on each reads one by one
	(Running star indexing and then mapping)
        """
        log('--->\nrunning STARUtil.run_star\n' +
            'params:\n{}'.format(json.dumps(input_params, indent=1)))

        alignment_objs = list()
        output_sam_file = ''
        upload_out = {}

        reads_info = reads['readsInfo']
        for rds in reads_info:
            rdsFiles = list()
            rdsName = input_obj_info['info'][1]
            ret_fwd = rds.get("file_fwd", None)
            if ret_fwd is not None:
                print("Done fetching FASTA file with name = {}".format(ret_fwd))
                rdsFiles.append(ret_fwd)
                if rds.get('file_rev', None) is not None:
                    rdsFiles.append(rds['file_rev'])

            star_ret = self._exec_star_pipeline(inpupt_params, rdsFiles, rdsName)
            if star_ret.get('star_output', None) is not None:
                if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
                    prefix = format(params[self.PARAM_IN_OUTFILE_PREFIX])
                    output_sam_file = '{}Aligned.out.sam'.format(prefix)
                else:
                    output_sam_file = 'Aligned.out.sam'

                output_sam_file = os.path.join(star_ret['star_output'], output_sam_file)

                #print("Uploading STAR output object...")
                # Upload the alignment
                upload_out = self.upload_STARalignment(input_params, reads, output_sam_file)
                alignment_ref = upload_out['ob_ref']
                alignObj = {'ref': alignment_ref,
                            'name': '{}_starAligned'.format(rds['name'])
                }
                alignment_objs.append({
                        'reads_ref': rds['object_ref'],
                        'AlignmentObj': alignObj
                })

            if ret_fwd is not None:
                os.remove(ret_fwd)
                if rds.get('file_rev', None) is not None:
                    os.remove(rds["file_rev"])
	# STEP 3: Generating report
        returnVal = {
                'output_info': {'upload_results': upload_out},
                'alignment_objs': alignment_objs,
        }

        report_info = self._generate_extended_report(alignment_set, input_params)
        report_out = {'report_info': report_info}
        returnVal.update(report_out)

        return returnVal


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
        ws = self.ws_client
        obj_ids = list()
        for ref in ref_list:
            obj_ids.append({"ref": ref})
        info = ws.get_object_info3({"objects": obj_ids})
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
        return 'release'

