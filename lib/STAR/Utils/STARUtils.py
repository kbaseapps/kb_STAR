
import time
import json
import os
import re
import uuid
import errno
import subprocess
import shutil
import sys
import zipfile
from pprint import pprint, pformat

from gff_utils import GFFUtils
from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
#from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils

from file_util import (
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
    GENOME_ANN_GTF = 'genome_annotation_GTF'
    #STAR_DATA = '/kb/module/testReads'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_OUTPUT_NAME = 'output_name'
    PARAM_IN_FASTA_REFS = 'genomeFastaFile_refs'
    PARAM_IN_FASTA_FILES = 'genomeFastaFiles'
    PARAM_IN_OUTFILE_PREFIX = 'outFileNamePrefix'
    PARAM_IN_STARMODE = 'runMode'
    PARAM_IN_THREADN = 'runThreadN'
    PARAM_IN_READS_REFS = 'readFilesIn_refs'
    PARAM_IN_READS_FILES = 'readFilesIn'

    INVALID_WS_OBJ_NAME_RE = re.compile('[^\\w\\|._-]')
    INVALID_WS_NAME_RE = re.compile('[^\\w:._-]')

    PARAM_IN_READS = 'sampleset_ref'
    PARAM_IN_GENOME = 'genome_ref'

    def __init__(self, config):
        self.config = config
        self.workspace_url = config['workspace-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.ws_client = Workspace(self.workspace_url)
        self.au = AssemblyUtil(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.scratch = config['scratch']
        self.working_dir = self.scratch
        self.gff_utils = GFFUtils(config)


    def _mkdir_p(self, dir):
        """
        _mkdir_p: make directory for given path
        """
        log('Creating a new dir: ' + dir)
        if not dir:
            return
        if not os.path.exists(dir):
            os.makedirs(dir)

    def _process_params(self, params):
        """
        _process_params:
                checks params passed to run_star method and set default values
        """
        log('Start validating run_star parameters')

        if (self.PARAM_IN_WS not in params or
                not params[self.PARAM_IN_WS]):
            raise ValueError(self.PARAM_IN_WS + ' parameter is required')
        if (self.PARAM_IN_OUTPUT_NAME not in params or
                not params[self.PARAM_IN_OUTPUT_NAME]):
            raise ValueError(self.PARAM_IN_OUTPUT_NAME + ' parameter is required')
        if params.get(self.PARAM_IN_STARMODE, None) is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'
	else:
            if params[self.PARAM_IN_STARMODE] == "genomeGenerate":
		if params.get(self.PARAM_IN_GENOME, None) is None:
                    raise ValueError(self.PARAM_IN_GENOME +
				' parameter is required for generating genome index')
                else:
                    genome_ann_file_path = os.join(self.STAR_IDX_DIR, self.GENOME_ANN_GTF)
                    params['sjdbGTFfile'] = self._get_genome_gtf_file(params[self.PARAM_IN_GENOME], genome_ann_file_path)

        if (params.get(self.PARAM_IN_STARMODE, None) is not None and
		params[self.PARAM_IN_STARMODE] != "genomeGenerate"):
            if params.get(self.PARAM_IN_READS, None) is None:
		raise ValueError(self.PARAM_IN_READS +
				' parameter is required for reads mapping')

        if params.get(self.PARAM_IN_THREADN, None) is not None:
            if not isinstance(params[self.PARAM_IN_THREADN], int):
                raise ValueError(self.PARAM_IN_HASH_THREADN + ' must be of type int')
	else:
             params[self.PARAM_IN_THREADN] = 1

        if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is None:
             params[self.PARAM_IN_OUTFILE_PREFIX] = 'STARoutput_'
	elif params[self.PARAM_IN_OUTFILE_PREFIX].find('/') != -1:
            raise ValueError(self.PARAM_IN_OUTFILE_PREFIX + ' cannot contain subfolder(s).')

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


    def _get_genome_gtf_file(self, gnm_ref, gtf_file_path):
        """
        Get data from genome object ref and return the GTF filename (with path)
        for STAR indexing and mapping.
        STAR uses the reference annotation to guide assembly and for creating alignment
        """
        if self.gff_utils.convert_genome_to_gtf(gnm_ref, gtf_file_path) == 0:
            return gtf_file_path
        else:
            return None


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

	# appending the standard optional inputs
        if params.get('sjdbGTFfile', None) is not None:
            idx_cmd.append('--sjdbGTFfile')
            idx_cmd.append(params['sjdbGTFfile'])
        if (params.get('sjdbOverhang', None) is not None
		and params['sjdbOverhang'] > 0):
            idx_cmd.append('--sjdbOverhang')
            idx_cmd.append(str(params['sjdbOverhang']))

        # STEP 2: appending the advanced optional inputs--TODO

        # STEP 3: return idx_cmd
        print ('STAR indexing CMD:')
        print ' '.join(idx_cmd)
        return idx_cmd

    def _construct_mapping_cmd(self, params):
	if params.get(self.PARAM_IN_STARMODE, None) is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'

        # STEP 1: set the working folder housing the STAR output results as well as the reads info
        star_out_dir = ''
	if params.get(self.STAR_OUT_DIR, None) is None:
            star_out_dir = self.scratch
	else:
            star_out_dir = params[self.STAR_OUT_DIR]

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

	if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
            mp_cmd.append('--' + self.PARAM_IN_OUTFILE_PREFIX)
            mp_cmd.append(os.path.join(star_out_dir, params[self.PARAM_IN_OUTFILE_PREFIX]))

        # STEP 3: appending the advanced optional inputs
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
        STAR_cmd = self._construct_indexing_cmd(params)

        p = subprocess.Popen(STAR_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()

        log('Return code: ' + str(retcode))
        if p.returncode != 0:
            raise ValueError('Error running STAR index generating, return code: ' + str(retcode) + '\n')

        return p.returncode

    def _exec_mapping(self, params):
        log('Running STAR mapping with params:\n' + pformat(params))
        STAR_cmd = self._construct_mapping_cmd(params)
        p = subprocess.Popen(STAR_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()
        log('Return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running STAR mapping, return code: ' + str(p.returncode) + '\n')

        return p.returncode


    def _exec_star(self, params):
        # build the parameters
        params_idx = {
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
		self.STAR_IDX_DIR: idxdir,
                'genomeFastaFiles': params[self.PARAM_IN_FASTA_FILES]
        }
        params_mp = {
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'readFilesIn': params[self.PARAM_IN_READS_FILES],
		self.STAR_IDX_DIR: idxdir,
		self.STAR_OUT_DIR: outdir,
		'outFileNamePrefix': params[self.PARAM_IN_OUTFILE_PREFIX]
        }

        if params.get('sjdbGTFfile', None) is not None:
            params_idx['sjdbGTFfile'] = params['sjdbGTFfile']
            params_mp['sjdbGTFfile'] = params['sjdbGTFfile']
        if params.get('sjdbOverhang', None) is not None :
            params_idx['sjdbOverhang'] = params['sjdbOverhang']
            params_mp['sjdbOverhang'] = params['sjdbOverhang']

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

        ret = 1
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
            ret = 1
            params_mp[self.PARAM_IN_STARMODE] = 'alignReads'
            try:
                ret = self._exec_mapping(params_mp)
                while( ret != 0 ):
                    time.sleep(1)
            except ValueError as emp:
                log('STAR mapping raised error:\n')
                pprint(emp)
            else:#no exception raised by STAR mapping and STAR returns 0, then move to saving and reporting  
                ret = {'STAR_idx': idxdir, 'STAR_output': outdir}
        return ret

    def upload_STARalignment(self, input_params, reads_info, alignment_file):
        """
        Uploads the alignment file + metadata.
        Returns the STAR alignment reference.
        """
        aligner_opts = dict()
        for k in input_params:
            aligner_opts[k] = str(input_params[k])
        pprint(reads_info)
        align_upload_params = {
            "destination_ref": "{}/{}".format(input_params['workspace_name'], input_params[self.PARAM_IN_OUTPUT_NAME]),
            "file_path": alignment_file,
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
        alignment_ref = ra_util.upload_alignment(align_upload_params)["obj_ref"]
        print("STAR alignment uploaded as object {}".format(alignment_ref))
        return alignment_ref

    # borrowed from kb_stringtie and modified for STAR
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

    def _generate_extended_report(self, obj_ref, params, star_ret):
        """
        generate_report: generate a summary STAR report, including index files and alignment output files
        """
        log('Generating summary report...')

        created_objects = list()
        created_objects.append({
            "ref": obj_ref,
            "description": "Reads {} aligned to Genome {}".format(params[self.PARAM_IN_READS], params[self.PARAM_IN_GENOME])
        })
        t0 = time.clock()
	index_files = self._get_output_file_list(self.STAR_IDX_DIR, star_ret['STAR_idx'])
        t1 = time.clock()
        output_files = self._get_output_file_list(self.STAR_OUT_DIR, star_ret['STAR_output'])
        t2 = time.clock()
        pprint( "Zipping index files used {} seconds ".format(t1-t0))
        pprint( "Zipping output files used {} seconds ".format(t2-t1))

        report_params = {
              'message': 'Created one alignment from the given sample set.',
              'workspace_name': params.get('workspace_name'),
              'objects_created': created_objects,
              'file_links': index_files + output_files,
              'direct_html_link_index': 0,
              'html_window_height': 0,#366,
              'summary_window_height': 0,#366,
              'report_object_name': 'kb_star_report_' + str(uuid.uuid4())
	}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        report_info = kbase_report_client.create_extended_report(report_params)

        t3 = time.clock()
        pprint( "Creating report used {} seconds ".format(t3-t2))
        return {'report_name': report_info['name'], 'report_ref': report_info['ref']}

    def _generate_report(self, obj_ref, params):
        """
        Creates a brief STAR report.
        """
        report_client = KBaseReport(self.callback_url, token=self.token)
        report_text = None
        created_objects = list()
        created_objects.append({
            "ref": obj_ref,
            "description": "Reads {} aligned to Genome {}".format(params[self.PARAM_IN_READS], params[self.PARAM_IN_GENOME])
        })

        report_text = "Created one alignment from the given sample set."
        report_info = report_client.create({
            "workspace_name": params[self.PARAM_IN_WS],
            "report": {
                "objects_created": created_objects,
                "text_message": report_text
            }
        })
        return {'report_name': report_info['name'], 'report_ref': report_info['ref']}

    def _convert_params(self, input_params):
        """
        Convert input parameters with KBase ref format into STAR parameters,
        and add the advanced options.
        """
	params = {
            'runMode': 'genomeGenerate',
            'runThreadN': input_params[self.PARAM_IN_THREADN],
            'outFileNamePrefix': input_params[self.PARAM_IN_OUTFILE_PREFIX]
	}

	# STEP 1: Converting refs to file locations in the scratch area
        smplset_ref = input_params.get(self.PARAM_IN_READS, None)
        reads_refs = list()
	if smplset_ref is not None:
            try:
		print("Fetching reads ref(s) from sampleset ref {}".format(smplset_ref))
		reads_refs = fetch_reads_refs_from_sampleset(smplset_ref, self.workspace_url, self.callback_url)
		print("Done fetching reads ref(s)!")
            except ValueError:
		print("Incorrect object type for fetching reads ref(s)!")
		raise

        readsInfo = list()
	readsFiles = list()
	for source_reads in reads_refs:
            try:
                print("Fetching FASTA file from reads reference {}".format(source_reads['ref']))
                ret_reads = fetch_reads_from_reference(source_reads['ref'], self.callback_url)
                ret_fwd = ret_reads.get("file_fwd", None)
		if ret_fwd is not None:
                    print("Done fetching FASTA file with name = {}".format(ret_fwd))
                    readsFiles.append(ret_reads["file_fwd"])
                    if ret_reads.get("file_rev", None) is not None:
                        readsFiles.append(ret_reads["file_rev"])
            except ValueError:
                print("Incorrect object type for fetching a FASTA file!")
                raise

            if ret_reads.get("file_fwd", None) is None:
                raise RuntimeError("FASTA file fetched from reads ref {} doesn't seem to exist!".format(source_reads['ref']))
            else:
                if source_reads.get("condition", None) is not None:
                    ret_reads["condition"] = source_reads["condition"]
                else:
                    ret_reads["condition"] = "N/A"
                readsInfo.append(ret_reads)

	params[self.PARAM_IN_READS_FILES] = readsFiles

        params[self.PARAM_IN_FASTA_FILES] = list()
        genome_ref = input_params.get(self.PARAM_IN_GENOME, None)
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
		params[self.PARAM_IN_FASTA_FILES].append(genome_fasta_file["path"])

        # STEP 2: Add advanced options from input_params to params
        sjdbGTFfile = input_params.get("sjdbGTFfile", None)
	if sjdbGTFfile is not None:
            params['sjdbGTFfile'] = sjdbGTFfile
            if input_params.get('sjdbOverhang', None) is not None :
                params['sjdbOverhang'] = input_params['sjdbOverhang']
            else:
                params['sjdbOverhang'] = 100

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

        return params


    def run_star(self, input_params):
        """
        run_star: run the STAR app
        """
        log('--->\nrunning STARUtil.run_star\n' +
            'params:\n{}'.format(json.dumps(input_params, indent=1)))

        # STEP 0: creating the directories for STAR
        outdir = os.path.join(self.scratch, self.STAR_OUT_DIR)
        self._mkdir_p(outdir)
        idxdir = os.path.join(self.scratch, self.STAR_IDX_DIR)
        self._mkdir_p(idxdir)

	# STEP 1: preprocessing the input parameters
        input_params = self._process_params(input_params)
        params = self._convert_params(input_params)

	# STEP 2: Running star
	star_ret = self._exec_star(params)

	# STEP 3: Uploading the alignment and generating report
        if not isinstance(star_ret, int):
            #print("Uploading STAR output object...")
            alignment_file = "{}Aligned.out.sam".format(input_params[self.PARAM_IN_OUTFILE_PREFIX])
            alignment_file = os.path.join(star_ret['STAR_output'], alignment_file)

            # Upload the alignment with ONLY the first reads_ref for now
            alignment_ref = self.upload_STARalignment(input_params, readsInfo[0], alignment_file)

            returnVal = {
                'output_folder': star_ret['STAR_output'],
                'alignment_ref': alignment_ref
            }

            #print("Creating STAR output report...")
            #report_out = self._generate_extended_report(alignment_ref, input_params, star_ret)
            report_out = self._generate_report(alignment_ref, input_params)

            returnVal.update(report_out)
        else:
            print("STAR failed with error!!!")
            returnVal = {
                'output_folder': 'star_raised an error',
                'alignment_ref': None
            }
        return returnVal

