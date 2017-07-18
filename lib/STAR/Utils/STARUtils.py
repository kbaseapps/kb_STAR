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

from STAR.Utils.Program_Runner import Program_Runner
from DataFileUtil.DataFileUtilClient import DataFileUtil
from Workspace.WorkspaceClient import Workspace as Workspace
from KBaseReport.KBaseReportClient import KBaseReport
#from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from KBParallel.KBParallelClient import KBParallel

from file_util import (
    valid_string,
    check_reference,
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
    PARAM_IN_GENOME = 'assembly_or_genome_ref'

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
        self.STAR_output = ''
        self.STAR_idx = ''
        self.star_runner = Program_Runner(self.STAR_BIN, self.scratch)
        self.parallel_runner = KBParallel(self.callback_url)


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
        if (params.get(self.PARAM_IN_OUTPUT_NAME, None) is None or
                not valid_string(params[self.PARAM_IN_OUTPUT_NAME])):
            raise ValueError("Parameter alignment output_name must be a valid Workspace object string, "
                      "not {}".format(params.get(self.PARAM_IN_OUTPUT_NAME, None)))
        if params.get(self.PARAM_IN_STARMODE, None) is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'
	else:
            if params[self.PARAM_IN_STARMODE] == "genomeGenerate":
		if params.get(self.PARAM_IN_GENOME, None) is None:
                    raise ValueError(self.PARAM_IN_GENOME +
				' parameter is required for generating genome index')
        params['sjdbGTFfile'] = self._get_genome_gtf_file(params[self.PARAM_IN_GENOME], self.STAR_idx)

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
            gfu_ret = gfu.genome_to_gff({'genome_ref': gnm_ref,
                                           'is_gtf': 1,
                                           'target_dir': gtf_file_dir})
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

        exitCode = self.star_runner.run(idx_cmd, self.scratch)

        return exitCode

    def _exec_mapping(self, params):
        log('Running STAR mapping with params:\n' + pformat(params))

        mp_cmd = self._construct_mapping_cmd(params)

        exitCode = self.star_runner.run(mp_cmd, self.scratch)

        return exitCode

    def _exec_star(self, params, rds_files, rds_name):
        # build the parameters
        params_idx = {
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
		self.STAR_IDX_DIR: self.STAR_idx,
                'genomeFastaFiles': params[self.PARAM_IN_FASTA_FILES]
        }
        aligndir = self.STAR_output
        if rds_name:
            aligndir = os.path.join(self.STAR_output, rds_name)
            self._mkdir_p(aligndir)
            print '\nSTAR output directory created: ' + aligndir
        params_mp = {
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'readFilesIn': rds_files,#params[self.PARAM_IN_READS_FILES],
		self.STAR_IDX_DIR: self.STAR_idx,
		'align_output': aligndir
        }

        if params.get('sjdbGTFfile', None) is not None:
            params_idx['sjdbGTFfile'] = params['sjdbGTFfile']
            params_mp['sjdbGTFfile'] = params['sjdbGTFfile']
        if params.get('sjdbOverhang', None) is not None :
            params_idx['sjdbOverhang'] = params['sjdbOverhang']
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
                ret = {'star_idx': self.STAR_idx, 'star_output': aligndir}
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
            "destination_ref": "{}/{}".format(input_params[self.PARAM_IN_WS], input_params[self.PARAM_IN_OUTPUT_NAME]),
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

    def _generate_extended_report(self, alignObj_refs, params):
        """
        generate_extended_report: generate a summary STAR report, including index files and alignment output files
        """
        log('Generating summary report...')

        created_objects = list()
        index_files = list()
        output_files = list()
        for oref in alignObj_refs:
            created_objects.append({
                "ref": oref,
                "description": "Reads {} aligned to Genome {}".format(params[self.PARAM_IN_READS], params[self.PARAM_IN_GENOME])
            })

        t0 = time.clock()
        #for out_dir in star_dirs:
            #index_files += self._get_output_file_list(self.STAR_IDX_DIR, out_dir['star_idx'])
        index_files = self._get_output_file_list(self.STAR_IDX_DIR, self.STAR_IDX_DIR)
        t1 = time.clock()
            #output_files += self._get_output_file_list(self.STAR_OUT_DIR, out_dir['star_output'])
        output_files = self._get_output_file_list(self.STAR_OUT_DIR, self.STAR_OUT_DIR)
        t2 = time.clock()
        pprint( "Zipping index files used {} seconds ".format(t1-t0))
        pprint( "Zipping output files used {} seconds ".format(t2-t1))

        report_params = {
              'message': 'Created a set of {} alignment(s) from the given sample set.'.format(len(alignObj_refs)),
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

        t3 = time.clock()
        pprint( "Creating report used {} seconds ".format(t3-t2))
        return {'report_name': report_info['name'], 'report_ref': report_info['ref']}

    def _generate_report(self, alignmentSet, params):
        """
        _generate_report: Creates a brief STAR report.
        """
        print("Creating STAR output report...in workspace " + params[self.PARAM_IN_WS])
        report_client = KBaseReport(self.callback_url, token=self.token)

        created_objects = list()
        for key, value in alignmentSet.iteritems():
            created_objects.append({
                'ref': key,
                'reads': value['readsName'],
                'alignment': value['name'],
                'description': 'Reads {} aligned to Genome {}'.format(value['readsName'], params[self.PARAM_IN_GENOME])
            })

        report_text = "Created one alignment set from the given reads set."
        report_info = report_client.create({
            'workspace_name': params[self.PARAM_IN_WS],
            "report": {
                "objects_created": created_objects,
                "text_message": report_text
            }
        })
        return {'report_name': report_info['name'], 'report_ref': report_info['ref']}


    def _get_reads_info(self, readsset_ref):
        reads_refs = list()
	if readsset_ref is not None:
            try:
		print("Fetching reads ref(s) from sampleset ref {}".format(readsset_ref))
		reads_refs = fetch_reads_refs_from_sampleset(readsset_ref, self.workspace_url, self.callback_url)
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
                if source_reads.get("condition", None) is not None:
                    ret_reads["condition"] = source_reads["condition"]
                else:
                    ret_reads["condition"] = None
                reads_info.append(ret_reads)
        return reads_info

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

    def _convert_params(self, input_params):
        """
        Convert input parameters with KBase ref format into STAR parameters,
        and add the advanced options.
        """
	params = {
            'runMode': 'genomeGenerate',
            'runThreadN': input_params[self.PARAM_IN_THREADN]
	}

	# STEP 1: Converting refs to file locations in the scratch area
        readsInfo = self._get_reads_info(input_params[self.PARAM_IN_READS])

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

        return {'input_parameters':params, 'reads_info': readsInfo}

    def _build_star_index(self, params):
        """
        Runs STAR in genomeGenerate mode to build the index files and directory for STAR mapping.
        It creates a directory as defined by self.STAR_IDX_DIR in the scratch area that houses the index files.
        """
        # check options and raise ValueError here as needed.
        idx_dir = self.STAR_idx
        try:
            os.mkdir(idx_dir)
        except OSError:
            print("Ignoring error for already existing {} directory".format(idx_dir))

        # build the indexing parameters
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

        return ret


    def _run_star_mapping(self, params, rds_files, rds_name):
        """
        Runs STAR in alignReads mode for STAR mapping.
        It creates a directory as defined by self.STAR_OUT_DIR with a subfolder named after the reads
        """
        # build the mapping parameters
        aligndir = self.STAR_output
        if rds_name:
            aligndir = os.path.join(self.STAR_output, rds_name)
            self._mkdir_p(aligndir)
            print '\nSTAR output directory created: ' + aligndir
        params_mp = {
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'readFilesIn': rds_files,#params[self.PARAM_IN_READS_FILES],
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
            ret = {'star_idx': self.STAR_idx, 'star_output': aligndir}

        return ret


    def run_single(self, reads_info, input_params):
        """
        Performs a single run of STAR against a single reads reference. The rest of the info
        is taken from the params dict - see the spec for details.
        """
        log('--->\nrunning STARUtil.run_single\n' +
            'params:\n{} on reads {}'.format(json.dumps(input_params, indent=1), reads_info))

        # 1. Get STAR index from genome.
        #    a. If it exists in cache, use that.
        #    b. Otherwise, build it
        idx_ret = self._build_star_index(input_params)
        alignments = dict()
        alignment_ref = None

        if idx_ret == 0:
            # 2. Fetch the reads file and make sure input params are correct.
            # if the reads ref came from a different sample set, then we need to drop that
            # reference inside the reads info object so it can be linked in the alignment-TODO
            #if reads_info["object_ref"] != input_params[self.PARAM_IN_READS]:
                #reads_info["readsset_ref"] = input_params[self.PARAM_IN_READS]
            # make sure condition info carries over if we have it

            if not "condition" in reads_info:
                reads_info["condition"] = input_params["condition"]

            rdsFiles = list()
            rdsName = ''
            ret_fwd = reads_info["file_fwd"]
            if ret_fwd is not None:
                rdsFiles.append(ret_fwd)
                rdsName = reads_info['file_name'].split('.')[0]
                if reads_info['file_rev'] is not None:
                    rdsFiles.append(rds['file_rev'])

            # 3. Finally all set, do the alignment and upload the output.
            star_mp_ret = self._run_star_mapping(input_params, rds_files, rds_name)
            if not isinstance(star_mp_ret, int):
                #print("Uploading STAR output object...")
                if input_params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
                    prefix = format(input_params[self.PARAM_IN_OUTFILE_PREFIX])
                    alignment_file = '{}Aligned.out.sam'.format(prefix)
                else:
                    alignment_file = 'Aligned.out.sam'

                alignment_file = os.path.join(star_ret['star_output'], alignment_file)

                # Upload the alignment
                alignment_ref = self.upload_STARalignment(input_params, reads_info, alignment_file)
                alignments[reads_info["object_ref"]] = {
                    "ref": alignment_ref,
                    "readsName": reads_info['file_name'],
                    "alignment_ame": input_params[self.PARAM_IN_OUTPUT_NAME]
                }

        return (alignments, alignment_ref)

    def run_batch(self, reads_infos, input_params):
        base_output_obj_name = input_params[self.PARAM_IN_OUTPUT_NAME]
        # build task list and send it to KBParallel
        tasks = list()
        for idx, rds_info in enumerate(reads_infos):
            single_param = dict(input_params)  # need a copy of the params
            single_param[self.PARAM_IN_OUTPUT_NAME] = "{}_{}".format(base_output_obj_name, idx)
            single_param["create_report"] = 0
            single_param[self.PARAM_IN_READS] = rds_info["object_ref"]
            if "condition" in rds_info:
                single_param["condition"] = rds_info["condition"]
            else:
                single_param["condition"] = "unspecified"

            tasks.append({
                "module_name": "kb_STAR",
                "function_name": "run_star",
                "version": "dev",
                "parameters": single_param
            })

        batch_run_params = {
            "tasks": tasks,
            "runner": "parallel",
            "concurrent_local_tasks": 3,
            "concurrent_njsw_tasks": 0,
            "max_retries": 2
        }
        parallel_runner = KBParallel(self.callback_url)
        results = parallel_runner.run_batch(batch_run_params)["results"]
        alignment_items = list()
        alignments = dict()
        for idx, result in enumerate(results):
            # idx of the result is the same as the idx of the inputs AND reads_infos
            if result["is_error"] != 0:
                raise RuntimeError("Failed a parallel run of STAR! {}".format(result["result_package"]["error"]))
            reads_ref = tasks[idx]["parameters"][self.PARAM_IN_READS]
            alignment_items.append({
                "ref": result["result_package"]["result"][0]["alignment_ref"],
                "label": reads_infos[idx].get(
                    "condition",
                    input_params.get("condition",
                               "unspecified condition"))
            })
            alignments[reads_ref] = {
                "ref": result["result_package"]["result"][0]["alignment_ref"],
                "name": tasks[idx]["parameters"][self.PARAM_IN_OUTPUT_NAME]
            }
        # build the final alignment set
        alignment_ref = self.upload_alignment_set(
            alignment_items, base_output_obj_name, input_params[self.PARAM_IN_WS]
        )
        return (alignments, alignment_ref)


    def run_star(self, input_params):
        """
        run_single: run the STAR app 
        """
        log('--->\nrunning STARUtil.run_star\n' +
            'params:\n{}'.format(json.dumps(input_params, indent=1)))

	# STEP 1: convert the input parameters (from refs to file paths, especially)
        params_ret = self._convert_params(input_params)
        params = params_ret.get('input_parameters', None)
        readsInfo = params_ret.get('reads_info', None)

	# STEP 2: Running star
        # Looping through for now, but later should implement the parallel processing here for all reads in readsInfo
        alignment_set = dict()
        star_out_dirs = list()
        for rds in readsInfo:
            rdsFiles = list()
            rdsName = ''
            ret_fwd = rds.get("file_fwd", None)
            if ret_fwd is not None:
                print("Done fetching FASTA file with name = {}".format(ret_fwd))
                rdsFiles.append(ret_fwd)
                rdsName = rds['file_name'].split('.')[0]
                if rds.get('file_rev', None) is not None:
                    rdsFiles.append(rds['file_rev'])

            star_ret = self._exec_star(params, rdsFiles, rdsName)
            if not isinstance(star_ret, int):
                #print("Uploading STAR output object...")
                if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
                    prefix = format(params[self.PARAM_IN_OUTFILE_PREFIX])
                    alignment_file = '{}Aligned.out.sam'.format(prefix)
                else:
                    alignment_file = 'Aligned.out.sam'

                alignment_file = os.path.join(star_ret['star_output'], alignment_file)

                # Upload the alignment
                alignment_ref = self.upload_STARalignment(input_params, rds, alignment_file)
                alignment_set[rds['object_ref']] = {
                        'ref': alignment_ref,
                        'readsName': rds['file_name'],
                        'name': params['output_name']
                }
                star_out_dirs.append(star_ret)

	# STEP 3: Generating report
        returnVal = {
                'alignment_ref': alignment_set
        }

        report_out = self._generate_extended_report(alignment_set, input_params)
        #report_out = self._generate_report(alignment_set, input_params)
        returnVal.update(report_out)

        return returnVal

