
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

#from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport
from ReadsUtils.ReadsUtilsClient import ReadsUtils
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
    STAR_IDX = '/kb/module/STAR_Genome_index/'
    #STAR_DATA = '/kb/module/testReads'
    PARAM_IN_WS = 'workspace_name'
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

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise


    def _process_params(self, params):
        """
        _process_params:
                checks params passed to run_star method and set default values
        """
        log('Start validating run_star parameters')

        if (self.PARAM_IN_WS not in params or
                not params[self.PARAM_IN_WS]):
            raise ValueError(self.PARAM_IN_WS + ' parameter is required')

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
             params[self.PARAM_IN_THREADN] = 1

        if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is None:
             params[self.PARAM_IN_OUTFILE_PREFIX] = 'STARoutput_'
	elif params[self.PARAM_IN_OUTFILE_PREFIX].find('/') != -1:
            raise ValueError(self.PARAM_IN_OUTFILE_PREFIX + ' cannot contain subfolder(s).')

	return params


    def _construct_indexing_cmd(self, params):
        # STEP 1: get the workspace and its folder where the genome indices are stored
        wsname = params[self.PARAM_IN_WS]

	# STEP 2: construct the command for running `STAR --runMode genomeGenerate...`
        idx_cmd = [self.STAR_BIN]
	idx_cmd.append('--genomeDir')
	idx_cmd.append(self.STAR_IDX)
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

        # appending the advanced optional inputs--TODO

        # STEP 3: return idx_cmd
        print ('STAR indexing CMD:')
        print ' '.join(idx_cmd)
        return idx_cmd

    def _construct_mapping_cmd(self, params):
        # STEP 1: get the working folder housing the STAR results as well as the reads info
        wsname = params[self.PARAM_IN_WS]
	if params.get(self.PARAM_IN_STARMODE, None) is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'

        star_out_dir = ''
	if params.get('star_output_dir', None) is None:
            star_out_dir = self.scratch
	else:
            star_out_dir = os.path.join(self.scratch, params['star_output_dir'])
        self._mkdir_p(star_out_dir)

        # STEP 2: construct the command for running STAR mapping
        mp_cmd = [self.STAR_BIN]
	mp_cmd.append('--genomeDir')
	mp_cmd.append(self.STAR_IDX)
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
        # appending the advanced optional inputs--TODO
        quant_modes = ["TranscriptomeSAM", "GeneCounts"]
        if (params.get('quantMode', None) is not None
                and params.get('quantMode', None) in quant_modes):
            mp_cmd.append('--quantMode')
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

        # STEP 3 return mp_cmd
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

    def __init__(self, config):
	self.workspace_url = config['workspace-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shock_url = config['shock-url']
        #self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.scratch = config['scratch']
        self.working_dir = self.scratch

    def _exec_star(self, params, star_outdir):
        outdir = os.path.join(self.scratch, star_outdir)
        self._mkdir_p(outdir)

        # build the parameters
        params_idx = {
                'workspace_name': params[self.PARAM_IN_WS],
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'genomeFastaFiles': params[self.PARAM_IN_FASTA_FILES]
        }

        if params.get('sjdbGTFfile', None) is not None:
            params_idx['sjdbGTFfile'] = params['sjdbGTFfile']
        if params.get('sjdbOverhang', None) is not None :
            params_idx['sjdbOverhang'] = params['sjdbOverhang']

        params_mp = {
                'workspace_name': params[self.PARAM_IN_WS],
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'readFilesIn': params[self.PARAM_IN_READS_FILES],
		'star_output_dir': outdir,
		'outFileNamePrefix': params[self.PARAM_IN_OUTFILE_PREFIX]
        }
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
                ret = outdir
        return ret

    def upload_STARalignment(self, input_params, reads_info, alignment_file):
        """
        Uploads the alignment file + metadata.
        Returns the STAR alignment reference.
        """
        aligner_opts = dict()
        for k in input_params:
            aligner_opts[k] = str(input_params[k])

        align_upload_params = {
            "destination_ref": "{}/{}".format(input_params['workspace_name'], input_params['alignment_name']),
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

    def _generate_extended_report(self, obj_ref, output_dir, params):
        """
        generate_report: generate summary report
        """
        log('Generating report')

        output_files = self._generate_output_file_list(output_dir)

        output_html_files = self._generate_html_report(output_dir,
                                                       params.get('assembly_ref'),
                                                       obj_ref,
                                                       params.get('out_header'))

        report_params = {
              'message': '',
              'workspace_name': params.get('workspace_name'),
              'file_links': output_files,
              'html_links': output_html_files,
              'direct_html_link_index': 0,
              'html_window_height': 266,
              'report_object_name': 'kb_star_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_report(self, params, reads_refs, alignments, alignment_set=None):
        """
        Builds and uploads the STAR report.
        """
        report_client = KBaseReport(self.callback_url)
        report_text = None
        created_objects = list()
        for k in alignments:
            created_objects.append({
                "ref": alignments[k],
                "description": "Reads {} aligned to Genome {}".format(k, params["genome_ref"])
            })
        if alignment_set is not None:
            created_objects.append({
                "ref": alignment_set,
                "description": "Set of all new alignments"
            })

        report_text = "Created {} alignments from the given alignment set.".format(len(alignments))
        report_info = report_client.create({
            "workspace_name": params["workspace_name"],
            "report": {
                "objects_created": created_objects,
                "text_message": report_text
            }
        })
        return {'report_name': report_info['name'], 'report_ref': report_info['ref']}


    def run_star(self, input_params):
        """
        run_star: run the STAR app
        """
        log('--->\nrunning STARUtil.run_star\n' +
            'params:\n{}'.format(json.dumps(input_params, indent=1)))

	#0. preprocessing the input parameters
        input_params = self._process_params(input_params)
	params = {
            'workspace_name': input_params[self.PARAM_IN_WS],
            'runMode': 'genomeGenerate',
            'runThreadN': input_params[self.PARAM_IN_THREADN],
            'outFileNamePrefix': input_params[self.PARAM_IN_OUTFILE_PREFIX]
	}

	#1. Converting refs to file locations in the scratch area
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
        readsNames = list()
	for source_reads in reads_refs:
            try:
                print("Fetching FASTA file from reads reference {}".format(source_reads['ref']))
                ret_reads = fetch_reads_from_reference(source_reads['ref'], self.callback_url)
                if ret_reads.get("file_fwd", None) is not None:
                    print("Done fetching FASTA file with name = {}".format(ret_reads.get("file_fwd", None)))
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
                if ret_reads.get("file_fwd", None) is not None:
                    readsFiles.append(ret_reads["file_fwd"])
                    readsNames.append(os.path.basename(ret_reads["file_fwd"]))
                    if ret_reads.get("file_rev", None) is not None:
                        readsFiles.append(ret_reads["file_rev"])

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

        sjdbGTFfile_ref = input_params.get("sjdbGTFfile_ref", None)
	if sjdbGTFfile_ref is not None:
            try:
		print("Fetching FASTA file from object {}".format(sjdbGTFfile_ref))
		sjdbGTFfile = fetch_fasta_from_object(sjdbGTFfile_ref, self.workspace_url, self.callback_url)
		print("Done fetching FASTA file! Path = {}".format(sjdbGTFfile.get("path", None)))
            except ValueError:
		print("Incorrect object type for fetching a FASTA file!")
		raise

            if sjdbGTFfile.get("path", None) is None:
		raise RuntimeError("FASTA file fetched from object {} doesn't seem to exist!".format(sjdbGTFfile_ref))
            else:
		params['sjdbGTFfile'] = sjdbGTFfile["path"]
                if input_params.get('sjdbOverhang', None) is not None :
                    params['sjdbOverhang'] = input_params['sjdbOverhang']

	#2. Running star
	star_out = self._exec_star(params, readsNames[0])

	#3. Uploading the alignment and generating report
        print("Uploading STAR output object and report...")
        alignment_file = "{}Aligned.out.sam".format(input_params[self.PARAM_IN_OUTFILE_PREFIX])
	alignment_file = os.path.join(star_out, alignment_file)

        # Upload the alignment with ONLY the first reads_ref for now
        input_params['alignment_name'] = "{}_Aligned".format(readsNames[0])
	alignment_ref = self.upload_STARalignment(input_params, readsInfo[0], alignment_file)

        reportVal = self._generate_report(alignment_ref, star_out, input_params)

        returnVal = {
            'output_folder': star_out,
            'alignment_ref': alignment_ref
        }

        returnVal.update(reportVal)

        return returnVal

