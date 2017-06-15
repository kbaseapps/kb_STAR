
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

from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils

from file_util import (
    fetch_reads_from_reference,
    fetch_reads_refs_from_sampleset,
    fetch_fasta_from_object
)

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class STARUtil:
    STAR_VERSION = 'STAR 2.5.3a'
    STAR_BIN = '/kb/deployment/bin/STAR'
    STAR_IDX = '/kb/module/STAR_genome_dir/'
    #STAR_DATA = '/kb/module/testReads'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_FASTA_FILES = 'genomeFastaFiles'
    PARAM_IN_OUTFILE_PREFIX = 'outFileNamePrefix'
    PARAM_IN_STARMODE = 'runMode'
    PARAM_IN_THREADN = 'runThreadN'
    PARAM_IN_READS_FILES = 'readFilesIn'
 
    INVALID_WS_OBJ_NAME_RE = re.compile('[^\\w\\|._-]')
    INVALID_WS_NAME_RE = re.compile('[^\\w:._-]')

    PARAM_IN_LIB = 'reads_ref'
    PARAM_IN_GENOME = 'genome_ref'
    PARAM_IN_ASSEMBLY = 'assembly_ref'

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
		if params.get(self.PARAM_IN_FASTA_FILES, None) is None:
		    raise ValueError(self.PARAM_IN_FASTA_FILES + 
				' parameter is required for generating genome index')
		if type(params[self.PARAM_IN_FASTA_FILES]) != list:
		    raise ValueError(self.PARAM_IN_FASTA_FILES + ' must be a list')
		if params.get(self.PARAM_IN_FASTA_FILES, None) is None:
		    raise ValueError('At least one FASTA file must be provided for generating genome index')

        if (params.get(self.PARAM_IN_STARMODE, None) is not None and 
		params[self.PARAM_IN_STARMODE] != "genomeGenerate"):
            if params.get(self.PARAM_IN_READS_FILES, None) is None:
		raise ValueError(self.PARAM_IN_READS_FILES + 
				' parameter is required for genome mapping')
            if type(params[self.PARAM_IN_READS_FILES]) != list:
		raise ValueError(self.PARAM_IN_READS_FILES + ' must be a list')
            if params.get(self.PARAM_IN_READS_FILES, None) is None:
		raise ValueError('At least one reads file must be provided for genome mapping')

        if params.get(self.PARAM_IN_THREADN, None) is not None:
            if not isinstance(params[self.PARAM_IN_THREADN], int):
                raise ValueError(self.PARAM_IN_HASH_THREADN + ' must be of type int')
	else:
	    params[self.PARAM_IN_THREADN] = 1

        if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is None:
	    params[self.PARAM_IN_OUTFILE_PREFIX] = 'STARoutput_'

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
	idx_cmd.append(params[self.PARAM_IN_THREADN])

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
        out_folder = params['out_folder']
        wsname = params[self.PARAM_IN_WS]

        # STEP 2: construct the command for running STAR mapping
        mp_cmd = [self.STAR_BIN]
	mp_cmd.append('--genomeDir')
	mp_cmd.append(out_folder)
	mp_cmd.append('--' + self.PARAM_IN_STARMODE)
	mp_cmd.append(params[self.PARAM_IN_STARMODE])
	mp_cmd.append('--' + self.PARAM_IN_THREADN)
	mp_cmd.append(params[self.PARAM_IN_THREADN])
        
	if params.get(self.PARAM_IN_READS_FILES, None) is not None:
            print('Input reads files:' + pformat(params[self.PARAM_IN_READS_FILES]))
	    mp_cmd.append('--' + self.PARAM_IN_READS_FILES)
            for reads_file in params[self.PARAM_IN_READS_FILES]:
	        mp_cmd.append(reads_file)
		readName, readsExtension = os.path.splitext(reads_file)
		if readsExtension == '.gz':
			mp_cmd.append('--readFilesCommand')
			mp_cmd.append('gunzip')
			mp_cmd.append('-c')

	if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
	    mp_cmd.append('--' + self.PARAM_IN_OUTFILE_PREFIX)
	    mp_cmd.append(os.path.join(out_folder, params[PARAM_IN_OUTFILE_PREFIX]))

        # appending the advanced optional inputs--TODO

        # STEP 3 return mp_cmd
        print ('STAR mapping CMD:')
        print ' '.join(mp_cmd)
        return mp_cmd

    def _exec_indexing(self, params):
        self.log('Running STAR index generating with params:\n' + pformat(params))
        STAR_cmd = self._construct_indexing_cmd(params)

        p = subprocess.Popen(STAR_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()

        self.log('Return code: ' + str(retcode))
        if p.returncode != 0:
            raise ValueError('Error running STAR index generating, return code: ' + str(retcode) + '\n')

        return p.returncode

    def _exec_mapping(self, params):
        self.log('Running STAR mapping with params:\n' + pformat(params))
        STAR_cmd = self._construct_mapping_cmd(params)
        p = subprocess.Popen(STAR_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()
        self.log('Return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running STAR mapping, return code: ' + str(p.returncode) + '\n')

        return p.returncode

    def __init__(self, config):
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
	self.workspace_url = config['workspace-url']
        self.working_dir = self.scratch

    def _exec_star(self, params):
        params = self._process_params(params)

        outdir = os.path.join(self.scratch, 'star_output_dir')
        self._mkdir_p(outdir)
        tmpdir = os.path.join(self.scratch, 'star_tmp_dir')
	self._mkdir_p(tmpdir)

        # build the parameters
        params_idx = {
                'workspace_name': params[self.PARAM_IN_WS],
                'runMode': params[self.PARAM_IN_STARMODE], #'genomeGenerate',
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
                'readsFilesIn': params[self.PARAM_IN_READS_FILES],
		'outFileNamePrefix': params[self.PARAM_IN_OUTFILE_PREFIX], 
                'out_folder': outdir
        }

        ret = 1
        try:
	    if(params[self.PARAM_IN_STARMODE]=='genomeGenerate'):
            	ret = self._exec_indexing(params_idx)
	    else:
		ret = 0
	    while( ret != 0 ):
                time.sleep(1)
        except ValueError as eindx:
            self.log('STAR genome indexing raised error:\n')
            print(eindx)
        else:#no exception raised by genome indexing and STAR returns 0, then run mapping
            ret = 1 
            try:
                ret = self._exec_mapping(params_mp)
                while( ret != 0 ):
                    time.sleep(1)
            except ValueError as eg:
                self.log('STAR mapping raised error:\n')
                print(eg)
            else:#no exception raised by STAR mapping and STAR returns 0, then move to saving and reporting  
                ret = outdir
        return ret


    def get_star_index(self, source_ref):
        """
        Builds or fetches the index file(s) as necessary, unpacks them in a directory.
        Returns a string representing the path and prefix of the index files.
        E.g. if there are a set of files like "foo.1.ht2", "foo.2.ht2", etc. all in the
        "my_reads" directory, this will return "my_reads/foo"
        """
        idx_prefix = self._fetch_star_index(source_ref, {})
        if idx_prefix:
            return idx_prefix
        else:
            return self._build_star_index(source_ref, {})

    def _fetch_star_index(self, source_ref, options):
        """
        Fetches STAR indexes from a remote location, if they're available.
        """
        # TODO: insert fetch code from file cache HERE.
        return None

    def _build_star_index(self, source_ref, options):
        """
        Runs star-build to build the index files and directory for use in STAR.
        It creates a directory called "kb_star_idx" and places the index files within.
        This also caches them, um, elsewhere so they can be found again.
        """
        # check options and raise ValueError here as needed.
        if source_ref is None:
            raise ValueError("Missing reference object needed to build a STAR index.")
        print("Building STAR index files for {}".format(source_ref))
        idx_dir = "kb_star_idx"
        idx_prefix = "kb_star_idx-" + str(uuid.uuid4())
        try:
            slef._mkdir_p(os.path.join(self.working_dir, idx_dir))
        except OSError:
            print("Ignoring error for already existing {} directory".format(idx_dir))
        try:
            print("Fetching FASTA file from object {}".format(source_ref))
            fasta_file = fetch_fasta_from_object(source_ref, self.workspace_url, self.callback_url)
            print("Done fetching FASTA file! Path = {}".format(fasta_file.get("path", None)))
        except ValueError:
            print("Incorrect object type for fetching a FASTA file!")
            raise

        fasta_path = fasta_file.get("path", None)
        if fasta_path is None:
            raise RuntimeError("FASTA file fetched from object {} doesn't seem to "
                               "exist!".format(source_ref))
        build_star_cmd = [
            "star-build",
            "-f",
            fasta_path
        ]
        if options.get('num_threads', None) is not None:
            build_star_cmd.extend(["-p", options['num_threads']])
        idx_prefix_path = os.path.join(self.working_dir, idx_dir, idx_prefix)
        build_star_cmd.append(idx_prefix_path)
        print("Executing build-star command: {}".format(build_star_cmd))
        p = subprocess.Popen(build_star_cmd, shell=False)
        ret_code = p.wait()
        if ret_code != 0:
            raise RuntimeError('Failed to generate STAR index files!')
        print("Done! STAR index files created with prefix {}".format(idx_prefix_path))
        return idx_prefix_path

    def upload_alignment(self, input_params, reads_info, alignment_file):
        """
        Uploads the alignment file + metadata.
        Returns the STAR alignment reference.
        """
        aligner_opts = dict()
        for k in input_params:
            aligner_opts[k] = str(input_params[k])

        align_upload_params = {
            "destination_ref": "{}/{}".format(input_params["workspace_name"], input_params["alignmentset_name"]),
            "file_path": alignment_file,
            "library_type": reads_info["style"],  # single or paired end,
            "condition": reads_info["condition"],
            "assembly_or_genome_ref": input_params["genome_ref"],
            "read_library_ref": reads_info["object_ref"],
            "aligned_using": "STAR",
            "aligner_version":self.STAR_VERSION,
            "aligner_opts": aligner_opts
        }

        if "sampleset_ref" in reads_info:
            align_upload_params["sampleset_ref"] = reads_info["sampleset_ref"]
        print("Uploading completed alignment")
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

	#1. Converting refs to file locations of the scratch area
        #obj_ref reads_ref;
	try:
	    print("Fetching FASTA file from object {}".format(reads_ref))
	    reads_fasta_file = fetch_fasta_from_object(reads_ref, self.workspace_url, self.callback_url)
	    print("Done fetching FASTA file! Path = {}".format(reads_fasta_file.get("path", None)))
	except ValueError:
	    print("Incorrect object type for fetching a FASTA file!")
	    raise

	if reads_fasta_file.get("path", None) is None:
	    raise RuntimeError("FASTA file fetched from object {} doesn't seem to "
		       "exist!".format(reads_ref))

        #obj_ref genome_ref;
	try:
	    print("Fetching FASTA file from object {}".format(genome_ref))
	    genome_fasta_file = fetch_fasta_from_object(genome_ref, self.workspace_url, self.callback_url)
	    print("Done fetching FASTA file! Path = {}".format(genome_fasta_file.get("path", None)))
	except ValueError:
	    print("Incorrect object type for fetching a FASTA file!")
	    raise

	if genome_fasta_file.get("path", None) is None:
	    raise RuntimeError("FASTA file fetched from object {} doesn't seem to "
		       "exist!".format(genome_ref))

	#obj_ref sjdbGTFfile_ref
	try:
	    print("Fetching FASTA file from object {}".format(sjdbGTFfile_ref))
	    sjdbGTFfile = fetch_fasta_from_object(sjdbGTFfile_ref, self.workspace_url, self.callback_url)
	    print("Done fetching FASTA file! Path = {}".format(sjdbGTFfile.get("path", None)))
	except ValueError:
	    print("Incorrect object type for fetching a FASTA file!")
	    raise

	if sjdbGTFfile.get("path", None) is None:
	    raise RuntimeError("FASTA file fetched from object {} doesn't seem to "
		       "exist!".format(sjdbGTFfile_ref))

	#list<obj_ref> genomeFastaFiles
	genomeFastaFiles = list()
	for source_ref in input_params['genomeFastaFiles']:
		try:
		    print("Fetching FASTA file from object {}".format(source_ref))
		    fasta_file = fetch_fasta_from_object(source_ref, self.workspace_url, self.callback_url)
		    print("Done fetching FASTA file! Path = {}".format(fasta_file.get("path", None)))
		except ValueError:
		    print("Incorrect object type for fetching a FASTA file!")
		    raise

		if fasta_file.get("path", None) is None:
		    raise RuntimeError("FASTA file fetched from object {} doesn't seem to "
			       "exist!".format(source_ref))

		genomeFastaFiles.append(fasta_file)

	#list<obj_ref> readsFilesIn
	readsFiles = list()
	for source_ref in input_params['readsFilesIn']:
		try:
		    print("Fetching FASTA file from object {}".format(source_ref))
		    reads_file = fetch_fasta_from_object(source_ref, self.workspace_url, self.callback_url)
		    print("Done fetching FASTA file! Path = {}".format(reads_file.get("path", None)))
		except ValueError:
		    print("Incorrect object type for fetching a FASTA file!")
		    raise

		if reads_file.get("path", None) is None:
		    raise RuntimeError("FASTA file fetched from object {} doesn't seem to "
			       "exist!".format(source_ref))

		readsFiles.append(reads_file)

	#2. Running star
	params = {
            'reads_fasta_file': reads_fasta_file,
            'genome_fasta_file': genome_fasta_file,
	    'runMode': input_params['runMode'],
	    'runThreadN': input_params['runThreadN'],
	    'genomeFastaFiles': genomeFastaFiles,
	    'sjdbGTFfile': sjdbGTFfile,
	    'sjdbOverhang': input_params['sjdbOverhang'],
	    'readFilesIn': readsFiles,
            'outFileNamePrefix': input_params['outFileNamePrefix']
	}
	star_out = self._exec_star(params)

        log('STAR result files have been saved to: {}'.format(star_out))
        log('STAR has generated files:\n{}'.format('\n'.join(os.listdir(star_out))))
	
	#alignment_file = os.path.join(self.working_dir, "{}.sam".format(output_file))
        print("Uploading STAR output object and report...")
        #alignment_ref = self.upload_alignment(input_params, reads, alignment_file)

        #reportVal = self._generate_report(alignment_ref, output_dir, input_params)

        returnVal = {
            'output_folder': star_out,
            'alignment_ref': alignment_ref
        }

        #returnVal.update(reportVal)

        return returnVal

