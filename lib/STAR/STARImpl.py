# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
from Bio import SeqIO
from pprint import pprint, pformat
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport

import re
import shutil
import subprocess
import numpy as np
from datetime import datetime
from pprint import pformat, pprint
import time
import uuid

from KBaseReport.baseclient import ServerError as _RepError
from kb_quast.kb_quastClient import kb_quast
from kb_quast.baseclient import ServerError as QUASTError
from ReadsUtils.ReadsUtilsClient import ReadsUtils 
from ReadsUtils.baseclient import ServerError
from Workspace.WorkspaceClient import Workspace as workspaceService
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
    GIT_COMMIT_HASH = "828fdeaf6369345768c448f3064522cbe55aeb95"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    STAR_BIN = '/kb/deployment/bin/STAR'
    STAR_DATA = '/kb/module/work/tmp'
    #STAR_DATA = '/kb/module/testReads'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_FASTA_FILES = 'genomeFastaFiles'
    PARAM_IN_OUTFILE_NAME = 'outFileNamePrefix'
    PARAM_IN_STARMODE = 'runMode'
    PARAM_IN_THREADN = 'runThreadN'
    PARAM_IN_READS_FILES = 'readFilesIn'
 
    INVALID_WS_OBJ_NAME_RE = re.compile('[^\\w\\|._-]')
    INVALID_WS_NAME_RE = re.compile('[^\\w:._-]')

    PARAM_IN_LIB = 'reads_ref'
    PARAM_IN_GENOME = 'genome_ref'
    PARAM_IN_ASSEMBLY = 'assembly_ref'

    def log(self, message, prefix_newline=False):
            print(('\n' if prefix_newline else '') +
                str(time.time()) + ': ' + str(message))

    def process_params(self, params):
        if (self.PARAM_IN_WS not in params or
                not params[self.PARAM_IN_WS]):
            raise ValueError(self.PARAM_IN_WS + ' parameter is required')

        if self.PARAM_IN_STARMODE not in params or params[self.PARAM_IN_STARMODE] is None:
            params[self.PARAM_IN_STARMODE] = 'alignReads'
        if self.PARAM_IN_STARMODE in params and params[self.PARAM_IN_STARMODE] == "genomeGenerate":
            if self.PARAM_IN_FASTA_FILES not in params:
		raise ValueError(self.PARAM_IN_FASTA_FILES + ' parameter is required for generating genome index')
            if type(params[self.PARAM_IN_FASTA_FILES]) != list:
		raise ValueError(self.PARAM_IN_FASTA_FILES + ' must be a list')
            if not params[self.PARAM_IN_FASTA_FILES]:
		raise ValueError('At least one FASTA file must be provided for generating genome index')

        if self.PARAM_IN_STARMODE in params and not params[self.PARAM_IN_STARMODE] == "genomeGenerate":
            if self.PARAM_IN_READS_FILES not in params:
		raise ValueError(self.PARAM_IN_READS_FILES + ' parameter is required for genome mapping')
            if type(params[self.PARAM_IN_READS_FILES]) != list:
		raise ValueError(self.PARAM_IN_READS_FILES + ' must be a list')
            if not params[self.PARAM_IN_READS_FILES]:
		raise ValueError('At least one reads file must be provided for genome mapping')

        if self.PARAM_IN_THREADN in params:
            if not isinstance(params[self.PARAM_IN_THREADN], int):
                raise ValueError(self.PARAM_IN_HASH_THREADN + ' must be of type int')
	else:
	    params[self.PARAM_IN_THREADN] = 1

        if (self.PARAM_IN_OUTFILE_NAME in params and
                params[self.PARAM_IN_OUTFILE_NAME]):
            if self.INVALID_WS_OBJ_NAME_RE.search(params[self.PARAM_IN_OUTFILE_NAME]):
            	raise ValueError('Invalid workspace object name ' +
                             params[self.PARAM_IN_OUTFILE_NAME])

    def construct_indexing_cmd(self, params):
        # STEP 1: get the working folder housing the STAR results as well as the reads info
        wsname = params[self.PARAM_IN_WS]
        out_folder = params['out_folder']

        # STEP 2: construct the command for running `STAR --runMode genomeGenerate...`
        indx_cmd = [self.STAR_BIN]
	indx_cmd.append('--genomeDir')
	indx_cmd.append(out_folder)
	indx_cmd.append('--' + self.PARAM_IN_STARMODE)
	indx_cmd.append('genomeGenerate')
	indx_cmd.append('--' + self.PARAM_IN_THREADN)
	indx_cmd.append(params[self.PARAM_IN_THREADN])

	if self.PARAM_IN_FASTA_FILES in params:
            print('Input fasta reads files:' + pformat(params[self.PARAM_IN_FASTA_FILES]))
	    indx_cmd.append('--' + self.PARAM_IN_FASTA_FILES)
            for fastaf in params[self.PARAM_IN_FASTA_FILES]:
	        indx_cmd.append(fastaf)

	# appending the standard optional inputs
        if 'sjdbGTFfile' in params and not (params['sjdbGTFfile'] is None):
            indx_cmd.append('--sjdbGTFfile')
            indx_cmd.append(params['sjdbGTFfile'])
        if 'sjdbOverhang' in params and not (params['sjdbOverhang'] is None) and params['sjdbOverhang'] > 0:
            indx_cmd.append('--sjdbOverhang')
            indx_cmd.append(str(params['sjdbOverhang']))
        # appending the advanced optional inputs--TODO
		
        # STEP 3: return indx_cmd
        print ('STAR indexing CMD:')
        print ' '.join(indx_cmd)
        return indx_cmd

    def construct_mapping_cmd(self, params):
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
        
	if self.PARAM_IN_READS_FILES in params:
            print('Input reads files:' + pformat(params[self.PARAM_IN_READS_FILES]))
	    mp_cmd.append('--' + self.PARAM_IN_READS_FILES)
            for readsf in params[self.PARAM_IN_READS_FILES]:
	        mp_cmd.append(readsf)

	mp_cmd.append('--' + self.PARAM_IN_OUTFILE_NAME)
	mp_cmd.append(params[PARAM_IN_OUTFILE_NAME])
        # appending the advanced optional inputs--TODO

        # STEP 3 return vg_cmd
        print ('STAR mapping CMD:')
        print ' '.join(mp_cmd)
        return mp_cmd

    def exec_indexing(self, params):
        self.log('Running STAR index generating with params:\n' + pformat(params))
        STAR_cmd = self.construct_indexing_cmd(params)

        p = subprocess.Popen(STAR_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()

        self.log('Return code: ' + str(retcode))
        if p.returncode != 0:
            raise ValueError('Error running STAR index generating, return code: ' + str(retcode) + '\n')

        return p.returncode

    def exec_mapping(self, params):
        self.log('Running STAR mapping with params:\n' + pformat(params))
        STAR_cmd = self.construct_mapping_cmd(params)
        p = subprocess.Popen(STAR_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()
        self.log('Return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running STAR mapping, return code: ' + str(p.returncode) + '\n')

        return p.returncode


    def exec_star(self, params):
        outdir = os.path.join(self.scratch, 'star_output_dir')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        tmpdir = os.path.join(self.scratch, 'star_tmp_dir')
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)

        # build the parameters
        params_indx = {
                'workspace_name': params[self.PARAM_IN_WS],
                'runMode': 'genomeGenerate',
		'runThreadN': params[self.PARAM_IN_THREADN],
                'genomeFastaFiles': params[self.PARAM_IN_FASTA_FILES],
                'out_folder': outdir
        }
        params_mp = {
                'workspace_name': params[self.PARAM_IN_WS],
                'runMode': params[self.PARAM_IN_STARMODE],
		'runThreadN': params[self.PARAM_IN_THREADN],
                'readsFilesIn': params[self.PARAM_IN_READS_FILES],
		'outFileNamePrefix': params[self.PARAM_IN_OUTFILE_NAME], 
                'out_folder': outdir
        }

        if 'sjdbGTFfile' in params and not (params['sjdbGTFfile'] is None):
            params_indx['sjdbGTFfile'] = params['sjdbGTFfile']
        if 'sjdbOverhang' in params and not (params['sjdbOverhang'] is None) :
            params_indx['sjdbOverhang'] = params['sjdbOverhang']

        ret = 1
        try:
            ret = self.exec_indexing(params_indx)
            while( ret != 0 ):
                time.sleep(1)
        except ValueError as eindx:
            self.log('STAR genome indexing raised error:\n')
            print(eindx)
        else:#no exception raised by genome indexing and STAR returns 0, then run mapping
            ret = 1 
            try:
                ret = self.exec_mapping(params_mp)
                while( ret != 0 ):
                    time.sleep(1)
            except ValueError as eg:
                self.log('STAR mapping raised error:\n')
                print(eg)
            else:#no exception raised by STAR mapping and STAR returns 0, then move to saving and reporting  
                ret = outdir
        return ret

    # adapted from
    # https://github.com/kbaseapps/kb_SPAdes/blob/master/lib/kb_SPAdes/kb_SPAdesImpl.py
    # which was adapted from
    # https://github.com/kbase/transform/blob/master/plugins/scripts/convert/trns_transform_KBaseFile_AssemblyFile_to_KBaseGenomes_ContigSet.py
    def load_stats(self, input_file_name):
        self.log('Starting conversion of FASTA to KBaseGenomeAnnotations.Assembly')
        self.log('Building Object.')
        if not os.path.isfile(input_file_name):
            raise Exception('The input file name {0} is not a file!'.format(input_file_name))
        with open(input_file_name, 'r') as input_file_handle:
            contig_id = None
            sequence_len = 0
            fasta_dict = dict()
            first_header_found = False
            # Pattern for replacing white space
            pattern = re.compile(r'\s+')
            for current_line in input_file_handle:
                if (current_line[0] == '>'):
                    # found a header line
                    # Wrap up previous fasta sequence
                    if not first_header_found:
                        first_header_found = True
                    else:
                        fasta_dict[contig_id] = sequence_len
                        sequence_len = 0
                    fasta_header = current_line.replace('>', '').strip()
                    try:
                        contig_id = fasta_header.strip().split(' ', 1)[0]
                    except:
                        contig_id = fasta_header.strip()
                else:
                    sequence_len += len(re.sub(pattern, '', current_line))
        # wrap up last fasta sequence
        if not first_header_found:
            raise Exception("There are no contigs in this file")
        else:
            fasta_dict[contig_id] = sequence_len
        return fasta_dict

    def generate_report(self, input_file_name, params, out_folder, wsname):
        self.log('Generating and saving report')

        fasta_stats = self.load_stats(input_file_name)
        lengths = [fasta_stats[contig_id] for contig_id in fasta_stats]

        assembly_ref = params[self.PARAM_IN_WS] + '/' + params[self.PARAM_IN_CS_NAME]

        report = ''
        report += 'Velvet results saved to: ' + wsname + '/' + out_folder + '\n'
        report += 'Assembly saved to: ' + assembly_ref + '\n'
        report += 'Assembled into ' + str(len(lengths)) + ' contigs.\n'
        report += 'Avg Length: ' + str(sum(lengths) / float(len(lengths))) + ' bp.\n'

        # compute a simple contig length distribution
        bins = 10
        counts, edges = np.histogram(lengths, bins)  # @UndefinedVariable
        report += 'Contig Length Distribution (# of contigs -- min to max ' + 'basepairs):\n'
        for c in range(bins):
            report += '   ' + str(counts[c]) + '\t--\t' + str(edges[c]) + ' to ' + str(edges[c + 1]) + ' bp\n'
        print('Running QUAST')
        kbq = kb_quast(self.callbackURL)
        quastret = kbq.run_QUAST({'files': [{'path': input_file_name,
                                             'label': params[self.PARAM_IN_CS_NAME]}]})
        print('Saving report')
        kbr = KBaseReport(self.callbackURL)
        report_info = kbr.create_extended_report(
            {'message': report,
             'objects_created': [{'ref': assembly_ref, 'description': 'Assembled contigs'}],
             'direct_html_link_index': 0,
             'html_links': [{'shock_id': quastret['shock_id'],
                             'name': 'report.html',
                             'label': 'QUAST report'}
                            ],
             'report_object_name': 'kb_star_report_' + str(uuid.uuid4()),
             'workspace_name': params[self.PARAM_IN_WS]
            })
        reportName = report_info['name']
        reportRef = report_info['ref']
        return reportName, reportRef

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        
        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']

        #END_CONSTRUCTOR
        pass


    def run_star(self, ctx, params):
        """
        The actual function is declared using 'funcdef' to specify the name
        and input/return arguments to the function.  For all typical KBase
        Apps that run in the Narrative, your function should have the 
        'authentication required' modifier.
        :param params: instance of type "STARParams" (Arguments for
           star_generate_indexes string reads_ref, assembly_ref and
           genome_ref: KBase style variable references string runMode:
           default: alignReads type of the run: alignReads => map reads
           genomeGenerate => generate genome files inputAlignmentsFromBAM =>
           input alignments from BAM. Presently only works with -outWigType
           and -bamRemoveDuplicates. liftOver => lift-over of GTF files
           (-sjdbGTFfile) between genome assemblies using chain file(s) from
           -genomeChainFiles. int runThreadN: default: 1 number of threads to
           run STAR list<string> genomeFastaFiles: path(s) to the fasta files
           with genomic sequences for genome generation. Only used if
           runMode==genomeGenerate.These files should be plain text FASTA
           files, they *cannot* be zipped. list<string> readFilesIn: default:
           Read1 Read2 paths to files that contain input read1 (and, if
           needed, read2) string sjdbGTFfile: default: -; path to the file
           with annotated transcripts in the standard GTF format int
           sjdbOverhang: default: 100; int>0: length of the donor/acceptor
           sequence on each side of the junctions, ideally = (ReadLength - 1)
           string outFileNamePrefix: you can change the file prefixes using
           --outFileNamePrefix /path/to/output/dir/prefix. By default, this
           parameter is ./, i.e. all output files are written in the current
           directory @optional sjdbGTFfile @optional sjdbOverhang) ->
           structure: parameter "reads_ref" of String, parameter
           "assembly_ref" of String, parameter "genome_ref" of String,
           parameter "workspace_name" of String, parameter "runMode" of
           String, parameter "runThreadN" of Long, parameter
           "genomeFastaFiles" of list of String, parameter "sjdbGTFfile" of
           String, parameter "sjdbOverhang" of Long, parameter "readFilesIn"
           of list of String, parameter "outFileNamePrefix" of String
        :returns: instance of type "STARResults" (Here is the definition of
           the output of the function.  The output can be used by other SDK
           modules which call your code, or the output visualizations in the
           Narrative.  'report_name' and 'report_ref' are special output
           fields- if defined, the Narrative can automatically render your
           Report.) -> structure: parameter "reads_alignment_ref" of String,
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_star
	# Step 1--generate genome indexes
        self.log('Running run_star with params:\n' + pformat(params))

        token = ctx['token']
        wsname = params[self.PARAM_IN_WS]
        self.process_params(params)
        star_out = self.exec_star(params)
        self.log('STAR final return: ' + str(star_out))

        # STEP 2: parse the output and save back to KBase, create report in the same time
        if isinstance(star_out, str) and star_out != '':
                output_align_reads = os.path.join(star_out, '*.Aligned.out.sam')
                # generate report from contigs.fa
                report_name, report_ref = self.generate_report(output_contigs, params, star_out, wsname)

                # STEP 3: contruct the output to send back
                output = {'report_name': report_name, 'report_ref': report_ref}
        else:
            output = {'report_name': 'STAR failed', 'report_ref': None}
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
