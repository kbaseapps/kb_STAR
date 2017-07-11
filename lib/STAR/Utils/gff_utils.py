import os
import logging
import uuid
from script_utils import runProgram, log
import errno
import glob
from datetime import datetime
from pprint import pprint
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from DataFileUtil.DataFileUtilClient import DataFileUtil
from DataFileUtil.baseclient import ServerError as DFUError
from Workspace.WorkspaceClient import Workspace as Workspace
import shutil


class GFFUtils:
    """
    A wrapper to the jhu.edu's gff utilities. See
    http://ccb.jhu.edu/software/stringtie/gff.shtml
    """

    def __init__(self, config, logger=None):
        self.config = config
        self.logger = logger
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = os.path.join(config['scratch'], 'gff_utils_' + str(uuid.uuid4()))
        self.token = os.environ['KB_AUTH_TOKEN']
        self.ws_url = config['workspace-url']
        self.ws_client = Workspace(self.ws_url)
        self.dfu = DataFileUtil(self.callback_url)
        self._mkdir_p(self.scratch)
        pass

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

    def convert_gff_to_gtf(self, gff_file_path, gtf_file_path):
        """
        RConverts the specified GFF3 file to a GTF2 format
        :param gff_file_path: File path to a reference genome in GFF3 format
        :param gtf_file_path: path to GTF2 output file
        :returns 0 if successful, else 1
        """
        log('Running gffread...', level=logging.INFO, logger=self.logger)

        # check if input file exists
        if not os.path.exists(gff_file_path):
            raise RuntimeError(None, 'Input gff file does not exist: ' + str(gff_file_path))

        try:
            gtf_args = " -E {0} -T -o {1}".format(gff_file_path, gtf_file_path)
            log("Executing: gffread {0}".format(gtf_args), level=logging.INFO, logger=self.logger)
            rval = runProgram(logger=self.logger, progName="gffread", argStr=gtf_args)

            if rval['stderr'] != '':
                return 1

        except Exception as ex:
            log("Error executing gffread {0}. {1}".format(gtf_args, ex.message),
                logging.ERROR, self.logger)
            return 1

        return 0

    def convert_genome_to_gff(self, genome_ref, target_dir):
        """
        Converts the specified kbase genome object to gff file format
        :param genome_ref: workspace reference to kbase genome object
        :param target_dir: directory to which to write the output gff file.
        :return: path to the gff file. The gff file has the same name as the specified genome object
        """
        gfu = GenomeFileUtil(self.callback_url)
        gff_file_path = gfu.genome_to_gff({'genome_ref': genome_ref,
                                           'target_dir': target_dir})['file_path']

        return gff_file_path

    def convert_genome_to_gtf(self, genome_ref, gtf_file_path):
        """
        Converts the specified kbase genome object to gtf file format
        :param genome_ref: workspace reference to kbase genome object
        :param gtf_file_path:  path to the output gtf file
        :return: 0 if successful, else 1.
        """
        gff_file_path = self.convert_genome_to_gff(genome_ref, self.scratch)
        return self.convert_gff_to_gtf(gff_file_path, gtf_file_path)


    def get_gtf_file(self, genome_ref):
        pprint('Fetching genome GTF file path for {}'.format(genome_ref))
        g_info = self.ws_client.get_object_info3({'objects':[{'ref':genome_ref}]})
        g_obj_info = g_info.get("infos", [[]])[0]
        if len(g_obj_info) == 0:
                raise RuntimeError("An error occurred while fetching type info from the Workspace. "
                                           "No information returned for reference {}".format(genome_ref))
        genome_name = g_obj_info[1]
        ws_id = g_obj_info[6]

        ws_gtf = genome_name + "_GTF_Annotation"

        self.logger.info('GTF file from genome_ref: ' + ws_gtf)

        gtf_ref = str(ws_id) + '/' + ws_gtf
        a_info = self.ws_client.get_object_info3({'objects': [{'ref': gtf_ref}]})
        a_obj_info = a_info.get("infos", [[]])[0]
        a_obj_ref = str(a_obj_info[6]) + '/' + str(a_obj_info[0])

        try:
            ret_obj = self.dfu.get_objects({'object_refs': [a_obj_ref]})['data']
        except DFUError as e:
            self.log('Logging stacktrace from workspace exception:\n' + e.data)
            raise

        self.logger.debug('============ Genome gtf object ========')
        for line in str(pprint(ret_obj)).split('\n'):
            self.logger.debug(line)
        self.logger.debug('=======  END of Genome gtf object ========')

        # set the output dir
        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds() * 1000)
        output_dir = os.path.join(self.scratch, 'download_gtf_' + str(timestamp))
        os.mkdir(output_dir)

        file_ret = self.dfu.shock_to_file({
            'shock_id': ret_obj[0]['data']['handle']['id'],
            'file_path': output_dir
        })
        gtf_files = glob.glob(output_dir + '/*.gtf')

        if len(gtf_files) > 0:
            return gtf_files[0]
        else:
            raise ValueError('GTF file not found in object {}'.format(ws_gtf))

