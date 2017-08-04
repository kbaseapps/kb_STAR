import json
import os
import re
import time
import copy
import sys
from pprint import pprint

from KBParallel.KBParallelClient import KBParallel
from STAR.Utils.STARUtils import STARUtils

from file_util import (
    extract_geneCount_matrix
)

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

class STAR_Aligner(object):

    def __init__(self, config, provenance, index_dir, output_dir):
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.parallel_runner = KBParallel(self.callback_url)
        self.star_utils = STARUtils(self.config, provenance)
        self.provenance = provenance
        self.star_idx_dir = index_dir
        self.star_out_dir = output_dir


    def run_star_indexing(self, input_params):
        """
        Runs STAR in genomeGenerate mode to build the index files and directory for STAR mapping.
        It creates a directory as defined by self.star_idx_dir in the scratch area that houses
        the index files.
        """
        genome_params = copy.deepcopy(input_params)

        # build the indexing parameters
        params_idx = self.star_utils._get_indexing_params(genome_params, self.star_idx_dir)

        ret = 1
        try:
            if genome_params[STARUtils.PARAM_IN_STARMODE]=='genomeGenerate':
                ret = self.star_utils._exec_indexing(params_idx)
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


    def run_star_mapping(self, params, rds_files, rds_name):
        """
        Runs STAR in alignReads mode for STAR mapping.
        It creates a directory as defined by self.star_out_dir with a subfolder named after the reads
        """
        params_mp = self.star_utils._get_mapping_params(
                        params, rds_files, rds_name, self.star_idx_dir, self.star_out_dir
                )

        retVal = {}
        params_mp[STARUtils.PARAM_IN_STARMODE] = 'alignReads'
        try:
            ret = self.star_utils._exec_mapping(params_mp)
            while( ret != 0 ):
                time.sleep(1)
        except ValueError as emp:
            log('STAR mapping raised error:\n')
            pprint(emp)
            retVal = {'star_idx': self.star_dx_dir, 'star_output': None}
        else:#no exception raised by STAR mapping and STAR returns 0, then move to saving and reporting  
            retVal = {'star_idx': self.star_idx_dir, 'star_output': params_mp.get('align_output')}

        return retVal


    def star_run_single(self, validated_params):
        """
        Performs a single run of STAR against a single reads reference. The rest of the info
        is taken from the params dict - see the spec for details.
        """
        log('--->\nrunning STARUtils.run_single\n' +
                'params:\n{}'.format(json.dumps(validated_params, indent=1)))

	# convert the input parameters (from refs to file paths, especially)
        params_ret = self.star_utils.convert_params(validated_params)
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
            bam_sort = ''
            prefix = ''
            if input_params.get('outSAMtype', None) == 'BAM':
                bam_sort = 'sortedByCoord'
            if input_params.get(STARUtils.PARAM_IN_OUTFILE_PREFIX, None) is not None:
                prefix = input_params[STARUtils.PARAM_IN_OUTFILE_PREFIX]
                output_bam_file = '{}Aligned.{}.out.bam'.format(prefix, bam_sort)
            else:
                output_bam_file = 'Aligned.{}.out.bam'.format(bam_sort)
            output_bam_file = os.path.join(star_mp_ret['star_output'], output_bam_file)

            #print("Uploading STAR output object...")
            # Upload the alignment
            upload_results = self.star_utils.upload_STARalignment(input_params, reads, output_bam_file)
            alignment_ref = upload_results['obj_ref']
            alignment_obj = {
                'ref': alignment_ref,
                'name': reads_ref['alignment_output_name']
            }
            alignment_objs.append({
                'reads_ref': reads_ref['ref'],
                'AlignmentObj': alignment_obj
            })

            singlerun_output_info['output_dir'] = star_mp_ret['star_output']
            singlerun_output_info['output_bam_file'] = output_bam_file

            singlerun_output_info['upload_results'] = upload_results

            #workspace_name = validated_params[STARUtils.PARAM_IN_WS]
            #expr_suffix = validated_params['expression_suffix']
            #gtf_file = validated_params['sjdbGTFfile']
            #expression_ref = self.star_utils._save_expression(
            #                    star_mp_ret['star_output'],
            #                    alignment_ref,
            #                    workspace_name,
            #                    gtf_file,
            #                    expr_suffix)

            if input_params.get("create_report", 0) == 1:
                report_info = self.star_utils.generate_report_for_single_run(singlerun_output_info, input_params)

        if ret_fwd is not None:
            os.remove(ret_fwd)
            if reads_info.get('file_rev', None) is not None:
                os.remove(reads_info["file_rev"])

        return {'alignmentset_ref': None,
                'output_directory': singlerun_output_info['output_dir'],
                'output_info': singlerun_output_info,
                'alignment_objs': alignment_objs,
                'report_name': report_info['name'],
                'report_ref': report_info['ref']
        }


    def star_run_batch(self, validated_params):
        ''' convert the input parameters (from refs to file paths, especially)'''
        reads_refs = self.star_utils.get_reads_refs(validated_params)

        # build task list and send it to KBParallel
        tasks = []
        for r in reads_refs:
            tasks.append(
                    self.star_utils.build_single_execution_task(
                        r['ref'], validated_params)
                    )

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

        batch_result = self.star_utils.process_batch_result(results, validated_params, reads_refs)
        batch_result['output_directory'] = self.star_out_dir

        return batch_result

