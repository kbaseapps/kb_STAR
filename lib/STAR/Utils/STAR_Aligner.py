import json
import os
import re
import time
import copy
import sys
from pprint import pprint

from KBParallel.KBParallelClient import KBParallel
from STAR.Utils.STARUtils import STARUtils
from kb_QualiMap.kb_QualiMapClient import kb_QualiMap
from SetAPI.SetAPIServiceClient import SetAPI

from file_util import (
    extract_geneCount_matrix
)

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

class STAR_Aligner(object):

    def __init__(self, config, provenance):
        self.config = config
        self.workspace_url = config['workspace-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.srv_wiz_url = config['srv-wiz-url']
        self.parallel_runner = KBParallel(self.callback_url)
        self.provenance = provenance
        self.star_utils = STARUtils(self.scratch,
                        self.workspace_url,
                        self.callback_url,
                        self.srv_wiz_url, provenance)
        self.set_api_client = SetAPI(self.srv_wiz_url, service_ver='dev')
        self.qualimap = kb_QualiMap(self.callback_url, service_ver='dev')
        self.star_idx_dir = None
        self.star_out_dir = None

        # from the provenance, extract out the version to run by exact hash if possible
        self.my_version = 'release'
        if len(provenance) > 0:
            if 'subactions' in provenance[0]:
                self.my_version = self.get_version_from_subactions('kb_STAR', provenance[0]['subactions'])
        print('Running kb_STAR version = ' + self.my_version)


    def run_align(self, params):
        # 0. create the star folders
        if self.star_idx_dir is None:
            (idx_dir, out_dir) = self.star_utils.create_star_dirs(self.scratch)
            self.star_idx_dir = idx_dir
            self.star_out_dir = out_dir

        # 1. validate & process the input parameters
        validated_params = self.star_utils.process_params(params)
        input_obj_info = self.star_utils.determine_input_info(validated_params)

	# 2. convert the inpddut parameters (from refs to file paths, especially)
        input_params = self.star_utils.convert_params(validated_params)

        returnVal = {
            "report_ref": None,
            "report_name": None
        }
        if input_obj_info['run_mode'] == 'single_library':
            returnVal = self.star_run_single(input_params)

        if input_obj_info['run_mode'] == 'sample_set':
            #returnVal = self.star_run_batch_parallel(input_params)
            returnVal = self.star_run_batch_sequential(input_params)

        return returnVal


    def star_run_single(self, input_params):
        """
        Performs a single run of STAR against a single reads reference. The rest of the info
        is taken from the params dict - see the spec for details.
        """
        log('--->\nrunning STAR_Aligner.star_run_single\n' +
                'params:\n{}'.format(json.dumps(input_params, indent=1)))

        # 0. get index
        self.get_index(input_params)

        # 1. Prepare for mapping
        rds = None
        reads_refs = input_params[STARUtils.SET_READS]
        for r in reads_refs:
            if r['ref'] == input_params[STARUtils.PARAM_IN_READS]:
                rds = r
                break

        reads_info = self.star_utils._get_reads_info(rds, input_params[STARUtils.PARAM_IN_READS])

        rds_name = rds['alignment_output_name'].replace(input_params['alignment_suffix'], '')

        alignment_objs = list()
        alignment_ref = None
        singlerun_output_info = {}
        report_info = {'name': None, 'ref': None}
        ret_val = None

        rds_files = list()
        ret_fwd = reads_info["file_fwd"]
        if ret_fwd is not None:
            rds_files.append(ret_fwd)
            if reads_info.get('file_rev', None) is not None:
                rds_files.append(reads_info['file_rev'])

        # 2. After all is set, do the alignment and upload the output.
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

            # Upload the alignment
            upload_results = self.star_utils.upload_STARalignment(
                                        input_params,
                                        rds,
                                        reads_info,
                                        output_bam_file)
            alignment_ref = upload_results['obj_ref']
            alignment_obj = {
                'ref': alignment_ref,
                'name': rds['alignment_output_name']
            }
            alignment_objs.append({
                'reads_ref': rds['ref'],
                'AlignmentObj': alignment_obj
            })

            singlerun_output_info['output_dir'] = star_mp_ret['star_output']
            singlerun_output_info['output_bam_file'] = output_bam_file
            singlerun_output_info['upload_results'] = upload_results

            if input_params.get("create_report", 0) == 1:
                report_info = self.star_utils.generate_report_for_single_run(singlerun_output_info, input_params)

            ret_val = {'alignmentset_ref': None,
                    'output_directory': singlerun_output_info['output_dir'],
                    'output_info': singlerun_output_info,
                    'alignment_objs': alignment_objs,
                    'report_name': report_info['name'],
                    'report_ref': report_info['ref']
            }
        else:
            ret_val = {'alignmentset_ref': None,
                    'output_directory': None,
                    'output_info': None,
                    'alignment_objs': None,
                    'report_name': None,
                    'report_ref': None
            }

        if ret_fwd is not None:
            os.remove(ret_fwd)
            if reads_info.get('file_rev', None) is not None:
                os.remove(reads_info["file_rev"])

        return ret_val


    def star_run_batch_sequential(self, input_params):
        """
        star_run_batch_sequential: running the STAR align by looping
        """
        log('--->\nrunning STAR_Aligner.star_run_batch_sequential\n' +
                'params:\n{}'.format(json.dumps(input_params, indent=1)))

        self.get_index(input_params)

        reads_refs = input_params[STARUtils.SET_READS]

        single_input_params = copy.deepcopy(input_params)

        # 1. Run the mapping one by one
        alignment_items = []
        alignment_objs = []
        rds_names = []
        for r in reads_refs:
            single_input_params[STARUtils.PARAM_IN_READS] = r['ref']
            single_input_params['create_report'] = 0
            single_ret = self.star_run_single(single_input_params)

            item = single_ret['alignment_objs'][0]
            a_obj = item['AlignmentObj']
            r_ref = item['reads_ref']
            alignment_objs.append(item)
            alignment_items.append({
                    'ref': a_obj['ref'],
                    'label': r.get(
                            'condition',
                            single_input_params.get('condition','unspecified'))
            })

            rds_names.append(r['alignment_output_name'].replace(single_input_params['alignment_suffix'], ''))

        # 2. Process all the results after mapping is done
        (batch_result, report_info) =self._batch_sequential_post_processing(
                                        alignment_items, rds_names, input_params)

        batch_result['output_directory'] = self.star_out_dir

        result = {'alignmentset_ref': result_obj_ref,
                'output_info': batch_result,
                'alignment_objs': alignment_objs,
                'report_name': report_info['name'],
                'report_ref': report_info['ref']
        }

        return result


    def _batch_sequential_post_processing(self, alignment_items, rds_names, params):
        '''
        process the mapping results of all the reads in the readsset_ref
        '''
        # 1. Save the alignment set
        set_name_map = self.star_utils.get_object_names([params[STARUtils.PARAM_IN_READS]])
        set_name = set_name_map[params[STARUtils.PARAM_IN_READS]]

        output_alignmentset_name = set_name + params['alignmentset_suffix']

        save_result = self.star_utils.upload_alignment_set(
                        alignment_items,
                        output_alignmentset_name,
                        params['output_workspace'])

        result_obj_ref = save_result['set_ref']

        # 2. Extract the ReadsPerGene counts if necessary
        output_dir = os.path.join(self.scratch, STARUtils.STAR_OUT_DIR)
        if params.get(STARUtils.PARAM_IN_OUTFILE_PREFIX, None) is not None:
            prefix = params[STARUtils.PARAM_IN_OUTFILE_PREFIX]
        else:
            prefix = ''

        gene_count_files = []
        if (params.get('quantMode', None) is not None
                    and (params['quantMode'] == 'Both'
                            or 'GeneCounts' in params['quantMode'])):
            for reads_name in rds_names:
                gene_count_files.append('{0}/{1}ReadsPerGene.out.tab'.format(reads_name, prefix))

            extract_geneCount_matrix(gene_count_files, output_dir)

        # 3. Reporting...
        report_info = {'name': None, 'ref': None}

        #run qualimap
        qualimap_report = self.qualimap.run_bamqc({'input_ref': result_obj_ref})
        qc_result_zip_info = qualimap_report['qc_result_zip_info']
        qc_result = [{'shock_id': qc_result_zip_info['shock_id'],
                      'name': qc_result_zip_info['index_html_file_name'],
                      'label': qc_result_zip_info['name']}]

        # create the report
        report_text = 'Ran on SampleSet or ReadsSet.\n\n'
        report_text += 'Created ReadsAlignmentSet: ' + str(output_alignmentset_name) + '\n\n'

        report_info = self.star_utils._generate_star_report(
                        result_obj_ref,
                        report_text,
                        qc_result,
                        params['output_workspace'],
                        output_dir)

        return (save_result, report_info)


    def star_run_batch_parallel(self, input_params):
        """
        star_run_batch_parallel: running the STAR align in batch parallelly
        """
        log('--->\nrunning STAR_Aligner.star_run_batch_parallel\n' +
                'params:\n{}'.format(json.dumps(input_params, indent=1)))

        reads_refs = input_params[STARUtils.SET_READS]
        # build task list and send it to KBParallel
        tasks = []
        for r in reads_refs:
            tasks.append(
                    self.build_single_execution_task(
                        r['ref'], input_params)
                    )

        batch_run_params = {'tasks': tasks,
                            'runner': 'parallel',
                            'max_retries': 2}

        if input_params.get('concurrent_local_tasks', None) is not None:
                batch_run_params['concurrent_local_tasks'] = input_params['concurrent_local_tasks']
        if input_params.get('concurrent_njsw_tasks', None) is not None:
                batch_run_params['concurrent_njsw_tasks'] = input_params['concurrent_njsw_tasks']

        results = self.parallel_runner.run_batch(batch_run_params)
        print('Batch run results=')
        pprint(results)

        batch_result = self.star_utils.process_batch_result(results, input_params, reads_refs)
        batch_result['output_directory'] = self.star_out_dir

        return batch_result


    def process_batch_result(self, batch_result, params, reads_refs):
        n_jobs = len(batch_result['results'])
        n_success = 0
        n_error = 0
        ran_locally = 0
        ran_njsw = 0

        set_name_map = self.star_utils.get_object_names([params[self.PARAM_IN_READS]])
        set_name = set_name_map[params[self.PARAM_IN_READS]]

        # reads alignment set items
        alignment_items = []
        alignment_objs = []
        rds_names = []

        for k in range(0, len(batch_result['results'])):
            reads_ref = reads_refs[k]
            rds_names.append(reads_ref['alignment_output_name'].replace(params['alignment_suffix'], ''))

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
        save_result = self.star_utils.upload_alignment_set(
                        alignment_items,
                        output_alignmentset_name,
                        params['output_workspace'])

        result_obj_ref = save_result['set_ref']

        output_dir = os.path.join(self.working_dir, STARUtils.STAR_OUT_DIR)
        if params.get(self.PARAM_IN_OUTFILE_PREFIX, None) is not None:
            prefix = params[self.PARAM_IN_OUTFILE_PREFIX]
        else:
            prefix = ''

        # Extract the ReadsPerGene counts if necessary
        gene_count_files = []
        if (params.get('quantMode', None) is not None
                    and (params['quantMode'] == 'Both'
                            or 'GeneCounts' in params['quantMode'])):
            for reads_name in rds_names:
                gene_count_files.append('{}/{}ReadsPerGene.out.tab'.format(reads_name, prefix))

            extract_geneCount_matrix(gene_count_files, output_dir)

        # Reporting...
        report_info = {'name': None, 'ref': None}

        #run qualimap
        qualimap_report = self.qualimap.run_bamqc({'input_ref': result_obj_ref})
        qc_result_zip_info = qualimap_report['qc_result_zip_info']
        qc_result = [{'shock_id': qc_result_zip_info['shock_id'],
                      'name': qc_result_zip_info['index_html_file_name'],
                      'label': qc_result_zip_info['name']}]

        # create the report
        report_text = 'Ran on SampleSet or ReadsSet.\n\n'
        report_text += 'Created ReadsAlignmentSet: ' + str(output_alignmentset_name) + '\n\n'
        report_text += 'Total ReadsLibraries = ' + str(n_jobs) + '\n'
        report_text += '        Successful runs = ' + str(n_success) + '\n'
        report_text += '            Failed runs = ' + str(n_error) + '\n'
        report_text += '       Ran on main node = ' + str(ran_locally) + '\n'
        report_text += '   Ran on remote worker = ' + str(ran_njsw) + '\n\n'

        report_info = self.star_utils._generate_star_report(
                        result_obj_ref,
                        report_text,
                        qc_result,
                        params['output_workspace'],
                        output_dir)

        result = {'alignmentset_ref': result_obj_ref,
                'output_info': batch_result,
                'alignment_objs': alignment_objs,
                'report_name': report_info['name'],
                'report_ref': report_info['ref']
        }

        return result

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

    def get_version_from_subactions(self, module_name, subactions):
        # go through each sub action looking for
        if not subactions:
            return 'dev' #'release'  # default to release if we can't find anything
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


    def run_star_indexing(self, input_params):
        """
        Runs STAR in genomeGenerate mode to build the index files and directory for STAR mapping.
        It creates a directory as defined by self.star_idx_dir in the scratch area that houses
        the index files.
        """
        ret_params = copy.deepcopy(input_params)
        ret_params[STARUtils.PARAM_IN_STARMODE] = 'genomeGenerate'

        # build the indexing parameters
        params_idx = self.star_utils._get_indexing_params(ret_params, self.star_idx_dir)

        ret = 1
        try:
            if ret_params[STARUtils.PARAM_IN_STARMODE]=='genomeGenerate':
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

        return (ret, params_idx[STARUtils.STAR_IDX_DIR])


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
            retVal = {'star_idx': self.star_idx_dir, 'star_output': None}
        else:#no exception raised by STAR mapping and STAR returns 0, then move to saving and reporting  
            retVal = {'star_idx': self.star_idx_dir, 'star_output': params_mp.get('align_output')}

        return retVal

    def get_index(self, input_params):
        '''
        get_index: generate the index if not yet existing
        '''
        gnm_ref = input_params[STARUtils.PARAM_IN_GENOME]
	if input_params.get('sjdbGTFfile', None) is None:
            input_params['sjdbGTFfile'] = self.star_utils._get_genome_gtf_file(
                                        gnm_ref,
                                        self.star_idx_dir)

        if not os.path.isfile(os.path.join(self.star_idx_dir, 'genomeParameters.txt')):
            # fetch genome fasta and GTF from refs to file location(s)
            input_params[STARUtils.PARAM_IN_FASTA_FILES] = self.star_utils._get_genome_fasta(gnm_ref)

            # generate the indices
            (idx_ret, idx_dir) = self.run_star_indexing(input_params)
            if idx_ret != 0:
                raise ValueError("Failed to generate genome indices, aborting...")
