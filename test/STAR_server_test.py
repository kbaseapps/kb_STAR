# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import os.path
import json  # noqa: F401
import time
import requests
import shutil

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from STAR.STARImpl import STAR
from STAR.STARServer import MethodContext
from STAR.authclient import KBaseAuth as _KBaseAuth

from biokbase.AbstractHandle.Client import AbstractHandle as HandleService
from ReadsUtils.baseclient import ServerError
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil

class STARTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('STAR'):
            cls.cfg[nameval[0]] = nameval[1]
            
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        
        # WARNING: don't call any logging methods on the context object,        
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'STAR',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = STAR(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        
        cls.readUtilsImpl = ReadsUtils(cls.callback_url, token=cls.token)
        cls.staged = {}

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_STAR_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    def load_fasta_file(self, filename, obj_name, contents):
        f = open(filename, 'w')
        f.write(contents)
        f.close()
        assemblyUtil = AssemblyUtil(self.callback_url)
        assembly_ref = assemblyUtil.save_assembly_from_fasta({'file': {'path': filename},
                                                              'workspace_name': self.getWsName(),
                                                              'assembly_name': obj_name
                                                              })
        return assembly_ref

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    # Uncomment to skip this test
    # @unittest.skip("skipped test_run_star")
    def test_run_star(self):
        # First load a test FASTA file as an KBase Assembly
        # get the test data
        pe_lib_info = self.getPairedEndLibInfo()
        pprint(pe_lib_info)

        # STAR parameters
        params = {
            'workspace_name': self.getWsName(),
            'output_contigset_name': 'STAR_test_contigset',
            'hash_length': 21,
            'read_libraries':[self.make_ref(pe_lib_info)],
            'min_contig_length': 500,
            'cov_cutoff': 5.2,
            'read_trkg': 0,
            'amos_file': 'yes',
            'exp_cov': 21.3,
            'ins_length': 400
        }

        # Second, call your implementation
        result = self.getImpl().run_star(self.getContext(), params)

        if not result[0]['report_ref'] is None:
                rep = self.wsClient.get_objects2({'objects': [{'ref': result[0]['report_ref']}]})['data'][0]
                print('REPORT object:')
                pprint(rep)

                self.assertEqual(rep['info'][1].rsplit('_', 1)[0], 'kb_velvet_report')
                self.assertEqual(rep['info'][2].split('-', 1)[0], 'KBaseReport.Report')
        else:
                print('Velvet failed!')

        # Second, call your implementation
        ret = self.getImpl().filter_contigs(self.getContext(),
                                            {'workspace_name': self.getWsName(),
                                             'assembly_input_ref': assembly_ref,
                                             'min_length': 10
                                             })
