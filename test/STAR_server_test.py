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

from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from GenomeFileUtil.GenomeFileUtilClient import GenomeFileUtil
from ReadsUtils.baseclient import ServerError
from ReadsUtils.ReadsUtilsClient import ReadsUtils

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

    
    def loadAssembly(self):
        if hasattr(self.__class__, 'assembly_ref'):
            return self.__class__.assembly_ref
        fasta_path = os.path.join(self.scratch, 'test.fna')
        shutil.copy(os.path.join('data', 'test.fna'), fasta_path)
        au = AssemblyUtil(self.callback_url)
        assembly_ref = au.save_assembly_from_fasta({'file': {'path': fasta_path},
                                                    'workspace_name': self.getWsName(),
                                                    'assembly_name': 'test_assembly'
                                                    })
        self.__class__.assembly_ref = assembly_ref
        return assembly_ref

    def loadGenome(self):
        if hasattr(self.__class__, 'genome_ref'):
            return self.__class__.genome_ref
        genbank_file_path = os.path.join(self.scratch, 'minimal.gbff')
        shutil.copy(os.path.join('data', 'minimal.gbff'), genbank_file_path)
        gfu = GenomeFileUtil(self.callback_url)
        genome_ref = gfu.genbank_to_genome({'file': {'path': genbank_file_path},
                                            'workspace_name': self.getWsName(),
                                            'genome_name': 'test_genome'
                                            })['genome_ref']
        self.__class__.genome_ref = genome_ref
        return genome_ref


    def loadPEReads(self):
        if hasattr(self.__class__, 'assembly_ref'):
            return self.__class__.assembly_ref
        fasta_path = os.path.join(self.scratch, 'test.fna')
        shutil.copy(os.path.join('data', 'test.fna'), fasta_path)
        au = AssemblyUtil(self.callback_url)
        assembly_ref = au.save_assembly_from_fasta({'file': {'path': fasta_path},
                                                    'workspace_name': self.getWsName(),
                                                    'assembly_name': 'test_assembly'
                                                    })
        self.__class__.assembly_ref = assembly_ref
        return assembly_ref

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
        params = {'command': 'STAR', 'options': ['--help']}
        self.getImpl().run_star(self.getContext(), params)

    # Uncomment to skip this test
    # @unittest.skip("skipped test_build_star_index_from_assembly")
    def test_build_star_index_from_assembly(self):

        # test build directly from an assembly, forget to add ws_for_cache so object will not be cached
        assembly_ref = self.loadAssembly()
        res = self.getImpl().get_star_index(self.getContext(), {'ref': assembly_ref})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 0)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 0)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_assembly')

        pprint(res)

        # do it again, and set ws_for_cache
        assembly_ref = self.loadAssembly()
        res = self.getImpl().get_star_index(self.getContext(), {'ref': assembly_ref,
                                                                   'ws_for_cache': self.getWsName()})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 0)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 1)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_assembly')

        pprint(res)

        # do it again, should retrieve from cache
        assembly_ref = self.loadAssembly()
        res = self.getImpl().get_star_index(self.getContext(), {'ref': assembly_ref})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 1)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 0)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_assembly')
        pprint(res)

    # Uncomment to skip this test
    # @unittest.skip("skipped test_build_star_index_from_genome")
    def test_build_star_index_from_genome(self):

        # finally, try it with a genome_ref instead
        genome_ref = self.loadGenome()
        res = self.getImpl().get_star_index(self.getContext(), {'ref': genome_ref})[0]
        self.assertIn('output_dir', res)
        self.assertIn('from_cache', res)
        self.assertEquals(res['from_cache'], 0)
        self.assertIn('pushed_to_cache', res)
        self.assertEquals(res['pushed_to_cache'], 0)
        self.assertIn('index_files_basename', res)
        self.assertEquals(res['index_files_basename'], 'test_genome_assembly')
        pprint(res)

    def loadSingleEndReads(self):
        if hasattr(self.__class__, 'se_reads_ref'):
            return self.__class__.se_reads_ref
        fq_path = os.path.join(self.scratch, 'reads_1_se.fq')
        shutil.copy(os.path.join('data', 'bt_test_data', 'reads_1.fq'), fq_path)

        ru = ReadsUtils(self.callback_url)
        se_reads_ref = ru.upload_reads({'fwd_file': fq_path,
                                        'wsname': self.getWsName(),
                                        'name': 'test_assembly',
                                        'sequencing_tech': 'artificial reads'})['obj_ref']
        self.__class__.se_reads_ref = se_reads_ref
        return se_reads_ref


    def loadPairedEndReads(self):
        if hasattr(self.__class__, 'pe_reads_ref'):
            return self.__class__.pe_reads_ref
        fq_path1 = os.path.join(self.scratch, 'reads_1.fq')
        shutil.copy(os.path.join('data', 'bt_test_data', 'reads_1.fq'), fq_path1)
        fq_path2 = os.path.join(self.scratch, 'reads_2.fq')
        shutil.copy(os.path.join('data', 'bt_test_data', 'reads_2.fq'), fq_path2)

        ru = ReadsUtils(self.callback_url)
        pe_reads_ref = ru.upload_reads({'fwd_file': fq_path1, 'rev_file': fq_path2,
                                        'wsname': self.getWsName(),
                                        'name': 'test_assembly',
                                        'sequencing_tech': 'artificial reads'})['obj_ref']
        self.__class__.pe_reads_ref = pe_reads_ref
        return pe_reads_ref




    def test_star_aligner(self):
        self.loadAssembly()
        self.loadSingleEndReads()
        self.loadPairedEndReads()
