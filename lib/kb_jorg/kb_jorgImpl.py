# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import json

from kb_jorg.Utils.JorgUtil import JorgUtil

#END_HEADER


class kb_jorg:
    '''
    Module Name:
    kb_jorg

    Module Description:
    A KBase module: kb_jorg
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "1.0.0"
    GIT_URL = "https://github.com/jungbluth/kb_jorg"
    GIT_COMMIT_HASH = "a0b04ac1561ebbed05152b2a9b4fdf7d8c019940"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        #END_CONSTRUCTOR
        pass

    def run_kb_jorg(self, ctx, params):
        """
        :param params: instance of type "JorgInputParams" (required
           params: assembly_ref: Genome assembly object reference
           binned_contig_name: BinnedContig object name and output file
           header workspace_name: the name of the workspace it gets saved to.
           reads_list: list of reads object
           (PairedEndLibrary/SingleEndLibrary) upon which CONCOCT will be run
           optional params: thread: number of threads; default 1
           min_contig_length: minimum contig length; default 1000
           contig_split_length: length to split long contigs; default 10000
           ref: https://github.com/BinPro/CONCOCT) -> structure: parameter
           "assembly_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "binned_contig_name" of String, parameter
           "workspace_name" of String, parameter "reads_list" of list of type
           "obj_ref" (An X/Y/Z style reference), parameter "thread" of Long,
           parameter "min_contig_length" of Long, parameter
           "contig_split_length" of Long
        :returns: instance of type "JorgResult" (result_folder: folder
           path that holds all files generated by run_kb_jorg report_name:
           report name generated by KBaseReport report_ref: report reference
           generated by KBaseReport) -> structure: parameter
           "result_directory" of String, parameter "binned_contig_obj_ref" of
           type "obj_ref" (An X/Y/Z style reference), parameter "report_name"
           of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN run_kb_jorg

        print('--->\nRunning kb_jorg.kb_jorg\nparams:')
        print(json.dumps(params, indent=1))

        for key, value in params.items():
            if isinstance(value, str):
                params[key] = value.strip()

        jorg_runner = JorgUtil(self.config)

        returnVal = jorg_runner.run_jorg_wf(params)

        #END run_kb_jorg

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method run_kb_jorg return value ' +
                             'returnVal is not type dict as required.')
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
