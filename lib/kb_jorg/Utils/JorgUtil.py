import errno
import json
import os
import subprocess
import sys
import time
import uuid
import zipfile
import copy

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.ReadsUtilsClient import ReadsUtils

from random import seed
from random import randint
# seed random number generator
seed(1)


# for future expansion
# from kb_jorg.BinningUtilities import BinningUtil as bu


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class JorgUtil:
    JORG_BASE_PATH = '/Jorg'
    JORG_RESULT_DIRECTORY = 'jorg_output_dir'
    MAPPING_THREADS = 16
    BBMAP_MEM = '30g'

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        self.ws_url = config['workspace-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.mgu = MetagenomeUtils(self.callback_url)

    def _validate_run_jorg_params(self, task_params):
        """
        _validate_run_jorg_params:
                validates params passed to run_jorg method
        """
        log('Start validating run_jorg params')

        # check for required parameters
        for p in ['assembly_ref', 'reads_file', 'workspace_name']:
            if p not in task_params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

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

    def _run_command(self, command):
        """
        _run_command: run command and print result
        """
        os.chdir(self.scratch)
        log('Start executing command:\n{}'.format(command))
        log('Command is running from:\n{}'.format(self.scratch))
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output, stderr = pipe.communicate()
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\n'.format(exitCode))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}\nStderr:\n{}'.format(exitCode, output, stderr)
            raise ValueError(error_msg)
            sys.exit(1)
        return (output, stderr)

    # this function has been customized to return read_type variable (interleaved vs single-end library)
    def stage_reads_file(self, reads_file):
        """
        stage_reads_file: download fastq file associated to reads to scratch area
                          and return result_file_path
        """

        log('Processing reads object list: {}'.format(reads_file))

        result_file_path = []
        read_type = []

        # getting from workspace and writing to scratch. The 'reads' dictionary now has file paths to scratch.
        reads = self.ru.download_reads({'read_libraries': reads_file, 'interleaved': None})['files']

        # reads_file is the list of file paths on workspace? (i.e. 12804/1/1).
        # "reads" is the hash of hashes where key is "12804/1/1" or in this case, read_obj and
        # "files" is the secondary key. The tertiary keys are "fwd" and "rev", as well as others.
        for read_obj in reads_file:
            files = reads[read_obj]['files']    # 'files' is dictionary where 'fwd' is key of file path on scratch.
            result_file_path.append(files['fwd'])
            read_type.append(files['type'])
            if 'rev' in files and files['rev'] is not None:
                result_file_path.append(files['rev'])

        return result_file_path, read_type

    def _get_contig_file(self, assembly_ref):
        """
        _get_contig_file: get contig file from GenomeAssembly object
        """
        contig_file = self.au.get_assembly_as_fasta({'ref': assembly_ref}).get('path')

        sys.stdout.flush()
        contig_file = self.dfu.unpack_file({'file_path': contig_file})['file_path']

        return contig_file
    #
    def retrieve_assembly(self, task_params):
        if os.path.exists(task_params['contig_file_path']):
            assembly = task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            assembly = self._get_contig_file(task_params['assembly_ref'])
        return assembly

    def deinterlace_raw_reads(self, fastq):
        fastq_forward = fastq.split('.fastq')[0] + "_forward.fastq"
        fastq_reverse = fastq.split('.fastq')[0] + "_reverse.fastq"
        command = 'deinterleave_fastq.sh < {} {} {}'.format(fastq, fastq_forward, fastq_reverse)
        try:
            self._run_command(command)
        except:
            raise Exception("Cannot deinterlace fastq file!")
        return (fastq_forward, fastq_reverse)

    def run_read_mapping_interleaved_pairs_mode(self, task_params, assembly, fastq, sam):
        read_mapping_tool = task_params['read_mapping_tool']
        log("running {} mapping in interleaved mode.".format(read_mapping_tool))
        random_seed_int = randint(0, 999999999)
        log("randomly selected seed (integer) used for read mapping is: {}".format(random_seed_int))
        if task_params['read_mapping_tool'] == 'bbmap':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'fast interleaved=true mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bowtie2_default':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            bt2index = os.path.basename(assembly) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 -x {} '.format(bt2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'bowtie2_very_sensitive':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            bt2index = os.path.basename(assembly) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 --very-sensitive -x {} '.format(bt2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'minimap2':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            command = 'minimap2 -ax sr -t {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} '.format(assembly)
            command += '{} '.format(fastq_forward)
            command += '{} > '.format(fastq_reverse)
            command += '{}'.format(sam)
        elif task_params['read_mapping_tool'] == 'hisat2':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            ht2index = os.path.basename(assembly) + '.ht2'
            command = 'hisat2-build {} '.format(assembly)
            command += '{} && '.format(ht2index)
            command += 'hisat2 -x {} '.format(ht2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '-S {} '.format(sam)
            command += '--seed {} '.format(random_seed_int)
            command += '--threads {}'.format(self.MAPPING_THREADS)
        log('running alignment command: {}'.format(command))
        out, err = self._run_command(command)
    #
    # def run_read_mapping_unpaired_mode(self, task_params, assembly, fastq, sam):
    #     read_mapping_tool = task_params['read_mapping_tool']
    #     log("running {} mapping in single-end (unpaired) mode.".format(read_mapping_tool))
    #     random_seed_int = randint(0, 999999999)
    #     log("randomly selected seed (integer) used for read mapping is: {}".format(random_seed_int))
    #     if task_params['read_mapping_tool'] == 'bbmap':
    #         log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
    #         command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
    #         command += 'threads={} '.format(self.MAPPING_THREADS)
    #         command += 'ref={} '.format(assembly)
    #         command += 'in={} '.format(fastq)
    #         command += 'out={} '.format(sam)
    #         command += 'fast interleaved=false mappedonly nodisk overwrite'
    #         # BBMap is deterministic without the deterministic flag if using single-ended reads
    #     elif task_params['read_mapping_tool'] == 'bowtie2_default':
    #         bt2index = os.path.basename(assembly) + '.bt2'
    #         command = 'bowtie2-build -f {} '.format(assembly)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '{} && '.format(bt2index)
    #         command += 'bowtie2 -x {} '.format(bt2index)
    #         command += '-U {} '.format(fastq)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '-S {}'.format(sam)
    #     elif task_params['read_mapping_tool'] == 'bowtie2_very_sensitive':
    #         bt2index = os.path.basename(assembly) + '.bt2'
    #         command = 'bowtie2-build -f {} '.format(assembly)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '{} && '.format(bt2index)
    #         command += 'bowtie2 --very-sensitive -x {} '.format(bt2index)
    #         command += '-U {} '.format(fastq)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '-S {}'.format(sam)
    #     elif task_params['read_mapping_tool'] == 'minimap2':
    #         command = 'minimap2 -ax sr -t {} '.format(self.MAPPING_THREADS)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '{} '.format(assembly)
    #         command += '{} > '.format(fastq)
    #         command += '{}'.format(sam)
    #     elif task_params['read_mapping_tool'] == 'hisat2':
    #         ht2index = os.path.basename(assembly) + '.ht2'
    #         command = 'hisat2-build {} '.format(assembly)
    #         command += '{} && '.format(ht2index)
    #         command += 'hisat2 -x {} '.format(ht2index)
    #         command += '-U {} '.format(fastq)
    #         command += '-S {} '.format(sam)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '--threads {}'.format(self.MAPPING_THREADS)
    #     log('running alignment command: {}'.format(command))
    #     out, err = self._run_command(command)

    def convert_sam_to_sorted_and_indexed_bam(self, sam):
        # create bam files from sam files
        sorted_bam = os.path.abspath(sam).split('.sam')[0] + "_sorted.bam"

        command = 'samtools view -F 0x04 -uS {} | '.format(sam)
        command += 'samtools sort - -o {}'.format(sorted_bam)

        log('running samtools command to generate sorted bam: {}'.format(command))
        self._run_command(command)

        # verify we got bams
        if not os.path.exists(sorted_bam):
            log('Failed to find bam file\n{}'.format(sorted_bam))
            sys.exit(1)
        elif(os.stat(sorted_bam).st_size == 0):
            log('Bam file is empty\n{}'.format(sorted_bam))
            sys.exit(1)

        # index the bam file
        command = 'samtools index {}'.format(sorted_bam)

        log('running samtools command to index sorted bam: {}'.format(command))
        self._run_command(command)

        return sorted_bam

    def generate_alignment_bams(self, task_params, assembly):
        """
            This function runs the selected read mapper and creates the
            sorted and indexed bam files from sam files using samtools.
        """

        reads_file = task_params['reads_file']

        (read_scratch_path, read_type) = self.stage_reads_file(reads_file)

        sorted_bam_file_list = []

        # list of reads files, can be 1 or more. assuming reads are either type unpaired or interleaved
        # will not handle unpaired forward and reverse reads input as seperate (non-interleaved) files

        for i in range(len(read_scratch_path)):
            fastq = read_scratch_path[i]
            fastq_type = read_type[i]

            sam = os.path.basename(fastq).split('.fastq')[0] + ".sam"
            sam = os.path.join(self.JORG_RESULT_DIRECTORY, sam)

            if fastq_type == 'interleaved':  # make sure working - needs tests
                log("Running interleaved read mapping mode")
                self.run_read_mapping_interleaved_pairs_mode(task_params, assembly, fastq, sam)
            else:  # running read mapping in single-end mode
                log("Running unpaired read mapping mode")
                self.run_read_mapping_unpaired_mode(task_params, assembly, fastq, sam)

            sorted_bam = self.convert_sam_to_sorted_and_indexed_bam(sam)

            sorted_bam_file_list.append(sorted_bam)

        return sorted_bam_file_list

    def generate_make_coverage_table_command(self, task_params, sorted_bam_file_list):
        # create the depth file for this bam
        #
        min_contig_length = task_params['min_contig_length']
        sorted_bam = task_params['sorted_bam']

        depth_file_path = os.path.join(self.scratch, str('jorg_depth.txt'))
        command = '/kb/module/lib/kb_jorg/bin/jgi_summarize_bam_contig_depths '
        command += '--outputDepth {} '.format(depth_file_path)
        command += '--minContigLength {} '.format(min_contig_length)
        command += '--minContigDepth 1 {}'.format(sorted_bam)

        log('running summarize_bam_contig_depths command: {}'.format(command))
        self._run_command(command)

        return depth_file_path

    def fix_generate_jorg_command_ui_bug(self, task_params):
        # needed to get checkbox for UI to work with string objects, for some
        # reason strings are converted to numerics when running inside KBase UI.
        parameter_high_contig_num = task_params['high_contig_num']
        parameter_single_end_reads = task_params['single_end_reads']

        if task_params['high_contig_num'] is 1:
            parameter_high_contig_num = '--high_contig_num yes'
        elif task_params['high_contig_num'] is 0:
            parameter_high_contig_num = '--high_contig_num no'

        if task_params['single_end_reads'] is 1:
            parameter_single_end_reads = '--single_end_reads yes'
        elif task_params['single_end_reads'] is 0:
            parameter_single_end_reads = '--single_end_reads no'

        return (parameter_high_contig_num, parameter_single_end_reads)

# *******************
    # def generate_jorg_coverage_table_from_bam(self, task_params):
    #     """
    #     generate_command: jorg generate coverage table
    #     """
    #     log("\n\nRunning generate_jorg_coverage_table_from_bam")
    #     command = 'python {}/scripts/jorg_coverage_table.py temp.bed '.format(self.JORG_BASE_PATH)
    #     command += '{}/*_sorted.bam > '.format(self.JORG_RESULT_DIRECTORY)
    #     command += '{}/coverage_table.tsv'.format(self.JORG_RESULT_DIRECTORY)
    #     log('Generated jorg generate coverage table from bam command: {}'.format(command))
    #
    #     self._run_command(command)
# *******************

    def generate_jorg_command(self, task_params):
        """
        generate_command: jorg
        """

        assembly_ref = task_params['assembly_ref']
        reads_file = task_params['reads_file']
        kmer_size = task_params['kmer_size']
        min_coverage = task_params['min_coverage']
        num_iterations = task_params['num_iterations']
        parameter_high_contig_num, parameter_single_end_reads = \
            self.fix_generate_jorg_command_ui_bug(task_params)

        log("\n\nRunning generate_jorg_command")
        command = '{}/jorg '.format(self.JORG_BASE_PATH)
        command += '--bin_fasta_file {} '.format(assembly_ref)
        command += '--reads_file {} '.format(reads_file)
        command += '--kmer_length {} '.format(kmer_size)
        command += '--min_coverage {} '.format(min_coverage)
        command += '--iterations {} '.format(num_iterations)
        command += ' {} '.format(parameter_high_contig_num)
        command += ' {}'.format(parameter_single_end_reads)
        log('Generated jorg command: {}'.format(command))

        self._run_command(command)


    def generate_output_file_list(self, result_directory):
        """
        generate_output_file_list: zip result files and generate file_links for report
        """
        log('Start packing result files')
        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'jorg_result.zip')

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:

            for dirname, subdirs, files in os.walk(result_directory):
                for file in files:
                    if (file.endswith('.sam') or
                        file.endswith('.bam') or
                        file.endswith('.bai') or
                       file.endswith('.summary')):
                            continue
                    if (dirname.endswith(self.BINNER_BIN_RESULT_DIR)):
                            continue
                    zip_file.write(os.path.join(dirname, file), file)
                if (dirname.endswith(self.BINNER_BIN_RESULT_DIR)):
                    baseDir = os.path.basename(dirname)
                    for file in files:
                        full = os.path.join(dirname, file)
                        zip_file.write(full, os.path.join(baseDir, file))

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'Files generated by JORG App'})

        return output_files

    def generate_html_report(self, result_directory, assembly_ref, binned_contig_obj_ref):
        """
        generate_html_report: generate html summary report
        """

        log('Start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        # get summary data from existing assembly object and bins_objects
        Summary_Table_Content = ''
        Overview_Content = ''
        (binned_contig_count, input_contig_count, total_bins_count) = \
            self.generate_overview_info(assembly_ref, binned_contig_obj_ref, result_directory)

        Overview_Content += '<p>Binned contigs: {}</p>'.format(binned_contig_count)
        Overview_Content += '<p>Input contigs: {}</p>'.format(input_contig_count)
        Overview_Content += '<p>Number of bins: {}</p>'.format(total_bins_count)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                report_template = report_template.replace('Summary_Table_Content',
                                                          Summary_Table_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for kb_jorg App'})
        return html_report

    def generate_overview_info(self, assembly_ref, binned_contig_obj_ref, result_directory):
        """
        _generate_overview_info: generate overview information from assembly and binnedcontig
        """

        # get assembly and binned_contig objects that already have some data populated in them
        assembly = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]
        binned_contig = self.dfu.get_objects({'object_refs': [binned_contig_obj_ref]})['data'][0]

        input_contig_count = assembly.get('data').get('num_contigs')
        binned_contig_count = 0
        total_bins_count = 0
        total_bins = binned_contig.get('data').get('bins')
        total_bins_count = len(total_bins)
        for bin in total_bins:
            binned_contig_count += len(bin.get('contigs'))

        return (binned_contig_count, input_contig_count, total_bins_count)

    def generate_report(self, binned_contig_obj_ref, task_params):
        """
        generate_report: generate summary report
        """
        log('Generating report')

        result_directory = os.path.join(self.scratch, "jorg_output_dir")

        task_params['result_directory'] = result_directory

        output_files = self.generate_output_file_list(task_params['result_directory'])

        output_html_files = self.generate_html_report(task_params['result_directory'],
                                                      task_params['assembly_ref'],
                                                      binned_contig_obj_ref)

        report_params = {
            'message': '',
            'workspace_name': task_params['workspace_name'],
            'file_links': output_files,
            'html_links': output_html_files,
            'direct_html_link_index': 0,
            'html_window_height': 266,
            'report_object_name': 'kb_jorg_report_' + str(uuid.uuid4())
        }

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def create_dict_from_depth_file(self, depth_file_path):
        # keep contig order (required by metabat2)
        depth_file_dict = {}
        with open(depth_file_path, 'r') as f:
            header = f.readline().rstrip().split("\t")
            # print('HEADER1 {}'.format(header))
            # map(str.strip, header)
            for line in f:
                # deal with cases were fastq name has spaces.Assume first
                # non white space word is unique and use this as ID.
                # line = line.rstrip()
                vals = line.rstrip().split("\t")
                if ' ' in vals[0]:
                    ID = vals[0].split()[0]
                else:
                    ID = vals[0]
                depth_file_dict[ID] = vals[1:]
            depth_file_dict['header'] = header
        return depth_file_dict

    def run_jorg(self, task_params):
        """
        run_jorg: jorg app

        required params:
            assembly_ref: Metagenome assembly object reference
            binned_contig_name: BinnedContig object name and output file header
            workspace_name: the name of the workspace it gets saved to.
            reads_file: list of reads object (PairedEndLibrary/SingleEndLibrary)
            upon which JORG will be run

        optional params:
            min_contig_length: minimum contig length; default 1000

            ref: https://github.com/BinPro/JORG/blob/develop/README.md
        """
        log('--->\nrunning JorgUtil.run_jorg\n' +
            'task_params:\n{}'.format(json.dumps(task_params, indent=1)))

        self._validate_run_jorg_params(task_params)

        # get assembly
        contig_file = self._get_contig_file(task_params['assembly_ref'])
        task_params['contig_file_path'] = contig_file

        # clean the assembly file so that there are no spaces in the fasta headers
        assembly = self.retrieve_assembly(task_params)
        task_params['contig_file_path'] = assembly

        # get reads
        (reads_list_file, read_type) = self.stage_reads_file(task_params['reads_file'])
        task_params['read_type'] = read_type
        task_params['reads_list_file'] = reads_list_file

        # prep result directory
        result_directory = os.path.join(self.scratch, self.JORG_RESULT_DIRECTORY)
        self._mkdir_p(result_directory)

        cwd = os.getcwd()
        log('changing working dir to {}'.format(result_directory))
        os.chdir(result_directory)



        # set up tasks for kbparallel to run alignments
        # this also submits run_alignments function in parallel
        # self.set_up_parallel_tasks(params)

        # run alignments, and update input contigs to use the clean file
        # this function has an internal loop to generate a sorted bam file for each input read file
        #
        # self.set_up_parallel_tasks(task_params)

        self.generate_alignment_bams(task_params, assembly)

        # not used right now
        depth_file_path = self.generate_make_coverage_table_command(task_params, sorted_bam_file_list)
        # depth_dict = self.create_dict_from_depth_file(depth_file_path)

        # run jorg prep, cut up fasta input
        self.generate_jorg_cut_up_fasta_command(task_params)

        # run jorg make coverage table from bam
        self.generate_jorg_coverage_table_from_bam(task_params)

        # run jorg prep and jorg
        self.generate_jorg_command(task_params)

        # run jorg post cluster merging command
        self.generate_jorg_post_clustering_merging_command(task_params)

        # run extract bins command
        self.generate_jorg_extract_fasta_bins_command(task_params)

        # run fasta renaming
        self.rename_and_standardize_bin_names(task_params)

        # revert fasta headers in bins
        self.revert_fasta_headers(task_params)

        self.make_binned_contig_summary_file_for_binning_apps(task_params)

        # file handling and management
        os.chdir(cwd)
        log('changing working dir to {}'.format(cwd))

        log('Saved result files to: {}'.format(result_directory))
        log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))

        # make new BinnedContig object and upload to KBase
        generate_binned_contig_param = {
            'file_directory': os.path.join(result_directory, self.BINNER_BIN_RESULT_DIR),
            'assembly_ref': task_params['assembly_ref'],
            'binned_contig_name': task_params['binned_contig_name'],
            'workspace_name': task_params['workspace_name']
        }

        binned_contig_obj_ref = \
            self.mgu.file_to_binned_contigs(generate_binned_contig_param).get('binned_contig_obj_ref')

        # generate report
        reportVal = self.generate_report(binned_contig_obj_ref, task_params)
        returnVal = {
            'result_directory': result_directory,
            'binned_contig_obj_ref': binned_contig_obj_ref
        }
        returnVal.update(reportVal)

        return returnVal