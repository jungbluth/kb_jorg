import errno
import json
import os
import subprocess
import sys
import time
import uuid
import zipfile
import copy
import shutil
import math
import glob
import pandas as pd

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.ReadsUtilsClient import ReadsUtils

from random import seed
from random import randint
# seed random number generator
seed(1)

from shutil import copyfile


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

        reads_file_check = isinstance(reads_file, list)
        if reads_file_check :
            log("Input reads_file detected as list. Great.")
        else:
            log("Input reads_file not a list, converting.")
            reads_file = [reads_file]


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


    def retrieve_assembly(self, task_params):
        if os.path.exists(task_params['contig_file_path']):
            assembly = task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            log("task_params['assembly_ref'] is {}".format(task_params['assembly_ref']))
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

    # not used right now because Jorg does not support single-end mode
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

        # sorted_bam_file_list = []

        # list of reads files, can be 1 or more. assuming reads are either type unpaired or interleaved
        # will not handle unpaired forward and reverse reads input as seperate (non-interleaved) files

        for i in range(len(read_scratch_path)):
            fastq = read_scratch_path[i]
            fastq_type = read_type[i]

            sam = os.path.basename(fastq).split('.fastq')[0] + ".sam"
            #sam = os.path.join(self.JORG_RESULT_DIRECTORY, sam)

            if fastq_type == 'interleaved':  # make sure working - needs tests
                log("Running interleaved read mapping mode")
                self.run_read_mapping_interleaved_pairs_mode(task_params, assembly, fastq, sam)
            else:  # running read mapping in single-end mode
                log("Running unpaired read mapping mode")
                self.run_read_mapping_unpaired_mode(task_params, assembly, fastq, sam)

            sorted_bam = self.convert_sam_to_sorted_and_indexed_bam(sam)

        return sorted_bam

    def generate_make_coverage_table_command(self, task_params, sorted_bam):
        # create the depth file for this bam
        #
        depth_file_path = os.path.join(self.scratch, str('depth.txt'))
        command = '/kb/module/lib/kb_jorg/bin/jgi_summarize_bam_contig_depths '
        command += '--outputDepth {} '.format(depth_file_path)
        command += '--minContigDepth 1 {}'.format(sorted_bam)

        log('running summarize_bam_contig_depths command: {}'.format(command))
        self._run_command(command)

        return depth_file_path

    def check_input_assembly_for_minimum_coverage(self, task_params):
        path_to_depth_file = os.path.join(self.scratch, str('depth.txt'))
        file1 = open(path_to_depth_file, 'r')
        lines = file1.readlines()
        detected_coverage_in_longest_contig = float(lines[1].split()[2])
        if detected_coverage_in_longest_contig >= float(task_params["min_coverage"]):
            log("The longest contig in the input assembly passes the input minimum coverage criteria (value: {}). Setting Jorg coverage parameter to 75% of this value.".format(detected_coverage_in_longest_contig))
            jorg_working_coverage = math.floor(float(lines[1].split()[2]) * 0.75)
            log("Working Jorg minimum coverage set to {}.".format(jorg_working_coverage))
        else:
            log("The longest contig in the input assembly has a coverage value of {}, which does not meet the minimum coverage criteria. Exiting Jorg.".format(float(lines[1].split()[2])))
            sys.exit(1)
        return jorg_working_coverage


    def fix_generate_jorg_command_ui_bug(self, task_params):
        # needed to get checkbox for UI to work with string objects, for some
        # reason strings are converted to numerics when running inside KBase UI.
        parameter_high_contig_num = task_params['high_contig_num']

        if task_params['high_contig_num'] is 1:
            parameter_high_contig_num = '--high_contig_num yes'
        elif task_params['high_contig_num'] is 0:
            parameter_high_contig_num = '--high_contig_num no'

        return parameter_high_contig_num


    def copy_required_jorg_input_files_to_cwd(self):
        manifest_template_file_source = "/Jorg/Example/manifest_template.conf"
        manifest_template_file_destination = os.path.join(self.scratch, str('manifest_template.conf'))
        shutil.move(manifest_template_file_source,manifest_template_file_destination)


## not working correctly in narrative
    def process_jorg_iteration_output_and_calc_stats(self):
        path_to_iterations_file = os.path.join(self.scratch, str("iterations.txt"))
        path_to_iterations_flat_file = os.path.join(self.scratch, str("iterations_flat.txt"))

        file1 = open(path_to_iterations_file, 'r')
        lines = file1.readlines()

        # used during debugging
        # log("start printing Jorg log")
        # datafile = glob.glob('Jorg*.log')[0]
        # z = open(datafile, "r")  # the a opens it in append mode
        # text = z.read()
        # log(text)
        # z.close()
        # log("end printing Jorg log")

        genome_num_fasta = []
        last_circle_check = []
        contig_name = []
        contig_length = []
        contig_GC_percent = []
        cumulative_length = []
        running_longest_single_fragment = 0
        running_longest_assembly_length = 0
        assembly_with_longest_single_fragment = 0
        assembly_with_longest_cumulative_assembly_length = 0
        i = -1
        j = -1
        with open(path_to_iterations_flat_file, 'a') as f:
            for line in lines:
                if line.startswith('Iteration'):
                    i += 1
                    genome_num_fasta.append(line.split()[1] + ".fasta")
                    last_circle_check.append(line.split()[1] + ".reduced")
                else:
                    if not line.startswith('contig_name'):
                        j += 1
                        contig_name.append(line.split()[0])
                        contig_length.append(line.split()[1])
                        contig_GC_percent.append(line.split()[2])
                        cumulative_length.append(line.split()[3])
                        if (int(contig_length[j]) > running_longest_single_fragment):
                            running_longest_single_fragment = int(contig_length[j])
                            assembly_with_longest_single_fragment = genome_num_fasta[i]
                        if (int(cumulative_length[j]) > running_longest_assembly_length):
                            running_longest_assembly_length = int(cumulative_length[j])
                            assembly_with_longest_cumulative_assembly_length = genome_num_fasta[i]
                        f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(genome_num_fasta[i],last_circle_check[i],contig_name[j],contig_length[j],contig_GC_percent[j],cumulative_length[j]))
        f.close()
        final_iteration_assembly = genome_num_fasta[len(genome_num_fasta)-1]
        log("running_longest_single_fragment {}, assembly_with_longest_single_fragment {}, assembly_with_longest_cumulative_assembly_length {}, final_iteration_assembly {}".format(running_longest_single_fragment, assembly_with_longest_single_fragment, assembly_with_longest_cumulative_assembly_length, final_iteration_assembly))
        return (running_longest_single_fragment, assembly_with_longest_single_fragment, assembly_with_longest_cumulative_assembly_length, final_iteration_assembly)


    def select_jorg_output_genome(self, task_params, running_longest_single_fragment, assembly_with_longest_single_fragment, assembly_with_longest_cumulative_assembly_length, final_iteration_assembly):
        if task_params["assembly_selection_criteria"] == "longest_single_fragment":
            output_jorg_assembly = assembly_with_longest_single_fragment
            log("Assembly {} selected as output assembly.".format(output_jorg_assembly))
        elif task_params["assembly_selection_criteria"] == "longest_single_fragment_filter":
            output_jorg_assembly = assembly_with_longest_single_fragment
            log("Assembly {} selected as output assembly.".format(output_jorg_assembly))
        elif task_params["assembly_selection_criteria"] == "longest_cumulative_assembly_length":
            output_jorg_assembly = assembly_with_longest_cumulative_assembly_length
            log("Assembly {} selected as output assembly.".format(output_jorg_assembly))
        elif task_params["assembly_selection_criteria"] == "longest_cumulative_assembly_length_filter":
            output_jorg_assembly = assembly_with_longest_cumulative_assembly_length
            log("Assembly {} selected as output assembly.".format(output_jorg_assembly))
        elif task_params["assembly_selection_criteria"] == "final_iteration_assembly":
            output_jorg_assembly = final_iteration_assembly
            log("Assembly {} selected as output assembly.".format(output_jorg_assembly))
        log("jorg assembly selected: {}".format(output_jorg_assembly))
        num_contigs = 0
        full_path_output_jorg_assembly = "Iterations/" + output_jorg_assembly
        file1 = open(full_path_output_jorg_assembly, 'r')
        lines = file1.readlines()
        for line in lines:
            if line.startswith('>'):
                num_contigs = num_contigs + 1
        output_jorg_assembly_name = os.path.basename(output_jorg_assembly).split(".fasta")[0]
        return output_jorg_assembly, output_jorg_assembly_name, num_contigs


    def uppercase_fastq_file(self, reads_file):
        output_fastq = str(reads_file) + "_uppercase.fastq"
        command = 'seqkit seq -u '
        command += '{} > '.format(reads_file)
        command += '{}'.format(output_fastq)
        log('uppercase_fastq_file: {}'.format(command))
        self._run_command(command)
        return output_fastq

    def clean_input_fasta(self, output_jorg_assembly):
        output_jorg_assembly_clean = str(output_jorg_assembly) + "_clean.fasta" # need to fix split command below, sloppy fix
        command = 'cut -d\' \' -f1 Iterations/{} > {}'.format(output_jorg_assembly, output_jorg_assembly_clean)
        log('clean_input_fasta: {}'.format(command))
        self._run_command(command)
        return output_jorg_assembly_clean

    def sort_fasta_by_length(self, input_fasta):
        output_fasta_sorted = input_fasta.rsplit('.', 1)[0] + "_sorted.fasta"
        command = 'seqkit sort '
        command += '-l {} '.format(input_fasta)
        command += '-r '.format(input_fasta)
        command += '> {} '.format(output_fasta_sorted)
        log('sort_fasta_by_length: {}'.format(command))
        self._run_command(command)
        return output_fasta_sorted

    def extract_mapping_tracks_from_bam(self, sorted_bam):
        output_sam = 'final.mapped.sam'
        command = 'bedtools genomecov -ibam '
        command += '{} '.format(sorted_bam)
        command += '-bg > circos_mapping_tracks.txt'
        log('extract_mapping_tracks_from_bam: {}'.format(command))
        self._run_command(command)
        file1 = open(os.path.abspath("circos_mapping_tracks.txt"), 'r')
        pandas_df = pd.read_table(file1, header=None)
        max_cov = round(pandas_df[3].max(),1)
        min_cov = round(pandas_df[3].min(),1)
        std_cov = round(pandas_df[3].std(),1)
        mean_cov = round(pandas_df[3].mean(),1)
        return max_cov, min_cov, std_cov, mean_cov

    def make_circos_karyotype_file(self, output_jorg_assembly_clean_sorted):
        from Bio import SeqIO
        count = 1
        path_to_circos_karyotype_file = "circos_karyotype.txt"
        with open(path_to_circos_karyotype_file, 'a') as f:
            for record in SeqIO.parse(output_jorg_assembly_clean_sorted, "fasta"):
                f.write("chr - {} {} 0 {} {}\n".format(record.id, count, len(record), record.id))
                count += 1
        f.close()

    def prep_circos_axis(self, max_cov):
        if max_cov < 30:
            max_cov = 30
        command = 'sed -i "s/^max.*/max   = {}/" /kb/module/lib/kb_jorg/circos/circos.conf '.format(max_cov)
        log('prep_circos_axis: {}'.format(command))
        self._run_command(command)

    def draw_circos_plot(self):
        command = 'circos -conf '
        command += '/kb/module/lib/kb_jorg/circos/circos.conf'
        log('draw_circos_plot: {}'.format(command))
        self._run_command(command)

    def make_circos_plot(self, task_params, reads_file, output_jorg_assembly):
        output_fastq = self.uppercase_fastq_file(reads_file)
        output_jorg_assembly_clean = self.clean_input_fasta(output_jorg_assembly)
        output_jorg_assembly_clean_sorted = self.sort_fasta_by_length(output_jorg_assembly_clean)
        sam = os.path.basename(output_fastq).split('.fastq')[0] + ".sam"
        self.run_read_mapping_interleaved_pairs_mode(task_params, output_jorg_assembly_clean_sorted, output_fastq, sam)
        sorted_bam = self.convert_sam_to_sorted_and_indexed_bam(sam)
        max_cov, min_cov, std_cov, mean_cov = self.extract_mapping_tracks_from_bam(sorted_bam)
        self.make_circos_karyotype_file(output_jorg_assembly_clean_sorted)
        self.prep_circos_axis(max_cov)
        self.draw_circos_plot()
        return output_jorg_assembly_clean_sorted, max_cov, min_cov, std_cov, mean_cov

    def run_circle_check_using_last(self, output_jorg_assembly):
        command = 'bash {}/lib/circle_check_using_last '.format(self.JORG_BASE_PATH)
        command += 'Iterations/{} '.format(output_jorg_assembly)
        self._run_command(command)
        path_to_last_output = output_jorg_assembly.split(".fasta")[0] + ".reduced"
        return path_to_last_output

    def process_last_output(self, task_params, path_to_last_output):
        file1 = open(path_to_last_output, 'r')
        remember_query = ""
        forward_circle_match = ""
        reverse_circle_match = ""
        lines = file1.readlines()
        circularized_contigs = []
        query_start = 0
        query_end = 0
        subject_start = 0
        subject_end = 0
        minimum_length_threshold = task_params['circle_min_overlap_length']
        for line in lines:
            if line.startswith('#'):
                pass
            else:
                if not remember_query == line.split()[0]: # clear this variable if different contig
                    remember_query = ""
                    longest_contig_length = 0
                remember_query = line.split()[0]
                if line.split()[0] == line.split()[1]: # if query and subject the same
                    if int(line.split()[14]) >= minimum_length_threshold: # if the match is above the minimum length threshold
                        if line.split()[2] == "100.00": # if a 100% match identified
                            if float(line.split()[10]) <= 10e-5: # if the expected value is significant
                                if line.split()[13] == line.split()[14]: # if length of query equals length of reference
                                    longest_contig_length = int(line.split()[14])
                                if forward_circle_match == "TRUE":
                                    if (query_start == line.split()[8]) and (query_end == line.split()[9]) and (subject_start == line.split()[6]) and (subject_end == line.split()[7]) and (int(subject_end) == int(longest_contig_length)) and (int(line.split()[7]) == int(longest_contig_length)):
                                        reverse_circle_match = "TRUE"
                                        circularized_contigs.append(line.split()[0])
                                elif line.split()[13] != line.split()[14]: # if length of query not equal to length of reference
                                    if (int(line.split()[7]) - int(line.split()[6])) == (int(line.split()[9]) - int(line.split()[8])): # if query start/end length equal to subject start/end length
                                        forward_circle_match = "TRUE"
                                        remember_query = line.split()[0]
                                        query_start = line.split()[6]
                                        query_end = line.split()[7]
                                        subject_start = line.split()[8]
                                        subject_end = line.split()[9]
        if len(circularized_contigs) >= 1:
            output_circle_text = "Yes! Contig(s): "
            for i in range(len(circularized_contigs)):
                output_circle_text = output_circle_text + str(circularized_contigs[i])
                if i != (len(circularized_contigs) - 1):
                    output_circle_text = output_circle_text + ","
        else:
            output_circle_text = "No"
        return circularized_contigs, output_circle_text

    def move_jorg_output_files_to_output_dir(self):
        dest = os.path.abspath(self.JORG_RESULT_DIRECTORY)
        files = os.listdir(os.path.abspath(self.scratch))
        for f in files:
            if (f.startswith("iterations") or f.startswith("Jorg") or f.startswith("manifest") or f.startswith("mira") or f.startswith("mirabait") or f.startswith("list") or f.startswith("depth") or f.startswith("circos") or f.endswith(".log") or f.endswith("sorted.fasta") or f.endswith("clean.fasta") or f.endswith(".reduced") or f.endswith(".tbl")):
                shutil.move(f, dest)

    def run_jorg_and_circos_workflow(self, task_params, jorg_working_coverage):
        """
        generate_command: jorg
        """

        assembly_ref = task_params['contig_file_path']
        reads_file = task_params['reads_list_file'][0]
        kmer_size = task_params['kmer_size']
        min_coverage = jorg_working_coverage
        num_iterations = task_params['num_iterations']
        parameter_high_contig_num = self.fix_generate_jorg_command_ui_bug(task_params)

        log("\n\nRunning run_jorg_and_circos_workflow")
        command = 'bash {}/jorg '.format(self.JORG_BASE_PATH)
        command += '--bin_fasta_file {} '.format(assembly_ref)
        command += '--reads_file {} '.format(reads_file)
        command += '--kmer_length {} '.format(kmer_size)
        command += '--min_coverage {} '.format(min_coverage)
        command += '--iterations {} '.format(num_iterations)
        command += '--runtime_cap 6 ' # runtime limit (in days) for running on KBase
        command += ' {}'.format(parameter_high_contig_num)
        log('Generated jorg command: {}'.format(command))
        self.copy_required_jorg_input_files_to_cwd()
        log("start running Jorg command")
        self._run_command(command)
        log("end running Jorg command")

        # process jorg output and calculate statistics
        running_longest_single_fragment, assembly_with_longest_single_fragment, assembly_with_longest_cumulative_assembly_length, final_iteration_assembly = self.process_jorg_iteration_output_and_calc_stats()

        #output_jorg_assembly = 'Iterations/1.fasta'
        output_jorg_assembly, output_jorg_assembly_name, num_contigs = self.select_jorg_output_genome(task_params, running_longest_single_fragment, assembly_with_longest_single_fragment, assembly_with_longest_cumulative_assembly_length, final_iteration_assembly)

        # run check for circularity
        path_to_last_output = self.run_circle_check_using_last(output_jorg_assembly)

        circularized_contigs, output_circle_text = self.process_last_output(task_params, path_to_last_output)

        # make circos plot
        output_jorg_assembly_clean_sorted, max_cov, min_cov, std_cov, mean_cov = self.make_circos_plot(task_params, reads_file, output_jorg_assembly)

        # move relevant files to output directory provided to user
        self.move_jorg_output_files_to_output_dir()

        return output_jorg_assembly_clean_sorted, output_jorg_assembly_name, num_contigs, output_circle_text, max_cov, min_cov, std_cov, mean_cov


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
                # for file in files:
                #     zip_file.write(os.path.join(dirname, file), file)
                if (dirname.endswith(self.JORG_RESULT_DIRECTORY)):
                    baseDir = os.path.basename(dirname)
                    for file in files:
                        full = os.path.join(dirname, file)
                        zip_file.write(full, os.path.join(baseDir, file))
        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'Files generated by Jorg App'})
        return output_files


    def generate_html_report(self, output_assembly_name, assembly_ref, assembly_stats):
        """
        generate_html_report: generate html summary report
        """
        log('Start generating html report')
        #html_report = list()
        output_directory = os.path.join(self.scratch, 'html_dir_' + str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        # get summary data from existing assembly object and bins_objects
        Summary_Table_Content = ''
        Overview_Content = ''

        # generate overview content

        # get png
        png_filename_l = [f for f in os.listdir(self.JORG_RESULT_DIRECTORY) if f.endswith('.png')]

        # Example
        # Overview_Content += '<p>Input contigs: {}</p>'.format(input_contig_count)
        Overview_Content += '<p>Iteration Selected: {}</p>'.format(assembly_stats['iteration'])
        Overview_Content += '<p>Number of contigs: {}</p>'.format(assembly_stats['num_contigs'])
        Overview_Content += '<p>Coverage (avg, sd, max, min): {}, {}, {}, {}</p>'.format(assembly_stats['mean_cov'],assembly_stats['std_cov'],assembly_stats['max_cov'],assembly_stats['min_cov'])
        Overview_Content += '<p>Circularized genome: {}</p>'.format(assembly_stats['circle_or_not'])
        for png_filename in png_filename_l:
            Overview_Content += '\n<embed src="{}" width="700px" height="700px">'.format(png_filename)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                report_template = report_template.replace('Summary_Table_Content',
                                                          Summary_Table_Content)
                result_file.write(report_template)

        # copy pdfs into html dir
        for png_filename in png_filename_l:
            shutil.copyfile(os.path.join(self.JORG_RESULT_DIRECTORY, png_filename), os.path.join(output_directory, png_filename))

        # save html dir to shock
        def dir_to_shock(dir_path, name, description):
            '''
            For regular directories or html directories
            name - for regular directories: the name of the flat (zip) file returned to ui
                   for html directories: the name of the html file
            '''
            dfu_fileToShock_ret = self.dfu.file_to_shock({
                'file_path': dir_path,
                'make_handle': 0,
                'pack': 'zip',
                })
            dir_shockInfo = {
                'shock_id': dfu_fileToShock_ret['shock_id'],
                'name': name,
                'description': description
                }
            return dir_shockInfo
        html_shockInfo = dir_to_shock(output_directory, 'report.html', 'HTML report for Jorg')
        """
        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for kb_concoct App'})
        return html_report
        """
        return [html_shockInfo]

    def generate_report(self, assembly_ref_obj, params, assembly_stats):
        """
        generate_report: generate summary report
        """
        log('Generating report')
        params['result_directory'] = self.JORG_RESULT_DIRECTORY
        output_files = self.generate_output_file_list(params['result_directory'])
        output_html_files = self.generate_html_report(params['result_directory'],params['output_assembly_name'],assembly_stats)
        report_params = {
              'message': '',
              'workspace_name': params.get('workspace_name'),
              'file_links': output_files,
              'html_links': output_html_files,
              'direct_html_link_index': 0,
              'html_window_height': 500,
              'report_object_name': 'kb_jorg_report_' + str(uuid.uuid4())}
        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)
        report_output = {'report_name': output['name'], 'report_ref': output['ref']}
        return report_output


    def run_jorg(self, task_params):
        """
        run_jorg: jorg app
        """
        log('--->\nrunning JorgUtil.run_jorg\n' +
            'task_params:\n{}'.format(json.dumps(task_params, indent=1)))

        # light validation on input parameters
        self._validate_run_jorg_params(task_params)

        # get assembly
        contig_file = self._get_contig_file(task_params['assembly_ref'])
        task_params['contig_file_path'] = contig_file

        # clean the assembly file so that there are no spaces in the fasta headers
        assembly = self.retrieve_assembly(task_params)
        task_params['contig_file_path'] = assembly

        # get reads
        (read_scratch_path, read_type) = self.stage_reads_file(task_params['reads_file'])
        task_params['read_type'] = read_type
        task_params['reads_list_file'] = read_scratch_path

        # prep result directory
        result_directory = os.path.join(self.scratch, self.JORG_RESULT_DIRECTORY)
        self._mkdir_p(result_directory)

        # map reads to determine input coverage
        sorted_bam = self.generate_alignment_bams(task_params, assembly)

        # extract depth information from bam files
        depth_file_path = self.generate_make_coverage_table_command(task_params, sorted_bam)

        # check to make sure input contigs have require coverage
        jorg_working_coverage = self.check_input_assembly_for_minimum_coverage(task_params)

        # run jorg and circos
        output_jorg_assembly_clean_sorted, output_jorg_assembly_name, num_contigs, output_circle_text, max_cov, min_cov, std_cov, mean_cov = self.run_jorg_and_circos_workflow(task_params, jorg_working_coverage)

        assembly_stats = {'iteration' : output_jorg_assembly_name, 'num_contigs' : num_contigs, 'mean_cov' : mean_cov, 'std_cov' : std_cov, 'min_cov' : min_cov, 'max_cov' : max_cov, 'circle_or_not' : output_circle_text}

        # file handling and management

        # log('Saved result files to: {}'.format(result_directory))
        # log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))

        dest = os.path.abspath(self.JORG_RESULT_DIRECTORY)

        assembly_ref_obj = self.au.save_assembly_from_fasta(
            {'file': {'path': dest + '/' + output_jorg_assembly_clean_sorted},
#            {'file': {'path': os.path.abspath(output_jorg_assembly_clean_sorted)},
             'workspace_name': task_params['workspace_name'],
             'assembly_name': task_params['output_assembly_name']
             })

        # generate report
        reportVal = self.generate_report(assembly_ref_obj, task_params, assembly_stats)
        returnVal = {
            'result_directory': result_directory,
            'assembly_obj_ref': assembly_ref_obj
        }

        returnVal.update(reportVal)
        return returnVal
