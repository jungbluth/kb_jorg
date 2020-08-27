/*
A KBase module: kb_jorg
*/

module kb_jorg {

    /* A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int boolean;

    /* An X/Y/Z style reference
    */
    typedef string obj_ref;

    /*
        required params:
        assembly_ref: Genome assembly object
        reads_file: reads object (PairedEndLibrary/SingleEndLibrary) upon which jorg will be run
        workspace_name: the name of the workspace it gets saved to.
        output_assembly_name: name of the output assembly

        optional params:
        read_mapping_tool: tool to use for read mapping
        kmer_size: specify kmer length for baiting
        min_coverage: minimum coverage value
        num_iterations: specify a number of iterations to use
        high_contig_num: are there more than 2500 contigs and is this expected
        single_end_reads: is the input library single-end type
        assembly_selection_criteria: criteria to select the final assembly output from Jorg
        circle_min_overlap_length: specify overlap length when checking for circularized contig
        inspect_all_contigs_using_last: check all contigs in the final selected assembly for circularization using LAST
        check_assembly_using_pilon: check for misassembled contigs using Pilon
        auto_stop_near_runtime_limit: limit number of iterations to avoid 7-day run limit associated errors on KBase cluster

        ref: https://github.com/jungbluth/jorg

    */
    typedef structure {
        obj_ref assembly_ref;
        string workspace_name;
        list<obj_ref> reads_file;
        string output_assembly_name;

        string read_mapping_tool;
        int kmer_size;
        int min_coverage;
        int num_iterations;
        string high_contig_num;
        string single_end_reads;
        string assembly_selection_criteria;
        int circle_min_overlap_length;
        string inspect_all_contigs_using_last;
        string check_assembly_using_pilon;
        string auto_stop_near_runtime_limit;

    } jorgInputParams;

    /*
        result_directory: folder path that holds all files generated by run_kb_jorg
        report_name: report name generated by KBaseReport
        report_ref: report reference generated by KBaseReport
    */
    typedef structure{
        string result_directory;
        obj_ref assembly_obj_ref;
        string report_name;
        string report_ref;
    }jorgResult;

    funcdef run_kb_jorg(jorgInputParams params)
        returns (jorgResult returnVal) authentication required;

};
