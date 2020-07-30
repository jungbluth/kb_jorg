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

        optional params:
        read_mapping_tool: tool to use for read mapping
        kmer_size: specify kmer length for baiting
        min_coverage: minimum coverage value
        num_iterations: specify a number of iterations to use
        ref: https://github.com/jungbluth/jorg

    */
    typedef structure {
        obj_ref assembly_ref;
        string workspace_name;
        obj_ref reads_file;

        string read_mapping_tool;
        int kmer_size;
        int min_coverage;
        int num_iterations;
        int num_iterations;
        int num_iterations;
        string high_contig_num;
        string single_end_reads;

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
