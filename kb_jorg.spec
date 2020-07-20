/*
A KBase module: kb_jorg
*/

module kb_jorg {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_kb_jorg(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
