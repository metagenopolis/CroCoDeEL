process buildSamplesList {
    input:
    file mgsProfiles

    output:
    file '*.csv'

    script:
    """
    python3 $params.project_dir/pythonScript/build_samples_list.py -i "$mgsProfiles"
    """
}

process crocodeelProcess{
    input:
    val(sample)

    output:
    path "*.txt"

    script:
    """
    python3 $params.project_dir/pythonScript/crocodeel.py -i $params.mgsProfilesPath -s $sample
    """
}

process writeResultsProcess{
    publishDir "${params.outputDirectoryPath}", mode: 'copy'

    input:
    file "contamination_results"

    output:
    file "*.txt"

    script:
    """
    cat contamination_results > contamination_results.txt
    """
}

process plotContaminationCasesProcess {
    publishDir "${params.outputDirectoryPath}", mode: 'copy'

    input:
    file "contamination_results"

    output:
    file "*.pdf"

    script:
    """
    python3 $params.project_dir/pythonScript/plot_contamination_cases.py -i $params.mgsProfilesPath -o $params.outputDirectoryPath
    """
}

process writeParametersFileProcess {
    publishDir "${params.outputDirectoryPath}", mode: 'copy'

    input:
    file "contamination_results"

    output:
    file "*.ini"

    script:
    """
    python3 $params.project_dir/pythonScript/write_parameters_file.py -i $params.mgsProfilesPath -o $params.outputDirectoryPath
    """
}

process writeSummaryReportProcess {
    publishDir "${params.outputDirectoryPath}", mode: 'copy'

    input:
    file "contamination_results"

    output:
    file "*.ini"

    script:
    """
    python3 $params.project_dir/pythonScript/write_summary_report.py -i $params.mgsProfilesPath -o $params.outputDirectoryPath
    """
}

process filterContaminationsByPlateProcess {
    publishDir "${params.outputDirectoryPath}", mode: 'copy'

    input:
    file "contamination_results"

    output:
    file "*.txt"

    script:
    """
    python3 $params.project_dir/pythonScript/filter_contamination_results_by_plate.py -o $params.outputDirectoryPath -p $params.plateMapPath
    """
}

process filterContaminationsIfLongitudinalSamplesProcess {
    publishDir "${params.outputDirectoryPath}", mode: 'copy'

    input:
    file "contamination_results"

    output:
    file "*.txt"

    script:
    """
    python3 $params.project_dir/pythonScript/filter_contamination_results_if_longitudinal_samples.py -o $params.outputDirectoryPath -m $params.longitudinalMetadataPath
    """
}

workflow CroCoDeEL {
    mgsProfiles = Channel.fromPath(params.mgsProfilesPath)
    contamination_results = buildSamplesList(mgsProfiles)| splitCsv(header:true)| map { row-> row.sample } | crocodeelProcess | \
    collectFile(name: 'results.txt', keepHeader:true, skip:1) | writeResultsProcess
    plotContaminationCasesProcess(contamination_results) 
    writeParametersFileProcess(contamination_results)
    writeSummaryReportProcess(contamination_results)
    if (params.plateMapPath != false & params.longitudinalMetadataPath == false){
        filterContaminationsByPlateProcess(contamination_results)
    }
    else if (params.longitudinalMetadataPath != false & params.plateMapPath == false){
        filterContaminationsIfLongitudinalSamplesProcess(contamination_results)
    }
    else if (params.longitudinalMetadataPath != false & params.plateMapPath != false){
        filterContaminationsByPlateProcess(contamination_results) | filterContaminationsIfLongitudinalSamplesProcess
    }
}