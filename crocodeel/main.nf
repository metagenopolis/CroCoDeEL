/**
Nextflow workflow for the detection of cross-sample contamination
**/

nextflow.enable.dsl=2

include {  CroCoDeEL } from './workflow/crocodeel'

log.info """\

            ======================================
             C R O C O D E E L    W O R K F L O W
            ======================================
Cross-sample Contamination Detection and Estimation of its Levels

Options:

    - mgsProfiles       : ${params.mgsProfilesPath}
    - outputDirectoryPath       : ${params.outputDirectoryPath}

 """

workflow {
     CroCoDeEL()
}