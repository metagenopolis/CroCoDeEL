#!/bin/bash

nextflow run main.nf \
    -w /export/mgps/home/lgoulet/scratch_global/tmp_nextflow_v3 \
    --mgsProfilesPath '/export/mgps/home/lgoulet/scratch_global/crocodeel/atraction_subset_mgs_profiles.tsv' \
    --outputDirectory '/export/mgps/home/lgoulet/scratch_global/crocodeel/results_crocodeel_subset_test' \
    --plateMapPath '/export/mgps/home/lgoulet/scratch_global/crocodeel/crocodeel/crocodeel/data/atraction_test/plate_map_test.txt' \
    --longitudinalMetadataPath '/export/mgps/home/lgoulet/scratch_global/crocodeel/crocodeel/crocodeel/data/atraction_test/longitudinal_test.txt'


