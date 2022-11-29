#!/bin/bash

i=1
for PARAM in participant cpus metadata
# folder_path=20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid2000-nb_clones2-nb_mut300
do i=$((i+1)) ; eval $PARAM=$(grep "^$PBS_ARRAYID," $INPUT_FILE | cut -d',' -f$i) ; done

export R_LIBS_USER='/data/users/ahamype1/BC_BILAT/dl_tools'
export PATH="/bioinfo/local/build/Centos/R/R-3.6.0/bin:/bioinfo/local/build/samtools-1.1/bin:$PATH"

R --file=./src/NGS_preprocess_server/script/superfreq/runSuperFreqOnDonor.R --args $participant $cpus $metadata

