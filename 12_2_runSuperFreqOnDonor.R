library(superFreq)
library(seqinr)

#take the INDIVIDUAL and cpu as input.
participant = as.character(commandArgs(TRUE)[1])
cpus = as.numeric(commandArgs(TRUE)[2])
metadata = as.character(commandArgs(TRUE)[3])

superFreq(
          metaDataFile = paste('/data/kdi_prod/project_result/1118/02.00/results/superFreq/', metadata, sep=''),  #all individuals in one metaData file.
          genome='hg19',
          normalDirectory='/data/kdi_prod/project_result/1118/02.00/results/superFreq/ReferenceNormals',
          Rdirectory='/data/kdi_prod/project_result/1118/02.00/results/superFreq/R',
          plotDirectory='/data/kdi_prod/project_result/1118/02.00/results/superFreq/plots',  #each individual will get a subdirectory in here
          reference='/data/kdi_prod/project_result/1118/02.00/results/superFreq/genomes/hg19.fa',
          mode='exome',
          cpus=cpus,
          captureRegions='/data/kdi_prod/project_result/1118/02.00/results/superFreq/captureRegions/SureSelect_Clinical_Research_Exome_Regions.bed',
          participants=participant  #run only on this individual
          )