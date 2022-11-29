
###########################################
# DNA sequencing analysis initial cohort #
###########################################

### step 1: pipeline of the bioinformatics platform to go from the fastq to
### processed BAM + mutect variant calling for cancer samples + GATK haplotype
### caller for the germline samples


### step 2: run SuperFreq v. 1.3.2, at the patient level and
### at the patient-side-level
for i in {1..6} ; do echo $j,P$i,2,superfreq_metadata_patient.csv >> results/superFreq/run_superfreq_patient.csv; j=$((j+1));  done
mkdir log/run_superfreq
qsub -t 11-18 -N run_superfreq -q batch -d $PWD -l walltime=12:00:00,mem=100gb,nodes=1:ppn=2 -o log/run_superfreq -e log/run_superfreq -v INPUT_FILE=results/superFreq/run_superfreq_patient.csv /data/kdi_prod/project_result/1118/02.00/src/NGS_preprocess_server/script/superfreq/run_superfreq.sh

for i in {1..6} ; do echo $j,P$i,2,superfreq_metadata_left.csv >> results/superFreq/run_superfreq_left.csv; j=$((j+1));  done
mkdir log/run_superfreq_left
qsub -t 11-18 -N run_superfreq -q batch -d $PWD -l walltime=12:00:00,mem=100gb,nodes=1:ppn=2 -o log/run_superfreq_left -e log/run_superfreq_left -v INPUT_FILE=results/superFreq/run_superfreq_left.csv /data/kdi_prod/project_result/1118/02.00/src/NGS_preprocess_server/script/superfreq/run_superfreq.sh

for i in {1..6} ; do echo $j,P$i,2,superfreq_metadata_left.csv >> results/superFreq/run_superfreq_left.csv; j=$((j+1));  done
mkdir log/run_superfreq_left
qsub -t 11-18 -N run_superfreq -q batch -d $PWD -l walltime=12:00:00,mem=100gb,nodes=1:ppn=2 -o log/run_superfreq_left -e log/run_superfreq_left -v INPUT_FILE=results/superFreq/run_superfreq_left.csv /data/kdi_prod/project_result/1118/02.00/src/NGS_preprocess_server/script/superfreq/run_superfreq.sh

# requires file the script run_superfreq.sh and runSuperFreqOnDonor.R



### step 3: manual review of common left and right mutations
### with IGV
# step 3A - extract bam