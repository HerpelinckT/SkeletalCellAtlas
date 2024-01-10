#!/bin/bash 

# A: Tim Herpelinck, 2023
# D: builds SCA cellranger pipeline by downloading required files, setting up the conda environment and building the reference genome.
# R: adapted from Sikkema et al., Nat. Med. (2023), https://rdcu.be/dvxX0; GitHub: https://github.com/LungCellAtlas/scRNAseq_pipelines

# Pipeline version
pipeline_version="1.0.0"


initialize_scRNA_seq() {

# DEFAULT PARAMETERS

# Ensembl release
ensrel="???"
# Genome
genomestring="GRCm?????"
# Species
species="mus_musculus"
# Expected cellranger version
expectedcrv="3.1.0"
# Resources for cellranger mkref
# number of cores
nthreads="20"
# Memory (GB)
memgb="12"
#script parts to run/skip:
download_files=true
install_env=false
build_ref_genome=true
download_ensembl_files=true

usage() {
    cat <<HELP_USAGE

# TBD!!!!!!!!!

HELP_USAGE
}

# TAKE IN ARGUMENTS

env_path="test"
work_dir="."






# CREATE LOG AND RECORD PARAMETERS

# store path to script location as script_dir
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# cd to workdir and store full path (no trailing slash)
cd $work_dir
work_dir=`pwd`

# define logfile
LOGFILE=$work_dir/"LOG_SCA_pipeline_setup.log"

# check if logfile exists
if [ -f ${LOGFILE} ]; then
    echo "ERROR: Log file ${LOGFILE} already exists. Remove to continue. Exiting..."
    exit 1
fi

# create log file
echo `date` > ${LOGFILE}
echo "Log file saved in: ${LOGFILE}"
# print pipeline version
echo "Skeletal Cell Atlas pipeline version: v${pipeline_version}" | tee -a ${LOGFILE}

# print selected parameter settings
echo "STEPS TO BE PERFORMED/SKIPPED:" | tee -a ${LOGFILE}

if [ "$download_files" == "true" ]; then
    echo "required files will be downloaded" | tee -a ${LOGFILE}
elif [ "$download_files" == "false" ]; then
    echo "skipping download of required files" | tee -a ${LOGFILE}
fi

if [ "$install_env" == "true" ]; then
    echo "conda environment will be installed" | tee -a ${LOGFILE}
elif [ "$install_env" == "false" ]; then
    echo "skipping installation of conda environment" | tee -a ${LOGFILE}
fi

if [ "build_ref_genome" == "true" ]; then
    echo "building of reference genome will be performed" | tee -a ${LOGFILE}
    if [ "$download_ensembl_files" == "false" ]; then
        echo "skipping download of ensembl files needed for reference genome" | tee -a ${LOGFILE}
    fi
elif [ "build_ref_genome" == "false" ]; then
    echo "skipping building of reference genome" | tee -a ${LOGFILE}
fi

# print parameters to console and logfile
echo "Params:" | tee -a ${LOGFILE}
echo "cellranger version: $expectedcrv" | tee -a ${LOGFILE}
echo "Ensembl release: $ensrel" | tee -a ${LOGFILE}
echo "Genome release: $genomestring" | tee -a ${LOGFILE}
echo "Species: $species" | tee -a ${LOGFILE}
echo "nthreads: $nthreads" | tee -a ${LOGFILE}
echo "memgb: $memgb" | tee -a ${LOGFILE}
echo "work directory: $work_dir" | tee -a ${LOGFILE}
echo "conda env directory: $env_path" | tee -a ${LOGFILE}

# request parameter confirmation
read -r -p "Are these parameters correct? [y/N]: " response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
then
    echo "Parameters confirmed, continuing." | tee -a ${LOGFILE}
else
    echo "Aborting." | tee -a ${LOGFILE}
    exit 1
fi


# CREATE CONDA ENVIRONMENT, if install_env=true

path_to_env="{$env_path}cr3-velocyto-scanpy"

if [ "$install_env" = true ]; then
    if [ -z "$env_path" ]; then
        echo "Error: Please provide a path for the Conda environment installation"
        usage
        exit 1
    fi

    echo "Creating conda environment in ${path_to_env}..." | tee -a ${LOGFILE}
    echo "Start time: `date`" | tee -a ${LOGFILE}

    mamba env create --prefix path_to_env -c $work_dir/conda-bld -c conda-forge -c bioconda -y cellranger=3.1.0=0 scanpy=1.4.4.post1=py_3  velocyto.py=0.17.17=py36hc1659b7_0 samtoools=1.10=h9402c20_2 conda=4.8.2=py36_0 nextflow=19.10 java-jdk=0.0.112 | tee -a ${LOGFILE}
    echo "Completed installation of conda environment cr3-velocyto-scanpy at $path_to_env"
    echo "End time: `date`" | tee -a ${LOGFILE}
fi

}

initialize_scRNA_seq "$@"