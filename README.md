# Exomeseq
Analysing sequence data from whole exome sequencing or targeted gene panel sequencing.

## Requirements

### Installed software
    * Ubuntu, or other compatible OS (tested with ubuntu 16.04)
    * Python 3 (should be installed by default)
    * Docker
    * Git
    * Unzip
    * Bgzip / tabix
    These can be installed by running the following commands:
    $ apt update
    $ apt install git unzip tabix docker.io python-all
    
### Hardware requirements
    * Storage: 256 GB
    * Memory (RAM): 8 GB
    * Processors (cores): 4
    
## Setup

The package consists of several python scripts that will

	1. Create all required folders and configuration files.
	2. Move any provided FASTQ files to their appropriate location.
	3. Prepare the configuration file with links to required reference files.
	4. Create and run scripts for analysing each sample. 

## Instructions:

    0. It is recommended to use a service like 'screen' or 'tmux' to detach the process. The terminal can then be closed without interrupting the program.
	1. Create a new project folder and move to this location
	$ mkdir projectA
	$ cd projectA
	 
    2. Copy all fastq files to the current folder. The name of the fastq files have to follow the following pattern: <sample_name>.1.fastq.gz
    $ cp <fastq_file1/2> <sample_name>.1.fastq.gz
    $ cp <fastq_file2/2> <sample_name>.2.fastq.gz
    
    3. Copy the bed file containing the target regions into the current folder.
    $ cp <target_regions>.bed .
    
    4. Download the pipeline code from GitHub:
    $ git clone https://github.com/si-medbif/exomeseq.git

#### The user can at this step choose one of two paths:
    
##### Fully automatic, no trimming of reads 
      
	5. Start the command to setup and run the pipeline 
	$ exomeseq/runme.sh
	
##### Half-automatic, user decides to trim reads or not

	5. Start the command to setup the pipeline 
	$ exomeseq/runsetup.sh
	
	6. Check the FastQC output (in "sample/FastQC_pre/") and then do either 6a or 6b:
	
	6a. Run the pipeline without trimming reads  
    $ exomeseq/run_analysis.sh
    
    6b. Run the pipeline after trimming reads
    $ exomeseq/run_qc_analysis.sh
	 
The combined report will be in the current folder when finished.

    $ ls -l full_report.txt
