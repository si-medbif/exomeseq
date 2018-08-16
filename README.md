# Exomeseq
Analysing sequence data from whole exome sequencing or targeted gene panel sequencing.

## Requirements
    * Ubuntu, or any other compatible OS
    * Python 3
    * Docker
    * Git
    
## Setup

The package consists of python scripts that will

	1. Create all required folders and configuration files.
	2. Move any provided FASTQ files to their appropriate location.
	3. Prepare the configuration file with links to required reference files.
	4. Create and run scripts for analysing each sample. 

###Instructions:

	1. Create a new project folder and move to this location
	$ mkdir projectA
	$ cd projectA
	 
    2. Copy all fastq files to the current folder. The name of the fastq files have to follow the following pattern: <sample_name>.1.fastq.gz
    $ cp <fastq_file1/2> <sample_name>.1.fastq.gz
    $ cp <fastq_file2/2> <sample_name>.2.fastq.gz
    
    3. Download the pipeline code from GitHub:
    $ git clone https://github.com/si-medbif/exomeseq.git
    
Fully automatic, no trimming of reads 
      
	4. Start the command to setup and run the pipeline 
	$ exomeseq/runme.sh
	
Half-automatic, user decides to trim reads or not

	4. Start the command to setup the pipeline 
	$ exomeseq/runsetup.sh
	
	5a. Run the pipeline without trimming reads  
    $ exomeseq/run_analysis.sh
    
    5b. Run the pipeline after trimming reads
    $ exomeseq/run_qc_analysis.sh
	 
The combined report will be in the current folder when finished.

    $ ls -l full_report.txt
    