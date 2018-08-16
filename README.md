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

Instructions:

	1. Create a new project folder and move to this location
	$ mkdir projectA
	$ cd projectA
	 
    2. Copy all fastq files to the current folder. The name of the fastq files have to follow the following pattern: <sample_name>.1.fastq.gz
    $ cp <fastq_file1/2> <sample_name>.1.fastq.gz
    $ cp <fastq_file2/2> <sample_name>.2.fastq.gz
    
    3. Download the pipeline code from GitHub:
    $ git clone https://github.com/si-medbif/exomeseq.git
      
	4. Start the command to setup and run the pipeline 
	$ exomeseq/runme.sh
	 
    5. The combined report will be in the current folder when finished.
    