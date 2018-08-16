# Exomeseq
Analysing sequence data from whole exome sequencing or targeted gene panel sequencing.

The package consists of python scripts that will

	1. Create a new folder based on the provided project name and parent folder.
	2. Prepare the configuration file with (dummy) links to required reference files.
	3. Copy/create any other files needed for the project into the target folder.

Instructions:

	1. Run setup.py. If the project folder already exists, setup.py will exit with an error message.
	$ ./setup.py <full-path-to-project-folder>  
	 
    1.1 Navigate to the newly created project folder
    $ cd <full-path-to-project-folder>  
     
    2. Place the fastq files in the fastq folder. The name of the fastq files have to follow the following pattern: <sample_name>.1.fastq.gz
    $ cp <fastq_file1/2> fastq/<sample_name>.1.fastq.gz
    $ cp <fastq_file2/2> fastq/<sample_name>.2.fastq.gz  
      
	2. [Optional] Add gene names to the file "genes.list" if there is a specific list of genes that should be included in the final report.
	  
	3. Add sample names (i.e. <sample_name>) to the file "samples.list".
	 
	4. Open the exome.cfg file and change the <CHANGEME> tag so the paths to the reference files are correct.
	   
	5. Run make_scripts.py. Use '-e' to specify if this is an exome data set.
	$ ./make_scripts.py -e <sample_name>  
	 
    6. Start analysis with:
    $ ./<sample_name>_GATK.sh
     
    7. Start annotation with:
    $ ./<sample_name>_annotation.sh
    
    8. Produce summarized report files with:
    $ ./parsevcf.py
    