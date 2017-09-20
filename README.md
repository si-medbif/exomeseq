# exomeseq
Analysing sequence data from whole exome sequencing or targeted gene panel sequencing.

The package consists of one python script that will

	1. Create a new folder based on the provided project name and parent folder.
	2. Prepare the configuration file with (dummy) links to required reference files.
	3. Copy/create any other files needed for the project into the target folder.

Instructions:

	1. Run the command (folder will be created):
	setup.py <full-path-to-project-folder>

	2. Add gene names to the file "genes.list" if there is a specific list of genes that should be included in the final report.
	3. Add sample names to the file "samples.list".
