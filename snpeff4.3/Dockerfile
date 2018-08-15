#################################################################
# Dockerfile
#
# Version:          1
# Software:         snpEff & SnpSift
# Software Version: 4.3
# Description:      Genomic variant annotations and functional effect prediction toolbox.
# Website:          http://snpeff.sourceforge.net/
# Tags:             Genomics
# Provides:         snpeff for GRCh38.86 and SnpSift
# Base Image:       biocontainers/biocontainers
# Build Cmd:        docker build --rm -t snpeff38:v1 snpeff4.3/.
# Run Cmd1:         docker run --rm -it snpeff38:v1 snpeff <options> <files>
# Run Cmd2:         docker run --rm -it snpeff38:v1 snpsift <options> <files>
#################################################################

# Set the base image to Ubuntu
FROM biocontainers/biocontainers

################## BEGIN INSTALLATION ######################

RUN cd && wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip && cd snpEff && java -jar snpEff.jar download GRCh38.86 && \
    cd ../ && rm snpEff_latest_core.zip && \
    bash -c 'echo -e "#!/bin/bash\njava -jar -Xmx4G /home/biodocker/snpEff/snpEff.jar \$@" > /home/biodocker/bin/snpeff' && chmod +x /home/biodocker/bin/snpeff && \
    bash -c 'echo -e "#!/bin/bash\njava -jar -Xmx4G /home/biodocker/snpEff/SnpSift.jar \$@" > /home/biodocker/bin/snpsift' && chmod +x /home/biodocker/bin/snpsift



ENV PATH=$PATH:/home/biocontainers/snpEff/scripts:/home/biocontainers/snpEff/scripts/gsa

# CHANGE WORKDIR TO /DATA
WORKDIR /data

# DEFINE DEFAULT COMMAND
CMD ["snpeff"]

##################### INSTALLATION END #####################
