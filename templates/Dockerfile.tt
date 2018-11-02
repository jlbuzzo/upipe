############################## HEADER ##########################################
#
# Dockerfile to make containers of Conda environments.
#
#
# Some information here...
################################################################################





################################################################################

# Inherit base image: debian:latest.
FROM continuumio/miniconda:latest

# Some metadata
MAINTAINER Jose Leonel L. Buzzo

WORKDIR /

# The default user at container's startup: lion.
RUN useradd -d /home/lion -ms /bin/bash -u 1541 -g root -p lion lion

USER root


# System wide installations. Already has curl, wget, bzip2, ca-certificates
# and git.
RUN apt-get update && \
	apt-get -y upgrade && \
	apt-get install -y \
		gcc \
		g++ \
		gfortran \
		make \
		parallel \
		autoconf \
		automake \
		libtool \
		libtbb-dev \
		libncurses5-dev \
		libbz2-dev \
		liblzma-dev \
		zlib1g-dev \
		pigz \
		tar \
		unzip \
		gnupg \
		gzip && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Adentrate the base work directory.
WORKDIR /home/lion

# Change to default user: lion.
#USER lion



############################## PROCEDURAL INSTALLATIONS ########################

# NOTE: These are the RECOMMENDED installation instructions!


# Create a Conda environment: ponga.
RUN conda create --name ponga
#	--override-channels -c iuc \
#	-c bioconda -c conda-forge \
#	-c defaults -c r python=2.7.13 \
#	bioconductor-deseq2=1.16.1 \
#	r-base=3.4.1 r-pheatmap \
#	bioconductor-vsn bioconductor-biocparallel=1.12.0 \
#	r-ggrepel \
#	star=2.5.3a \
#	bedtools=2.25.0 \
#	samtools=1.1 \
#	stringtie=1.3.3 \
#	igvtools=2.3.93 \
#	ucsc-genepredtobed ucsc-gtftogenepred ucsc-genepredtogtf ucsc-bedgraphtobigwig r-hexbin


# Another way to build ponga. Uncomment only if conda create fails.
#RUN ponga Build -s all


# Change to default user: lion.
USER lion

# Adjust user lion permissions and squire environment at startup.
ENV PATH=/opt/conda/etc/profile.d/:$PATH
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "source activate ponga" >> ~/.bashrc
RUN conda config --describe > ~/.condarc


# Another way to activate squire environment.
#RUN source activate ponga
#RUN conda config --set default_env ponga



############################## MAKEFILE PART ###################################

# This is a token for the Makefile to recognize container's environment.
RUN touch .ponga_switch



# Container's initial command.
CMD ["/bin/bash"]
