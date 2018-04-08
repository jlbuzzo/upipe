FROM debian:stretch-slim

ENV SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 \
		SAMTOOLS_PATH=/samtools-1.7

WORKDIR /

RUN apt-get update && \
		apt-get install -y \
			perl \
			parallel \
			pigz \
			wget \
			gcc \
			make \
			libncurses5-dev \
			autoconf \
			automake \
			libtool \
			zlib1g-dev \
			libbz2-dev \
			liblzma-dev \
			bzip2 \
			gzip && \
		apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget $SAMTOOLS_URL && tar xjf samtools-1.7.tar.bz2 && (cd $SAMTOOLS_PATH && make)

RUN cpan App::cpanminus
RUN cpanm Statistics::Basic Number::Format


COPY . /upipe

WORKDIR /upipe

ENV PATH=$SAMTOOLS_PATH:$PATH

CMD ["bash"]
