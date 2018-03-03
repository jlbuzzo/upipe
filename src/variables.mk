PIPELINE_NAME := RETRONATOR

SUFFIX := .bam

vpath %.pl /home/projects2/retrocopies/teste_leonel/library/src
SRC := /home/projects2/retrocopies/teste_leonel/library/src

INPUT_DIR ?= ../tst
OUTPUT_DIR ?= ./save
TEMP_PROCESS_DIR ?= $(OUTPUT_DIR)/temp/
DUMP_DIR ?= $(OUTPUT_DIR)/result/dump/


SORTED := TRUE
REF_ANNOTATION := /home/projects/databases/gencode/release26/gencode.v26.annotation.gtf
REP_ANNOTATION := /home/projects/local_collab/gh_RTC_simulation/make/library/annotation/rep.hg38.converted.bed
REP_ANNOTATION_manual := /home/projects/local_collab/gh_RTC_simulation/make/library/annotation/rep.hg38.converted.manual.bed
DISTANCE := 750000
MIN_DIST := 1000000

