BOOST = $(soft)/boost/1.44.0
LIBSEQ = $(soft)/libsequence/1.6.9
SAMTOOLS = $(prefix)/soft/samtools
BAMTOOLS = $(prefix)/soft/bamtools
GECO = $(prefix)/soft/geco

CC      = g++
CFLAGS  = -Wall -g -I$(BOOST)/include -I$(LIBSEQ)/include -I$(SAMTOOLS)/include -I$(BAMTOOLS)/include -I$(GECO)/include
LIBS    = -lbamtools -lsequence -lboost_system -lboost_regex -lboost_iostreams -lboost_filesystem -lboost_program_options -lgeco
LDFLAGS = -L$(BOOST)/lib -L$(LIBSEQ)/lib -L$(SAMTOOLS)/lib -L$(BAMTOOLS)/lib -L$(GECO)/lib
DEP     =
BIN_DIR = bin
BIN_NAME= test msConvert readGz cov_window bamDepth \
  sspStat sspDiff sspSetOG sspFilter sspSample sspSlide sspSelectInds sspConvert ssp2Vcf\
  bamPreprocess bamCheck viewRead \
  bamPickStretched bamPickOrphan bamPickReads bamPickBpReads \
  bam2Fastq bamISD bamStat \
  seqRec seqSplit alnStat
EXE     = $(addprefix $(BIN_DIR)/,$(BIN_NAME))
vpath %.h src 
vpath %.cc src
vpath %.o obj
all: $(BIN_DIR) $(EXE) 

bin/msConvert: msConvert.cc
	$(CC) $(CFLAGS) src/msConvert.cc $(LDFLAGS) $(LIBS) -o $@

bin/test: test.cc
	$(CC) $(CFLAGS) src/test.cc $(LDFLAGS) $(LIBS) -o $@

bin/readGz: readGz.cc
	$(CC) $(CFLAGS) src/readGz.cc $(LDFLAGS) $(LIBS) -o $@

bin/sspStat: sspStat.cc
	$(CC) $(CFLAGS) src/sspStat.cc $(LDFLAGS) $(LIBS) -o $@

bin/sspDiff: sspDiff.cc
	$(CC) $(CFLAGS) src/sspDiff.cc $(LDFLAGS) $(LIBS) -o $@

bin/sspSlide: sspSlide.cc
	$(CC) $(CFLAGS) src/sspSlide.cc $(LDFLAGS) $(LIBS) -o $@

bin/bamDepth: bamDepth.cc $(SAMTOOLS)/lib/libbam.a
	$(CC) $(CFLAGS) src/bamDepth.cc $(LDFLAGS) $(LIBS) -lbam -lm -lz -o $@

bin/sspFilter: sspFilter.cc ssp.o
	$(CC) $(CFLAGS) src/sspFilter.cc obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/sspSetOG: sspSetOG.cc obj/ssp.o
	$(CC) $(CFLAGS) src/sspSetOG.cc obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/sspSample: sspSample.cc ssp.o
	$(CC) $(CFLAGS) src/sspSample.cc obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/sspSelectInds: sspSelectInds.cc ssp.o
	$(CC) $(CFLAGS) src/sspSelectInds.cc obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/sspConvert: sspConvert.cc ssp.o
	$(CC) $(CFLAGS) src/sspConvert.cc obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/ssp2Vcf: ssp2Vcf.cc ssp.o
	$(CC) $(CFLAGS) src/ssp2Vcf.cc obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/cov_window: cov_window.cc read.o ssp.o
	$(CC) $(CFLAGS) src/cov_window.cc obj/read.o obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/seqRec: seqRec.h seqRec.cc geco_sequences.o read.o ssp.o
	$(CC) $(CFLAGS) src/seqRec.cc obj/geco_sequences.o obj/read.o obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/seqSplit: seqSplit.cc ssp.o
	$(CC) $(CFLAGS) src/seqSplit.cc obj/ssp.o $(LDFLAGS) $(LIBS) -o $@

bin/alnStat: alnStat.cc read.o
	$(CC) $(CFLAGS) src/alnStat.cc obj/read.o $(LDFLAGS) $(LIBS) -o $@

#high-level
bin/bamCheck: bamCheck.cc bam_utils.o
	$(CC) $(CFLAGS) src/bamCheck.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/viewRead: viewRead.cc bam_utils.o
	$(CC) $(CFLAGS) src/viewRead.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/bamPickStretched: bamPickStretched.cc bam_utils.o
	$(CC) $(CFLAGS) src/bamPickStretched.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/bamPickOrphan: bamPickOrphan.cc bam_utils.o
	$(CC) $(CFLAGS) src/bamPickOrphan.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/bamPickReads: bamPickReads.cc bam_utils.o
	$(CC) $(CFLAGS) src/bamPickReads.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/bamPickBpReads: bamPickBpReads.cc bam_utils.o
	$(CC) $(CFLAGS) src/bamPickBpReads.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/bamPreprocess: bamPreprocess.cc bam_utils.o read.o
	$(CC) $(CFLAGS) src/bamPreprocess.cc obj/bam_utils.o obj/read.o $(LDFLAGS) $(LIBS) -o $@

bin/bam2Fastq: bam2Fastq.cc bam_utils.o
	$(CC) $(CFLAGS) src/bam2Fastq.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/bamISD: bamISD.cc bam_utils.o
	$(CC) $(CFLAGS) src/bamISD.cc obj/bam_utils.o $(LDFLAGS) $(LIBS) -o $@

bin/bamStat: bamStat.cc bam_utils.o read.o
	$(CC) $(CFLAGS) src/bamStat.cc obj/bam_utils.o obj/read.o $(LDFLAGS) $(LIBS) -o $@

#low-level
obj/bam_utils.o: bam_utils.cc bam_utils.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $< -c -o $@

obj/geco_sequences.o: geco_sequences.cc geco_sequences.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $< -c -o $@

obj/ssp.o: ssp.cc ssp.h 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $< -c -o $@

obj/read.o: read.cc read.h 
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $< -c -o $@

obj/common.o: common.cc common.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $< -c -o $@

obj/init.o: init.cc init.h
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) $< -c -o $@

$(BIN_DIR):
	@echo Bin Dir is $(BIN_DIR)
	@if [ ! -d $(BIN_DIR) ]; then mkdir -p $(BIN_DIR); fi

clean:
	rm -f obj/*.o
	rm -f bin/* 

