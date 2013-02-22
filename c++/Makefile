.SUFFIXES:
.SUFFIXES: .c .o

CXX = g++
CC = $(CXX)

BOOST ?= $(soft_msi)/boost/1.47.0
LIBSEQ ?= $(soft_msi)/libsequence/1.7.3
SAMTOOLS = $(prefix)/soft/samtools
BAMTOOLS = $(prefix)/soft/bamtools
GECO = $(prefix)/soft/geco

DEBUG_LEVEL     = -g
EXTRA_CCFLAGS   = -Wall -m64
CXXFLAGS        = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS         = $(CXXFLAGS)

CPPFLAGS        = -I$(BOOST)/include -I$(LIBSEQ)/include -I$(SAMTOOLS)/include -I$(BAMTOOLS)/include -I$(GECO)/include
LDLIBS          = -lsequence -lboost_system -lboost_regex -lboost_iostreams -lboost_filesystem -lboost_program_options -lbamtools -lgeco
LDFLAGS         = -L$(BOOST)/lib -L$(LIBSEQ)/lib -L$(SAMTOOLS)/lib -L$(BAMTOOLS)/lib -L$(GECO)/lib

SRCEXT   = cc
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

BINNAMES = test msConvert readGz cov_window bamDepth blastTiling \
  sspStat sspDiff sspSetOG sspFilter sspSample sspSlide sspSelectInds sspConvert ssp2Vcf\
  bamPreprocess bamCheck viewRead \
  bamPickStretched bamPickOrphan bamPickReads bamPickBpReads \
  bam2Fastq bamISD bamStat \
  seqSplit alnStat \
#  seqRec
BINS    := $(addprefix $(BINDIR)/, $(BINNAMES))
SRCS    := $(wildcard $(SRCDIR)/*.$(SRCEXT))
OBJS    := $(patsubst $(SRCDIR)/%.$(SRCEXT), $(OBJDIR)/%.o, $(SRCS))


.PHONY: all clean

all: $(BINS)

$(OBJS): | $(OBJDIR)
$(OBJS): $(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(SRCDIR)/%.h

$(BINS): | $(BINDIR)
$(BINS): $(BINDIR)/%: $(SRCDIR)/%.$(SRCEXT)

BIN_SSP = sspFilter sspSelectInds sspSetOG sspSample sspConvert ssp2Vcf seqSplit
$(addprefix $(BINDIR)/, $(BIN_SSP)): $(OBJDIR)/ssp.o

BIN_BAM = bamCheck viewRead bamPickStretched bamPickOrphan bamPickReads bamPickBpReads bam2Fastq bamISD 
$(addprefix $(BINDIR)/, $(BIN_BAM)): $(OBJDIR)/bam_utils.o

$(BINDIR)/alnStat: $(OBJDIR)/read.o
$(BINDIR)/blastTiling: $(OBJDIR)/common.o
$(BINDIR)/bamDepth: $(SAMTOOLS)/lib/libbam.a -lbam -lm -lz
$(BINDIR)/cov_window: $(addprefix $(OBJDIR)/, read.o ssp.o)
$(BINDIR)/seqRec: $(addprefix $(OBJDIR)/, geco_sequences.o read.o ssp.o)
$(BINDIR)/bamPreprocess: $(addprefix $(OBJDIR)/, bam_utils.o read.o)
$(BINDIR)/bamStat: $(addprefix $(OBJDIR)/, bam_utils.o read.o)

$(BINDIR)/%:
	@echo "Linking $@..."
	@$(CC) $(CCFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $^ -o $@

$(OBJDIR)/%.o:
	@echo "Compiling $@..."
	@$(CC) $(CCFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) $< -c -o $@

$(BINDIR) $(OBJDIR):
	@mkdir -p $@

clean:
	$(RM) -r $(OBJDIR)
	$(RM) -r $(BINDIR)
