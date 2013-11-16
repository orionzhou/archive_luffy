.SUFFIXES:
.SUFFIXES: .c .o

CXX = g++
CC = $(CXX)

BOOST    ?= $(soft)/boost
LIBSEQ   ?= $(soft)/libsequence
SAMTOOLS ?= $(soft)/samtools
BAMTOOLS ?= $(soft)/bamtools
GECO     ?= $(soft)/geco

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

BINNAMES = test msConvert readGz blastTiling galTiling \
  sspSum sspStat sspDiff sspSetOG sspFilter sspSample sspSlide sspSubset sspConvert ssp2Vcf\
  bamViewRead bamPreprocess bamCheck  \
  bamPickStretched bamPickOrphan bamPickBpReads \
  bam2Fastq bamStat bamOrphanStat bamDeDup \
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

BIN_SSP = sspSum sspStat sspDiff sspFilter sspSubset sspSetOG sspSample sspConvert ssp2Vcf seqSplit
$(addprefix $(BINDIR)/, $(BIN_SSP)): $(OBJDIR)/ssp.o

BIN_BAM = bamCheck bamViewRead bamPickStretched bamPickOrphan bamPickBpReads bam2Fastq bamStat bamDeDup
$(addprefix $(BINDIR)/, $(BIN_BAM)): $(OBJDIR)/bam_utils.o

$(BINDIR)/alnStat: $(OBJDIR)/read.o
$(BINDIR)/blastTiling: $(OBJDIR)/location.o
$(BINDIR)/galTiling: $(OBJDIR)/location.o
#$(BINDIR)/bamDepth: $(SAMTOOLS)/lib/libbam.a -lbam -lm -lz
$(BINDIR)/cov_window: $(addprefix $(OBJDIR)/, read.o ssp.o)
$(BINDIR)/seqRec: $(addprefix $(OBJDIR)/, geco_sequences.o read.o ssp.o)
$(BINDIR)/bamPreprocess: $(addprefix $(OBJDIR)/, bam_utils.o read.o)
$(BINDIR)/bamOrphanStat: $(addprefix $(OBJDIR)/, bam_utils.o read.o)

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
