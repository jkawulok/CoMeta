all: cometa class2best genlist num_class tsk seq_gi2tax

BOOST_LIB = /software/local/libs/Boost/1.52.0-intel-13.0.1/lib
BOOST_H = /software/local/libs/Boost/1.52.0-intel-13.0.1/include

COMETA_ROOT_DIR = bin
COMETA_MAIN_DIR = src/_cometa
COMETA_seqGI2tax_DIR = src/_seq_gi2tax


CC     = g++
CFLAGS    =   -O3 -m64 -fopenmp  -I $(BOOST_H) -L $(BOOST_LIB)
CLINK    = -fopenmp -O3 

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@



cometa: $(COMETA_MAIN_DIR)/stdafx.o $(COMETA_MAIN_DIR)/timer.o $(COMETA_MAIN_DIR)/radix.o $(COMETA_MAIN_DIR)/fastq_reader.o \
		$(COMETA_MAIN_DIR)/kmer_reader.o $(COMETA_MAIN_DIR)/classifier.o $(COMETA_MAIN_DIR)/cometa.o $(COMETA_MAIN_DIR)/kmclass.o \
		$(COMETA_MAIN_DIR)/kmer_bin.o $(COMETA_MAIN_DIR)/meta_class.o $(COMETA_MAIN_DIR)/splitter.o $(COMETA_MAIN_DIR)/splitter_kmer.o \
		$(COMETA_MAIN_DIR)/splitter2class.o $(COMETA_MAIN_DIR)/wrapers.o $(COMETA_MAIN_DIR)/Tr_wrapers.o 
	$(CC) $(CLINK) -o $(COMETA_ROOT_DIR)/$@ $(COMETA_MAIN_DIR)/stdafx.o $(COMETA_MAIN_DIR)/timer.o $(COMETA_MAIN_DIR)/radix.o \
	$(COMETA_MAIN_DIR)/fastq_reader.o $(COMETA_MAIN_DIR)/kmer_reader.o $(COMETA_MAIN_DIR)/classifier.o $(COMETA_MAIN_DIR)/cometa.o \
	$(COMETA_MAIN_DIR)/kmclass.o $(COMETA_MAIN_DIR)/kmer_bin.o $(COMETA_MAIN_DIR)/meta_class.o $(COMETA_MAIN_DIR)/splitter.o \
	$(COMETA_MAIN_DIR)/splitter_kmer.o $(COMETA_MAIN_DIR)/splitter2class.o $(COMETA_MAIN_DIR)/wrapers.o $(COMETA_MAIN_DIR)/Tr_wrapers.o \
	$(COMETA_MAIN_DIR)/alibelf64.a \
	-L $(BOOST_LIB) -I $(BOOST_H) -lboost_system -lboost_filesystem -lboost_thread 



class2best: $(COMETA_MAIN_DIR)/class2best.o
	$(CC) $(CLINK) -o $(COMETA_ROOT_DIR)/$@ $(COMETA_MAIN_DIR)/class2best.o $(COMETA_MAIN_DIR)/stdafx.o  \
	-L $(BOOST_LIB) -I $(BOOST_H) -lboost_system -lboost_filesystem -lboost_thread


genlist: $(COMETA_MAIN_DIR)/gen_list.o
	$(CC) $(CLINK) -o $(COMETA_ROOT_DIR)/$@ $(COMETA_MAIN_DIR)/gen_list.o $(COMETA_MAIN_DIR)/stdafx.o  \
	-L $(BOOST_LIB) -I $(BOOST_H) -lboost_system -lboost_filesystem -lboost_thread


num_class: $(COMETA_MAIN_DIR)/num_class.o
	$(CC) $(CLINK) -o $(COMETA_ROOT_DIR)/$@ $(COMETA_MAIN_DIR)/num_class.o $(COMETA_MAIN_DIR)/stdafx.o  \
	-L $(BOOST_LIB) -I $(BOOST_H) -lboost_system -lboost_filesystem -lboost_thread



tsk:	$(COMETA_MAIN_DIR)/stdafx.o $(COMETA_MAIN_DIR)/kmer_reader.o $(COMETA_MAIN_DIR)/splitter_kmer.o $(COMETA_MAIN_DIR)/Tr_KMclass.o \
		$(COMETA_MAIN_DIR)/Tr_wrapers.o $(COMETA_MAIN_DIR)/TSK.o 
	$(CC) $(CLINK) -o $(COMETA_ROOT_DIR)/$@ $(COMETA_MAIN_DIR)/stdafx.o $(COMETA_MAIN_DIR)/kmer_reader.o $(COMETA_MAIN_DIR)/splitter_kmer.o \
	$(COMETA_MAIN_DIR)/Tr_KMclass.o $(COMETA_MAIN_DIR)/Tr_wrapers.o $(COMETA_MAIN_DIR)/TSK.o $(COMETA_MAIN_DIR)/alibelf64.a  \
	-L $(BOOST_LIB) -I $(BOOST_H) -lboost_system -lboost_filesystem -lboost_thread

	
seq_gi2tax:			$(COMETA_seqGI2tax_DIR)/seqGI2tax.o
	$(CC) $(CLINK) -o $(COMETA_ROOT_DIR)/$@ $(COMETA_seqGI2tax_DIR)/seqGI2tax.o \
	-L $(BOOST_LIB) -I $(BOOST_H) -lboost_system -lboost_filesystem -lboost_thread
	
	
clean:
	-rm $(COMETA_MAIN_DIR)/*.o
	-rm $(COMETA_seqGI2tax_DIR)/*.o


all: cometa class2best genlist num_class tsk seq_gi2tax



