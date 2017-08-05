# program name and object files
PROGRAM = test
OBJS = test.o

# define
CC = g++
#CFLAGS = -O2 -Wall
#LDFLAGS = -lm -lstdc++ -lpthread

CFLAGS=-c -Wall -std=c++0x -g3 -Ofast -msse2  -I. -I/usr/local/cuda/include -I/home/argon/caffe/include/ -I/home/argon/caffe/src/
#LDFLAGS= -L/usr/local/lib -L/usr/local/cuda/lib64 -L/usr/local/cuda/lib -lcudart -lcublas -lcurand-lglog -lgflags -lprotobuf -lleveldb -lsnappy -llmdb -lboost_system -lhdf5_hl -lhdf5 -lm -lopencv_core -lopencv_highgui -lopencv_imgproc -lboost_thread -lstdc++ -lcudnn -lcblas -latlas -L/home/argon/caffe/build/lib/ -lcaffe -lproto
LDFLAGS= -L/usr/local/lib -L/usr/local/cuda/lib64 -L/usr/local/cuda/lib -lcudart -lcublas -lcurand -lglog -lgflags -lprotobuf -lleveldb -lsnappy -llmdb -lboost_system -lhdf5_hl -lhdf5 -lm -lopencv_core -lopencv_highgui -lopencv_imgproc -lboost_thread -lstdc++ -lcblas -latlas -L/home/argon/caffe/build/lib/ -lcaffe

# target	 '$^' ... list of files.
$(PROGRAM): $(OBJS)
	$(CC) -o $(PROGRAM) $^ $(LDFLAGS)

# suffixe rule   '$<' ... top file name of list of files.
.cpp.o:
	$(CC) $(CFLAGS) -c $<

# delete target
.PHONY: clean
clean:
	$(RM) $(PROGRAM) $(OBJS)
