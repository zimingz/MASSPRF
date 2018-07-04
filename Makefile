CC=g++
CFLAGS = -std=c++0x -O3 -pthread -g

massprf: MASSprf.cpp PRFCluster.cpp base.cpp kfunc.cpp
		$(CC) $(CFLAGS) $^ -o bin/$@
