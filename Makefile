CC=g++
CFLAGS = -std=c++0x -O3 -pthread -g

massprf: MASSprf.cpp PRFCluster.cpp base.cpp
		$(CC) $(CFLAGS) $^ -o bin/$@
