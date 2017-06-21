CXXFLAGS=`pkg-config --cflags meep` `root-config --cflags` -O2 -g
LDFLAGS=`pkg-config --libs meep` `root-config --libs`


arthern_poc: arthern_poc.cc
	g++ $(CXXFLAGS) arthern_poc.cc -o arthern_poc $(LDFLAGS) 
