#Sean M. Law

include ../makefile.defs

OBJS = libmoletools.$(LIBEXT)

all:	$(BINDIR) libmoletools.$(LIBEXT)  larmorca

clean:
	rm -rf *.o
	rm -rf *stackdump
	rm -rf ../bin/*
	rm -rf *~

$(BINDIR):
	mkdir -p $(BINDIR)

libmoletools.$(LIBEXT):
	$(MAKE) -C ../lib libmoletools.$(LIBEXT)

larmorca: larmorca.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LIBFLAGS) $(INCFLAGS) -o $(BINDIR)/$@ $< $(LINKFLAGS)
	
%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(INCFLAGS) -c $< 
