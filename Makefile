CC = gcc

# Common headers and sources are given names for convenience
HDRS = hoalexer.h hoaparser.h simplehoa.h
SRCS = hoalexer.c hoaparser.c simplehoa.c

CFLAGS = -O3 -DNDEBUG
DBGFLAGS = -fsanitize=address -fno-omit-frame-pointer -g

.PHONY: clean

hoa2aig: $(SRCS) $(HDRS) hoa2aig.c aiger.c aiger.h
	$(CC) $(DBGFLAGS) -lm -o hoa2aig $(SRCS) aiger.c hoa2aig.c

# The parser is flex + bison based, everything is generated from
# hoa.l and hoa.y, the tokenizer and parser specifications
hoalexer.c: hoa.l hoaparser.c hoaparser.h
	flex --outfile=hoalexer.c --header-file=hoalexer.h hoa.l

hoaparser.c: hoa.y
	bison --defines --output=hoaparser.c hoa.y

parsertests: $(SRCS) $(HDRS) parsertests.c
	$(CC) $(DBGFLAGS) -o parsertests $(SRCS) parsertests.c
	ASAN_OPTIONS=detect_leaks=1
	cat examples/test1.ehoa | ./parsertests
	cat examples/test2.ehoa | ./parsertests
	cat examples/test3.ehoa | ./parsertests
	cat examples/aut1.ehoa | ./parsertests
	cat examples/aut2.ehoa | ./parsertests
	cat examples/aut3.ehoa | ./parsertests
	cat examples/aut4.ehoa | ./parsertests
	cat examples/aut5.ehoa | ./parsertests
	cat examples/aut6.ehoa | ./parsertests
	cat examples/aut7.ehoa | ./parsertests
	cat examples/aut8.ehoa | ./parsertests

clean:
	rm -f hoalexer.h hoalexer.c
	rm -f hoaparser.h hoaparser.c
	rm -f hoa2aig parsertests
