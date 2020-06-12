CC = clang

# Common headers and sources are given names for convenience
HDRS = hoalexer.h hoaparser.h simplehoa.h
SRCS = hoalexer.c hoaparser.c simplehoa.c

CFLAGS = -O3 -DNDEBUG
DBGFLAGS = -fsanitize=address -fno-omit-frame-pointer -g

hoa2aig: $(SRCS) $(HDRS) hoa2aig.c aiger/aiger.c aiger/aiger.h
	$(CC) $(CFLAGS) -lm -o hoa2aig $(SRCS) aiger/aiger.c hoa2aig.c

hoa2pg: $(SRCS) $(HDRS) hoa2pg.c
	$(CC) $(CFLAGS) -o hoa2pg $(SRCS) hoa2pg.c

hoacheck: $(SRCS) $(HDRS) hoacheck.c
	$(CC) $(CFLAGS) -o hoacheck $(SRCS) hoacheck.c

# The parser is flex + bison based, everything is generated from
# hoa.l and hoa.y, the tokenizer and parser specifications
hoalexer.c: hoa.l hoaparser.c hoaparser.h
	flex --outfile=hoalexer.c --header-file=hoalexer.h hoa.l

hoaparser.c: hoa.y
	bison --defines --output=hoaparser.c hoa.y

.PHONY: clean all

clean:
	rm -f hoalexer.h hoalexer.c
	rm -f hoaparser.h hoaparser.c
	rm -f hoa2aig hoacheck hoa2pg

all: hoa2aig hoa2pg
