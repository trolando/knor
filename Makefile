# Making the HOA2AIG file
HDRS = hoalexer.h hoaparser.h hoa.h
SRCS = hoalexer.c hoaparser.c

CFLAGS = -O3 -DNDEBUG

.PHONY: main tests clean

main: hoa2aig.c hoaparser.c hoalexer.c
	cc -o hoa2aig hoa2aig.c hoaparser.c hoalexer.c

# The parser is flex + bison based, everything is generated from
# hoa.l and hoa.y, the tokenizer and parser specifications
hoalexer.c: hoa.l hoa.y hoaparser.c hoaparser.h
	flex --outfile=hoalexer.c --header-file=hoalexer.h hoa.l

hoaparser.c: hoa.y
	bison --defines --output=hoaparser.c hoa.y

tests: tests.c $(SRCS) $(HDRS)
	clang -g -fsanitize=address tests.c $(SRCS) -o tests
	ASAN_OPTIONS=detect_leaks=1
	./tests

clean:
	rm hoalexer.h hoalexer.c
	rm hoaparser.h hoaparser.c
	rm hoa2aig tests
