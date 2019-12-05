# Making the HOA2AIG file
HDRS = hoalexer.h hoaparser.h
SRCS = hoalexer.c hoaparser.c

CFLAGS = -O3 -DNDEBUG

.PHONY: tests clean parser

# The parser is flex + bison based, everything is generated from
# hoa.l and hoa.y, the tokenizer and parser specifications
parser: hoa.y hoa.l
	bison --defines --output=hoaparser.c hoa.y
	flex --outfile=hoalexer.c --header-file=hoalexer.h hoa.l
	cc -o hoaparser hoaparser.c hoalexer.c

tests: tests.c $(SRCS) $(HDRS)
	clang -g -fsanitize=address tests.c $(SRCS) -o tests
	ASAN_OPTIONS=detect_leaks=1
	./tests

clean:
	rm hoalexer.c
	rm hoaparser.h hoaparser.c
	rm tests
