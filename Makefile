HDRS = hoalexer.h hoaparser.h simplehoa.h
SRCS = hoalexer.c hoaparser.c simplehoa.c hoa2aig.c

CFLAGS = -O3 -DNDEBUG
# With clang you can use the following instead
DBGFLAGS = -g -fsanitize=address

.PHONY: tests clean

hoa2aig: $(SRCS) $(HDRS)
	cc $(DBGFLAGS) -o hoa2aig $(SRCS)

# The parser is flex + bison based, everything is generated from
# hoa.l and hoa.y, the tokenizer and parser specifications
hoalexer.c: hoa.l hoaparser.c hoaparser.h
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
