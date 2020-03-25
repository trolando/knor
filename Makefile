CC = gcc
HDRS = hoalexer.h hoaparser.h simplehoa.h
SRCS = hoalexer.c hoaparser.c simplehoa.c hoa2aig.c

CFLAGS = -O3 -DNDEBUG
DBGFLAGS = -fsanitize=address -fno-omit-frame-pointer -g

.PHONY: clean

hoa2aig: $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) -o hoa2aig $(SRCS)

# The parser is flex + bison based, everything is generated from
# hoa.l and hoa.y, the tokenizer and parser specifications
hoalexer.c: hoa.l hoaparser.c hoaparser.h
	flex --outfile=hoalexer.c --header-file=hoalexer.h hoa.l

hoaparser.c: hoa.y
	bison --defines --output=hoaparser.c hoa.y

tests: $(SRCS) $(HDRS)
	$(CC) $(DBGFLAGS) -o tests $(SRCS)
	ASAN_OPTIONS=detect_leaks=1
	cat examples/test1.ehoa | ./tests
	cat examples/test2.ehoa | ./tests
	cat examples/test3.ehoa | ./tests
	cat examples/aut1.ehoa | ./tests
	cat examples/aut2.ehoa | ./tests
	cat examples/aut3.ehoa | ./tests
	cat examples/aut4.ehoa | ./tests
	cat examples/aut5.ehoa | ./tests
	cat examples/aut6.ehoa | ./tests
	cat examples/aut7.ehoa | ./tests
	cat examples/aut8.ehoa | ./tests

clean:
	rm -f hoalexer.h hoalexer.c
	rm -f hoaparser.h hoaparser.c
	rm -f hoa2aig tests
