SRCS = 
HDRS = 

CFLAGS = -O3 -DNDEBUG

.PHONY: tests clean all

all: hoa.y hoa.l
	yacc --defines --output=hoaparser.c hoa.y
	lex --outfile=hoalexer.c hoa.l
	cc -o hoaparser hoaparser.c hoalexer.c

tests: tests.c $(SRCS) $(HDRS)
	cc -g -fsanitize=address tests.c $(SRCS) -o tests
	ASAN_OPTIONS=detect_leaks=1
	./tests

clean:
	rm hoalexer.c
	rm hoaparser.h hoaparser.c
	rm tests
