FLAGS = -Wall -Wextra -g
DEPENDENCIES = fetch_seq.h align_seq.h

all: align

align: align.o fetch_seq.o align_seq.o
	gcc ${FLAGS} -o $@ $^

align.o: align.c fetch_seq.h align_seq.h
	gcc ${FLAGS} -c $<

fetch_seq.o: fetch_seq.c fetch_seq.h
	gcc ${FLAGS} -c $<

align_seq.o: align_seq.c align_seq.h
	gcc ${FLAGS} -c $<

clean:
	rm -f *.o align *.gch

