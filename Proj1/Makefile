PAR=spar
SEQ=sseq
TARGET=spar sseq

all: studentsseq.c studentspar.c
	@gcc studentspar.c -o spar -fopenmp -lm 
	@gcc studentsseq.c -o sseq -fopenmp -lm 

s:
	@./sseq

p:
	@./spar

clean:
	@rm $(TARGET)
