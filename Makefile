CC      = clang++
LD      = clang++
CFLAGS  = -Wall -O2 -I. -c
SUB     = submission

INCS = Image.h FFT.h 

OBJS = Image.o FFT.o main.o

all: fft

clean:
	rm -rf $(OBJS) fft $(SUB)

main.o: $(INCS) main.cxx
	$(CC) $(CFLAGS) main.cxx -o main.o

Image.o: $(INCS) Image.cxx
	$(CC) $(CFLAGS) Image.cxx -o Image.o

FFT.o: $(INCS) FFT.cxx
	$(CC) $(CFLAGS) FFT.cxx -o FFT.o

fft: $(OBJS)
	$(LD) $(OBJS) -o fft


submission: all
	mkdir -p $(SUB)
	./fft jobs.pgm $(SUB)/jobsLPF20.pgm -f 1 -pmax 0.2
	./fft jobs.pgm $(SUB)/jobsHPF20.pgm -f 2 -pmin 0.2
	./fft jobs.pgm $(SUB)/jobsBPF1005.pgm -f 3 -pmax 0.1 -pmin 0.05
	./fft jobs.pgm $(SUB)/jobsBPF2010.pgm -f 3 -pmax 0.2 -pmin 0.1
	./fft queenstower.pgm $(SUB)/queenstowerLPF20.pgm -f 1 -pmax 0.2
	./fft queenstower.pgm $(SUB)/queenstowerHPF20.pgm -f 2 -pmin 0.2
	./fft queenstower.pgm $(SUB)/queenstowerBPF1005.pgm -f 3 -pmax 0.1 -pmin 0.05
	./fft queenstower.pgm $(SUB)/queenstowerBPF2010.pgm -f 3 -pmax 0.2 -pmin 0.1
