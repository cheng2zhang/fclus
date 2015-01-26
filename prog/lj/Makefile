CC = icc -Wall -Wremarks -O2 -g
LM =
#CC = gcc -O3 -Wall -Wextra -g
#LM = -lm

CCd = gcc -Wall -Wextra -O0 -g
LMd = -lm

srcs = $(wildcard *.c)
bins = $(patsubst %.c,%, $(srcs))
bin2ds = $(patsubst %,%2d, $(bins))
bins_d = $(patsubst %,%_d, $(bins))
bin2ds_d = $(patsubst %,%2d_d, $(bins))
deps = $(wildcard *.h) Makefile

all: $(bins)

$(bins): % : %.c $(deps)
	$(CC) $< -o $@ $(LM)

$(bin2ds): %2d : %.c $(deps)
	$(CC) -DD=2 $< -o $@ $(LM)

$(bins_d): %_d : %.c $(deps)
	$(CCd) $< -o $@ $(LMd)

$(bins2ds_d): %2d_d : %.c $(deps)
	$(CCd) $< -o $@ $(LMd)

clean:
	rm -rf $(bins) $(bin2ds) $(bins_d) $(bins2ds_d) \
	  a.out *~ .*.un~ */*~ */.*.un~ *.dat *.pos *.log
	rstrip.py -Rlv

Dropbox: clean
	$(MAKE) -C .. $@

