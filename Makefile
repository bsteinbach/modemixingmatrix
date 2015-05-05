
modemix: modemix.c
	gcc -o modemix modemix.c -O2 -Wall -lm

clean:
	rm -f modemix
