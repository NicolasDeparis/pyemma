all: src/*
	gcc -fPIC -shared src/argsort.c -o libsort.so -lm
	mv -f libsort.so /usr/local/lib/
	gcc -fPIC -shared  src/kdtree.c -o  libkdtree.so  -lsort -lm
	mv -f libkdtree.so /usr/local/lib/
	cp src/*.h /usr/local/include/

local: src/*
	gcc -fPIC -shared src/argsort.c -o libsort.so -lm
	gcc -fPIC -shared  src/kdtree.c -o  libkdtree.so  -lsort -lm
