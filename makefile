
f := example

all: src

src: FORCE
	cd src; make; cd ..

FORCE:

run: 
	cd data/$(f); ../../bin/reuleaux $(f); cd ..

clean:
	cd src; rm -rf *.o;rm -rf *.d;cd ..
	rm -f lib/reuleaux.a
	rm -f bin/reuleaux
