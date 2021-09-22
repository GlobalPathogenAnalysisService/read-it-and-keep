MINIMAP2_VERSION = 2.22

readItAndKeep: ../minimap2-$(MINIMAP2_VERSION)/libminimap2.a
	g++ -Wall -O2 -I../minimap2-$(MINIMAP2_VERSION)/ -I../gzstream/  -pthread -O2 readItAndKeep.cpp ../minimap2-$(MINIMAP2_VERSION)/libminimap2.a ../gzstream/libgzstream.a -lz -lm -o readItAndKeep

../gzstream/libgzstream.a :
	$(MAKE) -C ../gzstream/

../minimap2-$(MINIMAP2_VERSION)/libminimap2.a:
	$(MAKE) -C ../minimap2-$(MINIMAP2_VERSION)

clean:
	$(RM) readItAndKeep
	$(MAKE) -C ../minimap2-$(MINIMAP2_VERSION) clean
	$(MAKE) -C ../gzstream/ clean

test:
	cd ../tests && pytest
