MINIMAP2_VERSION = 2.22

readItAndKeep: ../minimap2-$(MINIMAP2_VERSION)/libminimap2.a
	g++ -Wall -O2 -I../minimap2-$(MINIMAP2_VERSION)/  -pthread -O2 readItAndKeep.cpp ../minimap2-$(MINIMAP2_VERSION)/libminimap2.a -lz -lm -o readItAndKeep


../minimap2-$(MINIMAP2_VERSION)/libminimap2.a:
	$(MAKE) -C ../minimap2-$(MINIMAP2_VERSION)

clean:
	$(RM) readItAndKeep
	$(MAKE) -C ../minimap2-$(MINIMAP2_VERSION) clean

test:
	cd ../tests && pytest
