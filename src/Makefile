EXT_DIR = ext
MINIMAP2_VERSION = 2.22
MINIMAP2_DIR = $(EXT_DIR)/minimap2-$(MINIMAP2_VERSION)
GZSTREAM_DIR = $(EXT_DIR)/gzstream
LIBGZSTREAM_A = $(GZSTREAM_DIR)/libgzstream.a
LIBMINIMAP2_A = $(MINIMAP2_DIR)/libminimap2.a

readItAndKeep: $(LIBMINIMAP2_A) $(LIBGZSTREAM_A)
	g++ -std=c++11 -Wall -O2 -I$(MINIMAP2_DIR)/ -I$(GZSTREAM_DIR)  -pthread -O2 readItAndKeep.cpp $(LIBMINIMAP2_A) $(LIBGZSTREAM_A) -lz -lm -o readItAndKeep

$(LIBGZSTREAM_A) :
	$(MAKE) -C $(GZSTREAM_DIR)

$(LIBMINIMAP2_A) :
	$(MAKE) -C $(MINIMAP2_DIR)

clean:
	$(RM) readItAndKeep
	$(MAKE) -C $(MINIMAP2_DIR) clean
	$(RM) $(GZSTREAM_DIR)/*.o $(GZSTREAM_DIR)/*.a

test:
	cd ../tests && pytest
