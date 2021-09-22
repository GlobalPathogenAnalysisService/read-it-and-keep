EXT_DIR = ext
MINIMAP2_VERSION = 2.22
MINIMAP2_DIR = $(EXT_DIR)/minimap2-$(MINIMAP2_VERSION)
GZSTREAM_DIR = $(EXT_DIR)/gzstream

readItAndKeep: $(MINIMAP2_DIR)/libminimap2.a
	g++ -Wall -O2 -I$(MINIMAP2_DIR)/ -I$(GZSTREAM_DIR)  -pthread -O2 readItAndKeep.cpp $(MINIMAP2_DIR)/libminimap2.a $(GZSTREAM_DIR)/libgzstream.a -lz -lm -o readItAndKeep

$(GZSTREAM_DIR)/libgzstream.a :
	$(MAKE) -C $(GZSTREAM_DIR)

$(MINIMAP2_DIR)/libminimap2.a:
	$(MAKE) -C $(MINIMAP2_DIR)

clean:
	$(RM) readItAndKeep
	$(MAKE) -C $(MINIMAP2_DIR) clean
	$(MAKE) -C $(GZSTREAM_DIR) clean

test:
	cd ../tests && pytest
