static: landsat-index.static landsat-server.static

landsat-index.static: landsat-index.cpp lesser_landsat_scene.h projection.h rtree.hpp
	$(CXX) $(CXXFLAGS) $(GDAL_CFLAGS) \
         landsat-index.cpp $(LDFLAGS) -static-libstdc++ -static-libgcc -fopenmp -o $@
	strip $@
	chown "1000:1000" $@

landsat-server.static: server.c fullio.c fullio.h landsat-server.cpp pngwrite.c pngwrite.h
	$(CXX) $(CXXFLAGS) -fpie -pie -fstack-protector-strong -flto $(GDAL_CFLAGS) \
          -Wno-reorder -Wno-unused-parameter -Wno-write-strings -Wno-pointer-arith \
          server.c fullio.c pngwrite.c landsat-server.cpp $(LDFLAGS) -static-libstdc++ -static-libgcc -fopenmp -o $@
	strip $@
	chown "1000:1000" $@

clean-static:
	rm -f landsat-index.static landsat-server.static
