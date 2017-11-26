#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#define DEFAULT_INDEXFILE "/tmp/index.data"
#define DEFAULT_LIST_PREFIX "https://s3-us-west-2.amazonaws.com/landsat-pds/"
#define DEFAULT_READ_PREFIX "/vsicurl/https://s3-us-west-2.amazonaws.com/landsat-pds/"
#define RADIUS (6378137.0)
#define SMALL_TILE_SIZE (1<<6)
#define SMALL_TILE_SIZE2 (SMALL_TILE_SIZE * SMALL_TILE_SIZE)
#define STRING_BUFFER_SIZE (1<<10)
#define STRING_LEN (1<<10)
#define TILE_SIZE (1<<8)
#define TILE_SIZE2 (TILE_SIZE * TILE_SIZE)
#define WEBMERCATOR "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"

#endif
