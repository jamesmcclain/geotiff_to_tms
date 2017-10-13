#ifndef __FULLWRITE_H__
#define __FULLWRITE_H__

void fullwrite(int fd, const void * buffer, int bytes);
void fullread(int fd, void * buffer, int bytes);

#endif
