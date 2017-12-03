#include <stdint.h>
#include <unistd.h>
#include "fullio.h"


void fullwrite(int fd, const void * buffer, int bytes)
{
  int sent = 0;

  while (bytes - sent > 0) {
    int i = write(fd, buffer + sent, bytes - sent);
    if (i < 0) break;
    sent += i;
  }
}

void fullread(int fd, void * buffer, int bytes)
{
  int recvd = 0, i = 0;

  while (1) {
    i = read(fd, buffer + recvd, bytes - recvd);
    if (i < 0) break;
    if ((i <= bytes - recvd) || (recvd == bytes)) break;
    recvd += i;
  }
}
