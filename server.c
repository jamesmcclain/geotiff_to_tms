/*
 * Copyright (c) 2017, James McClain
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 * 3. All advertising materials mentioning features or use of this
 *    software must display the following acknowledgement: This product
 *    includes software developed by Dr. James W. McClain.
 * 4. Neither the names of the authors nor the names of the
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <netinet/in.h>
#include <netinet/ip.h>
#include <signal.h>
#include <sys/socket.h>

#include "ansi.h"
#include "load.h"

#define P (1<<5)


// http://localhost:8001/{z}/{x}/{y}.png
int main(int argc, const char ** argv)
{
  struct sockaddr_in sa;
  int fd1, fd2;
  uint8_t yes = 1;
  int p = P;

  sa.sin_family = AF_INET;
  sa.sin_addr.s_addr = htonl(INADDR_ANY);
  sa.sin_port = htons(8001);

  /* Number of servers to prefork */
  if (argc > 1)
    sscanf(argv[1], "%d", &p);

#if 0
  load(1);
  int fd = open("/tmp/tile.png", O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR);
  zxy(fd, 3, 5, 3, 1);
  /* zxy(fd, 10, 759, 448, 1); */
  fsync(fd); close(fd);
#endif

#if 1
  /* Create socket */
  if ((fd1 = socket(PF_INET, SOCK_STREAM, 0)) == -1) {
    fprintf(stderr, ANSI_COLOR_RED "'socket(2)' problem" ANSI_COLOR_RESET "\n");
    exit(-1);
  }
  setsockopt(fd1, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes));

  /* Bind */
  if (bind(fd1, (const struct sockaddr *)&sa, sizeof(sa)) == -1) {
    fprintf(stderr, ANSI_COLOR_RED "'bind(2)' problem" ANSI_COLOR_RESET "\n");
    exit(-1);
  }

  /* Listen */
  if (listen(fd1, 1<<10) == -1) {
    fprintf(stderr, ANSI_COLOR_RED "'listen(2)' problem" ANSI_COLOR_RESET "\n");
    close(fd1);
    exit(-1);
  }

  signal(SIGPIPE, SIG_IGN);
  for (int i = 0; (i < p-1) && fork(); ++i);

  /* Initialize backend */
  load(1);

  /* Handle requests */
  while (1) {
    char * twohundred =
      "HTTP/1.1 200 OK\r\n"
      "Content-Type: image/png\r\n\r\n";
    char buffer[1<<5];
    int z, x, y;

    if ((fd2 = accept(fd1, NULL, NULL)) == -1) {
      fprintf(stderr, ANSI_COLOR_RED "'accept(2)' problem" ANSI_COLOR_RESET "\n");
      close(fd1);
      close(fd2);
      exit(-1);
    }

    read(fd2, buffer, sizeof(buffer));
    sscanf(buffer, "GET /%d/%d/%d.png HTTP/1.1", &z, &x, &y);
    write(fd2, twohundred, strlen(twohundred));
    zxy(fd2, z, x, y, 1);

    fsync(fd2); shutdown(fd2, SHUT_RDWR); close(fd2);
  }

  fsync(fd1); shutdown(fd1, SHUT_RDWR); close(fd1);
#endif

  return 0;
}
