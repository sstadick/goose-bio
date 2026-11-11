#ifndef _WIN32
#define _FILE_OFFSET_BITS 64
#endif
#include "../include/faidx.h"
#ifdef _WIN32
#include <io.h>
#include <stdio.h>
int64_t goose_bio_faidx_seek(int32_t fd, int64_t offset) {
    return _lseeki64(fd, offset, SEEK_SET);
}
#else
#include <errno.h>
#include <unistd.h>
int64_t goose_bio_faidx_seek(int32_t fd, int64_t offset) {
    off_t result;
    do { result = lseek(fd, (off_t)offset, SEEK_SET); } while (result < 0 && errno == EINTR);
    return (int64_t)result;
}
#endif
