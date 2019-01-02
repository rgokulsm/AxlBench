#include <stdio.h>
#include <stdlib.h>

#define TYPE double

#define N 1048576
#define SCAN_BLOCK 16

#define SCAN_RADIX N/SCAN_BLOCK
#define BUCKETSIZE SCAN_BLOCK*SCAN_RADIX
