#include <stdio.h>
#include <stdlib.h>

#define FULLSTACK 1000000
#define EMPTYSTACK 0

typedef int qmsdata;
typedef enum{ false, true} boolean;

struct elem {
  qmsdata        d;
  struct elem *next;
  struct elem *prev;
};

typedef struct elem elem;

struct qmstack{
   int  cnt;
   elem *top;
   elem *bottom;
};

typedef struct qmstack qmstack;

void     qmsinit(qmstack *stk);
void     qmspush(qmsdata d, qmstack *stk);
qmsdata     qmstop(qmstack *stk);
qmsdata     qmspop(qmstack *stk);
boolean  qmsempty(const qmstack *stk);
boolean  qmsfull(const qmstack *stk);


