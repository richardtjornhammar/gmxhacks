#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#define MAXV            5000
#define EMPTY		0
#define FULL		30000
#define MAXL            30000
#define SLEN            100
#define	DIME		5
#define	DIM3		4
#define	SUP(X)	(X)*(X)

typedef enum boolean {false,true} boolean;

typedef struct dens{
  double D[20];
} dens;

typedef struct elem{
  double H[12][MAXL];
} elem;

typedef struct bund{
  int B[100];
} bund;

typedef struct CONTAINER {
  elem dat;
} CONTAINED;

typedef struct PART{
  bund dat;
} PARTITION;

typedef struct DEN{
  dens dat;
} DENSITY;

typedef CONTAINED *MAT;
typedef PARTITION *PAR;
typedef DENSITY   *DEN;

typedef struct complex {
        double real;
        double imag;
} complex;

typedef union listObject {
        double d;
        complex z;
        char c;
} listObject;

typedef struct ATOM {
	char		c;
	char		name[DIME+1];
	char		resi[DIME+1];
	char		ftyp[DIME*2];
	int		resnr;
	int		cgnr;
	float		pos[DIME+1];
	float		mass_ch[3];
	int		flag;
	int		qmmm;
	int		edge;
} ATOM;

typedef struct listLink {
//	enum {dubbel,komplex,kar} 	kind;
	ATOM data;
	struct listLink *next;
	struct listLink *prev;
} ELEMENT;

typedef ELEMENT	*LINK;

void mkListHash(LINK*, LINK );
LINK string_to_list(char*);
LINK str2list(char*);
LINK initAtom(void);
void addAtom(LINK);
LINK getAtom(LINK, int);
LINK getFirst(LINK);
void copyAtom(LINK, LINK);
void swapAtoms(LINK, LINK);
int iterCount(LINK);
void print_listc(LINK);
void printAtom(LINK);
void printAtomnr(LINK, int);
void printAtomList(LINK);
void printlink(LINK);
void printlinknr(LINK, int);
void concatenate(LINK, LINK);
void inAtomnr(LINK, LINK, int);
void insert(LINK, LINK, LINK);
void deletenr(LINK, int);
void delAtom(LINK,int);
void delete_list(LINK);
int listlength(LINK);
