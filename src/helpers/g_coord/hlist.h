#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#define MAXV            5000
#define EMPTY		0
#define FULL		40000
#define MAXL            40000
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

void search_top(LINK*, int, float*,int*, float*, 
		float*, MAT, float, PAR);
