#define DEBUG    0
#define PI       3.1415926535898
#define SIZE     3
#define DIM      3
#define MAXEL    1000
#define VERBO    1
#define MAXFILES 1000

typedef struct elem {
  double  R[DIM];
  int     nclust;
} elem;

typedef struct eleme {
  double R[DIM];
  double RAD;
  int NPOINTS;
  int IND[MAXEL];
} alam;

typedef struct CONTR {
  elem dat;
} CONTAINER;

typedef struct CONTD {
  alam dat;
} CONTAINED;

typedef CONTAINER *POINT;
typedef CONTAINED *CLUST;

void clusterQM(t_QMrec*);
void clusterMM(t_MMrec*,t_QMrec*);
void showCL(CLUST,FILE*);
void mergeCL(CLUST,CLUST);
void one_cluster(t_MMrec*,t_QMrec*);
