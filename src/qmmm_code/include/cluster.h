#define DEBUG    0
#define PI       3.1415926535898
#define SIZE     3
#define DIM      3
#define MAXEL    100000
#define VERBO    1
#define MAXFILES 10000

void printClustering(t_QMrec*, FILE*);
void sortClusters(t_QMrec*);
void qmmmPBC(t_QMrec*, rvec[2]);
void clusterAlgorithm(t_QMrec*, FILE*);
void makeDistance(t_QMrec*, FILE*,t_pbc*);
void clusterPositions(t_QMrec*);
void clustering(t_QMrec*, t_pbc*);
void one_cluster(t_MMrec*, t_QMrec*);
void clusterMM(t_MMrec*, t_QMrec*, t_pbc*);
double roundChrg(double);
int anint000(double r);
int anint001(double r);
int anint002(double r);
