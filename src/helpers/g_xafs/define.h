//	gcc list.c tools_*.c g_orient.c -o g_orient -lm
// gcc list.c tools_*.c g_exafs_new.c -o g_exafs -lm &&  cp g_exafs /home/richardt/bin/.
#define nm        0.000000001
#define k_b       1.3806503*nm*nm*nm*10000
#define e_q       1.60217646*nm*nm*0.1
#define epS       78.000000
#define Els       -2.837297
#define c_0       299792458.00
#define pi4       12.5663706143592
#define mu0       (pi4)*0.0000001
#define ep0       1/(c_0*c_0*mu0)
#define NA        6.02214179/nm/nm*100000
#define SQ(X)     (X)*(X)
#define CB(X)     (X)*SQ(X)
#define SGN(X)    ((X)>=0?1.0:-1.0)
#define RESX      13
#define RESY      10000
#define RHOw      33.427
#define MUw       0.0473

typedef struct vect{
  float* dat ;
  int i;
} vector;

typedef vector *pvec;

enum {
  etGRO=0, etPDB=1, etXYZ=2, etTOP=3, etITP=4, etDAT=5
};

enum {
  cmdROT=0,cmdMROT=1,cmdPUT=2,cmdTRANS=3,cmdERROR=4
};

void fatal(char*);
void pdb_create(LINK, LINK, FILE*, FILE*, float*);
void write_pdb(LINK, LINK, FILE*, FILE*, float*, int);
void write_pdb_tot(LINK, FILE*, float*, int, int);
void gro_create(LINK, LINK, FILE*, float*);
void trans_top(LINK, LINK, float*, int, int);
void merge_top(LINK, LINK, LINK, LINK, float*);
void assign_topology(LINK, LINK, LINK, FILE*, int);
LINK centerAtom(LINK, char*);
void init_itp(LINK, FILE*);
void gro_update(LINK, LINK, LINK, FILE*, float*, float*, int);
char* getAtomName(int);
float dotP(float*, float*);
int getAtomNr(char*);
void calc_dipole_rad(LINK, LINK, pvec, float*, int, int);
void calc_ndens(LINK*, LINK ,float , float, pvec, pvec, int, float*, pvec, int);
void calc_gID(pvec, float, float, float, pvec,int);
void calc_Polarisation(pvec, pvec, pvec, float,int);
void calc_simpPolarisation(LINK,LINK,float*,pvec,float,float,int,int);
void initPvec(pvec, int);
void update_topology(LINK, LINK, LINK, FILE*, int, float*, float*, int);
void printResults(float[RESX][RESY], FILE*, int);
void calc_densDIR(LINK[] , int, LINK, float , float , float[6] , float[3] , pvec, pvec, pvec, pvec, pvec, pvec, pvec, pvec,float[RESX][RESY], int, int, int,int,int);
void calc_cID(pvec, pvec, pvec, int);
void corner_search(LINK*, float*);
void calc_atomPOL(LINK*, int, LINK, float, float, float*, float*, pvec, pvec, float[RESX][RESY], int, int, int, int, int);
float rmPBC(LINK ,LINK , float[6]);
float rmPBCsol(LINK ,LINK , float[6]);
void scale(LINK,LINK,float);
