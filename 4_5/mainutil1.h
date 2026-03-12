#define gnr 128
#define gnz  gnr*2

#define iloop for(i=1;i<=gnr;i++)
#define i0loop for(i=0;i<=gnr;i++)
#define kloop for(k=1;k<=gnz;k++)
#define k0loop for(k=0;k<=gnz;k++)
#define ikloop iloop kloop
#define i0kloop i0loop kloop
#define ik0loop iloop k0loop
#define i0k0loop i0loop k0loop
#define iloopt for(i=1;i<=nrt;i++)
#define i0loopt for(i=0;i<=nrt;i++)
#define kloopt for(k=1;k<=nzt;k++)
#define k0loopt for(k=0;k<=nzt;k++)
#define ikloopt iloopt kloopt 
#define i0kloopt i0loopt kloopt 
#define ik0loopt iloopt k0loopt 



void print_last(double **phi_n, double **mu, double **u_n, 
                double **w_n, double **p_nmhf, double *byn);
void augmenc(double **phi);
double zero_extrap2( double a1, double a2 );
void augc(double **c, int nrt, int nzt);
void mat_comb3(double **a, double s1, double **b, double s2, double **c,
               double s3, double **d, 
               int nrl, int nrh, int ncl, int nch);

void augp(double **p, int nrt, int nzt);
double *dvector (long nl, long nh);


double **dmatrix(long nrl, long nrh, long ncl, long nch);


void free_dvector(double *v, long nl, long nh);


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);


void zero_vector ( double *a, int nl, int nh );



void veccopy (double *a, double *b, 
              int nla, int nha,
              int nlb, int nhb);

void pressure_update2(double **a, double **b, double **c);

void vecadd(double *a, double *b, double *c, int nl, int nr);


void vecsub(double *a, double *b, double *c, int nl, int nr);


void s_vec_add(double *a, double *b, double s,
              int nl, int nr);

void s_vec_mult(double *a, double fac, double *b, 
              int nl, int nr);


void veccomb( double *a, double s1, double *b, double s2, double *c, 
            int nl, int nr );



void vecmult(double *a, double *b, double *c, int nl, int nr );

void vecdiv(double *a, double *b, double *c, int nl, int nr );


double sum_vector( double *a, int nl, int nh );


double dot_product(double *a, double *b, int nl, int nr);


void mat_add(double **a, double **b, double **c, 
              int xl, int xr, int yl, int yr);



void mat_add2(double **a, double **b, double **c, 
              double **a2, double **b2, double **c2, 
              int xl, int xr, int yl, int yr);




void zero_matrix(double **a, int xl, int xr, int yl, int yr);


void mat_copy(double **a, double **b, 
              int xl, int xr, int yl, int yr);



void mat_copy2(double **a, double **b, 
              double **a2, double **b2, 
              int xl, int xr, int yl, int yr);




void mat_sub(double **a, double **b, double **c, 
            int nrl, int nrh, int ncl, int nch);

void mat_sub2(double **a, double **b, double **c, 
              double **a2, double **b2, double **c2, 
              int nrl, int nrh, int ncl, int nch);


void matmult(double **a, double **b, double **c, 
              int nrl, int nrh, int ncl, int nch);



void mat_comb(double **a, double s1, double **b, 
             double s2, double **c, 
             int nrl, int nrh, int ncl, int nch);


void s_mat_mult(double **a, double fac, double **b, 
              int nrl, int nrh, int ncl, int nch);




int s_mat_mult2(double **a, int arl, int arh, int acl, int ach,
                double fac, 
                double **b, int brl, int brh, int bcl, int bch);

double sum_matrix(double **a, 
                  int nrl, int nrh, int ncl, int nch);

double max_norm(double **a, double **b,
                int nrl, int nrh, int ncl, int nch);

double mat_normal2(double **a, 
               int nrl, int nrh, int ncl, int nch);

double mat_max(double **a, 
               int nrl, int nrh, int ncl, int nch);

void mat_mat_mult(double **a, double **b, double **c,
              int nrl, int nrh, int ncl, int nch);

int mat_mat_mult2(double **a, int arl, int arh, int acl, int ach,
                   double **b, int brl, int brh, int bcl, int bch,
                   double **c, int crl, int crh, int ccl, int cch);

void convert_mat_vec(int which_way, 
                double **a, int nrl, int nrh, int ncl, int nch,
                double *b, int nl, int nh);

void mat_vec_mult(double *a, double **b, double *c,
              int nrl, int nrh, int ncl, int nch);



void print_data(double **u, double **v, double **phi1, double **phi2,double **phi3,double **p);


void pressure_update(double **a);

void augmen_phi(double **phi, double **aug_phi);


void aug_uw(double **aug_u, double **aug_v, double **u, double **v, 
            int nxt, int nyt);
