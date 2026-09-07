
void initialization(double **u, double **w, double **p, double **phi1, double **phi2, double **phi3);

void cal_density(double **phi1, double **phi2,double **phi3, double den1, double den2,double den3, double **var_den);

void cal_vicosity(double **phi1, double **phi2,double **phi3, double den1, double den2,double den3, double **var_den);



void full_step(double **tu, double **tw, double **p, double **u, double **w);

void advection(double **adv_u, double **adv_w,double **u, double **w);
void temp_uw(double **tu, double **tw, double **u, double **w);
void Poisson_stage(double **tu, double **tw, double **p, double alpha);

void source_stage(double **tu, double **tw, double **divuw,
                  int nrt, int nzt, double alpha);
void Poisson(double **tu, double **tw, double **p);

void solve_Poisson_relaxation(double **u, double **f);
void relax(double **p, double **f,double **w,  int nrt, int nzt);

void source(double **tu, double **tw, double **divuw, int nrt, int nzt);
void div_uw(double **tu, double **tw, double **divuw, int nrt, int nzt);

void grad_p(double **p, double **dpdr, double **dpdz, int nrt, int nzt);

void laplace(double **p, double **lap_p, int nrt, int nzt);

void auguw(double **u, double **w, int nrt, int nzt);

void advection_c(double **u, double **w, double **c, double **adv_c);



void surface_tension(double **phi2, double **fr, double **fz);
void sf_force(double **phi, double **mu, double **fx, double **fy) ;
/****CH****/

void cahn(double **c_old, double **cc_old, double **c_new, double **mu, double theta, double **adv_phi);

void source_ch(double **sc, double **smu, double **c_old, double **cc_old, double theta, double **adv_c);

void relax_ch(double **c_new, double **Mo, double **mu_new, double **sc,
              double **smu, int nrt, int nzt);


double error(double **c_old, double **c_new, int nrt, int nzt);




void laplace_ch(double **a, double **lap_a, int nrt, int nzt);

/*******function********/
double mass_comp(double **phi);
double dfphi(double phi1);
void functiong(double **u, double **w, double **result);
void functionq(double **u, double **w, double **result);
/************/
