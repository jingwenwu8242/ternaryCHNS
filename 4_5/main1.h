
void initialization(double **u, double **w, double **p, double **phi1, double **phi2, double **phi3);

void cal_density(double **phi1, double **phi2,double **phi3, double den1, double den2,double den3, double **var_den);

void cal_vicosity(double **phi1, double **phi2,double **phi3, double den1, double den2,double den3, double **var_den);



void full_step(double **tu, double **tw, double **p, double **u, double **w);

void advection(double **adv_u, double **adv_w,double **u, double **w);
void temp_uw(double **tu, double **tw, double **u, double **w);
void Poisson(double **tu, double **tw, double **p);

void MG_Poisson(double **u, double **f);
void relax(double **p, double **f,double **w,  int nrt, int nzt);

void vcycle(double **uf, double **ff,double **wf, int nrf, int nzf, int ilevel);
void source(double **tu, double **tw, double **divuw, int nrt, int nzt);
void div_uw(double **tu, double **tw, double **divuw, int nrt, int nzt);

void grad_p(double **p, double **dpdr, double **dpdz, int nrt, int nzt);

void residual_den(double **r, double **u, double **f, double **den, int nrt, int nzt);

void laplace(double **p, double **lap_p, int nrt, int nzt);

void auguw(double **u, double **w, int nrt, int nzt);

void advection_c(double **u, double **w, double **c, double **adv_c);



void restrict1(double **cf, int nrf, int nzf,
              double **cc, int nrc, int nzc);

void prolong(double **cc, int nrc, int nzc,
             double **cf, int nrf, int nzf);




void surface_tension(double **phi2, double **fr, double **fz);
void sf_force(double **phi, double **mu, double **fx, double **fy) ;
/****CH****/

void cahn(double **c_old, double **cc_old, double **c_new, double **mu, double theta, double **adv_phi);

void source_ch(double **sc, double **smu, double **c_old, double **cc_old, double theta, double **adv_c);

void vcycle_ch(double **uf_new, double **Mo, double **wf_new, double **sc, double **smu, int nrf, int nzf, int ilevel);

void relax_ch(double **c_new, double **Mo, double **mu_new, double **sc,
              double **smu, int ilevel, int nrt, int nzt);
void restrict2(double **cf, double **cc, double **df, double **dc,
               int nrc, int nzc);


void prolong_ch(double **cc,  double **cf, double **dc, double **df,
                int nrc, int nzc);


void defect(double **duc, double **dwc,double **Mo, double **uf_new, double **wf_new,
            double **suf, double **swf, int nxf, int nyf);


double error(double **c_old, double **c_new, int nrt, int nzt);




void nonL(double **rc, double **rmu,double **Mo1, double **c_new, double **mu_new, int nrt, int nzt);


void laplace_Mo(double **Mo,double **a, double **lap_a, int nrt, int nzt);
void laplace_ch(double **a, double **lap_a, int nrt, int nzt);

/*******function********/
double mass_comp(double **phi);
double dfphi(double phi1);
void functiong(double **u, double **w, double **result);
void functionq(double **u, double **w, double **result);
/************/