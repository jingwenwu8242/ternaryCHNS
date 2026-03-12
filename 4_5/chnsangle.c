/*********************************

       Axisymmetric NSCH

variable density
variable viscosity

complex;
*********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "main1.h"
#include "mainutil1.h"
#define NR_END 1
#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)
int nr, nz, n_level, p_relax, c_relax, it;
double pi, dr, dz, dr2, dz2, rleft, rright, zleft, zright, h, dt, gam, Cahn, SS, lambda, **beta, **Mo1, **Mo2, **normphi3, **phi3, kappa,
    **worku, **workw, **workp, **adv_u, **adv_w, **adv_phi1, **adv_phi2, **muphi1, **omuphi1, **intmuphi1, **muphi2, **fr, **fz, **frphi1, **fzphi1, **frphi2, **fzphi2,
    vis1, vis2, vis3, mobil, sig, Re, **vis, rho1, rho2, rho3, velocity, **rho_v, gravity, Fr, Re, We, Pe, **intphi1, **intphi2, **intphi3, **intu, **intw,
    theta, theta1, theta2;

int main()
{
    extern int nr, nz, n_level, p_relax, c_relax, it;
    extern double pi, dr, dz, dr2, dz2, rleft, rright, zleft, zright, h, dt, gam, Cahn, SS, lambda, **beta, **Mo1, **Mo2, **normphi3, **phi3, kappa,
        **worku, **workw, **workp, **adv_u, **adv_w, **adv_phi1, **adv_phi2, **muphi1, **omuphi1, **intmuphi1, **muphi2, **fr, **fz, **frphi1, **fzphi1, **frphi2, **fzphi2,
        vis1, vis2, vis3, mobil, sig, **vis, rho1, rho2, rho3, velocity, **rho_v, gravity, Fr, Re, We, **intphi1, **intphi2, **intphi3, **intu, **intw,
        theta, theta1, theta2;

   
    int max_it, ns, i, k, count = 0;
    double **ou, **ow, **u, **w, **p, **np, **ophi1, **phi1, **nu, **nw, **nphi1, **ophi2, **phi2, **nphi2, **ophi3, T, mass1, mass2, res1, res2, red2, Pe, Cn;
 
    FILE *fu, *fw, *fp, *fphi1, *fphi2, *fphi3, *my, *mymass1, *mymass2;

    pi = 4.0 * atan(1.0);

    p_relax = 20;
    c_relax = 5;

    nr = gnr;
    nz = gnz;

    n_level = (int)(log(nr) / log(2) - 0.9);
    n_level -= 0;
    rleft = 0, rright = pi;
    zleft = 0.0, zright = 2 * pi;

    dr = rright / (double)nr;
    dz = zright / (double)nz;
    dr2 = pow(dr, 2);
    dz2 = pow(dz, 2);
    h = dr;

    /***********************/
    T = 1.0;
    dt = 0.2 * h * h;
    max_it = 40000;
    ns = (int)(max_it / 100 + 0.001);
    gam = 5.0 * h / (4.0 * sqrt(2.0) * atanh(0.9));
    Cahn = pow(gam, 2);

    rho1 = 1;
    rho2 = 1;
    rho3 = 1;
    vis1 = 1;
    vis2 = 1;
    vis3 = 1; // 
    Re = 1;
    We = 0.5;
    Pe = 0.001 / gam;
    mobil = 1 / Pe; 
    Cn = gam;
    sig = 1.0;

    lambda = 0.000;
    SS = 2.0;

    gravity = 0;
    Fr = 1;
    velocity = 0;
    theta = 0;
    theta1 = (180 - theta) * pi / 180;
    theta2 = theta * pi / 180.0; 
    kappa = 1e-8;
    /***********************/

    printf("nr = %d, nz = %d\n", nr, nz);
    printf("dt      = %f\n", dt);
    printf("max_it  = %d\n", max_it);
    printf("ns      = %d\n", ns);
    printf("n_level = %d\n", n_level);
    printf("sig      = %f\n", sig);
    printf("Cahn      = %f\n\n", Cahn);

 
    /********CH****/
    ophi1 = dmatrix(0, nr + 1, 0, nz + 1);
    phi1 = dmatrix(0, nr + 1, 0, nz + 1);
    nphi1 = dmatrix(0, nr + 1, 0, nz + 1);
    ophi2 = dmatrix(0, nr + 1, 0, nz + 1);
    phi2 = dmatrix(0, nr + 1, 0, nz + 1);
    nphi2 = dmatrix(0, nr + 1, 0, nz + 1);

    phi3 = dmatrix(0, nr + 1, 0, nz + 1);
    Mo1 = dmatrix(0, nr + 1, 0, nz + 1);
    Mo2 = dmatrix(0, nr + 1, 0, nz + 1);
    beta = dmatrix(1, nr, 1, nz);
    normphi3 = dmatrix(1, nr, 1, nz);

    muphi1 = dmatrix(0, nr + 1, 0, nz + 1);
    muphi2 = dmatrix(0, nr + 1, 0, nz + 1);
    omuphi1 = dmatrix(0, nr + 1, 0, nz + 1);
    intmuphi1 = dmatrix(0, nr + 1, 0, nz + 1);

    intphi1 = dmatrix(0, nr + 1, 0, nz + 1);
    intphi2 = dmatrix(0, nr + 1, 0, nz + 1);


    ou = dmatrix(-1, nr + 1, 0, nz + 1);
    u = dmatrix(-1, nr + 1, 0, nz + 1);
    nu = dmatrix(-1, nr + 1, 0, nz + 1);
    ow = dmatrix(0, nr + 1, -1, nz + 1);
    w = dmatrix(0, nr + 1, -1, nz + 1);
    nw = dmatrix(0, nr + 1, -1, nz + 1);
    p = dmatrix(0, nr + 1, 0, nz + 1);

    adv_u = dmatrix(0, nr, 1, nz);
    adv_w = dmatrix(1, nr, 0, nz);
    adv_phi1 = dmatrix(1, nr, 1, nz);
    adv_phi2 = dmatrix(1, nr, 1, nz);

    fr = dmatrix(0, nr, 1, nz);
    fz = dmatrix(1, nr, 0, nz);
    frphi1 = dmatrix(0, nr, 1, nz);
    fzphi1 = dmatrix(1, nr, 0, nz);
    frphi2 = dmatrix(0, nr, 1, nz);
    fzphi2 = dmatrix(1, nr, 0, nz);

    worku = dmatrix(0, nr, 1, nz); 
    workw = dmatrix(1, nr, 0, nz);
    workp = dmatrix(1, nr, 1, nz); 

    vis = dmatrix(0, nr + 1, 0, nz + 1);
    rho_v = dmatrix(0, nr + 1, 0, nz + 1);

    initialization(u, w, p, phi1, phi2, phi3);

    fu = fopen("E: \\u.m", "w");
    fw = fopen("E: \\w.m", "w");
    fp = fopen("E: \\p.m", "w");
    fphi1 = fopen("E: \\phi1.m", "w");
    fphi2 = fopen("E: \\phi2.m", "w");
    fphi3 = fopen("E: \\phi3.m", "w");
    fclose(fu);
    fclose(fw);
    fclose(fp);
    fclose(fphi1);
    fclose(fphi2);
    fclose(fphi3);

    print_data(u, w, phi1, phi2, phi3, p);

    mat_copy(nphi1, phi1, 1, nr, 1, nz);
    mat_copy(nphi2, phi2, 1, nr, 1, nz);
    mat_copy(ophi1, phi1, 1, nr, 1, nz);
    mat_copy(ophi2, phi2, 1, nr, 1, nz);
    mat_copy(nu, u, 0, nr, 1, nz);
    mat_copy(nw, w, 1, nr, 0, nz);
    mat_copy(ou, u, 0, nr, 1, nz);
    mat_copy(ow, w, 1, nr, 0, nz);

    augc(phi3, nr, nz);
    ikloop
    {

        Mo1[i][k] = 1.0 - phi3[i][k];
        normphi3[i][k] = sqrt(pow((phi3[i + 1][k] - phi3[i - 1][k]) / (2 * h), 2) + pow((phi3[i][k + 1] - phi3[i][k - 1]) / (2 * h), 2));
        if (phi3[i][k] > 0.5)
        {
            Mo2[i][k] = 0;
        }
        else
        {
            Mo2[i][k] = 1;
        }
    }
    augc(Mo1, nr, nz);
    augc(Mo2, nr, nz);



    for (it = 1; it <= max_it; it++)
    {

        ikloop
        {

            intphi1[i][k] = 2.0 * phi1[i][k] - ophi1[i][k];
            intphi2[i][k] = 2.0 * phi2[i][k] - ophi2[i][k];
            intmuphi1[i][k] = 2.0 * muphi1[i][k] - omuphi1[i][k];
        }
        augc(intphi1, nr, nz);
        augc(intphi2, nr, nz);


        cal_vicosity(intphi1, intphi2, phi3, vis1, vis2, vis3, vis);
        cal_density(intphi1, intphi2, phi3, rho1, rho2, rho3, rho_v);
        augc(vis, nr, nz);
        augc(rho_v, nr, nz);

        surface_tension(intphi1, frphi1, fzphi1);

        i0kloop
        {
            fr[i][k] = -frphi1[i][k] / We;

       
        }
        ik0loop
        {
            fz[i][k] = -fzphi1[i][k] / We;
        }

        full_step(nu, nw, p, u, w);

        advection_c(nu, nw, intphi1, adv_phi1);
        advection_c(nu, nw, intphi2, adv_phi2);


        ikloop
        {
            beta[i][k] = -1 / 3.0 * (dfphi(intphi1[i][k]) + dfphi(intphi2[i][k]) + dfphi(phi3[i][k]) + gam * intphi1[i][k] * (1 - intphi1[i][k]) * normphi3[i][k] * cos(theta1) / sqrt(2.0) + gam * intphi2[i][k] * (1 - intphi2[i][k]) * normphi3[i][k] * cos(theta2) / sqrt(2.0));
        }

        cahn(phi1, ophi1, nphi1, muphi1, theta1, adv_phi1);
        cahn(phi2, ophi2, nphi2, muphi2, theta2, adv_phi2);
        
        mat_copy(ophi1, phi1, 1, nr, 1, nz);
        mat_copy(ophi2, phi2, 1, nr, 1, nz);
        mat_copy(omuphi1, muphi1, 1, nr, 1, nz);
        mat_copy(phi1, nphi1, 1, nr, 1, nz);
        mat_copy(phi2, nphi2, 1, nr, 1, nz);

        mat_copy(ou, u, 0, nr, 1, nz);
        mat_copy(ow, w, 1, nr, 0, nz);
        mat_copy(u, nu, 0, nr, 1, nz);
        mat_copy(w, nw, 1, nr, 0, nz);



        if (it % ns == 0)
        {
            print_data(u, w, phi1, phi2, phi3, p);
            count++;
            printf("print out counts %d \n", count);

           
        }
    }


    return 0;
}

void initialization(double **u, double **w, double **p, double **phi1, double **phi2, double **phi3)
{
    extern int nr, nz;
    extern double h, gam;

    int i, k;
    double r, z, psi, R0, R1, alpha;
    R0 = 0.5;
    R1 = 0.2, alpha = 0.05;
    ikloop
    {

        r = ((double)i - 0.5) * h;
        z = ((double)k - 0.5) * h;


        phi3[i][k] = 0.5 + 0.5 * tanh(-(fabs(r + 0) - R1) / (0.1 * sqrt(2.0) * gam));
        phi1[i][k] = (1 - phi3[i][k]) * (0.5 + 0.5 * tanh((R0 - r + alpha * cos(z)) / ( 2*sqrt(2.0) * gam)));
        
        phi2[i][k] = 1.0 - phi1[i][k] - phi3[i][k];
    }

    zero_matrix(u, 0, nr + 1, 0, nz);
    zero_matrix(w, 0, nr, 0, nz + 1);
    zero_matrix(p, 0, nr + 1, 0, nz + 1);

    ik0loop
    {
        if (phi1[i][k] > 0.02)
        {
            w[i][k] = (0.5 * (phi1[i][k] + phi1[i][k + 1])) * velocity;
        }
       
    }
}

void cal_density(double **phi1, double **phi2, double **phi3, double den1, double den2, double den3, double **var_den)

{
    int i, k;

    ikloop
    {
    
        var_den[i][k] = 1; 
    }
}

void cal_vicosity(double **phi1, double **phi2, double **phi3, double den1, double den2, double den3, double **var_den)

{
    int i, k;
    ikloop
    {
     
        var_den[i][k] = 1; 
    }
}



void augp(double **p, int nrt, int nzt)
{ 

    int i, k;

    for (i = 1; i <= nrt; i++)
    {

        p[i][0] = p[i][1];
        p[i][nzt + 1] = p[i][nzt];
    }

    for (k = 0; k <= nzt + 1; k++)
    {

        p[0][k] = p[1][k];
        p[nrt + 1][k] = p[nrt][k];
    }
}

void augc(double **c, int nrt, int nzt) 
{                                       
    int i, k;

    for (k = 1; k <= nzt; k++)
    {
        c[0][k] = c[1][k];
        c[nrt + 1][k] = c[nrt][k];
    }

    for (i = 0; i <= nrt + 1; i++)
    {
        c[i][0] = c[i][1];
        c[i][nzt + 1] = c[i][nzt];
    }
}
void auguw(double **u, double **w, int nrt, int nzt)
{

    int i, k;
    extern double velocity;


    for (k = 1; k <= nzt; k++)
    {
     
        u[nrt][k] = 0;
        u[-1][k] = 0;
        u[0][k] = 0 - u[1][k];
        u[nrt + 1][k] = 0 - u[nrt - 1][k];
    }
  
    for (i = -1; i <= nrt + 1; i++)
    {
        u[i][0] = 0 - u[i][1];
        u[i][nzt + 1] = 0 - u[i][nzt];
    }

    for (k = 0; k <= nzt; k++)
    {
        w[0][k] = w[1][k];
        w[nrt + 1][k] = 0 - w[nrt][k];
    }

    for (i = 0; i <= nrt + 1; i++)
    {
        w[i][0] = w[i][nzt] = 0;
        w[i][-1] = 0 - w[i][1];
        w[i][nzt + 1] = 0 - w[i][nzt - 1];
    }
}


void full_step(double **tu, double **tw, double **p, double **u, double **w)
{

    extern double **adv_u, **adv_w;

    advection(adv_u, adv_w, u, w);

    temp_uw(tu, tw, u, w);
    auguw(tu, tw, nr, nz);

    Poisson(tu, tw, p);
}

void advection(double **adv_u, double **adv_w, double **u, double **w)
{
    extern int nr, nz;
    extern double h;

    int i, k;
    double r_c, r_p, r_m;

    auguw(u, w, nr, nz);
    i0kloop
    {
        r_c = ((double)i - 0.5 + 0.5) * dr;
        r_p = r_c + 0.5 * dr;
        r_m = r_c - 0.5 * dr;

        adv_u[i][k] = u[i][k] * (u[i + 1][k] - u[i - 1][k]) / (2 * h) + 0.25 * (r_p * w[i + 1][k] + r_m * w[i][k] + r_p * w[i + 1][k - 1] + r_m * w[i][k - 1]) * (u[i][k + 1] - u[i][k - 1]) / (2 * r_c * h);
    }
    ik0loop
    {
        r_c = ((double)i - 0.5) * dr;
        r_p = r_c + 0.5 * dr;
        r_m = r_c - 0.5 * dr;
    
        adv_w[i][k] = 0.25 * (r_p * (u[i][k + 1] + u[i][k]) + r_m * (u[i - 1][k + 1] + u[i - 1][k])) * (w[i + 1][k] - w[i - 1][k]) / (2 * r_c * h) + w[i][k] * (w[i][k + 1] - w[i][k - 1]) / (2 * h);
    }
}

void advection_c(double **u, double **w, double **phi1, double **adv_phi1) // B;
{
    extern int nr, nz;
    extern double h;

    int i, k;
    double r_c, r_p, r_m;

    augc(phi1, nr, nz);

    ikloop
    {
        r_c = ((double)i - 0.5) * h;
        r_p = r_c + 0.5 * h;
        r_m = r_c - 0.5 * h;
    
        adv_phi1[i][k] = 0.0;
        adv_phi1[i][k] += (r_p * u[i][k] * (phi1[i + 1][k] + phi1[i][k]) - r_m * u[i - 1][k] * (phi1[i][k] + phi1[i - 1][k])) / (2.0 * r_c * h);
        adv_phi1[i][k] += (w[i][k] * (phi1[i][k] + phi1[i][k + 1]) - w[i][k - 1] * (phi1[i][k] + phi1[i][k - 1])) / (2.0 * h);
    }
}


void temp_uw(double **tu, double **tw, double **u, double **w)
{
    extern int nr, nz;
    extern double dt, **phi3, kappa, dt;

    int i, k;
    double **temp1_u, **temp1_w, **temp1_g, **temp1_q;

    temp1_u = dmatrix(-1, nr + 1, 0, nz + 1);
    temp1_g = dmatrix(-1, nr + 1, 0, nz + 1);
    temp1_w = dmatrix(0, nr + 1, -1, nz + 1);
    temp1_q = dmatrix(0, nr + 1, -1, nz + 1);

    functiong(u, w, temp1_g);
    functionq(u, w, temp1_q);

    i0kloop
    {
        temp1_u[i][k] = u[i][k] / (1.0 + dt * 0.5 * (phi3[i + 1][k] + phi3[i][k]) / (0.5 * (rho_v[i + 1][k] + rho_v[i][k]) * kappa)) + dt * temp1_g[i][k];
    }

    ik0loop
    {
        temp1_w[i][k] = w[i][k] / (1.0 + dt * 0.5 * (phi3[i][k + 1] + phi3[i][k]) / (0.5 * (rho_v[i][k + 1] + rho_v[i][k]) * kappa)) + dt * temp1_q[i][k];
    }

    functiong(temp1_u, temp1_w, temp1_g);
    functionq(temp1_u, temp1_w, temp1_q);

    i0kloop
    {
        tu[i][k] = u[i][k] / (2 * (1.0 + dt * 0.5 * (phi3[i + 1][k] + phi3[i][k]) / (0.5 * (rho_v[i + 1][k] + rho_v[i][k]) * kappa))) + temp1_u[i][k] / (2 * (1.0 + dt * 0.5 * (phi3[i + 1][k] + phi3[i][k]) / (0.5 * (rho_v[i + 1][k] + rho_v[i][k]) * kappa))) + dt * temp1_g[i][k] / 2;
    }
    ik0loop
    {
        tw[i][k] = w[i][k] / (2 * (1.0 + dt * 0.5 * (phi3[i][k + 1] + phi3[i][k]) / (0.5 * (rho_v[i][k + 1] + rho_v[i][k]) * kappa))) + temp1_w[i][k] / (2 * (1.0 + dt * 0.5 * (phi3[i][k + 1] + phi3[i][k]) / (0.5 * (rho_v[i][k + 1] + rho_v[i][k]) * kappa))) + dt * temp1_q[i][k] / 2;
    }

    free_dmatrix(temp1_u, -1, nr + 1, 0, nz + 1);
    free_dmatrix(temp1_g, -1, nr + 1, 0, nz + 1);
    free_dmatrix(temp1_w, 0, nr + 1, -1, nz + 1);
    free_dmatrix(temp1_q, 0, nr + 1, -1, nz + 1);
}

void Poisson(double **tu, double **tw, double **p)
{
    extern int nr, nz;
    extern double **workp, **worku, **workw, **rho_v, **Mo2;

    int i, k;

    source(tu, tw, workp, nr, nz);

    MG_Poisson(p, workp);

    grad_p(p, worku, workw, nr, nz);
    i0kloop
    {
        tu[i][k] = tu[i][k] - dt * worku[i][k] / (0.5 * (rho_v[i + 1][k] + rho_v[i][k]));
        
    }
    ik0loop
    {
        tw[i][k] = tw[i][k] - dt * workw[i][k] / (0.5 * (rho_v[i][k + 1] + rho_v[i][k]));
        
    }
}
void source(double **tu, double **tw, double **divuw, int nrt, int nzt)
{
    extern double dt;

    int i, k;

    div_uw(tu, tw, divuw, nrt, nzt);
    ikloopt
    {
        divuw[i][k] = (divuw[i][k]) / dt;
    }
}

void div_uw(double **tu, double **tw, double **divuw, int nrt, int nzt)
{
    extern double rright;

    int i, k;
    double r_c, r_p, r_m, drt;

    drt = rright / (double)nrt;
    ikloopt
    {
        r_c = ((double)i - 0.5) * drt;
        r_p = r_c + 0.5 * drt;
        r_m = r_c - 0.5 * drt;
     
        divuw[i][k] = (r_p * tu[i][k] - r_m * tu[i - 1][k]) / (r_c * drt) + (tw[i][k] - tw[i][k - 1]) / drt;
    }
}

void MG_Poisson(double **u, double **f)
{
    extern int nr, nz;
    extern double **rho_v;

    int i, k, max_it = 1000, it_mg = 1;
    double tol = 1.0e-8, resid = 1.0, **work_t, **sor;

    work_t = dmatrix(1, nr, 1, nz);
    sor = dmatrix(1, nr, 1, nz);
    mat_copy(work_t, u, 1, nr, 1, nz);

    while (it_mg <= max_it && resid >= tol)
    {

        vcycle(u, f, rho_v, nr, nz, 1);
     
        pressure_update(u);

        ikloop
            sor[i][k] = work_t[i][k] - u[i][k];

        resid = mat_max(sor, 1, nr, 1, nz);

        mat_copy(work_t, u, 1, nr, 1, nz);
        it_mg++;
    }

    printf("Mac iteration = %d  residual = %16.14f \n", it_mg--, resid);
    free_dmatrix(work_t, 1, nr, 1, nz);
    free_dmatrix(sor, 1, nr, 1, nz);

    return;
}

void vcycle(double **uf, double **ff, double **wf, int nrf, int nzf, int ilevel)
{
    extern int n_level;

    relax(uf, ff, wf, nrf, nzf);

    if (ilevel < n_level)
    {

        int nrc, nzc;
        double **rf, **uc, **fc, **wc;

        nrc = nrf / 2;
        nzc = nzf / 2;

        rf = dmatrix(0, nrf + 1, 0, nzf + 1);
        uc = dmatrix(0, nrc + 1, 0, nzc + 1);
        fc = dmatrix(0, nrc + 1, 0, nzc + 1);
        wc = dmatrix(0, nrc + 1, 0, nzc + 1); 

        residual_den(rf, uf, ff, wf, nrf, nzf);

        restrict1(rf, nrf, nzf, fc, nrc, nzc);
        restrict1(wf, nrf, nzf, wc, nrc, nzc);

        augc(wc, nrc, nzc);

        zero_matrix(uc, 0, nrc + 1, 0, nzc + 1);

        vcycle(uc, fc, wc, nrc, nzc, ilevel + 1);

        prolong(uc, nrc, nzc, rf, nrf, nzf);
        mat_add(uf, uf, rf, 1, nrf, 1, nzf);

        relax(uf, ff, wf, nrf, nzf);

        free_dmatrix(rf, 0, nrf + 1, 0, nzf + 1);
        free_dmatrix(uc, 0, nrc + 1, 0, nzc + 1);
        free_dmatrix(fc, 0, nrc + 1, 0, nzc + 1);
        free_dmatrix(wc, 0, nrc + 1, 0, nzc + 1);
    }
}
void relax(double **p, double **f, double **w, int nrt, int nzt)
{
    extern int p_relax;
    extern double rright, gravity;

    int i, k, iter;
    double r_c, r_p, r_m, drt, drt2, a[4], coef, src, hr = 0.0;

    drt = rright / (double)nrt;
    drt2 = pow(drt, 2);

    for (iter = 1; iter <= p_relax; iter++)
    {
        augc(p, nrt, nzt);

        ikloopt
        {
            r_c = ((double)i - 0.5) * drt;
            r_p = r_c + 0.5 * drt;
            r_m = r_c - 0.5 * drt;

            a[0] = 2.0 / (w[i + 1][k] + w[i][k]);
            a[1] = 2.0 / (w[i][k] + w[i - 1][k]);
            a[2] = 2.0 / (w[i][k + 1] + w[i][k]);
            a[3] = 2.0 / (w[i][k] + w[i][k - 1]);

            src = f[i][k];
            coef = 0.0;

            src -= (a[0] * r_p * p[i + 1][k] + a[1] * r_m * p[i - 1][k]) / (r_c * drt2);
            coef -= (a[0] * r_p + a[1] * r_m) / (r_c * drt2);

            src -= (a[2] * p[i][k + 1] + a[3] * p[i][k - 1]) / drt2;
            coef -= (a[2] + a[3]) / drt2;

            p[i][k] = src / coef;
        }
    }
}

void grad_p(double **p, double **dpdr, double **dpdz, int nrt, int nzt)
{
    extern double rright, rho1, rho2, gravity;
    extern int nz;

    int i, k;
    double drt;

    drt = rright / (double)nrt;
    augc(p, nrt, nzt);
    i0kloopt
    {
        dpdr[i][k] = (p[i + 1][k] - p[i][k]) / drt;
    }

    ik0loopt
    {
        dpdz[i][k] = (p[i][k + 1] - p[i][k]) / drt;
    }
}


void residual_den(double **r, double **u, double **f, double **den, int nrt, int nzt) // r= u-f//
{
    int i, k;
    double **dpdr, **dpdz;

    dpdr = dmatrix(0, nrt, 1, nzt);
    dpdz = dmatrix(1, nrt, 0, nzt);

    grad_p(u, dpdr, dpdz, nrt, nzt); 

    i0kloopt
    {
        dpdr[i][k] = dpdr[i][k] / (0.5 * (den[i + 1][k] + den[i][k]));
    }

    ik0loopt
    {
        dpdz[i][k] = dpdz[i][k] / (0.5 * (den[i][k + 1] + den[i][k]));
    } 

    div_uw(dpdr, dpdz, r, nrt, nzt);  
    mat_sub(r, f, r, 1, nrt, 1, nzt); 

    free_dmatrix(dpdr, 0, nrt, 1, nzt);
    free_dmatrix(dpdz, 1, nrt, 0, nzt);
}

void sf_force(double **phi, double **mu, double **fx, double **fy) 
{
    extern int nr, nz;
    extern double h;

    int i, k;

    augc(phi, nr, nz);
    augc(mu, nr, nz);

    i0kloop
    {
        
        fx[i][k] = 0.5 * (phi[i + 1][k] - phi[i][k]) * (mu[i + 1][k] + mu[i][k]) / (h);
    }

    ik0loop
    {
        
        fy[i][k] = 0.5 * (phi[i][k + 1] - phi[i][k]) * (mu[i][k + 1] + mu[i][k]) / (h);
    }
}
void surface_tension(double **phi2, double **fr, double **fz)
{
    extern int nr, nz;
    extern double dr, gam, sig, gravity, **rho_v;

    int i, k;
    double adr, adz, rplus, rminus, rmain, curvature, fac, **half_u, **dpdr, **dpdz, **phi;

    fac = 6.0 * sqrt(2.0) * gam;

    half_u = dmatrix(-1, nr + 1, -1, nz + 1);
    dpdr = dmatrix(-1, nr + 1, -1, nz + 1);
    dpdz = dmatrix(-1, nr + 1, -1, nz + 1);
    phi = dmatrix(-1, nr + 2, -1, nz + 2);
    mat_copy(phi, phi2, 1, nr + 1, 1, nz + 1);

    augc(phi, nr, nz);

    for (i = 0; i < nr + 1; i++)
        for (k = 0; k < nz + 1; k++)
        {

            dpdr[i][k] = (phi[i + 1][k + 1] + phi[i + 1][k] - phi[i][k + 1] - phi[i][k]) / (2.0 * dr);
            dpdz[i][k] = (phi[i + 1][k + 1] - phi[i + 1][k] + phi[i][k + 1] - phi[i][k]) / (2.0 * dr);
            half_u[i][k] = sqrt(pow(dpdr[i][k], 2) + pow(dpdz[i][k], 2));
        }

    ikloop
    {

        if ((phi[i][k] < 0.98) && (phi[i][k] > 0.02))
        {

            rmain = ((double)i - 0.5) * dr;
            rplus = rmain + 0.5 * dr;
            rminus = rplus - dr;
     

            curvature = (rplus * (dpdr[i][k] / half_u[i][k] + dpdr[i][k - 1] / half_u[i][k - 1]) - rminus * (dpdr[i - 1][k] / half_u[i - 1][k] + dpdr[i - 1][k - 1] / half_u[i - 1][k - 1])) / (2.0 * rmain * dr)

                        + (dpdz[i][k] / half_u[i][k] - dpdz[i][k - 1] / half_u[i][k - 1] + dpdz[i - 1][k] / half_u[i - 1][k] - dpdz[i - 1][k - 1] / half_u[i - 1][k - 1]) / (2.0 * dr);

            adr = (dpdr[i][k] + dpdr[i - 1][k] + dpdr[i][k - 1] + dpdr[i - 1][k - 1]) / 4.0;
            adz = (dpdz[i][k] + dpdz[i - 1][k] + dpdz[i][k - 1] + dpdz[i - 1][k - 1]) / 4.0;

            fr[i][k] = fac * curvature * sqrt(adr * adr + adz * adz) * adr;

            fz[i][k] = fac * curvature * sqrt(adr * adr + adz * adz) * adz;
        }
        else
        {
            fr[i][k] = 0.0;
            fz[i][k] = 0.0;
        }
    }

    ikloop
    {

        if (i == nr)
        {
            dpdr[0][k] = 0.0;
            dpdr[nr][k] = 0.0;
        }
        else
            dpdr[i][k] = ((fr[i][k] + fr[i + 1][k]) * 0.5);

        if (k == nz)
        {
            dpdz[i][0] = 0.0;
            dpdz[i][nz] = 0.0;
        }
        else
            dpdz[i][k] = ((fz[i][k] + fz[i][k + 1])) * 0.5;
    }

    mat_copy(fr, dpdr, 0, nr, 1, nz);
    mat_copy(fz, dpdz, 1, nr, 0, nz);

    free_dmatrix(half_u, -1, nr + 1, -1, nz + 1);
    free_dmatrix(dpdr, -1, nr + 1, -1, nz + 1);

    free_dmatrix(dpdz, -1, nr + 1, -1, nz + 1);

    free_dmatrix(phi, -1, nr + 2, -1, nz + 2);

    return;
}


void cahn(double **c_old, double **cc_old, double **c_new, double **mu, double theta, double **adv_phi)
{
    extern int nr, nz;
    extern double **Mo1, **Mo2;

    int it_mg = 1, max_it_CH = 500;
    double resid = 1.0, tol = 1.0e-6, **sc, **smu, **ct;
    ct = dmatrix(1, nr, 1, nz);

    sc = dmatrix(1, nr, 1, nz);
    smu = dmatrix(1, nr, 1, nz);

    mat_copy(ct, c_old, 1, nr, 1, nz);

    source_ch(sc, smu, c_old, cc_old, theta, adv_phi); 

    while (it_mg <= max_it_CH && resid > tol)
    {

        vcycle_ch(c_new, Mo1, mu, sc, smu, nr, nz, 1);

        resid = error(ct, c_new, nr, nz);
        mat_copy(ct, c_new, 1, nr, 1, nz);

        it_mg++;
    }
    printf("cahn %16.14f   %d\n", resid, it_mg - 1);
    free_dmatrix(ct, 1, nr, 1, nz);
    free_dmatrix(sc, 1, nr, 1, nz);
    free_dmatrix(smu, 1, nr, 1, nz);
}


void source_ch(double **sc, double **smu, double **c_old, double **cc_old, double theta, double **adv_c)
{
    extern int nr, nz, it;
    extern double dt, h, mobil, gam, lambda, **Mo2, SS, **normphi3, **beta;
    int i, k, zero_norm, m, n;
    double **intc, rplus, rminus, rmain;
    intc = dmatrix(0, nr + 1, 0, nz + 1);

    ikloop
    {
        intc[i][k] = 2.0 * c_old[i][k] - cc_old[i][k];
    }
    augc(intc, nr, nz);

    ikloop
    {

        if (it == 1)
        {
            sc[i][k] = c_old[i][k] / dt - adv_c[i][k];
        }
        else
        {
            sc[i][k] = (4.0 * c_old[i][k] - cc_old[i][k]) / (2.0 * dt) - adv_c[i][k];
        }

        smu[i][k] = dfphi(intc[i][k]) + beta[i][k] - SS * intc[i][k] + gam * intc[i][k] * (1 - intc[i][k]) * normphi3[i][k] * cos(theta) / sqrt(2.0);
    }

    free_dmatrix(intc, 0, nr + 1, 0, nz + 1);
}

void vcycle_ch(double **uf_new, double **Mo1, double **wf_new, double **sc, double **smu, int nrf, int nzf, int ilevel)
{
    extern int n_level, it;

    relax_ch(uf_new, Mo1, wf_new, sc, smu, ilevel, nrf, nzf);

    if (ilevel < n_level)
    {

        int nrc, nzc;
        double **dcc, **dmuc, **uc_new, **wc_new,
            **oMo,
            **uc_def, **wc_def, **uf_def, **wf_def;

        nrc = nrf / 2;
        nzc = nzf / 2;
        uc_new = dmatrix(0, nrc + 1, 0, nzc + 1);
        wc_new = dmatrix(0, nrc + 1, 0, nzc + 1);
        dcc = dmatrix(1, nrc, 1, nzc);
        dmuc = dmatrix(1, nrc, 1, nzc);

        oMo = dmatrix(0, nrc + 1, 0, nzc + 1);

        uc_def = dmatrix(0, nrc + 1, 0, nzc + 1);
        wc_def = dmatrix(0, nrc + 1, 0, nzc + 1);
        uf_def = dmatrix(1, nrf, 1, nzf);
        wf_def = dmatrix(1, nrf, 1, nzf);

        defect(dcc, dmuc, Mo1, uf_new, wf_new, sc, smu, nrf, nzf);
        restrict1(Mo1, nrf, nzf, oMo, nrc, nrc);
        augc(oMo, nrc, nrc);

        zero_matrix(uc_def, 0, nrc + 1, 0, nzc + 1);
        zero_matrix(wc_def, 0, nrc + 1, 0, nzc + 1);

        vcycle_ch(uc_def, oMo, wc_def, dcc, dmuc, nrc, nzc, ilevel + 1);

        prolong_ch(uc_def, uf_def, wc_def, wf_def, nrc, nzc);

        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nrf, 1, nzf);

        relax_ch(uf_new, Mo1, wf_new, sc, smu, ilevel, nrf, nzf);

        free_dmatrix(dcc, 1, nrc, 1, nzc);
        free_dmatrix(dmuc, 1, nrc, 1, nzc);
        free_dmatrix(uc_new, 0, nrc + 1, 0, nzc + 1);
        free_dmatrix(wc_new, 0, nrc + 1, 0, nzc + 1);
        free_dmatrix(oMo, 0, nrc + 1, 0, nrc + 1);
        free_dmatrix(uf_def, 1, nrf, 1, nzf);
        free_dmatrix(wf_def, 1, nrf, 1, nzf);
        free_dmatrix(uc_def, 0, nrc + 1, 0, nzc + 1);
        free_dmatrix(wc_def, 0, nrc + 1, 0, nzc + 1);
    }
}

void relax_ch(double **c_new, double **Mo1, double **mu_new, double **sc,
              double **smu, int ilevel, int nrt, int nzt)
{

    extern int c_relax, it;
    extern double dt, Cahn, mobil, rright;

    int i, k, iter;
    double drt, drt2, a[4], f[2], det, r_p, r_m, r_c, xfac, yfac;

    drt = rright / (double)nrt;
    drt2 = pow(drt, 2);

    for (iter = 1; iter <= c_relax; iter++)
    {
        augc(c_new, nrt, nzt);
        augc(mu_new, nrt, nzt);

        ikloopt
        {
            r_c = ((double)i - 0.5) * drt;
            r_p = r_c + 0.5 * drt;
            r_m = r_c - 0.5 * drt;

            xfac = 0.5 * (Mo1[i + 1][k] + Mo1[i][k]) * r_p + 0.5 * (Mo1[i - 1][k] + Mo1[i][k]) * r_m;

            yfac = 0.5 * (Mo1[i][k + 1] + Mo1[i][k]) + 0.5 * (Mo1[i][k - 1] + Mo1[i][k]);
            if (it == 1)
            {
                a[0] = 1.0 / dt;
            }
            else
            {
                a[0] = 3.0 / (2.0 * dt);
            }

            a[1] = (xfac / r_c + yfac) * mobil / drt2;
            a[2] = -SS - Cahn * (xfac / r_c + yfac) / drt2;
            a[3] = 1.0;

            f[0] = sc[i][k];

            f[0] += 0.5 * (Mo1[i + 1][k] + Mo1[i][k]) * r_p * mu_new[i + 1][k] * mobil / (r_c * drt2);

            f[0] += 0.5 * (Mo1[i - 1][k] + Mo1[i][k]) * r_m * mu_new[i - 1][k] * mobil / (r_c * drt2);

            f[0] += 0.5 * (Mo1[i][k + 1] + Mo1[i][k]) * mu_new[i][k + 1] * mobil / (drt2);

            f[0] += 0.5 * (Mo1[i][k - 1] + Mo1[i][k]) * mu_new[i][k - 1] * mobil / (drt2);

            f[1] = smu[i][k];

            f[1] -= 0.5 * (Mo1[i + 1][k] + Mo1[i][k]) * Cahn * r_p * c_new[i + 1][k] / (r_c * drt2);

            f[1] -= 0.5 * (Mo1[i - 1][k] + Mo1[i][k]) * Cahn * r_m * c_new[i - 1][k] / (r_c * drt2);

            f[1] -= 0.5 * (Mo1[i][k + 1] + Mo1[i][k]) * Cahn * c_new[i][k + 1] / drt2;

            f[1] -= 0.5 * (Mo1[i][k - 1] + Mo1[i][k]) * Cahn * c_new[i][k - 1] / drt2;

            det = a[0] * a[3] - a[1] * a[2];

            c_new[i][k] = (a[3] * f[0] - a[1] * f[1]) / det;
            mu_new[i][k] = (-a[2] * f[0] + a[0] * f[1]) / det;
        }
    }
}
void defect(double **duc, double **dwc, double **Mo, double **uf_new, double **wf_new,
            double **suf, double **swf, int nxf, int nyf)
{
    double **ruf, **rwf, **rruf, **rrwf;

    ruf = dmatrix(1, nxf, 1, nyf);
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxf / 2, 1, nyf / 2);
    rrwf = dmatrix(1, nxf / 2, 1, nyf / 2);

    nonL(ruf, rwf, Mo, uf_new, wf_new, nxf, nyf);

    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf);

    restrict2(ruf, duc, rwf, dwc, nxf / 2, nyf / 2);

    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxf / 2, 1, nyf / 2);
    free_dmatrix(rrwf, 1, nxf / 2, 1, nyf / 2);
}

void nonL(double **rc, double **rmu, double **Mo1, double **c_new, double **mu_new, int nrt, int nzt)
{
    extern double dt, mobil, Cahn, SS;
    extern int it;

    int i, k;
    double **c_lap, **mu_lap;

    c_lap = dmatrix(1, nrt, 1, nzt);
    mu_lap = dmatrix(1, nrt, 1, nzt);

    laplace_Mo(Mo1, c_new, c_lap, nrt, nzt);
    laplace_Mo(Mo1, mu_new, mu_lap, nrt, nzt);
    ikloopt
    {
        if (it == 1)
        {
            rc[i][k] = c_new[i][k] / dt - mu_lap[i][k] * mobil;
        }
        else
        {
            rc[i][k] = 3.0 * c_new[i][k] / (2.0 * dt) - mu_lap[i][k] * mobil;
        }

        rmu[i][k] = mu_new[i][k] + Cahn * c_lap[i][k] - SS * c_new[i][k];
    }

    free_dmatrix(c_lap, 1, nrt, 1, nzt);
    free_dmatrix(mu_lap, 1, nrt, 1, nzt);
}

void laplace_ch(double **a, double **lap_a, int nrt, int nzt)
{
    extern double rright;

    int i, k;
    double r_i, r_iphf, r_imhf, drt, dzt,
        dadr_L, dadr_R, dadz_B, dadz_T;

    drt = rright / (double)nrt;
    dzt = drt;
    augc(a, nrt, nzt);
    ikloopt
    {
        r_i = ((double)i - 0.5) * drt;
        r_iphf = r_i + 0.5 * drt;
        r_imhf = r_i - 0.5 * drt;


        dadr_L = (a[i][k] - a[i - 1][k]) / drt;

        dadr_R = (a[i + 1][k] - a[i][k]) / drt;

        dadz_B = (a[i][k] - a[i][k - 1]) / dzt;

        dadz_T = (a[i][k + 1] - a[i][k]) / dzt;

        lap_a[i][k] = (r_iphf * dadr_R - r_imhf * dadr_L) / (r_i * drt) + (dadz_T - dadz_B) / dzt;
    }
}
void laplace_Mo(double **Mo, double **a, double **lap_a, int nrt, int nzt)
{
    extern double rright;

    int i, k;
    double r_i, r_iphf, r_imhf, drt, dzt,
        dadr_L, dadr_R, dadz_B, dadz_T;

    drt = rright / (double)nrt;
    dzt = drt;
    augc(a, nrt, nzt);
    ikloopt
    {
        r_i = ((double)i - 0.5) * drt;
        r_iphf = r_i + 0.5 * drt;
        r_imhf = r_i - 0.5 * drt;
 
        dadr_L = 0.5 * (Mo[i][k] + Mo[i - 1][k]) * (a[i][k] - a[i - 1][k]) / drt;

        dadr_R = 0.5 * (Mo[i][k] + Mo[i + 1][k]) * (a[i + 1][k] - a[i][k]) / drt;

        dadz_B = 0.5 * (Mo[i][k] + Mo[i][k - 1]) * (a[i][k] - a[i][k - 1]) / dzt;

        dadz_T = 0.5 * (Mo[i][k] + Mo[i][k + 1]) * (a[i][k + 1] - a[i][k]) / dzt;

        lap_a[i][k] = (r_iphf * dadr_R - r_imhf * dadr_L) / (r_i * drt) + (dadz_T - dadz_B) / dzt;
    }
}

void restrict1(double **cf, int nrf, int nzf,
               double **cc, int nrc, int nzc)
{
    extern double rright;

    int ir, iz;
    double drc, drf, rc_i, rf_2im1, rf_2i;

    drc = rright / (double)nrc;
    drf = 0.5 * drc;

    for (ir = 1; ir <= nrc; ir++)
    {

        rc_i = ((double)ir - 0.5) * drc;
        rf_2im1 = ((double)(2 * ir - 1) - 0.5) * drf;
        rf_2i = ((double)(2 * ir) - 0.5) * drf;


        for (iz = 1; iz <= nzc; iz++)
        {

            cc[ir][iz] = (rf_2im1 * (cf[2 * ir - 1][2 * iz - 1] + cf[2 * ir - 1][2 * iz]) + rf_2i * (cf[2 * ir][2 * iz - 1] + cf[2 * ir][2 * iz])) / (4.0 * rc_i);
        }
    }
}

void prolong(double **cc, int nrc, int nzc,
             double **cf, int nrf, int nzf)
{
    extern double rright;

    int ir, iz;
    double drc, drf, rc_i, rf_2im1, rf_2i;

    drc = rright / (double)nrc;
    drf = 0.5 * drc;

    for (ir = 1; ir <= nrc; ir++)
    {

        rc_i = ((double)ir - 0.5) * drc;
        rf_2im1 = ((double)(2 * ir - 1) - 0.5) * drf;
        rf_2i = ((double)(2 * ir) - 0.5) * drf;

        for (iz = 1; iz <= nzc; iz++)
        {

            cf[2 * ir - 1][2 * iz - 1] = rc_i / rf_2im1 * cc[ir][iz];
            cf[2 * ir - 1][2 * iz] = rc_i / rf_2im1 * cc[ir][iz];
            cf[2 * ir][2 * iz - 1] = rc_i / rf_2i * cc[ir][iz];
            cf[2 * ir][2 * iz] = rc_i / rf_2i * cc[ir][iz];
        }
    }
}

void restrict2(double **cf, double **cc, double **df, double **dc,
               int nrc, int nzc)
{
    extern double rright;

    int ir, iz;
    double drc, drf, rc_i, rf_2im1, rf_2i;

    drc = rright / (double)nrc;
    drf = 0.5 * drc;

    for (ir = 1; ir <= nrc; ir++)
    {

        rc_i = ((double)ir - 0.5) * drc;
        rf_2im1 = ((double)(2 * ir - 1) - 0.5) * drf;
        rf_2i = ((double)(2 * ir) - 0.5) * drf;


        for (iz = 1; iz <= nzc; iz++)
        {

            cc[ir][iz] = (rf_2im1 * (cf[2 * ir - 1][2 * iz - 1] + cf[2 * ir - 1][2 * iz]) + rf_2i * (cf[2 * ir][2 * iz - 1] + cf[2 * ir][2 * iz])) / (4.0 * rc_i);

            dc[ir][iz] = (rf_2im1 * (df[2 * ir - 1][2 * iz - 1] + df[2 * ir - 1][2 * iz]) + rf_2i * (df[2 * ir][2 * iz - 1] + df[2 * ir][2 * iz])) / (4.0 * rc_i);
        }
    }
}

void prolong_ch(double **cc, double **cf, double **dc, double **df,
                int nrc, int nzc)
{
    extern double rright;

    int ir, iz;
    double drc, drf, rc_i, rf_2im1, rf_2i;

    drc = rright / (double)nrc;
    drf = 0.5 * drc;

    for (ir = 1; ir <= nrc; ir++)
    {

        rc_i = ((double)ir - 0.5) * drc;
        rf_2im1 = ((double)(2 * ir - 1) - 0.5) * drf;
        rf_2i = ((double)(2 * ir) - 0.5) * drf;


        for (iz = 1; iz <= nzc; iz++)
        {

            cf[2 * ir - 1][2 * iz - 1] = rc_i / rf_2im1 * cc[ir][iz];
            cf[2 * ir - 1][2 * iz] = rc_i / rf_2im1 * cc[ir][iz];
            cf[2 * ir][2 * iz - 1] = rc_i / rf_2i * cc[ir][iz];
            cf[2 * ir][2 * iz] = rc_i / rf_2i * cc[ir][iz];

            df[2 * ir - 1][2 * iz - 1] = rc_i / rf_2im1 * dc[ir][iz];
            df[2 * ir - 1][2 * iz] = rc_i / rf_2im1 * dc[ir][iz];
            df[2 * ir][2 * iz - 1] = rc_i / rf_2i * dc[ir][iz];
            df[2 * ir][2 * iz] = rc_i / rf_2i * dc[ir][iz];
        }
    }
}

double error(double **c_old, double **c_new, int nrt, int nzt)
{
    double **uc, value;

    uc = dmatrix(1, nrt, 1, nzt);
    mat_sub(uc, c_new, c_old, 1, nrt, 1, nzt);

    value = mat_max(uc, 1, nrt, 1, nzt);
    free_dmatrix(uc, 1, nrt, 1, nzt);

    return value;
}



double *dvector(long nl, long nh)

{
    double *v;

    v = (double *)malloc((nh - nl + 1 + NR_END) * sizeof(double));
    memset(v, 0, (nh - nl + 1 + NR_END) * sizeof(double *));

    return v - nl + NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    double **m;
    long i, nrow = nrh - nrl + 1 + NR_END, ncol = nch - ncl + 1 + NR_END;

    m = (double **)malloc((nrow) * sizeof(double *));
    memset(m, 0, (nrow) * sizeof(double *));
    m += NR_END;
    m -= nrl;

    m[nrl] = (double *)malloc((nrow * ncol) * sizeof(double));
    memset(m[nrl], 0, (nrow * ncol) * sizeof(double));
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
    {
        m[i] = m[i - 1] + ncol;
    }

    return m;
}

void free_dvector(double *v, long nl, long nh)
{
    free(v + nl - NR_END);

    return;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free(m[nrl] + ncl - NR_END);
    free(m + nrl - NR_END);

    return;
}

void zero_vector(double *a, int nl, int nh)
{
    int i;

    for (i = nl; i <= nh; i++)
        a[i] = 0.0;

    return;
}

void veccopy(double *a, double *b,
             int nla, int nha,
             int nlb, int nhb)

{
    int i, idiff = nlb - nla;

    for (i = nla; i <= nha; i++)
        a[i] = b[i + idiff];

    return;
}

void vecadd(double *a, double *b, double *phi1, int nl, int nr)
{
    int i;

    for (i = nl; i <= nr; i++)
        a[i] = b[i] + phi1[i];

    return;
}

void vecsub(double *a, double *b, double *phi1, int nl, int nr)
{
    int i;

    for (i = nl; i <= nr; i++)
        a[i] = b[i] - phi1[i];

    return;
}

void s_vec_add(double *a, double *b, double s,
               int nl, int nr)
{
    int i;

    for (i = nl; i <= nr; i++)
        a[i] = b[i] + s;

    return;
}

void s_vec_mult(double *a, double fac, double *b,
                int nl, int nr)
{
    int i;

    for (i = nl; i <= nr; i++)
        a[i] = fac * b[i];

    return;
}

void veccomb(double *a, double s1, double *b, double s2, double *phi1,
             int nl, int nr)
{
    int i;

    for (i = nl; i <= nr; i++)
        a[i] = s1 * b[i] + s2 * phi1[i];

    return;
}

void vecmult(double *a, double *b, double *phi1, int nl, int nr)
{
    int i;

    for (i = nl; i <= nr; i++)
        a[i] = b[i] * phi1[i];
}

void vecdiv(double *a, double *b, double *phi1, int nl, int nr)
{
    int i;

    for (i = nl; i <= nr; i++)
        a[i] = b[i] / phi1[i];
}

double sum_vector(double *a, int nl, int nh)
{
    int i;
    double sum = 0.0;

    for (i = nl; i <= nh; i++)
        sum += a[i];

    return sum;
}

double dot_product(double *a, double *b, int nl, int nr)
{
    int i;
    double x = 0.0;

    for (i = nl; i <= nr; i++)
        x += a[i] * b[i];

    return x;
}


void mat_add(double **a, double **b, double **phi1,
             int xl, int xr, int yl, int yr)
{
    int i, k;

    for (i = xl; i <= xr; i++)
        for (k = yl; k <= yr; k++)
        {
            a[i][k] = b[i][k] + phi1[i][k];
        }

    return;
}

void mat_add2(double **a, double **b, double **phi1,
              double **a2, double **b2, double **c2,
              int xl, int xr, int yl, int yr)
{
    int i, k;

    for (i = xl; i <= xr; i++)
        for (k = yl; k <= yr; k++)
        {
            a[i][k] = b[i][k] + phi1[i][k];
            a2[i][k] = b2[i][k] + c2[i][k];
        }

    return;
}

void zero_matrix(double **a, int xl, int xr, int yl, int yr)
{
    int i, k;

    for (i = xl; i <= xr; i++)
        for (k = yl; k <= yr; k++)
        {

            a[i][k] = 0.0;
        }

    return;
}

void mat_copy(double **a, double **b,
              int xl, int xr, int yl, int yr)

{
    int i, k;

    for (i = xl; i <= xr; i++)
        for (k = yl; k <= yr; k++)

            a[i][k] = b[i][k];

    return;
}

void mat_copy2(double **a, double **b,
               double **a2, double **b2,
               int xl, int xr, int yl, int yr)

{
    int i, k;

    for (i = xl; i <= xr; i++)
        for (k = yl; k <= yr; k++)
        {

            a[i][k] = b[i][k];
            a2[i][k] = b2[i][k];
        }
    return;
}

void mat_sub(double **a, double **b, double **phi1,
             int nrl, int nrh, int ncl, int nch)
{
    int i, k;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
            a[i][k] = b[i][k] - phi1[i][k];

    return;
}

void mat_sub2(double **a, double **b, double **phi1,
              double **a2, double **b2, double **c2,
              int nrl, int nrh, int ncl, int nch)
{
    int i, k;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
        {
            a[i][k] = b[i][k] - phi1[i][k];
            a2[i][k] = b2[i][k] - c2[i][k];
        }

    return;
}

void matmult(double **a, double **b, double **phi1,
             int nrl, int nrh, int ncl, int nch)
{
    int i, k;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
            a[i][k] = b[i][k] * phi1[i][k];
}

void mat_comb(double **a, double s1, double **b,
              double s2, double **phi1,
              int nrl, int nrh, int ncl, int nch)
{
    int i, k;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
            a[i][k] = s1 * b[i][k] + s2 * phi1[i][k];

    return;
}

void mat_comb3(double **a, double s1, double **b, double s2, double **phi1,
               double s3, double **d,
               int nrl, int nrh, int ncl, int nch)
{
    int i, k;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
            a[i][k] = s1 * b[i][k] + s2 * phi1[i][k] + s3 * d[i][k];

    return;
}

void s_mat_mult(double **a, double fac, double **b,
                int nrl, int nrh, int ncl, int nch)

{
    int i, k;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
            a[i][k] = fac * b[i][k];

    return;
}

int s_mat_mult2(double **a, int arl, int arh, int acl, int ach,
                double fac,
                double **b, int brl, int brh, int bcl, int bch)
{
    int i, k;
    int abr, abc;

    if (((arh - arl) != (brh - brl)) ||
        ((ach - acl) != (bch - bcl)))
        return 1;

    abr = brl - arl;
    abc = bcl - acl;

    for (i = arl; i <= arh; i++)
        for (k = acl; k <= ach; k++)
            a[i][k] = fac * b[i + abr][k + abc];

    return 0;
}


double sum_matrix(double **a,
                  int nrl, int nrh, int ncl, int nch)
{
    int i, k;
    double sum = 0.0;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
            sum += a[i][k];

    return sum;
}


double max_norm(double **a, double **b,
                int nrl, int nrh, int ncl, int nch)
{
    int i, k;
    double x = 0.0;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
        {

            if (fabs(a[i][k] - b[i][k]) > x)
                x = fabs(a[i][k] - b[i][k]);
        }

    return x;
}

double mat_max(double **a,
               int nrl, int nrh, int ncl, int nch)
{
    int i, k;
    double x = 0.0;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
        {

            if (fabs(a[i][k]) > x)
                x = fabs(a[i][k]);
        }

    return x;
}

double mat_normal2(double **a,
                   int nrl, int nrh, int ncl, int nch)
{
    int i, k;
    double x = 0.0;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
        {

            x += pow(a[i][k], 2);
        }

    return (sqrt(x) / ((nrh - nrl + 1) * (nch - ncl + 1)));
}

void print_mat(FILE *ofptr, double **a, int xl, int xr, int yl, int yr)
{

    int i, k;
    for (i = xl; i <= xr; i++)
    {
        for (k = yl; k <= yr; k++)
        {

            fprintf(ofptr, " %f ", a[i][k]);
        }
        fprintf(ofptr, "\n ");
    }
    return;
}

void print_data(double **u, double **v, double **phi1, double **phi2, double **phi3, double **p)
{
    extern int nr, nz;

    char bufferphi1[2000], bufferphi2[2000], bufferphi3[2000], bufferu[2000], bufferv[2000], bufferp[2000];

    int i, k;

    FILE *fu, *fv, *fp, *fphi1, *fphi2, *fphi3;

    sprintf(bufferphi1, "E: \\phi1.m");
    sprintf(bufferphi2, "E: \\phi2.m");
    sprintf(bufferphi3, "E: \\phi3.m");
    sprintf(bufferu, "E: \\u.m");
    sprintf(bufferv, "E: \\w.m");
    sprintf(bufferp, "E: \\p.m");

    fphi1 = fopen(bufferphi1, "a");
    fphi2 = fopen(bufferphi2, "a");
    fphi3 = fopen(bufferphi3, "a");
    fu = fopen(bufferu, "a");
    fv = fopen(bufferv, "a");
    fp = fopen(bufferp, "a");

    for (i = 1; i <= nr; i++)
    {
        for (k = 1; k <= nz; k++)
        {

            fprintf(fu, " %f ", 0.5 * (u[i - 1][k] + u[i][k]));
            fprintf(fv, " %f ", 0.5 * (v[i][k] + v[i][k - 1]));

        }
        fprintf(fu, "\n ");
        fprintf(fv, "\n ");
    }

    print_mat(fphi1, phi1, 1, nr, 1, nz);
    print_mat(fphi2, phi2, 1, nr, 1, nz);
    print_mat(fphi3, phi3, 1, nr, 1, nz);
    print_mat(fp, p, 1, nr, 1, nz);

    fclose(fu);
    fclose(fv);
    fclose(fp);
    fclose(fphi1);
    fclose(fphi2);
    fclose(fphi3);

    return;
}

void pressure_update(double **a)
{

    extern int nr, nz;

    int i, k;
    double ave = 0.0;

    for (i = 1; i <= nr; i++)
        for (k = 1; k <= nz; k++)
        {

            ave = ave + a[i][k];
        }

    ave /= (nr + 0.0) * (nz + 0.0);

    for (i = 1; i <= nr; i++)
        for (k = 1; k <= nz; k++)
        {

            a[i][k] -= ave;
        }

    return;
}


double mass_comp(double **phi)
{
    extern int nr, nz;
    extern double h;
    int i, k;
    double sum1 = 0.0, r_c, sum2 = 0.0;

    for (k = 1; k <= nz; k++)
        for (i = 1; i <= nr; i++)
        {
            r_c = 1.0 * (i - 0.5) * h;

            sum1 += r_c * phi[i][k] * h * h;
            sum2 += r_c * h * h;
        }
    return (sum1 / sum2);
}

double dfphi(double phi1)
{
    double result;

    result = pow(phi1, 3) - 1.5 * pow(phi1, 2) + 0.5 * phi1;

    return result;
}

void functiong(double **u, double **w, double **result)
{
    extern int nr, nz;
    extern double h, Re, **rho_v, **vis, **fr, **fz, Fr, gravity, **adv_u, **adv_w, **phi3, kappa;

    int i, k;
    double r_c, r_p, r_m;

    auguw(u, w, nr, nz);
    i0kloop
    {
        r_c = ((double)i) * h;
        r_p = r_c + 0.5 * h;
        r_m = r_c - 0.5 * h;

        result[i][k] = (-adv_u[i][k] + fr[i][k] / (0.5 * (rho_v[i + 1][k] + rho_v[i][k])) + ((2.0 * r_p * vis[i + 1][k] * (u[i + 1][k] - u[i][k]) - 2.0 * r_m * vis[i][k] * (u[i][k] - u[i - 1][k])) / r_c - (vis[i + 1][k] + vis[i][k]) * u[i][k] * h * h / (r_c * r_c)

                                                                                             + 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (u[i][k + 1] - u[i][k]) - 0.25 * (vis[i][k - 1] + vis[i + 1][k - 1] + vis[i][k] + vis[i + 1][k]) * (u[i][k] - u[i][k - 1])

                                                                                             + 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (w[i + 1][k] - w[i][k]) - 0.25 * (vis[i][k - 1] + vis[i + 1][k - 1] + vis[i][k] + vis[i + 1][k]) * (w[i + 1][k - 1] - w[i][k - 1])) /
                                                                                                (h * h * Re * 0.5 * (rho_v[i + 1][k] + rho_v[i][k]))) /
                       (1.0 + dt * 0.5 * (phi3[i + 1][k] + phi3[i][k]) / (0.5 * (rho_v[i + 1][k] + rho_v[i][k]) * kappa));
    }
}
void functionq(double **u, double **w, double **result)
{
    extern int nr, nz;
    extern double h, Re, **rho_v, **vis, **fr, **fz, Fr, gravity, **adv_u, **adv_w, dt;

    int i, k;
    double r_c, r_p, r_m;

    auguw(u, w, nr, nz);
    ik0loop
    {
        r_c = ((double)i - 0.5) * h;
        r_p = r_c + 0.5 * h;
        r_m = r_c - 0.5 * h;

        result[i][k] = (-adv_w[i][k] + fz[i][k] / (0.5 * (rho_v[i][k + 1] + rho_v[i][k])) - 1 / Fr * gravity + ((r_p * 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (w[i + 1][k] - w[i][k]) - r_m * 0.25 * (vis[i - 1][k] + vis[i][k] + vis[i - 1][k + 1] + vis[i][k + 1]) * (w[i][k] - w[i - 1][k])) / r_c

                                                                                                                + 2.0 * vis[i][k + 1] * (w[i][k + 1] - w[i][k]) - 2.0 * vis[i][k] * (w[i][k] - w[i][k - 1])

                                                                                                                + (r_p * 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (u[i][k + 1] - u[i][k]) - r_m * 0.25 * (vis[i - 1][k] + vis[i][k] + vis[i - 1][k + 1] + vis[i][k + 1]) * (u[i - 1][k + 1] - u[i - 1][k])) / r_c) /
                                                                                                                   (h * h * Re * 0.5 * (rho_v[i][k + 1] + rho_v[i][k]))) /
                       (1.0 + dt * 0.5 * (phi3[i][k + 1] + phi3[i][k]) / (0.5 * (rho_v[i][k + 1] + rho_v[i][k]) * kappa));
    }
}
