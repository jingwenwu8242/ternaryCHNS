

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "main1.h"
#include "mainutil1.h"
#ifndef OUTPUT_DIR
#define OUTPUT_DIR "data6b/"
#endif
#define NR_END 1
int nr, nz, p_relax, c_relax, it;
double pi, dr, rright, h, dt, gam, Cahn, SS, **beta, **Mo1, **normphi3, **phi3, kappa,
    **worku, **workw, **workp, **adv_u, **adv_w, **adv_phi1, **adv_phi2, **muphi1, **muphi2, **fr, **fz, **frphi1, **fzphi1,
    vis1, vis2, vis3, mobil, sig, **vis, rho1, rho2, rho3, velocity, **rho_v, gravity, Fr, Re, We, **intphi1, **intphi2,
    theta, theta1, theta2, **midphi1;

int main()
{
    extern int nr, nz, p_relax, c_relax, it;
    extern double pi, dr, rright, h, dt, gam, Cahn, SS, **beta, **Mo1, **normphi3, **phi3, kappa,
        **worku, **workw, **workp, **adv_u, **adv_w, **adv_phi1, **adv_phi2, **muphi1, **muphi2, **fr, **fz, **frphi1, **fzphi1,
        vis1, vis2, vis3, mobil, sig, **vis, rho1, rho2, rho3, velocity, **rho_v, gravity, Fr, Re, We, **intphi1, **intphi2,
        theta, theta1, theta2, **midphi1;

    int max_it, ns, i, k, count = 0;
    double **u, **w, **p, **ophi1, **phi1, **nu, **nw, **nphi1, **ophi2, **phi2, **nphi2, mass1, Pe;

    clock_t start, end;
    double elapsed;
    start = clock();
    FILE *fu, *fw, *fp, *fphi1, *fphi2, *fphi3, *my, *mymass1;

    pi = 4.0 * atan(1.0);

    p_relax = 20;
    c_relax = 5;

    nr = gnr;
    nz = gnz;

    rright = 1;
    dr = rright / (double)nr;
    h = dr;
    dt = 0.000001;
    max_it = 14000;
    ns = (int)(max_it / 100 + 0.001);
    gam = 4.0 * h / (4.0 * sqrt(2.0) * atanh(0.9));

    Cahn = pow(gam, 2);

    rho1 = 1;
    rho2 = 1;
    rho3 = 100;
    vis1 = 1;
    vis2 = 1;
    vis3 = 100.0;
    Re = 1.0;
    We = 1;
    Pe = 1.0 / 100;
    mobil = 1/Pe;
    sig = 1.0;
    SS = 2.0;

    gravity = 0;
    Fr = 1;
    velocity = -0;
    theta = 60;
    theta1 = (180 - theta) * pi / 180;
    theta2 = theta * pi / 180.0;
    kappa = 1e-8;

    printf("nr = %d, nz = %d\n", nr, nz);
    printf("dt      = %f\n", dt);
    printf("max_it  = %d\n", max_it);
    printf("ns      = %d\n", ns);
    printf("sig      = %f\n", sig);
    printf("Cahn      = %f\n\n", Cahn);

    my = fopen(OUTPUT_DIR "remarks.m", "w");
    fprintf(my, "nr      = %d\n", nr);
    fprintf(my, "nz      = %d\n", nz);
    fprintf(my, "max_it  = %d\n", max_it);
    fprintf(my, "ns      = %d\n", ns);
    fprintf(my, "print times      = %d\n", max_it / ns + 1);
    fprintf(my, "dt      = %f\n", dt);
    fprintf(my, "h       = %f\n", h);
    fprintf(my, "vis1       = %f\n", vis1);
    fprintf(my, "vis2       = %f\n", vis2);
    fprintf(my, "rho1       = %f\n", rho1);
    fprintf(my, "rho2       = %f\n", rho2);
    fprintf(my, "epsilon_%d       = %f\n", (int)(gam * (4 * sqrt(2) * atanh(0.9)) / h + 0.5), gam);
    fprintf(my, "mobil       = %f\n", mobil);
    fprintf(my, "sig      = %f\n", sig);
    fprintf(my, "theta      = %f\n", theta);
    fclose(my);

    ophi1 = dmatrix(0, nr + 1, 0, nz + 1);
    phi1 = dmatrix(0, nr + 1, 0, nz + 1);
    nphi1 = dmatrix(0, nr + 1, 0, nz + 1);
    ophi2 = dmatrix(0, nr + 1, 0, nz + 1);
    phi2 = dmatrix(0, nr + 1, 0, nz + 1);
    nphi2 = dmatrix(0, nr + 1, 0, nz + 1);

    phi3 = dmatrix(0, nr + 1, 0, nz + 1);
    Mo1 = dmatrix(0, nr + 1, 0, nz + 1);
    beta = dmatrix(1, nr, 1, nz);
    normphi3 = dmatrix(1, nr, 1, nz);

    muphi1 = dmatrix(0, nr + 1, 0, nz + 1);
    muphi2 = dmatrix(0, nr + 1, 0, nz + 1);

    intphi1 = dmatrix(0, nr + 1, 0, nz + 1);
    intphi2 = dmatrix(0, nr + 1, 0, nz + 1);
    midphi1 = dmatrix(0, nr + 1, 0, nz + 1);

    u = dmatrix(-1, nr + 1, 0, nz + 1);
    nu = dmatrix(-1, nr + 1, 0, nz + 1);
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

    worku = dmatrix(0, nr, 1, nz);
    workw = dmatrix(1, nr, 0, nz);
    workp = dmatrix(1, nr, 1, nz);

    vis = dmatrix(0, nr + 1, 0, nz + 1);
    rho_v = dmatrix(0, nr + 1, 0, nz + 1);

    initialization(u, w, p, phi1, phi2, phi3);

    fu = fopen(OUTPUT_DIR "u.m", "w");
    fw = fopen(OUTPUT_DIR "w.m", "w");
    fp = fopen(OUTPUT_DIR "p.m", "w");
    fphi1 = fopen(OUTPUT_DIR "phi1.m", "w");
    fphi2 = fopen(OUTPUT_DIR "phi2.m", "w");
    fphi3 = fopen(OUTPUT_DIR "phi3.m", "w");
    fclose(fu);
    fclose(fw);
    fclose(fp);
    fclose(fphi1);
    fclose(fphi2);
    fclose(fphi3);
    mymass1 = fopen(OUTPUT_DIR "mass1.m", "w");

    mass1 = mass_comp(phi1);

    fprintf(mymass1, " %f\n", mass1);

    fclose(mymass1);

    print_data(u, w, phi1, phi2, phi3, p);

    mat_copy(nphi1, phi1, 1, nr, 1, nz);
    mat_copy(nphi2, phi2, 1, nr, 1, nz);
    mat_copy(ophi1, phi1, 1, nr, 1, nz);
    mat_copy(ophi2, phi2, 1, nr, 1, nz);
    mat_copy(nu, u, 0, nr, 1, nz);
    mat_copy(nw, w, 1, nr, 0, nz);

    augc(phi3, nr, nz);
    ikloop
    {

        Mo1[i][k] = 1.0 - phi3[i][k];
        normphi3[i][k] = sqrt(pow((phi3[i + 1][k] - phi3[i - 1][k]) / (2 * h), 2) + pow((phi3[i][k + 1] - phi3[i][k - 1]) / (2 * h), 2));
    }
    augc(Mo1, nr, nz);

    for (it = 1; it <= max_it; it++)
    {

        ikloop
        {

            intphi1[i][k] = 2.0 * phi1[i][k] - ophi1[i][k];
            intphi2[i][k] = 2.0 * phi2[i][k] - ophi2[i][k];

            midphi1[i][k] = 1.5 * phi1[i][k] - 0.5 * ophi1[i][k];
        }
        augc(intphi1, nr, nz);
        augc(intphi2, nr, nz);
        augc(midphi1, nr, nz);

        cal_vicosity(intphi1, intphi2, phi3, vis1, vis2, vis3, vis);
        cal_density(intphi1, intphi2, phi3, rho1, rho2, rho3, rho_v);
        augc(vis, nr, nz);
        augc(rho_v, nr, nz);

         surface_tension(midphi1, frphi1, fzphi1);

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
        mat_copy(phi1, nphi1, 1, nr, 1, nz);
        mat_copy(phi2, nphi2, 1, nr, 1, nz);

        mat_copy(u, nu, 0, nr, 1, nz);
        mat_copy(w, nw, 1, nr, 0, nz);

        printf("It=%d  mas=%f \n", it, mass_comp(phi1));

        if (it % ns == 0)
        {
            print_data(u, w, phi1, phi2, phi3, p);
            count++;
            printf("print out counts %d \n", count);

            end = clock();
            elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
            my = fopen(OUTPUT_DIR "remarks.m", "a");
            fprintf(my, "\n counts=    %d   iteration=    %d   CPU Time=   %f\n", count, it, elapsed);
            fclose(my);
            mass1 = mass_comp(phi1);

            mymass1 = fopen(OUTPUT_DIR "mass1.m", "a");

            fprintf(mymass1, " %f\n", mass1);

            fclose(mymass1);

        }
    }

    end = clock();
    elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time elapsed %f\n", elapsed);
    return 0;
}

void initialization(double **u, double **w, double **p, double **phi1, double **phi2, double **phi3)
{
    extern int nr, nz;
    extern double h, gam;

    int i, k;
    double r, z, psi, phi0;
    ikloop
    {

        r = ((double)i - 0.5) * h;
        z = ((double)k - 0.5) * h;

         phi0 = 0.5 + 0.5 * tanh((0.9 - sqrt(pow((r - 0), 2) + pow((z - 1.0), 2))) / (0.5 * 2 * sqrt(2.0) * gam));
        phi3[i][k] = (1 - phi0) * (0.5 + 0.5 * tanh(-fmax(fabs(r - 0) - 0.8, fabs(z + 0) - 0.25) / (0.5 * 2 * sqrt(2.0) * gam)));
        psi = (0.5 + 0.5 * tanh((0.20 - sqrt(pow((r - 0), 2) + pow((z - 0.25), 2))) / (2 * sqrt(2.0) * gam)));

        if (phi3[i][k] + psi > 1)
        {
            phi1[i][k] = 1 - phi3[i][k];
        }
        else
        {
            phi1[i][k] = psi;
        }

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
    (void)phi1;
    (void)phi2;
    (void)phi3;
    (void)den1;
    (void)den2;
    (void)den3;
    ikloop
    {
        var_den[i][k] = 1;
    }
}

void cal_vicosity(double **phi1, double **phi2, double **phi3, double den1, double den2, double den3, double **var_den)
{
    int i, k;
    (void)phi1;
    (void)phi2;
    (void)phi3;
    (void)den1;
    (void)den2;
    (void)den3;
    ikloop
    {
        var_den[i][k] = 1;
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

void full_step(double **tu, double **tw, double **p,
               double **u, double **w)
{
    extern double **adv_u, **adv_w;

    double gamma =
        1.0 - 1.0 / sqrt(2.0);

    advection(adv_u, adv_w, u, w);

    temp_uw(tu, tw, u, w);

    auguw(tu, tw, nr, nz);

    Poisson_stage(
        tu,
        tw,
        p,
        gamma * dt);

    auguw(tu, tw, nr, nz);
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

void advection_c(double **u, double **w, double **phi1, double **adv_phi1)
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

void temp_uw(double **tu, double **tw,
             double **u, double **w)
{
    extern int nr, nz;
    extern double dt, **phi3, kappa;
    extern double **adv_u, **adv_w;
    extern double **rho_v;

    int i, k;

    double gamma =
        1.0 - 1.0 / sqrt(2.0);

    double delta =
        1.0 - 1.0 / (2.0 * gamma);

    double phi_face;
    double rho_face;
    double lambda_face;
    double denom;

    double **temp1_u;
    double **temp1_w;

    double **g0;
    double **q0;

    double **g1;
    double **q1;

    double **p1;

    temp1_u =
        dmatrix(-1, nr + 1, 0, nz + 1);

    temp1_w =
        dmatrix(0, nr + 1, -1, nz + 1);

    g0 =
        dmatrix(-1, nr + 1, 0, nz + 1);

    q0 =
        dmatrix(0, nr + 1, -1, nz + 1);

    g1 =
        dmatrix(-1, nr + 1, 0, nz + 1);

    q1 =
        dmatrix(0, nr + 1, -1, nz + 1);

    p1 =
        dmatrix(0, nr + 1, 0, nz + 1);

    zero_matrix(
        p1,
        0, nr + 1,
        0, nz + 1);

    functiong(u, w, g0);
    functionq(u, w, q0);

    i0kloop
    {
        phi_face =
            0.5 *
            (phi3[i + 1][k] +
             phi3[i][k]);

        rho_face =
            0.5 *
            (rho_v[i + 1][k] +
             rho_v[i][k]);

        lambda_face =
            phi_face /
            (rho_face * kappa);

        denom =
            1.0 +
            gamma * dt * lambda_face;

        temp1_u[i][k] =
            (u[i][k] +
             gamma * dt * g0[i][k]) /
            denom;
    }

    ik0loop
    {
        phi_face =
            0.5 *
            (phi3[i][k + 1] +
             phi3[i][k]);

        rho_face =
            0.5 *
            (rho_v[i][k + 1] +
             rho_v[i][k]);

        lambda_face =
            phi_face /
            (rho_face * kappa);

        denom =
            1.0 +
            gamma * dt * lambda_face;

        temp1_w[i][k] =
            (w[i][k] +
             gamma * dt * q0[i][k]) /
            denom;
    }

    auguw(
        temp1_u,
        temp1_w,
        nr,
        nz);

    Poisson_stage(
        temp1_u,
        temp1_w,
        p1,
        gamma * dt);

    auguw(
        temp1_u,
        temp1_w,
        nr,
        nz);

    advection(
        adv_u,
        adv_w,
        temp1_u,
        temp1_w);

    functiong(
        temp1_u,
        temp1_w,
        g1);

    functionq(
        temp1_u,
        temp1_w,
        q1);

    i0kloop
    {
        phi_face =
            0.5 *
            (phi3[i + 1][k] +
             phi3[i][k]);

        rho_face =
            0.5 *
            (rho_v[i + 1][k] +
             rho_v[i][k]);

        lambda_face =
            phi_face /
            (rho_face * kappa);

        denom =
            1.0 +
            gamma * dt * lambda_face;

        tu[i][k] =
            (u[i][k]

             +

             dt *
                 (delta * g0[i][k]

                  +

                  (1.0 - delta) * g1[i][k])

             -

             (1.0 - gamma) / gamma *
                 (u[i][k]

                  +

                  gamma * dt * g0[i][k]

                  -

                  temp1_u[i][k])) /
            denom;
    }

    ik0loop
    {
        phi_face =
            0.5 *
            (phi3[i][k + 1] +
             phi3[i][k]);

        rho_face =
            0.5 *
            (rho_v[i][k + 1] +
             rho_v[i][k]);

        lambda_face =
            phi_face /
            (rho_face * kappa);

        denom =
            1.0 +
            gamma * dt * lambda_face;

        tw[i][k] =
            (w[i][k]

             +

             dt *
                 (delta * q0[i][k]

                  +

                  (1.0 - delta) * q1[i][k])

             -

             (1.0 - gamma) / gamma *
                 (w[i][k]

                  +

                  gamma * dt * q0[i][k]

                  -

                  temp1_w[i][k])) /
            denom;
    }

    free_dmatrix(
        temp1_u,
        -1, nr + 1,
        0, nz + 1);

    free_dmatrix(
        temp1_w,
        0, nr + 1,
        -1, nz + 1);

    free_dmatrix(
        g0,
        -1, nr + 1,
        0, nz + 1);

    free_dmatrix(
        q0,
        0, nr + 1,
        -1, nz + 1);

    free_dmatrix(
        g1,
        -1, nr + 1,
        0, nz + 1);

    free_dmatrix(
        q1,
        0, nr + 1,
        -1, nz + 1);

    free_dmatrix(
        p1,
        0, nr + 1,
        0, nz + 1);
}

void source_stage(double **tu,
                  double **tw,
                  double **divuw,
                  int nrt,
                  int nzt,
                  double alpha)
{
    int i, k;

    div_uw(
        tu,
        tw,
        divuw,
        nrt,
        nzt);

    ikloopt
    {
        divuw[i][k] =
            divuw[i][k] / alpha;
    }
}

void Poisson_stage(double **tu,
                   double **tw,
                   double **p,
                   double alpha)
{
    extern int nr, nz;

    extern double **workp;
    extern double **worku;
    extern double **workw;

    extern double **rho_v;
    extern double **phi3;

    extern double kappa;

    int i, k;

    double rho_face;
    double phi_face;
    double rho_eff_face;

    double **rho_eff;
    double **rho_save;

    rho_eff =
        dmatrix(
            0, nr + 1,
            0, nz + 1);

    ikloop
    {
        rho_eff[i][k] =
            rho_v[i][k] +
            alpha *
                phi3[i][k] /
                kappa;
    }

    augc(
        rho_eff,
        nr,
        nz);

    source_stage(
        tu,
        tw,
        workp,
        nr,
        nz,
        alpha);

    rho_save = rho_v;

    rho_v = rho_eff;

    solve_Poisson_relaxation(
        p,
        workp);

    rho_v = rho_save;

    grad_p(
        p,
        worku,
        workw,
        nr,
        nz);

    i0kloop
    {
        rho_face =
            0.5 *
            (rho_v[i + 1][k] +
             rho_v[i][k]);

        phi_face =
            0.5 *
            (phi3[i + 1][k] +
             phi3[i][k]);

        rho_eff_face =
            rho_face +
            alpha *
                phi_face /
                kappa;

        tu[i][k] =
            tu[i][k] -
            alpha *
                worku[i][k] /
                rho_eff_face;
    }

    ik0loop
    {
        rho_face =
            0.5 *
            (rho_v[i][k + 1] +
             rho_v[i][k]);

        phi_face =
            0.5 *
            (phi3[i][k + 1] +
             phi3[i][k]);

        rho_eff_face =
            rho_face +
            alpha *
                phi_face /
                kappa;

        tw[i][k] =
            tw[i][k] -
            alpha *
                workw[i][k] /
                rho_eff_face;
    }

    free_dmatrix(
        rho_eff,
        0, nr + 1,
        0, nz + 1);
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

void solve_Poisson_relaxation(double **u, double **f)
{
    extern int nr, nz;
    extern double **rho_v;

    int i, k, max_it = 1000, it_relax = 1;
    double tol = 1.0e-8, resid = 1.0, **work_t, **sor;

    work_t = dmatrix(1, nr, 1, nz);
    sor = dmatrix(1, nr, 1, nz);
    mat_copy(work_t, u, 1, nr, 1, nz);

    while (it_relax <= max_it && resid >= tol)
    {

        relax(u, f, rho_v, nr, nz);

        pressure_update(u);

        ikloop
            sor[i][k] = work_t[i][k] - u[i][k];

        resid = mat_max(sor, 1, nr, 1, nz);

        mat_copy(work_t, u, 1, nr, 1, nz);
        it_relax++;
    }

    printf("Pressure relaxation iteration = %d  residual = %16.14f \n", it_relax - 1, resid);
    free_dmatrix(work_t, 1, nr, 1, nz);
    free_dmatrix(sor, 1, nr, 1, nz);

    return;
}

void relax(double **p, double **f, double **w, int nrt, int nzt)
{
    extern int p_relax;
    extern double rright;

    int i, k, iter;
    double r_c, r_p, r_m, drt, drt2, a[4], coef, src;

    drt = rright / (double)nrt;
    drt2 = pow(drt, 2);

    for (iter = 1; iter <= p_relax; iter++)
    {

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

            if (i == 1)
            {
                src -= a[0] * r_p * p[i + 1][k] / (r_c * drt2);
                coef -= a[0] * (r_p) / (r_c * drt2);
            }

            else if (i == nrt)
            {
                src -= a[1] * r_m * p[i - 1][k] / (r_c * drt2);
                coef -= a[1] * r_m / (r_c * drt2);
            }

            else
            {
                src -= (a[0] * r_p * p[i + 1][k] + a[1] * r_m * p[i - 1][k]) / (r_c * drt2);
                coef -= (a[0] * r_p + a[1] * r_m) / (r_c * drt2);
            }

            if (k == 1)
            {
                src -= (a[2] * p[i][k + 1]) / drt2;
                coef -= a[2] / drt2;
            }

            else if (k == nzt)
            {
                src -= (a[3] * p[i][k - 1]) / drt2;
                coef -= a[3] / drt2;
            }

            else
            {
                src -= (a[2] * p[i][k + 1] + a[3] * p[i][k - 1]) / drt2;
                coef -= (a[2] + a[3]) / drt2;
            }

            p[i][k] = src / coef;
        }
    }
}

void grad_p(double **p, double **dpdr, double **dpdz, int nrt, int nzt)
{
    extern double rright;

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

void surface_tension(double **phi2, double **fr, double **fz)
{
    extern int nr, nz;
    extern double dr, gam, sig;

    int i, k;

    double rmain, rplus, rminus;
    double curvature;
    double adr, adz, grad_center;
    double fac;

    double **phi;
    double **dpdr;
    double **dpdz;
    double **grad_norm;
    double **cell_fr;
    double **cell_fz;

    const double tol = 1.0e-12;

    fac = 6.0 * sqrt(2.0) * gam * sig;

    phi = dmatrix(0, nr + 1, 0, nz + 1);

    dpdr = dmatrix(0, nr, 0, nz);
    dpdz = dmatrix(0, nr, 0, nz);
    grad_norm = dmatrix(0, nr, 0, nz);

    cell_fr = dmatrix(1, nr, 1, nz);
    cell_fz = dmatrix(1, nr, 1, nz);

    mat_copy(phi, phi2, 1, nr, 1, nz);
    augc(phi, nr, nz);

    for (i = 0; i <= nr; i++)
    {
        for (k = 0; k <= nz; k++)
        {
            dpdr[i][k] =
                (phi[i + 1][k + 1] + phi[i + 1][k] - phi[i][k + 1] - phi[i][k]) / (2.0 * dr);

            dpdz[i][k] =
                (phi[i + 1][k + 1] - phi[i + 1][k] + phi[i][k + 1] - phi[i][k]) / (2.0 * dr);

            grad_norm[i][k] =
                sqrt(dpdr[i][k] * dpdr[i][k] + dpdz[i][k] * dpdz[i][k]);
        }
    }

    ikloop
    {
        cell_fr[i][k] = 0.0;
        cell_fz[i][k] = 0.0;

        rmain =
            ((double)i - 0.5) * dr;

        rplus =
            rmain + 0.5 * dr;

        rminus =
            rmain - 0.5 * dr;

        if (grad_norm[i][k] < tol ||
            grad_norm[i][k - 1] < tol ||
            grad_norm[i - 1][k] < tol ||
            grad_norm[i - 1][k - 1] < tol)
        {
            curvature = 0.0;
        }
        else
        {

            curvature =
                (rplus *
                     (dpdr[i][k] / grad_norm[i][k]

                      +

                      dpdr[i][k - 1] / grad_norm[i][k - 1])

                 -

                 rminus *
                     (dpdr[i - 1][k] / grad_norm[i - 1][k]

                      +

                      dpdr[i - 1][k - 1] / grad_norm[i - 1][k - 1])) /
                    (2.0 * rmain * dr)

                +

                (dpdz[i][k] / grad_norm[i][k]

                 -

                 dpdz[i][k - 1] / grad_norm[i][k - 1]

                 +

                 dpdz[i - 1][k] / grad_norm[i - 1][k]

                 -

                 dpdz[i - 1][k - 1] / grad_norm[i - 1][k - 1]) /
                    (2.0 * dr);
        }

        adr =
            0.25 *
            (dpdr[i][k] + dpdr[i - 1][k] + dpdr[i][k - 1] + dpdr[i - 1][k - 1]);

        adz =
            0.25 *
            (dpdz[i][k] + dpdz[i - 1][k] + dpdz[i][k - 1] + dpdz[i - 1][k - 1]);

        grad_center =
            sqrt(adr * adr + adz * adz);

        cell_fr[i][k] =
            fac * curvature * grad_center * adr;

        cell_fz[i][k] =
            fac * curvature * grad_center * adz;
    }

    for (k = 1; k <= nz; k++)
    {
        fr[0][k] = 0.0;
        fr[nr][k] = 0.0;

        for (i = 1; i < nr; i++)
        {
            fr[i][k] =
                0.5 *
                (cell_fr[i][k] + cell_fr[i + 1][k]);
        }
    }

    for (i = 1; i <= nr; i++)
    {
        fz[i][0] = 0.0;
        fz[i][nz] = 0.0;

        for (k = 1; k < nz; k++)
        {
            fz[i][k] =
                0.5 *
                (cell_fz[i][k] + cell_fz[i][k + 1]);
        }
    }

    free_dmatrix(phi,
                 0, nr + 1,
                 0, nz + 1);

    free_dmatrix(dpdr,
                 0, nr,
                 0, nz);

    free_dmatrix(dpdz,
                 0, nr,
                 0, nz);

    free_dmatrix(grad_norm,
                 0, nr,
                 0, nz);

    free_dmatrix(cell_fr,
                 1, nr,
                 1, nz);

    free_dmatrix(cell_fz,
                 1, nr,
                 1, nz);
}

void cahn(double **c_old, double **cc_old, double **c_new, double **mu, double theta, double **adv_phi)
{
    extern int nr, nz;
    extern double **Mo1;

    int it_relax = 1, max_it_CH = 500;
    double resid = 1.0, tol = 1.0e-6, **sc, **smu, **ct;
    ct = dmatrix(1, nr, 1, nz);

    sc = dmatrix(1, nr, 1, nz);
    smu = dmatrix(1, nr, 1, nz);

    mat_copy(ct, c_old, 1, nr, 1, nz);

    source_ch(sc, smu, c_old, cc_old, theta, adv_phi);

    while (it_relax <= max_it_CH && resid > tol)
    {

        relax_ch(c_new, Mo1, mu, sc, smu, nr, nz);

        resid = error(ct, c_new, nr, nz);
        mat_copy(ct, c_new, 1, nr, 1, nz);

        it_relax++;
    }
    printf("cahn %16.14f   %d\n", resid, it_relax - 1);
    free_dmatrix(ct, 1, nr, 1, nz);
    free_dmatrix(sc, 1, nr, 1, nz);
    free_dmatrix(smu, 1, nr, 1, nz);
}

void source_ch(double **sc, double **smu, double **c_old, double **cc_old, double theta, double **adv_c)
{
    extern int nr, nz, it;
    extern double dt, gam, SS, **normphi3, **beta;
    int i, k;
    double **intc;
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

void relax_ch(double **c_new, double **Mo1, double **mu_new, double **sc,
              double **smu, int nrt, int nzt)
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

double error(double **c_old, double **c_new, int nrt, int nzt)
{
    double **uc, value;

    uc = dmatrix(1, nrt, 1, nzt);
    mat_sub(uc, c_new, c_old, 1, nrt, 1, nzt);

    value = mat_max(uc, 1, nrt, 1, nzt);
    free_dmatrix(uc, 1, nrt, 1, nzt);

    return value;
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

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free(m[nrl] + ncl - NR_END);
    free(m + nrl - NR_END);

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

void mat_sub(double **a, double **b, double **phi1,
             int nrl, int nrh, int ncl, int nch)
{
    int i, k;

    for (i = nrl; i <= nrh; i++)
        for (k = ncl; k <= nch; k++)
            a[i][k] = b[i][k] - phi1[i][k];

    return;
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

    sprintf(bufferphi1, OUTPUT_DIR "phi1.m");
    sprintf(bufferphi2, OUTPUT_DIR "phi2.m");
    sprintf(bufferphi3, OUTPUT_DIR "phi3.m");
    sprintf(bufferu, OUTPUT_DIR "u.m");
    sprintf(bufferv, OUTPUT_DIR "w.m");
    sprintf(bufferp, OUTPUT_DIR "p.m");

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
    extern double h, Re;
    extern double **rho_v, **vis, **fr;
    extern double **adv_u;

    int i, k;
    double r_c, r_p, r_m;
    double rho_face;
    double visc_term;

    auguw(u, w, nr, nz);

    i0kloop
    {
        r_c = ((double)i) * h;
        r_p = r_c + 0.5 * h;
        r_m = r_c - 0.5 * h;

        rho_face =
            0.5 * (rho_v[i + 1][k] + rho_v[i][k]);

        visc_term =
            ((2.0 * r_p * vis[i + 1][k] * (u[i + 1][k] - u[i][k]) - 2.0 * r_m * vis[i][k] * (u[i][k] - u[i - 1][k])) / r_c

             - (vis[i + 1][k] + vis[i][k]) * u[i][k] * h * h / (r_c * r_c)

             + 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (u[i][k + 1] - u[i][k])

             - 0.25 * (vis[i][k - 1] + vis[i + 1][k - 1] + vis[i][k] + vis[i + 1][k]) * (u[i][k] - u[i][k - 1])

             + 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (w[i + 1][k] - w[i][k])

             - 0.25 * (vis[i][k - 1] + vis[i + 1][k - 1] + vis[i][k] + vis[i + 1][k]) * (w[i + 1][k - 1] - w[i][k - 1])) /
            (h * h * Re * rho_face);

        result[i][k] =
            -adv_u[i][k] + fr[i][k] / rho_face + visc_term;
    }
}

void functionq(double **u, double **w, double **result)
{
    extern int nr, nz;
    extern double h, Re;
    extern double **rho_v, **vis, **fz;
    extern double Fr, gravity;
    extern double **adv_w;

    int i, k;
    double r_c, r_p, r_m;
    double rho_face;
    double visc_term;

    auguw(u, w, nr, nz);

    ik0loop
    {
        r_c = ((double)i - 0.5) * h;
        r_p = r_c + 0.5 * h;
        r_m = r_c - 0.5 * h;

        rho_face =
            0.5 * (rho_v[i][k + 1] + rho_v[i][k]);

        visc_term =
            ((
                 r_p * 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (w[i + 1][k] - w[i][k])

                 - r_m * 0.25 * (vis[i - 1][k] + vis[i][k] + vis[i - 1][k + 1] + vis[i][k + 1]) * (w[i][k] - w[i - 1][k])) /
                 r_c

             + 2.0 * vis[i][k + 1] * (w[i][k + 1] - w[i][k])

             - 2.0 * vis[i][k] * (w[i][k] - w[i][k - 1])

             + (r_p * 0.25 * (vis[i][k] + vis[i + 1][k] + vis[i][k + 1] + vis[i + 1][k + 1]) * (u[i][k + 1] - u[i][k])

                - r_m * 0.25 * (vis[i - 1][k] + vis[i][k] + vis[i - 1][k + 1] + vis[i][k + 1]) * (u[i - 1][k + 1] - u[i - 1][k])) /
                   r_c) /
            (h * h * Re * rho_face);

        result[i][k] =
            -adv_w[i][k] + fz[i][k] / rho_face - gravity / Fr + visc_term;
    }
}
