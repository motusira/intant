#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-12

static const double MID1[3][2] = {
    {0.95, 1.00},
    {1.05, 1.00},
    {1.10, 1.00},
};

static const double R1_TOMO[3][2] = {{1, 1}, {1, 1}, {1, 1}};
static const double R1_REGR[3][2] = {{1, 0}, {1, 0}, {1, 0}};

static const double MID2[3][3] = {
    {1.10, 0.90, 1.10},
    {1.40, 1.00, 0.80},
    {0.80, 1.40, 1.20},
};

static const char *outdir = ".";

static FILE *open_out(const char *name)
{
    char path[1024];
    snprintf(path, sizeof path, "%s/%s", outdir, name);
    FILE *f = fopen(path, "w");
    if (!f) {
        perror(path);
        exit(1);
    }
    return f;
}

static double det3(const double a[3][3])
{
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
         - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

static void inv3(const double a[3][3], double r[3][3])
{
    double d = det3(a);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            int r0 = (j + 1) % 3, r1 = (j + 2) % 3;
            int c0 = (i + 1) % 3, c1 = (i + 2) % 3;
            r[i][j] = (a[r0][c0] * a[r1][c1] - a[r0][c1] * a[r1][c0]) / d;
        }
}

static void svd32(const double a[3][2], double *smax, double *smin, double *gram_det)
{
    double g11 = 0, g12 = 0, g22 = 0;
    for (int i = 0; i < 3; i++) {
        g11 += a[i][0] * a[i][0];
        g12 += a[i][0] * a[i][1];
        g22 += a[i][1] * a[i][1];
    }
    double tr = g11 + g22, dt = g11 * g22 - g12 * g12;
    double disc = sqrt(fmax(tr * tr / 4 - dt, 0));
    *smax = sqrt(fmax(tr / 2 + disc, 0));
    *smin = sqrt(fmax(tr / 2 - disc, 0));
    *gram_det = dt;
}

static int row_feasible(const double mid[2], const double rad[2], double lambda)
{
    double lo2 = mid[1] - rad[1], hi2 = mid[1] + rad[1];
    double plo = fmin(lambda * lo2, lambda * hi2);
    double phi = fmax(lambda * lo2, lambda * hi2);
    return plo <= mid[0] + rad[0] + EPS && mid[0] - rad[0] <= phi + EPS;
}

static int a1_singular_for_lambda(const double R[3][2], double delta, double lambda)
{
    for (int i = 0; i < 3; i++) {
        double rad[2] = {delta * R[i][0], delta * R[i][1]};
        if (!row_feasible(MID1[i], rad, lambda))
            return 0;
    }
    return 1;
}

static double a1_delta_of_lambda(const double R[3][2], double lambda)
{
    double d = 0;
    for (int i = 0; i < 3; i++) {
        double den = R[i][0] + fabs(lambda) * R[i][1];
        double num = fabs(lambda * MID1[i][1] - MID1[i][0]);
        d = fmax(d, den > 0 ? num / den : (num > 0 ? INFINITY : 0));
    }
    return d;
}

static double a1_min_lambda(const double R[3][2], double lo, double hi)
{
    const double phi = (sqrt(5) - 1) / 2;
    double x1 = hi - phi * (hi - lo), x2 = lo + phi * (hi - lo);
    double f1 = a1_delta_of_lambda(R, x1), f2 = a1_delta_of_lambda(R, x2);
    for (int it = 0; it < 200; it++) {
        if (f1 < f2) {
            hi = x2; x2 = x1; f2 = f1;
            x1 = hi - phi * (hi - lo); f1 = a1_delta_of_lambda(R, x1);
        } else {
            lo = x1; x1 = x2; f1 = f2;
            x2 = lo + phi * (hi - lo); f2 = a1_delta_of_lambda(R, x2);
        }
    }
    return (lo + hi) / 2;
}

static int a1_singular_grid(const double R[3][2], double delta)
{
    for (int k = 0; k <= 200000; k++) {
        double lambda = 0.5 + 1.0 * k / 200000;
        if (a1_singular_for_lambda(R, delta, lambda))
            return 1;
    }
    return 0;
}

static double a1_bisect(const double R[3][2])
{
    double lo = 0, hi = 0.5;
    for (int it = 0; it < 40; it++) {
        double m = (lo + hi) / 2;
        if (a1_singular_grid(R, m)) hi = m; else lo = m;
    }
    return hi;
}

static void a1_build(const double R[3][2], double delta, double lambda, double A[3][2])
{
    for (int i = 0; i < 3; i++) {
        double m1 = MID1[i][0], m2 = MID1[i][1];
        double r1 = delta * R[i][0], r2 = delta * R[i][1];
        double res = lambda * m2 - m1;
        double den = r1 + fabs(lambda) * r2;
        double t = den > 0 ? res / den : 0;
        if (t > 1) t = 1;
        if (t < -1) t = -1;
        A[i][0] = m1 + t * r1;
        A[i][1] = m2 - t * r2 * (lambda >= 0 ? 1 : -1);
    }
}

static int a1_contains(const double R[3][2], double delta, const double A[3][2])
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 2; j++)
            if (fabs(A[i][j] - MID1[i][j]) > delta * R[i][j] + 1e-12)
                return 0;
    return 1;
}

static void print_m32(FILE *f, const char *name, const double A[3][2])
{
    fprintf(f, "%s =\n", name);
    for (int i = 0; i < 3; i++)
        fprintf(f, "    [ %12.9f  %12.9f ]\n", A[i][0], A[i][1]);
}

static void solve_a1(FILE *log, const char *tag, const char *title, const double R[3][2])
{
    fprintf(log, "==== A1 (%s) ====\n", title);
    {
        double smax, smin, gd;
        svd32(MID1, &smax, &smin, &gd);
        double rn = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 2; j++)
                rn += R[i][j] * R[i][j];
        fprintf(log, "mid A1: sigma_max = %.9f, sigma_min = %.9f, cond = %.3f, sigma_min/||R||_2 = %.9f\n",
                smax, smin, smax / smin, smin / sqrt(rn));
    }

    double lambda = a1_min_lambda(R, 0.5, 1.5);
    double delta = a1_delta_of_lambda(R, lambda);
    double delta_bis = a1_bisect(R);
    fprintf(log, "lambda*            = %.12f\n", lambda);
    fprintf(log, "delta* (minimax)   = %.12f\n", delta);
    fprintf(log, "delta* (bisection) = %.12f\n", delta_bis);
    fprintf(log, "singular(delta*-1e-6) = %d, singular(delta*+1e-9) = %d\n",
            a1_singular_grid(R, delta - 1e-6), a1_singular_for_lambda(R, delta + 1e-9, lambda));

    double A[3][2];
    a1_build(R, delta, lambda, A);
    print_m32(log, "A1'", A);
    double smax, smin, gd;
    svd32(A, &smax, &smin, &gd);
    fprintf(log, "A1' in A1: %s\n", a1_contains(R, delta, A) ? "yes" : "no");
    fprintf(log, "ratios a_i1/a_i2: %.12f %.12f %.12f\n",
            A[0][0] / A[0][1], A[1][0] / A[1][1], A[2][0] / A[2][1]);
    fprintf(log, "minors 2x2: %.3e %.3e %.3e\n",
            A[0][0] * A[1][1] - A[0][1] * A[1][0],
            A[0][0] * A[2][1] - A[0][1] * A[2][0],
            A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    fprintf(log, "det(A'^T A') = %.3e, sigma_max = %.9f, sigma_min = %.3e\n\n", gd, smax, smin);

    char name[256];
    snprintf(name, sizeof name, "a1_%s_delta_lambda.dat", tag);
    FILE *f = open_out(name);
    fprintf(f, "# lambda delta_min(lambda)\n");
    for (int k = 0; k <= 1000; k++) {
        double l = 0.85 + 0.35 * k / 1000;
        fprintf(f, "%.6f %.9f\n", l, a1_delta_of_lambda(R, l));
    }
    fclose(f);

    snprintf(name, sizeof name, "a1_%s_rows.dat", tag);
    f = open_out(name);
    fprintf(f, "# x=a_i2 y=a_i1 dx dy  x'=a'_i2 y'=a'_i1\n");
    for (int i = 0; i < 3; i++)
        fprintf(f, "%.9f %.9f %.9f %.9f %.9f %.9f\n", MID1[i][1], MID1[i][0],
                delta * R[i][1], delta * R[i][0], A[i][1], A[i][0]);
    fclose(f);

    snprintf(name, sizeof name, "a1_%s_params.gp", tag);
    f = open_out(name);
    fprintf(f, "delta_star = %.12f\nlambda_star = %.12f\n", delta, lambda);
    fclose(f);
}

static void a2_det_range(double delta, double *dmin, double *dmax)
{
    *dmin = INFINITY;
    *dmax = -INFINITY;
    for (int mask = 0; mask < 512; mask++) {
        double a[3][3];
        for (int k = 0; k < 9; k++)
            a[k / 3][k % 3] = MID2[k / 3][k % 3] + ((mask >> k) & 1 ? delta : -delta);
        double d = det3(a);
        if (d < *dmin) *dmin = d;
        if (d > *dmax) *dmax = d;
    }
}

static void solve_a2(FILE *log)
{
    fprintf(log, "==== A2 ====\n");
    double C[3][3];
    double d0 = det3(MID2);
    inv3(MID2, C);
    fprintf(log, "det(mid A2) = %.12f\n", d0);
    {
        double v[3] = {1, 1, 1}, nrm = 0;
        for (int it = 0; it < 500; it++) {
            double w[3] = {0, 0, 0}, u[3] = {0, 0, 0};
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    w[i] += C[i][j] * v[j];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    u[i] += C[j][i] * w[j];
            nrm = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
            for (int i = 0; i < 3; i++)
                v[i] = u[i] / nrm;
        }
        double smin = 1 / sqrt(nrm);
        fprintf(log, "mid A2: sigma_min = %.9f, sigma_min/||ee^T||_2 = %.9f\n", smin, smin / 3);
    }
    fprintf(log, "inv(mid A2) =\n");
    for (int i = 0; i < 3; i++)
        fprintf(log, "    [ %10.6f %10.6f %10.6f ]\n", C[i][0], C[i][1], C[i][2]);

    double best = -INFINITY;
    int by = 0, bz = 0;
    FILE *fv = open_out("a2_signs.dat");
    fprintf(fv, "# y1 y2 y3 z1 z2 z3 z^T C y\n");
    for (int ym = 0; ym < 8; ym++)
        for (int zm = 0; zm < 8; zm++) {
            double y[3], z[3], s = 0;
            for (int k = 0; k < 3; k++) {
                y[k] = (ym >> k) & 1 ? -1 : 1;
                z[k] = (zm >> k) & 1 ? -1 : 1;
            }
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    s += z[i] * C[i][j] * y[j];
            fprintf(fv, "%+.0f %+.0f %+.0f %+.0f %+.0f %+.0f %.9f\n",
                    y[0], y[1], y[2], z[0], z[1], z[2], s);
            if (s > best) {
                best = s;
                by = ym;
                bz = zm;
            }
        }
    fclose(fv);

    double y[3], z[3];
    for (int k = 0; k < 3; k++) {
        y[k] = (by >> k) & 1 ? -1 : 1;
        z[k] = (bz >> k) & 1 ? -1 : 1;
    }
    double delta = 1.0 / best;
    fprintf(log, "max z^T C y = %.12f, y = (%+.0f %+.0f %+.0f), z = (%+.0f %+.0f %+.0f)\n",
            best, y[0], y[1], y[2], z[0], z[1], z[2]);
    fprintf(log, "delta* (Rohn)      = %.12f\n", delta);

    double lo = 0, hi = 1;
    for (int it = 0; it < 60; it++) {
        double m = (lo + hi) / 2, dmin, dmax;
        a2_det_range(m, &dmin, &dmax);
        if (dmin <= 0 && dmax >= 0) hi = m; else lo = m;
    }
    fprintf(log, "delta* (bisection) = %.12f\n", hi);

    double A[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i][j] = MID2[i][j] - delta * y[i] * z[j];
    fprintf(log, "A2' =\n");
    int inside = 1;
    for (int i = 0; i < 3; i++) {
        fprintf(log, "    [ %12.9f %12.9f %12.9f ]\n", A[i][0], A[i][1], A[i][2]);
        for (int j = 0; j < 3; j++)
            if (fabs(A[i][j] - MID2[i][j]) > delta + 1e-12)
                inside = 0;
    }
    fprintf(log, "A2' in A2: %s\n", inside ? "yes" : "no");
    fprintf(log, "det A2' = %.3e\n", det3(A));

    double x[3] = {0, 0, 0}, r[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            x[i] += C[i][j] * y[j];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            r[i] += A[i][j] * x[j];
    fprintf(log, "kernel x = (%.6f %.6f %.6f), |A2' x| = %.3e\n\n",
            x[0], x[1], x[2], sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));

    FILE *f = open_out("a2_det_range.dat");
    fprintf(f, "# delta det_min det_max\n");
    for (int k = 0; k <= 400; k++) {
        double d = 2.5 * delta * k / 400, dmin, dmax;
        a2_det_range(d, &dmin, &dmax);
        fprintf(f, "%.9f %.9f %.9f\n", d, dmin, dmax);
    }
    fclose(f);

    f = open_out("a2_det_path.dat");
    fprintf(f, "# t det(mid - t y z^T)\n");
    for (int k = 0; k <= 400; k++) {
        double t = 2.5 * delta * k / 400, a[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                a[i][j] = MID2[i][j] - t * y[i] * z[j];
        fprintf(f, "%.9f %.9f\n", t, det3(a));
    }
    fclose(f);

    f = open_out("a2_params.gp");
    fprintf(f, "delta_star = %.12f\ndet_mid = %.12f\n", delta, d0);
    fclose(f);
}

int main(int argc, char **argv)
{
    if (argc > 1)
        outdir = argv[1];
    FILE *log = open_out("results.txt");
    solve_a1(log, "tomo", "томография, rad = delta * ones(3,2)", R1_TOMO);
    solve_a1(log, "regr", "регрессия, rad = delta * [1 0; 1 0; 1 0]", R1_REGR);
    solve_a2(log);
    fclose(log);
    return 0;
}
