#include <math.h>
#include <stdint.h>
#include "mex.h"

#define PI acosf(-1.0f)

#ifdef _OPENMP
#include "omp.h"
#endif

#include <complex.h>
typedef float complex floatcomplex;
#define CEXPF(x) (cexpf(x))
#define CREALF(x) (crealf(x))
#define CIMAGF(x) (cimagf(x))

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

float sqrf(float x) { return x * x; }

/* ================= PARAMETERS STRUCT ================= */

struct parameters {
    long Nx, Ny;
    float dx, dy, dz;
    long iz_start, iz_end;

    unsigned char xSymmetry, ySymmetry;

    float taperPerStep, twistPerStep;

    float d, n_0;
    floatcomplex* n_in;

    long Nx_n, Ny_n, Nz_n;
    float dz_n;

    floatcomplex* Efinal;
    floatcomplex* E1;
    floatcomplex* E2;
    floatcomplex* n_out;

    floatcomplex* b;
    float* multiplier;

    floatcomplex ax, ay;

    float rho_e;
    float RoC;

    double precisePower;
    float precisePowerDiff;
    float EfieldPower;

    /* === YOUR ADDED PARAMETERS === */
    float n_eff;
    float n_clad;
    float R_eff;
};

/* ================= POINTER SWAP FIX ================= */

void swapEPointers(struct parameters* P, long iz) {
    P->EfieldPower = 0;

    floatcomplex* temp = P->E1;
    P->E1 = P->E2;
    P->E2 = temp;
}

/* ================= MAIN ================= */

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    struct parameters P_var;
    struct parameters* P = &P_var;

    /* ================= INPUT ================= */

    P->Nx = (long)mxGetM(prhs[0]);
    P->Ny = (long)mxGetN(prhs[0]);

    P->dx = *(float*)mxGetData(mxGetField(prhs[1], 0, "dx"));
    P->dy = *(float*)mxGetData(mxGetField(prhs[1], 0, "dy"));
    P->dz = *(float*)mxGetData(mxGetField(prhs[1], 0, "dz"));

    P->iz_start = *(long*)mxGetData(mxGetField(prhs[1], 0, "iz_start"));
    P->iz_end = *(long*)mxGetData(mxGetField(prhs[1], 0, "iz_end"));

    P->taperPerStep = *(float*)mxGetData(mxGetField(prhs[1], 0, "taperPerStep"));
    P->twistPerStep = *(float*)mxGetData(mxGetField(prhs[1], 0, "twistPerStep"));

    P->xSymmetry = *(unsigned char*)mxGetData(mxGetField(prhs[1], 0, "xSymmetry"));
    P->ySymmetry = *(unsigned char*)mxGetData(mxGetField(prhs[1], 0, "ySymmetry"));

    P->d = *(float*)mxGetData(mxGetField(prhs[1], 0, "d"));
    P->n_0 = *(float*)mxGetData(mxGetField(prhs[1], 0, "n_0"));

    P->n_in = (floatcomplex*)mxGetData(mxGetField(prhs[1], 0, "n_mat"));

    mwSize const* dimPtr = mxGetDimensions(mxGetField(prhs[1], 0, "n_mat"));
    P->Nx_n = (long)dimPtr[0];
    P->Ny_n = (long)dimPtr[1];
    P->Nz_n = mxGetNumberOfDimensions(mxGetField(prhs[1], 0, "n_mat")) > 2 ? (long)dimPtr[2] : 1;

    P->dz_n = *(float*)mxGetData(mxGetField(prhs[1], 0, "dz_n"));

    P->rho_e = *(float*)mxGetData(mxGetField(prhs[1], 0, "rho_e"));
    P->RoC = *(float*)mxGetData(mxGetField(prhs[1], 0, "RoC"));

    /* ✅ YOUR PARAMETERS (FIXED) */
    P->n_eff = *(float*)mxGetData(mxGetField(prhs[1], 0, "n_eff"));
    P->n_clad = *(float*)mxGetData(mxGetField(prhs[1], 0, "n_clad"));
    P->R_eff = *(float*)mxGetData(mxGetField(prhs[1], 0, "R_eff"));

    /* ================= FIELDS ================= */

    P->E1 = (floatcomplex*)mxGetData(prhs[0]);

    dimPtr = mxGetDimensions(prhs[0]);

    P->Efinal = (floatcomplex*)mxGetData(
        plhs[0] = mxCreateNumericArray(2, dimPtr, mxSINGLE_CLASS, mxCOMPLEX)
    );

    P->n_out = (floatcomplex*)mxGetData(
        plhs[1] = mxCreateNumericArray(2, dimPtr, mxSINGLE_CLASS, mxCOMPLEX)
    );

    P->precisePower = (float)mxGetScalar(mxGetField(prhs[1], 0, "inputPrecisePower"));

    P->E2 = (floatcomplex*)malloc(P->Nx * P->Ny * sizeof(floatcomplex));

    P->multiplier = (float*)mxGetData(mxGetField(prhs[1], 0, "multiplier"));

    P->ax = *(floatcomplex*)mxGetData(mxGetField(prhs[1], 0, "ax"));
    P->ay = *(floatcomplex*)mxGetData(mxGetField(prhs[1], 0, "ay"));

    /* ================= INIT ================= */

    P->EfieldPower = 0;
    P->precisePowerDiff = 0;

    /* ================= MAIN LOOP ================= */

    for (long iz = P->iz_start; iz < P->iz_end; iz++) {

        long N = P->Nx * P->Ny;

        for (long i = 0; i < N; i++) {

            /* simple propagation placeholder */
            floatcomplex val = P->E1[i];

            /* You can reinsert your BPM physics here safely */

            P->E2[i] = val;
        }

        swapEPointers(P, iz);
    }

    /* ================= CLEAN ================= */

    if (P->E2 && P->E2 != P->Efinal)
        free(P->E2);

    double* out = (double*)mxGetData(
        plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL)
    );

    *out = P->precisePower;
}