#include "mechanism.H"
const int rmap[272] = {1,13,14,21,22,236,243,28,150,197,269,9,10,11,12,33,43,44,45,76,98,101,132,165,166,176,188,206,214,219,222,233,234,240,257,258,0,2,3,4,5,6,7,8,15,16,17,18,19,20,23,24,25,26,27,29,30,31,32,34,35,36,37,38,39,40,41,42,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,99,100,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,151,152,153,154,155,156,157,158,159,160,161,162,163,164,167,168,169,170,171,172,173,174,175,177,178,179,180,181,182,183,184,185,186,187,189,190,191,192,193,194,195,196,198,199,200,201,202,203,204,205,207,208,209,210,211,212,213,215,216,217,218,220,221,223,224,225,226,227,228,229,230,231,232,235,237,238,239,241,242,244,245,246,247,248,249,250,251,252,253,254,255,256,259,260,261,262,263,264,265,266,267,268,270,271};

// Returns 0-based map of reaction order
void GET_RMAP(int * _rmap)
{
for (int j=0; j<272; ++j) {
_rmap[j] = rmap[j];
}
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
const int ns[272] =
     {4,3,4,4,4,4,4,4,3,2,2,3,3,3,3,3,4,4,4,3,3,2,2,4,4,4,4,4,3,4,4,4,4,3,4,4,4,4,4,5,3,4,3,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,5,3,4,4,3,4,4,4,5,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,3,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,3,4,4,4,3,2,4,4,4,4,4,4,3,4,4,3,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,3,3,4,5,5,4,4,4,5,3,3,4,4,4,4,4,4,3,4,4,4,4,2,4,4,3,4,4,4,3,4,4,4,4,4,4,2,3,4,3,4,4,4,3,3,4,3,4,4,4,4,4,4,4,3,4,4,4,4,4,3,3,4,4,4,2,4,4,4,3,3,3,2,4,3};
const int kiv[1360] =
     {11,4,45,13,0,11,3,21,0,0,25,3,45,12,0,11,21,9,3,0,11,4,28,3,0,11,8,12,13,0,9,12,11,13,0,9,13,11,15,0,15,12,13,0,0,9,11,0,0,0,12,8,0,0,0,11,12,13,0,0,11,13,15,0,0,11,8,14,0,0,11,8,14,0,0,11,14,13,0,0,9,8,11,14,0,14,12,8,13,0,14,13,15,8,0,14,16,8,0,0,14,16,8,0,0,16,13,0,0,0,16,13,0,0,0,11,16,15,13,0,11,16,9,14,0,16,12,14,13,0,16,13,15,14,0,16,13,15,14,0,5,12,6,0,0,5,8,6,12,0,5,13,6,11,0,5,13,6,11,0,5,14,6,13,0,17,5,11,0,0,17,8,5,14,0,11,17,5,9,0,17,12,5,13,0,17,12,6,11,0,17,13,5,15,0,17,14,6,11,13,17,5,9,0,0,18,8,17,14,0,17,18,5,0,0,1,29,11,0,0,1,29,11,0,0,1,7,0,0,0,29,9,11,1,0,1,12,11,30,0,1,12,29,13,0,1,12,5,28,0,1,13,29,15,0,1,13,11,32,0,1,13,11,2,0,1,13,5,27,0,1,8,29,14,0,29,1,11,31,0,11,7,11,1,0,7,12,5,28,0,7,13,11,2,0,29,7,11,31,0,29,12,5,25,0,29,13,11,30,0,29,8,30,12,0,29,8,30,12,0,29,8,5,3,0,29,6,5,30,0,18,29,1,17,0,29,3,25,30,0,29,19,30,3,0,29,19,5,4,0,29,19,6,45,0,29,21,1,3,0,29,23,1,19,0,29,4,34,3,0,29,2,1,30,0,29,30,5,34,0,2,5,28,0,0,11,2,5,27,0,11,2,9,30,0,2,12,30,13,0,2,12,6,28,0,2,12,5,21,0,2,13,15,30,0,2,14,16,30,0,2,8,6,21,0,2,28,30,27,0,11,32,11,2,0,11,32,5,27,0,11,32,9,30,0,32,12,30,13,0,32,13,15,30,0,32,27,30,20,0,33,1,12,0,0,11,33,1,13,0,33,12,17,3,0,33,13,18,3,0,33,12,30,13,0,33,13,5,9,3,33,11,30,0,0,33,13,15,30,0,33,13,17,21,0,30,5,25,0,0,11,30,5,28,0,30,12,5,3,0,30,13,5,22,0,30,13,5,11,3,14,30,2,8,0,30,8,6,3,0,17,30,5,2,0,18,30,17,2,0,30,3,5,4,0,30,3,6,45,0,30,19,5,3,0,30,19,6,4,0,21,30,2,3,0,23,30,2,19,0,25,30,5,45,0,30,20,2,27,0,30,5,45,0,0,40,11,1,0,0,11,40,9,1,0,40,12,1,13,0,40,13,15,1,0,40,13,15,1,0,40,8,18,3,0,40,27,1,20,0,41,11,1,0,0,11,41,11,40,0,11,41,9,1,0,41,12,11,2,0,41,12,1,13,0,41,13,15,1,0,20,11,27,0,0,11,20,9,27,0,20,12,27,13,0,20,13,15,27,0,14,20,16,27,0,11,27,9,28,0,27,12,11,21,0,27,12,28,13,0,27,13,15,28,0,14,27,24,13,0,14,27,20,8,0,27,8,24,12,0,27,8,21,13,0,27,28,20,0,0,28,27,25,20,0,25,27,11,45,0,27,3,15,45,0,28,12,11,3,0,4,45,12,0,0,4,13,14,45,0,27,3,26,13,0,25,13,11,3,0,27,19,15,4,0,27,19,24,3,0,11,28,9,25,0,28,13,11,21,0,28,13,15,25,0,28,8,21,12,0,28,8,3,13,0,28,11,45,0,0,25,28,11,45,0,28,3,45,13,0,28,19,4,13,0,24,11,21,0,0,24,42,0,0,0,11,24,9,21,0,11,24,27,13,0,24,12,21,13,0,24,13,15,21,0,24,14,16,21,0,24,8,21,14,0,24,3,21,0,0,24,27,21,20,0,24,19,21,23,0,42,11,21,0,0,11,42,27,13,0,11,42,9,21,0,42,12,21,13,0,42,12,21,13,0,42,13,15,21,0,42,14,16,21,0,42,8,21,14,0,42,27,21,20,0,42,19,21,23,0,21,8,14,3,0,17,21,18,3,0,22,11,3,0,0,11,22,11,21,0,11,22,28,13,0,22,12,3,13,0,22,13,11,23,0,22,8,19,13,0,17,3,5,21,0,11,23,21,13,0,11,23,15,3,0,19,12,43,0,0,19,3,43,0,0,5,19,6,3,0,17,19,5,3,13,17,19,6,11,3,11,43,19,13,0,43,12,19,8,0,43,13,14,19,0,14,43,19,8,13,43,3,8,0,0,26,11,45,0,0,11,26,9,45,0,26,12,11,4,0,26,12,45,13,0,26,12,28,3,0,26,13,15,45,0,26,8,14,45,0,26,11,45,0,0,28,26,45,27,0,27,26,45,20,0,26,3,21,45,0,5,4,6,45,0,31,29,0,0,0,31,12,29,30,0,31,13,29,32,0,11,12,0,0,0,15,0,15,13,0,9,0,9,13,0,45,0,45,13,0,13,0,13,0,0,11,0,11,13,0,10,0,10,13,0,8,0,8,13,0,6,0,6,13,0,5,0,5,13,0,11,4,45,0,0,45,25,0,0,0,25,12,3,0,0,14,3,19,13,0,3,13,23,0,0,9,19,11,23,0,11,19,3,13,0,19,12,3,8,0,19,3,12,0,0,19,3,8,0,0,14,19,23,8,0,19,13,44,0,0,3,19,4,8,0,23,12,19,13,0,23,13,15,19,0,23,15,3,19,0,21,12,3,13,0,21,13,15,3,0,21,19,23,3,0,21,15,4,0,0,21,3,4,13,0,25,19,4,12,0,25,8,3,12,0,25,4,45,3,0,4,3,45,19,0,24,9,3,0,0,14,3,44,0,0,11,36,9,19,0,36,12,19,13,0,36,13,15,19,0,36,23,0,0,0,23,27,20,19,0,21,25,28,3,0,4,12,45,8,0,4,12,3,0,0,43,3,8,0,0,43,19,8,0,0,37,19,0,0,0,37,12,38,8,0,38,12,19,0,0};
const int nuv[1360] =
     {-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,2,0,0,0,-2,1,0,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-2,1,1,0,0,-1,2,0,0,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-2,2,1,0,0,-1,-1,1,1,0,-2,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,2,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,2,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,2,0,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-2,2,1,0,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,1,1,0,0,-2,2,1,0,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,2,0,0};
if (*i < 1) {
// Return max num species per reaction
*nspec = 5;
} else {
if (*i > 272) {
*nspec = -1;
} else {
*nspec = ns[*i-1];
for (int j=0; j<*nspec; ++j) {
ki[j] = kiv[(*i-1)*5 + j] + 1;
nu[j] = nuv[(*i-1)*5 + j];
}
}
}
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  q_f, amrex::Real *  q_r)
{
int id; // loop counter
amrex::Real c[52]; // temporary storage
amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (id = 0; id < 52; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, *T);

// convert to chemkin units
for (id = 0; id < 272; ++id) {
q_f[id] *= 1.0e-6;
q_r[id] *= 1.0e-6;
}
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void progressRateFR(amrex::Real *  q_f, amrex::Real *  q_r, amrex::Real *  sc, amrex::Real T)
{
const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T };// temperature cache
amrex::Real invT = 1.0 / tc[1];
// compute the Gibbs free energy
amrex::Real g_RT[52];
gibbs(g_RT, tc);

amrex::Real sc_qss[1];
amrex::Real kf_qss[0], qf_qss[0], qr_qss[0];
comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

return;
}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 12.011000; // C
awt[1] = 1.008000; // H
awt[2] = 15.999000; // O
awt[3] = 14.007000; // N
awt[4] = 39.950000; // Ar
awt[5] = 4.002602; // He
}

// get atomic weight for all elements
void CKAWT( amrex::Real *  awt)
{
atomicWeight(awt);
}

// Returns the elemental composition 
// of the speciesi (mdim is num of elements)
void CKNCF(int * ncf)
{
int id; // loop counter
int kd = 6; 
// Zero ncf
for (id = 0; id < kd * 52; ++ id) {
 ncf[id] = 0; 
}

// OH*
ncf[ 0 * kd + 1 ] = 1; // H
ncf[ 0 * kd + 2 ] = 1; // O

// HCN
ncf[ 1 * kd + 0 ] = 1; // C
ncf[ 1 * kd + 1 ] = 1; // H
ncf[ 1 * kd + 3 ] = 1; // N

// HNCO
ncf[ 2 * kd + 0 ] = 1; // C
ncf[ 2 * kd + 1 ] = 1; // H
ncf[ 2 * kd + 3 ] = 1; // N
ncf[ 2 * kd + 2 ] = 1; // O

// NO
ncf[ 3 * kd + 3 ] = 1; // N
ncf[ 3 * kd + 2 ] = 1; // O

// N2O
ncf[ 4 * kd + 3 ] = 2; // N
ncf[ 4 * kd + 2 ] = 1; // O

// CO
ncf[ 5 * kd + 0 ] = 1; // C
ncf[ 5 * kd + 2 ] = 1; // O

// CO2
ncf[ 6 * kd + 0 ] = 1; // C
ncf[ 6 * kd + 2 ] = 2; // O

// HNC
ncf[ 7 * kd + 0 ] = 1; // C
ncf[ 7 * kd + 1 ] = 1; // H
ncf[ 7 * kd + 3 ] = 1; // N

// O2
ncf[ 8 * kd + 2 ] = 2; // O

// H2
ncf[ 9 * kd + 1 ] = 2; // H

// AR
ncf[ 10 * kd + 4 ] = 1; // Ar

// H
ncf[ 11 * kd + 1 ] = 1; // H

// O
ncf[ 12 * kd + 2 ] = 1; // O

// OH
ncf[ 13 * kd + 1 ] = 1; // H
ncf[ 13 * kd + 2 ] = 1; // O

// HO2
ncf[ 14 * kd + 1 ] = 1; // H
ncf[ 14 * kd + 2 ] = 2; // O

// H2O
ncf[ 15 * kd + 1 ] = 2; // H
ncf[ 15 * kd + 2 ] = 1; // O

// H2O2
ncf[ 16 * kd + 1 ] = 2; // H
ncf[ 16 * kd + 2 ] = 2; // O

// HCO
ncf[ 17 * kd + 0 ] = 1; // C
ncf[ 17 * kd + 1 ] = 1; // H
ncf[ 17 * kd + 2 ] = 1; // O

// CH2O
ncf[ 18 * kd + 0 ] = 1; // C
ncf[ 18 * kd + 1 ] = 2; // H
ncf[ 18 * kd + 2 ] = 1; // O

// NO2
ncf[ 19 * kd + 3 ] = 1; // N
ncf[ 19 * kd + 2 ] = 2; // O

// NH3
ncf[ 20 * kd + 1 ] = 3; // H
ncf[ 20 * kd + 3 ] = 1; // N

// HNO
ncf[ 21 * kd + 1 ] = 1; // H
ncf[ 21 * kd + 3 ] = 1; // N
ncf[ 21 * kd + 2 ] = 1; // O

// HON
ncf[ 22 * kd + 1 ] = 1; // H
ncf[ 22 * kd + 3 ] = 1; // N
ncf[ 22 * kd + 2 ] = 1; // O

// HONO
ncf[ 23 * kd + 1 ] = 1; // H
ncf[ 23 * kd + 3 ] = 1; // N
ncf[ 23 * kd + 2 ] = 2; // O

// H2NO
ncf[ 24 * kd + 1 ] = 2; // H
ncf[ 24 * kd + 3 ] = 1; // N
ncf[ 24 * kd + 2 ] = 1; // O

// N
ncf[ 25 * kd + 3 ] = 1; // N

// NNH
ncf[ 26 * kd + 1 ] = 1; // H
ncf[ 26 * kd + 3 ] = 2; // N

// NH2
ncf[ 27 * kd + 1 ] = 2; // H
ncf[ 27 * kd + 3 ] = 1; // N

// NH
ncf[ 28 * kd + 1 ] = 1; // H
ncf[ 28 * kd + 3 ] = 1; // N

// CN
ncf[ 29 * kd + 0 ] = 1; // C
ncf[ 29 * kd + 3 ] = 1; // N

// NCO
ncf[ 30 * kd + 0 ] = 1; // C
ncf[ 30 * kd + 3 ] = 1; // N
ncf[ 30 * kd + 2 ] = 1; // O

// NCCN
ncf[ 31 * kd + 0 ] = 2; // C
ncf[ 31 * kd + 3 ] = 2; // N

// HOCN
ncf[ 32 * kd + 0 ] = 1; // C
ncf[ 32 * kd + 1 ] = 1; // H
ncf[ 32 * kd + 3 ] = 1; // N
ncf[ 32 * kd + 2 ] = 1; // O

// HCNO
ncf[ 33 * kd + 0 ] = 1; // C
ncf[ 33 * kd + 1 ] = 1; // H
ncf[ 33 * kd + 3 ] = 1; // N
ncf[ 33 * kd + 2 ] = 1; // O

// NCN
ncf[ 34 * kd + 0 ] = 1; // C
ncf[ 34 * kd + 3 ] = 2; // N

// HE
ncf[ 35 * kd + 5 ] = 1; // He

// HNO2
ncf[ 36 * kd + 1 ] = 1; // H
ncf[ 36 * kd + 3 ] = 1; // N
ncf[ 36 * kd + 2 ] = 2; // O

// N2O4
ncf[ 37 * kd + 3 ] = 2; // N
ncf[ 37 * kd + 2 ] = 4; // O

// N2O3
ncf[ 38 * kd + 3 ] = 2; // N
ncf[ 38 * kd + 2 ] = 3; // O

// CH
ncf[ 39 * kd + 0 ] = 1; // C
ncf[ 39 * kd + 1 ] = 1; // H

// H2CN
ncf[ 40 * kd + 0 ] = 1; // C
ncf[ 40 * kd + 1 ] = 2; // H
ncf[ 40 * kd + 3 ] = 1; // N

// HCNH
ncf[ 41 * kd + 0 ] = 1; // C
ncf[ 41 * kd + 1 ] = 2; // H
ncf[ 41 * kd + 3 ] = 1; // N

// HNOH
ncf[ 42 * kd + 1 ] = 2; // H
ncf[ 42 * kd + 3 ] = 1; // N
ncf[ 42 * kd + 2 ] = 1; // O

// NO3
ncf[ 43 * kd + 3 ] = 1; // N
ncf[ 43 * kd + 2 ] = 3; // O

// HONO2
ncf[ 44 * kd + 1 ] = 1; // H
ncf[ 44 * kd + 3 ] = 1; // N
ncf[ 44 * kd + 2 ] = 3; // O

// N2
ncf[ 45 * kd + 3 ] = 2; // N

// CH4
ncf[ 46 * kd + 0 ] = 1; // C
ncf[ 46 * kd + 1 ] = 4; // H

// C2H6
ncf[ 47 * kd + 0 ] = 2; // C
ncf[ 47 * kd + 1 ] = 6; // H

// N2H4
ncf[ 48 * kd + 1 ] = 4; // H
ncf[ 48 * kd + 3 ] = 2; // N

// N2H3
ncf[ 49 * kd + 1 ] = 3; // H
ncf[ 49 * kd + 3 ] = 2; // N

// N2H2
ncf[ 50 * kd + 1 ] = 2; // H
ncf[ 50 * kd + 3 ] = 2; // N

// H2NN
ncf[ 51 * kd + 1 ] = 2; // H
ncf[ 51 * kd + 3 ] = 2; // N

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(6);
ename[0] = "C";
ename[1] = "H";
ename[2] = "O";
ename[3] = "N";
ename[4] = "Ar";
ename[5] = "He";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(52);
kname[0] = "OH*";
kname[1] = "HCN";
kname[2] = "HNCO";
kname[3] = "NO";
kname[4] = "N2O";
kname[5] = "CO";
kname[6] = "CO2";
kname[7] = "HNC";
kname[8] = "O2";
kname[9] = "H2";
kname[10] = "AR";
kname[11] = "H";
kname[12] = "O";
kname[13] = "OH";
kname[14] = "HO2";
kname[15] = "H2O";
kname[16] = "H2O2";
kname[17] = "HCO";
kname[18] = "CH2O";
kname[19] = "NO2";
kname[20] = "NH3";
kname[21] = "HNO";
kname[22] = "HON";
kname[23] = "HONO";
kname[24] = "H2NO";
kname[25] = "N";
kname[26] = "NNH";
kname[27] = "NH2";
kname[28] = "NH";
kname[29] = "CN";
kname[30] = "NCO";
kname[31] = "NCCN";
kname[32] = "HOCN";
kname[33] = "HCNO";
kname[34] = "NCN";
kname[35] = "HE";
kname[36] = "HNO2";
kname[37] = "N2O4";
kname[38] = "N2O3";
kname[39] = "CH";
kname[40] = "H2CN";
kname[41] = "HCNH";
kname[42] = "HNOH";
kname[43] = "NO3";
kname[44] = "HONO2";
kname[45] = "N2";
kname[46] = "CH4";
kname[47] = "C2H6";
kname[48] = "N2H4";
kname[49] = "N2H3";
kname[50] = "N2H2";
kname[51] = "H2NN";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<53; k++) {
for (int l=0; l<53; l++) {
if(Jac[ 53 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<53; k++) {
for (int l=0; l<53; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 53 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the simplified (for preconditioning) system Jacobian
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<53; k++) {
for (int l=0; l<53; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 53 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}
}

nJdata[0] = nJdata_tmp;
}


// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 53;
int offset_col = nc * 53;
for (int k=0; k<53; k++) {
for (int l=0; l<53; l++) {
if(Jac[53*k + l] != 0.0) {
rowVals[nJdata_tmp] = l + offset_row; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
colPtrs[offset_col + (k + 1)] = nJdata_tmp;
}
}
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 53;
for (int l=0; l<53; l++) {
for (int k=0; k<53; k++) {
if(Jac[53*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 53;
for (int l=0; l<53; l++) {
for (int k=0; k<53; k++) {
if(Jac[53*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
rowPtrs[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 53;
for (int l=0; l<53; l++) {
for (int k=0; k<53; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[53*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 53;
for (int l=0; l<53; l++) {
for (int k=0; k<53; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[53*k + l] != 0.0) {
colVals[nJdata_tmp] = k + offset; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[offset + (l + 1)] = nJdata_tmp;
}
}
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU
// BASE 0
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<53; k++) {
for (int l=0; l<53; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 53*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[53*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 53*k + l;
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
colPtrs[k+1] = nJdata_tmp;
}
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
amrex::GpuArray<amrex::Real,2809> Jac = {0.0};
amrex::GpuArray<amrex::Real,52> conc = {0.0};
for (int n=0; n<52; n++) {
    conc[n] = 1.0/ 52.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<53; l++) {
for (int k=0; k<53; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[53*k + l] != 0.0) {
colVals[nJdata_tmp-1] = k+1; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
} else {
rowPtr[0] = 0;
int nJdata_tmp = 0;
for (int l=0; l<53; l++) {
for (int k=0; k<53; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[53*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
