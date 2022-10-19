#include "mechanism.H"
const int rmap[277] = {8,9,16,17,234,242,23,194,266,274,4,5,6,7,28,42,43,44,75,97,100,131,162,163,173,185,203,211,216,219,230,231,235,239,258,259,0,1,2,3,10,11,12,13,14,15,18,19,20,21,22,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,98,99,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,164,165,166,167,168,169,170,171,172,174,175,176,177,178,179,180,181,182,183,184,186,187,188,189,190,191,192,193,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,212,213,214,215,217,218,220,221,222,223,224,225,226,227,228,229,232,233,236,237,238,240,241,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,260,261,262,263,264,265,267,268,269,270,271,272,273,275,276};

// Returns 0-based map of reaction order
void GET_RMAP(int * _rmap)
{
for (int j=0; j<277; ++j) {
_rmap[j] = rmap[j];
}
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
const int ns[277] =
     {4,4,4,3,2,2,3,3,3,3,3,4,4,4,3,3,2,2,4,4,4,4,4,3,4,4,4,4,3,4,4,4,4,4,5,3,4,3,4,3,3,3,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,5,3,4,4,3,4,4,4,5,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,3,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,3,2,4,4,4,4,4,4,3,4,4,3,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,3,3,4,5,5,4,4,4,5,3,3,4,4,4,4,4,4,3,4,4,4,4,2,4,4,3,4,4,4,3,4,4,4,4,4,4,2,3,4,4,3,3,4,4,4,3,3,4,3,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,3,3,4,4,4,2,4,4,3,4,3,4,4,4,3,3,2,4,3};
const int kiv[1385] =
     {14,11,15,16,0,12,15,14,16,0,12,16,14,18,0,18,15,16,0,0,12,14,0,0,0,15,11,0,0,0,14,15,16,0,0,14,16,18,0,0,14,11,17,0,0,14,11,17,0,0,14,17,16,0,0,12,11,14,17,0,17,15,11,16,0,17,16,18,11,0,17,19,11,0,0,17,19,11,0,0,19,16,0,0,0,19,16,0,0,0,14,19,18,16,0,14,19,12,17,0,19,15,17,16,0,19,16,18,17,0,19,16,18,17,0,5,15,6,0,0,5,11,6,15,0,5,16,6,14,0,5,16,6,14,0,5,17,6,16,0,20,5,14,0,0,20,11,5,17,0,14,20,5,12,0,20,15,5,16,0,20,15,6,14,0,20,16,5,18,0,20,17,6,14,16,20,5,12,0,0,21,11,20,17,0,20,11,9,0,0,21,9,20,8,0,10,16,8,0,0,6,14,10,0,0,20,21,5,0,0,1,32,14,0,0,1,32,14,0,0,1,7,0,0,0,32,12,14,1,0,1,15,14,33,0,1,15,32,16,0,1,15,5,31,0,1,16,32,18,0,1,16,14,35,0,1,16,14,2,0,1,16,5,30,0,1,11,32,17,0,32,1,14,34,0,14,7,14,1,0,7,15,5,31,0,7,16,14,2,0,32,7,14,34,0,32,15,5,28,0,32,16,14,33,0,32,11,33,15,0,32,11,33,15,0,32,11,5,3,0,32,6,5,33,0,21,32,1,20,0,32,3,28,33,0,32,22,33,3,0,32,22,5,4,0,32,22,6,48,0,32,24,1,3,0,32,26,1,22,0,32,4,37,3,0,32,2,1,33,0,32,33,5,37,0,2,5,31,0,0,14,2,5,30,0,14,2,12,33,0,2,15,33,16,0,2,15,6,31,0,2,15,5,24,0,2,16,18,33,0,2,17,19,33,0,2,11,6,24,0,2,31,33,30,0,14,35,14,2,0,14,35,5,30,0,14,35,12,33,0,35,15,33,16,0,35,16,18,33,0,35,30,33,23,0,36,1,15,0,0,14,36,1,16,0,36,15,20,3,0,36,16,21,3,0,36,15,33,16,0,36,16,5,12,3,36,14,33,0,0,36,16,18,33,0,36,16,20,24,0,33,5,28,0,0,14,33,5,31,0,33,15,5,3,0,33,16,5,25,0,33,16,5,14,3,17,33,2,11,0,33,11,6,3,0,20,33,5,2,0,21,33,20,2,0,33,3,5,4,0,33,3,6,48,0,33,22,5,3,0,33,22,6,4,0,24,33,2,3,0,26,33,2,22,0,28,33,5,48,0,33,23,2,30,0,33,5,48,0,0,43,14,1,0,0,14,43,12,1,0,43,15,1,16,0,43,16,18,1,0,43,16,18,1,0,43,11,21,3,0,43,30,1,23,0,44,14,1,0,0,14,44,14,43,0,14,44,12,1,0,44,15,14,2,0,44,15,1,16,0,44,16,18,1,0,23,14,30,0,0,14,23,12,30,0,23,15,30,16,0,23,16,18,30,0,17,23,19,30,0,14,30,12,31,0,30,15,14,24,0,30,15,31,16,0,30,16,18,31,0,17,30,27,16,0,17,30,23,11,0,30,11,27,15,0,30,11,24,16,0,30,31,23,0,0,31,30,28,23,0,28,30,14,48,0,30,3,18,48,0,30,3,29,16,0,30,22,18,4,0,30,22,27,3,0,14,31,12,28,0,31,15,14,3,0,31,16,14,24,0,31,16,18,28,0,31,11,24,15,0,31,11,3,16,0,31,14,48,0,0,28,31,14,48,0,31,3,14,4,0,31,3,48,16,0,31,22,4,16,0,27,14,24,0,0,27,45,0,0,0,14,27,12,24,0,14,27,30,16,0,27,15,24,16,0,27,16,18,24,0,27,17,19,24,0,27,11,24,17,0,27,3,24,0,0,27,30,24,23,0,27,22,24,26,0,45,14,24,0,0,14,45,30,16,0,14,45,12,24,0,45,15,24,16,0,45,15,24,16,0,45,16,18,24,0,45,17,19,24,0,45,11,24,17,0,45,30,24,23,0,45,22,24,26,0,24,11,17,3,0,20,24,21,3,0,25,14,3,0,0,14,25,14,24,0,14,25,31,16,0,25,15,3,16,0,25,16,14,26,0,25,11,22,16,0,20,3,5,24,0,14,26,24,16,0,14,26,18,3,0,22,15,46,0,0,22,3,46,0,0,5,22,6,3,0,20,22,5,3,16,20,22,6,14,3,14,46,22,16,0,46,15,22,11,0,46,16,17,22,0,17,46,22,11,16,46,3,11,0,0,29,14,48,0,0,14,29,12,48,0,29,15,14,4,0,29,15,48,16,0,29,15,31,3,0,29,16,18,48,0,29,11,17,48,0,29,14,48,0,0,31,29,48,30,0,30,29,48,23,0,29,3,24,48,0,5,4,6,48,0,34,32,0,0,0,34,15,32,33,0,34,16,32,35,0,14,15,0,0,0,18,0,18,16,0,12,0,12,16,0,48,0,48,16,0,16,0,16,0,0,14,0,14,16,0,13,0,13,16,0,11,0,11,16,0,6,0,6,16,0,5,0,5,16,0,14,4,48,0,0,48,28,0,0,0,28,15,3,0,0,48,15,28,3,0,17,3,22,16,0,3,16,26,0,0,14,3,24,0,0,12,22,14,26,0,14,22,3,16,0,22,15,3,11,0,22,3,15,0,0,22,3,11,0,0,17,22,26,11,0,22,16,47,0,0,3,22,4,11,0,26,15,22,16,0,26,16,18,22,0,26,18,3,22,0,24,15,3,16,0,24,16,18,3,0,14,24,12,3,0,24,22,26,3,0,24,18,4,0,0,24,3,4,16,0,28,22,4,15,0,28,11,3,15,0,28,16,14,3,0,28,4,48,3,0,4,3,48,22,0,27,12,3,0,0,17,3,47,0,0,14,39,12,22,0,39,15,22,16,0,39,16,18,22,0,39,26,0,0,0,26,30,23,22,0,24,28,31,3,0,4,48,15,0,0,4,15,48,11,0,4,15,3,0,0,14,4,48,16,0,14,4,48,16,0,4,16,17,48,0,46,3,11,0,0,46,22,11,0,0,40,22,0,0,0,40,15,41,11,0,41,15,22,0,0};
const int nuv[1385] =
     {-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,2,0,0,0,-2,1,0,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-2,1,1,0,0,-1,2,0,0,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-2,2,1,0,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,0,0,-2,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,2,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,2,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,2,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,2,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,-1,1,1,1,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,1,-1,1,1,0,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,2,0,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-2,2,1,0,0,-1,-1,1,1,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-2,1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,0,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-1,-1,1,1,0,-1,-1,2,0,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,-1,1,1,0,-1,1,1,0,0,-2,2,1,0,0,-1,2,0,0,0,-1,-1,1,1,0,-1,-1,2,0,0};
if (*i < 1) {
// Return max num species per reaction
*nspec = 5;
} else {
if (*i > 277) {
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
amrex::Real c[54]; // temporary storage
amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (id = 0; id < 54; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, *T);

// convert to chemkin units
for (id = 0; id < 277; ++id) {
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
amrex::Real g_RT[54];
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
for (id = 0; id < kd * 54; ++ id) {
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

// HO2CHO
ncf[ 8 * kd + 0 ] = 1; // C
ncf[ 8 * kd + 1 ] = 2; // H
ncf[ 8 * kd + 2 ] = 3; // O

// O2CHO
ncf[ 9 * kd + 0 ] = 1; // C
ncf[ 9 * kd + 1 ] = 1; // H
ncf[ 9 * kd + 2 ] = 3; // O

// OCHO
ncf[ 10 * kd + 0 ] = 1; // C
ncf[ 10 * kd + 1 ] = 1; // H
ncf[ 10 * kd + 2 ] = 2; // O

// O2
ncf[ 11 * kd + 2 ] = 2; // O

// H2
ncf[ 12 * kd + 1 ] = 2; // H

// AR
ncf[ 13 * kd + 4 ] = 1; // Ar

// H
ncf[ 14 * kd + 1 ] = 1; // H

// O
ncf[ 15 * kd + 2 ] = 1; // O

// OH
ncf[ 16 * kd + 1 ] = 1; // H
ncf[ 16 * kd + 2 ] = 1; // O

// HO2
ncf[ 17 * kd + 1 ] = 1; // H
ncf[ 17 * kd + 2 ] = 2; // O

// H2O
ncf[ 18 * kd + 1 ] = 2; // H
ncf[ 18 * kd + 2 ] = 1; // O

// H2O2
ncf[ 19 * kd + 1 ] = 2; // H
ncf[ 19 * kd + 2 ] = 2; // O

// HCO
ncf[ 20 * kd + 0 ] = 1; // C
ncf[ 20 * kd + 1 ] = 1; // H
ncf[ 20 * kd + 2 ] = 1; // O

// CH2O
ncf[ 21 * kd + 0 ] = 1; // C
ncf[ 21 * kd + 1 ] = 2; // H
ncf[ 21 * kd + 2 ] = 1; // O

// NO2
ncf[ 22 * kd + 3 ] = 1; // N
ncf[ 22 * kd + 2 ] = 2; // O

// NH3
ncf[ 23 * kd + 1 ] = 3; // H
ncf[ 23 * kd + 3 ] = 1; // N

// HNO
ncf[ 24 * kd + 1 ] = 1; // H
ncf[ 24 * kd + 3 ] = 1; // N
ncf[ 24 * kd + 2 ] = 1; // O

// HON
ncf[ 25 * kd + 1 ] = 1; // H
ncf[ 25 * kd + 3 ] = 1; // N
ncf[ 25 * kd + 2 ] = 1; // O

// HONO
ncf[ 26 * kd + 1 ] = 1; // H
ncf[ 26 * kd + 3 ] = 1; // N
ncf[ 26 * kd + 2 ] = 2; // O

// H2NO
ncf[ 27 * kd + 1 ] = 2; // H
ncf[ 27 * kd + 3 ] = 1; // N
ncf[ 27 * kd + 2 ] = 1; // O

// N
ncf[ 28 * kd + 3 ] = 1; // N

// NNH
ncf[ 29 * kd + 1 ] = 1; // H
ncf[ 29 * kd + 3 ] = 2; // N

// NH2
ncf[ 30 * kd + 1 ] = 2; // H
ncf[ 30 * kd + 3 ] = 1; // N

// NH
ncf[ 31 * kd + 1 ] = 1; // H
ncf[ 31 * kd + 3 ] = 1; // N

// CN
ncf[ 32 * kd + 0 ] = 1; // C
ncf[ 32 * kd + 3 ] = 1; // N

// NCO
ncf[ 33 * kd + 0 ] = 1; // C
ncf[ 33 * kd + 3 ] = 1; // N
ncf[ 33 * kd + 2 ] = 1; // O

// NCCN
ncf[ 34 * kd + 0 ] = 2; // C
ncf[ 34 * kd + 3 ] = 2; // N

// HOCN
ncf[ 35 * kd + 0 ] = 1; // C
ncf[ 35 * kd + 1 ] = 1; // H
ncf[ 35 * kd + 3 ] = 1; // N
ncf[ 35 * kd + 2 ] = 1; // O

// HCNO
ncf[ 36 * kd + 0 ] = 1; // C
ncf[ 36 * kd + 1 ] = 1; // H
ncf[ 36 * kd + 3 ] = 1; // N
ncf[ 36 * kd + 2 ] = 1; // O

// NCN
ncf[ 37 * kd + 0 ] = 1; // C
ncf[ 37 * kd + 3 ] = 2; // N

// HE
ncf[ 38 * kd + 5 ] = 1; // He

// HNO2
ncf[ 39 * kd + 1 ] = 1; // H
ncf[ 39 * kd + 3 ] = 1; // N
ncf[ 39 * kd + 2 ] = 2; // O

// N2O4
ncf[ 40 * kd + 3 ] = 2; // N
ncf[ 40 * kd + 2 ] = 4; // O

// N2O3
ncf[ 41 * kd + 3 ] = 2; // N
ncf[ 41 * kd + 2 ] = 3; // O

// CH
ncf[ 42 * kd + 0 ] = 1; // C
ncf[ 42 * kd + 1 ] = 1; // H

// H2CN
ncf[ 43 * kd + 0 ] = 1; // C
ncf[ 43 * kd + 1 ] = 2; // H
ncf[ 43 * kd + 3 ] = 1; // N

// HCNH
ncf[ 44 * kd + 0 ] = 1; // C
ncf[ 44 * kd + 1 ] = 2; // H
ncf[ 44 * kd + 3 ] = 1; // N

// HNOH
ncf[ 45 * kd + 1 ] = 2; // H
ncf[ 45 * kd + 3 ] = 1; // N
ncf[ 45 * kd + 2 ] = 1; // O

// NO3
ncf[ 46 * kd + 3 ] = 1; // N
ncf[ 46 * kd + 2 ] = 3; // O

// HONO2
ncf[ 47 * kd + 1 ] = 1; // H
ncf[ 47 * kd + 3 ] = 1; // N
ncf[ 47 * kd + 2 ] = 3; // O

// N2
ncf[ 48 * kd + 3 ] = 2; // N

// CH4
ncf[ 49 * kd + 0 ] = 1; // C
ncf[ 49 * kd + 1 ] = 4; // H

// C2H6
ncf[ 50 * kd + 0 ] = 2; // C
ncf[ 50 * kd + 1 ] = 6; // H

// N2H4
ncf[ 51 * kd + 1 ] = 4; // H
ncf[ 51 * kd + 3 ] = 2; // N

// N2H3
ncf[ 52 * kd + 1 ] = 3; // H
ncf[ 52 * kd + 3 ] = 2; // N

// N2H2
ncf[ 53 * kd + 1 ] = 2; // H
ncf[ 53 * kd + 3 ] = 2; // N

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
kname.resize(54);
kname[0] = "OH*";
kname[1] = "HCN";
kname[2] = "HNCO";
kname[3] = "NO";
kname[4] = "N2O";
kname[5] = "CO";
kname[6] = "CO2";
kname[7] = "HNC";
kname[8] = "HO2CHO";
kname[9] = "O2CHO";
kname[10] = "OCHO";
kname[11] = "O2";
kname[12] = "H2";
kname[13] = "AR";
kname[14] = "H";
kname[15] = "O";
kname[16] = "OH";
kname[17] = "HO2";
kname[18] = "H2O";
kname[19] = "H2O2";
kname[20] = "HCO";
kname[21] = "CH2O";
kname[22] = "NO2";
kname[23] = "NH3";
kname[24] = "HNO";
kname[25] = "HON";
kname[26] = "HONO";
kname[27] = "H2NO";
kname[28] = "N";
kname[29] = "NNH";
kname[30] = "NH2";
kname[31] = "NH";
kname[32] = "CN";
kname[33] = "NCO";
kname[34] = "NCCN";
kname[35] = "HOCN";
kname[36] = "HCNO";
kname[37] = "NCN";
kname[38] = "HE";
kname[39] = "HNO2";
kname[40] = "N2O4";
kname[41] = "N2O3";
kname[42] = "CH";
kname[43] = "H2CN";
kname[44] = "HCNH";
kname[45] = "HNOH";
kname[46] = "NO3";
kname[47] = "HONO2";
kname[48] = "N2";
kname[49] = "CH4";
kname[50] = "C2H6";
kname[51] = "N2H4";
kname[52] = "N2H3";
kname[53] = "N2H2";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(Jac[ 55 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 55 * k + l] != 0.0){
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
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 55 * k + l] != 0.0){
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
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 55;
int offset_col = nc * 55;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if(Jac[55*k + l] != 0.0) {
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
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if(Jac[55*k + l] != 0.0) {
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
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if(Jac[55*k + l] != 0.0) {
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
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
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
int offset = nc * 55;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
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
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<55; k++) {
for (int l=0; l<55; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 55*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 55*k + l;
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
amrex::GpuArray<amrex::Real,3025> Jac = {0.0};
amrex::GpuArray<amrex::Real,54> conc = {0.0};
for (int n=0; n<54; n++) {
    conc[n] = 1.0/ 54.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
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
for (int l=0; l<55; l++) {
for (int k=0; k<55; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[55*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
