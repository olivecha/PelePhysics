#include "mechanism.H"
const int rmap[64] = {49,56,22,26,29,17,34,45,46,47,48,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,23,24,25,27,28,30,31,32,33,35,36,37,38,39,40,41,42,43,44,50,51,52,53,54,55,57,58,59,60,61,62,63};

// Returns 0-based map of reaction order
void GET_RMAP(int * _rmap)
{
for (int j=0; j<64; ++j) {
_rmap[j] = rmap[j];
}
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void CKINU(int * i, int * nspec, int * ki, int * nu)
{
const int ns[64] =
     {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,3,4,4,4,3,4,4,3,4,4,3,4,3,4,4,4,4,4,4,4,4,4,3,2,3,2,3,3,3,4,4,4,4,4,2,3,3,4,4,4,4,4};
const int kiv[256] =
     {18,6,14,0,14,2,0,6,14,7,5,0,5,17,3,14,17,6,5,0,17,7,5,13,17,7,9,14,17,2,13,6,17,0,5,1,17,0,18,7,5,16,3,17,16,6,5,13,16,2,19,6,16,7,9,17,14,16,5,18,16,0,9,18,16,0,15,7,12,5,16,0,5,12,3,16,12,6,16,7,12,7,9,16,15,2,8,18,15,5,18,0,5,15,3,18,15,6,5,1,15,7,9,18,5,0,13,0,5,13,3,0,13,7,9,0,1,18,6,0,5,1,18,7,5,1,18,7,1,6,0,0,1,7,8,18,11,0,6,0,8,0,11,7,5,11,0,7,11,6,0,2,19,6,13,7,19,2,13,8,19,8,10,13,5,2,6,7,3,6,5,7,3,7,5,9,9,6,7,0,5,3,0,0,5,7,9,0,6,2,0,0,5,6,7,0,5,2,8,0,5,8,7,0,5,8,3,2,5,8,9,6,8,6,2,7,8,7,9,2,8,7,9,2,7,10,0,0,8,10,2,0,8,10,2,0,5,10,3,8,5,10,9,7,10,7,9,8,10,7,9,8,10,6,8,7};
const int nuv[256] =
     {-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,2,1,-1,-1,1,1,-1,-1,1,1,-1,1,1,0,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,1,1,0,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,0,-1,-1,1,1,-1,-1,1,1,-1,1,1,0,-1,-1,1,1,-1,-1,1,1,-1,-1,2,0,-1,-1,1,1,-1,1,1,0,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,2,0,-2,1,0,0,-1,-1,1,0,-2,1,0,0,-1,-1,1,0,-1,-1,1,0,-1,-1,2,0,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-2,1,0,0,-2,1,1,0,-2,1,1,0,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1};
if (*i < 1) {
// Return max num species per reaction
*nspec = 4;
} else {
if (*i > 64) {
*nspec = -1;
} else {
*nspec = ns[*i-1];
for (int j=0; j<*nspec; ++j) {
ki[j] = kiv[(*i-1)*4 + j] + 1;
nu[j] = nuv[(*i-1)*4 + j];
}
}
}
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void CKKFKR(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  q_f, amrex::Real *  q_r)
{
int id; // loop counter
amrex::Real c[21]; // temporary storage
amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

// Compute conversion, see Eq 10
for (id = 0; id < 21; ++id) {
c[id] = x[id]*PORT;
}

// convert to chemkin units
progressRateFR(q_f, q_r, c, *T);

// convert to chemkin units
for (id = 0; id < 64; ++id) {
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
amrex::Real g_RT[21];
gibbs(g_RT, tc);

amrex::Real sc_qss[1];
amrex::Real kf_qss[0], qf_qss[0], qr_qss[0];
comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

return;
}

// save atomic weights into array
void atomicWeight(amrex::Real *  awt)
{
awt[0] = 14.007000; // N
awt[1] = 15.999000; // O
awt[2] = 1.008000; // H
awt[3] = 39.950000; // Ar
awt[4] = 4.002602; // He
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
int kd = 5; 
// Zero ncf
for (id = 0; id < kd * 21; ++ id) {
 ncf[id] = 0; 
}

// NO
ncf[ 0 * kd + 0 ] = 1; // N
ncf[ 0 * kd + 1 ] = 1; // O

// N2O
ncf[ 1 * kd + 0 ] = 2; // N
ncf[ 1 * kd + 1 ] = 1; // O

// O2
ncf[ 2 * kd + 1 ] = 2; // O

// H2
ncf[ 3 * kd + 2 ] = 2; // H

// AR
ncf[ 4 * kd + 3 ] = 1; // Ar

// H
ncf[ 5 * kd + 2 ] = 1; // H

// O
ncf[ 6 * kd + 1 ] = 1; // O

// OH
ncf[ 7 * kd + 2 ] = 1; // H
ncf[ 7 * kd + 1 ] = 1; // O

// HO2
ncf[ 8 * kd + 2 ] = 1; // H
ncf[ 8 * kd + 1 ] = 2; // O

// H2O
ncf[ 9 * kd + 2 ] = 2; // H
ncf[ 9 * kd + 1 ] = 1; // O

// H2O2
ncf[ 10 * kd + 2 ] = 2; // H
ncf[ 10 * kd + 1 ] = 2; // O

// NO2
ncf[ 11 * kd + 0 ] = 1; // N
ncf[ 11 * kd + 1 ] = 2; // O

// NH3
ncf[ 12 * kd + 2 ] = 3; // H
ncf[ 12 * kd + 0 ] = 1; // N

// HNO
ncf[ 13 * kd + 2 ] = 1; // H
ncf[ 13 * kd + 0 ] = 1; // N
ncf[ 13 * kd + 1 ] = 1; // O

// N
ncf[ 14 * kd + 0 ] = 1; // N

// N2H
ncf[ 15 * kd + 2 ] = 1; // H
ncf[ 15 * kd + 0 ] = 2; // N

// NH2
ncf[ 16 * kd + 2 ] = 2; // H
ncf[ 16 * kd + 0 ] = 1; // N

// NH
ncf[ 17 * kd + 2 ] = 1; // H
ncf[ 17 * kd + 0 ] = 1; // N

// N2
ncf[ 18 * kd + 0 ] = 2; // N

// H2NO
ncf[ 19 * kd + 2 ] = 2; // H
ncf[ 19 * kd + 0 ] = 1; // N
ncf[ 19 * kd + 1 ] = 1; // O

// HE
ncf[ 20 * kd + 4 ] = 1; // He

}

// Returns the vector of strings of element names
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
ename.resize(5);
ename[0] = "N";
ename[1] = "O";
ename[2] = "H";
ename[3] = "Ar";
ename[4] = "He";
}

// Returns the vector of strings of species names
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
kname.resize(21);
kname[0] = "NO";
kname[1] = "N2O";
kname[2] = "O2";
kname[3] = "H2";
kname[4] = "AR";
kname[5] = "H";
kname[6] = "O";
kname[7] = "OH";
kname[8] = "HO2";
kname[9] = "H2O";
kname[10] = "H2O2";
kname[11] = "NO2";
kname[12] = "NH3";
kname[13] = "HNO";
kname[14] = "N";
kname[15] = "N2H";
kname[16] = "NH2";
kname[17] = "NH";
kname[18] = "N2";
kname[19] = "H2NO";
kname[20] = "HE";
}

// compute the sparsity pattern of the chemistry Jacobian
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<22; k++) {
for (int l=0; l<22; l++) {
if(Jac[ 22 * k + l] != 0.0){
nJdata_tmp = nJdata_tmp + 1;
}
}
}

*nJdata = NCELLS * nJdata_tmp;
}



// compute the sparsity pattern of the system Jacobian
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<22; k++) {
for (int l=0; l<22; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 22 * k + l] != 0.0){
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
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

int nJdata_tmp = 0;
for (int k=0; k<22; k++) {
for (int l=0; l<22; l++) {
if(k == l){
nJdata_tmp = nJdata_tmp + 1;
} else {
if(Jac[ 22 * k + l] != 0.0){
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
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int nc=0; nc<NCELLS; nc++) {
int offset_row = nc * 22;
int offset_col = nc * 22;
for (int k=0; k<22; k++) {
for (int l=0; l<22; l++) {
if(Jac[22*k + l] != 0.0) {
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
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtrs[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 22;
for (int l=0; l<22; l++) {
for (int k=0; k<22; k++) {
if(Jac[22*k + l] != 0.0) {
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
int offset = nc * 22;
for (int l=0; l<22; l++) {
for (int k=0; k<22; k++) {
if(Jac[22*k + l] != 0.0) {
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
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int nc=0; nc<NCELLS; nc++) {
int offset = nc * 22;
for (int l=0; l<22; l++) {
for (int k=0; k<22; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1 + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[22*k + l] != 0.0) {
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
int offset = nc * 22;
for (int l=0; l<22; l++) {
for (int k=0; k<22; k++) {
if (k == l) {
colVals[nJdata_tmp] = l + offset; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[22*k + l] != 0.0) {
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
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

colPtrs[0] = 0;
int nJdata_tmp = 0;
for (int k=0; k<22; k++) {
for (int l=0; l<22; l++) {
if (k == l) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 22*k + l;
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[22*k + l] != 0.0) {
rowVals[nJdata_tmp] = l; 
indx[nJdata_tmp] = 22*k + l;
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
amrex::GpuArray<amrex::Real,484> Jac = {0.0};
amrex::GpuArray<amrex::Real,21> conc = {0.0};
for (int n=0; n<21; n++) {
    conc[n] = 1.0/ 21.000000 ;
}
aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

if (base == 1) {
rowPtr[0] = 1;
int nJdata_tmp = 1;
for (int l=0; l<22; l++) {
for (int k=0; k<22; k++) {
if (k == l) {
colVals[nJdata_tmp-1] = l+1; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[22*k + l] != 0.0) {
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
for (int l=0; l<22; l++) {
for (int k=0; k<22; k++) {
if (k == l) {
colVals[nJdata_tmp] = l; 
nJdata_tmp = nJdata_tmp + 1; 
} else {
if(Jac[22*k + l] != 0.0) {
colVals[nJdata_tmp] = k; 
nJdata_tmp = nJdata_tmp + 1; 
}
}
}
rowPtr[l+1] = nJdata_tmp;
}
}
}
