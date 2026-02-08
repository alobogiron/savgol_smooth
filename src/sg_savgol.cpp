#include "sg_savgol.h"

// ===== TABELA: polyorder 2 ou 3 (mesmos coeficientes para smoothing) =====
// Formato: { denom, w(|k|=M), ..., w(|k|=0), padding... }
static constexpr int32_t SG_P23_3_21[10][14] = {
  {    1,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 }, // 3
  {   35,   -3,   12,   17,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 }, // 5
  {   21,   -2,    3,    6,    7,    0,    0,    0,    0,    0,    0,    0,    0,    0 }, // 7
  {  231,  -21,   14,   39,   54,   59,    0,    0,    0,    0,    0,    0,    0,    0 }, // 9
  {  429,  -36,    9,   44,   69,   84,   89,    0,    0,    0,    0,    0,    0,    0 }, // 11
  {  143,  -11,    0,    9,   16,   21,   24,   25,    0,    0,    0,    0,    0,    0 }, // 13
  { 1105,  -78,  -13,   42,   87,  122,  147,  162,  167,    0,    0,    0,    0,    0 }, // 15
  {  323,  -21,   -6,    7,   18,   27,   34,   39,   42,   43,    0,    0,    0,    0 }, // 17
  { 2261, -136,  -51,   24,   89,  144,  189,  224,  249,  264,  269,    0,    0,    0 }, // 19
  { 3059, -171,  -76,    9,   84,  149,  204,  249,  284,  309,  324,  329,    0,    0 }, // 21
};

// ===== TABELA: polyorder 4 ou 5 (mesmos coeficientes para smoothing) =====
static constexpr int32_t SG_P45_3_21[10][14] = {
  {    1,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 }, // 3
  {    1,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0 }, // 5
  {  231,    5,  -30,   75,  131,    0,    0,    0,    0,    0,    0,    0,    0,    0 }, // 7
  {  429,   15,  -55,   30,  135,  179,    0,    0,    0,    0,    0,    0,    0,    0 }, // 9
  {  429,   18,  -45,  -10,   60,  120,  143,    0,    0,    0,    0,    0,    0,    0 }, // 11
  { 2431,  110, -198, -135,  110,  390,  600,  677,    0,    0,    0,    0,    0,    0 }, // 13
  { 46189, 2145, -2860, -2937, -165, 3755, 7500, 10125, 11063,    0,    0,    0,    0,    0 }, // 15
  { 4199,  195, -195, -260, -117,  135,  415,  660,  825,  883,    0,    0,    0,    0 }, // 17
  { 7429,  340, -255, -420, -290,   18,  405,  790, 1110, 1320, 1393,    0,    0,    0 }, // 19
  { 260015, 11628, -6460, -13005, -11220, -3940, 6378, 17655, 28190, 36660, 42120, 44003,    0,    0 }, // 21
};

static inline int clampi(int v, int lo, int hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

static inline int adjust_window(int L, int window) {
    if (L <= 0) return 0;
    if (window > L) window = L;
    if ((window & 1) == 0) window -= 1;
    if (window < 3) window = 3;
    if (window > 21) window = 21;
    return window;
}

bool sg_savgol_smooth_nearest(
    const float* x,
    float* y,
    size_t L,
    int window,
    int polyorder
) {
    if (!x || !y || L == 0) return false;

    window = adjust_window((int)L, window);
    if (window < 3 || window > 21 || (window & 1) == 0) return false;

    if (polyorder < 2) polyorder = 2;
    if (polyorder > 5) polyorder = 5;

    if (window < 4) polyorder = 2;

    const int row = (window - 3) / 2;
    const int32_t* c = (polyorder <= 3) ? SG_P23_3_21[row] : SG_P45_3_21[row];

    const int M = window / 2;
    const int last = (int)L - 1;
    const float denom = (float)c[0];

    for (int n = 0; n < (int)L; n++) {
        float acc = (float)c[1 + M] * x[n];

        for (int k = 1; k <= M; k++) {
            const float hk = (float)c[1 + (M - k)];
            const int im = clampi(n - k, 0, last);
            const int ip = clampi(n + k, 0, last);
            acc += hk * (x[im] + x[ip]);
        }

        y[n] = acc / denom;
    }

    return true;
}
