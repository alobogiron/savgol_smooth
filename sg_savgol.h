#pragma once
#include <cstdint>
#include <cstddef>

// SG Savitzky-Golay smoothing (deriv=0), mode=nearest (clamp)
// - x: entrada int32_t
// - y: saída float
// - L: tamanho do vetor
// - window: ímpar 3..21
// - polyorder: 2/3 usa tabela P23; 4/5 usa tabela P45
// Retorna false se parâmetros inválidos.
bool sg_savgol_smooth_nearest(
    const int32_t* x,
    float* y,
    size_t L,
    int window,
    int polyorder
);

