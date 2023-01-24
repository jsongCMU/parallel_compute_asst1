#include <stdio.h>
#include <algorithm>
#include <math.h>
#include "CMU418intrin.h"
#include "logger.h"
using namespace std;


void absSerial(float* values, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	if (x < 0) {
	    output[i] = -x;
	} else {
	    output[i] = x;
	}
    }
}

// implementation of absolute value using 15418 instrinsics
void absVector(float* values, float* output, int N) {
    __cmu418_vec_float x;
    __cmu418_vec_float result;
    __cmu418_vec_float zero = _cmu418_vset_float(0.f);
    __cmu418_mask maskAll, maskIsNegative, maskIsNotNegative;

    //  Note: Take a careful look at this loop indexing.  This example
    //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
    //  Why is that the case?
    for (int i=0; i<N; i+=VECTOR_WIDTH) {

	// All ones
	maskAll = _cmu418_init_ones();

	// All zeros
	maskIsNegative = _cmu418_init_ones(0);

	// Load vector of values from contiguous memory addresses
	_cmu418_vload_float(x, values+i, maskAll);               // x = values[i];

	// Set mask according to predicate
	_cmu418_vlt_float(maskIsNegative, x, zero, maskAll);     // if (x < 0) {

	// Execute instruction using mask ("if" clause)
	_cmu418_vsub_float(result, zero, x, maskIsNegative);      //   output[i] = -x;

	// Inverse maskIsNegative to generate "else" mask
	maskIsNotNegative = _cmu418_mask_not(maskIsNegative);     // } else {

	// Execute instruction ("else" clause)
	_cmu418_vload_float(result, values+i, maskIsNotNegative); //   output[i] = x; }

	// Write results back to memory
	_cmu418_vstore_float(output+i, result, maskAll);
    }
}

// Accepts an array of values and an array of exponents
// For each element, compute values[i]^exponents[i] and clamp value to
// 4.18.  Store result in outputs.
// Uses iterative squaring, so that total iterations is proportional
// to the log_2 of the exponent
void clampedExpSerial(float* values, int* exponents, float* output, int N) {
    for (int i=0; i<N; i++) {
	float x = values[i];
	float result = 1.f;
	int y = exponents[i];
	float xpower = x;
	while (y > 0) {
	    if (y & 0x1) {
			result *= xpower;
		}
	    xpower = xpower * xpower;
	    y >>= 1;
	}
	if (result > 4.18f) {
	    result = 4.18f;
	}
	output[i] = result;
    }
}

void clampedExpVector(float* values, int* exponents, float* output, int N) {
    __cmu418_vec_float x, result, xpower, vectorThreshold;
    __cmu418_vec_int y, vectorZero, vectorOne, vector_is_odd;
    __cmu418_mask maskAll, mask_yGTZ, mask_is_odd, mask_clamp;
    vectorZero = _cmu418_vset_int(0);
    vectorOne = _cmu418_vset_int(1);
    vectorThreshold = _cmu418_vset_float(4.18);
    for (int i=0; i<N; i+=VECTOR_WIDTH) {
        maskAll = _cmu418_init_ones(N-i);
        _cmu418_vload_float(x, values+i, maskAll); // x = values[i];
        _cmu418_vset_float(result, 1.0, maskAll); // result = 1.f;
        _cmu418_vload_int(y, exponents+i, maskAll); // y = exponents[i];
        _cmu418_vmove_float(xpower, x, maskAll); // xpower = x;
        _cmu418_vgt_int(mask_yGTZ, y, vectorZero, maskAll); // (y>0)
        while(_cmu418_cntbits(mask_yGTZ)) // while(y>0)
        {
            _cmu418_vbitand_int(vector_is_odd, y, vectorOne, maskAll); // y & 0x1
            _cmu418_vgt_int(mask_is_odd, vector_is_odd, vectorZero, maskAll); // Convert to mask
            _cmu418_vmult_float(result, result, xpower, mask_is_odd); // if(...) result *= xpower;
            _cmu418_vmult_float(xpower, xpower, xpower, maskAll); // xpower = xpower * xpower;
            _cmu418_vshiftright_int(y, y, vectorOne, maskAll); // y >>= 1;
            _cmu418_vgt_int(mask_yGTZ, y, vectorZero, maskAll); // (y>0)
        }
        _cmu418_vgt_float(mask_clamp, result, vectorThreshold, maskAll); // result > 4.18f
        _cmu418_vset_float(result, 4.18, mask_clamp); // if(...) result = 4.18f;
        _cmu418_vstore_float(output+i, result, maskAll); // output[i] = result;
    }
}


float arraySumSerial(float* values, int N) {
    float sum = 0;
    for (int i=0; i<N; i++) {
	sum += values[i];
    }

    return sum;
}

// Assume N % VECTOR_WIDTH == 0
// Assume VECTOR_WIDTH is a power of 2
float arraySumVector(float* values, int N) {
    __cmu418_vec_float sum;
    __cmu418_mask maskAll;
    maskAll = _cmu418_init_ones();
    _cmu418_vset_float(sum, 0.0, maskAll); // float sum = 0;
    for (int i=0; i<N; i+=VECTOR_WIDTH) { // for (int i=0; i<N; i++)
      __cmu418_vec_float temp_val;
      _cmu418_vload_float(temp_val, values+i, maskAll);
      _cmu418_vadd_float(sum, sum, temp_val, maskAll); // sum += values[i];
    }
    int counter = VECTOR_WIDTH;
    while(counter >>= 1){
      _cmu418_hadd_float(sum, sum);
      _cmu418_interleave_float(sum, sum);
    }
	return sum.value[0];
}
