#include <smmintrin.h> // For _mm_stream_load_si128
#include <emmintrin.h> // For _mm_mul_ps
#include <assert.h>
#include <stdint.h>

extern void saxpySerial(int N,
			float scale,
			float X[],
			float Y[],
			float result[]);


void saxpyStreaming(int N,
                    float scale,
                    float X[],
                    float Y[],
                    float result[])
{
    float scale_array[] = {scale, scale, scale, scale};
    // Float is 4 bytes / 32 bits
    // 128 bit operations are 4x larger, so can operate on 4 floats at once
    // Therefore only need a quarter of the usual number of iterations
    X-=4;
    Y-=4;
    result-=4;
    __m128 scale_4 = _mm_castsi128_ps(_mm_stream_load_si128((__m128i *) scale_array));
    for(int i = 0; i < N/4; i++)
    {
      // Load
      __m128 x_4 = _mm_castsi128_ps(_mm_stream_load_si128((__m128i *) (X+=4)));
      __m128 y_4 = _mm_castsi128_ps(_mm_stream_load_si128((__m128i *) (Y+=4)));
      // Compute
      __m128 res_4 = _mm_mul_ps(scale_4, x_4);
      res_4 = _mm_add_ps(res_4, y_4);
      // Store
      _mm_stream_ps((result+=4), res_4);
    }
}

