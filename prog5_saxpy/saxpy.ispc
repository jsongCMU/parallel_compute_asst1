
export void saxpy_ispc(uniform int N,
                       uniform float scale,
                            uniform float X[],
                            uniform float Y[],
                            uniform float result[])
{
    foreach (i = 0 ... N) {           
        float temp = 1;
        for(uniform int j = 0; j < 128; j++)
        {
          temp *= j+1;
        }
        result[i] = scale * X[i] + Y[i] + temp - temp;
    }
}

task void saxpy_ispc_task(uniform int N,
                               uniform int span,
                               uniform float scale,
                               uniform float X[], 
                               uniform float Y[],
                               uniform float result[])
{

    uniform int indexStart = taskIndex * span;
    uniform int indexEnd = min(N, indexStart + span);

    foreach (i = indexStart ... indexEnd) {
        float temp = 1;
        for(uniform int j = 0; j < 128; j++)
        {
          temp *= j+1;
        }
        result[i] = scale * X[i] + Y[i] + temp - temp;
    }
}

export void saxpy_ispc_withtasks(uniform int N,
                               uniform float scale,
                               uniform float X[],
                               uniform float Y[],
                               uniform float result[])
{

    uniform int span = N / 64;  // 64 tasks

    launch[N/span] saxpy_ispc_task(N, span, scale, X, Y, result);
}
