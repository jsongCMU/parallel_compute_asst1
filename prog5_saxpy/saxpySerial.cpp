
void saxpySerial(int N,
                       float scale,
                       float X[],
                       float Y[],
                       float result[])
{

    for (int i=0; i<N; i++) {
        float temp = 1;
        for(int j = 0; j < 512; j++)
        {
          temp *= j+1;
        }
        result[i] = scale*temp + X[i] + Y[i];
    }
}

