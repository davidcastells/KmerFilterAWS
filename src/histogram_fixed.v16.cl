#define KMER_CHUNKS     ((KMER_BINS)/32)



void readBigEndian512bits(__global unsigned char* restrict p, 
unsigned int* v15, unsigned int* v14, unsigned int* v13, unsigned int* v12,
unsigned int* v11, unsigned int* v10, unsigned int* v9, unsigned int* v8,
unsigned int* v7, unsigned int* v6, unsigned int* v5, unsigned int* v4,
unsigned int* v3, unsigned int* v2, unsigned int* v1, unsigned int* v0)
{
    *v15 = (p[0] << 24) | (p[1] << 16) | (p[2] << 8) | p[3];
    *v14 = (p[4] << 24) | (p[5] << 16) | (p[6] << 8) | p[7];
    *v13 = (p[8] << 24) | (p[9] << 16) | (p[10] << 8) | p[11];
    *v12 = (p[12] << 24) | (p[13] << 16) | (p[14] << 8) | p[15];
    *v11 = (p[16] << 24) | (p[17] << 16) | (p[18] << 8) | p[19];
    *v10 = (p[20] << 24) | (p[21] << 16) | (p[22] << 8) | p[23];
    *v9 = (p[24] << 24) | (p[25] << 16) | (p[26] << 8) | p[27];
    *v8 = (p[28] << 24) | (p[29] << 16) | (p[30] << 8) | p[31];
    *v7 = (p[32] << 24) | (p[33] << 16) | (p[34] << 8) | p[35];
    *v6 = (p[36] << 24) | (p[37] << 16) | (p[38] << 8) | p[39];
    *v5 = (p[40] << 24) | (p[41] << 16) | (p[42] << 8) | p[43];
    *v4 = (p[44] << 24) | (p[45] << 16) | (p[46] << 8) | p[47];
    *v3 = (p[48] << 24) | (p[49] << 16) | (p[50] << 8) | p[51];
    *v2 = (p[52] << 24) | (p[53] << 16) | (p[54] << 8) | p[55];
    *v1 = (p[56] << 24) | (p[57] << 16) | (p[58] << 8) | p[59];
    *v0 = (p[60] << 24) | (p[61] << 16) | (p[62] << 8) | p[63];
}



unsigned int get_kmer_index(unsigned int v15, unsigned int v14, unsigned int v13,  unsigned int v12,
unsigned int v11, unsigned int v10, unsigned int v9, unsigned int v8,
unsigned int v7, unsigned int v6, unsigned int v5, unsigned int v4,
unsigned int v3, unsigned int v2, unsigned int v1, unsigned int v0, int i)
{
unsigned int RET_MASK=(1<<10) -1;

unsigned int ret;
int b0 = 512  - (i)*2 -1 ; 		// number of the higher bit where the kmer begins
int b1 = 512 - (KMER_K+i)*2;	// number of the lower bit where the kmer begins
int wb0 = b0 % 32;		// b0 in 32 bits word
int wb1 = b1 % 32;		// b1 in 32 bits word
int k0 = b0 / 32;			// word that contains the b0
int k1 = b1 / 32;			// word that contains the b1
	unsigned int r0; unsigned int r1;

switch (k0)
{
case 15: r0 = v15; break;
case 14: r0 = v14; break;
case 13: r0 = v13; break;
case 12: r0 = v12; break;
case 11: r0 = v11; break;
case 10: r0 = v10; break;
case 9: r0 = v9; break;
case 8: r0 = v8; break;
case 7: r0 = v7; break;
case 6: r0 = v6; break;
case 5: r0 = v5; break;
case 4: r0 = v4; break;
case 3: r0 = v3; break;
case 2: r0 = v2; break;
case 1: r0 = v1; break;
case 0: r0 = v0; break;
}

switch (k1)
{
case 15: r1 = v15; break;
case 14: r1 = v14; break;
case 13: r1 = v13; break;
case 12: r1 = v12; break;
case 11: r1 = v11; break;
case 10: r1 = v10; break;
case 9: r1 = v9; break;
case 8: r1 = v8; break;
case 7: r1 = v7; break;
case 6: r1 = v6; break;
case 5: r1 = v5; break;
case 4: r1 = v4; break;
case 3: r1 = v3; break;
case 2: r1 = v2; break;
case 1: r1 = v1; break;
case 0: r1 = v0; break;
}


if (k1 == k0)
{
   ret = (r1 >> wb1);
}
else
{
   ret = (r1 >> wb1) | (r0 << (10-1 - wb0)) ; 
}

  return ret & RET_MASK;
}


unsigned int manhattanDistance32(unsigned int hp, unsigned int ht)
{
    unsigned int d = 0;
    
    #pragma unroll
    for (int k=0; k < 32; k++)
    {
        int bp = (hp >> k) & 0x1;
        int bt = (ht >> k) & 0x1;
        int dv = (bp & ~bt) ? 1:0; 
        d += dv;
    }
    
    return d;
}


unsigned int manhattanDistance(
        unsigned int hp0,
        unsigned int hp1,
        unsigned int hp2,
        unsigned int hp3,
        unsigned int hp4,
        unsigned int hp5,
        unsigned int hp6,
        unsigned int hp7,
        unsigned int hp8,
        unsigned int hp9,
        unsigned int hp10,
        unsigned int hp11,
        unsigned int hp12,
        unsigned int hp13,
        unsigned int hp14,
        unsigned int hp15,
        unsigned int hp16,
        unsigned int hp17,
        unsigned int hp18,
        unsigned int hp19,
        unsigned int hp20,
        unsigned int hp21,
        unsigned int hp22,
        unsigned int hp23,
        unsigned int hp24,
        unsigned int hp25,
        unsigned int hp26,
        unsigned int hp27,
        unsigned int hp28,
        unsigned int hp29,
        unsigned int hp30,
        unsigned int hp31,
        unsigned int ht0,
        unsigned int ht1,
        unsigned int ht2,
        unsigned int ht3,
        unsigned int ht4,
        unsigned int ht5,
        unsigned int ht6,
        unsigned int ht7,
        unsigned int ht8,
        unsigned int ht9,
        unsigned int ht10,
        unsigned int ht11,
        unsigned int ht12,
        unsigned int ht13,
        unsigned int ht14,
        unsigned int ht15,
        unsigned int ht16,
        unsigned int ht17,
        unsigned int ht18,
        unsigned int ht19,
        unsigned int ht20,
        unsigned int ht21,
        unsigned int ht22,
        unsigned int ht23,
        unsigned int ht24,
        unsigned int ht25,
        unsigned int ht26,
        unsigned int ht27,
        unsigned int ht28,
        unsigned int ht29,
        unsigned int ht30,
        unsigned int ht31
        )
{
    int d = 0;
    

    d += manhattanDistance32(hp0, ht0); 
    d += manhattanDistance32(hp1, ht1); 
    d += manhattanDistance32(hp2, ht2); 
    d += manhattanDistance32(hp3, ht3); 
    d += manhattanDistance32(hp4, ht4); 
    d += manhattanDistance32(hp5, ht5); 
    d += manhattanDistance32(hp6, ht6); 
    d += manhattanDistance32(hp7, ht7); 
    d += manhattanDistance32(hp8, ht8); 
    d += manhattanDistance32(hp9, ht9); 
    d += manhattanDistance32(hp10, ht10); 
    d += manhattanDistance32(hp11, ht11); 
    d += manhattanDistance32(hp12, ht12); 
    d += manhattanDistance32(hp13, ht13); 
    d += manhattanDistance32(hp14, ht14); 
    d += manhattanDistance32(hp15, ht15); 
    d += manhattanDistance32(hp16, ht16); 
    d += manhattanDistance32(hp17, ht17); 
    d += manhattanDistance32(hp18, ht18); 
    d += manhattanDistance32(hp19, ht19); 
    d += manhattanDistance32(hp20, ht20); 
    d += manhattanDistance32(hp21, ht21); 
    d += manhattanDistance32(hp22, ht22); 
    d += manhattanDistance32(hp23, ht23); 
    d += manhattanDistance32(hp24, ht24); 
    d += manhattanDistance32(hp25, ht25); 
    d += manhattanDistance32(hp26, ht26); 
    d += manhattanDistance32(hp27, ht27); 
    d += manhattanDistance32(hp28, ht28); 
    d += manhattanDistance32(hp29, ht29); 
    d += manhattanDistance32(hp30, ht30); 
    d += manhattanDistance32(hp31, ht31); 
    
    return (d+KMER_K-1)/KMER_K;
}



void computeHistogramPattern(__global unsigned char* restrict pattern , __global unsigned int* restrict patternIdx, unsigned int pi, 
        unsigned int* pH0,
        unsigned int* pH1,
        unsigned int* pH2,
        unsigned int* pH3,
        unsigned int* pH4,
        unsigned int* pH5,
        unsigned int* pH6,
        unsigned int* pH7,
        unsigned int* pH8,
        unsigned int* pH9,
        unsigned int* pH10,
        unsigned int* pH11,
        unsigned int* pH12,
        unsigned int* pH13,
        unsigned int* pH14,
        unsigned int* pH15,
        unsigned int* pH16,
        unsigned int* pH17,
        unsigned int* pH18,
        unsigned int* pH19,
        unsigned int* pH20,
        unsigned int* pH21,
        unsigned int* pH22,
        unsigned int* pH23,
        unsigned int* pH24,
        unsigned int* pH25,
        unsigned int* pH26,
        unsigned int* pH27,
        unsigned int* pH28,
        unsigned int* pH29,
        unsigned int* pH30,
        unsigned int* pH31
        )
{
    // in version 12 we assume a fixed value for the offset
    // 
    unsigned int offset = pi * 512 /  8;
    int len = 100;
    unsigned lenBytes = (100*2) / 8;
/*
    unsigned int offset = patternIdx[pi*INDEX_SIZE + 0];
    int len = patternIdx[pi*INDEX_SIZE+1];         // length in bases
    unsigned int lenBytes = alignedSequenceSize(patternIdx[pi*INDEX_SIZE+1]);
  */  
    
    *pH0 = 0;
    *pH1 = 0;
    *pH2 = 0;
    *pH3 = 0;
    *pH4 = 0;
    *pH5 = 0;
    *pH6 = 0;
    *pH7 = 0;
    *pH8 = 0;
    *pH9 = 0;
    *pH10 = 0;
    *pH11 = 0;
    *pH12 = 0;
    *pH13 = 0;
    *pH14 = 0;
    *pH15 = 0;
    *pH16 = 0;
    *pH17 = 0;
    *pH18 = 0;
    *pH19 = 0;
    *pH20 = 0;
    *pH21 = 0;
    *pH22 = 0;
    *pH23 = 0;
    *pH24 = 0;
    *pH25 = 0;
    *pH26 = 0;
    *pH27 = 0;
    *pH28 = 0;
    *pH29 = 0;
    *pH30 = 0;
    *pH31 = 0;


unsigned int v15;
unsigned int v14;
unsigned int v13;
unsigned int v12;
unsigned int v11;
unsigned int v10;
unsigned int v9;
unsigned int v8;
unsigned int v7;
unsigned int v6; 
unsigned int v5;
unsigned int v4;
unsigned int v3;
unsigned int v2;
unsigned int v1;
unsigned int v0;

    readBigEndian512bits(&pattern[offset], &v15,&v14,&v13,&v12,
&v11,&v10,&v9,&v8,
&v7,&v6,&v5,&v4,
&v3,&v2,&v1,&v0);

    int numBases = (512/2);


    #pragma unroll
    for (int i=0; i < (PATTERN_LEN-KMER_K+1); i++)
    {
        unsigned int kmer_index = get_kmer_index(v15,v14,v13,v12,v11,v10,v9,v8,v7,v6,v5,v4,v3,v2,v1,v0,i);

        int chunk = kmer_index / 32;
        int chunk_bit = kmer_index % 32;

        *pH0 |= (chunk == 0)? (1 << chunk_bit) : 0;
        *pH1 |= (chunk == 1)? (1 << chunk_bit) : 0;
        *pH2 |= (chunk == 2)? (1 << chunk_bit) : 0;
        *pH3 |= (chunk == 3)? (1 << chunk_bit) : 0;
        *pH4 |= (chunk == 4)? (1 << chunk_bit) : 0;
        *pH5 |= (chunk == 5)? (1 << chunk_bit) : 0;
        *pH6 |= (chunk == 6)? (1 << chunk_bit) : 0;
        *pH7 |= (chunk == 7)? (1 << chunk_bit) : 0;
        *pH8 |= (chunk == 8)? (1 << chunk_bit) : 0;
        *pH9 |= (chunk == 9)? (1 << chunk_bit) : 0;
        *pH10 |= (chunk == 10)? (1 << chunk_bit) : 0;
        *pH11 |= (chunk == 11)? (1 << chunk_bit) : 0;
        *pH12 |= (chunk == 12)? (1 << chunk_bit) : 0;
        *pH13 |= (chunk == 13)? (1 << chunk_bit) : 0;
        *pH14 |= (chunk == 14)? (1 << chunk_bit) : 0;
        *pH15 |= (chunk == 15)? (1 << chunk_bit) : 0;
        *pH16 |= (chunk == 16)? (1 << chunk_bit) : 0;
        *pH17 |= (chunk == 17)? (1 << chunk_bit) : 0;
        *pH18 |= (chunk == 18)? (1 << chunk_bit) : 0;
        *pH19 |= (chunk == 19)? (1 << chunk_bit) : 0;
        *pH20 |= (chunk == 20)? (1 << chunk_bit) : 0;
        *pH21 |= (chunk == 21)? (1 << chunk_bit) : 0;
        *pH22 |= (chunk == 22)? (1 << chunk_bit) : 0;
        *pH23 |= (chunk == 23)? (1 << chunk_bit) : 0;
        *pH24 |= (chunk == 24)? (1 << chunk_bit) : 0;
        *pH25 |= (chunk == 25)? (1 << chunk_bit) : 0;
        *pH26 |= (chunk == 26)? (1 << chunk_bit) : 0;
        *pH27 |= (chunk == 27)? (1 << chunk_bit) : 0;
        *pH28 |= (chunk == 28)? (1 << chunk_bit) : 0;
        *pH29 |= (chunk == 29)? (1 << chunk_bit) : 0;
        *pH30 |= (chunk == 30)? (1 << chunk_bit) : 0;
        *pH31 |= (chunk == 31)? (1 << chunk_bit) : 0;
    }
}

void computeHistogramText(__global unsigned char* restrict pattern , __global unsigned int* restrict patternIdx, unsigned int pi, 
        unsigned int* pH0,
        unsigned int* pH1,
        unsigned int* pH2,
        unsigned int* pH3,
        unsigned int* pH4,
        unsigned int* pH5,
        unsigned int* pH6,
        unsigned int* pH7,
        unsigned int* pH8,
        unsigned int* pH9,
        unsigned int* pH10,
        unsigned int* pH11,
        unsigned int* pH12,
        unsigned int* pH13,
        unsigned int* pH14,
        unsigned int* pH15,
        unsigned int* pH16,
        unsigned int* pH17,
        unsigned int* pH18,
        unsigned int* pH19,
        unsigned int* pH20,
        unsigned int* pH21,
        unsigned int* pH22,
        unsigned int* pH23,
        unsigned int* pH24,
        unsigned int* pH25,
        unsigned int* pH26,
        unsigned int* pH27,
        unsigned int* pH28,
        unsigned int* pH29,
        unsigned int* pH30,
        unsigned int* pH31
        )
{
    unsigned int offset = pi * 512 / 8;
    int len = 140;
    unsigned int lenBytes = (140 * 2) / 8;
/*
    unsigned int offset = patternIdx[pi*INDEX_SIZE + 0];
    int len = patternIdx[pi*INDEX_SIZE+1];         // length in bases
    unsigned int lenBytes = alignedSequenceSize(patternIdx[pi*INDEX_SIZE+1]);
 */   
    
    *pH0 = 0;
    *pH1 = 0;
    *pH2 = 0;
    *pH3 = 0;
    *pH4 = 0;
    *pH5 = 0;
    *pH6 = 0;
    *pH7 = 0;
    *pH8 = 0;
    *pH9 = 0;
    *pH10 = 0;
    *pH11 = 0;
    *pH12 = 0;
    *pH13 = 0;
    *pH14 = 0;
    *pH15 = 0;
    *pH16 = 0;
    *pH17 = 0;
    *pH18 = 0;
    *pH19 = 0;
    *pH20 = 0;
    *pH21 = 0;
    *pH22 = 0;
    *pH23 = 0;
    *pH24 = 0;
    *pH25 = 0;
    *pH26 = 0;
    *pH27 = 0;
    *pH28 = 0;
    *pH29 = 0;
    *pH30 = 0;
    *pH31 = 0;


    unsigned int v15;
unsigned int v14;
unsigned int v13;
unsigned int v12;
unsigned int v11;
unsigned int v10;
unsigned int v9;
unsigned int v8;
unsigned int v7;
unsigned int v6; 
unsigned int v5;
unsigned int v4;
unsigned int v3;
unsigned int v2;
unsigned int v1;
unsigned int v0;

 readBigEndian512bits(&pattern[offset], &v15,&v14,&v13,&v12,
&v11,&v10,&v9,&v8,
&v7,&v6,&v5,&v4,
&v3,&v2,&v1,&v0);

    #pragma unroll
    for (int i=0; i < (TEXT_LEN-KMER_K+1); i++)
    {
        unsigned int kmer_index = get_kmer_index(v15,v14,v13,v12,v11,v10,v9,v8,v7,v6,v5,v4,v3,v2,v1,v0,i);

        int chunk = kmer_index / 32;
        int chunk_bit = kmer_index % 32;

        *pH0 |= (chunk == 0)? (1 << chunk_bit) : 0;
        *pH1 |= (chunk == 1)? (1 << chunk_bit) : 0;
        *pH2 |= (chunk == 2)? (1 << chunk_bit) : 0;
        *pH3 |= (chunk == 3)? (1 << chunk_bit) : 0;
        *pH4 |= (chunk == 4)? (1 << chunk_bit) : 0;
        *pH5 |= (chunk == 5)? (1 << chunk_bit) : 0;
        *pH6 |= (chunk == 6)? (1 << chunk_bit) : 0;
        *pH7 |= (chunk == 7)? (1 << chunk_bit) : 0;
        *pH8 |= (chunk == 8)? (1 << chunk_bit) : 0;
        *pH9 |= (chunk == 9)? (1 << chunk_bit) : 0;
        *pH10 |= (chunk == 10)? (1 << chunk_bit) : 0;
        *pH11 |= (chunk == 11)? (1 << chunk_bit) : 0;
        *pH12 |= (chunk == 12)? (1 << chunk_bit) : 0;
        *pH13 |= (chunk == 13)? (1 << chunk_bit) : 0;
        *pH14 |= (chunk == 14)? (1 << chunk_bit) : 0;
        *pH15 |= (chunk == 15)? (1 << chunk_bit) : 0;
        *pH16 |= (chunk == 16)? (1 << chunk_bit) : 0;
        *pH17 |= (chunk == 17)? (1 << chunk_bit) : 0;
        *pH18 |= (chunk == 18)? (1 << chunk_bit) : 0;
        *pH19 |= (chunk == 19)? (1 << chunk_bit) : 0;
        *pH20 |= (chunk == 20)? (1 << chunk_bit) : 0;
        *pH21 |= (chunk == 21)? (1 << chunk_bit) : 0;
        *pH22 |= (chunk == 22)? (1 << chunk_bit) : 0;
        *pH23 |= (chunk == 23)? (1 << chunk_bit) : 0;
        *pH24 |= (chunk == 24)? (1 << chunk_bit) : 0;
        *pH25 |= (chunk == 25)? (1 << chunk_bit) : 0;
        *pH26 |= (chunk == 26)? (1 << chunk_bit) : 0;
        *pH27 |= (chunk == 27)? (1 << chunk_bit) : 0;
        *pH28 |= (chunk == 28)? (1 << chunk_bit) : 0;
        *pH29 |= (chunk == 29)? (1 << chunk_bit) : 0;
        *pH30 |= (chunk == 30)? (1 << chunk_bit) : 0;
        *pH31 |= (chunk == 31)? (1 << chunk_bit) : 0;

    }
}


/**
 * 
 * @param pattern
 * @param patternIdx
 * @param text
 * @param textIdx
 * @param workload
 * @param wi
 */
unsigned int computeTask(__global unsigned char* restrict pattern , __global unsigned int* restrict patternIdx, unsigned int pi ,
                         __global unsigned char* restrict text, __global unsigned int* restrict textIdx, unsigned int ti)
{
    
    unsigned int hp0;
    unsigned int hp1;
    unsigned int hp2;
    unsigned int hp3;
    unsigned int hp4;
    unsigned int hp5;
    unsigned int hp6;
    unsigned int hp7;
    unsigned int hp8;
    unsigned int hp9;
    unsigned int hp10;
    unsigned int hp11;
    unsigned int hp12;
    unsigned int hp13;
    unsigned int hp14;
    unsigned int hp15;
    unsigned int hp16;
    unsigned int hp17;
    unsigned int hp18;
    unsigned int hp19;
    unsigned int hp20;
    unsigned int hp21;
    unsigned int hp22;
    unsigned int hp23;
    unsigned int hp24;
    unsigned int hp25;
    unsigned int hp26;
    unsigned int hp27;
    unsigned int hp28;
    unsigned int hp29;
    unsigned int hp30;
    unsigned int hp31;
    unsigned int ht0;
    unsigned int ht1;
    unsigned int ht2;
    unsigned int ht3;
    unsigned int ht4;
    unsigned int ht5;
    unsigned int ht6;
    unsigned int ht7;
    unsigned int ht8;
    unsigned int ht9;
    unsigned int ht10;
    unsigned int ht11;
    unsigned int ht12;
    unsigned int ht13;
    unsigned int ht14;
    unsigned int ht15;
    unsigned int ht16;
    unsigned int ht17;
    unsigned int ht18;
    unsigned int ht19;
    unsigned int ht20;
    unsigned int ht21;
    unsigned int ht22;
    unsigned int ht23;
    unsigned int ht24;
    unsigned int ht25;
    unsigned int ht26;
    unsigned int ht27;
    unsigned int ht28;
    unsigned int ht29;
    unsigned int ht30;
    unsigned int ht31;
    
    computeHistogramPattern(pattern, patternIdx, pi, 
            &hp0, &hp1, &hp2, &hp3, &hp4, &hp5, &hp6, &hp7, &hp8, &hp9, 
            &hp10, &hp11, &hp12, &hp13, &hp14, &hp15, &hp16, &hp17, &hp18, &hp19, 
            &hp20, &hp21, &hp22, &hp23, &hp24, &hp25, &hp26, &hp27, &hp28, &hp29, 
            &hp30, &hp31);
    computeHistogramText(text, textIdx, ti, 
            &ht0, &ht1, &ht2, &ht3, &ht4, &ht5, &ht6, &ht7, &ht8, &ht9, 
            &ht10, &ht11, &ht12, &ht13, &ht14, &ht15, &ht16, &ht17, &ht18, &ht19, 
            &ht20, &ht21, &ht22, &ht23, &ht24, &ht25, &ht26, &ht27, &ht28, &ht29, 
            &ht30, &ht31);
    
    unsigned int d =  manhattanDistance(
            hp0, hp1, hp2, hp3, hp4, hp5, hp6, hp7, hp8, hp9,
            hp10, hp11, hp12, hp13, hp14, hp15, hp16, hp17, hp18, hp19,
            hp20, hp21, hp22, hp23, hp24, hp25, hp26, hp27, hp28, hp29,
            hp30, hp31,
            ht0, ht1, ht2, ht3, ht4, ht5, ht6, ht7, ht8, ht9,
            ht10, ht11, ht12, ht13, ht14, ht15, ht16, ht17, ht18, ht19,
            ht20, ht21, ht22, ht23, ht24, ht25, ht26, ht27, ht28, ht29,
            ht30, ht31);
    
    return d;
}
