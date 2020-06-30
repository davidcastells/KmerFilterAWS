//#include "ihc_apint.h"

#define WORKLOAD_TASK_SIZE  3
#define INDEX_SIZE          2
#define BASE_SIZE           2
#define LOAD_BASES_ALIGNMENT_BITS   512

#define KMER_K          5
#define KMER_K_BITS     (KMER_K*BASE_SIZE)
#define KMER_BINS       (1 << KMER_K_BITS)

#define PATTERN_LEN      100
#define TEXT_LEN         140 

unsigned int alignedSequenceSize(int bases);

#include "histogram_fixed.v16.cl"

/**
 * Compute the number of bytes required to store the number of bases, considering
 * that we require memory alignment
 * @param bases
 * @return 
 */
unsigned int alignedSequenceSize(int bases)
{
    unsigned int lenBits = bases * BASE_SIZE;         // length in bases
    unsigned int lenAlignedUnits = ((lenBits  + (LOAD_BASES_ALIGNMENT_BITS-1)) / LOAD_BASES_ALIGNMENT_BITS) * LOAD_BASES_ALIGNMENT_BITS; // round up

    unsigned int lenAlignedUnitsBytes = lenAlignedUnits / 8;
    return lenAlignedUnitsBytes;
}

/**
 * from version 11 we assume that pattern index, and text index is the same worload index  
 */
void doWorkloadTask(__global unsigned char* restrict pattern , __global unsigned int* restrict patternIdx, 
                    __global unsigned char* restrict text, __global unsigned int* restrict textIdx,
                    /*__global*/ unsigned int* workload, unsigned int wi, unsigned int li)
{
    unsigned int pi = wi; // workload[wi*WORKLOAD_TASK_SIZE+0];
    unsigned int ti = wi; // workload[wi*WORKLOAD_TASK_SIZE+1];

    
    unsigned int d = computeTask(pattern, patternIdx, pi, text, textIdx, ti);
    
#ifdef FPGA_DEBUG
    printf("[FPGA] pi=%d  ti=%d ", pi, ti);
    printf(" task %d = %d\n", wi, d);
#endif
    
    //workload[wi*WORKLOAD_TASK_SIZE+2] = d;
    workload[li] = d;
}


__kernel void kmer(__global unsigned char* restrict pattern ,
		   __global unsigned int* restrict patternIdx, 
                   __global unsigned char* restrict text,
		   __global unsigned int* restrict textIdx,
                   __global unsigned int* workload, 
		   unsigned int workloadLength)
{
	#define WORKLOAD_CHUNK 1024*16
	unsigned int workload_result[WORKLOAD_CHUNK];
	
    for (int i=0; i < workloadLength; /*i++*/)
    {
    	int base_i = i;
        int li;
	
    	// compute to local memory
    	for (li=0; (li < WORKLOAD_CHUNK) && (i < workloadLength); li++, i++)
	{
	   doWorkloadTask(pattern, patternIdx, text, textIdx, workload_result, i, li);
	}
		
	// transfer the results back to the main table
	for (int ti=0; ti < li; ti++)
		workload[(base_i + ti)*WORKLOAD_TASK_SIZE+2] = workload_result[ti];
     }
}
