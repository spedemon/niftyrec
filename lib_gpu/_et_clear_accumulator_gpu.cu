
#include "_et_clear_accumulator_gpu.h"

#define BLOCK 256

void et_clear_accumulator_gpu(float **d_accumulator, nifti_image *accumulator)
{
	cudaMemset((void*)*d_accumulator,0.0f,accumulator->nx*accumulator->ny*accumulator->nz*sizeof(float));
}				

