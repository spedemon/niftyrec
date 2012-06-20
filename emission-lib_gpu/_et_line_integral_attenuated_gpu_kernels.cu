/*
 *  _et_line_integrals_attenuated_gpu_kernels.cu
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_line_integral_attenuated_gpu.h"

__device__ __constant__ int3 c_ImageSize;


__global__ void et_line_integral_attenuated_gpu_kernel(float *g_activity, float *g_attenuation, float *g_sinogram)
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.y;
	if(tid<pixelNumber){
		unsigned int index=tid;
                float sum_attenuation=0.0f;
		float sum_activity=0.0f;
//                float background_activity=1.0f;
		for(unsigned int z=0; z<c_ImageSize.z; z++)
                        {
                        sum_attenuation += g_attenuation[index];
			sum_activity    += g_activity[index]*exp(-sum_attenuation);
			index += pixelNumber;
			}
//                sum_activity += background_activity*exp(-sum_attenuation);
		g_sinogram[tid]=sum_activity;
	}
	return; 	
}





