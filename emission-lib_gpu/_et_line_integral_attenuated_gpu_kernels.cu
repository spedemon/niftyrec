#include "_et_line_integral_attenuated_gpu.h"

__device__ __constant__ int3 c_ImageSize;


__global__ void et_line_integral_attenuated_gpu_kernel(float *g_activity, float *g_attenuation, float *g_sinogram)
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.z;
	if(tid<pixelNumber){
		unsigned int index=tid;
                float sum_attenuation=0.0f;
		float sum_activity=0.0f;
		for(unsigned int y=0; y<c_ImageSize.y; y++)
                        {
                        sum_attenuation += g_attenuation[index];
			sum_activity    += g_activity[index]*exp(-sum_attenuation);
			index += pixelNumber;
			}
		g_sinogram[tid]=sum_activity;
		//g_sinogram[tid]=g_activity[tid+128*128*50];
	}
	return; 	
}



