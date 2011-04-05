
#include "_et_line_backproject_attenuated_gpu.h"

__device__ __constant__ int3 c_backprojection_size;


__global__ void et_line_backproject_attenuated_gpu_kernel(float *g_sinogram, float *g_backprojection, float *g_attenuation)
{
	__shared__ float s_sino[BLOCK];
	
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_backprojection_size.x * c_backprojection_size.z;
	if(tid<pixelNumber){
		//load sinogram to shared mem
		s_sino[threadIdx.x] = g_sinogram[tid];
		
		//backproject
                float sum_attenuation = 0.0f;
		unsigned int index = tid;
		for(unsigned int y=0; y < c_backprojection_size.y; y++){
                        sum_attenuation += g_attenuation[index];
			g_backprojection[index] = s_sino[threadIdx.x]*exp(-sum_attenuation);
                        //g_backprojection[index] = g_sinogram[tid];
			index += pixelNumber;
		}
	}
	return;  
}


