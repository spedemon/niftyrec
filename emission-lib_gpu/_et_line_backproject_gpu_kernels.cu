
#include "_et_line_backproject_gpu.h"

__device__ __constant__ int3 c_backprojection_size;


__global__ void et_line_backproject_gpu_kernel(float *g_sinogram, float *g_backprojection)
{
	__shared__ float s_sino[BLOCK];
	
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_backprojection_size.x * c_backprojection_size.z;
	if(tid<pixelNumber){
		//load sinogram to shared mem
		s_sino[threadIdx.x] = g_sinogram[tid];
		
		//backproject
		unsigned int index = tid;
		for(unsigned int y=0; y < c_backprojection_size.y; y++){
			g_backprojection[index] = s_sino[threadIdx.x];
                        //g_backprojection[index] = g_sinogram[tid];
			index += pixelNumber;
		}
	}
	return;  
}


