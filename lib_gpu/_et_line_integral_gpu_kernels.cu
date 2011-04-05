#include "_et_line_integral_gpu.h"

__device__ __constant__ int3 c_ImageSize;

/*
__global__ void et_line_integral_gpu_kernel(float *g_activity, float *g_sinogram)
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
//	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.y;
	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.z;
	if(tid<pixelNumber){
		unsigned int index=tid;
		float sum=0.0f;
		for(unsigned int y=0; y<c_ImageSize.y; y++){
			sum += g_activity[index];
			//index += pixelNumber;
			index += c_ImageSize.x;
			
//			if(tid==0)
//			    g_sinogram[y] = g_activity[y*10];
//			if(tid==11)
//			    g_sinogram[y+10] = g_activity[y*10+11];			    
		}
		g_sinogram[tid]=sum;
		//g_sinogram[tid]=g_activity[index];
	}
	return;  
*/


__global__ void et_line_integral_gpu_kernel(float *g_activity, float *g_sinogram)
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.z;
	if(tid<pixelNumber){
		unsigned int index=tid;
		float sum=0.0f;
		for(unsigned int y=0; y<c_ImageSize.y; y++)
                        {
			sum += g_activity[index];
			index += pixelNumber;
			}
		g_sinogram[tid]=sum;
		//g_sinogram[tid]=g_activity[tid+128*128*50];
	}
	return; 	


/*
__global__ void et_line_integral_gpu_kernel(float *g_activity, float *g_sinogram)
{
	const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int pixelNumber = c_ImageSize.x * c_ImageSize.z;
	if(tid<pixelNumber){
		unsigned int index= (tid / c_ImageSize.z)*(c_ImageSize.z*c_ImageSize.y-c_ImageSize.z)+tid;
		float sum=0.0f;
                for(unsigned int y=0; y<c_ImageSize.y; y++)
                        {
                        sum   += g_activity[index];
                        index += c_ImageSize.z;
                        }
		g_sinogram[tid]=sum;
	}
	return; 
*/
	
/*	extern __shared__ float line[];
	unsigned int t      = threadIdx.x;
	unsigned int bx     = blockIdx.x;
	unsigned int by     = blockIdx.y;
	unsigned int bdim   = blockDim.x;
	unsigned int gdim_x = gridDim.x;
	unsigned int gdim_y = gridDim.y;
	unsigned int bdim_h = bdim/2;
	unsigned int index  = 0;

	// perform first level of reduction, while reading from global memory and writing to shared memory
	index = by*bdim + bx*gdim_y; 
	line[t] = g_activity[t + index];
	__syncthreads();
	index += bdim_h;
	line[t] = line[t] + g_activity[t + index];
	__syncthreads();

	for (index=bdim_h; index>32; index>>=1)
	{
	    if (t < index)
	        line[t] += line[t + index];
	    __syncthreads();
	}
	if (t < 32)
	{
	    line[t] += line[t + 32];
	    line[t] += line[t + 16];
	    line[t] += line[t + 8];
	    line[t] += line[t + 4];
	    line[t] += line[t + 2];
	    line[t] += line[t + 1];
	}

	// write result for this block to global mem
	index = by*bdim + bx*gdim_y; 
	if (t == 0) g_sinogram[index] = line[0]; */	
}


