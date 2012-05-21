#include "_et_line_integral_gpu.h"
#include "_et_line_integral_gpu_kernels.cu"

#define BLOCK 256

void et_line_integral_gpu(float **d_activity, float **d_sinogram, int cam, nifti_image *img)
{
	int3 imageSize = make_int3(img->dim[1],img->dim[2],img->dim[3]);
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(c_ImageSize,&imageSize,sizeof(int3)));
	
	const unsigned int Grid = (unsigned int)ceil(img->dim[1]*img->dim[2]/(float)BLOCK);
	dim3 B1(BLOCK,1,1);
	dim3 G1(Grid,1,1);
	
	float *currentCamPointer = (*d_sinogram) + cam * img->dim[1] * img->dim[2] ;
	
	et_line_integral_gpu_kernel <<<G1,B1>>> (*d_activity, currentCamPointer);

	CUDA_SAFE_CALL(cudaThreadSynchronize());
}


