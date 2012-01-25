
#include "_et_convolveSeparable2D_gpu_kernels.cu"
#include "_et_common.h"


////////////////////////////////////////////////////////////////////////////////
// 2D Convolution
////////////////////////////////////////////////////////////////////////////////

int et_convolveSeparable2D_gpu(float **d_data, int *data_size, float **d_kernel_separated, int *kernel_size, float **d_result)
{
    int status = 1;
    const int dataH = data_size[0];
    const int dataW = data_size[1];
    const int kernelRadius = (kernel_size[0]-1)/2;

    const int n_slices = data_size[2];
    const int data_slice_size = dataH * dataW;

    float *d_Data, *d_Kernel_separated, *d_Result, *d_Buffer;

    cutilSafeCall( cudaMalloc((void **)&d_Buffer , dataW * dataH * sizeof(float)) );

    //Convolve slices one by one
    for (int slice=0; slice<n_slices; slice++)
        {
        //Determine slice pointer
        d_Data = (*d_data) + slice * data_slice_size; 
        d_Result = (*d_result) + slice * data_slice_size;
        d_Kernel_separated = (*d_kernel_separated) + slice * 2 * kernel_size[0];

        //Convolve  //FIXME: make it up to 10% faster by using texture memory for the kernels 
        setConvolutionKernel(d_Kernel_separated,kernelRadius);
        convolutionRowsGPU(d_Buffer,d_Data,dataW,dataH,kernelRadius);
        convolutionColumnsGPU(d_Result,d_Buffer,dataW,dataH,kernelRadius);

//        cutilSafeCall( cutilDeviceSynchronize() );
        }
    cutilSafeCall( cudaFree(d_Buffer) );

    status = 0;
    return status;
}



