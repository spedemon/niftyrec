
#include "_et_convolveSeparable2D.h"


extern "C" void convolutionRow(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR)
{
    for(int y = 0; y < imageH; y++)
        for(int x = 0; x < imageW; x++){
            float sum = 0;
            for(int k = -kernelR; k <= kernelR; k++){
                int d = x + k;
                if(d >= 0 && d < imageW)
                    sum += h_Src[y * imageW + d] * h_Kernel[kernelR - k];
            }
            h_Dst[y * imageW + x] = sum;
        }
}


extern "C" void convolutionColumn(float *h_Dst,float *h_Src,float *h_Kernel,int imageW,int imageH,int kernelR)
{
    for(int y = 0; y < imageH; y++)
        for(int x = 0; x < imageW; x++){
            float sum = 0;
            for(int k = -kernelR; k <= kernelR; k++){
                int d = y + k;
                if(d >= 0 && d < imageH)
                    sum += h_Src[d * imageW + x] * h_Kernel[kernelR - k];
            }
            h_Dst[y * imageW + x] = sum;
        }
}


int et_convolveSeparable2D(nifti_image* inputImage, float *kernelSeparated, int kernelLengthH, int kernelLengthW, nifti_image *outputImage, float background)
{
    int imageH = inputImage->ny;
    int imageW = inputImage->nx;
    int n_slices = inputImage->nz;
    int kernelRadiusH = (kernelLengthH-1)/2;
    int kernelRadiusW = (kernelLengthW-1)/2;

    int image_slice_size = imageH * imageW;
    int kernel_slice_size = kernelLengthW + kernelLengthH;

    float *buffer = (float*) malloc(imageH*imageW*sizeof(float));

    float *input, *output, *kernel_row, *kernel_column;

    //Convolve slices one by one
    for (int slice=0; slice<n_slices; slice++)
        {
        // determine slice pointers
        input = (float*)inputImage->data + slice * image_slice_size;
        output = (float*)outputImage->data + slice * image_slice_size; 
        kernel_row = kernelSeparated + slice * kernel_slice_size;
        kernel_column = kernelSeparated + slice * kernel_slice_size + kernelLengthW;

        // convolve slice
        convolutionRow(buffer,input,kernel_row,imageW,imageH,kernelRadiusW);
        convolutionColumn(output,buffer,kernel_column,imageW,imageH,kernelRadiusH);
        }

    free(buffer);

    return 0;
}


