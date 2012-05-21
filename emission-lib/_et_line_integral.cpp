
#include "_et_line_integral.h"
#ifdef _OPENMP
#include "omp.h"
#endif

void et_line_integral(nifti_image *inputImage, nifti_image *sinoImage, int cam)
{
    float *sino_data = (float *) (sinoImage->data) + cam*inputImage->nx*inputImage->ny ;
    float *input_data = (float *) (inputImage->data);
    memset((void*) sino_data,0,inputImage->nx*inputImage->ny*sizeof(float));
    for(int y=0; y<inputImage->ny; y++) {
//#pragma omp parallel for
        for(int x=0; x<inputImage->nx; x++) {
            for(int z=0; z<inputImage->nz; z++) {
                sino_data[x + y*inputImage->nx] = sino_data[x + y*inputImage->nx] + input_data[x + y*inputImage->nx + z*inputImage->nx*inputImage->ny];
            }
        }
    }

}


