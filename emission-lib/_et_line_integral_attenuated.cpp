
#include "_et_line_integral_attenuated.h"
#ifdef _OPENMP
#include "omp.h"
#endif

void et_line_integral_attenuated(nifti_image *inputImage, nifti_image *attenuationImage, nifti_image *sinoImage, int cam) //FIXME test attenuation
{
    float *sino_data = (float *) (sinoImage->data) + cam*inputImage->nx*inputImage->ny ;
    float *input_data = (float *) (inputImage->data);
    float *attenuation_data = (float *) (attenuationImage->data);
    float sum_attenuation;
    int index;
    memset((void*) sino_data,0,inputImage->nx*inputImage->ny*sizeof(float));
    for(int y=0; y<inputImage->ny; y++) {
//#pragma omp parallel for
        for(int x=0; x<inputImage->nx; x++) {
            sum_attenuation = 0;
            for(int z=0; z<inputImage->nz; z++) {
                index = x + y*inputImage->nx + z*inputImage->nx*inputImage->ny; 
                sum_attenuation += attenuation_data[index]; 
                sino_data[x + y*inputImage->nx] = sino_data[x + y*inputImage->nx] + input_data[index]*exp(-sum_attenuation); 
            }
        }
    }

}


