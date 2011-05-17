
#include "_et_accumulate.h"

void et_accumulate(nifti_image *inputImage, nifti_image *accumulatorImage)
{
    float *input_data = (float*) inputImage->data;
    float *accumulator_data = (float*) accumulatorImage->data;
    for(int i=0; i<accumulatorImage->nvox; i++)
        {
        accumulator_data[i] += input_data[i];
        }
}
