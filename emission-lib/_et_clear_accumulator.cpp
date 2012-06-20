/*
 *  _et_clear_accumulator.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#include "_et_clear_accumulator.h"

void et_clear_accumulator(nifti_image *accumulator)
{
    float *data = (float*) accumulator->data;
    //for(int i=0; i<accumulator->nvox; i++)
    //    data[i] = 0.0f;
    memset(data, 0, accumulator->nvox*sizeof(float));
}
