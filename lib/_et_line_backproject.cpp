
#include "_et_line_backproject.h"

void et_line_backproject(nifti_image *sinogramImage, nifti_image *bkprImage, int cam)
{
    float *sino_data  = (float *) (sinogramImage->data) + cam*bkprImage->nx*bkprImage->nz ;
    float *bkpr_data = (float *)  (bkprImage->data);

    for(int z=0; z<bkprImage->nz; z++) {
        for(int y=0; y<bkprImage->ny; y++) {
            for(int x=0; x<bkprImage->nx; x++) {
                bkpr_data[x + y*bkprImage->nx + z*bkprImage->nx*bkprImage->ny] = sino_data[x + z*bkprImage->ny];
            }
        }
    }

}

