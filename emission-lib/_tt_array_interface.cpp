
#include "_tt_array_interface.h"
#include "_et_common.h"

extern "C" int tt_array_project(float *attenuation, int *attenuation_size, float *projection, int *projection_size, float *image_origin, float *detector_origin, float *detector_shape, float *psf, int *psf_size, float background, int GPU)
{
	int status = 1;
        int n_projections;
        int no_psf = 0;

        n_projections = projection_size[2];
        fprintf_verbose("tt_array_project: number of projections: %d\n", n_projections);

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        /* Check consistency of input */
        //Attenuation must be of size [NxNxm]; m>=2
        if (attenuation_size[0] != attenuation_size[1] || attenuation_size[2]<2)
                {
                fprintf_verbose("tt_array_project: 3D attenuation must be of size [N,N,m]; m>=2.\n");
                return status;
                }
        //Size of projection must be consistent with attenuation size
        if (projection_size[0] != attenuation_size[0] || projection_size[1] != attenuation_size[2] || projection_size[2] != n_projections) 
                {
                fprintf_verbose("tt_array_project: 3D projection must be of size [N,m,n_projections] for attenuation of size [N,N,m] and 'n_projections' projections.\n");
                return status;
                }
        //Size of psf must be odd and consistent with attenuation size
        if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=attenuation_size[2])
                    {
                    fprintf_verbose("tt_array_project: 3D psf must be of size [h,k,m] for attenuation of size [N,N,m]; h,k odd.\n");
                    return status;
                    }
                }

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = attenuation_size[0];
	dim[2]    = attenuation_size[1];
	dim[3]    = attenuation_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1; 
	nifti_image *attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        attenuationImage->data = (float *)(attenuation);    
	
	// Allocate the result nifti image
	dim[1] = projection_size[0];
	dim[2] = projection_size[1];
	dim[3] = projection_size[2];	   

        nifti_image *projectionImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        projectionImage->data = (float *)(projection);

	//Allocate Point Spread Function nifti image
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float *)(psf);
            }

        //Do projection
        #ifdef _USE_CUDA
        if (GPU)
            status = tt_project_gpu(attenuationImage, projectionImage, psfImage, n_projections, image_origin, detector_origin, detector_shape, background);
        else
            status = tt_project(attenuationImage, projectionImage, psfImage, n_projections, image_origin, detector_origin, detector_shape, background);
        #else
            if (GPU)
                fprintf_verbose( "tt_array_project: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = tt_project(attenuationImage, projectionImage, psfImage, n_projections, image_origin, detector_origin, detector_shape, background);
        #endif

	//Free
	if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
	if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
	(void)nifti_free_extensions( attenuationImage ) ;
	free(attenuationImage) ;
	
	if( projectionImage->fname != NULL ) free(projectionImage->fname) ;
	if( projectionImage->iname != NULL ) free(projectionImage->iname) ;
	(void)nifti_free_extensions( projectionImage ) ;
	free(projectionImage) ;

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

	return status;
}

