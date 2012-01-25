
#include "_et_array_interface.h"
#include "_et_common.h"

extern "C" int et_array_affine(float *image_ptr, int *image_size, float *transformed_image_ptr, int *transformed_image_size, float *affine_ptr, int *affine_size, float background, int GPU)
{
	int status;

	// Allocate source image 
        int dim[8];
	dim[0]    = 3;
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *sourceImage;
        sourceImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sourceImage->data = static_cast<void *>(image_ptr);
	
	// Allocate the result image
	nifti_image *transformedImage = nifti_copy_nim_info(sourceImage);
        transformedImage->data = static_cast<void *>(transformed_image_ptr);

        // Allocate and initialise Affine transformation matrix
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        //FIXME

	// Rotate
	#ifdef _USE_CUDA
	if(GPU)
	    status = et_affine_gpu(sourceImage, transformedImage, affineTransformation, background);
	else
	    status = et_affine(sourceImage, transformedImage, affineTransformation, background);
      	#else
      	    status = et_affine(sourceImage, transformedImage, affineTransformation, background);
      	#endif
      	
	//Free
	if( sourceImage->fname != NULL ) free(sourceImage->fname) ;
	if( sourceImage->iname != NULL ) free(sourceImage->iname) ;
	(void)nifti_free_extensions( sourceImage ) ;
	free(sourceImage) ;
	
	if( transformedImage->fname != NULL ) free(transformedImage->fname) ;
	if( transformedImage->iname != NULL ) free(transformedImage->iname) ;
	(void)nifti_free_extensions( transformedImage ) ;
	free(transformedImage) ;

        free(affineTransformation);

	return status;
}



extern "C" int et_array_rotate(float *image, int *size, float *rotated_image, float *angles, float *centers, float background, int GPU)
{
	int status;

	// Allocate source image 
        int dim[8];
	dim[0]    = 3;
	dim[1]    = size[0];
	dim[2]    = size[1];
	dim[3]    = size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *sourceImage;
        sourceImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sourceImage->data = static_cast<void *>(image);
	
        // Allocate the result image
	nifti_image *resultImage = nifti_copy_nim_info(sourceImage);
        resultImage->data = static_cast<void *>(rotated_image);

        // Rotate
        #ifdef _USE_CUDA
        if(GPU)
            status = et_rotate_gpu(sourceImage, resultImage, angles[0],angles[1],angles[2], centers[0], centers[1], centers[2], background);
        else
            status = et_rotate(sourceImage, resultImage, angles[0],angles[1],angles[2], centers[0], centers[1], centers[2], background);
      	#else
      	    status = et_rotate(sourceImage, resultImage, angles[0],angles[1],angles[2], centers[0], centers[1], centers[2], background);
        #endif
      	
	//Free
	if( sourceImage->fname != NULL ) free(sourceImage->fname) ;
	if( sourceImage->iname != NULL ) free(sourceImage->iname) ;
	(void)nifti_free_extensions( sourceImage ) ;
	free(sourceImage) ;
	
	if( resultImage->fname != NULL ) free(resultImage->fname) ;
	if( resultImage->iname != NULL ) free(resultImage->iname) ;
	(void)nifti_free_extensions( resultImage ) ;
	free(resultImage) ;

	return status;
}



extern "C" int et_array_project(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU)
{
	int status = 1;
	int dims;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        float *cameras_array;

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

	// 2D or 3D?
        dims = 3;
	if (activity_size[2] == 1)
            dims = 2;

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3))
            {
            fprintf_verbose("et_array_project: Incorrect size of cameras %d %d. 'Cameras' must be either [n_cameras x 3] or [n_cameras x 1].\n",cameras_size[0],cameras_size[1]);
            return status;
            }
        if (dims==2)
            //Activity must be of size [NxN]
            {
            if (activity_size[0] != activity_size[1])
                {
                fprintf_verbose("et_array_project: 2D activity must be of size [N,N].\n");
                return status;
                }
            //Size of sinogram must be consistent with activity size
            if (sinogram_size[0] != activity_size[0] || sinogram_size[1] != n_cameras) 
                {
                fprintf_verbose("et_array_project: 2D sinogram must be of size [N,n_cameras] for activity of size [N,N] and 'n_cameras' cameras.\n");
                return status;
                }
            //Size of psf must be odd
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1)
                    {
                    fprintf_verbose("et_array_project: 2D psf must be of size [h,k]; h,k odd.\n");
                    return status;
                    }
                }
            }
        if (dims==3)
            //Activity must be of size [NxNxm]; m>=2
            {
            if (activity_size[0] != activity_size[1] || activity_size[2]<2)
                {
                fprintf_verbose("et_array_project: 3D activity must be of size [N,N,m]; m>=2.\n");
                return status;
                }
            //Size of sinogram must be consistent with activity size
            if (sinogram_size[0] != activity_size[0] || sinogram_size[1] != activity_size[2] || sinogram_size[2] != n_cameras) 
                {
                fprintf_verbose("et_array_project: 3D sinogram must be of size [N,m,n_cameras] for activity of size [N,N,m] and 'n_cameras' cameras.\n");
                return status;
                }
            //Size of psf must be odd and consistent with activity size
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=activity_size[2])
                    {
                    fprintf_verbose("et_array_project: 3D psf must be of size [h,k,m] for activity of size [N,N,m]; h,k odd.\n");
                    return status;
                    }
                }
            }

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (n_cameras_axis == 3)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[0*n_cameras+cam] = cameras[cam];
            }

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = activity_size[0];
	dim[2]    = activity_size[1];
	dim[3]    = activity_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	//fprintf_verbose("\nS: %d %d %d",activity_size[0],activity_size[1],activity_size[2]);
	nifti_image *activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        activityImage->data = (float *)(activity);

        // Allocate attenuation nifti image
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);        
            }

	// Allocate the result nifti image
	//fprintf_verbose( "\nN CAMERAS: %d ",n_cameras);
	if (dims == 2)
	   {
	   dim[1] = activity_size[0];
	   dim[2] = n_cameras;
	   dim[3] = 1;
	   }
	else
	   {
	   dim[1] = activity_size[0];
	   dim[2] = activity_size[2];
	   dim[3] = n_cameras;	   
	   }
        nifti_image *sinogramImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sinogramImage->data = (float *)(sinogram);

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
            status = et_project_gpu(activityImage, sinogramImage, psfImage, attenuationImage, cameras_array, n_cameras, background, background_attenuation);
        else
            status = et_project(activityImage, sinogramImage, psfImage, attenuationImage, cameras_array, n_cameras, background, background_attenuation);
        #else
            if (GPU)
                fprintf_verbose( "et_array_project: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_project(activityImage, sinogramImage, psfImage, attenuationImage, cameras_array, n_cameras, background, background_attenuation);
        #endif

	//Free
	if( activityImage->fname != NULL ) free(activityImage->fname) ;
	if( activityImage->iname != NULL ) free(activityImage->iname) ;
	(void)nifti_free_extensions( activityImage ) ;
	free(activityImage) ;

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }
	
	if( sinogramImage->fname != NULL ) free(sinogramImage->fname) ;
	if( sinogramImage->iname != NULL ) free(sinogramImage->iname) ;
	(void)nifti_free_extensions( sinogramImage ) ;
	free(sinogramImage) ;

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);

	return status;
}





extern "C" int et_array_backproject(float *sino, int *sino_size, float *bkpr, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU)
{
	int status;
	int dims;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        float *cameras_array;

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

	// 2D or 3D?
        dims = 3;
	if (sino_size[2] == 1)
            dims = 2;

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3))
            {
            fprintf_verbose("et_array_backproject: 'Cameras' must be either [n_cameras x 3] or [n_cameras x 1]\n");
            return status;
            }
        if (dims==2)
            //Sino must be of size [Nxn_cameras]
            {
            if (sino_size[1] != n_cameras)
                {
                fprintf_verbose("et_array_backproject: 2D sinogram must be of size [N,n_cameras].\n");
                return status;
                }
            //Size of psf must be odd
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]!=sino_size[0])
                    {
                    fprintf_verbose("et_array_backproject: 2D psf must be of size [h,N]; h odd.\n");
                    return status;
                    }
                }
            }
        if (dims==3)
            //Sino must be of size [Nxmxn_cameras]; m>=2
            {
            if (sino_size[2] != n_cameras)
                {
                fprintf_verbose("et_array_backproject: 3D sino must be of size [N,m,n_cameras].\n");
                return status;
                }
            //Size of psf must be odd and consistent with activity size
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=sino_size[0])
                    {
                    fprintf_verbose("et_array_backproject: 3D psf must be of size [h,k,N] for activity of size [N,N,m]; h,k odd.\n");
                    return status;
                    }
                }
            }

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (n_cameras_axis == 3)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[0*n_cameras+cam] = cameras[cam];
            }


	// Allocate backprojection (result) image 
        int dim[8];
	dim[0]    = 3;
	dim[1]    = bkpr_size[0];
	dim[2]    = bkpr_size[1];
	dim[3]    = bkpr_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *bkprImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        bkprImage->data = (float*) bkpr;
	
        // Allocate attenuation image 
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float*) attenuation;
            }

	// Allocate the sinogram (input) image
	dim[3]    = n_cameras;
        nifti_image *sinoImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        sinoImage->data = (float*) sino;

	//Allocate Point Spread Function
        nifti_image *psfImage = NULL;
        if (!no_psf)
            {
            dim[1] = psf_size[0];
            dim[2] = psf_size[1];
            dim[3] = psf_size[2];
            psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            psfImage->data = (float*) psf;
            }
	//Backproject
	#ifdef _USE_CUDA
	if(GPU)
	    status = et_backproject_gpu(sinoImage, bkprImage, psfImage, attenuationImage, cameras_array, n_cameras, background, background_attenuation);
	else
	    status = et_backproject(sinoImage, bkprImage, psfImage, attenuationImage, cameras_array, n_cameras, background, background_attenuation);	    
	#else
	    status = et_backproject(sinoImage, bkprImage, psfImage, attenuationImage, cameras_array, n_cameras, background, background_attenuation);	    
	#endif
	
	//Free (free nifti images but not their data arrays)
	if( sinoImage->fname != NULL ) free(sinoImage->fname) ;
	if( sinoImage->iname != NULL ) free(sinoImage->iname) ;
	(void)nifti_free_extensions( sinoImage ) ;
	free(sinoImage) ;
	
	if( bkprImage->fname != NULL ) free(bkprImage->fname) ;
	if( bkprImage->iname != NULL ) free(bkprImage->iname) ;
	(void)nifti_free_extensions( bkprImage ) ;
	free(bkprImage) ;

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);

	return status;
}


extern "C" int et_mlem_spect(float *sinogram_data, int size_x, int size_y, int n_cameras, float firstcamera, float lastcamera, int iterations, int use_psf, int use_ddpsf, int psf_size_x, int psf_size_y, float *psf_data, int use_attenuation, float *attenuation_data, float *activity_data, int GPU)
{
    return 0;
}

extern "C" int et_array_fisher_grid(float *activity_ptr, int *activity_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation, int *attenuation_size, float epsilon, float background, float background_attenuation, int GPU)
{
        int from_projection = 0;
	int status = 1;
	int dims;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        float *cameras_array;

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

	// 2D or 3D?
        dims = 3;
	if (activity_size[2] == 1)
            dims = 2;

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3))
            {
            fprintf_verbose("et_array_project: Incorrect size of cameras %d %d. 'Cameras' must be either [n_cameras x 3] or [n_cameras x 1].\n",cameras_size[0],cameras_size[1]);
            return status;
            }
        if (dims==2)
            //Activity must be of size [NxN]
            {
            if (activity_size[0] != activity_size[1])
                {
                fprintf_verbose("et_array_project: 2D activity must be of size [N,N].\n");
                return status;
                }
            //Size of psf must be odd
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1)
                    {
                    fprintf_verbose("et_array_project: 2D psf must be of size [h,k]; h,k odd.\n");
                    return status;
                    }
                }
            }
        if (dims==3)
            //Activity must be of size [NxNxm]; m>=2
            {
            if (activity_size[0] != activity_size[1] || activity_size[2]<2)
                {
                fprintf_verbose("et_array_project: 3D activity must be of size [N,N,m]; m>=2.\n");
                return status;
                }
            //Size of psf must be odd and consistent with activity size
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=activity_size[2])
                    {
                    fprintf_verbose("et_array_project: 3D psf must be of size [h,k,m] for activity of size [N,N,m]; h,k odd.\n");
                    return status;
                    }
                }
            }

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (n_cameras_axis == 3)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[0*n_cameras+cam] = cameras[cam];
            }

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = activity_size[0];
	dim[2]    = activity_size[1];
	dim[3]    = activity_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	//fprintf_verbose("\nS: %d %d %d",activity_size[0],activity_size[1],activity_size[2]);
	nifti_image *activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        activityImage->data = (float *)(activity_ptr);

        // Allocate grid nifti image
        nifti_image *gridImage = NULL;
        gridImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        gridImage->data = (float *)(grid_ptr);        

        // Allocate attenuation nifti image
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);        
            }

	// Allocate the result Fisher nifti image
        dim[1] = fisher_size[0];
        dim[2] = fisher_size[1];
        dim[3] = 1;
        nifti_image *fisherImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        fisherImage->data = (float *)(fisher_ptr);

        nifti_image *fisherpriorImage;
        if (fisher_prior_ptr!=NULL)
            {
            fisherpriorImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            fisherpriorImage->data = (float *)(fisher_prior_ptr);
            }

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

        //Compute Fisher Information Matrix
        #ifdef _USE_CUDA
        if (GPU)
            status = et_fisher_grid_gpu(from_projection, activityImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, n_cameras, epsilon, background, background_attenuation); 
        else
            status = et_fisher_grid(from_projection, activityImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, n_cameras, epsilon, background, background_attenuation); 
        #else
            if (GPU)
                fprintf_verbose( "et_array_project: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_fisher_grid(from_projection, activityImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, n_cameras, epsilon, background, background_attenuation); 
        #endif

	//Free
	if( activityImage->fname != NULL ) free(activityImage->fname) ;
	if( activityImage->iname != NULL ) free(activityImage->iname) ;
	(void)nifti_free_extensions( activityImage ) ;
	free(activityImage) ;

	(void)nifti_free_extensions( gridImage ) ;
	free(gridImage) ;

	(void)nifti_free_extensions( fisherImage ) ;
	free(fisherImage) ;

        if (fisher_prior_ptr!=NULL)
            {         
            (void)nifti_free_extensions( fisherpriorImage ) ;
            free(fisherpriorImage) ;
            }

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);

	return status;
}


extern "C" int et_array_fisher_grid_projection(float *sinogram_ptr, int *sinogram_size, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation, int *attenuation_size, float epsilon, float background, float background_attenuation, int GPU)
{
        int from_projection = 1;
	int status = 1;
	int dims;
        int n_cameras;
        int n_cameras_axis;
        int no_psf = 0;
        int no_attenuation = 0;
        float *cameras_array;

        n_cameras = cameras_size[0];
        n_cameras_axis = cameras_size[1];

	// 2D or 3D?
        dims = 3;
	if (bkpr_size[2] == 1)
            dims = 2;

        //PSF or not?
        if (psf_size[0] == 0 && psf_size[1] == 0 && psf_size[2] == 0)
            no_psf = 1;

        //attenuation or not?
        if (attenuation_size[0] == 0 && attenuation_size[1] == 0 && attenuation_size[2] == 0)
            no_attenuation = 1;

        /* Check consistency of input */
        // Cameras must specify all 3 axis of rotation (3D array) or can be a 1D array if rotation is only along z axis.
        if (!(n_cameras_axis == 1 || n_cameras_axis == 3))
            {
            fprintf_verbose("et_array_project: Incorrect size of cameras %d %d. 'Cameras' must be either [n_cameras x 3] or [n_cameras x 1].\n",cameras_size[0],cameras_size[1]);
            return status;
            }
        if (dims==2)
            //Activity must be of size [NxN]
            {
            if (bkpr_size[0] != bkpr_size[1])
                {
                fprintf_verbose("et_array_project: 2D activity must be of size [N,N].\n");
                return status;
                }
            //Size of psf must be odd
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1)
                    {
                    fprintf_verbose("et_array_project: 2D psf must be of size [h,k]; h,k odd.\n");
                    return status;
                    }
                }
            }
        if (dims==3)
            //Activity must be of size [NxNxm]; m>=2
            {
            if (bkpr_size[0] != bkpr_size[1] || bkpr_size[2]<2)
                {
                fprintf_verbose("et_array_project: 3D activity must be of size [N,N,m]; m>=2.\n");
                return status;
                }
            //Size of psf must be odd and consistent with activity size
            if (!no_psf)
                {
                if (psf_size[0]%2!=1 || psf_size[1]%2!=1 || psf_size[2]!=bkpr_size[2])
                    {
                    fprintf_verbose("et_array_project: 3D psf must be of size [h,k,m] for activity of size [N,N,m]; h,k odd.\n");
                    return status;
                    }
                }
            }

        // Allocate array for cameras
        cameras_array = (float *)malloc(n_cameras*3*sizeof(float));
        if (n_cameras_axis == 3)
            memcpy((void*) cameras_array, (void*) cameras, n_cameras*3*sizeof(float));
        if (n_cameras_axis == 1)
            {
            memset(cameras_array, 0, n_cameras*3*sizeof(float));
            for (int cam=0; cam<n_cameras; cam++)
                cameras_array[0*n_cameras+cam] = cameras[cam];
            }

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = sinogram_size[0];
	dim[2]    = sinogram_size[1];
	dim[3]    = sinogram_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *projectionImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        projectionImage->data = (float *)(sinogram_ptr);

        // Allocate grid nifti image
	dim[1]    = bkpr_size[0];
	dim[2]    = bkpr_size[1];
	dim[3]    = bkpr_size[2];
        nifti_image *gridImage = NULL;
        gridImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        gridImage->data = (float *)(grid_ptr);        

        // Allocate attenuation nifti image
        nifti_image *attenuationImage = NULL;
        if(!no_attenuation)
            {
            attenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            attenuationImage->data = (float *)(attenuation);        
            }

	// Allocate the result Fisher nifti image
        dim[1] = fisher_size[0];
        dim[2] = fisher_size[1];
        dim[3] = 1;
        nifti_image *fisherImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        fisherImage->data = (float *)(fisher_ptr);

        nifti_image *fisherpriorImage;
        if (fisher_prior_ptr!=NULL)
            {
            fisherpriorImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            fisherpriorImage->data = (float *)(fisher_prior_ptr);
            }

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

        //Compute Fisher Information Matrix
        #ifdef _USE_CUDA
        if (GPU)
            status = et_fisher_grid_gpu(from_projection, projectionImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, n_cameras, epsilon, background, background_attenuation); 
        else
            status = et_fisher_grid(from_projection, projectionImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, n_cameras, epsilon, background, background_attenuation); 
        #else
            if (GPU)
                fprintf(stderr,"et_array_project: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_fisher_grid(from_projection, projectionImage, gridImage, fisherImage, fisherpriorImage, psfImage, attenuationImage, cameras_array, n_cameras, epsilon, background, background_attenuation); 
        #endif

	//Free
	if( projectionImage->fname != NULL ) free(projectionImage->fname) ;
	if( projectionImage->iname != NULL ) free(projectionImage->iname) ;
	(void)nifti_free_extensions( projectionImage ) ;
	free(projectionImage) ;

	(void)nifti_free_extensions( gridImage ) ;
	free(gridImage) ;

	(void)nifti_free_extensions( fisherImage ) ;
	free(fisherImage) ;

        if (fisher_prior_ptr!=NULL)
            {         
            (void)nifti_free_extensions( fisherpriorImage ) ;
            free(fisherpriorImage) ;
            }

        if(!no_attenuation)
            {
            if( attenuationImage->fname != NULL ) free(attenuationImage->fname) ;
            if( attenuationImage->iname != NULL ) free(attenuationImage->iname) ;
            (void)nifti_free_extensions( attenuationImage ) ;
            free(attenuationImage) ;
            }

        if (!no_psf)
            {
            if( psfImage->fname != NULL ) free(psfImage->fname) ;
            if( psfImage->iname != NULL ) free(psfImage->iname) ;
            (void)nifti_free_extensions( psfImage ) ;
            free(psfImage) ;
            }

        free(cameras_array);

	return status;
}


extern "C" int et_array_convolve(float *image, int *image_size, float *out, int *out_size, float *psf, int *psf_size, int GPU)
{
        int status = 1;

	// Allocate source nifti image 
        int dim[8];
	dim[0]    = 3;//dims;            FIXME: bug in cudaCommon_transferNiftiToArrayOnDevice for 2D nifti images
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	//fprintf_verbose("\nS: %d %d %d",image_size[0],image_size[1],image_size[2]);
	nifti_image *imageImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        imageImage->data = (float *)(image);
	
	// Allocate the result nifti image
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];

        nifti_image *outImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        outImage->data = (float *)(out);

	//Allocate Point Spread Function nifti image
	//fprintf_verbose("\nPSF: %d %d %d",psf_size[0],psf_size[1],psf_size[2]);
	dim[1] = psf_size[0];
	dim[2] = psf_size[1];
	dim[3] = psf_size[2];
	nifti_image *psfImage;
        psfImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        psfImage->data = (float *)(psf);

        //Do convolution
        #ifdef _USE_CUDA
        if (GPU)
            status = et_convolve_gpu(imageImage, outImage, psfImage);
        else
            status = et_convolve(imageImage, outImage, psfImage);
        #else
            if (GPU)
                fprintf_verbose( "et_array_convolve: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
            status = et_convolve(imageImage, outImage, psfImage);
        #endif

	//Free
	if( imageImage->fname != NULL ) free(imageImage->fname) ;
	if( imageImage->iname != NULL ) free(imageImage->iname) ;
	(void)nifti_free_extensions( imageImage ) ;
	free(imageImage) ;
	
	if( outImage->fname != NULL ) free(outImage->fname) ;
	if( outImage->iname != NULL ) free(outImage->iname) ;
	(void)nifti_free_extensions( outImage ) ;
	free(outImage) ;

	if( psfImage->fname != NULL ) free(psfImage->fname) ;
	if( psfImage->iname != NULL ) free(psfImage->iname) ;
	(void)nifti_free_extensions( psfImage ) ;
	free(psfImage) ;

	return status;
}


/*
extern "C" int et_array_joint_histogram(float *matrix_A, float *matrix_B, int *joint_histogram, int matrix_dimensions, int *matrix_size, int histogram_size, float min_A, float max_A, float min_B, float max_B, int GPU)
{
	int status;

	// Allocate matrix_A and matrix_B nifti
        int dim[8] = {1,1,1,1,1,1,1,1};
	dim[0]    = matrix_dimensions;
	for(int i=1; i<=matrix_dimensions; i++)
		dim[i] = matrix_size[i-1];
	
	nifti_image *matrix_A_Image = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        matrix_A_Image->data = (float*) matrix_A;

	nifti_image *matrix_B_Image = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        matrix_B_Image->data = (float*) matrix_B;
                
	// Allocate the joint histogram
	for(int i=0; i<=8; i++)
		dim[i] = 1;
	dim[0]    = 2;
	dim[1]    = histogram_size;
	dim[2]    = histogram_size;
        nifti_image *joint_histogram_Image = nifti_make_new_nim(dim, NIFTI_TYPE_INT32, false);
        joint_histogram_Image->data = static_cast<void *>(joint_histogram);
        
	//Compute joint histogram
	#ifdef _USE_CUDA
	if(GPU)
	    status = et_joint_histogram_gpu(matrix_A_Image, matrix_B_Image, joint_histogram_Image, min_A, max_A, min_B, max_B);
	else
	    status = et_joint_histogram(matrix_A_Image, matrix_B_Image, joint_histogram_Image, min_A, max_A, min_B, max_B);
	#else
	    status = et_joint_histogram(matrix_A_Image, matrix_B_Image, joint_histogram_Image, min_A, max_A, min_B, max_B);
	#endif
	
	//Free
	if( matrix_A_Image->fname != NULL ) free(matrix_A_Image->fname) ;
	if( matrix_A_Image->iname != NULL ) free(matrix_A_Image->iname) ;
	(void)nifti_free_extensions( matrix_A_Image ) ;
	free(matrix_A_Image) ;
	
	if( matrix_B_Image->fname != NULL ) free(matrix_B_Image->fname) ;
	if( matrix_B_Image->iname != NULL ) free(matrix_B_Image->iname) ;
	(void)nifti_free_extensions( matrix_B_Image ) ;
	free(matrix_B_Image) ;

	if( joint_histogram_Image->fname != NULL ) free(joint_histogram_Image->fname) ;
	if( joint_histogram_Image->iname != NULL ) free(joint_histogram_Image->iname) ;
	(void)nifti_free_extensions( joint_histogram_Image ) ;
	free(joint_histogram_Image) ;

	return status;
}
*/


extern "C" int et_array_list_gpus(int *gpu_count, int *gpus_info_array)
{
        int status = 1;
	#ifdef _USE_CUDA
        status = et_list_gpus(gpu_count, gpus_info_array);
	#else
        gpu_count[0] = 0;
        status = 0;
	#endif
        return status;
}



extern "C" int et_array_set_gpu(int id)
{
        int status = 1;
	#ifdef _USE_CUDA
        status = et_set_gpu(id);
	#else
        status = 1;
	#endif
        return status;
}



