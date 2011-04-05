
#include "_et.h"
#include "_et_common.h"
#include<stdio.h>



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   CPU   ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int et_affine(nifti_image *sourceImage, nifti_image *transformedImage, mat44 *affineTransformation, float background)
{
    int status = 1;
    /* Allocate the deformation Field image */
    nifti_image *positionFieldImage = nifti_copy_nim_info(sourceImage);
    positionFieldImage->dim[0]=positionFieldImage->ndim=5;
    positionFieldImage->dim[1]=positionFieldImage->nx=sourceImage->nx;
    positionFieldImage->dim[2]=positionFieldImage->ny=sourceImage->ny;
    positionFieldImage->dim[3]=positionFieldImage->nz=sourceImage->nz;
    positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
    positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
    positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
    positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
    positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
    positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
    positionFieldImage->nbyper = sizeof(float);
    positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
	
    /* Apply affine */
    reg_affine_positionField(       affineTransformation,
                                    sourceImage,
                                    positionFieldImage );
    /* Resample the source image */
    reg_resampleSourceImage<float>( sourceImage,
                                    sourceImage,
                                    transformedImage,
                                    positionFieldImage,
                                    NULL,
                                    1,
                                    background);
    nifti_image_free(positionFieldImage);
    status = 0;
    return status;
}



int et_rotate(nifti_image *sourceImage, nifti_image *resultImage, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, float background)
{
	int status = 1;
	
	//Create transformation matrix
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
	et_create_rotation_matrix(affineTransformation, theta_x, theta_y, theta_z, center_x, center_y, center_z);
	
	//Apply affine transformation
	status = et_affine(sourceImage, resultImage, affineTransformation, background);

	//Free
	free(affineTransformation);

	return status;
}



int et_project(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
{
        /* Check consistency of input */
        //...
        //fprintf_verbose("\n ET_PROJECT: %d %d %d, %d %d %d, %d %d %d",activityImage->nx,activityImage->ny,activityImage->nz,sinoImage->nx,sinoImage->ny,sinoImage->nz,psfImage->nx,psfImage->ny,psfImage->nz);
	/* Allocate the deformation Field image */
	nifti_image *positionFieldImage = nifti_copy_nim_info(activityImage);
	positionFieldImage->dim[0]=positionFieldImage->ndim=5;
	positionFieldImage->dim[1]=positionFieldImage->nx=activityImage->nx;
	positionFieldImage->dim[2]=positionFieldImage->ny=activityImage->ny;
	positionFieldImage->dim[3]=positionFieldImage->nz=activityImage->nz;
	positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
	positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
	positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
	positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
	positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
	positionFieldImage->nbyper = sizeof(float);
	positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
	
	/* Allocate arrays */
        int dim[8];
	dim[0]    = 3;
	dim[1]    = activityImage->dim[1];
	dim[2]    = activityImage->dim[2];
	dim[3]    = activityImage->dim[3];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *rotatedImage;
        rotatedImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        rotatedImage->data = (int *)malloc(activityImage->nvox*sizeof(int));	

	/* Define centers of rotation */
	float center_x = ((float)(activityImage->nx - 1)) / 2.0;
	float center_y = ((float)(activityImage->ny - 1)) / 2.0;
	float center_z = ((float)(activityImage->nz - 1)) / 2.0;
		
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
	
	for(int cam=0; cam<n_cameras; cam++){
		// Apply affine //
                fprintf_verbose( "et_project: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
		et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z);
		reg_affine_positionField(	affineTransformation,
						activityImage,
						positionFieldImage );
		// Resample the source image //
		reg_resampleSourceImage<float>(	activityImage,
						activityImage,
						rotatedImage,
						positionFieldImage,
						NULL,
						1,
						background );	

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    {
                    //et_convolveFFT2D(         rotatedImage,
                    //                          image_size;
                    //                          psfImage,
                    //                          psf_size,
                    //                          rotatedImage);
                    }
	
		// Integrate along lines //
		et_line_integral(		rotatedImage,
						sinoImage,
						cam );
	}

	/*Free*/
	nifti_image_free(rotatedImage);
	nifti_image_free(positionFieldImage);
	free(affineTransformation);

	return 0;
}



int et_backproject(nifti_image *sinogramImage, nifti_image *accumulatorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
{
        /* Check consistency of input */
        //...

	/* Allocate the deformation Field image */
	nifti_image *positionFieldImage = nifti_copy_nim_info(accumulatorImage);
	positionFieldImage->dim[0]=positionFieldImage->ndim=5;
	positionFieldImage->dim[1]=positionFieldImage->nx=accumulatorImage->nx;
	positionFieldImage->dim[2]=positionFieldImage->ny=accumulatorImage->ny;
	positionFieldImage->dim[3]=positionFieldImage->nz=accumulatorImage->nz;
	positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
	positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
	positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
	positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
	positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
	positionFieldImage->nbyper = sizeof(float);
	positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
	
	/* Allocate arrays */
        int dim[8];
	dim[0]    = 3;
	dim[1]    = accumulatorImage->dim[1];
	dim[2]    = accumulatorImage->dim[2];
	dim[3]    = accumulatorImage->dim[3];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *rotatedImage;
        rotatedImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        rotatedImage->data = (int *)malloc(accumulatorImage->nvox*sizeof(int));	

	nifti_image *temp_backprojectionImage;
        temp_backprojectionImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        temp_backprojectionImage->data = (int *)malloc(accumulatorImage->nvox*sizeof(int));	        

	/* Define centers of rotation */
	float center_x = ((float)(accumulatorImage->nx - 1)) / 2.0;
	float center_y = ((float)(accumulatorImage->ny - 1)) / 2.0;
	float center_z = ((float)(accumulatorImage->nz - 1)) / 2.0;

	/* Clear accumulator */
	et_clear_accumulator(accumulatorImage);	
	
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
	
	for(int cam=0; cam<n_cameras; cam++){

		/* Line Backproject */
		et_line_backproject(		sinogramImage,
						temp_backprojectionImage,
						cam );
		
		/* Rotate */
		et_create_rotation_matrix(	affineTransformation,
						-cameras[0*n_cameras+cam],
						-cameras[1*n_cameras+cam],
						-cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z);
						
		reg_affine_positionField(	affineTransformation,
						accumulatorImage,
						positionFieldImage);

		reg_resampleSourceImage(	temp_backprojectionImage,
						temp_backprojectionImage,
						rotatedImage,
						positionFieldImage,
						NULL,
						1,
						background);
		
		//fprintf_verbose("\n>> %d %d %d %d ",accumulatorImage->nx, accumulatorImage->ny, accumulatorImage->nz, accumulatorImage->nvox);
		
		/* Accumulate */
		et_accumulate(			rotatedImage,
						accumulatorImage );
	}

        //for (int i=1; i<accumulatorImage->nvox; i++)
        //    ((float*)accumulatorImage->data)[i] = ((float*)rotatedImage->data)[i];
	/*Free*/
	nifti_image_free(rotatedImage);
	nifti_image_free(temp_backprojectionImage);
	nifti_image_free(positionFieldImage);
	free(affineTransformation);

	return 0;
}




int et_project_backproject(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_alpha, float *cameras_beta, float *cameras_gamma)
{
    return 0;
}



int et_convolve(nifti_image *inImage, nifti_image *outImage, nifti_image *psfImage)
{
    int status = 1;
    return status;
}


/*
int et_joint_histogram(nifti_image *matrix_A_Image, nifti_image *matrix_B_Image, nifti_image *joint_histogram_Image, float min_A, float max_A, float min_B, float max_B)
{
    return 0;
}
*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   GPU   ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _USE_CUDA

int et_affine_gpu(nifti_image *sourceImage, nifti_image *resultImage, mat44 *affineTransformation, float background)
{
	/* initialise the cuda arrays */
	cudaArray *sourceImageArray_d;
	float     *resultImageArray_d;
	float4    *positionFieldImageArray_d;
	int       *mask_d;

	/* Allocate the deformation Field image */
	nifti_image *positionFieldImage = nifti_copy_nim_info(sourceImage);
	positionFieldImage->dim[0]=positionFieldImage->ndim=5;
	positionFieldImage->dim[1]=positionFieldImage->nx=sourceImage->nx;
	positionFieldImage->dim[2]=positionFieldImage->ny=sourceImage->ny;
	positionFieldImage->dim[3]=positionFieldImage->nz=sourceImage->nz;
	positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
	positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
	positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
	positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
	positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
	positionFieldImage->nbyper = sizeof(float);
	positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
	
	/* Allocate arrays on the device */
	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&resultImageArray_d, resultImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, resultImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, resultImage->dim)) return 1;

	/* Transfer data from the host to the device */
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;
	int *mask_h=(int *)malloc(resultImage->nvox*sizeof(int));
	for(int i=0; i<resultImage->nvox; i++) mask_h[i]=i;
	CUDA_SAFE_CALL(cudaMemcpy(mask_d, mask_h, resultImage->nvox*sizeof(int), cudaMemcpyHostToDevice));
	free(mask_h);

	/* Apply affine */
	reg_affine_positionField_gpu(	affineTransformation,
					resultImage,
					&positionFieldImageArray_d);
	/* Resample the source image */
	reg_resampleSourceImage_gpu(	resultImage,
					sourceImage,
					&resultImageArray_d,
					&sourceImageArray_d,
					&positionFieldImageArray_d,
					&mask_d,
					resultImage->nvox,
					background);
	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(resultImage, &resultImageArray_d)) return 1;
	/*Free*/
	cudaCommon_free(&sourceImageArray_d);
	cudaCommon_free((void **)&resultImageArray_d);
	cudaCommon_free((void **)&mask_d);
	cudaCommon_free((void **)&positionFieldImageArray_d);
	nifti_image_free(positionFieldImage);

	return 0;
}



int et_rotate_gpu(nifti_image *sourceImage, nifti_image *resultImage, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, float background)
{
	int status = 1;
	
	//Create transformation matrix
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
	et_create_rotation_matrix(affineTransformation, theta_x, theta_y, theta_z, center_x, center_y, center_z);
	
	//Apply affine transformation
	status = et_affine_gpu(sourceImage, resultImage, affineTransformation, background);

	//Free
	free(affineTransformation);

	return status;
}



int et_project_gpu(nifti_image *activity, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
{
	/* initialise the cuda arrays */
	cudaArray *activityArray_d;               //stores input activity, makes use of fetch unit
        cudaArray *attenuationArray_d;            //stores input attenuation coefficients, makes use of fetch unit
	float     *sinoArray_d;                   //stores sinogram (output)
	float     *rotatedArray_d;                //stores activity aligned with current camera
        float     *rotatedAttenuationArray_d;     //stores attenuation coefficients aligned with current camera
	float4    *positionFieldImageArray_d;     //stores position field for rotation of activity and attenuation
	int       *mask_d;                        //binary mask that defines active voxels (typically all active)
        float     *psfArray_d;                    //stores point spread function
        int       psf_size[3];
        int       image_size[3];

        /* Check consistency of input */
        //..
	
	/* Allocate arrays on the device */
	if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activity->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, activity->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, activity->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, activity->dim)) return 1;

	/* Transfer data from the host to the device */
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activity)) return 1;
	int *mask_h=(int *)malloc(activity->nvox*sizeof(int));
	for(int i=0; i<activity->nvox; i++) mask_h[i]=i;
	CUDA_SAFE_CALL(cudaMemcpy(mask_d, mask_h, activity->nvox*sizeof(int), cudaMemcpyHostToDevice));
	free(mask_h);
	
	/* Define centers of rotation */
	float center_x = ((float)(activity->nx - 1)) / 2.0;
	float center_y = ((float)(activity->ny - 1)) / 2.0;
	float center_z = ((float)(activity->nz - 1)) / 2.0;
		
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) return 1;
            image_size[0] = activity->dim[1];
            image_size[1] = activity->dim[2];
            image_size[2] = activity->dim[3];
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            }
        if (attenuationImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) return 1;
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) return 1;
            }

	for(unsigned int cam=0; cam<n_cameras; cam++){

		// Apply affine //
		et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z);
		reg_affine_positionField_gpu(	affineTransformation,
						activity,
						&positionFieldImageArray_d);

		// Resample the source image //
		reg_resampleSourceImage_gpu(	activity,
						activity,
						&rotatedArray_d,
						&activityArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						activity->nvox,
						background);

                // Resample the attenuation map //
                if (attenuationImage != NULL)
		    reg_resampleSourceImage_gpu(	attenuationImage,
						attenuationImage,
						&rotatedAttenuationArray_d,
						&attenuationArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						attenuationImage->nvox,
						background_attenuation);

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    et_convolveFFT2D_gpu(       &rotatedArray_d,
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &rotatedArray_d);

		// Integrate along lines //
                if (attenuationImage != NULL)
                    et_line_integral_attenuated_gpu(	&rotatedArray_d,
						&rotatedAttenuationArray_d, 
						&sinoArray_d,
						cam,
						activity);
                else
                    {
		    et_line_integral_gpu(	&rotatedAttenuationArray_d,
						&sinoArray_d,
						cam,
						activity);
                    }
	}


	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(sinoImage, &sinoArray_d)) return 1;

	/*Free*/
	cudaCommon_free(&activityArray_d);
	cudaCommon_free((void **)&rotatedArray_d);
	cudaCommon_free((void **)&sinoArray_d);
	cudaCommon_free((void **)&mask_d);
	cudaCommon_free((void **)&positionFieldImageArray_d);
        if (psfImage != NULL)
            cudaCommon_free((void **)&psfArray_d);
	free(affineTransformation);
        if (attenuationImage != NULL)
            {
            cudaCommon_free(&attenuationArray_d);
            cudaCommon_free((void **)&rotatedAttenuationArray_d);
            }
	return 0;
}


int et_backproject_gpu(nifti_image *sinoImage, nifti_image *accumulatorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
{
	/* initialise the cuda arrays */
	cudaArray *backprojectionArray_d;
	cudaArray *attenuationArray_d;
	float     *temp_backprojection_d;
	float     *sinoArray_d;
	float     *rotatedArray_d;
	float     *rotatedAttenuationArray_d;
        float     *attenuationPlaneArray_d;
	float     *accumulatorArray_d;
	float4    *positionFieldImageArray_d;
	int       *mask_d;
        float     *psfArray_d;
        int       psf_size[3];
        int       image_size[3];

	/* Allocate the deformation Field image */
	nifti_image *positionFieldImage = nifti_copy_nim_info(accumulatorImage);
	positionFieldImage->dim[0]=positionFieldImage->ndim=5;
	positionFieldImage->dim[1]=positionFieldImage->nx = accumulatorImage->nx;
	positionFieldImage->dim[2]=positionFieldImage->ny = accumulatorImage->ny;
	positionFieldImage->dim[3]=positionFieldImage->nz = accumulatorImage->nz;
	positionFieldImage->dim[4]=positionFieldImage->nt = 1; positionFieldImage->pixdim[4]=positionFieldImage->dt = 1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu = 3; positionFieldImage->pixdim[5]=positionFieldImage->du = 1.0;
	positionFieldImage->dim[6]=positionFieldImage->nv = 1; positionFieldImage->pixdim[6]=positionFieldImage->dv = 1.0;
	positionFieldImage->dim[7]=positionFieldImage->nw = 1; positionFieldImage->pixdim[7]=positionFieldImage->dw = 1.0;
	positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
	positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
	positionFieldImage->nbyper = sizeof(float);
	positionFieldImage->data=NULL;
	
	/* Allocate arrays on the device */

        cudaChannelFormatDesc backprojectionArray_d_chdesc = cudaCreateChannelDesc<float>();
        cudaExtent backprojectionArray_d_extent;
        backprojectionArray_d_extent.width  = accumulatorImage->nx;
        backprojectionArray_d_extent.height = accumulatorImage->ny;
        backprojectionArray_d_extent.depth  = accumulatorImage->nz;
        cudaError_t cuda_status1 = cudaMalloc3DArray(&backprojectionArray_d, &backprojectionArray_d_chdesc, backprojectionArray_d_extent);

        if(cudaCommon_allocateArrayToDevice<float>(&temp_backprojection_d, accumulatorImage->dim)) return 1;

        if (attenuationImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) return 1;
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) return 1;	
//            if(cudaCommon_allocateArrayToDevice<float>(&attenuationPlaneArray_d, attenuationImage->dim)) return 1;	//FIXME size
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) return 1;
            }

	if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, accumulatorImage->dim)) return 1;	
	if(cudaCommon_allocateArrayToDevice<float>(&accumulatorArray_d, accumulatorImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, accumulatorImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, accumulatorImage->dim)) return 1;

	/* Transfer data from the host to the device */
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sinoArray_d,sinoImage)) return 1;
	int *mask_h=(int *)malloc(accumulatorImage->nvox*sizeof(int));
	for(int i=0; i<accumulatorImage->nvox; i++) mask_h[i]=i;
	CUDA_SAFE_CALL(cudaMemcpy(mask_d, mask_h, accumulatorImage->nvox*sizeof(int), cudaMemcpyHostToDevice));
	free(mask_h);

	/* Define centers of rotation */
	float center_x = ((float)(accumulatorImage->nx - 1)) / 2.0;
	float center_y = ((float)(accumulatorImage->ny - 1)) / 2.0;
	float center_z = ((float)(accumulatorImage->nz - 1)) / 2.0;

	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));

	/* Clear accumulator */
	et_clear_accumulator_gpu(		&accumulatorArray_d,
						accumulatorImage );

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) return 1;
            image_size[0] = accumulatorImage->dim[1];
            image_size[1] = accumulatorImage->dim[2];
            image_size[2] = accumulatorImage->dim[3];
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            }


		
	for(int cam=0; cam<n_cameras; cam++){

                // Rotate attenuation //
                if (attenuationImage != NULL)
                    et_create_rotation_matrix(	affineTransformation,
						cameras[0*n_cameras+cam],
						cameras[1*n_cameras+cam],
						cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z);
                    reg_affine_positionField_gpu(affineTransformation,
						attenuationImage,
						&positionFieldImageArray_d);
		    reg_resampleSourceImage_gpu(attenuationImage,
						attenuationImage,
						&rotatedAttenuationArray_d,
						&attenuationArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						attenuationImage->nvox,
						background_attenuation);

		// Line Backproject //
                if (attenuationImage != NULL)
                    {
                    et_line_backproject_attenuated_gpu(	&sinoArray_d,
						&temp_backprojection_d,
						&rotatedAttenuationArray_d,
                                                cam,
						accumulatorImage);
                    }
                else
                    et_line_backproject_gpu(	&sinoArray_d,
						&temp_backprojection_d,
						cam,
						accumulatorImage);

                // Copy to texture bound memory (for rotation) //
                cudaError_t cuda_status;
                cudaMemcpy3DParms p = {0};
              
                p.srcPtr.ptr        = temp_backprojection_d;
                p.srcPtr.pitch      = 0;
                p.srcPtr.xsize      = accumulatorImage->nx;
                p.srcPtr.ysize      = accumulatorImage->ny;
                p.dstArray          = backprojectionArray_d;
                p.extent.width      = accumulatorImage->nx;
                p.extent.height     = accumulatorImage->ny;
                p.extent.depth      = accumulatorImage->nz;
                p.kind              = cudaMemcpyDeviceToDevice;
                cuda_status         = cudaMemcpy3D(&p);

		// Rotate backprojection //
		et_create_rotation_matrix(	affineTransformation,
						-cameras[0*n_cameras+cam],
						-cameras[1*n_cameras+cam],
						-cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z);
		reg_affine_positionField_gpu(	affineTransformation,
						accumulatorImage,
						&positionFieldImageArray_d);
		reg_resampleSourceImage_gpu(	accumulatorImage,
						accumulatorImage,
						&rotatedArray_d,
						&backprojectionArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						accumulatorImage->nvox,
						background);

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    {
                    et_convolveFFT2D_gpu(       &rotatedArray_d,
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                    }

		// Accumulate //
		et_accumulate_gpu(		&rotatedArray_d,
						&accumulatorArray_d,
						accumulatorImage );
	}

	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(accumulatorImage, &accumulatorArray_d)) return 1; 

	/*Free*/
        if (attenuationImage != NULL)
            {
            cudaCommon_free(&backprojectionArray_d);
            cudaCommon_free(&attenuationArray_d);
            cudaCommon_free((void **)&rotatedAttenuationArray_d);
//         cudaCommon_free((void **)&attenuationPlaneArray_d);
            }
	cudaCommon_free((void **)&rotatedArray_d);
	cudaCommon_free((void **)&sinoArray_d);
	cudaCommon_free((void **)&accumulatorArray_d);
	cudaCommon_free((void **)&mask_d);
	cudaCommon_free((void **)&positionFieldImageArray_d);
	cudaCommon_free((void **)&temp_backprojection_d);
        if (psfImage != NULL)
            cudaCommon_free((void **)&psfArray_d);
	nifti_image_free(positionFieldImage);
	free(affineTransformation);

	return 0;
}


int et_project_backproject_gpu(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_theta_x, float *cameras_theta_y, float *cameras_theta_z)
{
	int status = 1;
	//status = rotate();
	return status;
}



int et_convolve_gpu(nifti_image *inImage, nifti_image *outImage, nifti_image *psfImage)
{
    int status = 1;

    int data_size[2];
    int kernel_size[2];

    float *d_Input;
    float *d_Kernel;
    float *d_Result;

    fprintf_verbose("et_convolve_gpu - Image size: %d %d %d\n",inImage->dim[1],inImage->dim[2],inImage->dim[3]);

    fprintf_verbose("et_convolve_gpu - Allocating memory...\n");
    if(cudaCommon_allocateArrayToDevice<float>(&d_Input, inImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float>(&d_Kernel, psfImage->dim)) return 1;
    if(cudaCommon_allocateArrayToDevice<float>(&d_Result, outImage->dim)) return 1;

    fprintf_verbose("et_convolve_gpu - Uploading to GPU device...\n");
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&d_Input, inImage)) return 1;
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&d_Kernel, psfImage)) return 1;
    if(cudaCommon_transferNiftiToArrayOnDevice<float>(&d_Result, outImage)) return 1;

    fprintf_verbose("et_convolve_gpu - Performing 2D convolution...\n");
    data_size[0] = inImage->dim[1];
    data_size[1] = inImage->dim[2];
    data_size[2] = inImage->dim[3];
    kernel_size[0] = psfImage->dim[1];
    kernel_size[1] = psfImage->dim[2];
    kernel_size[2] = psfImage->dim[3];
    
    et_convolveFFT2D_gpu(&d_Input, data_size, &d_Kernel, kernel_size, &d_Result);

    fprintf_verbose("et_convolve_gpu - Reading back GPU FFT results...\n");
    //cutilSafeCall( cudaMemcpy(h_Result, d_PaddedData, fftH * fftW * sizeof(float), cudaMemcpyDeviceToHost) );
    if(cudaCommon_transferFromDeviceToNifti(outImage, &d_Result)) return 1; 

    fprintf_verbose("et_convolve_gpu - Freeign GPU memory...\n");    
    cutilSafeCall( cudaFree(d_Input) );
    cutilSafeCall( cudaFree(d_Kernel) );
    cutilSafeCall( cudaFree(d_Result) );
    fprintf_verbose("et_convolve_gpu - Freeign GPU memory...\n");   

    status = 0;
    return status;
}


/*
int et_joint_histogram_gpu(nifti_image *matrix_A_Image, nifti_image *matrix_B_Image, nifti_image *joint_histogram_Image, float min_A, float max_A, float min_B, float max_B)
{
	float *matrix_A_d;
	float *matrix_B_d;
	int *joint_histogram_d;
	
	// Allocate arrays on device 
	CUDA_SAFE_CALL(cudaMalloc((void **)&matrix_A_d,  matrix_A_Image->nvox*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&matrix_B_d,  matrix_B_Image->nvox*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&joint_histogram_d,  joint_histogram_Image->nvox*sizeof(int)));
        cudaMemset((void*)joint_histogram_d,0,joint_histogram_Image->nvox*sizeof(int));
        
	// Transfer data from the host to the device 
	CUDA_SAFE_CALL(cudaMemcpy(matrix_A_d, matrix_A_Image->data, matrix_A_Image->nvox*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(matrix_B_d, matrix_B_Image->data, matrix_B_Image->nvox*sizeof(float), cudaMemcpyHostToDevice));

	// Compute joint histogram 
	et_joint_histogram_gpu(		&matrix_A_d,
					&matrix_B_d,
					&joint_histogram_d,
					matrix_A_Image->nvox,
					joint_histogram_Image->nx,
					min_A,
					max_A,
					min_B,
					max_B );
	
	// Transfer joint histogram back to host 
	fprintf_verbose( "\nJH: %d %d %d %f %f %f %f",joint_histogram_Image->nx, joint_histogram_Image->nvox, matrix_A_Image->nvox, min_A, max_A, min_B, max_B);
	CUDA_SAFE_CALL(cudaMemcpy(joint_histogram_Image->data, joint_histogram_d, joint_histogram_Image->nvox*sizeof(int), cudaMemcpyDeviceToHost));
	
	// Free arrays on device 
	cudaCommon_free((void **)&matrix_A_d);
	cudaCommon_free((void **)&matrix_B_d);
	cudaCommon_free((void **)&joint_histogram_d);
	
	return 0;
}
*/


int et_list_gpus(int *device_count_out, int *devices)
{
    /* Initialise the cuda card */
    int status = 1;
    struct cudaDeviceProp deviceProp;
    int device_count = 0;
    int multiprocessors = 0;
    int clock = 0;
    int gflops = 0;
    int globalmem = 0;

    cudaGetDeviceCount( &device_count );

    int device_count_max = device_count;
    if (device_count_max > MAX_DEVICES)
        device_count_max = MAX_DEVICES;

    for (int dev=0; dev<device_count_max; dev++)
        {
        cudaGetDeviceProperties(&deviceProp, dev);
        multiprocessors = deviceProp.multiProcessorCount;
        clock = deviceProp.clockRate;
        gflops = multiprocessors * clock;
        globalmem = (int)floor(deviceProp.totalGlobalMem/1000000.0);
        //fprintf_verbose("\nDevice %d: %d MP, %d GHz, %d GFlops, %d Mb",dev,multiprocessors,clock,gflops,globalmem);
        devices[SIZE_OF_INFO*dev+0] = dev;
        devices[SIZE_OF_INFO*dev+1] = gflops;
        devices[SIZE_OF_INFO*dev+2] = multiprocessors;
        devices[SIZE_OF_INFO*dev+3] = clock;
        devices[SIZE_OF_INFO*dev+4] = globalmem;
        }
    device_count_out[0] = device_count;
    status = 0;
    return status;
}


int et_set_gpu(int id)
{
    int status = 1;
    struct cudaDeviceProp deviceProp;

    CUDA_SAFE_CALL(cudaSetDevice( id ));
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, id ));
    if (deviceProp.major < 1)
        {
        printf("ERROR\tThe specified graphical card does not exist.\n");
        status = 1;
	}
    else
        status = 0;
    return status;
}



#endif





