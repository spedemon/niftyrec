
#include "_tt.h"
#include "_et_common.h"


int tt_project(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background)
{
    int status = 1;
    status = 0;
    return status;
}


#ifdef _USE_CUDA
int tt_project_gpu(nifti_image *attenuation, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background)
{
	/* initialise the cuda arrays */
        cudaArray *attenuationArray_d;            //stores input attenuation coefficients, makes use of fetch unit
	float     *projectionArray_d;             //stores projections (output)
        float     *perspectiveAttenuationArray_d; //stores attenuation coefficients aligned with current camera
	float4    *positionFieldImageArray_d;     //stores position field for rotation of attenuation
	int       *mask_d;                        //binary mask that defines active voxels (typically all active)
        float     *psfArray_d;                    //stores point spread function
        int       psf_size[3];
        int       image_size[3];

	/* Allocate arrays on the device */
	if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuation->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&projectionArray_d, projectionImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&perspectiveAttenuationArray_d, attenuation->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, attenuation->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, attenuation->dim)) return 1;

	/* Transfer data from the host to the device */
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuation)) return 1;
	int *mask_h=(int *)malloc(attenuation->nvox*sizeof(int));
	for(int i=0; i<attenuation->nvox; i++) mask_h[i]=i;
	CUDA_SAFE_CALL(cudaMemcpy(mask_d, mask_h, attenuation->nvox*sizeof(int), cudaMemcpyHostToDevice));
	free(mask_h);

	/* Alloc transformation matrix */
	mat44 *transformationMatrix = (mat44 *)calloc(1,sizeof(mat44));

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) return 1;
            image_size[0] = attenuation->dim[1];
            image_size[1] = attenuation->dim[2];
            image_size[2] = attenuation->dim[3];
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            }

	for(unsigned int projection=0; projection<n_projections; projection++){

		// Compute deformation field //

//		et_create_rotation_matrix(transformationMatrix, 0,0,0,0,0,0);

//		reg_affine_positionField_gpu(	transformationMatrix,
//						attenuation,
//						&positionFieldImageArray_d);

                tt_perspective_positionField_gpu( attenuation,
                                                  image_origin, 
                                                  detector_origin,
                                                  detector_shape,
                                                  &positionFieldImageArray_d);
fprintf(stderr, "iO: %f %f %f\n dO: %f %f %f\n dS: %f %f \n",image_origin[0],image_origin[1],image_origin[2],detector_origin[0],detector_origin[1],detector_origin[2], detector_shape[0], detector_shape[1]);

		// Resample attenuation image //
		reg_resampleSourceImage_gpu(	attenuation,
						attenuation,
						&perspectiveAttenuationArray_d,
						&attenuationArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						attenuation->nvox,
						background);

		// Correct for non-uniform voxel size //

                // Apply Depth Dependent Point Spread Function to attenuation //
//                if (psfImage != NULL)
//                    {
//                    et_convolveFFT2D_gpu(       &perspectiveArray_d,
//                                                image_size,
//                                                &psfArray_d,
//                                                psf_size,
//                                                &perspectiveArray_d);
//                    }

//	float *temp_h=(float *)malloc(attenuation->nvox*sizeof(float));
//	for(int i=0; i<attenuation->nvox; i++) temp_h[i]=(float)9;
//	CUDA_SAFE_CALL(cudaMemcpy(perspectiveAttenuationArray_d, temp_h, attenuation->nvox*sizeof(float), cudaMemcpyHostToDevice));

		// Integrate along lines //
                et_line_integral_gpu(		&perspectiveAttenuationArray_d,
						&projectionArray_d,
						projection,
						attenuation);
	}

	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(projectionImage, &projectionArray_d)) return 1;

	/*Free*/
	cudaCommon_free(&attenuationArray_d);
	cudaCommon_free((void **)&perspectiveAttenuationArray_d);
	cudaCommon_free((void **)&projectionArray_d);
	cudaCommon_free((void **)&mask_d);
	cudaCommon_free((void **)&positionFieldImageArray_d);
	free(transformationMatrix);
        if (psfImage != NULL)
            cudaCommon_free((void **)&psfArray_d);
	return 0;
}
#endif

