/*
 *  _et.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et.h"
#include "_et_common.h"
#include <stdio.h>
#include <time.h>

#define max(a,b)	(((a) > (b)) ? (a) : (b))
#define min(a,b)	(((a) < (b)) ? (a) : (b))


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   CPU   ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p)
{
  return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
           ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}


int et_is_block_multiple(int size)
{
    if (size % (ET_BLOCK_SIZE*ET_BLOCK_SIZE) == 0)
        return 1;
    return 0;
}
int et_get_block_size(void)
{
    return ET_BLOCK_SIZE*ET_BLOCK_SIZE;
}


//! Affine transformation of nifti_image 
/*!
  \param *sourceImage the source image to be transformed. 
  \param *transformedImage the transformed image. 
  \param *affineTransformation the [4x4] transformed matrix. 
  \param background the background value when resampling the transformed image. 
*/
int et_affine(nifti_image *sourceImage, nifti_image *transformedImage, mat44 *affineTransformation, float background)
{
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
    return 0;
}


//! Rotate a nifti_image in 3D 
/*!
  \param *sourceImage the source image to be transformed. 
  \param *resultImage the rotated image. 
  \param theta_x the rotation angle around x axis in radians. 
  \param theta_y the rotation angle around y axis in radians. 
  \param theta_z the rotation angle around z axis in radians. 
  \param center_x the center of rotation along x.  
  \param center_y the center of rotation along y. 
  \param center_z the center of rotation along z. 
  \param background the background value when resampling the transformed image. 
*/
int et_rotate(nifti_image *sourceImage, nifti_image *resultImage, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, float background)
{
        int status; 
	//Create transformation matrix
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
	et_create_rotation_matrix(affineTransformation, theta_x, theta_y, theta_z, center_x, center_y, center_z, XYZ_ROTATION);
	
	//Apply affine transformation
	status = et_affine(sourceImage, resultImage, affineTransformation, background);

	//Free
	free(affineTransformation);

	return status;
}



//! Projection for Emission Imaging
/*!
  \param *activityImage the activity (or its estimate). NULL for attenuation and background activity only. 
  \param *sinoImage the photon counts in projection space. 
  \param *psfImage the depth-dependent point spread function, NULL for no PSF. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_project(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{
        int separable_psf = 0;
        int psf_size[3];

        /* Check consistency of input */
        nifti_image *referenceImage;              // this image holds information about image size and voxel size (activity or attenuation might not be defined (NULL pointers) ) 
        if (activityImage==NULL && attenuationImage==NULL)
            {
            fprintf(stderr, "et_project: Error - define at least one between activityImage and attenuationImage. \n");
            return niftyrec_error_parameters; 
            }
        else if (attenuationImage==NULL)
            referenceImage=activityImage;
        else
            referenceImage=attenuationImage; 

	/* Allocate the deformation Field image */
        alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS); 
	nifti_image *positionFieldImage = nifti_copy_nim_info(referenceImage);
	positionFieldImage->dim[0]=positionFieldImage->ndim=5;
	positionFieldImage->dim[1]=positionFieldImage->nx=referenceImage->nx;
	positionFieldImage->dim[2]=positionFieldImage->ny=referenceImage->ny;
	positionFieldImage->dim[3]=positionFieldImage->nz=referenceImage->nz;
	positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
	positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
	positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
	positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
	positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
	positionFieldImage->nbyper = sizeof(float);
	positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper);
        if (positionFieldImage->data==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
        alloc_record_add(memory_record,(void*)positionFieldImage,ALLOCTYPE_NIFTI);
	
	/* Allocate arrays */
        int dim[8];
	dim[0]    = 3;
	dim[1]    = referenceImage->dim[1];
	dim[2]    = referenceImage->dim[2];
	dim[3]    = referenceImage->dim[3];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *rotatedImage=NULL;
        if (activityImage!=NULL)
            {
            rotatedImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedImage->data = (float *)malloc(referenceImage->nvox*sizeof(float));	
            if (rotatedImage->data==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
            alloc_record_add(memory_record,(void*)rotatedImage,ALLOCTYPE_NIFTI);
            }

	nifti_image *rotatedAttenuationImage=NULL;
        if (attenuationImage != NULL)
            {
            rotatedAttenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedAttenuationImage->data = (float *)malloc(referenceImage->nvox*sizeof(float));
            if (rotatedAttenuationImage->data==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
            alloc_record_add(memory_record,(void*)rotatedAttenuationImage,ALLOCTYPE_NIFTI);	
            }

	/* Define centers of rotation */
	float center_x = ((float)(referenceImage->nx - 1)) / 2.0;
	float center_y = ((float)(referenceImage->ny - 1)) / 2.0;
	float center_z = ((float)(referenceImage->nz - 1)) / 2.0;
		
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        if (affineTransformation==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
        alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_GUEST);

        /* Decide whether to use FFT or separate the convolution */
        float *psfSeparated=NULL;
        float psf_norm;
	if(psfImage!=NULL)
            {
            if (psfImage->nx <= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1) 
                {
                separable_psf=1;
                psf_size[0] = psfImage->dim[1];
                psf_size[1] = psfImage->dim[2];
                psf_size[2] = psfImage->dim[3];
                psfSeparated = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
                if (psfSeparated==NULL) {alloc_record_destroy(memory_record); return niftyrec_error_alloccpu;}; 
                alloc_record_add(memory_record,(void*)psfSeparated,ALLOCTYPE_GUEST);
                for (int n=0; n<psf_size[2];n++) 
                    {
                    psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                    for (int i=0;i<psf_size[0];i++) 
                        {
                        psfSeparated[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                        psfSeparated[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                        }
                    }
                }
            }

        /* Project */
	for(int cam=0; cam<n_cameras; cam++){
		// Apply affine //
                fprintf_verbose( "et_project: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
		et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z, XYZ_ROTATION);
		reg_affine_positionField(	affineTransformation,
						referenceImage,
						positionFieldImage );
		// Resample the source image //
                if (activityImage != NULL)
                    reg_resampleSourceImage<float>(activityImage,
						activityImage,
						rotatedImage,
						positionFieldImage,
						NULL,
						1,
						background );	

                // Resample the attenuation map //
                if (attenuationImage != NULL)
                    {
                    reg_resampleSourceImage<float>(attenuationImage,
						attenuationImage,
						rotatedAttenuationImage,
						positionFieldImage,
						NULL,
						1,
						background_attenuation );	                    
                    }

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL && activityImage != NULL)
                    {
                    if (separable_psf)
                        et_convolveSeparable2D( rotatedImage,
                                                psfSeparated,
                                                psf_size[0],
                                                psf_size[1],
                                                rotatedImage, 
                                                background );
                    else
                        et_convolve2D(          rotatedImage,
                                                psfImage,
                                                rotatedImage, 
                                                background );
                    }

		// Integrate along lines //
                et_line_integral_attenuated(    rotatedImage, 
                                                rotatedAttenuationImage, 
                                                sinoImage, 
                                                cam, 
                                                background); 
	}

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* sino_data = (float*) sinoImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<sinoImage->nvox; i++)
                {
                if (sino_data[i] < 0)
                    sino_data[i] = 0;
                }
            }

	/* Deallocate memory */
        return alloc_record_destroy(memory_record); 
}



//! Back-projection for Emission Imaging
/*!
  \param *sinogramImage the data to be back-projected in projection space. 
  \param *backprojectionImage the output backprojection. 
  \param *psfImage the depth-dependent point spread function, NULL for no point spread function. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_backproject(nifti_image *sinogramImage, nifti_image *backprojectionImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{

        int separable_psf = 0;
        int psf_size[3];

	/* Allocate the deformation Field image */
        alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);

	nifti_image *positionFieldImage = nifti_copy_nim_info(backprojectionImage);
	positionFieldImage->dim[0]=positionFieldImage->ndim=5;
	positionFieldImage->dim[1]=positionFieldImage->nx=backprojectionImage->nx;
	positionFieldImage->dim[2]=positionFieldImage->ny=backprojectionImage->ny;
	positionFieldImage->dim[3]=positionFieldImage->nz=backprojectionImage->nz;
	positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu=3;positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
	positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
	positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
	positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
	positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
	positionFieldImage->nbyper = sizeof(float); 
	positionFieldImage->data = (void *)calloc(positionFieldImage->nvox, positionFieldImage->nbyper); 
        if (positionFieldImage->data==NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,positionFieldImage,ALLOCTYPE_NIFTI);
	
	/* Allocate arrays */
        int dim[8];
	dim[0]    = 3;
	dim[1]    = backprojectionImage->dim[1];
	dim[2]    = backprojectionImage->dim[2];
	dim[3]    = backprojectionImage->dim[3];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *rotatedImage;
        rotatedImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        rotatedImage->data = (int *)malloc(backprojectionImage->nvox*sizeof(int));	
        if (rotatedImage->data == NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,rotatedImage,ALLOCTYPE_NIFTI); 

	nifti_image *temp_backprojectionImage;
        temp_backprojectionImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        temp_backprojectionImage->data = (int *)malloc(backprojectionImage->nvox*sizeof(int));	        
        if (temp_backprojectionImage->data == NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,temp_backprojectionImage,ALLOCTYPE_NIFTI); 

	nifti_image *rotatedAttenuationImage;
        if (attenuationImage != NULL)
            {
            rotatedAttenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedAttenuationImage->data = (int *)malloc(attenuationImage->nvox*sizeof(int));	
            if (rotatedAttenuationImage->data == NULL) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_alloccpu; }
            alloc_record_add(memory_record,rotatedAttenuationImage,ALLOCTYPE_NIFTI); 
            }

	/* Define centers of rotation */
	float center_x = ((float)(backprojectionImage->nx - 1)) / 2.0;
	float center_y = ((float)(backprojectionImage->ny - 1)) / 2.0;
	float center_z = ((float)(backprojectionImage->nz - 1)) / 2.0;

	/* Clear accumulator */
	et_clear_accumulator(backprojectionImage);	
	
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        if (affineTransformation == NULL) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_alloccpu; }
        alloc_record_add(memory_record,affineTransformation,ALLOCTYPE_GUEST); 

        /* Decide whether to use FFT or separate the convolution */
        float *psfSeparated=NULL;
        float psf_norm;
	if(psfImage!=NULL)
            {
            if (psfImage->nx <= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1) 
                {
                separable_psf=1;
                psf_size[0] = psfImage->dim[1];
                psf_size[1] = psfImage->dim[2];
                psf_size[2] = psfImage->dim[3];
                psfSeparated = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
                if (psfSeparated == NULL) {
                    alloc_record_destroy(memory_record); 
                    return niftyrec_error_alloccpu; }
                alloc_record_add(memory_record,psfSeparated,ALLOCTYPE_GUEST); 
                for (int n=0; n<psf_size[2];n++) 
                    {
                    psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                    for (int i=0;i<psf_size[0];i++) 
                        {
                        psfSeparated[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                        psfSeparated[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                        }
                    }
                }
            }       

        /* Backproject */
	for(int cam=0; cam<n_cameras; cam++){
                /* Rotate attenuation */

                if (attenuationImage != NULL)                
                    {
                    et_create_rotation_matrix(	affineTransformation,
						cameras[0*n_cameras+cam],
						cameras[1*n_cameras+cam],
						cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z,
						XYZ_ROTATION);
                    reg_affine_positionField(	affineTransformation,
						attenuationImage,
						positionFieldImage);
                    reg_resampleSourceImage(	attenuationImage,
						attenuationImage,
						rotatedAttenuationImage,
						positionFieldImage,
						NULL,
						1,
						background_attenuation );
                    }

		/* Line Backproject */
                if (attenuationImage != NULL)
                    {
                    et_line_backproject_attenuated(    sinogramImage,
                                                temp_backprojectionImage, 
                                                rotatedAttenuationImage, 
                                                cam );
                    }
                else
                    {
                    et_line_backproject(        sinogramImage,
                                                temp_backprojectionImage,
                                                cam );
                    }

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    {
                    if (separable_psf)
                        et_convolveSeparable2D( temp_backprojectionImage,
                                                psfSeparated,
                                                psf_size[0],
                                                psf_size[1],
                                                temp_backprojectionImage, 
                                                0.0f );
                    else
                        et_convolve2D(          temp_backprojectionImage,
                                                psfImage,
                                                temp_backprojectionImage, 
                                                0.0f );
                    }

		/* Rotate backprojection */
		et_create_rotation_matrix(	affineTransformation,
						-cameras[0*n_cameras+cam],
						-cameras[1*n_cameras+cam],
						-cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z,
						ZYX_ROTATION);
						
		reg_affine_positionField(	affineTransformation,
						backprojectionImage,
						positionFieldImage);

		reg_resampleSourceImage(	temp_backprojectionImage,
						temp_backprojectionImage,
						rotatedImage,
						positionFieldImage,
						NULL,
						1,
						background  );

		/* Accumulate */
		et_accumulate(			rotatedImage,
						backprojectionImage );
	}

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* accumulator_data = (float*) backprojectionImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<backprojectionImage->nvox; i++)
                {
                if (accumulator_data[i] < 0)
                   accumulator_data[i] = 0;
                }
            }
	/*Free*/
        return alloc_record_destroy(memory_record); 
}



//! Fisher Information Matrix of a grid of voxels, Emission Imaging
/*!
  \param from_projection whether the input image is a projection image (1) or activity image (0)
  \param *inputImage input image: projection image or activity image. 
  \param *gridImage grid image (same size as the activity), indexes from 1 to N_points at grid locations, 0 elsewhere. 
  \param *fisherImage output Fisher Information Matrix
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_fisher_grid(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *fisherpriorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
{
    int status = 0;
    int psf_size_x = psfImage->nx;
    int psf_size_y = psfImage->ny;
    int psf_size_semi_x = (psf_size_x-1)/2;
    int psf_size_semi_y = (psf_size_y-1)/2;
    float *fisher_matrix = (float *) fisherImage->data;
    float *fisher_matrix_prior=NULL;
    if (fisherpriorImage!=NULL)
        fisher_matrix_prior = (float *) fisherpriorImage->data;

    // 1) Project object and invert the sinogram elements
    int dim[8];
    dim[0] = 3;
    dim[1] = gridImage->dim[1];
    dim[2] = gridImage->dim[3];
    dim[3] = n_cameras;	  
    dim[4] = 1;	  
    dim[5] = 1;	  
    dim[6] = 1;	  
    dim[7] = 1;	  
    nifti_image *invsinogramImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
    float *invsinogram;
    if (from_projection==0)
        {
        invsinogram = (float*) malloc(dim[1]*dim[2]*dim[3]*sizeof(float));
        invsinogramImage->data = (float *)(invsinogram);
        status = et_project(inputImage, invsinogramImage, psfImage, attenuationImage, cameras, n_cameras, background, background_attenuation, 1);
        if (status)
            {
            fprintf(stderr,"'et_fisher_grid': error while calculating projection\n");
            return status;
            } 
        }
    else
        {
        invsinogramImage = inputImage;
        invsinogram = (float*)invsinogramImage->data;
        }

//    if (epsilon<=eps) epsilon=eps;
//    for (int i=0; i<invsinogramImage->nvox; i++)
//        invsinogram[i]=1/(invsinogram[i]+epsilon);

    // 2) Create (3xN) matrix of coordinates of the grid points; these will be rotated according to each camera position. 
    int n_grid_elements = fisherImage->dim[1];
    int *grid_coords  = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    for (int i=0; i<n_grid_elements*3; i++)
        grid_coords[i]=-1;
    int n=0;
    float *grid_data = (float*)gridImage->data;
    for (int x=0; x<gridImage->nx; x++) {
        for (int y=0; y<gridImage->ny; y++) {
            for (int z=0; z<gridImage->nz; z++) {
                if (grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] !=0) { //FIXME: make sure this works for non cubic images
                    n = grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] - 1;
                    if (grid_coords[n*3]!=-1)
                        return ET_ERROR_BADGRID;     // this happens if there are two elements of the grid with the same index
                    if (n > n_grid_elements) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is bigger than the numbe of elements 
                    if (n < 0) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is negative
                    grid_coords[n*3]=x;
                    grid_coords[n*3+1]=y;
                    grid_coords[n*3+2]=z;
                    }
                }
            }
        }
    // 3) For each camera position, update the FIM
    for (int i=0; i<n_grid_elements*n_grid_elements; i++)
        fisher_matrix[i]=0;
    if (fisher_matrix_prior!=NULL)
        {
        for (int i=0; i<n_grid_elements*n_grid_elements; i++)
            fisher_matrix_prior[i]=0;
        }

    float center_x = ((float)(gridImage->nx - 1)) / 2.0;
    float center_y = ((float)(gridImage->ny - 1)) / 2.0;
    float center_z = ((float)(gridImage->nz - 1)) / 2.0;
		
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    int *grid_coords_rotated = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    float position[3]; float position_rotated[3];
    
    for (int cam=0; cam<n_cameras; cam++)
        {
        float *invsino_data = (float *) (invsinogramImage->data) + cam*gridImage->nx*gridImage->nz ;
        // 3a) rotate the grid coordinates
        et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z, XYZ_ROTATION);
        for (int n=0; n<n_grid_elements; n++)
            {
            position[0]=grid_coords[3*n]; position[1]=grid_coords[3*n+1]; position[2]=grid_coords[3*n+2]; 
            reg_mat44_mul(affineTransformation, position, position_rotated);
            grid_coords_rotated[3*n] = position_rotated[0]; grid_coords_rotated[3*n+1] = position_rotated[1]; grid_coords_rotated[3*n+2] = position_rotated[2]; //implicit rounding (NN interpolation)
            }
        // 3b) for each pair, compute the FIM element (for the current camera, then it all sums up)
        int bbox0[2];
        int bbox1[2];
        int bbox_size[2];
        int Z;
        int psf_i_x, psf_i_y, psf_j_x, psf_j_y;
        float *PSF_i; float *PSF_j;
        float fisher_element;
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords_rotated[3*i];
            int i_y = grid_coords_rotated[3*i+1];
            int i_z = grid_coords_rotated[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                int j_x = grid_coords_rotated[3*j];
                int j_y = grid_coords_rotated[3*j+1];
                int j_z = grid_coords_rotated[3*j+2];
                // 3b1) compute bounding box
                int B_x = psf_size_x - abs(j_x-i_x);     // Bounding box size
                int B_y = psf_size_y - abs(j_y-i_y);
                int p_x = max(i_x,j_x)-psf_size_semi_x;  // start of bounding box in projection space
                int p_y = max(i_y,j_y)-psf_size_semi_y;
                // 3b2) compute coordinates of the PSFs
                if (B_x>0 && B_y>0 && i_z>=0 && j_z>=0 && i_z<gridImage->nz && j_z<gridImage->nz)  // check if the PSFs overlap and if the plane if within the field of view 
                    {
                    if ((j_x >= i_x) && (j_y > i_y))
                        {
                        psf_i_x = j_x - i_x;
                        psf_i_y = j_y - i_y;
                        psf_j_x = 0;
                        psf_j_y = 0;
                        }
                    else if ((j_x < i_x) && (j_y >= i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = j_y-i_y;
                        psf_j_x = i_x-j_x;
                        psf_j_y = 0;
                        }
                    else if ((j_x <= i_x) && (j_y < i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = 0;
                        psf_j_x = i_x-j_x;
                        psf_j_y = i_y-j_y;
                        }
                    else
                        {
                        psf_i_x = j_x-i_x;
                        psf_i_y = 0;
                        psf_j_x = 0;
                        psf_j_y = i_y-j_y;
                        }         
                    // 3b3) update pointers to the PSFs (function of z)
                    PSF_i = (float*) (psfImage->data) + i_z * psf_size_x*psf_size_y; 
                    PSF_j = (float*) (psfImage->data) + j_z * psf_size_x*psf_size_y; 
                    // 3b4) update the Fisher Information matrix
                    for (int x=0; x<B_x; x++)
                        {
                        for (int y=0; y<B_y; y++)
                            {
                            // check if the point is within the projection space
                            int proj_x = p_x+x;
                            int proj_y = p_y+y;
                            if (proj_x>=0 && proj_y>=0 && proj_x<invsinogramImage->nx && proj_y<invsinogramImage->ny)
                                {
                                fisher_element = fisher_matrix[i*n_grid_elements+j] + PSF_j[(psf_j_y+y)*psf_size_x+(psf_j_x+x)]*PSF_i[(psf_i_y+y)*psf_size_x+(psf_i_x+x)]*invsino_data[proj_y*invsinogramImage->nx + proj_x];
                                //invsinogram[cam*invsinogramImage->nx*invsinogramImage->nx+proj_y*invsinogramImage->nx + proj_x];
                                fisher_matrix[i*n_grid_elements+j] = fisher_element;
                                }
                            }
                        }
                    }
                }
            }
        }

    // 3c) Fisher Information of the prior
    if (fisher_matrix_prior!=NULL)
        {
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords[3*i];
            int i_y = grid_coords[3*i+1];
            int i_z = grid_coords[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                int j_x = grid_coords[3*j];
                int j_y = grid_coords[3*j+1];
                int j_z = grid_coords[3*j+2];
                if (abs(i_x-j_x)<=1 && abs(i_y-j_y)<=1 && abs(i_z-j_z)<=1)
                    fisher_matrix_prior[i*n_grid_elements+j] = 1;
//                float dist = sqrt((i_x-j_x)^2 + (i_y-j_y)^2 + (i_z-j_z)^2);
//                fisher_matrix_prior[i*n_grid_elements+j] = dist;
                }
            }
        }

    // 4) Fill matrix (the other half)
    for (int i=0; i<n_grid_elements; i++)
        for (int j=i+1; j<n_grid_elements; j++)
            fisher_matrix[j*n_grid_elements+i] = fisher_matrix[i*n_grid_elements+j];

    if (fisher_matrix_prior!=NULL)
        for (int i=0; i<n_grid_elements; i++)
            for (int j=i+1; j<n_grid_elements; j++)
                fisher_matrix_prior[j*n_grid_elements+i] = fisher_matrix_prior[i*n_grid_elements+j];


    // Free
    free(grid_coords);
    free(affineTransformation);
    free(grid_coords_rotated);
    if (from_projection==0)
        {
        free(invsinogram);
        (void)nifti_free_extensions( invsinogramImage) ;
        free(invsinogramImage) ;
        }
    return status;
}



//! Gradient of the attenuation map. Use this in order to estimate the attenuation map from the emission data. 
/*!
  \param *gradientImage outpup gradient in voxel space. 
  \param *sinoImage input sinogram.  
  \param *activityImage input activity (estimate). 
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_gradient_attenuation(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values) 
{
    return 1;
}



//! Convolve a stack of 2D images. 
/*!
  \param *inImage input stack of images. 
  \param *outImage output convolved stack of images. 
  \param *psfImage convolution kernel. 
*/
int et_convolve(nifti_image *inImage, nifti_image *outImage, nifti_image *kernelImage)
{
    int status = 1;
    return status;
}


int et_project_partial(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *partialsumImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{
return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   GPU   ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _USE_CUDA

//! Affine transformation of nifti_image on GPU
/*!
  \param *sourceImage the source image to be transformed. 
  \param *resultImage the transformed image. 
  \param *affineTransformation the [4x4] transformed matrix. 
  \param background the background value when resampling the transformed image. 
*/
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
	cudaMemcpy(mask_d, mask_h, resultImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
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


//! Rotate a nifti_image in 3D on the GPU
/*!
  \param *sourceImage the source image to be transformed. 
  \param *resultImage the rotated image. 
  \param theta_x the rotation angle around x axis in radians. 
  \param theta_y the rotation angle around y axis in radians. 
  \param theta_z the rotation angle around z axis in radians. 
  \param center_x the center of rotation along x.  
  \param center_y the center of rotation along y. 
  \param center_z the center of rotation along z. 
  \param background the background value when resampling the transformed image. 
*/
int et_rotate_gpu(nifti_image *sourceImage, nifti_image *resultImage, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z, float background)
{
	int status = 1;
	
	//Create transformation matrix
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
	et_create_rotation_matrix(affineTransformation, theta_x, theta_y, theta_z, center_x, center_y, center_z, XYZ_ROTATION);
	
	//Apply affine transformation
	status = et_affine_gpu(sourceImage, resultImage, affineTransformation, background);

	//Free
	free(affineTransformation);

	return status;
}


//! Projection for Emission Imaging, on GPU
/*!
  \param *activityImage the activity (or its estimate). NULL for attenuation and background activity only. 
  \param *sinoImage the photon counts in projection space. 
  \param *psfImage the depth-dependent point spread function, NULL for no PSF. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_project_gpu(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{
	/* initialise the cuda arrays */
	cudaArray *activityArray_d=NULL;               //stores input activity, makes use of fetch unit
        cudaArray *attenuationArray_d=NULL;            //stores input attenuation coefficients, makes use of fetch unit
	float     *sinoArray_d=NULL;                   //stores sinogram (output)
	float     *rotatedArray_d=NULL;                //stores activity aligned with current camera
        float     *rotatedAttenuationArray_d=NULL;     //stores attenuation coefficients aligned with current camera
	float4    *positionFieldImageArray_d=NULL;     //stores position field for rotation of activity and attenuation
	int       *mask_d=NULL;                        //binary mask that defines active voxels (typically all active)
        float     *psfArray_d=NULL;                    //stores point spread function
        float     *psfSeparatedArray_d=NULL;           //stores point spread function
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

        /* Check consistency of input */
        nifti_image *referenceImage;              // this image holds information about image size and voxel size (activity or attenuation might not be defined (NULL pointers) ) 
        if (activityImage==NULL && attenuationImage==NULL)
            {
            fprintf(stderr, "et_project_gpu: Error - define at least one between activityImage and attenuationImage. \n");
            return niftyrec_error_parameters; 
            }
        else if (attenuationImage==NULL)
            referenceImage=activityImage;
        else
            referenceImage=attenuationImage; 
	
	/* Allocate arrays on the device and transfer data to the device */

        alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);         
        // Activity
        if (activityImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;} 
            alloc_record_add(memory_record,(void*)activityArray_d,ALLOCTYPE_CUDA_ARRAY);
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, activityImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;}  
            alloc_record_add(memory_record,(void*)rotatedArray_d,ALLOCTYPE_CUDA);  
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_transfergpu;} 
            }

        // Singoram 
	if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;} 
        alloc_record_add(memory_record,(void*)sinoArray_d,ALLOCTYPE_CUDA);
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, referenceImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;} 
        alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA);
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, referenceImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;} 
        alloc_record_add(memory_record,(void*)mask_d,ALLOCTYPE_CUDA);

	// Mask 
	int *mask_h=(int *)malloc(referenceImage->nvox*sizeof(int));
	for(int i=0; i<referenceImage->nvox; i++) mask_h[i]=i;
	cudaMemcpy(mask_d, mask_h, referenceImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
	free(mask_h);
	
	/* Define centers of rotation */
	float center_x = ((float)(referenceImage->nx - 1)) / 2.0;
	float center_y = ((float)(referenceImage->ny - 1)) / 2.0;
	float center_z = ((float)(referenceImage->nz - 1)) / 2.0;
		
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_GUEST);

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;} 
            alloc_record_add(memory_record,(void*)psfArray_d,ALLOCTYPE_CUDA);
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_transfergpu;} 
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = referenceImage->dim[1];
            image_size[1] = referenceImage->dim[2];
            image_size[2] = referenceImage->dim[3];

            if (psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            if(cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)) != cudaSuccess) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)psfSeparatedArray_d,ALLOCTYPE_CUDA);

            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice);
            free(psfSeparatedArray_h);
            }

        if (attenuationImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) {
                 alloc_record_destroy(memory_record); 
                 return niftyrec_error_transfergpu;} 
            alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY);
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) {
                 alloc_record_destroy(memory_record); 
                 return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)rotatedAttenuationArray_d,ALLOCTYPE_CUDA);
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) {
                 alloc_record_destroy(memory_record);
                 return niftyrec_error_transfergpu;}
            }

	for(unsigned int cam=0; cam<n_cameras; cam++){
                fprintf_verbose( "et_project: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
		// Apply affine //
		et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z, XYZ_ROTATION);
		reg_affine_positionField_gpu(	affineTransformation,
						referenceImage,
						&positionFieldImageArray_d);

		// Resample the activity image //
                if (activityImage != NULL)
		    reg_resampleSourceImage_gpu(activityImage,
						activityImage,
						&rotatedArray_d,
						&activityArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						activityImage->nvox,
						background);

                // Resample the attenuation map //
                if (attenuationImage != NULL)
		    reg_resampleSourceImage_gpu(attenuationImage,
						attenuationImage,
						&rotatedAttenuationArray_d,
						&attenuationArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						attenuationImage->nvox,
						background_attenuation);

                // Apply Depth Dependent Point Spread Function //
                if ((psfImage != NULL) && (activityImage!=NULL))
                    {
                    if (separable_psf)
                        {
                        int status = et_convolveSeparable2D_gpu(
                                                &rotatedArray_d, 
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                        if (status)
                            {
                            alloc_record_destroy(memory_record);
                            return niftyrec_error_kernel;
                            }
                        }
                    else
                        et_convolveFFT2D_gpu(   &rotatedArray_d, 
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                    }

		// Integrate along lines //
                if ((activityImage!=NULL) && (attenuationImage != NULL))
                    et_line_integral_attenuated_gpu(	rotatedArray_d,
						rotatedAttenuationArray_d, 
						sinoArray_d,
                                                NULL,
						cam,
						referenceImage,
                                                background);
                else if ((activityImage!=NULL) && (attenuationImage == NULL))
		    et_line_integral_attenuated_gpu(	rotatedArray_d,
                                                NULL,
						sinoArray_d,
                                                NULL,
						cam,
						referenceImage,
                                                background);
                else if ((activityImage==NULL) && (attenuationImage != NULL))
                    et_line_integral_attenuated_gpu(NULL,
						rotatedAttenuationArray_d, 
						sinoArray_d,
                                                NULL,
						cam,
						referenceImage,
                                                background);
                else 
                    {
                    alloc_record_destroy(memory_record);
                    return niftyrec_error_parameters;
                    }
	}

	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(sinoImage, &sinoArray_d)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;}

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* sino_data = (float*) sinoImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<sinoImage->nvox; i++)
                if (sino_data[i] < 0)
                    sino_data[i] = 0;
            }

	/*Free*/
        int status = alloc_record_destroy(memory_record); 
        return status; 
}



//! Back-projection for Emission Imaging, on GPU.
/*!
  \param *sinogramImage the data to be back-projected in projection space. 
  \param *backprojectionImage the output backprojection. 
  \param *psfImage the depth-dependent point spread function, NULL for no point spread function. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_backproject_gpu(nifti_image *sinoImage, nifti_image *backprojectionImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
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
        float     *psfSeparatedArray_d;
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

	/* Allocate the deformation Field image */
        alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);  
	
	/* Allocate arrays on the device */

        cudaChannelFormatDesc backprojectionArray_d_chdesc = cudaCreateChannelDesc<float>();
        cudaExtent backprojectionArray_d_extent;
        backprojectionArray_d_extent.width  = backprojectionImage->nx;
        backprojectionArray_d_extent.height = backprojectionImage->ny;
        backprojectionArray_d_extent.depth  = backprojectionImage->nz;
        if (cudaMalloc3DArray(&backprojectionArray_d, &backprojectionArray_d_chdesc, backprojectionArray_d_extent) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)backprojectionArray_d,ALLOCTYPE_CUDA_ARRAY);

        
        //allocate backprojection
        cudaPitchedPtr temp_backprojection_pitched; 
        cudaExtent temp_backprojection_extent = make_cudaExtent(sizeof(float)*backprojectionImage->nx,backprojectionImage->ny,backprojectionImage->nz); 	
        cudaError_t cuda_error = cudaMalloc3D(&temp_backprojection_pitched, temp_backprojection_extent); 	
        if(cuda_error != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        temp_backprojection_d = (float*) temp_backprojection_pitched.ptr;
        alloc_record_add(memory_record,(void*)temp_backprojection_d,ALLOCTYPE_CUDA);

        
        if (attenuationImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim) != cudaSuccess) {
                alloc_record_destroy(memory_record);
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY);
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim) != cudaSuccess) {
                alloc_record_destroy(memory_record);
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)rotatedAttenuationArray_d,ALLOCTYPE_CUDA);
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage) != cudaSuccess) {
                alloc_record_destroy(memory_record);
                return niftyrec_error_transfergpu;}
            }

	if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)sinoArray_d,ALLOCTYPE_CUDA);
	if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, backprojectionImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)rotatedArray_d,ALLOCTYPE_CUDA);	
	if(cudaCommon_allocateArrayToDevice<float>(&accumulatorArray_d, backprojectionImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)accumulatorArray_d,ALLOCTYPE_CUDA);
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, backprojectionImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA);
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, backprojectionImage->dim) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_allocgpu;}
        alloc_record_add(memory_record,(void*)mask_d,ALLOCTYPE_CUDA);

	/* Transfer data from the host to the device */
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sinoArray_d,sinoImage)) return 1;
	int *mask_h=(int *)malloc(backprojectionImage->nvox*sizeof(int)); 
        if (mask_h==NULL) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_alloccpu;}
	for(int i=0; i<backprojectionImage->nvox; i++) mask_h[i]=i;
	if (cudaMemcpy(mask_d, mask_h, backprojectionImage->nvox*sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_transfergpu;}
	free(mask_h);

	/* Define centers of rotation */
	float center_x = ((float)(backprojectionImage->nx - 1)) / 2.0;
	float center_y = ((float)(backprojectionImage->ny - 1)) / 2.0;
	float center_z = ((float)(backprojectionImage->nz - 1)) / 2.0;

	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        if (affineTransformation==NULL) {
            alloc_record_destroy(memory_record);
            return niftyrec_error_alloccpu;}
        alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_GUEST);

	/* Clear accumulator */
	et_clear_accumulator_gpu(&accumulatorArray_d,backprojectionImage );

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim) != cudaSuccess) {
                alloc_record_destroy(memory_record);
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)psfArray_d,ALLOCTYPE_CUDA); 
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage) != cudaSuccess) {
                alloc_record_destroy(memory_record);
                return niftyrec_error_transfergpu;}
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = backprojectionImage->dim[1];
            image_size[1] = backprojectionImage->dim[2];
            image_size[2] = backprojectionImage->dim[3];

            if (psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            if (cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)) != cudaSuccess) {
                alloc_record_destroy(memory_record);
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)psfSeparatedArray_d,ALLOCTYPE_CUDA); 
            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            if (cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) {
                alloc_record_destroy(memory_record);
                free(psfSeparatedArray_h);
                return niftyrec_error_transfergpu;}
            free(psfSeparatedArray_h);
            }

	for(int cam=0; cam<n_cameras; cam++){
                // Rotate attenuation //
                if (attenuationImage != NULL)
                    {
                    et_create_rotation_matrix(	affineTransformation,
						cameras[0*n_cameras+cam],
						cameras[1*n_cameras+cam],
						cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z,
						XYZ_ROTATION);
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
                    }
		// Line Backproject //
                if (attenuationImage != NULL)
                    {
                    et_line_backproject_attenuated_gpu(	&sinoArray_d,
						&temp_backprojection_d,
						&rotatedAttenuationArray_d,
                                                cam,
						backprojectionImage);
                    }
                else
                    {
                    et_line_backproject_gpu(	&sinoArray_d,
						&temp_backprojection_d,
						cam,
						backprojectionImage);
                    }

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    {
                    if (separable_psf)
                        et_convolveSeparable2D_gpu( &temp_backprojection_d,
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
                    else
                        et_convolveFFT2D_gpu(   &temp_backprojection_d,
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
                    }

                // Copy to texture bound memory (for rotation) //
                cudaError_t cuda_status;
                cudaExtent volumeSize = make_cudaExtent(backprojectionImage->nx, backprojectionImage->ny, backprojectionImage->nz);
                cudaMemcpy3DParms copyparms={0};
                copyparms.extent = volumeSize;
                copyparms.dstArray = backprojectionArray_d;
                copyparms.kind = cudaMemcpyDeviceToDevice; 
                copyparms.srcPtr = temp_backprojection_pitched;
                cuda_status = cudaMemcpy3D(&copyparms);

		        if (cuda_status != cudaSuccess)
		        {
				    fprintf(stderr, "Error copying to texture bound memory: %s\n",cudaGetErrorString(cuda_status));
				    return 1;
                }
		
		// Rotate backprojection //
		et_create_rotation_matrix(	affineTransformation,
						-cameras[0*n_cameras+cam],
						-cameras[1*n_cameras+cam],
						-cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z,
						ZYX_ROTATION);
		reg_affine_positionField_gpu(	affineTransformation,
						backprojectionImage,
						&positionFieldImageArray_d);
		reg_resampleSourceImage_gpu(	backprojectionImage,
						backprojectionImage,
						&rotatedArray_d,
						&backprojectionArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						backprojectionImage->nvox,
						background);

		// Accumulate //
		et_accumulate_gpu(&rotatedArray_d,
						&accumulatorArray_d,
						backprojectionImage );
	}

	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(backprojectionImage, &accumulatorArray_d)) return 1; 

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* accumulator_data = (float*) backprojectionImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<backprojectionImage->nvox; i++)
                if (accumulator_data[i] < 0)
                    accumulator_data[i] = 0;
            }

	/*Free*/
        int status = alloc_record_destroy(memory_record); 
        return status; 
}



//! Fisher Information Matrix of a grid of voxels, Emission Imaging, on GPU. 
/*!
  \param from_projection whether the input image is a projection image (1) or activity image (0)
  \param *inputImage input image: projection image or activity image. 
  \param *gridImage grid image (same size as the activity), indexes from 1 to N_points at grid locations, 0 elsewhere. 
  \param *fisherImage output Fisher Information Matrix
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_fisher_grid_gpu(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *fisherpriorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
{
    int status = 0;
    int psf_size_x = psfImage->nx;
    int psf_size_y = psfImage->ny;
    int psf_size_semi_x = (psf_size_x-1)/2;
    int psf_size_semi_y = (psf_size_y-1)/2;
    float *fisher_matrix = (float *) fisherImage->data;
    float *fisher_matrix_prior = NULL;
    if (fisherpriorImage!=NULL)
        fisher_matrix_prior = (float *) fisherpriorImage->data;

    // 1) Project object and invert the sinogram elements
    int dim[8];
    dim[0] = 3;
    dim[1] = gridImage->dim[1];
    dim[2] = gridImage->dim[3];
    dim[3] = n_cameras;	  
    dim[4] = 1;	  
    dim[5] = 1;	  
    dim[6] = 1;	  
    dim[7] = 1;	  
    nifti_image *invsinogramImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
    float *invsinogram;
    if (from_projection==0)
        {
        invsinogram = (float*) malloc(dim[1]*dim[2]*dim[3]*sizeof(float));
        invsinogramImage->data = (float *)(invsinogram);
        status = et_project_gpu(inputImage, invsinogramImage, psfImage, attenuationImage, cameras, n_cameras, background, background_attenuation, 1);
        if (status)
            {
            fprintf_verbose("'et_fisher_grid': error while calculating projection\n");
            return status;
            } 
        }
    else
        {
        invsinogramImage = inputImage;
        invsinogram = (float*)invsinogramImage->data;
        }

//    if (epsilon<=eps) epsilon=eps;
//    for (int i=0; i<invsinogramImage->nvox; i++)
//        invsinogram[i]=1/(invsinogram[i]+epsilon);

    // 2) Create (3xN) matrix of coordinates of the grid points; these will be rotated according to each camera position. 
    int n_grid_elements = fisherImage->dim[1];
    int *grid_coords  = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    for (int i=0; i<n_grid_elements*3; i++)
        grid_coords[i]=-1;
    int n=0;
    float *grid_data = (float*)gridImage->data;
    for (int x=0; x<gridImage->nx; x++) {
        for (int y=0; y<gridImage->ny; y++) {
            for (int z=0; z<gridImage->nz; z++) {
                if (grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] !=0) { //FIXME: make sure this works for non cubic images
                    n = grid_data[z*gridImage->ny*gridImage->nx+y*gridImage->nx+x] - 1;
                    if (grid_coords[n*3]!=-1)
                        return ET_ERROR_BADGRID;     // this happens if there are two elements of the grid with the same index
                    if (n > n_grid_elements) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is bigger than the numbe of elements 
                    if (n < 0) 
                        return ET_ERROR_BADGRID;     // this happens if at least one of the elements of the grid is negative
                    grid_coords[n*3]=x;
                    grid_coords[n*3+1]=y;
                    grid_coords[n*3+2]=z;
                    }
                }
            }
        }
    // 3) For each camera position, update the FIM
    for (int i=0; i<n_grid_elements*n_grid_elements; i++)
        fisher_matrix[i]=0;
    if (fisher_matrix_prior!=NULL)
        for (int i=0; i<n_grid_elements*n_grid_elements; i++)
            fisher_matrix_prior[i]=0; 
   

    float center_x = ((float)(gridImage->nx - 1)) / 2.0;
    float center_y = ((float)(gridImage->ny - 1)) / 2.0;
    float center_z = ((float)(gridImage->nz - 1)) / 2.0;
		
    mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
    int *grid_coords_rotated = (int*) malloc(n_grid_elements*sizeof(int)*3); 
    float position[3]; float position_rotated[3];
    
    for (int cam=0; cam<n_cameras; cam++)
        {
        float *invsino_data = (float *) (invsinogramImage->data) + cam*gridImage->nx*gridImage->nz ;
        // 3a) rotate the grid coordinates
        et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z, XYZ_ROTATION);
        for (int n=0; n<n_grid_elements; n++)
            {
            position[0]=grid_coords[3*n]; position[1]=grid_coords[3*n+1]; position[2]=grid_coords[3*n+2]; 
            reg_mat44_mul(affineTransformation, position, position_rotated);
            grid_coords_rotated[3*n] = position_rotated[0]; grid_coords_rotated[3*n+1] = position_rotated[1]; grid_coords_rotated[3*n+2] = position_rotated[2]; //implicit rounding (NN interpolation)
            }
        // 3b) for each pair, compute the FIM element (for the current camera, then it all sums up)
        int bbox0[2];
        int bbox1[2];
        int bbox_size[2];
        int Z;
        int psf_i_x, psf_i_y, psf_j_x, psf_j_y;
        float *PSF_i; float *PSF_j;
        float fisher_element;
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords_rotated[3*i];
            int i_y = grid_coords_rotated[3*i+1];
            int i_z = grid_coords_rotated[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                int j_x = grid_coords_rotated[3*j];
                int j_y = grid_coords_rotated[3*j+1];
                int j_z = grid_coords_rotated[3*j+2];
                // 3b1) compute bounding box
                int B_x = psf_size_x - abs(j_x-i_x);     // Bounding box size
                int B_y = psf_size_y - abs(j_y-i_y);
                int p_x = max(i_x,j_x)-psf_size_semi_x;  // start of bounding box in projection space
                int p_y = max(i_y,j_y)-psf_size_semi_y;
                // 3b2) compute coordinates of the PSFs
                if (B_x>0 && B_y>0 && i_z>=0 && j_z>=0 && i_z<gridImage->nz && j_z<gridImage->nz)  // check if the PSFs overlap and if the plane if within the field of view 
                    {
                    if ((j_x >= i_x) && (j_y > i_y))
                        {
                        psf_i_x = j_x - i_x;
                        psf_i_y = j_y - i_y;
                        psf_j_x = 0;
                        psf_j_y = 0;
                        }
                    else if ((j_x < i_x) && (j_y >= i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = j_y-i_y;
                        psf_j_x = i_x-j_x;
                        psf_j_y = 0;
                        }
                    else if ((j_x <= i_x) && (j_y < i_y))
                        {
                        psf_i_x = 0;
                        psf_i_y = 0;
                        psf_j_x = i_x-j_x;
                        psf_j_y = i_y-j_y;
                        }
                    else
                        {
                        psf_i_x = j_x-i_x;
                        psf_i_y = 0;
                        psf_j_x = 0;
                        psf_j_y = i_y-j_y;
                        }         
                    // 3b3) update pointers to the PSFs (function of z)
                    PSF_i = (float*) (psfImage->data) + i_z * psf_size_x*psf_size_y; 
                    PSF_j = (float*) (psfImage->data) + j_z * psf_size_x*psf_size_y; 
                    // 3b4) update the Fisher Information matrix
                    for (int x=0; x<B_x; x++)
                        {
                        for (int y=0; y<B_y; y++)
                            {
                            // check if the point is within the projection space
                            int proj_x = p_x+x;
                            int proj_y = p_y+y;
                            if (proj_x>=0 && proj_y>=0 && proj_x<invsinogramImage->nx && proj_y<invsinogramImage->ny)
                                {
                                fisher_element = fisher_matrix[i*n_grid_elements+j] + PSF_j[(psf_j_y+y)*psf_size_x+(psf_j_x+x)]*PSF_i[(psf_i_y+y)*psf_size_x+(psf_i_x+x)]*invsino_data[proj_y*invsinogramImage->nx + proj_x];
                                //invsinogram[cam*invsinogramImage->nx*invsinogramImage->nx+proj_y*invsinogramImage->nx + proj_x];
                                fisher_matrix[i*n_grid_elements+j] = fisher_element;
                                }
                            }
                        }
                    }
                }
            }
        }

    // 3c) Fisher Information of the prior
    if (fisher_matrix_prior!=NULL)
        {
        for (int i=0; i<n_grid_elements; i++)
            {
            int i_x = grid_coords[3*i];
            int i_y = grid_coords[3*i+1];
            int i_z = grid_coords[3*i+2];
            for (int j=i; j<n_grid_elements; j++)
                {
                if (i!=j)
                    {
                    int j_x = grid_coords[3*j];
                    int j_y = grid_coords[3*j+1];
                    int j_z = grid_coords[3*j+2];
                    if (abs(i_x-j_x)<=3 && abs(i_y-j_y)<=3 && abs(i_z-j_z)<=3)
                        fisher_matrix_prior[i*n_grid_elements+j] = 1;
                    }
//                float dist = sqrt((i_x-j_x)^2 + (i_y-j_y)^2 + (i_z-j_z)^2);
//                fisher_matrix_prior[i*n_grid_elements+j] = dist;
                }
            }
        }

    // 4) Fill matrix (the other half)
    for (int i=0; i<n_grid_elements; i++)
        for (int j=i+1; j<n_grid_elements; j++)
            fisher_matrix[j*n_grid_elements+i] = fisher_matrix[i*n_grid_elements+j];

    if (fisher_matrix_prior!=NULL)
        for (int i=0; i<n_grid_elements; i++)
            for (int j=i+1; j<n_grid_elements; j++)
                fisher_matrix_prior[j*n_grid_elements+i] = fisher_matrix_prior[i*n_grid_elements+j];

    // Free
    free(grid_coords);
    free(affineTransformation);
    free(grid_coords_rotated);
    if (from_projection==0)
        {
        free(invsinogram);
        (void)nifti_free_extensions( invsinogramImage) ;
        free(invsinogramImage) ;
        }
    return status;
} 



//! Gradient of the attenuation map. Use this in order to estimate the attenuation map from the emission data. GPU version. 
/*!
  \param *gradientImage outpup gradient in voxel space. 
  \param *sinoImage input sinogram.  
  \param *activityImage input activity (estimate). 
  \param *psfImage Point Spread Function
  \param *attenuationImage attenuation map, save size as the activity. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_gradient_attenuation_gpu(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values) 
{
	/* initialise the cuda arrays */
	cudaArray *activityArray_d;               //stores input activity, makes use of fetch unit
	cudaArray *backprojectionArray_d;
	cudaArray *attenuationArray_d;
	float     *temp_backprojection_d;
	float     *sinoArray_d;
	float     *rotatedArray_d;
	float     *rotatedAttenuationArray_d;
        float     *attenuationPlaneArray_d;
	float     *gradientArray_d;
	float4    *positionFieldImageArray_d;
	int       *mask_d;
        float     *psfArray_d;
        float     *psfSeparatedArray_d;
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

	/* Allocate the deformation Field image */
	nifti_image *positionFieldImage = nifti_copy_nim_info(gradientImage);
	positionFieldImage->dim[0]=positionFieldImage->ndim=5;
	positionFieldImage->dim[1]=positionFieldImage->nx = gradientImage->nx;
	positionFieldImage->dim[2]=positionFieldImage->ny = gradientImage->ny;
	positionFieldImage->dim[3]=positionFieldImage->nz = gradientImage->nz;
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
        backprojectionArray_d_extent.width  = gradientImage->nx;
        backprojectionArray_d_extent.height = gradientImage->ny;
        backprojectionArray_d_extent.depth  = gradientImage->nz;
        cudaError_t cuda_status1 = cudaMalloc3DArray(&backprojectionArray_d, &backprojectionArray_d_chdesc, backprojectionArray_d_extent);

	if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) return 1;
        if(cudaCommon_allocateArrayToDevice<float>(&temp_backprojection_d, gradientImage->dim)) return 1;

        if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) return 1;
        if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) return 1;	
        if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) return 1;

	if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, gradientImage->dim)) return 1;	
	if(cudaCommon_allocateArrayToDevice<float>(&gradientArray_d, gradientImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, gradientImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, gradientImage->dim)) return 1;

	/* Transfer data from the host to the device */
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sinoArray_d,sinoImage)) return 1;
	int *mask_h=(int *)malloc(gradientImage->nvox*sizeof(int));
	for(int i=0; i<gradientImage->nvox; i++) mask_h[i]=i;
	cudaMemcpy(mask_d, mask_h, gradientImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
	free(mask_h);

	/* Define centers of rotation */
	float center_x = ((float)(gradientImage->nx - 1)) / 2.0;
	float center_y = ((float)(gradientImage->ny - 1)) / 2.0;
	float center_z = ((float)(gradientImage->nz - 1)) / 2.0;

	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));

	/* Clear gradient */
	et_clear_accumulator_gpu(		&gradientArray_d,
						gradientImage );
        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) return 1;
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = gradientImage->dim[1];
            image_size[1] = gradientImage->dim[2];
            image_size[2] = gradientImage->dim[3];

            if (psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice);
            free(psfSeparatedArray_h);
            }

	for(int cam=0; cam<n_cameras; cam++){
                // Rotate attenuation and activity //
                et_create_rotation_matrix(      affineTransformation,
						cameras[0*n_cameras+cam],
						cameras[1*n_cameras+cam],
						cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z,
						XYZ_ROTATION);
                reg_affine_positionField_gpu(   affineTransformation,
						attenuationImage,
						&positionFieldImageArray_d);
                reg_resampleSourceImage_gpu(    attenuationImage,
						attenuationImage,
						&rotatedAttenuationArray_d,
						&attenuationArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						attenuationImage->nvox,
						background_attenuation);
                reg_resampleSourceImage_gpu(    activityImage,
						activityImage,
						&rotatedArray_d,
						&activityArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						activityImage->nvox,
						background);

		// Line Backproject, compute gradient of the likelihood with respect of the attenuation coefficients //
                et_attenuation_gradient_gpu(    &rotatedArray_d, 
                                                &sinoArray_d, 
						&temp_backprojection_d, 
						&rotatedAttenuationArray_d, 
                                                cam, 
						gradientImage); 


                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
                    {
                    if (separable_psf)
                        et_convolveSeparable2D_gpu( &temp_backprojection_d,
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
                    else
                        et_convolveFFT2D_gpu(   &temp_backprojection_d,
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &temp_backprojection_d);
                    }

                // Copy to texture bound memory (for rotation) //
                cudaError_t cuda_status;
                cudaMemcpy3DParms p = {0};
              
                p.srcPtr.ptr        = temp_backprojection_d;
                p.srcPtr.pitch      = 0;
                p.srcPtr.xsize      = gradientImage->nx;
                p.srcPtr.ysize      = gradientImage->ny;
                p.dstArray          = backprojectionArray_d;
                p.extent.width      = gradientImage->nx;
                p.extent.height     = gradientImage->ny;
                p.extent.depth      = gradientImage->nz;
                p.kind              = cudaMemcpyDeviceToDevice;
                cuda_status         = cudaMemcpy3D(&p);
		        if (cuda_status != cudaSuccess)
				    {
			        fprintf(stderr, "Error copying to texture bound memory: %s\n",cudaGetErrorString(cuda_status));
			        return 1;
		            }
		// Rotate backprojection //
		et_create_rotation_matrix(	affineTransformation,
						-cameras[0*n_cameras+cam],
						-cameras[1*n_cameras+cam],
						-cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z,
						ZYX_ROTATION);
		reg_affine_positionField_gpu(	affineTransformation,
						gradientImage,
						&positionFieldImageArray_d);
		reg_resampleSourceImage_gpu(	gradientImage,
						gradientImage,
						&rotatedArray_d,
						&backprojectionArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						gradientImage->nvox,
						background);

		// Accumulate //
		et_accumulate_gpu(		&rotatedArray_d,
						&gradientArray_d,
						gradientImage );
	}

	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(gradientImage, &gradientArray_d)) return 1; 

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* gradient_data = (float*) gradientImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<gradientImage->nvox; i++)
                if (gradient_data[i] < 0)
                    gradient_data[i] = 0;
            }

	/*Free*/
	cudaCommon_free(&activityArray_d);
        cudaCommon_free(&backprojectionArray_d);

        cudaCommon_free(&attenuationArray_d);
        cudaCommon_free((void **)&rotatedAttenuationArray_d);

	cudaCommon_free((void **)&rotatedArray_d);
	cudaCommon_free((void **)&sinoArray_d);
	cudaCommon_free((void **)&gradientArray_d);
	cudaCommon_free((void **)&mask_d);
	cudaCommon_free((void **)&positionFieldImageArray_d);
	cudaCommon_free((void **)&temp_backprojection_d);
        if (psfImage != NULL)
            {
            cudaCommon_free((void **)&psfArray_d);
            if (separable_psf)
                cudaCommon_free((void **)&psfSeparatedArray_d);
            }
	nifti_image_free(positionFieldImage);
	free(affineTransformation);

	return 0;
}



//! Convolve a stack of 2D images on the GPU. 
/*!
  \param *inImage input stack of images. 
  \param *outImage output convolved stack of images. 
  \param *psfImage convolution kernel. 
*/
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



//! List NVIDIA CUDA compatible GPU's installed in the system. 
/*!
  \param *device_count_out output, number of installed GPUs. 
  \param *devices outoput, GPU devices compute capability and ID's. See et_list_gpus_mex.  
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


//! Set GPU to be used by NiftyRec. 
/*!
  \param id GPU ID. 
*/
int et_set_gpu(int id)
{
    int status = 1;
    struct cudaDeviceProp deviceProp;

    cudaSetDevice( id );
    cudaGetDeviceProperties(&deviceProp, id );
    if (deviceProp.major < 1)
        {
        printf("ERROR - The specified graphical card does not exist.\n");
        status = 1;
	}
    else
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
	cudaMalloc((void **)&matrix_A_d,  matrix_A_Image->nvox*sizeof(float));
	cudaMalloc((void **)&matrix_B_d,  matrix_B_Image->nvox*sizeof(float));
	cudaMalloc((void **)&joint_histogram_d,  joint_histogram_Image->nvox*sizeof(int));
        cudaMemset((void*)joint_histogram_d,0,joint_histogram_Image->nvox*sizeof(int));
        
	// Transfer data from the host to the device 
	cudaMemcpy(matrix_A_d, matrix_A_Image->data, matrix_A_Image->nvox*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(matrix_B_d, matrix_B_Image->data, matrix_B_Image->nvox*sizeof(float), cudaMemcpyHostToDevice);

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
	cudaMemcpy(joint_histogram_Image->data, joint_histogram_d, joint_histogram_Image->nvox*sizeof(int), cudaMemcpyDeviceToHost);
	
	// Free arrays on device 
	cudaCommon_free((void **)&matrix_A_d);
	cudaCommon_free((void **)&matrix_B_d);
	cudaCommon_free((void **)&joint_histogram_d);
	
	return 0;
}
*/



//! Projection for Emission Imaging, on GPU
/*!
  \param *activityImage the activity (or its estimate). NULL for attenuation and background activity only. 
  \param *sinoImage the photon counts in projection space. 
  \param *psfImage the depth-dependent point spread function, NULL for no PSF. 
  \param *attenuationImage the attenuation map, NULL for no attenuation. 
  \param *cameras [n_camerasx3] array of camera orientations in radians. 
  \param n_cameras number of projections (camera positions). 
  \param background the activity background (used when activity is rotated and resampled). 
  \param background_attenuation the attenuation background (used when the attenuation map is rotated and resampled). 
*/
int et_project_partial_gpu(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *partialsumImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation, int truncate_negative_values)
{
	/* initialise the cuda arrays */
	cudaArray *activityArray_d=NULL;               //stores input activity, makes use of fetch unit
        cudaArray *attenuationArray_d=NULL;            //stores input attenuation coefficients, makes use of fetch unit
	float     *sinoArray_d=NULL;                   //stores sinogram (output)
	float     *partialsumArray_d=NULL;             //stores sinogram (output)
	float     *rotatedArray_d=NULL;                //stores activity aligned with current camera
        float     *rotatedAttenuationArray_d=NULL;     //stores attenuation coefficients aligned with current camera
	float4    *positionFieldImageArray_d=NULL;     //stores position field for rotation of activity and attenuation
	int       *mask_d=NULL;                        //binary mask that defines active voxels (typically all active)
        float     *psfArray_d=NULL;                    //stores point spread function
        float     *psfSeparatedArray_d=NULL;           //stores point spread function
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

        /* Check consistency of input */
        nifti_image *referenceImage;              // this image holds information about image size and voxel size (activity or attenuation might not be defined (NULL pointers) ) 
        if (activityImage==NULL && attenuationImage==NULL)
            {
            fprintf(stderr, "et_project_partial_gpu: Error - define at least one between activityImage and attenuationImage. \n");
            return niftyrec_error_parameters; 
            }
        else if (attenuationImage==NULL)
            referenceImage=activityImage;
        else
            referenceImage=attenuationImage; 
	
	/* Allocate arrays on the device and transfer data to the device */

        alloc_record *memory_record = alloc_record_create(RECORD_MAXELEMENTS);         
        // Activity
        if (activityImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;} 
            alloc_record_add(memory_record,(void*)activityArray_d,ALLOCTYPE_CUDA_ARRAY);
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, activityImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;}  
            alloc_record_add(memory_record,(void*)rotatedArray_d,ALLOCTYPE_CUDA);  
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_transfergpu;} 
            }
        // Partial sum
        if(cudaCommon_allocateArrayToDevice<float>(&partialsumArray_d, referenceImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;} 
        alloc_record_add(memory_record,(void*)partialsumArray_d,ALLOCTYPE_CUDA);

        // Singoram 
	if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;} 
        alloc_record_add(memory_record,(void*)sinoArray_d,ALLOCTYPE_CUDA);
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, referenceImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_allocgpu;} 
        alloc_record_add(memory_record,(void*)positionFieldImageArray_d,ALLOCTYPE_CUDA);
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, referenceImage->dim)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;} 
        alloc_record_add(memory_record,(void*)mask_d,ALLOCTYPE_CUDA);

	// Mask 
	int *mask_h=(int *)malloc(referenceImage->nvox*sizeof(int));
	for(int i=0; i<referenceImage->nvox; i++) mask_h[i]=i;
	cudaMemcpy(mask_d, mask_h, referenceImage->nvox*sizeof(int), cudaMemcpyHostToDevice);
	free(mask_h);
	
	/* Define centers of rotation */
	float center_x = ((float)(referenceImage->nx - 1)) / 2.0;
	float center_y = ((float)(referenceImage->ny - 1)) / 2.0;
	float center_z = ((float)(referenceImage->nz - 1)) / 2.0;
		
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));
        alloc_record_add(memory_record,(void*)affineTransformation,ALLOCTYPE_GUEST);

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;} 
            alloc_record_add(memory_record,(void*)psfArray_d,ALLOCTYPE_CUDA);
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_transfergpu;} 
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = referenceImage->dim[1];
            image_size[1] = referenceImage->dim[2];
            image_size[2] = referenceImage->dim[3];

            if (psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            if(cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)) != cudaSuccess) {
                alloc_record_destroy(memory_record); 
                return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)psfSeparatedArray_d,ALLOCTYPE_CUDA);

            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice);
            free(psfSeparatedArray_h);
            }

        if (attenuationImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) {
                 alloc_record_destroy(memory_record); 
                 return niftyrec_error_transfergpu;} 
            alloc_record_add(memory_record,(void*)attenuationArray_d,ALLOCTYPE_CUDA_ARRAY);
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) {
                 alloc_record_destroy(memory_record); 
                 return niftyrec_error_allocgpu;}
            alloc_record_add(memory_record,(void*)rotatedAttenuationArray_d,ALLOCTYPE_CUDA);
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) {
                 alloc_record_destroy(memory_record);
                 return niftyrec_error_transfergpu;}
            }

	for(unsigned int cam=0; cam<n_cameras; cam++){
                fprintf_verbose( "et_project: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
		// Apply affine //
		et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z, XYZ_ROTATION);
		reg_affine_positionField_gpu(	affineTransformation,
						referenceImage,
						&positionFieldImageArray_d);

		// Resample the activity image //
                if (activityImage != NULL)
		    reg_resampleSourceImage_gpu(activityImage,
						activityImage,
						&rotatedArray_d,
						&activityArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						activityImage->nvox,
						background);

                // Resample the attenuation map //
                if (attenuationImage != NULL)
		    reg_resampleSourceImage_gpu(attenuationImage,
						attenuationImage,
						&rotatedAttenuationArray_d,
						&attenuationArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						attenuationImage->nvox,
						background_attenuation);

                // Apply Depth Dependent Point Spread Function //
                if ((psfImage != NULL) && (activityImage!=NULL))
                    {
                    if (separable_psf)
                        {
                        int status = et_convolveSeparable2D_gpu(
                                                &rotatedArray_d, 
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                        if (status)
                            {
                            alloc_record_destroy(memory_record);
                            return niftyrec_error_kernel;
                            }
                        }
                    else
                        et_convolveFFT2D_gpu(   &rotatedArray_d, 
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                    }

		// Integrate along lines //
                if ((activityImage!=NULL) && (attenuationImage != NULL))
                    et_line_integral_attenuated_gpu(	rotatedArray_d,
						rotatedAttenuationArray_d, 
						sinoArray_d,
                                                partialsumArray_d,
						cam,
						referenceImage,
                                                background);
                else if ((activityImage!=NULL) && (attenuationImage == NULL))
		    et_line_integral_attenuated_gpu(	rotatedArray_d,
                                                NULL,
						sinoArray_d,
                                                partialsumArray_d,
						cam,
						referenceImage,
                                                background);
                else if ((activityImage==NULL) && (attenuationImage != NULL))
                    et_line_integral_attenuated_gpu(NULL,
						rotatedAttenuationArray_d, 
						sinoArray_d,
                                                partialsumArray_d,
						cam,
						referenceImage,
                                                background);
                else 
                    {
                    alloc_record_destroy(memory_record);
                    return niftyrec_error_parameters;
                    }
                if (cudaMemcpy(((float*) partialsumImage->data)+referenceImage->nx*referenceImage->ny*referenceImage->nz*cam, partialsumArray_d, referenceImage->nx*referenceImage->ny*referenceImage->nz*sizeof(float), cudaMemcpyDeviceToHost) != cudaSuccess) {
                    alloc_record_destroy(memory_record); 
                    return niftyrec_error_transfergpu;}
	}

	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(sinoImage, &sinoArray_d)) {
            alloc_record_destroy(memory_record); 
            return niftyrec_error_transfergpu;}

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* sino_data = (float*) sinoImage->data;
        if (truncate_negative_values)
            {
            for (int i=0; i<sinoImage->nvox; i++)
                if (sino_data[i] < 0)
                    sino_data[i] = 0;
            }

	/*Free*/
        return alloc_record_destroy(memory_record); 
}



//! Reset GPU
int et_reset_gpu()
{
//    #ifdef _CUDA_4
    cudaDeviceReset();
    return niftyrec_success;
//    #else
//    return niftyrec_error_nogpubuilt;  
//    #endif
}

#endif





