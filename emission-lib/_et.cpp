/*
 *  _et.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  Centre for Medical Image Computing (CMIC)
 *  University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et.h"
#include "_et_common.h"
#include<stdio.h>
#ifdef _OPENMP
#include "omp.h"
#endif

#define max(a,b)	(((a) > (b)) ? (a) : (b))
#define min(a,b)	(((a) < (b)) ? (a) : (b))

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////   CPU   ///////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
#ifdef _OPENMP
    omp_set_dynamic( 0 );
    omp_set_num_threads( 8 );
#endif

        int separable_psf = 0;
        int psf_size[3];

        /* Check consistency of input */
        //...

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
        rotatedImage->data = (float *)malloc(activityImage->nvox*sizeof(float));	

	nifti_image *rotatedAttenuationImage;
        if (attenuationImage != NULL)
            {
            rotatedAttenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedAttenuationImage->data = (float *)malloc(activityImage->nvox*sizeof(float));	
            }

	/* Define centers of rotation */
	float center_x = ((float)(activityImage->nx - 1)) / 2.0;
	float center_y = ((float)(activityImage->ny - 1)) / 2.0;
	float center_z = ((float)(activityImage->nz - 1)) / 2.0;
		
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));

        /* Decide whether to use FFT or separate the convolution */
        float *psfSeparated=NULL;
        float psf_norm;
	if(psfImage!=NULL)
            {
            if(1) //if (psfImage->nx <= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1) 
                {
                separable_psf=1;
                psf_size[0] = psfImage->dim[1];
                psf_size[1] = psfImage->dim[2];
                psf_size[2] = psfImage->dim[3];
                psfSeparated = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
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

                // Resample the attenuation map //
                if (attenuationImage != NULL)
                    {
                    reg_resampleSourceImage<float>(    attenuationImage,
						attenuationImage,
						rotatedAttenuationImage,
						positionFieldImage,
						NULL,
						1,
						background_attenuation );	                    
                    }

                // Apply Depth Dependent Point Spread Function //
                if (psfImage != NULL)
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
                if (attenuationImage != NULL)
                    {
                    et_line_integral_attenuated(rotatedImage, 
                                                rotatedAttenuationImage, 
                                                sinoImage, 
                                                cam );
                    }
                else
                    {
                    et_line_integral(           rotatedImage,
                                                sinoImage,
                                                cam );
                    }

	}

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* sino_data = (float*) sinoImage->data;
//#pragma omp parallel for 
        for (int i=0; i<sinoImage->nvox; i++)
            {
            if (sino_data[i] < 0)
                sino_data[i] = 0;
            }

	/* Deallocate memory */
	nifti_image_free(rotatedImage);
        if (attenuationImage != NULL)
            nifti_image_free(rotatedAttenuationImage);
        if (psfImage != NULL)
            if (separable_psf)
                free(psfSeparated);
	nifti_image_free(positionFieldImage);
	free(affineTransformation);

	return 0;
}



int et_backproject(nifti_image *sinogramImage, nifti_image *accumulatorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
{
#ifdef _OPENMP
    omp_set_num_threads( 8 );
#endif
        int separable_psf = 0;
        int psf_size[3];

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

	nifti_image *rotatedAttenuationImage;
        if (attenuationImage != NULL)
            {
            rotatedAttenuationImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
            rotatedAttenuationImage->data = (int *)malloc(attenuationImage->nvox*sizeof(int));	
            }

	/* Define centers of rotation */
	float center_x = ((float)(accumulatorImage->nx - 1)) / 2.0;
	float center_y = ((float)(accumulatorImage->ny - 1)) / 2.0;
	float center_z = ((float)(accumulatorImage->nz - 1)) / 2.0;

	/* Clear accumulator */
	et_clear_accumulator(accumulatorImage);	
	
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));

        /* Decide whether to use FFT or separate the convolution */
        float *psfSeparated=NULL;
        float psf_norm;
	if(psfImage!=NULL)
            {
            if(1) //if (psfImage->nx <= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1) 
                {
                separable_psf=1;
                psf_size[0] = psfImage->dim[1];
                psf_size[1] = psfImage->dim[2];
                psf_size[2] = psfImage->dim[3];
                psfSeparated = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
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
						center_z);
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
						background  );

		/* Accumulate */
		et_accumulate(			rotatedImage,
						accumulatorImage );
	}

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* accumulator_data = (float*) accumulatorImage->data;
//#pragma omp parallel for 
        for (int i=0; i<accumulatorImage->nvox; i++)
            {
            if (accumulator_data[i] < 0)
               accumulator_data[i] = 0;
            }
	/*Free*/
	nifti_image_free(rotatedImage);
        if (attenuationImage != NULL)        
            nifti_image_free(rotatedAttenuationImage);
        if (psfImage != NULL)
            if (separable_psf)
                free(psfSeparated);
	nifti_image_free(temp_backprojectionImage);
	nifti_image_free(positionFieldImage);
	free(affineTransformation);

	return 0;
}



int et_fisher_grid(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *fisherpriorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float epsilon, float background, float background_attenuation)
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

    if (epsilon<=eps) epsilon=eps;
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
        status = et_project(inputImage, invsinogramImage, psfImage, attenuationImage, cameras, n_cameras, background, background_attenuation);
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
        et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z);
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



int et_project_backproject(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_alpha, float *cameras_beta, float *cameras_gamma)
{
    return 1;
}


int et_gradient_attenuation(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation) 
{
    return 1;
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



int et_project_gpu(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation)
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
        float     *psfSeparatedArray_d;           //stores point spread function
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

        /* Check consistency of input */
        //..
	
	/* Allocate arrays on the device */
	if(cudaCommon_allocateArrayToDevice<float>(&activityArray_d, activityImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&sinoArray_d, sinoImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float>(&rotatedArray_d, activityImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<float4>(&positionFieldImageArray_d, activityImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice<int>(&mask_d, activityImage->dim)) return 1;

	/* Transfer data from the host to the device */
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&activityArray_d,activityImage)) return 1;
	int *mask_h=(int *)malloc(activityImage->nvox*sizeof(int));
	for(int i=0; i<activityImage->nvox; i++) mask_h[i]=i;
	CUDA_SAFE_CALL(cudaMemcpy(mask_d, mask_h, activityImage->nvox*sizeof(int), cudaMemcpyHostToDevice));
	free(mask_h);
	
	/* Define centers of rotation */
	float center_x = ((float)(activityImage->nx - 1)) / 2.0;
	float center_y = ((float)(activityImage->ny - 1)) / 2.0;
	float center_z = ((float)(activityImage->nz - 1)) / 2.0;
		
	/* Alloc transformation matrix */
	mat44 *affineTransformation = (mat44 *)calloc(1,sizeof(mat44));

        /* Allocate and initialize kernel for DDPSF */
        if (psfImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&psfArray_d, psfImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&psfArray_d, psfImage)) return 1;
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = activityImage->dim[1];
            image_size[1] = activityImage->dim[2];
            image_size[2] = activityImage->dim[3];

            if (psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            CUDA_SAFE_CALL(cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)));
            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            CUDA_SAFE_CALL(cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice));
            free(psfSeparatedArray_h);
            }

        if (attenuationImage != NULL)
            {
            if(cudaCommon_allocateArrayToDevice<float>(&attenuationArray_d, attenuationImage->dim)) return 1;
            if(cudaCommon_allocateArrayToDevice<float>(&rotatedAttenuationArray_d, attenuationImage->dim)) return 1;
            if(cudaCommon_transferNiftiToArrayOnDevice<float>(&attenuationArray_d,attenuationImage)) return 1;
            }

	for(unsigned int cam=0; cam<n_cameras; cam++){
                fprintf_verbose( "et_project: Rotation: %f  %f  %f  \n",cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam]);
		// Apply affine //
		et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z);
		reg_affine_positionField_gpu(	affineTransformation,
						activityImage,
						&positionFieldImageArray_d);

		// Resample the source image //
		reg_resampleSourceImage_gpu(	activityImage,
						activityImage,
						&rotatedArray_d,
						&activityArray_d,
						&positionFieldImageArray_d,
						&mask_d,
						activityImage->nvox,
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
                    {
                    if (separable_psf)
                        et_convolveSeparable2D_gpu(&rotatedArray_d, 
                                                image_size,
                                                &psfSeparatedArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                    else
                        et_convolveFFT2D_gpu(   &rotatedArray_d, 
                                                image_size,
                                                &psfArray_d,
                                                psf_size,
                                                &rotatedArray_d);
                    }

		// Integrate along lines //
                if (attenuationImage != NULL)
                    {
                    et_line_integral_attenuated_gpu(	&rotatedArray_d,
						&rotatedAttenuationArray_d, 
						&sinoArray_d,
						cam,
						activityImage);
                    }
                else
                    {
		    et_line_integral_gpu(	&rotatedArray_d,
						&sinoArray_d,
						cam,
						activityImage);
                    }
	}


	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(sinoImage, &sinoArray_d)) return 1;

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* sino_data = (float*) sinoImage->data;
        for (int i=0; i<sinoImage->nvox; i++)
            if (sino_data[i] < 0)
                sino_data[i] = 0;

	/*Free*/
	cudaCommon_free(&activityArray_d);
	cudaCommon_free((void **)&rotatedArray_d);
	cudaCommon_free((void **)&sinoArray_d);
	cudaCommon_free((void **)&mask_d);
	cudaCommon_free((void **)&positionFieldImageArray_d);
        if (psfImage != NULL)
            {
            cudaCommon_free((void **)&psfArray_d);
            if (separable_psf)
                cudaCommon_free((void **)&psfSeparatedArray_d);
            }
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
        float     *psfSeparatedArray_d;
        int       psf_size[3];
        int       image_size[3];
        int       separable_psf=0;

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
            psf_size[0] = psfImage->dim[1];
            psf_size[1] = psfImage->dim[2];
            psf_size[2] = psfImage->dim[3];
            image_size[0] = accumulatorImage->dim[1];
            image_size[1] = accumulatorImage->dim[2];
            image_size[2] = accumulatorImage->dim[3];

            if (psf_size[0]<= (MAX_SEPARABLE_KERNEL_RADIUS*2)+1)
                separable_psf=1;
            }

        if (separable_psf)
            {
            CUDA_SAFE_CALL(cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)));
            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            CUDA_SAFE_CALL(cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice));
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
                    }
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
                    {
                    et_line_backproject_gpu(	&sinoArray_d,
						&temp_backprojection_d,
						cam,
						accumulatorImage);
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

		// Accumulate //
		et_accumulate_gpu(		&rotatedArray_d,
						&accumulatorArray_d,
						accumulatorImage );
	}

	/* Transfer result back to host */
	if(cudaCommon_transferFromDeviceToNifti(accumulatorImage, &accumulatorArray_d)) return 1; 

        /* Truncate negative values: small negative values may be found due to FFT and IFFT */
        float* accumulator_data = (float*) accumulatorImage->data;
        for (int i=0; i<accumulatorImage->nvox; i++)
            if (accumulator_data[i] < 0)
                accumulator_data[i] = 0;

	/*Free*/
        cudaCommon_free(&backprojectionArray_d);
        if (attenuationImage != NULL)
            {
            cudaCommon_free(&attenuationArray_d);
            cudaCommon_free((void **)&rotatedAttenuationArray_d);
            }
	cudaCommon_free((void **)&rotatedArray_d);
	cudaCommon_free((void **)&sinoArray_d);
	cudaCommon_free((void **)&accumulatorArray_d);
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


int et_project_backproject_gpu(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_theta_x, float *cameras_theta_y, float *cameras_theta_z)
{
	int status = 1;
	//status = rotate();
	return status;
}


int et_fisher_grid_gpu(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *fisherpriorImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float epsilon, float background, float background_attenuation)
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

    if (epsilon<=eps) epsilon=eps;
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
        status = et_project_gpu(inputImage, invsinogramImage, psfImage, attenuationImage, cameras, n_cameras, background, background_attenuation);
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
        et_create_rotation_matrix(affineTransformation, cameras[0*n_cameras+cam], cameras[1*n_cameras+cam], cameras[2*n_cameras+cam], center_x, center_y, center_z);
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


int et_gradient_attenuation_gpu(nifti_image *gradientImage, nifti_image *sinoImage, nifti_image *activityImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation) 
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
	CUDA_SAFE_CALL(cudaMemcpy(mask_d, mask_h, gradientImage->nvox*sizeof(int), cudaMemcpyHostToDevice));
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
            CUDA_SAFE_CALL(cudaMalloc((void **)&psfSeparatedArray_d, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float)));
            float *psfSeparatedArray_h = (float*) malloc((psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float));
            float psf_norm;
            for (int n=0; n<psf_size[2];n++) {
                psf_norm = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + (psf_size[0]-1)/2];
                for (int i=0;i<psf_size[0];i++) {
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 * psf_size[0] + i];
                    psfSeparatedArray_h[(psf_size[0]+psf_size[1])*n + psf_size[0] + i] = ((float*)psfImage->data)[psf_size[0]*psf_size[1]*n + (psf_size[0]-1)/2 + i * psf_size[0]] / psf_norm;
                    }
                }
            CUDA_SAFE_CALL(cudaMemcpy(psfSeparatedArray_d, psfSeparatedArray_h, (psf_size[0]+psf_size[1])*psf_size[2]*sizeof(float), cudaMemcpyHostToDevice));
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
						center_z);
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

		// Rotate backprojection //
		et_create_rotation_matrix(	affineTransformation,
						-cameras[0*n_cameras+cam],
						-cameras[1*n_cameras+cam],
						-cameras[2*n_cameras+cam],
						center_x,
						center_y, 
						center_z);
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
        for (int i=0; i<gradientImage->nvox; i++)
            if (gradient_data[i] < 0)
                gradient_data[i] = 0;

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





