
#include "_reg.h"

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_reg_bspline.h"
#include "_reg_mutualinformation.h"
#include "_reg_ssd.h"
#include "_reg_tools.h"
#include "float.h"

#ifdef _USE_CUDA
	#include "_reg_cudaCommon.h"
	#include "_reg_resampling_gpu.h"
	#include "_reg_affineTransformation_gpu.h"
	#include "_reg_bspline_gpu.h"
	#include "_reg_mutualinformation_gpu.h"
	#include "_reg_tools_gpu.h"
#endif

#define JH_PW_APPROX 2


int reg_gradient_NMI_nodes(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointImage, nifti_image* nodeGradientImage, int binning)
{
	int status = 1;
	return status;
}


int reg_gradient_NMI_nodes_gpu(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointImage, nifti_image* nodeGradientImage, int binning)
{
	float *targetImageArray_d=NULL;
	cudaArray *sourceImageArray_d=NULL;
	float4 *controlPointImageArray_d=NULL;
	float4 *positionFieldImageArray_d=NULL;
	int *targetMask_d=NULL;
	float4 *resultGradientArray_d=NULL;
	float4 *voxelNMIGradientArray_d=NULL;
	float4 *nodeNMIGradientArray_d=NULL;
	float *resultImageArray_d=NULL;
	float *logJointHistogram_d=NULL;

        float *probaJointHistogram = (float *)malloc(binning*(binning+2)*sizeof(float));
        float *logJointHistogram = (float *)malloc(binning*(binning+2)*sizeof(float));
	float *entropies = (float *)malloc(4*sizeof(float));
	double *entropies_double = (double *)malloc(4*sizeof(double));	

	int smoothingRadius[3];
	smoothingRadius[0] = (int)floor( 2.0*controlPointImage->dx/targetImage->dx );
	smoothingRadius[1] = (int)floor( 2.0*controlPointImage->dy/targetImage->dy );
	smoothingRadius[2] = (int)floor( 2.0*controlPointImage->dz/targetImage->dz );

	int activeVoxelNumber = targetImage->nvox;

	if(cudaCommon_allocateArrayToDevice<float>(&targetImageArray_d, targetImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&targetImageArray_d, targetImage)) return 1;
fprintf(stderr,"\n1.\n");
	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;
fprintf(stderr,"\n1..\n");
	if(cudaCommon_allocateArrayToDevice<float>(&resultImageArray_d, targetImage->dim)) return 1;
fprintf(stderr,"\n1...\n");
        fprintf(stderr,"\nControl point image dim: %d %d %d %d %d %d\n",controlPointImage->dim[0],controlPointImage->dim[1],controlPointImage->dim[2],controlPointImage->dim[3],controlPointImage->dim[4],controlPointImage->dim[5]);
	if(cudaCommon_allocateArrayToDevice<float4>(&controlPointImageArray_d, controlPointImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float4>(&controlPointImageArray_d,controlPointImage)) return 1;
fprintf(stderr,"\n1....\n");
	CUDA_SAFE_CALL(cudaMalloc((void **)&logJointHistogram_d, binning*(binning+2)*sizeof(float)));
	if(cudaCommon_allocateArrayToDevice<float4>(&voxelNMIGradientArray_d, sourceImage->dim)) return 1; //result
	if(cudaCommon_allocateArrayToDevice<float4>(&nodeNMIGradientArray_d, controlPointImage->dim)) return 1;
fprintf(stderr,"\n1.....\n");
CUDA_SAFE_CALL(cudaMalloc((void **)&targetMask_d, activeVoxelNumber*sizeof(int)));
int *mask_h=(int *)malloc(activeVoxelNumber*sizeof(int));
for(int i=0; i<activeVoxelNumber; i++) mask_h[i]=i;
CUDA_SAFE_CALL(cudaMemcpy(targetMask_d, mask_h, activeVoxelNumber*sizeof(int), cudaMemcpyHostToDevice));
fprintf(stderr,"\n1......\n");
	CUDA_SAFE_CALL(cudaMalloc((void **)&positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&resultGradientArray_d, activeVoxelNumber*sizeof(float4)));
fprintf(stderr,"\n3\n");
	nifti_image *tempImage=NULL;
	tempImage = nifti_copy_nim_info(sourceImage);
	tempImage->datatype = NIFTI_TYPE_FLOAT32;
	tempImage->nbyper = sizeof(float);
	tempImage->data = (float *)malloc(sourceImage->nvox*sizeof(float));
fprintf(stderr,"\n4\n");
	// generate the position field //
	reg_bspline_gpu(			controlPointImage,
						targetImage,
						&controlPointImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber);
fprintf(stderr,"\n5\n");
	// Resample the source image //
	reg_resampleSourceImage_gpu(		targetImage, //result
						sourceImage,
						&resultImageArray_d,
						&sourceImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber,
						0);
fprintf(stderr,"\n6\n");
	// The result image is transfered back to the host //
	if(cudaCommon_transferFromDeviceToNifti(tempImage, &resultImageArray_d)) return 1;

	/* Create joint histogram and transfer to the device*/
	reg_getEntropies<float>(		targetImage,
						tempImage, //resultImage
						2,
						binning,
						probaJointHistogram,
						logJointHistogram,
						entropies,
						mask_h);

//for(int i=0; i<binning*(binning+2); i++)
//    fprintf(stderr,"%f ",logJointHistogram[i]);
//fprintf(stderr,"\n");

	entropies_double[0]=entropies[0];entropies_double[1]=entropies[1];entropies_double[2]=entropies[2];entropies_double[3]=entropies[3];
fprintf(stderr,"\nBinning: %d  Entropies: %f %f %f %f",binning, entropies_double[0],entropies_double[1],entropies_double[2],entropies_double[3]);

	// Tranfer the histogram to the device //
	CUDA_SAFE_CALL(cudaMemcpy(logJointHistogram_d, logJointHistogram, binning*(binning+2)*sizeof(float), cudaMemcpyHostToDevice));

	// NMI Gradient //
fprintf(stderr,"%d\n",targetImage->nvox);
	reg_getSourceImageGradient_gpu(						targetImage,
										sourceImage,
										&sourceImageArray_d,
										&positionFieldImageArray_d,
										&resultGradientArray_d,
										activeVoxelNumber);

	reg_getVoxelBasedNMIGradientUsingPW_gpu(				targetImage,
										tempImage, 
										&targetImageArray_d,     // V
										&resultImageArray_d,     // V
										&resultGradientArray_d,  // V
										&logJointHistogram_d,    // V
										&voxelNMIGradientArray_d,
										&targetMask_d,           // V
										activeVoxelNumber,       //
										entropies_double,        //
										binning);                //

fprintf(stderr,"\ngradientImage: %d %d %d\n",nodeGradientImage->nx,nodeGradientImage->ny,nodeGradientImage->nz);
fprintf(stderr,"sourceImage: %d %d %d\n",sourceImage->nx,sourceImage->ny,sourceImage->nz);
fprintf(stderr,"Smoothing radius: %d %d %d\n", smoothingRadius[0],smoothingRadius[1],smoothingRadius[2]);
fprintf(stderr,"Image size: %d %d %d\n", sourceImage->nx, sourceImage->ny, sourceImage->nz);

//	reg_smoothImageForCubicSpline_gpu(					sourceImage,
//										&voxelNMIGradientArray_d,
//										smoothingRadius);

	reg_voxelCentric2NodeCentric_gpu(					targetImage,
										controlPointImage,
										&voxelNMIGradientArray_d,
										&nodeNMIGradientArray_d);

fprintf(stderr,"Voxel number: %d\n",activeVoxelNumber);
fprintf(stderr,"Nodes size:   %d\n",nodeGradientImage->nvox);

CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, nodeNMIGradientArray_d, nodeGradientImage->nvox*sizeof(float), cudaMemcpyDeviceToHost));
//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, logJointHistogram_d, binning*(binning+2)*sizeof(float), cudaMemcpyDeviceToHost));
//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, voxelNMIGradientArray_d, activeVoxelNumber*4*sizeof(float), cudaMemcpyDeviceToHost));

//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, resultGradientArray_d, activeVoxelNumber*4*sizeof(float), cudaMemcpyDeviceToHost));
//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, positionFieldImageArray_d, activeVoxelNumber*4*sizeof(float), cudaMemcpyDeviceToHost));
//CUDA_SAFE_CALL(cudaMemcpy(nodeGradientImage->data, targetMask_d, activeVoxelNumber*1*sizeof(float), cudaMemcpyDeviceToHost));

	cudaCommon_free( (void **)&targetImageArray_d );
	cudaCommon_free( &sourceImageArray_d );
	cudaCommon_free( (void **)&controlPointImageArray_d );
	cudaCommon_free( (void **)&resultImageArray_d );
	cudaCommon_free( (void **)&positionFieldImageArray_d );
	CUDA_SAFE_CALL(cudaFree(targetMask_d));
	cudaCommon_free((void **)&resultGradientArray_d);
	cudaCommon_free((void **)&voxelNMIGradientArray_d);
	cudaCommon_free((void **)&nodeNMIGradientArray_d);
	cudaCommon_free((void **)&logJointHistogram_d);

        free(probaJointHistogram);
        free(logJointHistogram);
	free(entropies);
	free(entropies_double);
        free(mask_h);

        nifti_image_free(tempImage);

fprintf(stderr,"_reg Done \n\n");
	return 0;
}




int reg_gradient_voxel_to_nodes(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage)
{
	int status = 1;
	return status;
}



int reg_gradient_voxel_to_nodes_gpu(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage)
{
	int status = 1;
	float4 *voxelNMIGradientArray_d=NULL;
	float4 *nodeNMIGradientArray_d=NULL;

	if(cudaCommon_allocateArrayToDevice(&voxelNMIGradientArray_d, gradientImage->dim)) return 1;
	if(cudaCommon_allocateArrayToDevice(&nodeNMIGradientArray_d, controlPointImage->dim)) return 1;	

	int smoothingRadius[3];
	smoothingRadius[0] = (int)floor( 2.0*controlPointImage->dx/gradientImage->dx );
	smoothingRadius[1] = (int)floor( 2.0*controlPointImage->dy/gradientImage->dy );
	smoothingRadius[2] = (int)floor( 2.0*controlPointImage->dz/gradientImage->dz );

fprintf(stderr, "\n Image delta: %f %f %f",gradientImage->dx,gradientImage->dy,gradientImage->dz);
fprintf(stderr, "\n CP    delta: %f %f %f\n",controlPointImage->dx,controlPointImage->dy,controlPointImage->dz);

	if(cudaCommon_transferNiftiToArrayOnDevice<float4> (&voxelNMIGradientArray_d,gradientImage)) return 1;

	reg_smoothImageForCubicSpline_gpu(					gradientImage,
										&voxelNMIGradientArray_d,
										smoothingRadius);

	reg_voxelCentric2NodeCentric_gpu(					gradientImage,
										controlPointImage,
										&voxelNMIGradientArray_d,
										&nodeNMIGradientArray_d);

	CUDA_SAFE_CALL(cudaMemcpy(cpGradientImage->data, nodeNMIGradientArray_d, cpGradientImage->nvox*sizeof(float), cudaMemcpyDeviceToHost));

	cudaCommon_free((void **)&voxelNMIGradientArray_d);
	cudaCommon_free((void **)&nodeNMIGradientArray_d);
	status = 0;
	return status;
}



int reg_resample_spline(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage)
{
	int status = 1;
	return status;
}


int reg_resample_spline_gpu(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage)
{
	int status = 1;

	cudaArray *sourceImageArray_d=NULL;
	float4 *controlPointImageArray_d=NULL;
	float4 *positionFieldImageArray_d=NULL;
	int *targetMask_d=NULL;

	float *resultImageArray_d=NULL;

	int *targetMask;
	targetMask = (int *)malloc(sourceImage->nvox*sizeof(int));
	int activeVoxelNumber=0;
	for(unsigned int i=0; i<sourceImage->nvox; i++)
		targetMask[i]=i;
	activeVoxelNumber=sourceImage->nvox;

        int binning = 128; //FIXME make it an argument

	if(cudaCommon_allocateArrayToDevice<float>(&sourceImageArray_d, sourceImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float>(&sourceImageArray_d,sourceImage)) return 1;

	if(cudaCommon_allocateArrayToDevice<float>(&resultImageArray_d, sourceImage->dim)) return 1;

	if(cudaCommon_allocateArrayToDevice<float4>(&controlPointImageArray_d, controlPointImage->dim)) return 1;
	if(cudaCommon_transferNiftiToArrayOnDevice<float4>(&controlPointImageArray_d,controlPointImage)) return 1;

	// Index of the active voxel is stored
	int *targetMask_h; CUDA_SAFE_CALL(cudaMallocHost((void **)&targetMask_h, activeVoxelNumber*sizeof(int)));
	int *targetMask_h_ptr = &targetMask_h[0];
	for(unsigned int i=0;i<sourceImage->nvox;i++)
		{
		if(targetMask[i]!=-1) 
			*targetMask_h_ptr++=i;
		}
	CUDA_SAFE_CALL(cudaMalloc((void **)&targetMask_d, activeVoxelNumber*sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(targetMask_d, targetMask_h, activeVoxelNumber*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaFreeHost(targetMask_h));

	CUDA_SAFE_CALL(cudaMalloc((void **)&positionFieldImageArray_d, activeVoxelNumber*sizeof(float4)));

	// generate the position field //
	reg_bspline_gpu(			controlPointImage,
						sourceImage,  //targetImage
						&controlPointImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber);

	// Resample the source image //
	reg_resampleSourceImage_gpu(		resultImage,
						sourceImage,
						&resultImageArray_d,
						&sourceImageArray_d,
						&positionFieldImageArray_d,
						&targetMask_d,
						activeVoxelNumber,
						0);

	if(cudaCommon_transferFromDeviceToNifti(resultImage, &resultImageArray_d)) return 1;

	cudaCommon_free( &sourceImageArray_d );
	cudaCommon_free( (void **)&controlPointImageArray_d );
	cudaCommon_free( (void **)&resultImageArray_d );
	cudaCommon_free( (void **)&positionFieldImageArray_d );
	CUDA_SAFE_CALL(cudaFree(targetMask_d));

	status = 0;
	return status;
}



nifti_image* reg_initialize_control_points(int control_point_size[], float gridSpacing[])
{
	int dim_cpp[8];
	dim_cpp[0]=5;
	dim_cpp[1]=control_point_size[0];
	dim_cpp[2]=control_point_size[1];
	dim_cpp[3]=control_point_size[2];
	dim_cpp[5]=4;
	dim_cpp[4]=dim_cpp[6]=dim_cpp[7]=1;
	nifti_image* controlPointImage = nifti_make_new_nim(dim_cpp, NIFTI_TYPE_FLOAT32, true);
	controlPointImage->cal_min=0;
	controlPointImage->cal_max=0;
	controlPointImage->pixdim[0]=1.0f;
	controlPointImage->pixdim[1]=controlPointImage->dx=gridSpacing[0];
	controlPointImage->pixdim[2]=controlPointImage->dy=gridSpacing[1];
	controlPointImage->pixdim[3]=controlPointImage->dz=gridSpacing[2];
	controlPointImage->pixdim[4]=controlPointImage->dt=1.0f;
	controlPointImage->pixdim[5]=controlPointImage->du=1.0f;
	controlPointImage->pixdim[6]=controlPointImage->dv=1.0f;
	controlPointImage->pixdim[7]=controlPointImage->dw=1.0f;
	controlPointImage->qform_code=1;
        controlPointImage->quatern_b=0.f;
        controlPointImage->quatern_c=0.f;
        controlPointImage->quatern_d=0.f;
        controlPointImage->qfac=1.f;
	controlPointImage->qoffset_x = -gridSpacing[0];
	controlPointImage->qoffset_y = -gridSpacing[1];
        controlPointImage->qoffset_z = -gridSpacing[2];

        controlPointImage->qto_xyz = nifti_quatern_to_mat44(0.f, 0.f, 0.f,
		-gridSpacing[0], -gridSpacing[1], -gridSpacing[2],
		gridSpacing[0], gridSpacing[1], gridSpacing[2], 1.f);
        controlPointImage->qto_ijk = nifti_mat44_inverse(controlPointImage->qto_xyz);
	controlPointImage->sform_code=0;
	return controlPointImage;
}


nifti_image *reg_initialize_image(int image_size[], float image_spacing[])
{
	int dim[8];
	dim[0]    = 3;
	dim[1]    = image_size[0];
	dim[2]    = image_size[1];
	dim[3]    = image_size[2];
	dim[4]    = 1;
	dim[5]    = 1;
	dim[6]    = 1;
	dim[7]    = 1;
	nifti_image *image = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
//	image->cal_min=0;
//	image->cal_max=0;
	image->pixdim[0]=1.0f;
	image->pixdim[1]=image->dx=image_spacing[0];
	image->pixdim[2]=image->dy=image_spacing[1];
	image->pixdim[3]=image->dz=image_spacing[2];
	image->pixdim[4]=image->dt=1.0f;
	image->pixdim[5]=image->du=1.0f;
	image->pixdim[6]=image->dv=1.0f;
	image->pixdim[7]=image->dw=1.0f;
	image->qform_code=1;
        image->quatern_b=0.f;
        image->quatern_c=0.f;
        image->quatern_d=0.f;
        image->qfac=1.f;
	image->qoffset_x = 0.f;
	image->qoffset_y = 0.f;
        image->qoffset_z = 0.f;

        image->qto_xyz = nifti_quatern_to_mat44(0.f, 0.f, 0.f,
		0, 0, 0,
		image_spacing[0], image_spacing[1], image_spacing[2], 1.f);
fprintf(stderr,"\nTransformation matrix: %f %f %f %f",image->qto_xyz.m[0][0],image->qto_xyz.m[0][1],image->qto_xyz.m[0][2],image->qto_xyz.m[0][3]);
fprintf(stderr,"\nTransformation matrix: %f %f %f %f",image->qto_xyz.m[1][0],image->qto_xyz.m[1][1],image->qto_xyz.m[1][2],image->qto_xyz.m[1][3]);
fprintf(stderr,"\nTransformation matrix: %f %f %f %f",image->qto_xyz.m[2][0],image->qto_xyz.m[2][1],image->qto_xyz.m[2][2],image->qto_xyz.m[2][3]);
fprintf(stderr,"\nTransformation matrix: %f %f %f %f",image->qto_xyz.m[3][0],image->qto_xyz.m[3][1],image->qto_xyz.m[3][2],image->qto_xyz.m[3][3]);
        image->qto_ijk = nifti_mat44_inverse(image->qto_xyz);
	image->sform_code=0;
	return image;	
}


nifti_image *reg_initialize_deformation_field(int image_size[], float image_spacing[])
{
	nifti_image *positionFieldImage = reg_initialize_image(image_size, image_spacing);
        positionFieldImage->dim[0]=positionFieldImage->ndim=5;
        positionFieldImage->dim[4]=positionFieldImage->nt=1;positionFieldImage->pixdim[4]=positionFieldImage->dt=1.0;
	positionFieldImage->dim[5]=positionFieldImage->nu=4;
        positionFieldImage->pixdim[5]=positionFieldImage->du=1.0;
        positionFieldImage->dim[6]=positionFieldImage->nv=1;positionFieldImage->pixdim[6]=positionFieldImage->dv=1.0;
        positionFieldImage->dim[7]=positionFieldImage->nw=1;positionFieldImage->pixdim[7]=positionFieldImage->dw=1.0;
        positionFieldImage->nvox=positionFieldImage->nx*positionFieldImage->ny*positionFieldImage->nz*positionFieldImage->nt*positionFieldImage->nu;
        positionFieldImage->datatype = NIFTI_TYPE_FLOAT32;
        positionFieldImage->nbyper = sizeof(float);
	return positionFieldImage;
}



int free_nifti_image_except_data(nifti_image *image)
{
	if( image->fname != NULL ) free(image->fname) ;
	if( image->iname != NULL ) free(image->iname) ;
	(void)nifti_free_extensions( image ) ;
	free(image) ;
	return 0;
}

