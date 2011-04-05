
#include "_reg_array_interface.h"
#include "_reg_affineTransformation.h"
#include "_reg_bspline.h"

int reg_array_compute_control_point_size(int size[], int imageSize[], float gridSpacing[])
{
	size[0]=floor(imageSize[0]/gridSpacing[0])+4;
	size[1]=floor(imageSize[1]/gridSpacing[1])+4;
	size[2]=floor(imageSize[2]/gridSpacing[2])+4;
	return 0;
}


int reg_array_bspline_initialiseControlPointGridWithAffine(float *controlPoints, float *affineTransformation, int imageSize[], float gridSpacing[])
{
	mat44 *affineTransformation_mat44 = (mat44 *)calloc(1,sizeof(mat44));
	// Initialize affine transform matrix
	
	int k=0;
	for (unsigned int i=0; i<4; i++) {
		for (unsigned int j=0; j<4; j++) {				
			affineTransformation_mat44->m[i][j] = affineTransformation[k];
			k++;
		}
	}

	int control_point_size[3];
	reg_array_compute_control_point_size(control_point_size, imageSize, gridSpacing);
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, gridSpacing);
	controlPointImage->data = (float*) controlPoints; //FIXME: verify that nifti_make_new_nim() does not alloc ->data

	if(reg_bspline_initialiseControlPointGridWithAffine(affineTransformation_mat44, controlPointImage)) return 1;

	//Free
	free_nifti_image_except_data(controlPointImage);

	return 0;
}



int reg_array_gradient_NMI_nodes(float *target, float *source, float *gradient, int *image_size, float *control_points, float gridSpacing[], int binning, int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing);
	sourceImage->data = (float *)(source);
	
	nifti_image *targetImage = reg_initialize_image(image_size, image_spacing);
	targetImage->data = (float *)(target);

	//Create control points nifti image
	int control_point_size[3];
	reg_array_compute_control_point_size(control_point_size, image_size, gridSpacing);
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, gridSpacing);
	controlPointImage->data = (float*) (control_points); //FIXME: verify that nifti_make_new_nim() does not alloc ->data

	//Create control points gradient nifti image
	nifti_image *gradientImage=NULL;
	gradientImage = nifti_copy_nim_info(controlPointImage);
	gradientImage->datatype = NIFTI_TYPE_FLOAT32;
	gradientImage->nbyper = sizeof(float);
	gradientImage->data = (float*) (gradient);

	//Compute NMI gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_gradient_NMI_nodes_gpu(targetImage, sourceImage, controlPointImage, gradientImage, binning);
	else
		status = reg_gradient_NMI_nodes(targetImage, sourceImage, controlPointImage, gradientImage, binning);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_NMI_nodes_array: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
	status = reg_gradient_NMI_nodes(targetImage, sourceImage, controlPointImage, gradientImage, binning);
	#endif

	//Free
	free_nifti_image_except_data(targetImage);
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(controlPointImage);
	free_nifti_image_except_data(gradientImage);

	return status;
}


int reg_array_gradient_voxel_to_nodes(float *control_points_gradient, float *gradient, float *control_points, int image_size[], int cp_size[], float grid_spacing[], int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *gradientImage = reg_initialize_deformation_field(image_size, image_spacing);
	gradientImage->data = (float *)(gradient);

	//Create control points nifti image
	int control_point_size[3];
	reg_array_compute_control_point_size(control_point_size, image_size, grid_spacing);
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, grid_spacing);
	controlPointImage->data = (float*) (control_points); //FIXME: verify that nifti_make_new_nim() does not alloc ->data

	//Create control points gradient nifti image
	nifti_image *cpGradientImage=NULL;
	cpGradientImage = nifti_copy_nim_info(controlPointImage);
	cpGradientImage->datatype = NIFTI_TYPE_FLOAT32;
	cpGradientImage->nbyper = sizeof(float);
	cpGradientImage->data = (float*) control_points_gradient;

	//Compute gradient
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_gradient_voxel_to_nodes_gpu(cpGradientImage, gradientImage, controlPointImage);
	else
		status = reg_gradient_voxel_to_nodes(cpGradientImage, gradientImage, controlPointImage);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_NMI_nodes_array: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
	status = reg_gradient_voxel_to_nodes(cpGradientImage, gradientImage, controlPointImage);
	#endif

	//Free
	free_nifti_image_except_data(gradientImage);
	free_nifti_image_except_data(cpGradientImage);
	free_nifti_image_except_data(controlPointImage);

	return status;
}



int reg_array_resample_spline(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_point_size[], float spacing[], int GPU)
{
	int status = 1;

	// Allocate source and target nifti image 
	float image_spacing[3] = {1,1,1};
	nifti_image *sourceImage = reg_initialize_image(image_size, image_spacing);
	sourceImage->data = (float *)(image_ptr);

	nifti_image *resultImage = reg_initialize_image(image_size, image_spacing);
	resultImage->data = (float *)(outimage_ptr);

	//Create control points nifti image
	nifti_image *controlPointImage = reg_initialize_control_points(control_point_size, spacing);
	controlPointImage->data = (float*) nodes_ptr; //FIXME: verify that nifti_make_new_nim() does not alloc ->data

	//Resample image
	#ifdef _USE_CUDA
	if (GPU)
		status = reg_resample_spline_gpu(resultImage, sourceImage, controlPointImage);
	else
		status = reg_resample_spline(resultImage, sourceImage, controlPointImage);
	#else
	if (GPU)
		fprintf(stderr, "reg_gradient_NMI_nodes_array: No GPU support. In order to activate GPU acceleration please configure with GPU flag and compile.");
	status = reg_resample_spline(resultImage, sourceImage, controlPointImage);
	#endif	

	//Free 
	free_nifti_image_except_data(sourceImage);
	free_nifti_image_except_data(resultImage);
	free_nifti_image_except_data(controlPointImage);

	return status;
}





