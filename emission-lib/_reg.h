

#include "_reg_tools.h"

int reg_gradient_NMI_nodes(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointsImage, nifti_image* gradientImage, int binning);
int reg_gradient_NMI_nodes_gpu(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointsImage, nifti_image* gradientImage, int binning);

int reg_gradient_voxel_to_nodes(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage);
int reg_gradient_voxel_to_nodes_gpu(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage);

int reg_resample_spline(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage);
int reg_resample_spline_gpu(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage);

nifti_image* reg_initialize_control_points(int control_point_size[], float gridSpacing[]);
nifti_image *reg_initialize_image(int image_size[], float image_spacing[]);
nifti_image *reg_initialize_deformation_field(int image_size[], float image_spacing[]);
int free_nifti_image_except_data(nifti_image *image);

