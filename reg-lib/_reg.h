

#include "_reg_tools.h"

nifti_image* reg_initialize_control_points(int control_point_size[], float gridSpacing[]);
nifti_image *reg_initialize_image(int image_size[], float image_spacing[]);
nifti_image *reg_initialize_deformation_field(int image_size[], float image_spacing[]);
int free_nifti_image_except_data(nifti_image *image);

int reg_gradient_NMI_nodes(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointsImage, nifti_image* gradientImage, int binning);
int reg_gradient_voxel_to_nodes(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage);
int reg_resample_spline(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage);
int reg_image_gradient(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage);
int reg_ssd_gradient(nifti_image *ssdGradientImage, nifti_image *targetImage, nifti_image *sourceImage, nifti_image *gradientImage, int smoothingRadius[]);
int reg_gaussian_smooth(nifti_image *niftiImage, float sigma);
int reg_scale_amplitude(nifti_image *niftiImage, float min_value, float max_value);
int reg_gradient_jacobian_determinant(nifti_image *controlPointImage, nifti_image *nodeGradientImage, int image_size[], float image_spacing[], float weight);
int reg_gradient_bending_energy(nifti_image *controlPointImage, nifti_image *nodeGradientImage, int image_size[], float image_spacing[], float weight);

#ifdef _USE_CUDA
int reg_image_gradient_gpu(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage);
int reg_gradient_NMI_nodes_gpu(nifti_image* targetImage, nifti_image* sourceImage, nifti_image* controlPointsImage, nifti_image* gradientImage, int binning);
int reg_gradient_voxel_to_nodes_gpu(nifti_image *cpGradientImage, nifti_image *gradientImage, nifti_image *controlPointImage);
int reg_resample_spline_gpu(nifti_image *resultImage, nifti_image *sourceImage, nifti_image *controlPointImage);
int reg_ssd_gradient_gpu(nifti_image *ssdGradientImage, nifti_image *targetImage, nifti_image *sourceImage, nifti_image *gradientImage, int smoothingRadius[]);
int reg_gaussian_smooth_gpu(nifti_image *niftiImage, float sigma);
int reg_scale_amplitude_gpu(nifti_image *niftiImage, float min_value, float max_value);
int reg_gradient_jacobian_determinant_gpu(nifti_image *controlPointImage, nifti_image *nodeGradientImage, int image_size[], float image_spacing[], float weight);
int reg_gradient_bending_energy_gpu(nifti_image *controlPointImage, nifti_image *nodeGradientImage, int image_size[], float image_spacing[], float weight);
#endif

