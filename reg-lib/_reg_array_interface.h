

#include "_reg.h"

int reg_array_compute_control_point_size(int size[], int imageSize[], float gridSpacing[]);
int reg_array_bspline_initialiseControlPointGridWithAffine(float *controlPoints, float *affineTransformation, int imageSize[], float gridSpacing[]);
int reg_array_gradient_NMI_nodes(float *target, float *source, float *gradient, int *image_size, float *control_points, float gridSpacing[], int binning, int enable_gpu);
int reg_array_gradient_voxel_to_nodes(float *cp_gradient_ptr, float *gradient_ptr, float *cp_ptr, int image_size[], int cp_size[], float grid_spacing[], int enable_gpu);
int reg_array_resample_spline(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_points_size[], float *spacing, int enable_gpu);
int reg_array_image_gradient(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_points_size[], float *spacing, int enable_gpu);
int reg_array_ssd_gradient(float *ssd_gradient_ptr, float *target_ptr, float *source_ptr, float *gradient_ptr, int image_size[], int smoothing_radius[], int enable_gpu);
int reg_array_gaussian_smooth(float *image_ptr, int image_size[], float smoothing_sigma, int enable_gpu); 
int reg_array_scale_amplitude(float *image_ptr, int image_size[], float min_value, float max_value, int enable_gpu); 
int reg_array_gradient_jacobian_determinant(float *nodes_gradient_ptr, float *control_points_ptr, int image_size[], int cp_size[], float cp_spacing[], float weight, int GPU);
int reg_array_gradient_bending_energy_gpu(float *nodes_gradient_ptr, float *control_points_ptr, int image_size[], int cp_size[], float cp_spacing[], float weight, int GPU);
