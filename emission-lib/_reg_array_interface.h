

#include "_reg.h"

int reg_array_compute_control_point_size(int size[], int imageSize[], float gridSpacing[]);
int reg_array_bspline_initialiseControlPointGridWithAffine(float *controlPoints, float *affineTransformation, int imageSize[], float gridSpacing[]);
int reg_array_gradient_NMI_nodes(float *target, float *source, float *gradient, int *image_size, float *control_points, float gridSpacing[], int binning, int GPU);
int reg_array_gradient_voxel_to_nodes(float *cp_gradient_ptr, float *gradient_ptr, float *cp_ptr, int image_size[], int cp_size[], float grid_spacing[], int GPU);
int reg_array_resample_spline(float *outimage_ptr, float *image_ptr, float *nodes_ptr, int image_size[], int control_points_size[], float spacing[], int enable_gpu);
