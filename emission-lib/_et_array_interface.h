/*
 *  _et_array_interface.h
 *  
 *  NiftyRec
 *  Stefano Pedemonte, Oct. 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


/**
* Plain C-types interface to NiftyRec
*/

#include "_et.h"

extern "C" int et_array_affine(float *image_ptr, int *image_size, float *transformed_image_ptr, int *transformed_image_size, float *affine_ptr, int *affine_size, float background, int GPU);
extern "C" int et_array_rotate(float *image_ptr, int *size_ptr, float *rotated_image_ptr, float *angles_ptr, float *centers_ptr, float background, int GPU);
extern "C" int et_array_project(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU, int truncate_negative_values);
extern "C" int et_array_backproject(float *sino, int *sino_size, float *bkpr, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU, int truncate_negative_values); 
extern "C" int et_array_calculate_size_psf(unsigned int *psf_size_x, unsigned int *psf_size_y, float fwhm_pixels_dist0, float sensitivity0, float dist0, float fwhm_pixels_dist1, float sensitivity1, float dist1); 
extern "C" int et_array_make_psf(float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float fwhm_pixels_dist0, float sensitivity0, float dist0, float fwhm_pixels_dist1, float sensitivity1, float dist1, unsigned int N_psf_planes); 
extern "C" int et_array_make_cameras(float *cameras_data, float firstcamera_deg, float lastcamera_deg, unsigned int n_cameras, unsigned int rotation_axis);
extern "C" int et_array_osem_spect(float *activity_data, unsigned int size_x, unsigned int size_y, unsigned int subset_order, float *sinogram_data, int n_cameras, float firstcamera, float lastcamera, unsigned int rotation_axis, int iterations, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);
extern "C" int et_array_osem_step(float *activity_data, int size_x, int size_y, unsigned int subset_order, float *sinogram_data, int n_cameras, float *cameras_data, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);
extern "C" int et_array_mlem_spect(float *activity_data, unsigned int size_x, unsigned int size_y, float *sinogram_data, int n_cameras, float firstcamera, float lastcamera, unsigned int rotation_axis, int iterations, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);
extern "C" int et_array_mlem_step(float *activity_data, int size_x, int size_y, float *sinogram_data, int n_cameras, float *cameras_data, int use_attenuation, float *attenuation_data, int use_psf, float *psf_data, unsigned int psf_size_x, unsigned int psf_size_y, float background, float background_attenuation, float epsilon, int use_gpu);
//extern "C" int et_array_joint_histogram(float *matrix_A, float *matrix_B, int *joint_histogram, int matrix_dimensions, int *matrix_size, int histogram_size, float min_A, float max_A, float min_B, float max_B, int GPU);
extern "C" int et_array_convolve(float *image_ptr, int *image_size, float *out_ptr, int *out_size, float *psf_ptr, int *psf_size, int GPU);
extern "C" int et_array_list_gpus(int *gpu_count, int *gpus_info_array);
extern "C" int et_array_set_gpu(int id);
extern "C" int et_array_fisher_grid(float *activity_ptr, int *activity_size, float *cameras_ptr, int *cameras_size, float *psf_ptr, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation_ptr, int *attenuation_size, float epsilon, float background, float background_attenuation, int enable_gpu);
extern "C" int et_array_fisher_grid_projection(float *sinogram_ptr, int *sinogram_size, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *grid_ptr, float *fisher_ptr, float *fisher_prior_ptr, int *fisher_size, float *attenuation, int *attenuation_size, float epsilon, float background, float background_attenuation, int GPU);
extern "C" int et_array_gradient_attenuation(float *sino, int *sino_size, float *activity, int *activity_size, float *gradient, int *gradient_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU, int truncate_negative_values);
extern "C" int et_array_project_partial(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *partialsum, int *partialsum_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU, int truncate_negative_values, int do_rotate_partial);
extern "C" int et_array_isinstalled();
extern "C" int et_array_reset_gpu();
extern "C" int et_array_histogram_weighted(float *inputdata_ptr, float *weights_ptr, float *histogram_ptr, int N, int N_classes, int N_bins, double min_value, double max_value); 

