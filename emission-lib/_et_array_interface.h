
#include "_et.h"

extern "C" int et_array_affine(float *image_ptr, int *image_size, float *transformed_image_ptr, int *transformed_image_size, float *affine_ptr, int *affine_size, float background, int GPU);
extern "C" int et_array_rotate(float *image_ptr, int *size_ptr, float *rotated_image_ptr, float *angles_ptr, float *centers_ptr, float background, int GPU);
extern "C" int et_array_project(float *activity, int *activity_size, float *sinogram, int *sinogram_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU);
extern "C" int et_array_backproject(float *sino, int *sino_size, float *bkpr, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *attenuation, int *attenuation_size, float background, float background_attenuation, int GPU); 
extern "C" int et_mlem_spect(float *sinogram_data, int size_x, int size_y, int n_cameras, float firstcamera, float lastcamera, int iterations, int use_psf, int use_ddpsf, int psf_size_x, int psf_size_y, float *psf_data, int use_attenuation, float *attenuation_data, float *activity_data, int GPU);
//extern "C" int et_array_joint_histogram(float *matrix_A, float *matrix_B, int *joint_histogram, int matrix_dimensions, int *matrix_size, int histogram_size, float min_A, float max_A, float min_B, float max_B, int GPU);
extern "C" int et_array_convolve(float *image_ptr, int *image_size, float *out_ptr, int *out_size, float *psf_ptr, int *psf_size, int GPU);
extern "C" int et_array_list_gpus(int *gpu_count, int *gpus_info_array);
extern "C" int et_array_set_gpu(int id);
extern "C" int et_array_fisher_grid(float *activity_ptr, int *activity_size, float *cameras_ptr, int *cameras_size, float *psf_ptr, int *psf_size, float *grid_ptr, float *fisher_ptr, int *fisher_size, float *attenuation_ptr, int *attenuation_size, float epsilon, float background, float background_attenuation, int enable_gpu);
extern "C" int et_array_fisher_grid_projection(float *sinogram_ptr, int *sinogram_size, int *bkpr_size, float *cameras, int *cameras_size, float *psf, int *psf_size, float *grid_ptr, float *fisher_ptr, int *fisher_size, float *attenuation, int *attenuation_size, float epsilon, float background, float background_attenuation, int GPU);
