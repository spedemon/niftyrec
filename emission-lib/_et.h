
#define SIZE_OF_INFO 5
#define MAX_DEVICES 16

#define eps 0.000000000001f

#include "_reg_tools.h"

#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_et_line_integral.h"
#include "_et_line_integral_attenuated.h"
#include "_et_accumulate.h"
#include "_et_line_backproject.h"
#include "_et_line_backproject_attenuated.h"
#include "_et_clear_accumulator.h"
#include "_et_convolve2D.h"
//#include "_et_joint_histogram.h"

#define ET_ERROR_BADGRID 2

#ifdef _USE_CUDA
#include "_reg_cudaCommon.h"
#include "_reg_resampling_gpu.h"
#include "_reg_affineTransformation_gpu.h"
#include "_et_line_integral_gpu.h"
#include "_et_line_integral_attenuated_gpu.h"
#include "_et_accumulate_gpu.h"
#include "_et_line_backproject_gpu.h"
#include "_et_line_backproject_attenuated_gpu.h"
#include "_et_clear_accumulator_gpu.h"
#include "_et_convolveFFT2D_gpu.h"
#include "_et_convolveSeparable2D_gpu.h"
//#include "_et_joint_histogram_gpu.h"

int et_rotate_gpu(nifti_image *sourceImage, nifti_image *resultImage, float alpha, float beta, float gamma, float center_x, float center_y, float center_z, float background);
int et_project_gpu(nifti_image *activity, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float backgroundAttenuation);
int et_backproject_gpu(nifti_image *sinogram, nifti_image *accumulator, nifti_image *psf, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation);
int et_joint_histogram_gpu(nifti_image *matrix_A_Image, nifti_image *matrix_B_Image, nifti_image *joint_histogram_Image, float min_A, float max_A, float min_B, float max_B);
int et_project_backproject_gpu(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_alpha, float *cameras_beta, float *cameras_gamma);
int et_affine_gpu(nifti_image *sourceImage, nifti_image *resultImage, mat44 *affineTransformation, float background);
int et_convolve_gpu(nifti_image *inImage, nifti_image *outImage, nifti_image *psfImage);
int et_fisher_grid_gpu(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *priorfisherImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras_array, int n_cameras, float epsilon, float background, float background_attenuation); 
int et_list_gpus(int *device_count_out, int *devices);
int et_set_gpu(int id);
#endif

int et_rotate(nifti_image *sourceImage, nifti_image *resultImage, float alpha, float beta, float gamma, float center_x, float center_y, float center_z, float background);
int et_project(nifti_image *activityImage, nifti_image *sinoImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation);
int et_backproject(nifti_image *sinogram, nifti_image *accumulator, nifti_image *psf, nifti_image *attenuationImage, float *cameras, int n_cameras, float background, float background_attenuation);
int et_joint_histogram(nifti_image *matrix_A_Image, nifti_image *matrix_B_Image, nifti_image *joint_histogram_Image, float min_A, float max_A, float min_B, float max_B);
int et_project_backproject(nifti_image *activity, nifti_image *sino, nifti_image *psf, int n_cameras, float *cameras_alpha, float *cameras_beta, float *cameras_gamma);
int et_affine(nifti_image *sourceImage, nifti_image *transformedImage, mat44 *affineTransformation, float background);
int et_convolve(nifti_image *inImage, nifti_image *outImage, nifti_image *psfImage);
int et_fisher_grid(int from_projection, nifti_image *inputImage, nifti_image *gridImage, nifti_image *fisherImage, nifti_image *priorfisherImage, nifti_image *psfImage, nifti_image *attenuationImage, float *cameras, int n_cameras, float epsilon, float background, float background_attenuation);



