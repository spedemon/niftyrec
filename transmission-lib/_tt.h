

#include "_reg_tools.h"
#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_et_line_integral.h"
#include "_et_accumulate.h"
#include "_et_line_backproject.h"
#include "_et_clear_accumulator.h"
#include "_tt_project_ray.h"
#include "_tt_common.h"
#include "_tt_backproject_ray.h"
//#include "_et_joint_histogram.h"

#ifdef _USE_CUDA
#include "_reg_cudaCommon.h"
#include "_reg_resampling_gpu.h"
#include "_reg_affineTransformation_gpu.h"
#include "_et_line_integral_gpu.h"
#include "_et_accumulate_gpu.h"
#include "_et_line_backproject_gpu.h"
#include "_et_clear_accumulator_gpu.h"
#include "_et_convolveFFT2D_gpu.h"
//#include "_tt_perspective_gpu.h"
#include "_tt_project_ray_gpu.h"
#include "_tt_backproject_ray_gpu.h"
//#include "_et_joint_histogram_gpu.h"
#endif




#ifdef _USE_CUDA
int tt_project_gpu(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background);
#endif

int tt_project(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background);

int tt_backproject_ray(float h_projections[], u_int_2 detector_pixels, u_int n_projections, float out_backprojection[], float_2 detector_size[], float_3 detector_transl[], float_3 detector_rotat[], float_3 source_pos[], u_int_3 volume_voxels, float_3 volume_size, float t_step, int interpolation, int use_gpu);

int tt_project_ray(VolumeType h_volume[], u_int_3 volume_voxels, float out_projections[], u_int n_projections, float_2 detector_scale[], float_3 detector_transl[], float_3 detector_rotat[], u_int_2 detector_pixels, float_3 source_position[], float_3 volume_size, float t_step, int use_gpu);


