

#include "_reg_tools.h"
#include "_reg_resampling.h"
#include "_reg_affineTransformation.h"
#include "_et_line_integral.h"
#include "_et_accumulate.h"
#include "_et_line_backproject.h"
#include "_et_clear_accumulator.h"
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
#include "_tt_perspective_gpu.h"
//#include "_et_joint_histogram_gpu.h"

int tt_project_gpu(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background);
#endif

int tt_project(nifti_image *attenuationImage, nifti_image *projectionImage, nifti_image *psfImage, int n_projections, float *image_origin, float *detector_origin, float *detector_shape, float background);

