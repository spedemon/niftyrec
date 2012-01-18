
#ifndef _TT_BACKPROJECT_RAY_GPU_H
#define _TT_BACKPROJECT_RAY_GPU_H

#include <_tt_common.h>
// Utilities and System includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cutil_inline.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
#include <sys/time.h>
#include <cutil_math.h>

extern "C" int set_inViewMatrix(float *invViewMatrix, float_2 detector_scale, float_3 detector_transl, float_3 detector_rotat);
extern "C" int tt_backproject_cpu(float *out_backprojection, float *current_projection, float *invViewMatrix, uint2 detectorPixels, float3 sourcePosition, uint3 volumeVoxels, float3 volumeSize, float t_step, int interpolation);
extern "C" void tt_backproject_ray_kernel(dim3 gridSize, dim3 blockSize, float *d_projection, float *d_output, uint2 detectorPixels, float3 sourcePosition, uint3 volumeVoxels, float3 volumeSize, float t_step, int interpolation);
extern "C" void copyInvViewMatrix_bk(float *invViewMatrix, size_t sizeofMatrix);

//extern "C" int tt_backproject_ray_array(float h_projections[], u_int_2 detector_pixels, u_int n_projections, float out_backprojection[], float_2 detector_size[], float_3 detector_transl[], float_3 detector_rotat[], float_3 source_pos[], u_int_3 volume_voxels, float_3 volume_size, float t_step, int interpolation, int use_gpu);

#endif
