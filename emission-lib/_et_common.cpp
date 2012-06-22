/*
 *  _et_common.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_et_common.h"


extern "C" int fprintf_verbose (const char *__restrict __format, ...)
{
    va_list args;
    va_start(args,__format);  
    int return_status = 1;

    #ifdef _VERBOSE
    return_status = vfprintf(stderr, __format, args);
    #else
    return_status = 0;
    #endif

     va_end(args);
     return return_status;
 }





int et_create_rotation_matrix(mat44 *transformationMatrix, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z)
{
	int status = 0;
	float s_theta_x, c_theta_x, s_theta_y, c_theta_y, s_theta_z, c_theta_z;
	mat44 *rotation_x = (mat44 *)calloc(1,sizeof(mat44));
	mat44 *rotation_y = (mat44 *)calloc(1,sizeof(mat44));
	mat44 *rotation_z = (mat44 *)calloc(1,sizeof(mat44));
	
	// Initialize affine transform matrix
	s_theta_x = sin(theta_x);
	c_theta_x = cos(theta_x);
	s_theta_y  = sin(theta_y);
	c_theta_y  = cos(theta_y);
	s_theta_z = sin(theta_z);
	c_theta_z = cos(theta_z);
	 
	rotation_z->m[0][0] = c_theta_z;
	rotation_z->m[0][1] = -s_theta_z;
	rotation_z->m[0][2] = 0.0;
	rotation_z->m[0][3] = - center_x*c_theta_z + center_y*s_theta_z + center_x;
	rotation_z->m[1][0] = s_theta_z;
	rotation_z->m[1][1] = c_theta_z;
	rotation_z->m[1][2] = 0.0;
	rotation_z->m[1][3] = - center_x*s_theta_z - center_y*c_theta_z + center_y;
	rotation_z->m[2][0] = 0.0;
	rotation_z->m[2][1] = 0.0;
	rotation_z->m[2][2] = 1.0;
	rotation_z->m[2][3] = 0.0;
	rotation_z->m[3][0] = 0.0;
	rotation_z->m[3][1] = 0.0;
	rotation_z->m[3][2] = 0.0;
	rotation_z->m[3][3] = 1.0;	

	rotation_y->m[0][0] = c_theta_y;
	rotation_y->m[0][1] = 0.0;
	rotation_y->m[0][2] = s_theta_y;
	rotation_y->m[0][3] = - c_theta_y*center_x - s_theta_y*center_z + center_x;
	rotation_y->m[1][0] = 0.0;
	rotation_y->m[1][1] = 1.0;
	rotation_y->m[1][2] = 0.0;
	rotation_y->m[1][3] = 0.0;
	rotation_y->m[2][0] = -s_theta_y;
	rotation_y->m[2][1] = 0.0;
	rotation_y->m[2][2] = c_theta_y;
	rotation_y->m[2][3] = s_theta_y*center_x - c_theta_y*center_z + center_z;
	rotation_y->m[3][0] = 0.0;
	rotation_y->m[3][1] = 0.0;
	rotation_y->m[3][2] = 0.0;
	rotation_y->m[3][3] = 1.0;	
	
	rotation_x->m[0][0] = 1.0;
	rotation_x->m[0][1] = 0.0;
	rotation_x->m[0][2] = 0.0;
	rotation_x->m[0][3] = 0.0;
	rotation_x->m[1][0] = 0.0;
	rotation_x->m[1][1] = c_theta_x;
	rotation_x->m[1][2] = -s_theta_x;
	rotation_x->m[1][3] = -c_theta_x*center_y + s_theta_x*center_z + center_y;
	rotation_x->m[2][0] = 0.0;
	rotation_x->m[2][1] = s_theta_x;
	rotation_x->m[2][2] = c_theta_x;
	rotation_x->m[2][3] = -s_theta_x*center_y - c_theta_x*center_z + center_z;
	rotation_x->m[3][0] = 0.0;
	rotation_x->m[3][1] = 0.0;
	rotation_x->m[3][2] = 0.0;
	rotation_x->m[3][3] = 1.0;	
	
	*transformationMatrix = reg_mat44_mul(rotation_y, rotation_x);
	*transformationMatrix = reg_mat44_mul(rotation_z, transformationMatrix);
	
	free(rotation_x);
	free(rotation_y);
	free(rotation_z);

	return status;
}




extern "C" char *niftyrec_error_msg(int status)
{

    switch (status)
        {
        case niftyrec_success:
            sprintf(_niftyrec_msg,"niftyrec: OK.\n");
            break;
        case niftyrec_error_unspecified:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - unspecified error.\n");
            break;
        case niftyrec_error_parameters:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - inconsisten paramenters.\n");
            break;
        case niftyrec_error_kernel:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - gpu kernel execution error.\n");
            break;
        case niftyrec_error_allocgpu:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - GPU memory allocation.\n");
            break;
        case niftyrec_error_alloccpu:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - cpu memory allocation.\n");
            break;
        case niftyrec_error_transfergpu:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - transfering memory from guest to GPU or vice versa.\n");
            break;
        case niftyrec_error_nogpubuilt:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - niftyrec was built without GPU support. Please rebuild or install binary built with GPU support in order to use the GPU.\n");
            break;
        case niftyrec_error_nogpuinstalled:
            sprintf(_niftyrec_msg,"niftyrec: ERROR - no GPUs installed in the system.\n");
            break;
        default: 
            sprintf(_niftyrec_msg,"niftyrec: unknown status.\n");
            break;
        }
    return _niftyrec_msg;
}


