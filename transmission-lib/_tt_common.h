
#ifndef _TT_COMMON_H
#define _TT_COMMON_H

#include <math.h>

struct mat44{                   /** 4x4 matrix struct **/
  float m[4][4] ;
};
typedef mat44 mat44;

struct mat33{                   /** 3x3 matrix struct **/
  float m[3][3] ;
};
typedef mat33 mat33;

extern "C" mat44 reg_mat44_mul(mat44 *A, mat44 *B);
extern "C" int create_rotation_matrix44(mat44 *transformationMatrix, float theta_x, float theta_y, float theta_z, float center_x, float center_y, float center_z);

#endif
