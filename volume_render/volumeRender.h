
#include <GL/glew.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include <cutil_inline.h>
#include <cutil_gl_inline.h>

#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>
#include <cuda_gl_interop.h>

#include <volumeRender_kernel.h>

typedef unsigned int uint;
typedef unsigned char uchar;

extern "C" void run_gui(int volume_x, int volume_y, int volume_z, int camera_x, int camera_y);

extern "C" int WindowDump(bool single);
extern "C" int ScreenShot(void);


#define MAX_EPSILON_ERROR 5.00f
#define THRESHOLD         0.30f
