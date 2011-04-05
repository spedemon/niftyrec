
from ctypes import cdll, CDLL, pydll, PyDLL, CFUNCTYPE
from ctypes.util import find_library
from ctypes import string_at, byref, c_int, c_long, c_size_t, c_char_p, c_double, c_void_p, c_float
from ctypes import Structure, pointer, cast, addressof, c_int, c_double, c_float
from ctypes import POINTER as P
from numpy import zeros, asarray, float32, int32

lib_name = "lib_et_array_interface"

MAX_DEVICES = 16
SIZE_INFO = 5

L = None
for ext in ['so','dylib','dll']:
    try:
        L = CDLL(lib_name + "." + ext)
    except OSError,e:
        pass
if not L:
    raise ImportError("Cannot find "+lib_name+" shared object.")



def project(activity, cameras, psf, attenuation, gpu=1, background=0.0, background_attenuation=0.0):
  """Project activity to multiple cameras"""
  gpu        = int32(gpu)
  background = float32(background)
  background_attenuation = float32(background_attenuation)

  cameras=cameras.transpose()

  #Convert all matrices to float32
  if not(activity.dtype==float32):
      activity = asarray(activity,dtype=float32)
  if not(attenuation.dtype==float32):
      attenuation = asarray(attenuation,dtype=float32)
  if not(psf.dtype==float32):
      psf = asarray(psf,dtype=float32)
  if not(cameras.dtype==float32):
      cameras = asarray(cameras,dtype=float32)

  L.et_array_project.restype  = c_int
  L.et_array_project.argtypes = [P(c_float), P(c_int), P(c_float), P(c_int), P(c_float), P(c_int), P(c_float), P(c_int), P(c_float), P(c_int), c_float, c_float, c_int]

  activity_p      = activity.ctypes.data_as(P(c_float))
  activity_size   = asarray(activity.shape,dtype=int32)
  activity_size_p = activity_size.ctypes.data_as(P(c_int))
  cameras_p      = cameras.ctypes.data_as(P(c_float))
  cameras_size   = asarray(cameras.shape,dtype=int32)
  cameras_size_p = cameras_size.ctypes.data_as(P(c_int))
  psf_p      = psf.ctypes.data_as(P(c_float))
  psf_size   = asarray(psf.shape,dtype=int32)
  psf_size_p = psf_size.ctypes.data_as(P(c_int))
  attenuation_p      = attenuation.ctypes.data_as(P(c_float))
  attenuation_size   = asarray(attenuation.shape,dtype=int32)
  attenuation_size_p = attenuation_size.ctypes.data_as(P(c_int))
  projection_size   = [activity_size[0],activity_size[2],cameras_size[0]] 
  projection        = zeros(projection_size,dtype=float32)
  projection_p      = projection.ctypes.data_as(P(c_float))
  projection_size   = asarray(projection.shape,dtype=int32)
  projection_size_p = projection_size.ctypes.data_as(P(c_int))

  status = L.et_array_project(activity_p, activity_size_p, projection_p, projection_size_p, cameras_p, cameras_size_p, psf_p, psf_size_p, attenuation_p, attenuation_size_p, background, background_attenuation, gpu)
  if status:
      projection = None
  return projection



def backproject(projections, cameras, psf, attenuation, gpu=1, background=0, background_attenuation=0.0):
  """Backproject from multiple cameras to body space"""
  gpu        = int32(gpu)
  background = float32(background)
  background_attenuation = float32(background_attenuation)

  cameras = cameras.transpose()  

  #Convert all matrices to float32
  if not(projections.dtype==float32):
      projections = asarray(projections,dtype=float32)
  if not(attenuation.dtype==float32):
      attenuation = asarray(attenuation,dtype=float32)
  if not(psf.dtype==float32):
      psf = asarray(psf,dtype=float32)
  if not(cameras.dtype==float32):
      cameras = asarray(cameras,dtype=float32)

  L.et_array_backproject.restype  = c_int
  L.et_array_backproject.argtypes = [P(c_float), P(c_int), P(c_float), P(c_int), P(c_float), P(c_int), P(c_float), P(c_int), P(c_float), P(c_int), c_float, c_float, c_int]
  
  projections_p      = projections.ctypes.data_as(P(c_float))
  projections_size   = asarray(projections.shape,dtype=int32)
  projections_size_p = projections_size.ctypes.data_as(P(c_int))
  cameras_p          = cameras.ctypes.data_as(P(c_float))
  cameras_size       = asarray(cameras.shape,dtype=int32)
  cameras_size_p     = cameras_size.ctypes.data_as(P(c_int))
  psf_p              = psf.ctypes.data_as(P(c_float))
  psf_size           = asarray(psf.shape,dtype=int32)
  psf_size_p         = psf_size.ctypes.data_as(P(c_int))
  attenuation_p      = attenuation.ctypes.data_as(P(c_float))
  attenuation_size   = asarray(attenuation.shape,dtype=int32)
  attenuation_size_p = attenuation_size.ctypes.data_as(P(c_int))
  bkprojection_size   = [projections_size[0],projections_size[0],projections_size[1]]
  bkprojection        = zeros(bkprojection_size,dtype=float32)
  bkprojection_p      = bkprojection.ctypes.data_as(P(c_float))
  bkprojection_size   = asarray(bkprojection.shape,dtype=int32)
  bkprojection_size_p = bkprojection_size.ctypes.data_as(P(c_int))

  status = L.et_array_backproject(projections_p, projections_size_p, bkprojection_p, bkprojection_size_p, cameras_p, cameras_size_p, psf_p, psf_size_p, attenuation_p, attenuation_size_p, background, background_attenuation, gpu)
  if status:
      bkprojection = None
  return bkprojection
  pass



def set_gpu(gpu_id):
  """Set a GPU device with given ID. Use in conjunction with get_gpu_list()"""
  L.et_array_set_gpu.restype  = c_int
  L.et_array_set_gpu.argtypes = [c_int]
  return L.et_array_set_gpu(gpu_id)



def list_gpus():
  """Get list of GPU devices"""
  L.et_array_list_gpus.restype  = c_int
  L.et_array_list_gpus.argtypes = [P(c_int),P(c_int)]
  gpu_count   = asarray(0,dtype=int32)
  gpu_count_p = gpu_count.ctypes.data_as(P(c_int))
  gpus_info    = zeros([1,MAX_DEVICES*SIZE_INFO],dtype=int32)
  gpus_info_p  = gpus_info.ctypes.data_as(P(c_int))
  status = L.et_array_list_gpus(gpu_count_p, gpus_info_p)
  if status:
      return None
  gpu_count = int(gpu_count)
  gpus_info = gpus_info[0,0:SIZE_INFO*gpu_count]
  gpus_info = gpus_info.tolist()
  return [gpu_count, gpus_info]



#def affine():
#  """3D affine transformation of 3D image on GPU"""
#  pass



#def rotate():
#  """3D rotation of 3D image on GPU"""
#  pass



#def convolve():
#  """2D convolution of 3D omages 'slice by slice' on GPU"""
#  pass



