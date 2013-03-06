
from ctypes import cdll, CDLL, pydll, PyDLL, CFUNCTYPE
from ctypes.util import find_library
from ctypes import string_at, byref, c_int, c_long, c_size_t, c_char_p, c_double, c_void_p, c_float
from ctypes import Structure, pointer, cast, addressof, c_int, c_double, c_float
from ctypes import POINTER as P
from numpy import zeros, asarray, float32, int32, ones, single, linspace, round, pi
from thread import start_new_thread
from time import sleep
from numpy.random import rand
from scipy.signal import convolve

lib_name = "lib_et_array_interface"

MAX_DEVICES = 16
SIZE_INFO = 5
ET_BACKGROUND_ACTIVITY = 0.0
ET_BACKGROUND_ATTENUATION = 0.0
ET_EPSILON = 0.000001
NIFTYREC_WEBSITE = "http://niftyrec.sourceforge.net"


L = None
for extension in ['so','dylib','dll']:
    try:
        L = CDLL(lib_name + "." + extension)
    except OSError as e:
        pass
if not L:
    raise(ImportError("Cannot find NiftyRec libraries ("+lib_name+")."))
    print("Please make sure that: ")
    print(" 1) NiftyRec is installed. If not, download from "+NIFTYREC_WEBSITE+" and install.")
    print(" 2) The path to the NiftyRec libraries is included in the system libraries path: ")
    print("    - LINUX: type the following at the terminal before you launch python: ")
    print("               export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/pathtoniftyreclibraries/ ")
    print("               (where '/pathtoniftyreclibraries/' is typically '/usr/local/lib/' or '/home/yourusername/niftyrec/lib/')")
    print("    - MAC OSx: type the following at the terminal before you launch python: ")
    print("               export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/pathtoniftyreclibraries/ ")
    print("               (where '/pathtoniftyreclibraries/' is typically '/usr/local/lib/' or '/Users/yourusername/niftyrec/lib/')")
    print("    - WINDOWS: Open the control panel, click on classic view, system, advanced ")
    print("               1) click on 'Start'")
    print("               2) click on 'Control Panel' and switch to 'Classic View' on the left panel.  ")
    print("               3) Click on 'System' to open the 'system properties' window. ")
    print("               4) Select the 'Advanced' tab in 'System properties' and click on 'Environment Variables'. ")
    print("               5) Scroll down the list under 'System Variables' and find the 'Path' variable. ")
    print("               6) Click on 'Path' to highlight it and then on the 'Edit' button. A new window 'Edit System Variables' should appear. ")
    print("               7) Click on the field next to 'Variable Value' and using the right arrow '--->' on the keyboard, reach to the end of the string. Then type ';C:\pathtoniftyreclibraries\ ' ")
    print("               8) (where '/pathtoniftyreclibraries/' is typically '\Program Files\NiftyRec\NiftyRec\lib\ '.) ")
    print("               9) Click 'OK' then 'OK' again. Then click 'OK' once again to close the 'system properties' window.  ")
    print("               10) Relaunch Python ")



def et_project(activity, cameras, psf=None, attenuation=None, gpu=1, background=0.0, background_attenuation=0.0):
  """Project activity to multiple cameras"""
  if gpu and not (et_is_block_multiple(activity.shape[0]) and et_is_block_multiple(activity.shape[1]) and et_is_block_multiple(activity.shape[2])):
      print "Error et_project() - With GPU acceleration enabled, the volume size has to be multiple of ",et_get_block_size()
      return None

  gpu        = int32(gpu)
  background = float32(background)
  background_attenuation = float32(background_attenuation)

  cameras=cameras.transpose()

  #Convert all matrices to float32
  if not(activity.dtype==float32):
      activity = asarray(activity,dtype=float32)
  if attenuation != None:
      if not(attenuation.dtype==float32):
          attenuation = asarray(attenuation,dtype=float32)
  if psf != None:
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
  if psf != None:
      psf_p      = psf.ctypes.data_as(P(c_float))
      psf_size   = asarray(psf.shape,dtype=int32)
      psf_size_p = psf_size.ctypes.data_as(P(c_int))
  else:
      psf_p      = None
      psf_size   = asarray(zeros([1,3]),dtype=int32)
      psf_size_p = psf_size.ctypes.data_as(P(c_int))
  if attenuation != None:
      attenuation_p      = attenuation.ctypes.data_as(P(c_float))
      attenuation_size   = asarray(attenuation.shape,dtype=int32)
      attenuation_size_p = attenuation_size.ctypes.data_as(P(c_int))
  else:
      attenuation_p      = None
      attenuation_size   = asarray(zeros([1,3]),dtype=int32)
      attenuation_size_p = attenuation_size.ctypes.data_as(P(c_int))      
  projection_size   = [activity_size[0],activity_size[1],cameras_size[0]]
  projection        = zeros(projection_size,dtype=float32)
  projection_p      = projection.ctypes.data_as(P(c_float))
  projection_size   = asarray(projection.shape,dtype=int32)
  projection_size_p = projection_size.ctypes.data_as(P(c_int))
  
  status = L.et_array_project(activity_p, activity_size_p, projection_p, projection_size_p, cameras_p, cameras_size_p, psf_p, psf_size_p, attenuation_p, attenuation_size_p, background, background_attenuation, gpu)
  if status:
      projection = None
  return projection



def et_backproject(projections, cameras, psf=None, attenuation=None, gpu=1, background=0, background_attenuation=0.0):
  """Backproject from multiple cameras to body space"""
  if gpu and not (et_is_block_multiple(projections.shape[0]) and et_is_block_multiple(projections.shape[1])):
      print "Error et_backproject() - With GPU acceleration enabled, the volume size has to be multiple of ",et_get_block_size()
      return None

  gpu        = int32(gpu)
  background = float32(background)
  background_attenuation = float32(background_attenuation)

  cameras = cameras.transpose()  

  #Convert all matrices to float32
  if not(projections.dtype==float32):
      projections = asarray(projections,dtype=float32)
  if attenuation!=None:
      if not(attenuation.dtype==float32):
          attenuation = asarray(attenuation,dtype=float32)
  if psf!=None:
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
  if psf != None:
      psf_p      = psf.ctypes.data_as(P(c_float))
      psf_size   = asarray(psf.shape,dtype=int32)
      psf_size_p = psf_size.ctypes.data_as(P(c_int))
  else:
      psf_p      = None
      psf_size   = asarray(zeros([1,3]),dtype=int32)
      psf_size_p = psf_size.ctypes.data_as(P(c_int))
  if attenuation != None:
      attenuation_p      = attenuation.ctypes.data_as(P(c_float))
      attenuation_size   = asarray(attenuation.shape,dtype=int32)
      attenuation_size_p = attenuation_size.ctypes.data_as(P(c_int))
  else:
      attenuation_p      = None
      attenuation_size   = asarray(zeros([1,3]),dtype=int32)
      attenuation_size_p = attenuation_size.ctypes.data_as(P(c_int))     
  bkprojection_size   = [projections_size[0],projections_size[1],projections_size[0]]
  bkprojection        = zeros(bkprojection_size,dtype=float32)
  bkprojection_p      = bkprojection.ctypes.data_as(P(c_float))
  bkprojection_size   = asarray(bkprojection.shape,dtype=int32)
  bkprojection_size_p = bkprojection_size.ctypes.data_as(P(c_int))
  
  status = L.et_array_backproject(projections_p, projections_size_p, bkprojection_p, bkprojection_size_p, cameras_p, cameras_size_p, psf_p, psf_size_p, attenuation_p, attenuation_size_p, background, background_attenuation, gpu)
  if status:
      bkprojection = None
  return bkprojection
  pass



def et_set_gpu(gpu_id):
  """Set a GPU device with given ID. Use in conjunction with get_gpu_list()"""
  L.et_array_set_gpu.restype  = c_int
  L.et_array_set_gpu.argtypes = [c_int]
  return L.et_array_set_gpu(gpu_id)



def et_list_gpus():
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



def et_is_block_multiple(integer_value): 
    L.et_array_is_block_multiple.restype  = c_int
    L.et_array_is_block_multiple.argtypes = [c_int]
    integer_value = int32(integer_value)
    status = L.et_array_is_block_multiple(integer_value)
    return status 


def et_get_block_size():
    L.et_array_get_block_size.restype = c_int
    L.et_array_get_block_size.argtypes = []
    return L.et_array_get_block_size()


def et_osmapem_step(subset_order, activity_old, sinogram, cameras, attenuation, psf, beta_osl, gradient_prior, use_gpu, background_activity=ET_BACKGROUND_ACTIVITY, background_attenuation=ET_BACKGROUND_ATTENUATION, epsilon=ET_EPSILON):
    """Maximum A Posteriori Expectation Maximisation - One Step Late"""
    N1 = activity_old.shape[0]
    N2 = activity_old.shape[1]
    N3 = activity_old.shape[2]
    N_cameras = cameras.shape[1]
    N_cameras_sub = int32(round(N_cameras/subset_order))
    cameras_indexes = int32(round(0.5 + (N_cameras-0.5) * rand(N_cameras_sub))-1)
#    cameras_sub = cameras[:,cameras_indexes]
    cameras_sub = zeros((3,N_cameras_sub))
    cameras_sub[:,:] = cameras[:,cameras_indexes]
    #compute sensitivity for the subset
    normalization = backproject(ones((N1,N2,N_cameras_sub),dtype=single), cameras_sub, psf, attenuation, use_gpu, background_activity, background_attenuation)
    normalization = (normalization - beta_osl * gradient_prior)
    normalization = normalization*(normalization>0) + epsilon
    #project and backproject
    proj = project(activity_old, cameras_sub, psf, attenuation, use_gpu, background_activity, background_attenuation)
    proj = proj*(proj>0) + epsilon
    update = backproject(sinogram[cameras_indexes,:,:].reshape(N1,N2,N_cameras_sub) / proj, cameras_sub, psf, attenuation, use_gpu, background_activity, background_attenuation)
    update = update*(update>0)
    activity_new = activity_old * update;
    activity_new = activity_new / (normalization)
    return activity_new





class Reconstructor:
    def __init__(self,volume_size):
        self._reconstructing = False
        self.volume_size = volume_size
        self.clear_activity()
        self.N_cameras   = 0
        self.cameras     = None
        self.psf         = None
        self.attenuation = zeros(self.volume_size,dtype=single)
        self.sinogram    = None
        self.use_gpu     = True
        self.background_activity = single(BG_ACTIVITY)
        self.background_attenuation = single(BG_ATTENUATION)
        self.epsilon     = EPSILON
        
    def set_callback_status(self,function):
        self.callback_status = function

    def set_callback_updateactivity(self,function):
        self.callback_updateactivity = function



    def set_activity(self,activity):
        self.activity = single(activity)

    def clear_activity(self):
        self.activity = ones(self.volume_size,dtype=single)

    def set_attenuation(self,attenuation):
        self.attenuation = single(attenuation)

    def set_sinogram(self,sinogram):
        if (sinogram.shape[1] != self.volume_size[0]) or (sinogram.shape[2] != self.volume_size[1]):
#        if (sinogram.shape[0] != self.volume_size[0]) or (sinogram.shape[1] != self.volume_size[1]):
            print("Sinogram size not consistent with activity eatimate, need reinterpolation")
            print("Sinogram size: "+str(sinogram.shape))
            print("Activity size: "+str(self.volume_size))
            return
        self.sinogram = single(sinogram)

    def set_psf_matrix(self,psf):
        self.psf = single(psf)

    def set_psf_two_points(self,psf0, plane0, psf1, plane1):
        #FIXME
        self.psf = single(ones((3,3,self.volume_size[2])))

    def set_cameras_matrix(self,cameras):
        self.cameras = single(cameras)
        self.N_cameras = cameras.shape[1]

    def set_cameras_equispaced(self,first_camera, last_camera, N_cameras, axis=1):
        first_camera = first_camera*pi/180.0
        last_camera  = last_camera*pi/180.0
        cameras = zeros((3,N_cameras))
        cameras[axis,:] = linspace(first_camera, last_camera, N_cameras)
        self.N_cameras = N_cameras
        self.cameras = single(cameras)        

    def set_use_gpu(self,enable):
        self.use_gpu = int(enable)

    def set_background_activity(self,bg):
        self.background_activity = single(bg)

    def set_background_attenuation(self,bg):
        self.background_attenuation = single(bg)

    def set_epsilon(self, epsilon):
        self.epsilon = single(epsilon)



    def has_all_parameters(self):
        if self.N_cameras == 0:
            return False
        if self.cameras == None:
            return False
        if self.sinogram == None: 
            return False
        return True


    def is_reconstructing(self):
        return self._reconstructing


    def reconstruct(self, method, parameters):
        if self.is_reconstructing():
            print("Reconstruction running")
            return
        print("Reconstruction started")
        if not self.has_all_parameters(): 
            print("Reconstructor does not have all necessary parameters")
            return 
        self._stop=False
        start_new_thread(self._reconstruct,(method,parameters))
#        self._reconstruct(method,parameters)

    def _reconstruct(self, method, parameters):
        if method=="osem":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
                #gobject.idle_add(self.callback_status,"%d"%((100.0*step)/steps)+"% OSEM reconstruction",(1.0*step)/steps)
                #gobject.idle_add(self.callback_updateactivity, self.activity )
                print('activity:'+str(self.activity.shape)+str(self.activity.dtype))
                print('sinogram:'+str(self.sinogram.shape)+str(self.sinogram.dtype))
                print('cameras: '+str(self.cameras.shape)+str(self.cameras.dtype))
                if self.psf!=None:
                    print 'psf:     ', self.psf.shape, self.psf.dtype
                else: 
                    print 'psf:     None'
                if self.attenuation!=None:
                    print 'attenuation:',self.attenuation.shape, self.attenuation.dtype
                else: 
                    print 'attenuation:None'
                print('use gpu: ',str(self.use_gpu))
                beta_osl=0
                gradient_prior=0
                self.activity = osmapem_step(subset_order, self.activity, self.sinogram, self.cameras, self.attenuation, 
                self.psf, beta_osl, gradient_prior, self.use_gpu, self.background_activity, self.background_attenuation, self.epsilon)
                #gobject.idle_add(self.callback_updateactivity, self.activity)
                if self._stop:
                    break
        if method=="quadratic_MRF":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
                #gobject.idle_add(self.callback_status,"%d"%((100.0*step)/steps)+"% TV OSEM reconstruction",(1.0*step)/steps)
                #gobject.idle_add(self.callback_updateactivity, self.activity )
                print('activity:'+str(self.activity.shape)+str(self.activity.dtype))
                print('sinogram:'+str(self.sinogram.shape)+str(self.sinogram.dtype))
                print('cameras: '+str(self.cameras.shape)+str(self.cameras.dtype))
                if self.psf!=None:
                    print 'psf:     ', self.psf.shape, self.psf.dtype
                else: 
                    print 'psf:     None'
                if self.attenuation!=None:
                    print 'attenuation:',self.attenuation.shape, self.attenuation.dtype
                else: 
                    print 'attenuation:None'
                print('use gpu: ',str(self.use_gpu))
                beta_osl=parameters['beta']/1000.0
                kernel = -1*ones((3,3,3))
                kernel[1,1,1] = 26
                gradient_prior=convolve(self.activity,kernel,'same')
                self.activity = osmapem_step(subset_order, self.activity, self.sinogram, self.cameras, self.attenuation, 
                self.psf, beta_osl, gradient_prior, self.use_gpu, self.background_activity, self.background_attenuation, self.epsilon)
                #gobject.idle_add(self.callback_updateactivity, self.activity)
                if self._stop:
                    break
        elif method=="je":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
                #gobject.idle_add(self.callback_status,"OSEM reconstruction",(1.0*step)/steps)
                #gobject.idle_add(self.callback_updateactivity, self.activity )
                sleep(0.1)
                if self._stop:
                    break
        elif method=="ccje":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
                #gobject.idle_add(self.callback_status,"OSEM reconstruction",(1.0*step)/steps)
                #gobject.idle_add(self.callback_updateactivity, self.activity )
                sleep(0.1)
                if self._stop:
                    break
        print("Reconstruction done")
        self._reconstructing = False
        sleep(0.1)
        self.callback_status(" ",0.0)       

    def stop(self):
        if self.is_reconstructing():
            self._stop=True


