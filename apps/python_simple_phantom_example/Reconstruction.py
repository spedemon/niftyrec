
from thread import start_new_thread
from time import sleep
from numpy import zeros, ones, single, linspace, int32, round, pi, asarray
from numpy.random import rand
try: 
    from scipy.signal import convolve
except ImportError:
    HAVE_SCIPY = False
else:
    HAVE_SCIPY = True
import gobject

try:
	from NiftyRec.NiftyRec import et_project as project
	from NiftyRec.NiftyRec import et_backproject as backproject
except:
    print 'Error importing NiftyRec !'
    HAS_NIFTYREC = Falsec
    project = None
    backproject = None
else:
    HAS_NIFTYREC = True

BG_ACTIVITY = 0.0
BG_ATTENUATION = 0.0
EPSILON = 0.001


def osmapem_step(subset_order, activity_old, sinogram, cameras, attenuation, psf, beta_osl, gradient_prior, use_gpu, background_activity=BG_ACTIVITY, background_attenuation=BG_ATTENUATION, epsilon=EPSILON, verbose=False):

    N1 = activity_old.shape[0]
    N2 = activity_old.shape[1]
    N3 = activity_old.shape[2]
    N_cameras = cameras.shape[1]
    N_cameras_sub = int32(round(N_cameras/subset_order))
    cameras_indexes = int32(round(0.5 + (N_cameras-0.5) * rand(N_cameras_sub))-1)
    #cameras_sub = cameras[:,cameras_indexes]

    cameras_sub = zeros((3,N_cameras_sub))
    cameras_sub[:,:] = cameras[:,cameras_indexes]
#    print cameras
#    cameras_sub=cameras 

    #compute sensitivity for the subset
    normalization = backproject(ones((N1,N2,N_cameras_sub),dtype=single), cameras_sub, psf, attenuation, use_gpu, background_activity, background_attenuation)
    normalization = (normalization - beta_osl * gradient_prior)
    normalization = normalization*(normalization>0) + epsilon

    #project and backproject
    proj = project(activity_old, cameras_sub, psf, attenuation, use_gpu, background_activity, background_attenuation)
    proj = proj*(proj>0) + epsilon
    
    update = backproject(sinogram[cameras_indexes,:,:].reshape(N1,N2,N_cameras_sub) / proj, cameras_sub, psf, attenuation, use_gpu, background_activity, background_attenuation)
    update = update*(update>0)

    if verbose:
    	print "proj max",proj.max()
    	print "proj min",proj.min()
    	print "update max",update.max()
    	print "update min",update.min()
    	print "normalization max",normalization.max()
    	print "normalization min",normalization.min()
    	
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
        self.attenuation = None
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
#        if (sinogram.shape[1] != self.volume_size[0]) or (sinogram.shape[2] != self.volume_size[1]):
#        if (sinogram.shape[0] != self.volume_size[0]) or (sinogram.shape[1] != self.volume_size[1]):
#            print "Sinogram size not consistent with activity estimate, need reinterpolation"
#            print "Sinogram size: ",sinogram.shape
#            print "Activity size: ",self.volume_size
#            return
        self.sinogram = single(sinogram)
	return

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


    def reconstruct(self, method, parameters,verbose=False):
        if self.is_reconstructing():
            print "Reconstruction running"
            return
        print "Reconstruction started"
        if not self.has_all_parameters(): 
            print "Reconstructor does not have all necessary parameters"
            return 
        self._stop=False
        start_new_thread(self._reconstruct,(method,parameters,verbose))
#        self._reconstruct(method,parameters)

    def _reconstruct(self, method, parameters,verbose=False):
        if method=="osem":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
                gobject.idle_add(self.callback_status,"%d"%((100.0*step)/steps)+"% OSEM reconstruction",(1.0*step)/steps)
                gobject.idle_add(self.callback_updateactivity, self.activity )
		if verbose:	
		    print '\n----------------------------------\nStep %i of %i' % (step+1,steps)
                    print 'activity:', self.activity.shape, self.activity.dtype
                    print 'sinogram:', self.sinogram.shape, self.sinogram.dtype
                    print 'cameras: ', self.cameras.shape, self.cameras.dtype
                    if self.psf!=None: 
                        print 'psf:     ', self.psf.shape, self.psf.dtype
                    else:
                        print 'psf: None'
                    if self.attenuation!=None: 
                        print 'attenuation:',self.attenuation.shape, self.attenuation.dtype
                    else:
                        print 'attenuation: None'
                    print 'use gpu: ',self.use_gpu
                beta_osl=0
                gradient_prior=0
                self.activity = osmapem_step(subset_order, self.activity, self.sinogram, self.cameras, self.attenuation,\
                			     self.psf, beta_osl, gradient_prior, self.use_gpu, self.background_activity, self.background_attenuation, self.epsilon, verbose)
                gobject.idle_add(self.callback_updateactivity, self.activity)
                if self._stop:
                    break
        if method=="tv":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
                gobject.idle_add(self.callback_status,"%d"%((100.0*step)/steps)+"% TV OSEM reconstruction",(1.0*step)/steps)
                gobject.idle_add(self.callback_updateactivity, self.activity )
		if verbose:
		    print '\n----------------------------------\nStep %i of %i' % (step+1,steps)
                    print 'activity:', self.activity.shape, self.activity.dtype
                    print 'sinogram:', self.sinogram.shape, self.sinogram.dtype
                    print 'cameras: ', self.cameras.shape, self.cameras.dtype
                    if self.psf!=None: 
                        print 'psf:     ', self.psf.shape, self.psf.dtype
                    else:
                        print 'psf: None'
                    if self.attenuation!=None: 
                        print 'attenuation:',self.attenuation.shape, self.attenuation.dtype
                    else:
                        print 'attenuation: None'
                    print 'use gpu: ',self.use_gpu
                beta_osl=parameters['beta']/1000.0
                kernel = -1*ones((3,3,3))
                kernel[1,1,1] = 26
                if HAVE_SCIPY:
                    gradient_prior=convolve(self.activity,kernel,'same')
                else:
                    print "\nPlease install Scipy! (cannot find scipy.signal)\n"
                    exit(-1)
                self.activity = osmapem_step(subset_order, self.activity, self.sinogram, self.cameras, self.attenuation,\
					     self.psf, beta_osl, gradient_prior, self.use_gpu, self.background_activity, self.background_attenuation, self.epsilon, verbose)
                gobject.idle_add(self.callback_updateactivity, self.activity)
                if self._stop:
                    break
        elif method=="je":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
		if verbose: print '\n----------------------------------\nStep %i of %i' % (step+1,steps)
                gobject.idle_add(self.callback_status,"OSEM reconstruction",(1.0*step)/steps)
                gobject.idle_add(self.callback_updateactivity, self.activity )
                sleep(0.1)
                if self._stop:
                    break
        elif method=="ccje":
            steps = parameters['steps']
            subset_order = parameters['subset_order']
            self._reconstructing = True
            for step in range(steps):
		if verbose: print '\n----------------------------------\nStep %i of %i' % (step+1,steps)
                gobject.idle_add(self.callback_status,"OSEM reconstruction",(1.0*step)/steps)
                gobject.idle_add(self.callback_updateactivity, self.activity )
                sleep(0.1)
                if self._stop:
                    break
        print "Reconstruction done"
        self._reconstructing = False
        sleep(0.1)
        self.callback_status(" ",0.0)       

    def stop(self):
        if self.is_reconstructing():
            self._stop=True

