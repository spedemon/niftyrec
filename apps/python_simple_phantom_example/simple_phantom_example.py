#!/usr/bin/env python
#
# NiftyRec Python Example Program
#
# D. Duke
# Energy Systems Division
# Argonne National Laboratory, Illinois, USA
#
# Supply a phantom via an 8-bit image, project it to a sinogram and then try and reconstruct it using et_project and et_backproject.
# Demonstrates use of the python interface to NiftyRec.
#
# NOTE: This program requires PIL (Python Imaging Library) to load a TIFF file and display it.
#
# Last Update: March 5, 2013.
#
##########################################################################################################################################

# Modules from your Python distribution
import os, time, sys
try:
	import numpy
except ImportError:
	print '\nPlease install the NumPy module!'
	os._exit(1)	

# NiftyRec Modules
try:
	from NiftyRec.NiftyRec import et_project as project
	from NiftyRec.NiftyRec import et_list_gpus as list_gpus
except ImportError:
	print '\nPlease install the Python NiftyRec module! It can be found in NiftyRec-x.x.x/Python/NiftyRec'
	os._exit(1)	

# Local Modules
try:
	from Reconstruction import Reconstructor
	from simple_phantom_image_handling import *
except ImportError:
	print '\nCan\t find Reconstruction.py and/or simple_phantom_image_handling.py. Please ensure they are in your cwd.'
	os._exit(1)	

##########################################################################################################################################
# Configure the settings below appropriately for your system.

verbose=True					# Show additional output while processing.
use_the_GPU=False				# Set to TRUE if you have a CUDA-enabled card and have build NiftyRec to use it.
input_phantom="phantom.tif"			# Source phantom image.
theta0=0. ; theta1=360.; n_cameras=100		# Number of equispaced virtual projections of the phantom
method='osem'  					# Methods available in Python are osem, tv, je, ccje.
params={'subset_order':8,'steps':30,'beta':0.}	# Parameters that methods may need in NiftyRec.
padFactor=1.0					# How much zero padding around the image as a fraction of image size. Must be >= 1.
resize_display=5				# Zoom the small image when showing it onscreen, so we can see it better.
psf_one=numpy.ones((1,1))			# Unity point spread function

### BEGIN MAIN ###########################################################################################################################

print '\nSimple Phantom Example for NiftyRec\n'

# Start a timer
start_time=time.time()

# Load image
print 'Loading phantom...'
i,N=image2array(input_phantom)		# Load phantom image (size N*n, N>n)
M=1 					# One image frame loaded (M = no. of slices)
i=numpy.true_divide(i,numpy.max(i)) 	# Normalize image
i=makeSquareImg(i)			# Square off the image to N*N
i,N=padImg(i,int(N*padFactor)) 		# Additionally pad image by some ratio of N.

# Scale-up image if number of camera angles is larger
if n_cameras>N:
	N=n_cameras
	i=rescaleImAr(i,(N,N))

# Create Reconstructor
volume_size = (N,M,N)
print 'The test phantom volume is of size',volume_size
r=Reconstructor(volume_size)

# Camera angles
r.set_cameras_equispaced(theta0,theta1,n_cameras, axis=1)  # r.cameras is in radians!

# If the number of camera angles is less than the image size, pad out with theta=0 instances, they won't add more information.
dummy_views_to_add=0
original_n_camera_angles=n_cameras
if N>n_cameras:
	dummy_views_to_add=N-n_cameras
	r.set_cameras_matrix(numpy.hstack((r.cameras,numpy.zeros((3,dummy_views_to_add)))))

# Apply PSF to every view
psf_mat=numpy.zeros((psf_one.shape[0],psf_one.shape[1],N))
for k in range(N): psf_mat[:,:,k]=psf_one
r.set_psf_matrix( psf_mat )
print 'The point spread function is of size',psf_mat.shape

# Check if CUDA is to be used , and show some diagnostic info.
r.set_use_gpu(use_the_GPU)
if use_the_GPU:	CPUGPU='GPU'
else:		CPUGPU='CPU'
print 'Using GPU?',r.use_gpu
if use_the_GPU: print 'GPU Info: ',list_gpus()

# Use NiftyRec to make the phantom's sinogram that we will then use to try and reconstruct the phantom.
print 'Making sinogram of phantom on the %s...' % CPUGPU
input_phantom_array=numpy.zeros(volume_size)				      	  # Make an empty volume
input_phantom_array[ 0:i.shape[0], 0, 0:i.shape[1] ]=i.astype(numpy.float32)  	  # Put image array into volume's first slice
NRsino=project(input_phantom_array, r.cameras, r.psf, r.attenuation, use_the_GPU) # Run et_project
print 'Size of NiftyRec\'s sinogram:',NRsino.shape

# If we had to add any 'dummy' views to pad the solution, zero out the sinogram for those views.
# et_project will obviously calculate the true projection at theta=0, and we want to remove that superfluous information.
if dummy_views_to_add>0: NRsino[original_n_camera_angles:,:,:]=-1

# Load sinogram and set callbacks
r.set_sinogram(NRsino)
r.set_callback_status(callback_status_handler)
r.set_callback_updateactivity(callback_updateactivity_handler)

# Check if we're ready to go, else quit.
print 'All paremeters set?',r.has_all_parameters(),'\n\n'
if not r.has_all_parameters(): os._exit(1)

# Begin reconstruction
sys.stdout.write( 'Reconstructing using `%s\' method on the %s...' % (method,CPUGPU) )
r.reconstruct(method,params,verbose)  # remove underscore to background-thread

# Waiting loop
time.sleep(1)			# r._reconstructing may not go high immediately.
while r._reconstructing:
	time.sleep(1)		# Don't hammer on the CPU while we wait!
	sys.stdout.flush()

# Stop the clock
elapsed_time = time.time() - start_time

# Make NumPy arrays of the first slice of the volume for the input phantom & reconstruction.
slice_input=input_phantom_array[:,0,:]
slice_output=r.activity[:,0,:]

# Determine sum square error in the reconstruction of the phantom
reconstruction_sum_sq_err = numpy.true_divide(numpy.sum((slice_input-slice_output)**2),numpy.sum(slice_input**2))
print 'Fraction error of reconstruction: %f\n' % reconstruction_sum_sq_err
print 'Time elapsed: %i seconds (on the %s)' % (numpy.ceil(elapsed_time),CPUGPU)

# Display sinogram, input & output images using PIL.
sinogram_slice_imageArray=numpy.transpose(NRsino.reshape(r.N_cameras,N))
displayImg(sinogram_slice_imageArray,resize=resize_display,title='NiftyRec sinogram')
displayImg(slice_output,resize=resize_display,title='NiftyRec reconstruction')
displayImg(slice_input,resize=resize_display,title='Original data')
