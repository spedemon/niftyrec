# Image-handling and Callback Subroutines for NiftyRec Python Example Program
#
# D. Duke
# Energy Systems Division
# Argonne National Laboratory, Illinois, USA

import os

try:
	import numpy
except ImportError:
	print("\nPlease install the NumPy module!")
	os._exit(1)	

try:
	from PIL import Image,ImageDraw
except ImportError:
	print("\nPlease install the Python Imaging Library!")
	os._exit(1)

# Very simple callback handlers so we can see onscreen what happened
def callback_status_handler(*kwargs):
	print 'Status handler callback:', kwargs
def callback_updateactivity_handler(*kwargs):
	print("Updatea-ctivity handler callback:"+str(kwargs))

# Load monochrome image into numpy array using PIL
def image2array(file):
	Im=Image.open(file).convert('I')
	ImAr=numpy.array(Im.getdata(),'float').reshape(Im.size[1],Im.size[0])
	L=numpy.max(ImAr.shape)
	return ImAr,L

# Pad the image so it is square.
def makeSquareImg(im):
	if im.shape[0] != im.shape[1]:
		return im
	else:
		L=numpy.max(im.shape)
		newIm=numpy.zeros((L,L))
		newIm[0:im.shape[0],0:im.shape[1]]=im[:,:]
		return newIm

# Grow the image to some new size by padding with zeros - centering.
def padImg(imAr,newL):
	if newL == numpy.max(imAr.shape):
		return imAr, newL
	dx=(newL-imAr.shape[0])/2
	dy=(newL-imAr.shape[0])/2
	if (dx<=0)|(dy<=0):
		print('Image would be cropped!')
		return imAr, newL
	newIm=numpy.zeros((imAr.shape[0]+dx*2,imAr.shape[1]+dy*2))
	newIm[dx:dx+imAr.shape[0],dy:dy+imAr.shape[1]]=imAr
	return newIm, newL

# Use PIL to rescale an image array using nearest-neighbour method, 16 bit rasterizing.
def rescaleImAr(ImAr,newSz):
	m=numpy.min(ImAr)
	M=numpy.max(ImAr)
	AdjImAr=(ImAr-m)*65535./float(M-m) # Go to 16 bit int
	Im=Image.fromarray(AdjImAr.astype(numpy.uint16),'I;16').convert('I')
	Im=Im.resize(newSz,Image.NEAREST)
	OutImAr=numpy.array(Im.getdata(),'float').reshape(Im.size[1],Im.size[0])
	return (OutImAr*float(M-m) / 65535.)+m # Go back to floating point

# Display monochrome image array (currently set to do 8 bit downsampling, many image viewers do not like 16 bit images.)
def displayImg(imAr,manualContrast=(0,0),resize=1,title=''):
	# Adjust contrast
	if manualContrast==(0,0):
		m=numpy.min(imAr)
		M=numpy.max(imAr)
	else:
		m,M=manualContrast
	AdjImAr=(imAr-m)*255./float(M-m)

        # Convert from Numpy array to PIL Image
        h=imAr.shape[1]
        w=imAr.shape[0]
        Im=Image.fromarray(AdjImAr.astype(numpy.uint16),'I;16').convert('I')
        
	# Up scale
	if resize>1: Im=Im.resize((AdjImAr.shape[0]*resize,AdjImAr.shape[1]*resize),Image.NEAREST)

	# Title
	if title!='':
		draw=ImageDraw.Draw(Im)
		draw.text((5,5),title,fill=(255))

        print('displaying frame of size %i x %i' % (h,w))
        Im.show()
        return

# Renormalize an image array from 0->1 - image array should be float type first!
def normalizeImg(imAr):
	iMax=numpy.max(imAr)
	iMin=numpy.min(imAr)
	new_im=numpy.true_divide(imAr-iMin,iMax-iMin)
	return new_im
