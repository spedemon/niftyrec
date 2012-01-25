function image = et_backproject(sinogram, cameras, attenuation, psf, use_gpu, background, background_attenuation)

%ET_BACKPROJECT
%    Backprojection function for Emission Tomographic reconstruction
%
%Description
%    Function for backprojection of activity from detector space.
%
%    IMAGE = ET_BACKPROJECT(SINOGRAM, CAMERAS, ATTENUATION, PSF, USE_GPU, BACKGROUND, BACKGROUND_ATTENUATION)
%
%    SINOGRAM is a 2D or 3D sinogram.
%
%    CAMERAS specifies camera positions and it can be expressed in two forms: 
%    a matrix of size [n,3] representing angular position of each camera 
%    with respect of x,y,z axis; or a column vector of length n where for each 
%    camera, only rotation along z axis is specified. This is the most common 
%    situation in PET and SPECT where the gantry rotates along a single axis.
%
%    ATTENUATION specifies attenuation coefficients. 
%
%    PSF is a Depth-Dependent Point Spread Function.
%
%    USE_GPU is optional and it enables GPU acceleration if a compatible GPU 
%    device is installed in the system. By default use_gpu is set to 0 (disabled).
%
%    BACKGROUND is the value the background is set to when performing rotation.
%    It defaults to 0.
%
%    BACKGROUND_ATTENUATION is the value the attenuation background is set 
%    to when performing rotation. It defaults to 0.
%
%GPU acceleration
%    If a CUDA compatible Grahpics Processing Unit (GPU) is installed, 
%    the projection algorithm can take advantage of it. Set use_gpu parameter
%    to 1 to enable GPU acceleration. If GPU acceleration is not available, 
%    the value of the parameter is uninfluential.
%
%Algorithm notes
%    Rotation based projection algorithm with trilinear interpolation.
%    FFT-based Depth-Dependent Point Spread Function.
%
%Reference
%    Pedemonte, Bousse, Erlandsson, Modat, Arridge, Hutton, Ourselin, 
%    "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%Example
%    N = 128;
%    n_cameras = 100;
%    use_gpu = 1;
%    background = 0;
%    background_attenuation=0;
%    sinogram = ones(N,N,n_cameras);
%    attenuation = zeros(N,N,N);
%    PSF = ones(11,11,N);
%    cameras = [0:pi/n_cameras:pi-pi/n_cameras]';
%    image = et_backproject(sinogram,cameras,attenuation,PSF,use_gpu,background,background_attenuation);
%
%See also
%    ET_PROJECT, ET_MAPEM_STEP, ET_MLEM_DEMO
%    ET_LIST_GPUS, ET_SET_GPU
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL.
%Gower Street, London, UK


if not(exist('attenuation'))
    attenuation = 0;
end

if not(exist('psf'))
    psf = 0;
end

if not(exist('use_gpu'))
    use_gpu = 0;
end

if not(exist('background'))
    background = 0;
end

if not(exist('background_attenuation'))
    background_attenuation = 0;
end

image = et_backproject_mex(sinogram, cameras, attenuation, psf, use_gpu, background, background_attenuation);

