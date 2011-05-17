function sinogram = et_project(activity, cameras, attenuation, psf, use_gpu, background, background_attenuation)
%ET_PROJECT
%    Projection function for Emission Tomographic reconstruction
%
%Description
%    Function for projection of activity into detector space.
%
%    SINOGRAM = ET_PROJECT(ACTIVITY, CAMERAS, PSF, USE_GPU, BACKGROUND)
%
%    ATIVITY is a 2D or 3D matrix of activity.
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
%    Depth-Dependent Point Spread Function is applyed in the frequency domain.
%
%Reference
%    Pedemonte, Bousse, Erlandsson, Modat, Arridge, Hutton, Ourselin, 
%    "GPU Accelerated Rotation-Based Emission Tomography Reconstruction", NSS/MIC 2010
%
%Example
%   N = 128;
%   use_gpu = 1;
%   activity = ones(N,N,N);
%   attenuationi = zeros(N,N,N);
%   PSF = ones(3,3,N);
%   cameras = [0:pi/100:pi];
%   sinogram = et_project(activity,cameras,attenuation,PSF,use_gpu);
%
%See also
%   ET_BACKPROJECT, ET_MLEM_RECONSTRUCT,
%   ET_LIST_GPUS, ET_SET_GPU
%
% 
%Stefano Pedemonte
%Copyright 2009-2010 CMIC-UCL
%Gower Street, London, UK


if not(exist('attenuation'))
    attenuation = 0;
end

if not(exist('background_attenuation'))
    background_attenuation = 0;
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

sinogram = et_project_mex(activity, cameras, attenuation, psf, use_gpu, background, background_attenuation);

