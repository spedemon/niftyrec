function [sinogram,partial] = et_project_partial(activity, cameras, attenuation, psf, use_gpu, background, background_attenuation, truncate_negative_values)
%ET_PROJECT
%    Experimental - Projector for Emission Tomographic reconstruction. Returns
%    partial line integrals. 
%
%Description
%    Function for projection of activity into detector space.
%
%    [SINOGRAM,PARTIAL] = ET_PROJECT(ACTIVITY, CAMERAS, ATTENUATION, PSF, USE_GPU, BACKGROUND, BACKGROUND_ATTENUATION, TRUNCATE_NEGATIVE_VALUES)
%
%    SINOGRAM is the projection 
%
%    PARTIAL are the partial line intergrals. It's a volume of size
%    [Nx,Ny,Nz,Ncameras] where Nx,Ny,Nz is the size of the activity and
%    attenuation volumes, Ncameras is the number of camera positions. 
%    Each 3D volume [:,:,:,i] contains the partial line integrals for the
%    i-th camera. 
%
%    ATIVITY is a 2D or 3D matrix of activity.
%
%    CAMERAS specifies camera positions and it can be expressed in two forms: 
%    a matrix of size [n,3] representing angular position of each camera 
%    with respect of x,y,z axis; or a column vector of length n where for each 
%    camera, only rotation along z axis is specified. This is the most common 
%    situation in PET and SPECT where the gantry rotates along a single axis.
%
%    ATTENUATION specifies attenuation coefficients. It must be the same size as ACTIVITY.
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
%    TRUNCATE_NEGATIVE_VALUES defaults to 1. If 1 then negative results are truncated to 0. 
%    If 0 there is no truncation. 
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
%   N = 128;
%   use_gpu = 1;
%   activity = ones(N,N,N);
%   attenuation = zeros(N,N,N);
%   PSF = ones(7,7,N);
%   cameras = [0:pi/100:pi]';
%   sinogram = et_project(activity,cameras,attenuation,PSF,use_gpu);
%
%See also
%   ET_BACKPROJECT, ET_MAPEM_STEP, ET_MLEM_DEMO
%   ET_LIST_GPUS, ET_SET_GPU
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
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

if not(exist('truncate_negative_values'))
    truncate_negative_values = 1;
end

[sinogram,partial] = et_project_partial_mex(activity, cameras, attenuation, psf, use_gpu, background, background_attenuation, truncate_negative_values);
