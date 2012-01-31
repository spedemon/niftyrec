
% ET_OSEM_DEMO
%     NiftyRec Demo: MLEM SPECT reconstruction 
%
%See also
%   ET_MLEM_DEMO, ET_MAPEM_MRF_DEMO, ET_MAPEM_STEP
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N          = 128;
N_cameras  = 120;
cameras    = linspace(0,2*pi,N_cameras)';
psf        = ones(5,5,N);
N_counts   = 50e6;

iter_mlem    = 50;
subset_order = 32;
GPU          = 1;

%% Simulate SPECT scan 
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
phantom = et_spherical_phantom(N,N,N,N/8,30,10,N/4,N/3,N/2) .* mask;
attenuation = 0;
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram = poissrnd(ideal_sinogram);


%% Reconstruction:

disp('Reconstructing..');
activity = ones(N,N,N);
for i=1:iter_mlem
    fprintf('\nMLEM step: %d',i);
    activity = et_osmapem_step(subset_order, activity, sinogram, cameras, attenuation, psf, 0, 0, GPU, 0, 0.0001);
    imagesc(activity(:,:,floor(N/2))); colormap gray; axis square; pause(0.2)
end
disp('Done');

