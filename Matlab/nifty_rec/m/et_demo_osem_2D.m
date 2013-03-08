
% ET_DEMO_OSEM_2D
%     NiftyRec Demo: Ordered Subsets Expectation Maximisation (OSEM) SPECT
%     reconstruction - 2D.
%
%See also
%   ET_DEMO_MLEM, ET_DEMO_MAPEM_MRF, ET_OSMAPEM_STEP
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N          = 128;
N_cameras  = 120;
cameras = zeros(N_cameras,3);
cameras(:,2)=(pi/180)*(0:180/N_cameras:180-180/N_cameras);
psf        = ones(5,5,N);
N_counts   = 150e6/N;

iter_osem    = 100;
subset_order = 8;
GPU          = 0;

%% Simulate SPECT scan 
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,1,N,N*0.5,1,0,(N+1)/2,1,(N+1)/2);
phantom = mask .* et_spiral_phantom(N,1,N); 
attenuation = 0;
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram = et_poissrnd(ideal_sinogram);
figure(1); imagesc(squeeze(sinogram)); colormap gray; axis equal tight; 

%% Reconstruction:
disp('Reconstructing..');
activity = ones(N,1,N);
for i=1:iter_osem
    fprintf('\nMLEM step: %d',i);
    activity = mask.*et_osmapem_step(subset_order, activity, sinogram, cameras, attenuation, psf, 0, 0, GPU, 0, 0.0001);
    figure(2); imagesc(squeeze(activity(:,1,:))); colormap gray; axis square; pause(0.2)
end
disp('Done');

if GPU
    et_reset_gpu();
end

