
% ET_DEMO_OSEM
%     NiftyRec Demo: Ordered Subsets Expectation Maximisation (OSEM) SPECT 
%     reconstruction - 3D.
%
%See also
%   ET_DEMO_MLEM, ET_DEMO_MAPEM_MRF, ET_OSMAPEM_STEP
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

phantom_type = 0;  % 0 for 'brain FDG PET'; 1 for 'sphere in uniform background'

%% Simulate SPECT scan 
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
if phantom_type == 0
    phantom = et_load_nifti('activity_128.nii'); 
    phantom = phantom.img .* mask;
else
    phantom = et_spherical_phantom(N,N,N,N/8,30,10,N/4,N/3,N/2) .* mask;
end
attenuation = 0;
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram = et_poissrnd(ideal_sinogram);


%% Reconstruction:

disp('Reconstructing..');
hFig = figure(); set(hFig, 'Position', get(0,'ScreenSize')/2); 
activity = ones(N,N,N);
for i=1:iter_mlem
    fprintf('\nMLEM step: %d',i);
    activity = et_osmapem_step(subset_order, activity, sinogram, cameras, attenuation, psf, 0, 0, GPU, 0, 0.0001);
    scale = 400;
    subplot(1,3,1); image(scale*flipud(squeeze(mask(:,:,floor(N/2)).*activity(:,:,floor(N/2)))')); colormap gray; axis square off tight; 
    subplot(1,3,2); image(scale*flipud(squeeze(mask(:,floor(N/2),:).*activity(:,floor(N/2),:))')); colormap gray; axis square off tight; 
    subplot(1,3,3); image(scale*flipud(squeeze(mask(floor(N/2),:,:).*activity(floor(N/2),:,:))')); colormap gray; axis square off tight; pause(0.05);
end


%% Cleanup 

if GPU
    et_reset_gpu();
end

disp('Done');
