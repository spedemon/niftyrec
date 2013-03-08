
% ET_DEMO_MLEM
%     NiftyRec Demo: MLEM SPECT reconstruction 
%
%See also
%   ET_DEMO_OSEM, ET_DEMO_MAPEM_MRF, ET_MAPEM_STEP
%
% 
%Stefano Pedemonte
%Copyright 2009-2013 CMIC-UCL
%Gower Street, London, UK


%% Parameters
N          = 128;
N_cameras  = 120;
cameras    = linspace(0,2*pi,N_cameras)';
psf        = ones(5,5,N);
N_counts   = 50e6;

iter_mlem  = 30;
GPU        = 1;

phantom_type   = 3;  % 0 for 'brain FDG PET'; 1 for sphere in uniform background; 2 for spiral

%% Initialise plot
figure(1); subplot(2,3,1); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,2); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,3); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,4); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,5); image(zeros(N,N)); colormap gray; axis square off tight; 
figure(1); subplot(2,3,6); image(zeros(N,N)); colormap gray; axis square off tight; pause(0.1); 

%% Load or generate phantom 
disp('Creating synthetic sinogram..');
mask = et_spherical_phantom(N,N,N,N*0.45,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
if phantom_type == 0
    phantom = et_load_nifti('activity_128.nii'); 
    phantom = phantom.img .* mask;
elseif phantom_type == 1
    phantom = et_spherical_phantom(N,N,N,N/8,60,20,N/4,N/2,N/3) .* mask;
else
    phantom = et_spiral_phantom(N,N,N,60,20,8,2,10,40,2) .* mask;
end
attenuation = 0;
figure(1); subplot(2,3,1); imagesc(squeeze(phantom(:,floor(N/2),:))); colormap gray; axis square off tight;  pause(0.1); 

%% Simulate SPECT scan 
ideal_sinogram = et_project(phantom, cameras, attenuation, psf, GPU);
ideal_sinogram = ideal_sinogram/sum(ideal_sinogram(:))*N_counts;
sinogram = et_poissrnd(ideal_sinogram);

disp('Visualising sinogram..');
for i=1:N_projections  
    figure(1); subplot(2,3,2); image(0.3*squeeze(sinogram(:,:,i))'); colormap gray; axis square off tight;  pause(0.1);     
    if i<N/2
        figure(1); subplot(2,3,3); image(0.3*squeeze(sinogram(:,i,:))); colormap gray; axis square off tight;  pause(0.1); 
    end
end

%% Compute normalization volume
disp('Computing normalization..');
norm = et_backproject(ones(N,N,length(cameras)), cameras, attenuation, psf, GPU) ;

%% Reconstruction:
disp('Reconstructing..');
activity = ones(N,N,N);
for i=1:iter_mlem
    fprintf('\nMLEM step: %d',i);
    activity = et_mapem_step(activity, norm, sinogram, cameras, attenuation, psf, 0, 0, GPU, 0, 0.0001);
    scale = 400;
    figure(1); subplot(2,3,4); image(scale*squeeze(mask(:,floor(N/2),:).*activity(:,floor(N/2),:))); colormap gray; axis square off tight;     
    figure(1); subplot(2,3,5); image(scale*squeeze(mask(floor(N/2),:,:).*activity(floor(N/2),:,:))); colormap gray; axis square off tight;
    figure(1); subplot(2,3,6); image(scale*squeeze(mask(:,:,floor(N/2)).*activity(:,:,floor(N/2))')); colormap gray; axis square off tight; pause(0.05);     
end



%% Cleanup 
if GPU
    et_reset_gpu();
end

disp('MLEM Done');

