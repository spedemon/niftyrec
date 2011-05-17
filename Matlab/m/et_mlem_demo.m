
%ET_MLEM
%Demo MLEM SPECT reconstruction 


%% Parameters
cameras    = [0:360/120:360-360/120]'*pi/180;
psf        = ones(15,15,N);
N          = 128;

iter_mlem = 30;
GPU        = 1;

disp('Loading sinogram..');
%sinogram   = load_nii('TOMO_I123_EM001_DS.img');
%sinogram   = double(sinogram.img);
%for i=1:120
%    sinogram(:,:,i) = sinogram(:,:,i)';
%end

sinogram = poissrnd( abs(et_project(et_spherical_phantom(N,N,N,N/8,100,0,N/4,N/3,N/2), cameras, psf, GPU)));

%% Reconstruction:

%Create normalization volume
disp('Creating normalization..');
norm = et_backproject(ones(N,N,length(cameras)), cameras, psf, GPU) ;
activity = et_reorder_activity_in(ones(N,N,N));

%Reconstruction
disp('Reconstructing..');
for i=1:iter_mlem
    fprintf('\nMLEM step: %d',i);
    activity = et_mapem_step(activity, norm, sinogram, cameras, psf, 0, 0, GPU, 0, 0.0001);
    a = et_reorder_activity_out(activity);
    imagesc(a(:,:,floor(N/4))); colormap gray; axis square; pause(0.2)
end
activity = et_reorder_activity_out(activity);

disp('Done');

