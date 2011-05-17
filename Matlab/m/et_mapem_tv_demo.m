
%ET_MAPEM_TV_DEMO
%Demo Maximum A Posteriori (MAPEM) SPECT reconstruction with Total Variation prior


%% Parameters
cameras    = [0:360/120:360-360/120]'*pi/180;
psf        = 0;
N          = 128;
beta       = 10;

iter_mapem = 200;
GPU        = 1;

disp('Loading sinogram..');
%sinogram   = load_nii('TOMO_I123_EM001_DS.img');
%sinogram   = double(sinogram.img);
%for i=1:120
%    sinogram(:,:,i) = sinogram(:,:,i)';
%end
sinogram = poissrnd(et_project(et_spherical_phantom(N,N,N,N/8,100,0.1), cameras, psf, GPU));

%% Reconstruction:

%Create normalization volume
disp('Creating normalization..');
norm = et_backproject(ones(N,N,length(cameras)), cameras, psf, GPU) ;
activity = ones(N,N,N);

%Create kernel for Total Variation gradient
k1 = [-1,-1,-1;-1,-1,-1;-1,-1,-1];
k2 = [-1,-1,-1;-1,26,-1;-1,-1,-1];
k3 = [-1,-1,-1;-1,-1,-1;-1,-1,-1];
kernel = zeros(3,3,3);
kernel(:,:,1) = k1; kernel(:,:,2) = k2; kernel(:,:,3) = k3;

%Reconstruction
disp('Reconstructing..');
for i=1:iter_mapem
    fprintf('\nMLEM step: %d',i);
    prior_gradient = - convn(activity,kernel,'same');
    activity = et_mapem_step(activity, norm, sinogram, cameras, psf, beta, prior_gradient, GPU);
    imagesc(activity(:,:,N/2)); colormap gray; axis square; pause(0.1)
end

disp('Done');

