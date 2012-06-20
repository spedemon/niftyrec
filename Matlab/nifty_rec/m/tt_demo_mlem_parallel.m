

% TT_DEMO_MLEM_PARALLEL
%     NiftyRec Demo: iterative reconstruction transmission tomography with parallel rays 
%
%See also
%   TT_DEMO_PROJECT
%
% 
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


n_iter = 1000;
N = 128; 
N_cameras = 120; 
GPU = 1; 

load attenuation_128.mat
mask = et_spherical_phantom(N,N,N,N/2,1,0,(N+1)/2,(N+1)/2,(N+1)/2);
a = mask.*et_rotate(attenuation_128,[pi/2,0,0],[64.5,64.5,64.5],1,0)/100000;
cameras = zeros(N_cameras,3);
cameras(:,2)=(pi/180)*(0:360/N_cameras:360-360/N_cameras);

sino = et_project(zeros(N,N,N), cameras, a, 0, GPU);
min_invsino = min(1./sino(:)); 
max_invsino = max(1./sino(:)); 
image((-min_invsino+1./sino(:,:,1))/max_invsino*64); colormap gray; axis tight off equal; 
for i=1:120
    image((-min_invsino+1./sino(:,:,i))/max_invsino*64); colormap gray; axis tight off equal; pause(0.1)
end

B = et_backproject(sino, cameras, 0, 0, GPU); 
attenuation = 0.01*ones(N,N,N); 
for i =1:n_iter
    fprintf('iter %d \n',i);
    update = (et_backproject(et_project(zeros(N,N,N), cameras, attenuation, 0, GPU), cameras, 0, 0, GPU) + 0.0001) ./ (B + 0.0001);
    attenuation = mask.*(attenuation.*update);
    image(3000*attenuation(:,:,64)); axis equal tight off; pause(0.1); 
end
    
