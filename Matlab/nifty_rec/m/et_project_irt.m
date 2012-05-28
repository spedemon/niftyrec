function sinogram = et_project_irt(activity, cameras, attenuation, psf)
%ET_PROJECT_IRT
%    Projection function for Emission Tomographic reconstruction
%
%Description
%    It behaves like et_project, but uses IRT. 
%    Projects activity into detector space.
%
%    SINOGRAM = ET_PROJECT_IRT(ACTIVITY, CAMERAS, ATTENUATION, PSF)
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
%Reference
%    IRT reconstruction toolbox. Fessler
%
%Example
%   N = 128;
%   use_gpu = 1;
%   activity = ones(N,N,N);
%   attenuation = zeros(N,N,N);
%   PSF = ones(7,7,N);
%   cameras = [0:pi/100:pi]';
%   sinogram = et_project_irt(activity,cameras,attenuation,PSF);
%
%See also
%   ET_PROJECT, ET_BACKPROJECT_IRT, ET_MLEM_DEMO, ET_OSEM_DEMO
%
%
%Stefano Pedemonte
%Copyright 2009-2012 CMIC-UCL
%Gower Street, London, UK


if not(exist('attenuation'))
    attenuation = 0;
end

N_x = size(activity,1);
N_y = size(activity,3);
N_z = size(activity,2);
N_projections = size(cameras,1); 
if size(cameras,2)==3
    rotation_angle = 180/pi*(cameras(N_projections,2)-cameras(1,2)); 
    initial_rotation_angle = 180/pi*cameras(1,2);
else
    rotation_angle = 180/pi*(cameras(N_projections)-cameras(1)); 
    initial_rotation_angle = 180/pi*cameras(1);
end

activity2 = zeros(N_x,N_y,N_z);
for i=1:N_z
    activity2(:,:,i)=activity(:,i,:);
end

f.chat = 0;
f.option = {};
f.test = '3s';
dx=1;
ig = image_geom('nx', N_x, 'ny', N_y, 'nz', N_z, 'dx', dx, 'dz', dx);
f.psfs = '/tmp/t,psfs.fld';
f.blur = ',fft'; 
fld_write(f.psfs, psf)
f.na = N_projections;
f.dir = test_dir;
if attenuation==0
    f.mumap = '-';
else
    f.mumap = [f.dir 'mumap.fld'];
    fld_write(f.mumap, mumap); 
end
f.sfilter = 1; 
f.sys_type = '3s@%g,%g,%g,%g,%g,%d%s@%s@%s@-%d,%d,%d';
f.sys_type = sprintf(f.sys_type, ig.dx, ig.dx, ig.dz, rotation_angle,initial_rotation_angle, f.sfilter, f.blur, f.mumap, f.psfs, ig.nx, ig.nz, f.na); 
f.option = {f.option{:}, 'chat', f.chat, 'checkmask', 0};
G = Gtomo3(f.sys_type, ig.mask, ig.nx, ig.ny, ig.nz, f.option{:}, 'nthread', jf('ncore'));
sinogram = G*activity2;



