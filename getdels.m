%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute finite differences dv/dt,dx/dt,dy/dt,dz/dt
function [deldvdt,deldvdx,deldvdy,deldvdz] = getdels( v,n,vectordimension,xsize,ysize,zsize,tsize,dx,dy,dz,dt)

deldvdt = zeros(vectordimension,xsize,ysize,zsize,tsize);
deldvdx = zeros(vectordimension,xsize,ysize,zsize,tsize);
deldvdy = zeros(vectordimension,xsize,ysize,zsize,tsize);
deldvdz = zeros(vectordimension,xsize,ysize,zsize,tsize);

%  dv/dt
if (tsize>2)
    deldvdt(:,:,:,:,1) =  (v(:,:,:,:,2) - v(:,:,:,:,1))/dt; %forward difference for border 
    deldvdt(:,:,:,:,2:end-1) = (v(:,:,:,:,3:end) - v(:,:,:,:,1:end-2))/(2*dt);  %centreal difference
    deldvdt(:,:,:,:,end) = (v(:,:,:,:,end) - v(:,:,:,:,end-1))/dt; %forward difference for border 
end

% dv/dx
deldvdx(:,1,:,:,:) = (v(:,2,:,:,:) - v(:,1,:,:,:))/dx; %forward difference for border 
deldvdx(:,2:end-1,:,:,:) = (v(:,3:end,:,:,:) - v(:,1:end-2,:,:,:))/(2*dx); %centreal difference
deldvdx(:,end,:,:,:) = (v(:,end,:,:,:) - v(:,end-1,:,:,:))/dx; %forward difference for border 

% dv/dy
deldvdy(:,:,1,:,:) = (v(:,:,2,:,:) - v(:,:,1,:,:))/dy; %forward difference for border 
deldvdy(:,:,2:end-1,:,:) = (v(:,:,3:end,:,:) - v(:,:,1:end-2,:,:))/(2*dy); %centreal difference
deldvdy(:,:,end,:,:) = (v(:,:,end,:,:) - v(:,:,end-1,:,:))/dy; %forward difference for border 

% dv/dz
if (zsize>2)
    deldvdz(:,:,:,1,:) = (v(:,:,:,2,:) - v(:,:,:,1,:))/dz; %forward difference for border 
    deldvdz(:,:,:,2:end-1,:) = (v(:,:,:,3:end,:) - v(:,:,:,1:end-2,:))/(2*dz); %centreal difference
    deldvdz(:,:,:,end,:) = (v(:,:,:,end,:) - v(:,:,:,end-1,:))/dz; %forward difference for border 
end

deldvdx = reshape(deldvdx(:,:,:,:,:),[vectordimension,n]);
deldvdy = reshape(deldvdy(:,:,:,:,:),[vectordimension,n]);
deldvdz = reshape(deldvdz(:,:,:,:,:),[vectordimension,n]);

end

