%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean up
clear all;
close all;

%% parameters to be setup
%% example: data set with preset parameters 
example =  'fourcenterAnalytic'
%example =  'Cylinder2D';
%example =  'ocean';

%% load data
if strcmp(example, 'ocean')
    % observer field optimization parameters
    kWeight = 1; % weight of K matrix 
    lambda = 1;  % weight of D matrix 
    mu = 0.05;   % weight of N matrix 

    % optimization algorithm parameters
    iterations = 500; % number of iterations 
    usePreconditioner = true;
    tolerance = 1e-6; % optimization target
    
    domainV = LoadOceanData('./data/ocean/Ocean_geostrophic_velocity.mat');  
    resampledDomain = getResampledDomain2D(domainV, domainV.xsize, domainV.ysize, domainV.tsize);    
    v = resampledDomain.grid;
    v(isnan(v))=0;
    v(3,:,:,:,:)=0; 
elseif strcmp(example, 'fourcenterAnalytic')
    % observer field optimization parameters
    kWeight = 1; % weight of K matrix 
    lambda = 1;  % weight of D matrix 
    mu = 0.01;   % weight of N matrix 

    % optimization algorithm parameters
    iterations = 200; % number of iterations 
    usePreconditioner = true;
    tolerance = 1e-6; % optimization target
    
    % data specifics 
    % define domain
    xsize = 32;
    ysize = xsize;
    tsize = ysize;
    
    domainV = getAnalytic2D('fourcenter');
    domainU = getAnalytic2D('centerCW');
    pushForwardAnalytic= getAnalytic2D('CW');
    
    domainVMinusUinT.x = @(x,y,t) domainV.x(pushForwardAnalytic.x(x,y,t),pushForwardAnalytic.y(x,y,t)) - domainU.x(pushForwardAnalytic.x(x,y,t),pushForwardAnalytic.y(x,y,t));
    domainVMinusUinT.y = @(x,y,t) domainV.y(pushForwardAnalytic.x(x,y,t),pushForwardAnalytic.y(x,y,t)) - domainU.y(pushForwardAnalytic.x(x,y,t),pushForwardAnalytic.y(x,y,t));
    
    domainVMinusUin0.x = @(x,y,t) pushForwardAnalytic.x(domainVMinusUinT.x(x,y,t), domainVMinusUinT.y(x,y,t), -t);
    domainVMinusUin0.y = @(x,y,t) pushForwardAnalytic.y(domainVMinusUinT.x(x,y,t), domainVMinusUinT.y(x,y,t), -t);
    
    domainAnalytic = getAnalytic2D('empty');
    domainAnalytic.u = @(xgrid,ygrid,tgrid)arrayfun(domainVMinusUin0.x, xgrid, ygrid, tgrid);
    domainAnalytic.v = @(xgrid,ygrid,tgrid)arrayfun(domainVMinusUin0.y, xgrid, ygrid, tgrid);
    
    domainV = getResampledDomain2D(domainAnalytic, xsize,ysize,tsize);
    v = domainV.grid;
    v(3,:,:,:,:)=0; 
    
elseif strcmp(example, 'Cylinder2D')
    
    % observer field optimization parameters
    kWeight = 1; % weight of K matrix 
    lambda = 1;  % weight of D matrix 
    mu = 0.05;   % weight of N matrix 
    
    % optimization algorithm parameters
    %use these settings for a very low residual (very high computation times ~10hours)
    %iterations = 15000; % number of iterations 
    %tolerance = 1e-6; % optimization target
    
    %use these settings for a medium residual (high computation times ~5hours)
    iterations = 15000; % number of iterations 
    tolerance = 1e-4; % optimization target
    
    %use these settings for faster approximation with higher residual
    %iterations = 1500; % number of iterations 
    %tolerance = 1e-4; % optimization target
    
    usePreconditioner = true;
    
    % data specifics 
    domainV = LoadRawBinaryData('./data/Cylinder2D/Cylinder2D.mat', 'data/Cylinder2D/Cylinder2D.raw', true);
    v = zeros(3, domainV.xsize, domainV.ysize, domainV.zsize, domainV.tsize);
    v(:,:,:,:,:) = domainV.data;
end

% grid sizes
xsize = domainV.w;
ysize = domainV.h;
zsize = domainV.d;
tsize = domainV.timeSteps;

xmax= domainV.xmax;
xmin= domainV.xmin;
ymax= domainV.ymax;
ymin= domainV.ymin;
zmax= domainV.zmax;
zmin= domainV.zmin;
tmax= domainV.tmax;
tmin= domainV.tmin;
spaceunit= 1;
timeunit= 1;
vectordimension = 3;
datatype='float';
    
deltaT = (tmax-tmin)/tsize; 
deltaX = (xmax-xmin)/xsize;
deltaY = (xmax-xmin)/ysize;
deltaZ = (xmax-xmin)/zsize;

slicesize = xsize*ysize;
volumesize = slicesize*zsize;
n = volumesize * tsize;

%% K = 0.5*(del+delTranspose)
% K ... Killing "operator"
% K*U ... Killing energy 
ticK = tic;

% compute K 
Kopt = getK( n, xsize, zsize, slicesize, volumesize, deltaX, deltaY, deltaZ );
Kopt = kWeight * Kopt;

sec = toc(ticK);
disp([num2str(sec),'s: K'])

%% D = -dudt+delVU-delUV
% compute D 
ticD = tic;
[deldvdt,deldvdx,deldvdy,deldvdz] = getdels(v,n,3,xsize,ysize,zsize,tsize,deltaX, deltaY, deltaZ, deltaT);
Dopt = getD(v,deldvdx,deldvdy,deldvdz,n,xsize,zsize,tsize,slicesize,volumesize,deltaX,deltaY,deltaZ,deltaT);

sec = toc(ticD);
disp([num2str(sec),'s: D'])

%% compute N
ticN = tic;
Nopt = sparse(1:3*n,1:3*n,-1*ones(1,3*n),3*n,3*n);

sec = toc(ticN);
disp([num2str(sec),'s: N'])


%% AT A u = AT b
ticAb = tic;

% left side
Aopt = sparse([Kopt;lambda*Dopt;mu*Nopt]);

AT = transpose(Aopt);
A = AT*Aopt;


% right-hand side
b1 = zeros( 9*n, 1 );
b2 = -lambda*deldvdt(:);
b3 = -mu * v(:);

b = [b1;b2;b3];
b = AT*b;

sec = toc(ticAb);
disp([num2str(sec),'s: A'])

clearvars AT
spparms('spumoni',2)

%% optimization
ticOptimization = tic;

if ~usePreconditioner
    [u,flag,relres,iter,resvec]=cgs(A,b,tolerance,iterations);
else
    alpha = max(sum(abs(A),2)./diag(A))-2;
    L = ichol(A,struct('type','ict','droptol',1e-4,'diagcomp',alpha,'michol','on'));
    [u,flag,relres,iter,resvec]=pcg(A,b,tolerance,iterations,L,L');
end

sec = toc(ticOptimization);
disp(['residual (|A*u-b|): ', num2str(sum(abs(A*u-b)))])

disp([num2str(sec),'s: optimization'])

residual = sum(abs(A*u-b));

%% save result as raw file
mkdir(['./data/',example]);
signature = ['_iterations', num2str(iterations),'_lambda', num2str(lambda), '_mu', num2str(mu), '_residual' ,num2str(residual)];
% save one mat file
save(['./data/',example,'/',example,signature,'.mat'],'xmin','xmax','ymin','ymax','zmin','zmax','tmin','tmax','spaceunit','timeunit','xsize','ysize','zsize','tsize','vectordimension','datatype');
% save _u.raw and _v.raw files

if isfield(domainV,'files')
    N = xsize*ysize*zsize;
    for k = 1:length(domainV.files)
        fileName = domainV.files(k).name;
        fprintf(1, 'Now saving %s\n', fileName);
        fileID = fopen(strcat('./data/',example,'/',fileName,'_',signature,'_u','.raw'),'wb');     
        startIndex = 1 + N*3*(k-1);
        endIndex   = N*3*k ;
        fwrite(fileID,u(startIndex:endIndex),datatype);
        fclose(fileID);
    end
else
    fileID = fopen(strcat('./data/',example,'/',example,signature,'_u','.raw'),'wb');    
    if (vectordimension <3)
        uu = zeros(numel(u)*2/3,1);
        uu(1:2:end) = u(1:3:end);
        uu(2:2:end) = u(2:3:end);
        fwrite(fileID,uu,datatype);
    else
        fwrite(fileID,u,datatype);
    end
    fclose(fileID);
end

fileID = fopen(strcat('./data/',example,'/',example,signature,'_v','.raw'),'wb');    
if (vectordimension <3)
    vv = zeros(numel(v)*2/3,1);
    vv(1:2:end) = v(1:3:end);
    vv(2:2:end) = v(2:3:end);
    fwrite(fileID,vv,datatype);
else
    fwrite(fileID,v,datatype);
end
fclose(fileID);