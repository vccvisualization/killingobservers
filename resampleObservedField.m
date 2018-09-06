%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function grid = resampleObservedField(domainVin, domainUin, R, newSizeX, newSizeY)
% resampleObservedField Resamples a 2D time-dependent vector field v for the observer field u at observation time R
%
% grid = resampleObservedField(domainVin, domainUin, R, newSizeX, newSizeY)
%
% domainVin contains v as a domain object
% domainUin contains u as a domain object
% newSizeX .. the x-size of the resampled grid
% newSizeY .. the y-size of the resampled grid
%
% grid .. the resulting resampled vector field
%
% a domain object obj must have the following members:
% obj.xmin .. minimum physical x 
% obj.xmax .. maximum physical x 
% obj.ymin .. minimum physical y 
% obj.ymax .. maximum physical y 
% obj.tmin .. minimum physical t 
% obj.tmax .. maximum physical t
% obj.xsize .. size of the grid in x 
% obj.ysize .. size of the grid in y
% obj.tsize .. size of the grid in t
% obj.u .. function handle @(x,y,t) that evaluates the first vector component of the field at position x,y,t
% obj.v .. function handle @(x,y,t) that evaluates the second vector component of the field at position x,y,t

    % vector fields given as domain objects
    global domainU; 
    global domainV;
    global gui;
      
    domainV = domainVin;
    domainU = domainUin;
    
    % sanity checking
    if R < domainV.tmin || R > domainV.tmax
        disp 'ERROR: parameter R is outside the timeframe of the input data'
        return;
    end
    
    if domainV.xmax ~= domainU.xmax || domainV.xmin ~= domainU.xmin || domainV.ymax ~= domainU.ymax || domainV.ymin ~= domainU.ymin || domainV.zmax ~= domainU.zmax || domainV.zmin ~= domainU.zmin || domainV.tmax ~= domainU.tmax || domainV.tmin ~= domainU.tmin
        disp 'ERROR: all dimensions of domainU and domainV must match!'
        return;
    end
    
    if domainV.xmax < domainV.xmin || domainV.ymax < domainV.ymin || domainV.zmax < domainV.zmin || domainV.tmax  < domainV.tmin
        disp 'ERROR: all minimum dimensions must be smaller or equal the maximum dimension'
        return;
    end

    % integration settings
    global pathLineSettings
    pathLineSettings.cols = 5;
    pathLineSettings.rows = 5;
    pathLineSettings.iterations = 2;
  
    gui.resampleRows = newSizeY;
    gui.resampleCols = newSizeX;
    
    grid = resampleObservedFieldVolume(R, domainV.tsize);
end

function grid = resampleObservedFieldVolume(r, steps)
    global domainV;
    global domainU;
    global gui;
    global pathLineSettings;
    
    localResampleRows = gui.resampleRows;
    localResampleCols = gui.resampleCols;
    localTmin = domainV.tmin;
    localTmax = domainV.tmax;
    localIterations = pathLineSettings.iterations;
    
    tic
    
    %span multiple threads for parallel resampling
    parfor indexT = 1:steps
        %localResult = zeros(3, resampleRows,resampleCols,1);
        t =  localTmin + (indexT-1)*(localTmax-localTmin)/(steps-1);
        %disp(['progress ',num2str(100*(indexT-1)/(steps-1)), '% (',num2str(indexT), ' of ',num2str(steps),' slices)']);
        localResult{indexT} = resampleObservedSlice(domainV, domainU, r, t, localResampleRows, localResampleCols, localIterations);
    end 
    
    disp 'parfor finished';
    toc
    tic
    
    %collect results into one array
    grid = zeros(3, localResampleRows,localResampleCols,1,steps);
    for indexT = 1:steps
        grid(1:2,:,:,:,indexT) = localResult{indexT};
    end

    grid(3,:,:,:,:) = 1;
    grid = permute(grid,[1,3,2,4,5]);
    disp 'copying finished';
    toc
end



%% compute timestep T for fixed time R
% R denotes the fixed observer time 
% T is the observed time slice
function grid = resampleObservedSlice(domainV, domainU, R, T, resampleRows, resampleCols, iterations)
    
    %% find lookup positions 
    lookupPosInT = RtoT(domainU, resampleCols, resampleRows, R, T, iterations);

    %% lookup v-u in slice T
    tslice = ones(resampleRows,resampleCols)*T;
    VminusU_x = arrayfun(@(x,y,t)(domainV.u(x,y,t)-domainU.u(x,y,t)), permute(lookupPosInT(1,:),[2,1]), permute(lookupPosInT(2,:),[2,1]), tslice(:)); 
    VminusU_y = arrayfun(@(x,y,t)(domainV.v(x,y,t)-domainU.v(x,y,t)), permute(lookupPosInT(1,:),[2,1]), permute(lookupPosInT(2,:),[2,1]), tslice(:));

    
    % interleave VminusU_x with VminusU_y
    VminusU(1,:) = VminusU_x;
    VminusU(2,:) = VminusU_y;

    %% transform (v-u)
    %% pre-compute integration lookup table (T to R)
    preIntegratedTtoR = TtoR(domainU, resampleCols, resampleRows, R, T, iterations);
    % define interpolation grid for pre-computed integration
    xcoords = linspace(domainV.xmin,domainV.xmax, resampleCols);
    ycoords = linspace(domainV.ymin,domainV.ymax, resampleRows);
    
    % get function handle for gridded interpolation
    preIntegratedTtoRX = griddedInterpolant({ycoords,xcoords},permute(preIntegratedTtoR(1,:,:),[2,3,1]),'linear', 'linear');
    preIntegratedTtoRY = griddedInterpolant({ycoords,xcoords},permute(preIntegratedTtoR(2,:,:),[2,3,1]),'linear', 'linear');

    % get FTLE positions (the eps depends on the dimensions of the input grids only not on the resample grid dimensions)
    epsilonX = (domainV.xmax-domainV.xmin)/ domainV.xsize;
    epsilonY = (domainV.ymax-domainV.ymin)/ domainV.ysize;
    
    epsX = ones(resampleRows,resampleCols)*epsilonX;
    epsY = ones(resampleRows,resampleCols)*epsilonY;

    posEastInT(1,:,:)  = permute(lookupPosInT(1,:,:),[2,3,1])+epsX;
    posEastInT(2,:,:)  = permute(lookupPosInT(2,:,:),[2,3,1]);
    posWestInT(1,:,:)  = permute(lookupPosInT(1,:,:),[2,3,1])-epsX;
    posWestInT(2,:,:)  = permute(lookupPosInT(2,:,:),[2,3,1]);
    posNorthInT(1,:,:) = permute(lookupPosInT(1,:,:),[2,3,1]);
    posNorthInT(2,:,:) = permute(lookupPosInT(2,:,:),[2,3,1])+epsY;
    posSouthInT(1,:,:) = permute(lookupPosInT(1,:,:),[2,3,1]);
    posSouthInT(2,:,:) = permute(lookupPosInT(2,:,:),[2,3,1])-epsY;

    % the gridded interpolant expects a y,x lookup (because y is mapped to rows and x to columns)
    posEastXinR = preIntegratedTtoRX(permute(posEastInT(2,:,:),[2,3,1]),permute(posEastInT(1,:,:),[2,3,1]));
    posEastYinR = preIntegratedTtoRY(permute(posEastInT(2,:,:),[2,3,1]),permute(posEastInT(1,:,:),[2,3,1]));
    posWestXinR = preIntegratedTtoRX(permute(posWestInT(2,:,:),[2,3,1]),permute(posWestInT(1,:,:),[2,3,1]));
    posWestYinR = preIntegratedTtoRY(permute(posWestInT(2,:,:),[2,3,1]),permute(posWestInT(1,:,:),[2,3,1]));
    posNorthXinR = preIntegratedTtoRX(permute(posNorthInT(2,:,:),[2,3,1]),permute(posNorthInT(1,:,:),[2,3,1]));
    posNorthYinR = preIntegratedTtoRY(permute(posNorthInT(2,:,:),[2,3,1]),permute(posNorthInT(1,:,:),[2,3,1]));
    posSouthXinR = preIntegratedTtoRX(permute(posSouthInT(2,:,:),[2,3,1]),permute(posSouthInT(1,:,:),[2,3,1]));
    posSouthYinR = preIntegratedTtoRY(permute(posSouthInT(2,:,:),[2,3,1]),permute(posSouthInT(1,:,:),[2,3,1]));

    e1x = (posEastXinR-posWestXinR) /   (2*epsilonX);
    e1y = (posEastYinR-posWestYinR) /   (2*epsilonY);

    e2x = (posNorthXinR-posSouthXinR) / (2*epsilonX);
    e2y = (posNorthYinR-posSouthYinR) / (2*epsilonY);

    % interleave e1x and e2x for upper row of transformation matrix 
    exArray(1,:) = e1x(:);
    exArray(2,:) = e2x(:);
    % interleave e1y and e2y for lower row of transformation matrix 
    eyArray(1,:) = e1y(:);
    eyArray(2,:) = e2y(:);

    % transform (v-u) at each point to the new local coordinate system 
    % pointwise multiplication of upper row (matrix multiplication upper row)
    vnew(1,:) = exArray(:) .* VminusU(:);
    % pointwise multiplication of lower row (matrix multiplication lower row)
    vnew(2,:) = eyArray(:) .* VminusU(:);

    % add even and odd columns to reduce the local 2x2 matrix to a 2x1 vector (per element)
    result = vnew(:,1:2:end) + vnew(:,2:2:end);
   
    % reshape to a 2 x h x w vector field
    grid = reshape(result(:), 2,resampleRows,resampleCols);
    
end


function posAtR = TtoR(domainU, w, h, R, T, iterations)
    posAtR = RtoT(domainU, w, h, T, R, iterations);
end

function posAtT = RtoT(domainU, w, h, R, T, iterations)
    [p_t, num_points, num_lines] = timeSlicedPathlines(domainU, w, h, R, T, iterations);
    posArray = [p_t(end, :, 1);p_t(end, :, 2)];
    posAtT = reshape(posArray,2,h,w);
end

function [sliceTx, sliceTy] = lookupV(domainV, w, h, resamplePosX, resamplePosY, t)
    tslice = ones(w,h)*t;
    sliceTx = arrayfun(domainV.u, resamplePosX(:), resamplePosY(:), tslice(:));
    sliceTx = reshape(sliceTx, h, w);
    sliceTy = arrayfun(domainV.v, resamplePosX(:), resamplePosY(:), tslice(:));
    sliceTy = reshape(sliceTy, h, w);
end 


function [p_t, num_points, num_lines, tspan] = timeSlicedPathlines(domain, w, h, t1, t2, iterations)
    pathLineX = linspace(domain.xmin, domain.xmax, w);  
    pathLineY = linspace(domain.ymin, domain.ymax, h);   
    [pathLineStartX,pathLineStartY] = meshgrid(pathLineX, pathLineY);
    if t1==t2
        % return the set of starting positions
        num_points = 1;
        num_lines = w*h;
        p_t = zeros(num_points,num_lines,2);
        p_t(1,:,1) = pathLineStartX(:);
        p_t(1,:,2) = pathLineStartY(:);
        return;
    end

    disp 'WARNING number of iterations for the integration is independent of actual distance'
    tspan = linspace(t1,t2,iterations);
    
    options = odeset('RelTol',1e-4,'AbsTol',1e-4); % ODE solver options
    [~,F] = ode45(@evaluateFlow,tspan,[pathLineStartX(:);pathLineStartY(:)],options,domain.u,domain.v);

    xp_t = F(:,1:end/2);
    yp_t = F(:,end/2+1:end);

    [num_points, num_lines] = size(xp_t);
    p_t = zeros(num_points,num_lines,2);
    p_t(:,:,1) = xp_t;
    p_t(:,:,2) = yp_t;
end


%% evaluateFlow function called by ode45
% input arguments:
%   t: time
%   y: a vector of initial conditions - x and y positions are concatenated
%	u,v: function handles to evaluate the two vector components
% output arguments:
%   dy: concatenated velocity vector
function dy = evaluateFlow( t, y, u, v )
    Np = numel(y)/2;
    dy = zeros(2*Np,1);
	dy(1:Np,1) = u( y(1:Np,1),y(Np+1:2*Np,1),t*ones(Np,1) );
	dy(Np+1:2*Np,1) = v( y(1:Np,1),y(Np+1:2*Np,1),t*ones(Np,1) );
end
