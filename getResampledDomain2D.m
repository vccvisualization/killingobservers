%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% resample 2D domain on a grid of size w x h x timeSteps
function resampledDomain = getResampledDomain2D(domain, w, h, timeSteps)
    resampledDomain = domain;
    resampledDomain.w = w;
    resampledDomain.h = h;
    resampledDomain.d = 1;
    resampledDomain.timeSteps = timeSteps;
    resampledDomain.spaceunit = 1;
    resampledDomain.timeunit = 1;
    
    % support other common variable names
    resampledDomain.xsize = resampledDomain.w;
    resampledDomain.ysize = resampledDomain.h;
    resampledDomain.zsize = resampledDomain.d;
    resampledDomain.tsize = resampledDomain.timeSteps;
    
    x = linspace(domain.xmin, domain.xmax, w);        
    y = linspace(domain.ymin, domain.ymax, h);     

    [yi,xi] = meshgrid(y,x);

    for t = 1:timeSteps
        tPhysical = resampledDomain.tmin;
        if timeSteps > 1
            tPhysical = resampledDomain.tmin + (resampledDomain.tmax-resampledDomain.tmin) / (timeSteps-1) * (t-1);
        end
        timegrid(1:resampledDomain.w,1:resampledDomain.h) = tPhysical;
        resampledDomain.resampledU(:,:,1,t) = arrayfun(domain.u, xi, yi, timegrid);
        resampledDomain.resampledV(:,:,1,t) = arrayfun(domain.v, xi, yi, timegrid);
    end
    resampledDomain.grid(1,:) = resampledDomain.resampledU(:);
    resampledDomain.grid(2,:) = resampledDomain.resampledV(:);
    resampledDomain.grid = reshape(resampledDomain.grid, 2, w, h, 1, timeSteps);
end