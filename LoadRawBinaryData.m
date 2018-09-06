%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function expects a filename of a '.mat' file that describes 
% the contents of the raw binary file of the same name.
% the raw file has the extension 'raw'.
% returns a domain-object containing the physical dimensions
% and gridded interpolant functions for the 'u' and 'v' component
function domain = LoadRawBinaryData(headerfilename, filename, adjustDomain)
    [domain.data,domain.w,domain.h,domain.d,domain.timeSteps, domain.xmax, domain.xmin, domain.ymax, domain.ymin, domain.zmax, domain.zmin, domain.tmax, domain.tmin, domain.spaceunit, domain.timeunit, domain.datatype, domain.vectordimension] = loadfile(headerfilename, filename);
    
    % support other common variable names
    domain.xsize = domain.w;
    domain.ysize = domain.h;
    domain.zsize = domain.d;
    domain.tsize = domain.timeSteps;
    
    if adjustDomain
        % shifting domain boundaries to make integration numerically more stable
        domainxspan = domain.xmax - domain.xmin;
        domainyspan = domain.ymax - domain.ymin;
        domainzspan = domain.zmax - domain.zmin;
    
        domain.xmin = domainxspan * -0.5;
        domain.xmax = domainxspan *  0.5;
        domain.ymin = domainyspan * -0.5;
        domain.ymax = domainyspan *  0.5;
        domain.zmin = domainzspan * -0.5;
        domain.zmax = domainzspan *  0.5;
        domain.tmax = domain.tmax-domain.tmin;
        domain.tmin = 0;
    end

    
        
    x = linspace(domain.xmin, domain.xmax, domain.w);        
    y = linspace(domain.ymin, domain.ymax, domain.h);     
    t = linspace(domain.tmin, domain.tmax, domain.timeSteps);     
    %disp 'WARNING: extrapolation is enabled for debugging'
    u_interp = griddedInterpolant({x',y',t'},permute(domain.data(1,:,:,:,:), [2,3,5,1,4]), 'linear','linear');
    v_interp = griddedInterpolant({x',y',t'},permute(domain.data(2,:,:,:,:), [2,3,5,1,4]), 'linear','linear');
    domain.u = @(x,y,t) u_interp(x,y,t);
    domain.v = @(x,y,t) v_interp(x,y,t);
end