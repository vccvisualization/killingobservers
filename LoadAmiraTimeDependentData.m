%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function expects a filename of a '.am' file for dim=2 and a path with
% '.am' files for each timestep for dim=3
% returns a domain-object containing the physical dimensions
% the data on the grid 
% and gridded interpolant functions for the 'x' and 'y' component for 2D
% grids and for 'x', 'y', and 'z' components for 3D grids.
% for dim=2 the filename needs to be a '*.am' file with time being the
% z-axis. for dim=3 the filename is a path that contains one file per
% timestep alphabetically ordered. 
% tmin and tmax need to be specified for dim=3 only
function domain = LoadAmiraTimeDependentData(filename, dim, tmin, tmax)
    if dim==2
        % we expect a time dependent field stored in a 3D volume
        % the x and y dimensions are preserved 
        % the z dimension is interpreted as time
        [domain.header, domain.data] = LoadAmiraFile(filename);
        domain.w = domain.header.xSize;
        domain.h = domain.header.ySize;
        domain.d = 1;
        
        domain.timeSteps = domain.header.zSize;
        % support other common variable names
        domain.xsize = domain.w;
        domain.ysize = domain.h;
        domain.zsize = domain.d;
        domain.tsize = domain.timeSteps;
        
        % physical dimensions
        domain.xmax = domain.header.BoundingBox.xMax;
        domain.xmin = domain.header.BoundingBox.xMin; 
        domain.ymax = domain.header.BoundingBox.yMax;
        domain.ymin = domain.header.BoundingBox.yMin;
        domain.zmax = 0;
        domain.zmin = 0;
        domain.tmax = domain.header.BoundingBox.zMax;
        domain.tmin = domain.header.BoundingBox.zMin;
        domain.spaceunit = 1; 
        domain.timeunit = 1;
        domain.datatype = 'float';
        domain.vectordimension = 3;

        x = linspace(domain.xmin, domain.xmax, domain.w);        
        y = linspace(domain.ymin, domain.ymax, domain.h);     
        t = linspace(domain.tmin, domain.tmax, domain.timeSteps);     

        u_interp = griddedInterpolant({x',y',t'},permute(domain.data(1,:,:,:), [2,3,4,1]), 'linear','none');
        v_interp = griddedInterpolant({x',y',t'},permute(domain.data(2,:,:,:), [2,3,4,1]), 'linear','none');
        domain.u = @(x,y,t) u_interp(x,y,t);
        domain.v = @(x,y,t) v_interp(x,y,t);
    elseif dim == 3
        % we expect a 3D time dependent volume
        % each timestep is saved in a separate file in the given directory
        directory = filename;
        files = dir(fullfile(directory,'*.am')); 
        % dummy element
        domain.data = zeros(4);
        domain.files = files;
        for k = 1:length(files)
            baseFileName = files(k).name;
            fullFileName = fullfile(directory, baseFileName);
            fprintf(1, 'Now reading %s\n', fullFileName);
            [domainTemp.header, domainTemp.data] = LoadAmiraFile(fullFileName);
            if k == 1
                % do the preallocation
                domain.data = zeros(3,domainTemp.header.xSize,domainTemp.header.ySize,domainTemp.header.zSize,length(files));
            end
            domain.data(:,:,:,:,k) = domainTemp.data;
        end
        domain.w = domainTemp.header.xSize;
        domain.h = domainTemp.header.ySize;
        domain.d = domainTemp.header.zSize;
        domain.timeSteps = length(files);
        domain.xmax = domainTemp.header.BoundingBox.xMax;
        domain.xmin = domainTemp.header.BoundingBox.xMin; 
        domain.ymax = domainTemp.header.BoundingBox.yMax;
        domain.ymin = domainTemp.header.BoundingBox.yMin;
        domain.zmax = domainTemp.header.BoundingBox.zMax;
        domain.zmin = domainTemp.header.BoundingBox.zMin;

        
        domain.tmax = tmax;
        domain.tmin = tmin;
        domain.spaceunit = 1; 
        domain.timeunit = 1;
        domain.datatype = 'float';
        domain.vectordimension = 3;

        x = linspace(domain.xmin, domain.xmax, domain.w);        
        y = linspace(domain.ymin, domain.ymax, domain.h);     
        z = linspace(domain.zmin, domain.zmax, domain.d);     
        t = linspace(domain.tmin, domain.tmax, domain.timeSteps);     

        u_interp = griddedInterpolant({x',y',z',t'},permute(domain.data(1,:,:,:,:), [2,3,4,5,1]), 'linear','none');
        v_interp = griddedInterpolant({x',y',z',t'},permute(domain.data(2,:,:,:,:), [2,3,4,5,1]), 'linear','none');
        w_interp = griddedInterpolant({x',y',z',t'},permute(domain.data(3,:,:,:,:), [2,3,4,5,1]), 'linear','none');
        domain.x = @(x,y,z,t) u_interp(x,y,z,t);
        domain.y = @(x,y,z,t) v_interp(x,y,z,t);
        domain.z = @(x,y,z,t) w_interp(x,y,z,t);
    end
end