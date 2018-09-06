%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% raw binary file loading
% the raw file has the extension 'raw'.
% returns a domain-object containing the physical dimensions
% and gridded interpolant functions for the 'u' and 'v' component
function domain = LoadRawSyntheticData(vectordimension, xsize,ysize,zsize,tsize,filename)
    fId = fopen(filename);
    dataset = reshape(fread(fId,'float'),[vectordimension,xsize,ysize,zsize,tsize]);
    fclose(fId);
    domain.data = dataset(1:vectordimension,:,:,:,:);
    if (vectordimension == 2 ) 
        domain.data(3,:,:,:,:)=0; 
    end    
    
    domain.w = xsize;
    domain.h = ysize;
    domain.d = zsize;
    domain.timeSteps = tsize;
    domain.xmax = xsize-1;
    domain.xmin = 0;
    domain.ymax = ysize-1;
    domain.ymin = 0;
    domain.zmax = zsize-1;
    domain.zmin = 0;
    domain.tmax = tsize-1;
    domain.tmin = 0;
    domain.spaceunit = 1;
    domain.timeunit = 1;
    domain.datatype = 'float';
    domain.vectordimension =3;
    
    x = linspace(domain.xmin, domain.xmax, domain.w);        
    y = linspace(domain.ymin, domain.ymax, domain.h);     
    t = linspace(domain.tmin, domain.tmax, domain.timeSteps);     
    
    u_interp = griddedInterpolant({x',y',t'},permute(domain.data(1,:,:,:,:), [2,3,5,1,4]), 'linear','none');
    v_interp = griddedInterpolant({x',y',t'},permute(domain.data(2,:,:,:,:), [2,3,5,1,4]), 'linear','none');
    domain.u = @(x,y,t) u_interp(x,y,t);
    domain.v = @(x,y,t) v_interp(x,y,t);
end