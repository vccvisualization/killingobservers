%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function domain = getAnalytic2D(name)
if strcmp(name, 'empty')
    domain.xmin = -2;
    domain.xmax = 2;
    domain.ymin = -2;
    domain.ymax = 2;
    domain.zmin = 0;
    domain.zmax = 0;
    domain.tmin = 0;
    domain.tmax = 2*pi;

elseif strcmp(name, 'twocenter')
    domain.xmin = -2;
    domain.xmax = 2;
    domain.ymin = -2;
    domain.ymax = 2;
    domain.zmin = 0;
    domain.zmax = 0;
    domain.tmin = 0;
    domain.tmax = 2*pi;
    
    domain.x = @(x,y,t) 5/2 * (-y*(4*y*y+4*x*x+1));
    domain.y = @(x,y,t) 5/2 * ( x*(4*y*y+4*x*x-1));

elseif strcmp(name, 'fourcenter')
    domain.xmin = -2;
    domain.xmax = 2;
    domain.ymin = -2;
    domain.ymax = 2;
    domain.zmin = 0;
    domain.zmax = 0;
    domain.tmin = 0;
    domain.tmax = 2*pi;
    
    domain.x = @(x,y,t) -x*(2*y*y-1)*exp(-x*x-y*y);
    domain.y = @(x,y,t) y*(2*x*x-1)*exp(-x*x-y*y);
    
    
elseif strcmp(name, 'centerCW')
    domain.xmin = -2;
    domain.xmax = 2;
    domain.ymin = -2;
    domain.ymax = 2;
    domain.zmin = 0;
    domain.zmax = 0;
    domain.tmin = 0;
    domain.tmax = 2*pi;
    domain.x = @(x,y,t) -y;
    domain.y = @(x,y,t) x;
elseif strcmp(name, 'centerCCW')
    domain.xmin = -2;
    domain.xmax = 2;
    domain.ymin = -2;
    domain.ymax = 2;
    domain.zmin = 0;
    domain.zmax = 0;
    domain.tmin = 0;
    domain.tmax = 2*pi;
    domain.x = @(x,y,t) y;
    domain.y = @(x,y,t) -x;

elseif strcmp(name, 'CCW')
    domain.xmin = -2;
    domain.xmax = 2;
    domain.ymin = -2;
    domain.ymax = 2;
    domain.zmin = 0;
    domain.zmax = 0;
    domain.tmin = 0;
    domain.tmax = 2*pi;
   
    domain.x = @(x,y,t) x*cos(-t)+y*-sin(-t);
    domain.y = @(x,y,t) x*sin(-t)+y*cos(-t);
  
elseif strcmp(name, 'CW')
    domain.xmin = -2;
    domain.xmax = 2;
    domain.ymin = -2;
    domain.ymax = 2;
    domain.zmin = 0;
    domain.zmax = 0;
    domain.tmin = 0;
    domain.tmax = 2*pi;
   
    domain.x = @(x,y,t) x*cos(t)+y*-sin(t);
    domain.y = @(x,y,t) x*sin(t)+y*cos(t);
end
    
if ~strcmp(name, 'empty')
    domain.u = domain.x;
    domain.v = domain.y;
end
end