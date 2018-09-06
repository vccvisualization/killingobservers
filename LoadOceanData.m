%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function domain = LoadOceanData(filename)

dataV = load(filename);

domain.tmin = dataV.time(1);
domain.tmax = domain.tmin+91;

% domain sizes
domain.w = 390;  
domain.h = 210; 
domain.d = 1;
domain.timeSteps = 14;

% support other common variable names
domain.xsize = domain.w;
domain.ysize = domain.h;
domain.zsize = domain.d;
domain.tsize = domain.timeSteps;
    
% actual domain of ocean data set: 
% lon: -12..29 lat: -45..-19

% actual grid size: 
% w: 165, h: 105, t: 14

% paper parameters of Haller 2016: 
domain.xmin = -4;
domain.xmax = 9;
domain.ymin = -35;
domain.ymax = -28;
domain.zmin = 1;
domain.zmax = 1;

u_interp = griddedInterpolant({dataV.lon,dataV.lat,dataV.time},permute(dataV.UT,[2,1,3]),'linear','none');
v_interp = griddedInterpolant({dataV.lon,dataV.lat,dataV.time},permute(dataV.VT,[2,1,3]),'linear','none');

domain.u = @(x,y,t) u_interp(x,y,t);
domain.v = @(x,y,t) v_interp(x,y,t);

