%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

example =  'fourcenterAnalytic'
%example =  'Cylinder2D';
%example =  'ocean';

R = 0;
if strcmp(example, 'ocean')
    % ocean dataset
    % physical time range of data: 733010..733010+91
    Rslice = 0; % range 0..13
    R = 733010 + (91/13) * Rslice;
elseif strcmp(example, 'fourcenterAnalytic')
    % physical time range of data: 0..2*pi
    Rslice = 0; % range 0..31
    R=(2*pi/31)*Rslice;
elseif strcmp(example, 'Cylinder2D')
    % physical time range of data: 0..8
    Rslice = 0; % range 0..1000
    R=(8/1000)*Rslice;
end

%% pick file with GUI
[file,path] = uigetfile('*.mat','open header file','header.mat');
filename = [path,file];

% expects a filename that will be used to construct 3 associated filenames.
% the headerfile file [filename, '.mat'] holds the information of the grid
% data.
% the file [filename, '_v.raw'] holds the binary data of the field v.
% the file [filename, '_u.raw'] holds the binary data of the field u.
filenameBase = filename;
[filepath, name, ext] = fileparts(filename);
if strcmp(ext,'.mat')
    % strip off extension
    filenameBase = [filepath, '/', name];
end
filename = [filenameBase,'.mat'];
filenameV = [filenameBase,'_v.raw'];
filenameU = [filenameBase,'_u.raw'];


  %% load data and init vector fields
domainV = LoadRawBinaryData(filename, filenameV, false);
domainU = LoadRawBinaryData(filename, filenameU, false);


resampleSizeX = domainV.w;
resampleSizeY = domainV.h;
takeTime = tic;
%% start resampling script
grid = resampleObservedField(domainV, domainU, R, resampleSizeX, resampleSizeY);

sec = toc(takeTime);

signature = [filenameV,'_r',num2str(R)];
	% set the z component of the velocity

fileID = fopen(strcat(signature,'.raw'),'wb');
% save the permutated (dimensions 2 and 3 are switched to match input file format) grid as binary raw
disp(['saving file ', strcat(signature,'.raw')])
fwrite(fileID,grid,domainV.datatype);
fclose(fileID);
xmin = domainV.xmin;
xmax = domainV.xmax;
ymin = domainV.ymin;
ymax = domainV.ymax;
zmin = domainV.zmin;
zmax = domainV.zmax;
tmin = domainV.tmin;
tmax = domainV.tmax;
spaceunit = domainV.spaceunit;
timeunit = domainV.timeunit;
xsize = domainV.w;
ysize = domainV.h;
zsize = 1;
tsize = domainV.timeSteps;
vectordimension = domainV.vectordimension;
datatype = domainV.datatype;
save([signature,'.mat'],'xmin','xmax','ymin','ymax','zmin','zmax','tmin','tmax','spaceunit','timeunit','xsize','ysize','zsize','tsize','vectordimension','datatype');


