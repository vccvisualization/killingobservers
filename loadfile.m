%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [v,xsize,ysize,zsize,tsize, xmax, xmin, ymax, ymin, zmax, zmin, tmax, tmin, spaceunit, timeunit, datatype, vectordimension] = loadfile(headerfilename, filename)
    load(headerfilename);
    fId = fopen(filename);
    dataset = reshape(fread(fId,datatype),[vectordimension,xsize,ysize,zsize,tsize]);
    fclose(fId);
    v = dataset(1:vectordimension,:,:,:,:);
    if (vectordimension < 3 ) 
        v(3,:,:,:,:)=0; 
    end
end