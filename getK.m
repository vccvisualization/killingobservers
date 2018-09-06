%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute the E matrix -> del(u) + transpose(del(u))
function [ K ] = getK( n, xsize, zsize, slicesize, volumesize, dx, dy, dz)

if (zsize > 2)
    vecx = zeros(36*n,1);
    vecy = zeros(36*n,1);
    vecz = zeros(36*n,1);
else
    vecx = zeros(24*n,1);
    vecy = zeros(24*n,1);
    vecz = zeros(24*n,1);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  0.5 * del u
counter = 1;
% compute 0.5*d(u)/d(x)
for i=0:n-1 
     modulo = mod(i, xsize); 
     if (modulo == 0)   
            vecx(counter:counter+5) = [ 9*i+1, 9*i+1, 9*i+2, 9*i+2, 9*i+3, 9*i+3 ]; %border
            vecy(counter:counter+5) = [ 3*i+1, 3*i+4, 3*i+2, 3*i+5, 3*i+3, 3*i+6 ];
            vecz(counter:counter+5) = [ -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx ];
            counter=counter+6;
     else if (modulo == xsize-1)
            vecx(counter:counter+5) = [ 9*i+1, 9*i+1, 9*i+2, 9*i+2, 9*i+3, 9*i+3 ]; %border
            vecy(counter:counter+5) = [ 3*i-2, 3*i+1, 3*i-1, 3*i+2, 3*i, 3*i+3 ];
            vecz(counter:counter+5) = [ -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx ];
            counter=counter+6;
         else
            vecx(counter:counter+5) = [ 9*i+1, 9*i+1, 9*i+2, 9*i+2, 9*i+3, 9*i+3 ];
            vecy(counter:counter+5) = [ 3*i-2, 3*i+4, 3*i-1, 3*i+5, 3*i, 3*i+6 ];
            vecz(counter:counter+5) = [ -0.25/dx, 0.25/dx, -0.25/dx, 0.25/dx, -0.25/dx, 0.25/dx ];
            counter=counter+6;
         end
    end
end


% compute 0.5*d(u)/d(y)
for i=0:n-1 
   modulo = mod(i, slicesize); 
   if (modulo < xsize )    
        vecx(counter:counter+5) = [ 9*i+4, 9*i+4,  9*i+5,  9*i+5,  9*i+6,  9*i+6 ]; %border
        vecy(counter:counter+5) = [ 3*i+1, 3*(i+xsize)+1, 3*i+2, 3*(i+xsize)+2, 3*i+3, 3*(i+xsize)+3 ];
        vecz(counter:counter+5) = [ -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy ];
        counter=counter+6;
   else if ( modulo > slicesize - xsize - 1 )
           vecx(counter:counter+5) = [ 9*i+4, 9*i+4,  9*i+5,  9*i+5,  9*i+6,  9*i+6 ]; %border
           vecy(counter:counter+5) = [ 3*(i-xsize)+1, 3*i+1, 3*(i-xsize)+2, 3*i+2, 3*(i-xsize)+3, 3*i+3 ];
           vecz(counter:counter+5) = [ -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy];
           counter=counter+6;   
       else
            vecx(counter:counter+5) = [ 9*i+4, 9*i+4,  9*i+5,  9*i+5,  9*i+6,  9*i+6 ]; 
            vecy(counter:counter+5) = [ 3*(i-xsize)+1, 3*(i+xsize)+1, 3*(i-xsize)+2, 3*(i+xsize)+2, 3*(i-xsize)+3, 3*(i+xsize)+3 ];
            vecz(counter:counter+5) = [ -0.25/dy, 0.25/dy, -0.25/dy, 0.25/dy, -0.25/dy, 0.25/dy ];
            counter=counter+6;           
       end
   end
end

% compute 0.5*d(u)/d(z)
if (zsize > 2)  %%compute z gradient only for 3D data
    for i=0:n-1 
        modulo = mod(i, volumesize); 
        if ( modulo < slicesize )
            vecx(counter:counter+5) = [ 9*i+7, 9*i+7,  9*i+8,  9*i+8,  9*i+9,  9*i+9 ]; %border
            vecy(counter:counter+5) = [ 3*i+1, 3*(i+slicesize)+1, 3*i+2, 3*(i+slicesize)+2, 3*i+3, 3*(i+slicesize)+3 ];
            vecz(counter:counter+5) = [ -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz];
            counter=counter+6;        
        else if (modulo > volumesize - slicesize - 1 )
            vecx(counter:counter+5) = [ 9*i+7, 9*i+7,  9*i+8,  9*i+8,  9*i+9,  9*i+9 ]; %border
            vecy(counter:counter+5) = [ 3*(i-slicesize)+1, 3*i+1, 3*(i-slicesize)+2, 3*i+2, 3*(i-slicesize)+3, 3*i+3 ];
            vecz(counter:counter+5) = [ -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz ];
            counter=counter+6;         
        else 
            vecx(counter:counter+5) = [ 9*i+7, 9*i+7,  9*i+8,  9*i+8,  9*i+9,  9*i+9 ];
            vecy(counter:counter+5) = [ 3*(i-slicesize)+1, 3*(i+slicesize)+1, 3*(i-slicesize)+2, 3*(i+slicesize)+2, 3*(i-slicesize)+3, 3*(i+slicesize)+3 ];
            vecz(counter:counter+5) = [ -0.25/dz, 0.25/dz, -0.25/dz, 0.25/dz, -0.25/dz, 0.25/dz ];
            counter=counter+6;
            end
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ADD SYMETRIC PART transposed 0.5 * (del u)

% compute 0.5*d(u)/d(x)

for i=0:n-1
     modulo = mod(i, xsize); 
     if (modulo == 0)   
            vecx(counter:counter+5) = [ 9*i+1, 9*i+1, 9*i+4, 9*i+4, 9*i+7, 9*i+7 ]; %border
            vecy(counter:counter+5) = [ 3*i+1, 3*i+4, 3*i+2, 3*i+5, 3*i+3, 3*i+6 ];
            vecz(counter:counter+5) = [ -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx];
            counter=counter+6;
     else if (modulo == xsize-1)
            vecx(counter:counter+5) = [ 9*i+1, 9*i+1, 9*i+4, 9*i+4, 9*i+7, 9*i+7 ]; %border
            vecy(counter:counter+5) = [ 3*i-2, 3*i+1, 3*i-1, 3*i+2, 3*i, 3*i+3 ];
            vecz(counter:counter+5) = [ -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx ];
            counter=counter+6;
         else
            vecx(counter:counter+5) = [ 9*i+1, 9*i+1, 9*i+4, 9*i+4, 9*i+7, 9*i+7 ];
            vecy(counter:counter+5) = [ 3*i-2, 3*i+4, 3*i-1, 3*i+5, 3*i, 3*i+6 ];
            vecz(counter:counter+5) = [ -0.25/dx, 0.25/dx, -0.25/dx, 0.25/dx, -0.25/dx, 0.25/dx ];
            counter=counter+6;
         end
    end
end

% compute 0.5*d(u)/d(y)

for i=0:n-1 
   modulo = mod(i, slicesize); 
   if (modulo < xsize )    
        vecx(counter:counter+5) = [ 9*i+2, 9*i+2,  9*i+5,  9*i+5,  9*i+8,  9*i+8 ]; %border
        vecy(counter:counter+5) = [ 3*i+1, 3*(i+xsize)+1, 3*i+2, 3*(i+xsize)+2, 3*i+3, 3*(i+xsize)+3 ];
        vecz(counter:counter+5) = [ -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy ];
        counter=counter+6;
   else if ( modulo > slicesize - xsize - 1 )
           vecx(counter:counter+5) = [ 9*i+2 9*i+2,  9*i+5,  9*i+5,  9*i+8,  9*i+8 ]; %border
           vecy(counter:counter+5) = [ 3*(i-xsize)+1, 3*i+1, 3*(i-xsize)+2, 3*i+2, 3*(i-xsize)+3, 3*i+3 ];
           vecz(counter:counter+5) = [ -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy];
           counter=counter+6;   
       else
            vecx(counter:counter+5) = [ 9*i+2, 9*i+2,  9*i+5,  9*i+5,  9*i+8,  9*i+8 ];
            vecy(counter:counter+5) = [ 3*(i-xsize)+1, 3*(i+xsize)+1, 3*(i-xsize)+2, 3*(i+xsize)+2, 3*(i-xsize)+3, 3*(i+xsize)+3 ];
            vecz(counter:counter+5) = [ -0.25/dy, 0.25/dy, -0.25/dy, 0.25/dy, -0.25/dy, 0.25/dy ];
            counter=counter+6;           
       end
   end
end

% compute 0.5*d(u)/d(z)

if (zsize > 2)  %%compute z gradient only for 3D data
    for i=0:n-1 
        modulo = mod(i, volumesize); 
        if ( modulo < slicesize )
            vecx(counter:counter+5) = [ 9*i+3, 9*i+3,  9*i+6,  9*i+6,  9*i+9,  9*i+9 ]; %border
            vecy(counter:counter+5) = [ 3*i+1, 3*(i+slicesize)+1, 3*i+2, 3*(i+slicesize)+2, 3*i+3, 3*(i+slicesize)+3 ];
            vecz(counter:counter+5) = [ -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz];
            counter=counter+6;        
        else if (modulo > volumesize - slicesize - 1 )
            vecx(counter:counter+5) = [ 9*i+3, 9*i+3,  9*i+6,  9*i+6,  9*i+9,  9*i+9 ]; %border
            vecy(counter:counter+5) = [ 3*(i-slicesize)+1, 3*i+1, 3*(i-slicesize)+2, 3*i+2, 3*(i-slicesize)+3, 3*i+3 ];
            vecz(counter:counter+5) = [ -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz ];
            counter=counter+6;         
        else 
            vecx(counter:counter+5) = [ 9*i+3, 9*i+3,  9*i+6,  9*i+6,  9*i+9,  9*i+9 ];
            vecy(counter:counter+5) = [ 3*(i-slicesize)+1, 3*(i+slicesize)+1, 3*(i-slicesize)+2, 3*(i+slicesize)+2, 3*(i-slicesize)+3, 3*(i+slicesize)+3 ];
            vecz(counter:counter+5) = [ -0.25/dz, 0.25/dz, -0.25/dz, 0.25/dz, -0.25/dz, 0.25/dz ];
            counter=counter+6;
            end
        end
    end
end

% build sparse matrix
K = sparse( vecx, vecy, vecz, 9*n, 3*n );

end

