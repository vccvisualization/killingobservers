%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2018 King Abdullah University Of Science and Technology 
%
% Contact: 
% Peter Rautek peter.rautek@kaust.edu.sa
% Matej Mlejnek matej.mlejnek@kaust.edu.sa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ D ] = getD( v, deldvdx, deldvdy, deldvdz, n, xsize, zsize, tsize, slicesize, volumesize, dx, dy, dz, dt )

%add del v(u) component
vecx = [1:3*n,1:3*n,1:3*n];
vecy = zeros(9*n,1);
vecy(1:3:end) = [1:3:3*n-2,2:3:3*n-1,3:3:3*n];
vecy(2:3:end) = vecy(1:3:end);
vecy(3:3:end) = vecy(1:3:end);
vecz = [deldvdx(:),deldvdy(:),deldvdz(:)];
vecz = vecz(:);
counterdelvu = 9*n+1;

counter = counterdelvu;

% add -del u(v) component
for i=0:n-1 
     modulo = mod(i, xsize); 
     if (modulo == 0)   
           vecy(counter:counter+5) = [ 3*i+1, 3*i+4, 3*i+2, 3*i+5, 3*i+3, 3*i+6]; %border
           vecz(counter:counter+5) = -[ -1/dx, 1/dx, -1/dx, 1/dx, -1/dx, 1/dx ];
     else if (modulo == xsize-1)
           vecy(counter:counter+5)= [ 3*i-2, 3*i+1, 3*i-1, 3*i+2, 3*i, 3*i+3]; %border
           vecz(counter:counter+5) = -[ -1/dx, 1/dx, -1/dx, 1/dx, -1/dx, 1/dx ];
         else
           vecy(counter:counter+5) = [ 3*i-2, 3*i+4, 3*i-1, 3*i+5, 3*i, 3*i+6 ];
           vecz(counter:counter+5) = -[ -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx, -0.5/dx, 0.5/dx ];
         end
     end
     vecx(counter:counter+5) = [ 3*i+1, 3*i+1, 3*i+2, 3*i+2, 3*i+3, 3*i+3];
     counter=counter+6;
end

vecz(counterdelvu:end) = vecz(counterdelvu:end) .* transpose([repelem(v(1:3:end),6)]);

for i=0:n-1 
   modulo = mod(i, slicesize); 
   if (modulo < xsize ) %first row border
            vecy(counter:counter+5) = [ 3*i+1, 3*(i+xsize)+1, 3*i+2, 3*(i+xsize)+2, 3*i+3, 3*(i+xsize)+3 ]; 
            vecz(counter:counter+5) = -[ -1/dy, 1/dy, -1/dy, 1/dy, -1/dy, 1/dy ];
   else if ( modulo > slicesize - xsize - 1 ) %last row border
            vecy(counter:counter+5) = [ 3*(i-xsize)+1, 3*i+1, 3*(i-xsize)+2, 3*i+2, 3*(i-xsize)+3, 3*i+3  ]; 
            vecz(counter:counter+5) = -[ -1/dy, 1/dy, -1/dy, 1/dy, -1/dy, 1/dy ];
       else
            vecy(counter:counter+5) = [ 3*(i-xsize)+1, 3*(i+xsize)+1, 3*(i-xsize)+2, 3*(i+xsize)+2, 3*(i-xsize)+3, 3*(i+xsize)+3 ];
            vecz(counter:counter+5) = -[ -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy, -0.5/dy, 0.5/dy ];
       end
   end
   vecx(counter:counter+5) = [ 3*i+1, 3*i+1, 3*i+2, 3*i+2, 3*i+3, 3*i+3 ];
   counter=counter+6;
end

vecz(counterdelvu+6*n:end) = vecz(counterdelvu+6*n:end) .* transpose([repelem(v(2:3:end),6)]);

if (zsize > 2)  %%compute z gradient only for 3D data
    for i=0:n-1 
        modulo = mod(i, volumesize); 
        if ( modulo < slicesize ) %first slice border
            vecy(counter:counter+5) = [ 3*i+1, 3*(i+slicesize)+1, 3*i+2, 3*(i+slicesize)+2, 3*i+3, 3*(i+slicesize)+3 ]; 
            vecz(counter:counter+5) = -[ -1/dz, 1/dz, -1/dz, 1/dz, -1/dz, 1/dz ];
        else if (modulo > volumesize - slicesize - 1 ) %last slice border
            vecy(counter:counter+5) = [ 3*(i-slicesize)+1, 3*i+1, 3*(i-slicesize)+2, 3*i+2, 3*(i-slicesize)+3, 3*i+3]; 
            vecz(counter:counter+5) = -[ -1/dz, 1/dz, -1/dz, 1/dz -1/dz, 1/dz ];
            else 
            vecy(counter:counter+5) = [ 3*(i-slicesize)+1, 3*(i+slicesize)+1, 3*(i-slicesize)+2, 3*(i+slicesize)+2, 3*(i-slicesize)+3, 3*(i+slicesize)+3 ];
            vecz(counter:counter+5) = -[ -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz, -0.5/dz, 0.5/dz  ];
            end
        end
        vecx(counter:counter+5) = [ 3*i+1, 3*i+1, 3*i+2, 3*i+2, 3*i+3, 3*i+3 ];
        counter=counter+6;           
    end
    vecz(counterdelvu+12*n:end) = vecz(counterdelvu+12*n:end) .* transpose([repelem(v(3:3:end),6)]);
end

% add -du/dt component
if (tsize > 2)  %%compute t gradient only for time-dependent data
    for i=1:3*n 
        if ( i < 3*volumesize + 1  ) %first volume border
             vecx(counter:counter+1) = [ i, i ];
             vecy(counter:counter+1) = [ i, i+3*volumesize ]; 
             vecz(counter:counter+1) = -[ -1/dt, 1/dt ];
             counter=counter+2;
        else if (i > 3*volumesize*tsize - 3*volumesize )  %last volume border
             vecx(counter:counter+1) = [ i, i ];
             vecy(counter:counter+1) = [ i-3*volumesize, i ];
             vecz(counter:counter+1) = -[ -1/dt, 1/dt ];
             counter=counter+2;
            else 
                vecx(counter:counter+1) =  [ i, i ];
                vecy(counter:counter+1) =  [ i-3*volumesize, i+3*volumesize ];
                vecz(counter:counter+1) = -[ -0.5/dt, 0.5/dt ];
                counter=counter+2;
            end
        end
    end
end

D = sparse(vecx, vecy, vecz, 3*n, 3*n);

end

