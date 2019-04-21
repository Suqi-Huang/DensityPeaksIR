function [ classInitial ] = singularFind( obj, rho, delta)
%SINGULARFIND Find the density peaks in rho-delta space.
% rho--density of each pixel, a vector with the size of [m*n, 1].
% delta--delta-space of each pixel, a vector with the size of [m*n, 1].
% classInitial--the sign vector with the size of [m*n, 1], the positions 
% with non-zero values represent density peaks.

nLength = length(rho);
classInitial = zeros(nLength, 1);

%% find the singular points with product
product = rho.*delta;
[productSort, index] = sort(product, 'descend');
cl = 1;
% threashold = max(product) / obj.rhoDeltaThd;  % set the threashold
threashold = productSort(obj.numSeeds);
for i = 1 : nLength
    if ( product(index(i)) >= threashold )
        classInitial( index(i) ) = cl;
        cl = cl+1;
    end
end

end

