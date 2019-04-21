function [ rho, delta ] = iterationElection( obj, rhoMat )
%iterationElection Compute the density and delta-distance of each pixel.
% rhoMat--density map (as input image in this paper), a matrix with a size of [m,n].
% rho--the density of each pixel, a vector with the size of [m*n, 1].
% delta--the delta-distance of each pixel, a vector with the size of [m*n, 1].

%%  Initialization output result
[nRows, nColumns] = size(rhoMat);
rho = reshape(rhoMat, nRows*nColumns, 1);   % density
delta = zeros(nRows*nColumns, 1);   % close distance
neighbor = zeros(nRows*nColumns, 1);    % neighbor

%%  Initialization the intermediate value
% real position in feature space, x and y coodinate.
[numRows, numColumns] = size(rhoMat);
pxVector = (1 : numRows)';
pyVector = 1 : numColumns;
pxMat = repmat( pxVector , 1, numColumns);
pyMat = repmat(pyVector, numRows, 1);
% neighbor is probably zero. so let minimum of rhoMat is 
% realmin('double'). 
rhoMat = rhoMat + realmin('double');

%% circulation and calculate
while( numRows > 1 || numColumns > 1 )
    % Padding zeros boundaries.
    rhoPadMat = zeros(numRows+2, numColumns+2);
    pxPadMat = zeros(numRows+2, numColumns+2);
    pyPadMat = zeros(numRows+2, numColumns+2);
    signPadMat = zeros(numRows+2, numColumns+2);
    % Initialize
    rhoPadMat( 2 : numRows+1 , 2 : numColumns+1 ) = rhoMat;
    pxPadMat( 2 : numRows+1 , 2 : numColumns+1 ) = pxMat;
    pyPadMat(  2 : numRows+1 , 2 : numColumns+1 ) = pyMat;

    % Finding the neighbor of every pixel.
    % Sliding 3*3 window and compare and compute distance.
    for i = 2 : numRows+1
        for j = 2 : numColumns+1
            pxC = pxPadMat(i, j);   % center position
            pyC = pyPadMat(i, j);
            if (pxC ~= 0 && pyC ~= 0)   % find the cluster point. 
                rhoPatch = rhoPadMat( i-1:i+1 , j-1:j+1 );
                rhoVector = rhoPatch(:);
                pxPatch = pxPadMat( i-1:i+1 , j-1:j+1 );
                pyPatch = pyPadMat( i-1:i+1 , j-1:j+1 );
                pxVector = pxPatch(:);
                pyVector = pyPatch(:);
                maxRho = max( rhoVector );
                rhoPoint = rhoVector(5);
                
                if ( rhoPoint < maxRho )    % centre point not local maximum.
                    candidateIndex = ( rhoVector > rhoPoint );
                    pxCandidate = pxVector(candidateIndex);
                    pyCandidate = pyVector(candidateIndex);
                    distCandidate = sqrt( ( pxCandidate - pxC ).^2 + ...
                        ( pyCandidate - pyC ).^2 );
                    rhoCandidate = rhoVector(candidateIndex);
                    attraction = rhoCandidate ./ distCandidate;
                    [~, pos] = max(attraction);               
                    pxN = pxCandidate(pos);
                    pyN = pyCandidate(pos);
                    neighbor( (pyC-1)*nRows + pxC ) = ...
                                (pyN-1)*nRows + pxN;    %   set neighbor.
                    delta( (pyC-1)*nRows + pxC ) = distCandidate(pos);
                else                     % rhoPoint is maximum
                    maxIndex = ( rhoVector == maxRho );
                    nMax = length( maxIndex );
                    if (nMax == 1)
                        % sign for locl maximum.
                        signPadMat(i, j) = 1;     
                    else
                        [~, maxPos] = sort(rhoVector, 'descend');
                        if( maxPos(1) == 5 )
                            signPadMat(i, j) = 1;
                        else
                            pxN = pxVector( maxPos(1) );
                            pyN = pyVector( maxPos(1) );
                            neighbor( (pyC-1)*nRows + pxC ) = ...
                                (pyN-1)*nRows + pxN;    %   set neighbor.
                            delta( (pyC-1)*nRows + pxC ) = ...
                                sqrt( (pxN-pxC)^2+(pyN-pyC)^2 );
                        end
                    end
                end
                
            end     % end if(pxC ~= 0 || pyC ~= 0)

        end
    end

    % shrink
    numRowsU =ceil( numRows / 2);
    numColumnsU = ceil( numColumns / 2 );
    % upper layer parameter
    rhoMatU = zeros(numRowsU, numColumnsU);
    pxMatU = zeros(numRowsU, numColumnsU);
    pyMatU = zeros(numRowsU, numColumnsU);

    signPadMat = signPadMat > 0;
    % current layer signMat
    signMat = signPadMat( 2:numRows+1, 2:numColumns+1 );
    % upper layer signMat
    % matShrink is shrinking function.
    [signMatU, seqOrder] = matShrink( obj, signMat );
    % ’˚¿ÌÀ≥–Ú
    rhoExtract = rhoMat(signMat);
    rhoMatU(signMatU) = rhoExtract(seqOrder);
    pxExtract = pxMat(signMat);
    pxMatU(signMatU) = pxExtract(seqOrder);
    pyExtract = pyMat(signMat);
    pyMatU(signMatU) = pyExtract(seqOrder);

    rhoMat = rhoMatU;
    pxMat = pxMatU;
    pyMat = pyMatU;

    [numRows, numColumns] = size(rhoMat);
end

delta( (pyMat-1)*nRows + pxMat ) = sqrt( nRows^2 + ...
    nColumns^2 );
% the maximum neighbor is itself.
neighbor( (pyMat-1)*nRows + pxMat ) = ...
    (pyMat-1)*nRows + pxMat;    % the maximum neighbor

end

