classdef DensityPeaksIR < handle
    %CLUSTERREGGROW Summary of this class goes here
    %   This is the class of the proposed method.
    
    properties
        numSeeds = 20;
        tarRate = 0.01*0.15;
        thdQuatile = 3;
    end
    
    methods
         
        function [tarPos, tarCon] = finalDetect(obj, rhoMat)
            m = size(rhoMat, 1);
            [rho, delta] = iterationElection( obj, rhoMat );
            [ classInitial ] = singularFind( obj, rho, delta );
            singularIndex = find( classInitial ~=  0 );
            classCenterRows = mod( singularIndex, m );
            classCenterRows(classCenterRows == 0) = m;
            classCenterCols = ceil( singularIndex / m );
            seedPos = [ classCenterCols, classCenterRows ];
            grayJump = regionGrow( obj, rhoMat, seedPos );
            confidence = confidenceCal( obj, grayJump );
            posIndex = confidence > obj.thdQuatile;
            tarPos = seedPos(posIndex, :);
            tarCon = confidence(posIndex);
        end
        
        [ rho, delta ] = iterationElection( obj, rhoMat );
        [ classInitial ] = singularFind( obj, rho, delta );
        [ grayJump ] = regionGrow( obj, rhoMat, seedPos );
        [ conf ] = confidenceCal( obj, inputVector );
        
    end
    
    methods( Access = protected )
        [ outputMat, outputSeq ] = matShrink( obj, inputMat );
        [ areaSize ] = calRefSize( obj, rows, cols );
    end
end



