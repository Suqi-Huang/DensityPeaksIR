function [ outputMat, outputSeq ] = matShrink( ~, inputMat )
%MATSHRINK Summary of this function goes here
%   One step in the calculation of delta-distance.

[ nRowsIn, nColumnsIn ] = size(inputMat);
nRowsOut = ceil( nRowsIn / 2 );
nColumnsOut = ceil( nColumnsIn /2 );
outputMat = zeros(nRowsOut, nColumnsOut);

[rowsIn, columnsIn] = find( inputMat > 0 );
rowsOut = ceil( rowsIn / 2 );
columnsOut = ceil( columnsIn / 2 );
p = rowsOut + (columnsOut-1)*nRowsOut;
outputMat(p) = 1;
outputMat = outputMat > 0;

rowVector = (1 : nRowsIn)'; 
columnVector = 1 : nColumnsIn;
rowMat = repmat(rowVector, 1, nColumnsIn);
columnMat = repmat(columnVector, nRowsIn, 1);

pxIn = rowMat( inputMat );
pyIn = columnMat( inputMat );
pxOut = ceil( pxIn /2 );
pyOut = ceil( pyIn / 2 );
feature= (pyOut-1) * nRowsOut + pxOut;
[~, outputSeq] = sort( feature, 'ascend' );

end

