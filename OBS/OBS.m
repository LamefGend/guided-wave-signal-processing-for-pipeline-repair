%% Optimal Baseline Subtraction
% From 2 lists, Data and Baseline, returns the baseline member that has the
% lower RMS error. (not tested to very long baseline lists)
%
% Syntax:
%                                       Residual = OBS ( Data , Baseline );
%   [ Residual, BestBase, RMS_ERR, ErrorSorted ] = OBS ( Data , Baseline );
%
% Inputs: 
%   Data     : List of GW data, in the form M x N x O, where M is the
%              members, N is the channels or decomposition order and M 
%              is the time samples.
%              This list is the data of interest.
%   Baseline : List of GW data, in the form M x N x O, where M is the
%              members, N is the channels or decomposition order and M 
%              is the time samples.
%              This list is the baseline.
%
% Outputs:
%   Residual    : <Description>
%   BestIndexes : <Description>
%   RMS_ERR     : <Description>
%   ErrorSorted : <Description>
%
% See also
%
% Author  : Henrique Tormen Haan de Oliveira
% Date    : 2021
% Version : 1.0
%
% Reference BibTex
% no ref
function [ Residual, BestIndexes, RMS_ERR, ErrorSorted ] = OBS ( Data , Baseline )
%% Redimensioning the data 
D = permute( Data    , [1 4 2 3] ); %keep the damage   list in the first dimension
B = permute( Baseline, [4 1 2 3] ); %keep the baseline list in the first dimension

%% Matching sizes
D = repmat( D, [     1      size(B,2) 1 1 ]);
B = repmat( B, [ size(D,1)      1     1 1 ]);

%% Calculating the residuals to all combinations
Res = D - B;

%% Calculating the differnce RMS 
RMS_ERR = sqrt(sum(sum(Res.^2,4)/size(D,4),3)/size(D,3));

%% Sorting and picking
[ ~, ErrorSorted ] = sort(RMS_ERR,2);
BestIndexes = ErrorSorted ( : , 1);

%% Redo the calculation using the best fit
Residual = Data - Baseline(BestIndexes,:,:);

end