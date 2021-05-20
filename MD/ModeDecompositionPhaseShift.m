%% ModeDecompositionPhaseShift
% Mode decomposition "simple signal processing" coded using as base the
% Mr. Lowe, MJS 1998 paper.
%
% Syntax: 
%   Decomposition = ModeDecompositionPhaseShift ( TimeTraces );
%
% Inputs:
%   TimeTraces  : Scalar acoustic data, 2-dimensional matrix ( N x M ). N 
%                 is the number of channels around a pipe and M is the time
%                 sample count.
%
% Outputs:
%   OrderVector   : Stores the circunferencial order vector ( 1 x M ) where
%                   M is the time sample count.
%   Decomposition : Scalar acoustic data, 2-dimensional matrix ( O x M )
%                   where O is the Order and M is the time samples count.
%
% See also
%
% Author  : Henrique Tormen Haan de Oliveira
% Date    : 2021
% Version : 1.0
%
% Reference BibTex
% @article{lowe1998mode,
%   title={The mode conversion of a guided wave by a part-circumferential notch in a pipe},
%   author={Lowe, MJS and Alleyne, DN and Cawley, P},
%   year={1998}
% }
function [ Decomposition, OrderVector ] = ModeDecompositionPhaseShift ( TimeTraces )
N   = size( TimeTraces, 1 );
pad = size( TimeTraces, 2 );

U_t = fft( TimeTraces, pad, 2 );

Decomposition( 1, : ) = mean( TimeTraces, 1 );
for i = 2 : N / 2
    ang = ( ( 0 : N-1 ) / N )' * 2*pi *( i - 1 ); % 1x32
    
    p = repmat( exp( -1i*ang ),[ 1 size( U_t, 2 ) ] );
    U_tt = U_t .* p;
    
    u_tt = ifft( U_tt, [], 2, 'symmetric' );
    
    Decomposition( i, : ) = mean( u_tt, 1 );
    OrderVector = 0 : ( size( Decomposition, 1 ) - 1 );
end
end

% Extracted from reference
% pg 4 - Section 4 Rsults of Reflections Study
% "Some simple signal processing was necessary in order to
% determine the amplitude of each of the refiected modes. An
% identical methodology was applied to both the experimental
% and the finite element results. For the refiection of the axially
% symmetric L(0,2) mode, the 16 individual signals from the
% transducers ( or 32 signals from the nodes) were simply added.
% The resulting signal was thus exactly as if the transducers were
% wired together, as reported in previous studies ofL(0,2) refiection.
% For the other two modes, a phase delay of NB!21f was
% added to each signal before summing them. N is the circumferential
% order number and (} is the angular distance from the
% centre of the notch. Thus a separate processing calculation was
% performed in order to extract the amplitude of each of the three
% modes from the multiple transducer ( or node) records. Since
% the signals were rather narrow band, the processing could reasonably
% have been performed directly on the raw time records.
% However, for better accuracy, the calculations were performed
% in the frequency domain, and then inverted to give a processed
% time-domain record for each mode."
