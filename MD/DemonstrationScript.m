clear all
close all
clc

INPUT = 2;

switch INPUT
    case 1, load('../Data/FEM/Single Emitter/Data_1.mat')
    case 2, load('../Data/FEM/Single Emitter/Data_1.mat')
    case 3, load('../Data/FEM/Single Emitter/Data_1.mat')
end

% Constructing time and position vectors
TimeVectorSec   = (0:(size(Data,2)-1)) / SamplingFrequencyHz;

aux = Positions(1,:);
[ ~, Radius, ~ ] = cart2pol ( aux(:,1), aux(:,2), aux(:,3) );
PerimeterMeters = (0:(size(Data,1)-1))/size(Data,1) * pi * 2 * Radius;

% Surfing Raw Data
figure( 1 );
colormap jet;
surf( TimeVectorSec, PerimeterMeters, Data );
shading flat;
view( 2 );
axis tight;
ylabel( 'Arc Position [ m ]' );
xlabel( 'Time [ s ]' );

% Mode Order Decomposition 
[ DecomposedData, Order ] = ModeDecompositionPhaseShift ( Data );

Envelope_dec_data = abs( hilbert( DecomposedData' ) )';

% Surfing Decomposed Data
figure( 2 );
colormap jet;
surf( TimeVectorSec, Order, DecomposedData );
shading flat;
view( 2 );
axis tight;
ylabel( 'Arc Position [ m ]' );
xlabel( 'Time [ s ]' );

% Plot Decomposed Data
figure( 3 );
subplot 211
plot( TimeVectorSec, 20*log10( Envelope_dec_data(1:4,:) ) );
xlabel( 'Time [ s ]' );
ylim([-90 0])
legend('Order = 0','Order = 1','Order = 2','Order = 3');
subplot 212
plot( TimeVectorSec, Envelope_dec_data(1:4,:) );
xlabel( 'Time [ s ]' );
legend('Order = 0','Order = 1','Order = 2','Order = 3');


openfig('..\Data\Dispersion\DispersionFig.fig')

