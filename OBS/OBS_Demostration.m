%% LOADING DATA STRUCTURES
clear all
close all
clc

count = 0;
for i = 1:16
    count = count + 1;
    base_raw(count,1) =  load(sprintf('../Data/FEM/Elasticity Modulus Study/Baseline_%d.mat',i));
end
count = 0;
for i = 1
    count = count + 1;
    data_raw(count,1) = load(sprintf('../Data/FEM/Elasticity Modulus Study/Data_Damage_%d.mat',i));
end
MatPropBase   = load('../Data/FEM/Elasticity Modulus Study/MatInputs_base.mat');
MatPropDamage = load('../Data/FEM/Elasticity Modulus Study/MatInputs_damage.mat');
clear count i
%% EXTRACT DATA FROM DATA STRUCTURES
% Extract the matrices form the structures and concatenate them into a
% single array. the size is in the form of M x N x O, where M is the members,
% N is the channels or decomposition order and M is the time samples.
clear base
for i = 1 : size(base_raw,1)
    base(i,:,:) = base_raw(i,1).Data;
end
for i = 1 : size(data_raw,1)
    data(i,:,:) = data_raw(i,1).Data;
end
%% PLOT INPUT CHANNEL AVERAGES
figure(1);clf;
subplot 211;hold all
for i = 1 : size(base,1)
    plot(squeeze(mean(base(i,:,:),2)) , 'DisplayName', sprintf('Modulus: %.1f GPa',MatPropBase.Eo(i)/1e9))
end
xlabel('Time samples')
ylabel('Normalized Amplitude')
grid on
title('Baseline Data')
legend show

subplot 212;hold all
for i = 1 : size(data,1)
    plot(squeeze(mean(data(i,:,:),2)) , 'DisplayName', sprintf('Modulus: %.1f GPa',MatPropDamage.Eo(i)/1e9))
end
xlabel('Time samples')
ylabel('Normalized Amplitude')
grid on
title('Damage Data')
%% MODE EXTRACTION
for i = 1 : size(base,1)
    aux = shiftdim(base(i,:,:),1);
    aux = ModeDecompositionPhaseShift ( aux );
    base_decomp ( i, :, : ) = aux;
end
for i = 1 : size(data,1)
    aux = shiftdim(data(i,:,:),1);
    aux = ModeDecompositionPhaseShift ( aux );
    data_decomp ( i, :, : ) = aux;
end
if 0
    %% Excldude data ( The model (FEM) containing the damage was simulated with the same material properties as the baseline model #5, so the difference will be as close as possible from the perfect substraction)
    base_decomp([2 3 4 5],:,:)=[];
end
%% PLOT DECOMPOSITION
% Select which circunferential order to use in the following steps: 1-T01; 2-F12; 3-F22 ...
iMode = 1;

figure(2);clf;
subplot 211;hold all;
for i = 1 : size(base_decomp,1)
    plot(squeeze(base_decomp(i,iMode,:)))
end
subplot 212;hold all
for i = 1 : size(data_decomp,1)
    plot(squeeze(data_decomp(i,iMode,:)))
end
%% APPLIES OBS FUNCTION
[ ~, BestBase, ~, ErrSort ] = OBS ( data_decomp(:,iMode,:) , base_decomp(:,iMode,:) );

BestBase = ErrSort(5);

Residual = data_decomp - base_decomp(BestBase,:,:);

%% PLOT AFTER OBS
iData = 1;
iMode = 1;

figure(3);clf
    subplot 211; hold all
        plot( squeeze( data_decomp(    iData       , iMode, : ) ), '.-', 'DisplayName', 'Data'       );
        plot( squeeze( base_decomp( BestBase(iData), iMode, : ) ), '.-', 'DisplayName', 'BestBase'   );
        plot( squeeze( Residual   (    iData       , iMode, : ) ), '.-', 'DisplayName', 'Difference' );
            grid on
            xlabel('Time samples')
            ylabel('Normalized Amplitude')
            title('OBS Result (linear)')
    subplot 212; hold all
        plot( 20*log10( abs( hilbert( squeeze( data_decomp(    iData       , iMode, : ) ) ) ) ), '-', 'DisplayName', 'Data' );
        plot( 20*log10( abs( hilbert( squeeze( base_decomp( BestBase(iData), iMode, : ) ) ) ) ), '-', 'DisplayName', 'BestBase' );
        plot( 20*log10( abs( hilbert( squeeze( Residual   (    iData       , iMode, : ) ) ) ) ), '-', 'DisplayName', 'Difference' );
            grid on
            xlabel('Time samples')
            ylabel('Amplitude 0-refernce')
            title('OBS Result (log)')
%% Applies OS
[Rst, E, beta, beta_min]  = OS( data_decomp(    iData       , iMode, : ), base_decomp( BestBase(iData), iMode, : ));

% %% PLOT AFTER OS
figure(4);clf
    subplot 211;hold all;
        plot( 20*log10( abs( hilbert( squeeze( data_decomp(    iData       , iMode, : ) ) ) ) ), '.-', 'DisplayName', 'T(0,1) Damage' );
        plot( 20*log10( abs( hilbert( squeeze( base_decomp( BestBase(iData), iMode, : ) ) ) ) ), '.-', 'DisplayName', 'T(0,1) Base'   );
        plot( 20*log10( abs( hilbert( squeeze( Residual   (    iData       , iMode, : ) ) ) ) ), '.-', 'DisplayName', 'T(0,1) Difference after OBS' );
        plot( 20*log10( abs( hilbert( squeeze( Rst        (      1         ,   1  , : ) ) ) ) ), '.-', 'DisplayName', 'T(0,1) Diff after OS' );
            grid on
            legend show
            xlabel('Time samples')
            ylabel('Amplitude 0-refernce')
    subplot 223;hold all;
        plot( beta    , E                , '.-'  );
        stem( beta_min, E(beta==beta_min), 'ro-' );
            grid on
            xlim([min(beta) max(beta)])
            xlabel('beta parameter');
            ylabel('RMS Error');