%% Optimal Stretch
% <Description>
%
% Syntax:
%                         Residue = OS( Damage, Base2Str );
%   [Residue, E, beta, DBeta_min] = OS( Damage, Base2Str );
%
% Inputs: 
%   Damage   : <Description>
%   Base2Str : <Description>
%
% Outputs:
%    Residue   : <Description>
%    E         : <Description>
%    beta      : <Description>
%    DBeta_min : <Description>
%
% See also
%
% Author  : Henrique Tormen Haan de Oliveira
% Date    : 2021
% Version : 1.0
%
% Reference BibTex
% no ref
function [Residue, E, beta, DBeta_min] = OS( Damage, Base2Str )
DEBUG = 0;

MODE = 'spline';
% MODE = 'linear';

if DEBUG
    %%
    clearvars -except Damage Base2Str DEBUG MODE
end

%% Input Shift
[Damage2  , Nshift] = shiftdim(Damage);
Base2Str2           = shiftdim(Base2Str);

%%
% %% TIME CUT
pini     = 1;
pfin     = numel( Damage2 );

% %% BETA DEF
beta_step = 0.005;
beta = [0:beta_step:.05 (.05+2*beta_step):2*beta_step:.2 (.2+4*beta_step):4*beta_step:.8]/100;
beta = cat(2,-1*flip(beta(2:end),2),beta);

% %% INI
fa    = 2;
pad   = 2^nextpow2(size(Damage2,1));
dw    = fa/pad;%[metade]
w_pre = ( 0 : 1 : (pad/2 - 1) )';
w     = w_pre * dw;

% FFT + CUT IN HALF
dfft_full  = fft( Damage2  , pad, 1); 
bfft_full  = fft( Base2Str2, pad, 1); 

dfft = dfft_full( 1:end/2, 1 );
bfft = bfft_full( 1:end/2, 1 );

% WINDOW
dfft = dfft .* tukeywin(size(dfft,1),0.025);
bfft = bfft .* tukeywin(size(bfft,1),0.025);

for i = 1 : numel(beta)
    %%
    dw_str = dw * (1 + beta(i));
    w_str  = w_pre * dw_str;

    bfft_str = interp1( w_str, bfft, w, MODE );

    rfft = dfft - bfft_str;
    if any(isnan(rfft))
        rfft = rfft(~isnan(rfft));
        size(rfft)
    end
    
    E( i, 1 ) = sqrt( mean( rfft .* conj(rfft) ) );
    
    if DEBUG && 1
        %%
        clf
        subplot 211;hold all
            plot(abs(bfft    ),'.-','Color',[1 1 1]*.2)
            plot(abs(bfft_str),'.-')
            plot(abs(dfft    ),'.-')
            plot(abs(rfft    ),'.-')
                xlim([0 400])
        ha = subplot(2,2,3);hold all
            plot(beta,'Visible','off')
            plot(beta(1:i),E,'r.')
            plot(beta(i),E(i))
        
        drawnow
        pause(.005)
    end
end
%%

[~,i_mins] = min(E,[],1);
DBeta_min  = beta(i_mins);


if DEBUG
    %%
    clf
    hold all
    plot(beta        ,E        ,'.-','DisplayName','Error curve')
    plot(beta(i_mins),E(i_mins),'ro','DisplayName','Min Error')
    ylabel('Error')
    xlabel('beta parameter')
    legend show
    grid on
    ylim([0 max(E)])
end


%%
dfft = dfft_full( 1:end/2, 1 );
bfft = bfft_full( 1:end/2, 1 );

% WINDOW
dfft = dfft .* tukeywin(size(dfft,1),0.025);
bfft = bfft .* tukeywin(size(bfft,1),0.025);

dw_str = dw * (1 + DBeta_min);
w_str  = w_pre * dw_str;

bfft_str = interp1( w_str, bfft, w, MODE );
% bfft_str = interp1( w_str, bfft, w, 'linear' );

rfft = dfft - bfft_str;
% if any(isnan(rfft))
%     error('1')
% end
rfft(isnan(rfft)) = 0;

Err = sqrt( mean( rfft .* conj(rfft) ) );

rfft2 = cat( 1, rfft    , 0, flip( rfft    (2:end,1) ,1 ) );

r = ifft( rfft2   , [], 1, 'symmetric' );

Residue = r ( 1 : numel(Damage2) , 1);

if DEBUG
    %%
    b_str = cat( 1, bfft_str, 0, flip( bfft_str(2:end,1) ,1 ) );
    dfft2 = cat( 1, dfft, 0, flip( dfft(2:end,1) ,1 ) );
    
    b_ = ifft( bfft_str, [], 1, 'symmetric' );
    d2 = ifft( dfft2   , [], 1, 'symmetric' );
    
    clf;
    subplot 211,hold all;
        plot(Damage2)
        plot(d2)
        plot(r)
    subplot 212;hold all;
        plot(abs(dfft))
        plot(abs(bfft))
        plot(abs(bfft_str))
        plot(abs(rfft))
    
end

%% Output Shift
Residue = shiftdim(Residue, - Nshift);
end