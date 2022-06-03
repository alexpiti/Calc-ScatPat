function [E,the,phi] = calc_ScatPat( freq , MetaSurf , Illumin , Method )

% FUNCTION [E,the,phi] = calc_ScatPat( freq , MetaSurf , Illumin , Method )
%
% Calculates the scattering pattern produced on the upper hemisphere by 
% the illumination of given abstracted flat metasurface with an EM field 
% at given frequency. Currently supporting only plane and spherical wave
% illumination, and a rectangular square-cell metasurface of Mx-by-Ny
% cells. The configuration of the metasurface is given in terms of complex
% (amplitude+phase) reflection coefficients of each cell, assumed uncoupled
% from the rest of the cells, obviously corresponding to the given 
% frequency, polarization and incidence direction.
%
% This function implements the Huygens-Frensel diffraction principle.
% https://en.wikipedia.org/wiki/Huygens%E2%80%93Fresnel_principle
%
% ===== Inputs =====
% Check the "nargin==0" block for values of these arguments
%  * frequency [Hz] scalar // CW frequency for this pattenern
%  * MetaSurf -- structure defining the MS properties, with fields:
%     .duc [m] scalar // cell-width (square)
%     .NumCellsXY [.] 1x2 vector // no of cells in Ny/rows and Mx/columns
%     .Refl_Ampli [.]  Ny-by-Mx matrix // reflection-coeff amplitude profile 
%     .Refl_Phase [rad] Ny-by-Mx matrix // reflection-coeff phase profile
%     .Cell_ScatPat_Exp [.] scalar // exponent in expression: cos(theta)^n
%  * Illumin -- structure defining the illumination/incident wave properties
%     .Type [.] Illumination/source type: 1=Plane wave, else=Spherical
%     .DOA [deg] 1x2 vector // Direction-of-Arrival for plane wave [theta,phi]
%     .xyz [m] 1x3 vector // coordinates to the spherical wave center
%  * Method -- structure with misc parameters for the calculation
%     .SPAR [deg] scalar // Scattering-Pattern Angle-Resolution
%     .focusDir [deg] 1x2 vector // "focus" calc on a specific direction
%     .angspan [deg] scalar // if set to zero, we scan entire hemisphere,
%          else theta,phi are scanned at +/-(angspan/2) around focusDir
%
% ===== Outputs =====
%  * E -- scattering pattern (E-field complex amplitude)
%  * [the,phi] -- matrices from meshgrid, with the directions corresponding
%       to the scattering patter. They are: theta=0:SPAR:90, phi=0:SPAR:360
%
% Alexandros Pitilakis | Thessaloniki Greece | June 2022
%
% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis

% Function plots also the MS configuration. Set this to zero to disable.
Plot_MS_config = 1; 

% Constants
c0 = 3e8; % [m/s] Light speed in vacuum

% ------------------------------------------
% Define test params (used when F5'ing here)
% ------------------------------------------
if nargin == 0
    
    clc; close all;

    freq = 1e12; % [Hz] Operating frequncy

    MetaSurf.duc = c0/freq/8; % [m] cell-width (square)
    MetaSurf.NumCellsXY = [ 40 40 ]; % [.] number of cells in each dimension
    MetaSurf.Refl_Ampli = ones(  MetaSurf.NumCellsXY(1) , MetaSurf.NumCellsXY(2) );
    MetaSurf.Refl_Phase = 2*pi*rand( MetaSurf.NumCellsXY(1) , MetaSurf.NumCellsXY(2) );  % rad
    MetaSurf.Cell_ScatPat_Exp = 0; % [.] exponent in expression: cos(theta)^n

    Illumin.Type = 1; % Illumination/source type: 1=Plane, else=Spherical
    Illumin.DOA = [ 45 , 0 ]; % [deg] theta+phi angles (plane waves)
    Illumin.xyz = [ 0 , 0 , c0/freq*5 ]; % [m] sphere center (spherical)

    Method.SPAR = 3; % [deg] % SPAR: Scattering-Pattern Angle-Resolution
    Method.focusDir = [45,180]; % [deg](theta,phi) to "focus" calc on
    Method.angspan = 0; % [deg] if set to zero, we scan entire hemisphere
            
end

% ------------------------------------------
% Preparations
% ------------------------------------------

% Retrieve params from input
duc = MetaSurf.duc; % [m] unit-cell size (square)
N = MetaSurf.NumCellsXY(1); % y-dim (up/down, #-of-row)
M = MetaSurf.NumCellsXY(2); % x-dim (left/right, #-of-column)
Rp_mn = MetaSurf.Refl_Phase; % [rad] phase profile
Ra_mn = MetaSurf.Refl_Ampli; % [.] amplitude profile

ST = Illumin.Type; % Source/Illumination type: 0=Plane, 1=Spherical
thei = Illumin.DOA(1); % [deg] theta incident/illum (for plane)
phii = Illumin.DOA(2); % [deg] phi incident/illum (for plane)
xyz = Illumin.xyz; % [m] sphere center (for spherical)

SPAR = Method.SPAR; % Scat-Pat Ang-Resolution
thef = Method.focusDir(1); % [deg] lobe-focus theta
phif = Method.focusDir(2); % [deg] lobe-focus phi
angspan = Method.angspan; % [deg] ang-span on lobe-focus-direction
cosexp = MetaSurf.Cell_ScatPat_Exp; % [.] exponent in expression: cos(theta)^n


% Some derived params
wl = c0 / freq; % [m] free-space wavelength
k0 = 2*pi/wl; % [rad/m] free-space wavenumber
[mm,nn] = meshgrid(1:M,1:N); % index-array in MS rows/columns
dxk = k0 * duc; % unit cell size, x-dimension*kappa (d*k0=2*pi * d/wl)
dyk = k0 * duc; % unit cell size, y-dimension*kappa
xyzks = k0*xyz; % [.] spher-source xyz, normed to kappa0


% ------------------------------------------
% Plot MS configuration (phase & ampl)
% ------------------------------------------
if Plot_MS_config == 1
    
    % Phase-profile across the MS
    subplot(1,2,1)
    imagesc(fliplr(Rp_mn*180/pi)); set(gca,'YDir','Normal')
    axis equal tight; 
    xlabel('x-dim ("M") #cells'); ylabel('y-dim ("N") #cells');
    title( '\angle\Phi_{nm} (deg)');
    set(gca,'XMinorTick','on','YMinorTick','on')
    colormap(jet);
    caxis([0 360])
    hcb1 = colorbar;
    set(hcb1, 'YTick', 0:45:360 ); 
    
    % Uneven amplitude plot across the MS
    subplot(1,2,2)
    imagesc(fliplr(Ra_mn)); axis equal tight; set(gca,'YDir','Normal')
    hcb2 = colorbar;    
    xlabel('x-dim ("M") #cells'); ylabel('y-dim ("N") #cells');
    title( '|r_{nm}|');
    set(gca,'XMinorTick','on','YMinorTick','on')
    caxis([0 1])

end

% ------------------------------------------
% Source configuration
% ------------------------------------------

% Direction of source (phi,theta), as "seen" from each unit-cell
if ST == 1
    % For Plane-Wave (PW) incidence, theta & phi are constant across the MS
    %   e.g. for normal inc.: theta=0 (and phi=X, where X="indifference")
    theSmn = zeros(N,M) + thei*pi/180; % (rad)
    phiSmn = zeros(N,M) + phii*pi/180; % (rad)
else
    % For Spherical-Wave (SW, point source), theta & phi are rigorously
    % calculated from geometric formulas.
    
    % 0. Calc unit-cell coords
    [xuc,yuc] = meshgrid( (-(M-1)/2:(M-1)/2)*dxk , (-(N-1)/2:(N-1)/2)*dyk ) ;
    zuc = zeros(size(xuc));
    % 1. Transform source coords
    xsrc = xyzks(1); %(xyzks(1)-0.5)*(M-1)/2*dxk;
    ysrc = xyzks(2); %(xyzks(2)-0.5)*(N-1)/2*dyk;
    zsrc = xyzks(3);
    % 2. cartesians dists of source from each unit-cell center
    xd_sfuc = xsrc - xuc;
    yd_sfuc = ysrc - yuc;
    zd_sfuc = zsrc - zuc;
    % 3. transform to sphericals
    dist_suc = sqrt(xd_sfuc.^2 + yd_sfuc.^2 + zd_sfuc.^2);
    theSmn = atan(sqrt(xd_sfuc.^2 + yd_sfuc.^2)./zd_sfuc); % (rad)
    phiSmn = atan2( yd_sfuc, xd_sfuc ); % (rad)
    
end


% Source (input) wave-front amplitude & phase on each unit-cell. Amplitude
% is usually assumed ==1, as inc. wave travels in lossless medium (air),
% comes from farfield and thus illuminates the metasurface homogeneously.
% Phase depends on the wave "shape" (plane, spherical, gaussian etc).
if ST == 1
    % * For normal PW incidence: A=1 and phase=const (e.g. zero)
    % * For oblique PW incidence in x0z plane: phase = k0*dx*sin(theta)
    Ia_mn = ones(N,M);
    Ip_mn = ( dxk*(mm-1).*cos(phiSmn) - dyk*(nn-1).*sin(phiSmn) ).*sin(theSmn) ;
else
    % Point Source:
    Ia_mn = ones(N,M);
    Ip_mn = dist_suc; % Because phase=beta*dz-->k0*distance
end


% =========================================================================
% Huygens-Fresnel Principle for Scattering Pattern ("hemisphere")
% =========================================================================

% Angle-resolution of Scatt-Patt:
dsp = SPAR; % degree-step-phi
dst = SPAR; % degree-step-theta

% Directions on "hemisphere" space. If angspace>0, it only calculates
% around a specific direction (e.g. around the approx reflection-lobe) 
% to cut-down on simulation times.
if angspan > 0
    phi1 = pi/180*( (-angspan/2:dsp:angspan/2) + phif );
    the1 = pi/180*( (-angspan/2:dst:angspan/2) + thef );
    
    if any( the1 < 0 )
        in = the1<0;
        phi1( in ) = phi1(in)+180;
        the1 = abs(the1);
    end

else
    phi1 = (0:dsp:360) *pi/180;
    the1 = (0:dst:90)  *pi/180;
end

% Directions in the "hemisphere", gridded
[phi,the]=meshgrid(phi1,the1);

% Scat-Patt formation: 
E = 0*the; % this will hold the pattern (complex)

% We scan all the cells in the MS, and sum-up their contribution.
cc = 0;
fprintf( ' -- Scatt-Patt calc. progress:\n' );
for m1 = 1:M
    for n1 = 1:N
        
        cc = cc+1;
        ProgBarString( cc, M*N )
        
        % "Interfere" the complex-valued illumination/incidence profile
        % with the cell reflection profile.
       auxFix = cos( theSmn(n1,m1) ).^cosexp ...
                .*Ia_mn(n1,m1).*exp(1j*Ip_mn(n1,m1)) ...
                .*Ra_mn(n1,m1).*exp(1j*Rp_mn(n1,m1));        
        
        % Add the contributions of source-interaction w/ each unit-cell
        % *** The two exp's here are the "zeta_mn" in Hamidreza PhD
        auxi = auxFix...
            .*exp( +1j*dxk*m1.*sin(the).*cos(phi))...
            .*exp( -1j*dyk*n1.*sin(the).*sin(phi));
        
        % Add to pattern
        E = E + cos( the ).^cosexp .* auxi;
        
    end
end
fprintf( 'Done.\n' );

% Test plot
if nargin == 0
    
    try
        figure;
        Unorm = abs(E).^2/max(abs(E(:)).^2);
        plot_3D_Pattern( the*180/pi, phi*180/pi, Unorm , [3 6] )
    catch
        error( ' Get "plot_3D_Pattern" from GitHub to plot the pattern' )
    end    
    
end


end

% #########################################################################
% Misc Routines
% #########################################################################

function ProgBarString( ni , N , L )

% Prints a progress-bar string in the CMD window. Meant for use in for
% loops, e.g. with ni=3 and N=20 inputs will print a 15% filled bar that
% spans the entire command window width. Optional argument L controls the
% lenght of the bar, e.g., if L=20 then the bar width will be 20 monospaced
% characters on the CMD window.
%
% Alexandros Pitilakis, Thessaloniki/Greece, 2020


% Test inputs
if nargin == 0
    ni = 3;
    N = 20;
end

% Get CMD window character width (L) and scale ni by it
if nargin < 3
    L = max( get( 0 , 'CommandWindowSize' ) ) - 4; 
end
n = round( ni*L/N );

% Prepare progress-bar string
str = [ '|' , repmat(' ',[1,L]) , '|' ];
if n==1,
    str(2) = '>';    
elseif n<L    
    str(2:(2+n-1)) = '-';
    str(2+n) = '>';
elseif n==L
    str = [ '|' , repmat('=',[1,L]) , '|' ];
end

% Add text with completion percentage in the middle of the string
numstr = sprintf( ' %3d %% ' , round(ni/N*100) );
iis = round(L/2)-3;
length(numstr);
str( iis:(iis+length(numstr)-1) ) = numstr;

% Vocalize in CMD
if nargout == 0
    if ni > 1 % for ni=2,3,... you gotta backspace (delete) previous str
        fprintf( repmat( '\b' , [1,L+3] ) );
    end
    fprintf( '%s\n' , str );
end

end
















