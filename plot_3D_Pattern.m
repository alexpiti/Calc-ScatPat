function plot_3D_Pattern( theta_deg, phi_deg, U_Lin, PlotType )

% FUNCTION plot_3D_Pattern( theta_deg, phi_deg, U_Lin, PlotType )
%
% Plots 3D-pattern U(theta,phi) where:
%  * theta (deg): inclination angle, typically 0...180 or 0...90
%  * phi (deg): azimuthg angle, typically 0...360
%  * U: normalized intensity (W/sterrad), linear scale
% The 2D-arrays (theta,phi,U) are of the same size.
%
% PlotType is an integer index that controls how the data is plotted:
%   1 : |E| = sqrt(U), E-field amplitude, Linear, 3D-spherical 
%   2 : P = |E|^2 = (U), power or intensity, Linear, 3D-spherical
%   3 : P_dB = U_dB, power or intensity, LOG-scale, 3D-spherical [Default]
%   4 : |E| = sqrt(U), E-field amplitude, Linear, 3D-cylindrical
%   5 : P = |E|^2 = (U), power or intensity, Linear, 3D-cylindrical
%   6 : P_dB = U_dB, power or intensity, LOG-scale, 3D-cylindrical
%
% PlotType can be a 2D-array, in which case subplots are generated, e.g.,
% if PlotType = [ 2 3 ; 5 6 ] we get a 2-by-2 panel, with corresp. plots.
%
% You can edit the "plot3Dparams" struct below, for various ploting
% alternatives.
%
% MIT License | Copyright (c) 2022 Alexandros Pitilakis, Thessaloniki/Greece

% -------------------------------------------------------------------------
%% Plot Params (Default values)
% -------------------------------------------------------------------------
global plot3Dparams

% General stuff
plot3Dparams.dB_Lim = 30;
plot3Dparams.ra = 1.2; % length of normalized x/y/z-axis vector
plot3Dparams.cbo = 1; % colorbar on ?
plot3Dparams.View = [50 30]; % view [el,az] angles
plot3Dparams.ScaleColors = 0; % scale colors for uniformity acorss |E|,P,dB?

% Spherical-coord-system params
plot3Dparams.Sph.xyz = 1; % plot xyz axes ?
plot3Dparams.Sph.Sphere = 0.2; % plot Sphere ? 0=No, >0: Boldness (0.1=Barely Vis)
plot3Dparams.Sph.Circles = 1; % plot cut-circles ?

% Cylindrical-coord-system params
plot3Dparams.Cyl.Flat_or_Surf = 1; % 1=flat (pcolor-like), 2=surf
plot3Dparams.Cyl.Circles_degStep = 15; % deg-step for concentric-circles
plot3Dparams.Cyl.Rays_degStep = 30; % deg=step for radial-rays
plot3Dparams.Cyl.Circles.Annotate = 1; % annotate theta=90 ?
plot3Dparams.Cyl.Rays.Annotate = 1; % annotate phi=0 ?
plot3Dparams.Cyl.xyz = 1; % plot xy axes ?
plot3Dparams.Cyl.CapSideWalls = 1 ; % plot "cap" and "SWs" of cylinder ?
plot3Dparams.Cyl.gridCol = 0.2*[1 1 1]; % cyl-grid line-color

% -------------------------------------------------------------------------
%% Test inputs
% -------------------------------------------------------------------------
if nargin == 0
    close all; clc;
        
    phi = 0:1:360; % [deg]
    the = 0:1:90; % [deg]
    [phi_deg,theta_deg] = meshgrid( phi,the ); % [deg]
    
    % Test plot for a 2D array of isotropic scatterers (so that we need
    % only check the Array Factor, AF). The array is set-up to point
    % towards theta=30deg and phi=-120.
    t_max = 30;
    p_max = -120;
    
    % Uniform 2D-array parameters
    Nx = 5; % [.] Number of elements of grid in x-dimenson
    Ny = Nx; % [.] Number of elements of grid in y-dimenson
    kdx = 2*pi * 0.5; % [rad] kappa*dx, where dx: distance in x-dim
    kdy = kdx;        % [rad] kappa*dy, where dy: distance in y-dim
    bx = -kdx*sind(t_max)*cosd(p_max); % [rad] Phase difference in x-dim, from electrical feed
    by = -kdy*sind(t_max)*sind(p_max); % [rad] Phase difference in y-dim, from electrical feed
    
    % Relative phase-difference (from electrical feed and geometric spacing)
    psix = kdx * sind(theta_deg) .* cosd(phi_deg) + bx; % [rad] x-dimension
    psiy = kdy * sind(theta_deg) .* sind(phi_deg) + by; % [rad] y-dimension
    
    % Array Factor in 1D (x- and y-dimensions)
    AF1Dx = sin( 0.5*Nx*psix ) ./ sin( 0.5*psix ); % [.]
    AF1Dy = sin( 0.5*Ny*psiy ) ./ sin( 0.5*psiy ); % [.]
    
    % Correct cases where division 0/0 --> "NaN" ("Not A Number") error
    ix = isnan(AF1Dx);
    iy = isnan(AF1Dy);
    AF1Dx( ix ) = Nx * cos( 0.5*Nx*psix(ix) ) ./ cos( 0.5*psix(ix) ) ;
    AF1Dy( iy ) = Ny * cos( 0.5*Ny*psiy(iy) ) ./ cos( 0.5*psiy(iy) ) ;
    
    % Array Factor in 2D (power)
    AF2D_power  = abs( abs(AF1Dx) .* abs(AF1Dy) ).^2;
    
    % Normalize
    AF2D_power_n = AF2D_power / max( AF2D_power(:) );
    
    %AF2D_power_n = sind(theta_deg).^2;
    
    % Calculate Directivity
    Un = AF2D_power_n;
    D = 4*pi*Un / trapz( phi , trapz( the,Un.*sind(theta_deg) ) ) / (pi/180)^2;
    D0max = max(D(:));
    fprintf(' AF max-directivity is %4.2f [dBi]\n' , 10*log10(D0max) );
    
    U_Lin = D;
    
    
    % Dipole wl/2
    %U_Lin = sind(theta_deg).^2;
    
end
if nargin == 0
    PlotType = [ 1 2 3 ; 4 5 6 ];
    %PlotType = [3 6];
end
if nargin == 3
    PlotType = 3;
    plot3Dparams.cbo = 1; % colorbar on
    plot3Dparams.ScaleColors = 1; % linear scaling of color-axis
    plot3Dparams.dB_Lim = 30; % min-max is 30dB
end

% -------------------------------------------------------------------------
%% Error Checks
% -------------------------------------------------------------------------

% Check theta and phi
if any(size(theta_deg)==1) || any(size(phi_deg)==1)
    error(' ## (theta,phi) must be meshgrid-produced 2D arrays' )
end
if max(theta_deg(:))<90 || max(phi_deg(:))<90
    %error(' ## (theta,phi) must be in DEGREES')
end

% Check U
if any( U_Lin < 0 ) 
    error(' ## Make sure U is in linear power-scale' )
end
if any( ~isreal( U_Lin ) )
    error(' ## Make sure U is real' );
end
if any( size(U_Lin) ~= size(theta_deg) )
    error(' ## Make sure U is same size as theta_deg and phi_deg' );
end
if max(U_Lin(:))~=1
    disp(' ** Normalizing U_Lin for max==1 (or 0dB)' );
    U_Lin = U_Lin/max(U_Lin(:));
end

% -------------------------------------------------------------------------
%% Plots
% -------------------------------------------------------------------------
for kp = 1:numel(PlotType)
    
    % Get proper subplot index
    [I,J] = ind2sub(size(PlotType),kp);
    kpSP = sub2ind(size(PlotType'),J,I);
    
    % Prepare subplot
    hsp = subplot(size(PlotType,1),size(PlotType,2),kpSP);
    
    %hsp = axes;
    %set( hsp, 'Position' , [0.05+(J-1)*0.30, 0.45-(I-1)*0.45, 0.35, 0.45] );
    hold on;
    
    % Prepare data, assuming U_Lin is intensity (units of power), Linear
    switch mod( -1+PlotType(kp),3 )+1 % PlotType(kp)
        case 1,
            r = sqrt(U_Lin);
            titstr = 'Ampl. |E| (Lin)';
        case 2,
            r = U_Lin;
            titstr = 'Power |E|^2 (Lin.)';
        case 3,
            r = 10*log10( U_Lin );            
            titstr = 'Log (dB)';
    end
    
    % Call appropriate plot function (see below), for either 3D-spherical
    % or 3D-cylindrical style
    if PlotType(kp) <=3
        plot_3D_Spherical(theta_deg,phi_deg,r,PlotType(kp)) 
        title( [titstr, ' | Spher.'] );
    else
        plot_3D_Cylindrical(theta_deg,phi_deg,r,PlotType(kp))
        title( [titstr, ' | Cyl.'] );
    end
    
    
end

end

% #########################################################################
%% Plot 3D cylindrical
% #########################################################################
function plot_3D_Cylindrical( t,p,r,PT )

global plot3Dparams

gridCol = plot3Dparams.Cyl.gridCol;

% Preps
ro = r; % store originally supplied r (need for color-scaling)
if (PT-3) == 3  % Check for dB and rescale in [0,1], so that it's plottable:
    R = plot3Dparams.dB_Lim;
    r = (max(r,-R) + R)/R;    
end

% Convert from cylindrical (t,p,r) to cartesian (x,y,z) coordinates:
rhoc = t; % [theta]
phic = p; % [phi]
x = rhoc .* cosd( phic );
y = rhoc .* sind( phic );
z = r;

% If color-scaling was requested
if plot3Dparams.ScaleColors==1 && plot3Dparams.Cyl.Flat_or_Surf==2 % Surf
    if (PT-3) == 3 % dB
        col_scale = 10.^(ro/10);
    elseif (PT-3) == 1 % |E|
        col_scale = ro.^2;
    else
        col_scale = r;
    end    
else % Unscaled or Flat
    col_scale = r;
end
FlatOut = (plot3Dparams.Cyl.Flat_or_Surf==1); % If Flat, this is ==1

% >>>> Plot using SURFACE function <<<<
surf(x,y,z*(1-FlatOut),col_scale);
shading flat
caxis([0 1])

% === Various other plot "cosmetics" ==
xlabel( '\theta cos( \phi )' );
ylabel( '\theta sin( \phi )' );
hold on;
axis tight off;
set(gca,'DataAspectRatio',[1 1 1/90])

% Set plot-viewing angle
switch plot3Dparams.Cyl.Flat_or_Surf    
    case 1, 
        view(2);
    case 2,    
        view( plot3Dparams.View )
end

% If colorbar is on
if plot3Dparams.cbo == 1
    hcb = colorbar;
    if (PT-3) == 3 % dB-plot
        yTick = get(hcb,'YTick');
        yTick_dBn = -(1-yTick)*plot3Dparams.dB_Lim;
        set(hcb,'YTickLabel',num2str(yTick_dBn(:)))
    end
end

% -- plot concentric circles (equal-theta) lines
if plot3Dparams.Cyl.Circles_degStep > 0
    rho = max(x(:));
    a = 0:1:360;
    rho_step = plot3Dparams.Cyl.Circles_degStep; % [deg]
    for myrho = rho:-rho_step:rho_step
        plot3( myrho*cosd(a), myrho*sind(a), 0*a , 'Color' , gridCol );   
    end
    
    % -- annotate concentric circles
    if plot3Dparams.Cyl.Circles.Annotate == 1
       text( 0.75*rho , -0.75*rho , 0 ,...
           sprintf('\\theta=%d^o',round(rho)) ,'FontSize' , 8 );
    end   
    
end

% -- plot radial (equal-phi) lines
if plot3Dparams.Cyl.Rays_degStep > 0
    rho = max(x(:));
    phi_step = plot3Dparams.Cyl.Rays_degStep; % [deg]
    for myphi = 0:phi_step:180-phi_step
        plot3( rho*cosd(myphi+[0 180]), rho*sind(myphi+[0 180]), [0 0] , ...
            'Color' , gridCol );   
    end
    
    % -- annotate radial rays
    if plot3Dparams.Cyl.Rays.Annotate == 1
       text( 1.1*rho , 0 , 0 , '\phi=0' , ...
           'FontSize' , 8 , 'Color' , 'r');
    end
    
end

% -- plot cylinder "cap" and "sidewalls"
if plot3Dparams.Cyl.CapSideWalls == 1 && plot3Dparams.Cyl.Flat_or_Surf == 2
    
    %cap rays
    rho = max(x(:));
    phi_step = plot3Dparams.Cyl.Rays_degStep; % [deg]
    for myphi = 0:phi_step:180-phi_step
        plot3( rho*cosd(myphi+[0 180]), rho*sind(myphi+[0 180]), [1 1] ,...
            'Color' , gridCol );   
    end
    
    %cap circles
    a = 0:1:360;
    rho_step = plot3Dparams.Cyl.Circles_degStep; % [deg]
    for myrho = rho:-rho_step:rho_step
        plot3( myrho*cosd(a), myrho*sind(a), 1+0*a ,...
            'Color' , gridCol );
    end
    
    %vertical beams
    for myphi = 0:phi_step:360-phi_step
        plot3( rho*cosd(myphi)*[1 1], rho*sind(myphi)*[1 1], [0 1],...
            'Color' , gridCol );
    end
    
end

% -- plot xyz "axes"
if plot3Dparams.Cyl.xyz == 1
    ra = (plot3Dparams.ra)*max(t(:));
    a = 0:1:360;
    lw = 2;
    plot3( [0 ra],[0, 0],[0 0], 'r' , 'LineWidth' , lw )
    plot3( [0 0],[0, ra],[0 0], 'g' , 'LineWidth' , lw )
    plot3( max(t(:))*cosd(a), max(t(:))*sind(a), 0*a, 'b' , 'LineWidth' , lw );
    
    if plot3Dparams.Cyl.CapSideWalls == 1 && plot3Dparams.Cyl.Flat_or_Surf == 2
        ra = max(t(:));
        plot3( ra*cosd(a), ra*sind(a), 1+0*a, 'b' , 'LineWidth' , lw );
        plot3( [ra ra -ra -ra 0],[0 0 0 0 0],[0 1 1 0 0], 'r' , 'LineWidth' , lw )
        plot3( [0 0 0 0 0],[ra ra -ra -ra 0],[0 1 1 0 0], 'g' , 'LineWidth' , lw )
        
    end
    
end

end

% #########################################################################
%% Plot 3D spherical
% #########################################################################
function plot_3D_Spherical( t,p,r,PT )

global plot3Dparams

% Check for dB and scale it in [0,1], so that it's plottable:
if PT==3
    ro = r; % store originally-supplied r (needed for color-scaling)
    R = plot3Dparams.dB_Lim;
    r = (max(r,-R) + R)/R;
end

% Convert from spherical (t,p,r) to cartesian (x,y,z) coordinates
x = r .* sind(t) .* cosd(p);
y = r .* sind(t) .* sind(p);
z = r .* cosd(t);

% If color-scaling was requiest
if plot3Dparams.ScaleColors==1
    if PT==3 % dB
        col_scale = 10.^(ro/10);
    elseif PT==1 % |E| (linear)
        col_scale = r.^2;
    else % default P=U=|E|^2 (linear)
        col_scale = r;
    end
else
    col_scale = r;
end

% >>>> Plot using SURFACE function <<<<
surf(x,y,z,col_scale);
shading flat 
caxis([0 1])

% === Now various other plot "cosmetics" ==
hold on;
view( plot3Dparams.View )

% If colorbar is on
if plot3Dparams.cbo == 1
    hcb = colorbar;
    if PT == 3 % dB-plot
        yTick = get(hcb,'YTick');
        yTick_dBn = -(1-yTick)*plot3Dparams.dB_Lim;
        set(hcb,'YTickLabel',num2str(yTick_dBn(:)))
    end
end

% --- plot xyz axes
if plot3Dparams.Sph.xyz == 1;
    ra = plot3Dparams.ra;
    lw = 2;
    plot3( [0 ra],[0 0],[0 0], 'r' , 'LineWidth' , lw ); % x-axis
    plot3( [0 0],[0 ra],[0 0], 'g' , 'LineWidth' , lw ); % y-axis
    plot3( [0 0],[0 0],[0 ra], 'b' , 'LineWidth' , lw ); % z-axis    
end

% -- plot a "sphere" of isotropic (max-U) radiation for reference
if plot3Dparams.Sph.Sphere > 0;
    [xx,yy,zz] = sphere(360/15); % sphere with 15-deg spacing in its (phi,theta)
    if max(t(:)) == 90
        zz = max(zz,0); % crop z<0 part if upper-hemisphere only
    end
    hs = surf(xx,yy,zz);
    set(hs,'FaceColor','None','EdgeColor',0*[1 1 1],...
        'EdgeAlpha',plot3Dparams.Sph.Sphere )    
end

% --- plot "cut-circles" along the three planes xy (blue), xz (green), zy (red)
if plot3Dparams.Sph.Circles == 1;
    a = 0:1:360;
    b = 0:1:(180+180*(max(t(:))==180));
    plot3( 1*cosd(a) , 1*sind(a) , 0*a , 'b' );
    plot3( 0*b , 1*cosd(b) , 1*sind(b) , 'r' );
    plot3( 1*cosd(b) , 0*b , 1*sind(b) , 'g' );
end

axis equal tight off;



end


