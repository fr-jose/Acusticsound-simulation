% ******************************************************
% FFP
% ******************************************************
clear all

% *** set up Munk SSP ***
zs  = 1000.0;			% source depth
f   = 50; omega = 2 * pi * f	% frequency in Hertz
eps = 0.00737; c0  = 1500;

d   = 5000;			% bottom depth
nz  = 500;			% number of finite-difference points
h   = d / nz; h2  = h * h;	% mesh spacing
z   = linspace( 0, d, nz );	% grid coordinates

x = 2 * ( z - 1300 ) / 1300;
c = c0 * ( 1 + eps * ( x - 1 + exp( -x ) ) );

plot( z, c ); view( 90, 90 );
xlabel( 'Depth (m)' ); ylabel( 'Sound Speed (m/s)' )

% ******************************************************
% set up finite difference matrix
% ******************************************************

v = omega * omega ./ ( c .* c );

%D = sparse( 1:nz, 1:nz  , -2*ones( 1, nz   ) / h2 + v, nz, nz, nz   );
%E = sparse( 2:nz, 1:nz-1,    ones( 1, nz-1 ) / h2,     nz, nz, nz-1 );
%A = D + E + E';
D = -2*ones( nz, 1 ) / h2 + v';
E =    ones( nz, 1 ) / h2;
A = spdiags([E D E], -1:1, nz, nz);

% ******************************************************
% Solve for successive wavenumbers
% ******************************************************

nk = 512;				% # of k points, for FFT = 2^n

kmin = omega / 1550.0;
kmax = omega / 1500.0;
deltak = ( kmax - kmin ) / ( nk - 1 );
alpha = deltak;			% offset into complex plane

k = linspace( kmin, kmax, nk ) - i * alpha * ones( 1, nk );

% computer range vector
nr     = nk;
rmax   = 2.0 * pi / deltak;
deltar = rmax / nr;
r = linspace( deltar, rmax, nr );

isd = zs / d * nz;		% index of source depth
e = zeros( nz, 1 );
e( isd ) = 1.0 / h;		% point source forcing

% solve the system repeatedly

for ik = 1:nk
   if ( mod( ik, 50 ) == 0 )
       ik % display every 50th step
   end
   D = spdiags( k( ik )^2 * ones( nz, 1 ), 0, nz, nz );
   B = A - D;
   %g( :, ik ) = B \ e;
   if ( ik == 1 )
       g( :, ik ) = B \e;
   else
       %M=ones( size( B ) );
   g( :, ik ) = symmlq( B, e );  %0.001, [], squeeze( g(:, ik-1) ) );
end
end

% ******************************************************
% FFT and display results
% ******************************************************

% --- subsample G

takez = 1:5:nz;
zt = z( takez );
gt = g( takez, : );

figure
pcolor( real( k ), zt, abs( gt ) ); ...
shading flat; colormap( gray ); colorbar; view( 0, -90 );
xlabel( 'Wavenumber (1/m)' ); ylabel( 'Depth (m)' );
title('Green''s function')

figure
plot( real( k ), abs( gt( 8, : ) ) ); ...
xlabel( 'Wavenumber (1/m)' ); ylabel( 'Magnitude' );
title('Green''s function')

% --- FFT with sqrt( k ) weight to get pressure
gt = full( gt * spdiags( sqrt( k ), 0, nk ) );
p = fft( gt' );

% --- range factor
p =  ( deltak / sqrt( 2.0 * pi ) ) * p' * ...
     spdiags( exp( alpha * r ) ./ sqrt( r ), 0, nr );
tl = full( 20.0 * log10( abs( p ) ) );

% subsample and plot

taker = 1:nr;
rt = r( taker );
tlt = tl( :, taker );

figure;
pcolor( rt, zt, tlt ); ...
axis( [ 0 100000 0 5000 ] ); caxis( [ -100 -60 ] ); ...
shading flat; colormap( gray ); colorbar; view( 0, -90 );
xlabel( 'Range (m)' ); ylabel( 'Depth (m)' );
title('FFP intensity')

figure
plot( rt, tlt( 8, : ) ); ...
xlabel( 'Range (m)' ); ylabel( 'Intensity' );
title('FFP intensity')
axis( [ 0 100000 -100 -60 ] );
