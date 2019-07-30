function K = errFnKernelmorlet( t1,t2,M )

Kerr = (erf((M - (t1+t2))/(2)) + erf(( M + (t1+t2))/(2)));

% Just the integral over the kerenl...
% K = (1/2)*(a*d)*(sqrt(pi))*exp(-(((t1-t2)/(2*a*d))^2))*Kerr;

% Including all constant factors except complex exp

K = (1/(2*M)).*exp(-(((t1-t2)/(2)).^2)).*Kerr;

K = K.*exp(-1i*2*pi*(t1-t2));
