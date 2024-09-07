% calculate epsilon at given frequency (wavelength) for given fitting
% x - fitting coefficients
% lb, ub - lower and upper values for fitting coefficients
% w - wavelength (if arg='lambda') or frequency (otherwise)
% terms_num, terms_sz - see definition in fitting.m
function eps=calculate_eps(x,w,terms_num,terms_sz,arg,lb,ub)

if (strcmp(arg,'lambda'))
    w=1/w; % move to frequency representation
end

% moving x coefficients inside (lb,ub) interval
xx = move_inside_period(x,lb,ub);

j=1;
eps=xx(j);
j=j+terms_sz(1);

for i=1:terms_num(2) % Debye
    deps=xx(j);
    gamma=xx(j+1);
    eps=eps+deps/(1-1.i*gamma*w);
    j=j+terms_sz(2);
end
for i=1:terms_num(3) % Drude
    deps=xx(j);
    w0=xx(j+1);
    gamma=xx(j+2);
    eps=eps+deps*w0^2/(-2.i*gamma*w-w^2);
    j=j+terms_sz(3);
end
for i=1:terms_num(4) % Lorentz
    deps=xx(j);
    w0=xx(j+1);
    gamma=xx(j+2);
    eps=eps+deps*w0^2/(w0^2-2.i*gamma*w-w^2);
    j=j+terms_sz(4);
end
for i=1:terms_num(5) % modified Lorentz
    deps=xx(j);
    w0=xx(j+1);
    gamma=xx(j+2);
    gamma_num=xx(j+3);    
    eps=eps+deps*(w0^2-1.i*gamma_num*w)/(w0^2-2.i*gamma*w-w^2);
%    eps=eps-(deps*w0*w0-1.i*gamma_num*w)/(w*w+2.i*gamma*w-w0*w0);        
%    eps = eps + deps*w0*(exp(1.i*gamma_num)/(w0-w-1.i*gamma) + exp(-1.i*gamma_num)/(w0+w+1.i*gamma));
    j=j+terms_sz(5);
end
