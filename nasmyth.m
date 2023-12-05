function [phi,varargout] = nasmyth(varargin)
% NASMYTH Nasmyth universal shear spectrum.
%
% There are 4 basic forms of the function:
% 1) [phi,k] = nasmyth(e,nu,N)
% 2) phi     = nasmyth(e,nu,k)
% 3) [phi,k] = nasmyth(0,N);
% 4) phi     = nasmyth(0,k);
%
% Form 1)
% Returns the Nasmyth spectrum for the dissipation rate e (W/kg)
% and viscosity nu (m^2/s).  The length of the returned spectrum 
% phi is N, k is the wavenumber in cpm. 
% Default values are nu = 1e-6 and N = 1000, if nu or N are omitted.
% Example calls: [phi,k] = nasmyth(1e-7,1.2e-6,512)
%                [phi,k] = nasmyth(1e-7,1.2e-6)
%                [phi,k] = nasmyth(1e-7)
%
% Form 2)
% Same as form 1) except that the Nasmyth spectrum is evaluated at 
% the wavenumbers given in the input vector k (in cpm).
% Example calls: phi = nasmyth(1e-7,1.2e-6,logspace(-1,3,512));
%
% Form 3)
% Returns the non-dimensional Nasmyth spectrum (G2 spectrum) of length
% N points. The wavenumber is k=k'/ks [where k' is in cpm (see Oakey 1982)]
% runs from 1e-4 to 1.
%
% Form 4)
% Same as 3 except that the non-dimensional spectrum is evaluated at
% the wavenumbers given in the input vector k (in k'/ks).
% Example calls: phi = nasmyth(0,logspace(-3,0,512));
%
% Note: For forms 1) and 2), the dissipation rate can be a vector, 
% e.g. e = [1e-7 1e-6 1e-5],  in which case phi is a matrix whose 
% columns contain the scaled Nasmyth spectra for the elements in e. 
%
% The form of the spectrum is computed from Lueck's (1995) fit to the 
% Nasmyth points listed by Oakey (1982).
% 
% References: Oakey, N. S., 1982: J. Phys. Ocean., 12, 256-271.
%             Lueck, R. G., 1995: personal communication.
%             Wolk et al., 2000: J. Atmos. Ocean. Tech., submitted.

% Development history:
% 1st Version By D. Huang, SEOS, UVic.
% 1996/02/05, Rolf Lueck, CEOR UVic: G2's formula fitted
% 1997/07/01, Fabian Wolk, SEOS UVic: allow input for vector epsilon
% 2000/07/28, Fabian Wolk, Alec Electronics: non-dimensional output (forms 3 and 4)
% 2000/09/28, Fabian Wolk, Alec Electronics: wavenumber input (form 2)

% Included in TurbTools with permission of Rolf Lueck.
% Fabian Wolk 2001/02/03
% 2001 Alec Electronics Co., Ltd.
% #Revision 0.0: 2001/08/09#

% argument checking
error(nargchk(1,3,nargin));
[scaled,e,nu,N,k] = checkArgs(varargin,nargin);

if scaled  % forms 1) and 2)
   e = e(:)';
   Ne = length(e);
   ks = (e./nu.^3).^(1/4); % Kolmogorov wavenumber(s)    
   ks = ks(ones(N,1),:);   
   if isempty(k) % form 1)
      x = logspace(-4,0,N)';
      x = x(:, ones(1,Ne)); 
   else          % form 2)
      k = k(:,ones(1,Ne));
      x = k./ks;
   end
   G2 = 8.05*x.^(1/3)./(1+(20*x).^(3.7)); % Lueck's fit
   k = x.*ks;    
   e = e(ones(N,1),:);   
   phi = e.^(3/4) * nu^(-1/4) .* G2;
   varargout{1} = k;            
else    % forms 3) and 4)
   if isempty(k) % form 3)
      k = logspace(-4,0,N)';  % k = k_hat/k_s, as in Oakey 1982.
      varargout{1} = k;
   end
   phi = 8.05*k.^(1/3)./(1+(20*k).^(3.7)); % Lueck's fit
end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [scaled,e,nu,N,k] = checkArgs(arg,narg);
% Helper function for Nasmyth.m

% the default values:
scaled = 0;
e  = 1e-6;
nu = 1e-6;
N  = 1000;
k  = [];


if any(arg{1} ~= 0) % it's form 1) or 2)
   scaled = 1;
   if narg == 3
      e = arg{1};
      nu = arg{2};
      N = arg{3};
   elseif narg == 2
      e = arg{1};
      nu = arg{2};
   elseif narg == 1
      e = arg{1};
   end
else                % it is form 3) or 4)
   scaled = 0;
   N = arg{2};
end

if length(N) > 1 % last argument is vector means it's a wavenumber vector
   if all(size(N)>1)
      error('Sorry, can''t have matrix wavenumber input.');
   else
      k = N(:);
      N = length(k);
   end
end



   
   
