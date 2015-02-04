function W = kernel(r, h, der)
%KERNEL  Calculate the kernel value, its gradient, or Hessian.
%    The specified location 'r' must be a column vector of size D
%    (one of 1, 2, or 3). The support of the kernel is a D-dimensional
%    sphere of radius 2*h.
%    W = kernel(r, h, 0)  returns the (scalar) value of the kernel
%    W = kernel(r, h, 1)  returns grad(W), a row vector of size D
%    W = kernel(r, h, 2)  returns Hess(W), a matrix of size DxD.

% Problem dimension
dim = length(r);

% Set the scaling factor as appropriate for the
% current dimension (so that W is a partition of
% unity in the continuous case).
switch dim
    case 1
        factor = 6 * h;
    case 2
        factor = 2.8 * pi * h^2;
    case 3
        factor = 4 * pi * h^3;
end

q = norm(r)/h;

% Special case: q = 0
if q < 1e-6
   switch der
       case 0
           W = 4;
       case 1
           W = zeros(1, dim);
       case 2
           W = -(12/h^2) * eye(dim);
   end
   
   W = W / factor;
   
   return
end

% Get the value and derivatives of the cubic function.
[Wq, Wq_p, Wq_pp] = cubic_function(q);

% Assemble the output quantity, depending on the requested
% derivative level.
switch der
    case 0
        W = Wq;
    case 1
        W = Wq_p / (q * h^2) * r';        
    case 2
        W = Wq_p / (q * h^2) * eye(dim) + ...
            (Wq_pp - Wq_p / q) / (q^2 * h^4) * (r * r');
end


% Scale the output quantity.
W = W / factor;

return

% -------------------------------------

function [Wq, Wq_p, Wq_pp] = cubic_function(q)

if (q >= 2)
    Wq = 0;
    Wq_p = 0;
    Wq_pp = 0;
elseif (q >= 1)
    Wq = (2-q)^3;
    Wq_p = -3 * (2-q)^2;
    Wq_pp = 6 * (2-q);
else
    Wq = (2-q)^3 - 4*(1-q)^3;
    Wq_p = -3*(2-q)^2 + 12*(1-q)^2;
    Wq_pp = 6*(2-q) -24*(1-q);
end

return


  

