function W = kernel(r, h, der)
%KERNEL  Calculate the kernel value, its gradient, or Hessian.
%    The specified location 'r' must be a 2D column vector. 
%    The support of the kernel is a circle of radius 2*h.
%    W = kernel(r, h, 0)  returns the (scalar) value of the kernel.
%    W = kernel(r, h, 1)  returns grad(W), a 2D row vector.
%    W = kernel(r, h, 2)  returns Hess(W), a 2x2 matrix.

% Set the scaling factor as appropriate for a 2D problem (so that W is
% a partition of unity in the continuous case).
factor = 2.8 * pi * h^2;

q = norm(r)/h;

% Special case: q = 0
if q < 1e-6
   switch der
       case 0
           W = 4;
       case 1
           W = [0 0];
       case 2
           W = -(12/h^2) * eye(2);
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
        W = Wq_p / (q * h^2) * eye(2) + ...
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


  

