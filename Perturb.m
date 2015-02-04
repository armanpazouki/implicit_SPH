function out = Perturb(vec3, thresh)
out = vec3 + thresh.*(2 * rand(size(vec3)) - 1);
