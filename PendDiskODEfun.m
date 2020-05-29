function dZdT = PendDiskODEfun(t,Z,Torques,Forces,return_lambda) 

n = length(Z);
if nargin < 5
    return_lambda = false;
elseif return_lambda
    n = n+4;
end
dZdT = zeros(n,1);

q = Z(1:6);
dq = Z(7:12);

dZdT(1:6) = dq; % Velocities
[A,Q] = AQ_PendDisk_num([q(:);dq(:)]',Torques,Forces);

ddqL = A\Q;
dZdT(7:12) = ddqL(1:6);
if return_lambda
    dZdT(13:n) = ddqL(7:10);
end