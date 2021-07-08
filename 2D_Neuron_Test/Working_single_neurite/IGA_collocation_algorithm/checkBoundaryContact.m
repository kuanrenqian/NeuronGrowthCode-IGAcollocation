function [state] = checkBoundaryContact(phi)

state = 0;

len = length(phi);
M = sqrt(len);
phi = reshape(phi,M,M);

threshould = 0.1;

if ( sum(sum(phi(1:5,:))) >threshould || ...
        sum(sum(phi(end-5:end,:))) >threshould || ...
        sum(sum(phi(:,1:5))) >threshould || ...
        sum(sum(phi(:,end-5:end))) >threshould)
    sum(sum(phi(1:5,:)))
    sum(sum(phi(end-5:end,:)))
    sum(sum(phi(:,1:5)))
    sum(sum(phi(:,end-5:end)))
    state = 1;
end
