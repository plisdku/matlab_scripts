function x=cummax(x)
% x = cummax(x)
% by Pierre Grinspan (copied from a Matlab forum)
n=size(x,1);
if n==1
elseif n==2, x(2,:)=max(x);
else
x(2:2:n,:)=cummax(max(x(1:2:n-1,:),x(2:2:n,:)));
x(3:2:n,:)=max(x(3:2:n,:),x(2:2:n-1,:));
end


