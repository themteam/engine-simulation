function [res] = valODE45_Pression( theta,P,V)
global gamma

res=zeros(length(P),1);
for i=1:length(P)
    res(i)=P(i)*V(i)^(gamma);
end 
ref= P(1)*V(1)^(gamma);
plot(theta,res,'b*-',theta,ref,'ro-');
end

