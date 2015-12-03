% %Verification de la dérivée
% for i=1:length(teta)-1
%     res(i)=V(i+1)-dvdta(i)*(teta(i+1)-teta(i));
% end
% subplot(3,1,1);
% plot(teta,V);
% grid on
% subplot(3,1,2)
% plot(teta(1:end-1),res);
% grid on
% subplot(3,1,3);
% plot(teta,dvdta);
% grid on

% figure;plot(theta,y(:,3))
% 
% valODE45_Pression(theta,y(:,3),V,y0(3))
% 
% dP=y(:,3);
% P=zeros(length(dP),1);
% P(1)=P0;
% 
% for i=2:length(theta)
%     P(i)=P(i-1)+dP(i-1)*(theta(i)-theta(i-1));
% end
% plot(P)