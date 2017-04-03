clear all
clc
% 
% A = [1 1 1 1 1 0 0 0 0 0 1; -3 -1 1 2 2 1 1 1 1 0 -1; 0 -1 -3 -2 -3 0 0 0 0 0 -2; 0 0 -1 -1 0 0 0 0 0 1 0]
% I = eye(11)
% 
% for i = 1:11
%     for j = 1:4
%         AA(j,i) = A(j,i);
%     end
%     for j = 5:15
%         AA(j,i) = I(j-4,i);
%     end
% end
% AA
% rref(AA)

Name = {'rho','mu','k','L','W','H','D','T','P'}
Unit = {'kg/m3','Pas','W/mk','m','m','m','m','K','Pa'}
[d,f] = unit2si(Unit)