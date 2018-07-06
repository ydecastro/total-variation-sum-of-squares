k = input('relaxation order=');
mpol xp
mup = meas(xp);
mpol xm
mum = meas(xm);
% build problem data
p = (0:k)'; % powers
B = zeros(size(p));
Xp = [-3/4,1/2]; Xm = [1/8];
for i = 1:length(Xp), B = B + Xp(i).^p; end  % positive atoms
for i = 1:length(Xm), B = B - Xm(i).^p; end % negative atoms
% build semidefinite relaxation
P = msdp(min(mass(mup)+mass(mum)),...
       mom(xp.^p)-mom(xm.^p)==B,...
       (xp+1)*(xp-1)<=0, (xm+1)*(xm-1)<=0);
% semidefinite program in SeDuMi format
[A,b,c,K,b0,s] = msedumi(P);
% solve semidefinite program
[x,y] = sedumi(A,b,c,K);
% extract moment matrices
z=c-A'*y;
Mm=reshape(z(K.f+K.l+(1:K.s(1)^2)),K.s(1),K.s(1));
Mp=reshape(z(K.f+K.l+K.s(1)^2+(1:K.s(2)^2)),K.s(2),K.s(2));
% lower bound on the objective function
disp(['lower bound=' num2str(b0+s*c'*x)]);
% try to extract the atoms
try
 mext(Mm,1,k) % positive atoms
 mext(Mp,1,k) % negative atoms
end
% display polynomial certificate
zp=rot90(x(1:K.f),-1);
X=linspace(-1,1,1e3);
close all
plot(X,polyval(zp,X),'k','linewidth',3)
hold on
plot([-1 1],[1 1],':k');plot([-1 1],[-1 -1],':k');
for i = 1:length(Xp), plot(Xp(i),1,'.r','markersize',30); end
for i = 1:length(Xm), plot(Xm(i),-1,'.b','markersize',30); end
axis([-1 1 -1.1 1.3])
