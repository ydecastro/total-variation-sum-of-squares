% bivariate example
mpol xp 2
mup = meas(xp);
mpol xm 2
mum = meas(xm);
% build problem data
d = 12; % degree
p = genpow(3,d); p = p(:,2:end); % powers
b = zeros(size(p,1),1);
Xp = [-1/2 1/2; 1/2 -1/2; 1/2 1/2; 0 0]; Xm = [0 -1/2; 1/2 0];
for i = 1:size(Xp,1), b = b + Xp(i,1).^p(:,1) .* Xp(i,2).^p(:,2); end % positive atoms
for i = 1:size(Xm,1), b = b - Xm(i,1).^p(:,1) .* Xm(i,2).^p(:,2); end % negative atoms
% build semidefinite relaxation
k = 12; % order of the relaxation
ME = [mass(mup)-mass(mum)==b(1)];
for i = 2:size(p,1)
 ME = [ME; mom(xp(1)^p(i,1).*xp(2)^p(i,2))-mom(xm(1)^p(i,1).*xm(2)^p(i,2))==b(i)];
end
MC = [xp(1)^2<=1, xp(2)^2<=1, xm(1)^2<=1, xm(2)^2<=1];
P = msdp(min(mass(mup)+mass(mum)),ME,MC,k);
% semidefinite program in SeDuMi format
[A,b,c,K,b0,s] = msedumi(P);
% solve semidefinite program
[x,y] = sedumi(A,b,c,K);
% extract moment matrices
z=c-A'*y;
Mp=reshape(z(K.f+K.l+(1:K.s(1)^2)),K.s(1),K.s(1));
Mm=reshape(z(K.f+K.l+K.s(1)^2+(1:K.s(2)^2)),K.s(2),K.s(2));
% lower bound on the objective function
disp(['lower bound=' num2str(b0+s*c'*x)]);
% try to extract the atoms
try
 mext(Mp,2,2*k) % positive atoms
 mext(Mm,2,2*k) % negative atoms
end
% display polynomial certificate
close all
[X1,X2]=meshgrid(linspace(-1,1,1e3));
certif=x(1:K.f)'*mmon(xp,d);
X3=eval(vectorize(certif,'X1','X2'));
mesh(X1,X2,X3,'FaceColor','interp','EdgeColor','none','FaceLighting','phong');
view(-15,15)
axis vis3d
camlight left
xlabel x_1
ylabel x_2
colormap spring
hold on
for i = 1:size(Xp,1), plot3(Xp(i,1),Xp(i,2),1,'.r','markersize',30); end
for i = 1:size(Xm,1), plot3(Xm(i,1),Xm(i,2),-1,'.b','markersize',30); end

