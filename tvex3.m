rand('state',31400) % random seed fixed so that the example is the one in the paper
Xp = []; Xm = [];
for i = 1:5 % random points on the sphere
 xs = randn(1,3); xs = xs/norm(xs); Xp = [Xp;xs]; % positive atoms
 xs = randn(1,3); xs = xs/norm(xs); Xm = [Xm;xs]; % negative atoms
end
%%
minimal_separation          = pi;
for i=1:4
    for j=(i+1):5
        minimal_separation  = min([minimal_separation ...
                                acos(Xp(i)*Xp(j)) acos(Xm(i)*Xm(j)) ...
                                acos(Xm(i)*Xp(j)) acos(Xp(i)*Xm(j))]);
    end
end
degree_separation_condition  = ceil(2/minimal_separation);
fprintf('Separation condition requires at least degree %i\n', degree_separation_condition);

%%
mpol xp 3
mup = meas(xp);
mpol xm 3
mum = meas(xm);
k = input('relaxation order=');
% build problem data
p = genpow(4,k); p = p(:,2:end); % powers
b = zeros(size(p,1),1);
for i = 1:size(Xp,1), b = b + Xp(i,1).^p(:,1) .* Xp(i,2).^p(:,2) .* Xp(i,3).^p(:,3); end % positive atoms
for i = 1:size(Xm,1), b = b - Xm(i,1).^p(:,1) .* Xm(i,2).^p(:,2) .* Xm(i,3).^p(:,3); end % negative atoms
ME = [mass(mup)-mass(mum)==b(1)];
for i = 2:size(p,1)
 ME = [ME; mom(xp(1)^p(i,1).*xp(2)^p(i,2).*xp(3)^p(i,3))-mom(xm(1)^p(i,1).*xm(2)^p(i,2).*xm(3)^p(i,3))==b(i)];
end
MC = [xp(1)^2+xp(2)^2+xp(3)^2==1; xm(1)^2+xm(2)^2+xm(3)^2==1];
P = msdp(min(mass(mup)+mass(mum)),ME,MC);
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
 mext(Mp,3,k) % positive atoms
 mext(Mm,3,k) % negative atoms
end
%% display polynomial certificate on the sphere
close all
[X1,X2,X3]=sphere(500);
certif=x(1:size(p,1))'*mmon(xp,k);
X4=eval(vectorize(certif,'X1','X2','X3'));
h=surf(X1,X2,X3,X4,'EdgeColor','none');%,'FaceLighting','phong');
axis vis3d
camlight left
xlabel x_1
ylabel x_2
zlabel x_3
view(45,45)
hold on
for i = 1:size(Xp,1), plot3(Xp(i,1),Xp(i,2),Xp(i,3),'.r','markersize',20); end
for i = 1:size(Xm,1), plot3(Xm(i,1),Xm(i,2),Xm(i,3),'.b','markersize',20); end
colorbar

