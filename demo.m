
clear;
%% ---------------------------------------------------
% Centralized PDFPO^2 to get the exact minimizer x*

d = 50;
p = 500;

%% Generate data
sample = 20;
mu1 = zeros([1,p]);
sigma1 = zeros([p,p]);
for i=1:p
    for j=1:p
        sigma1(i,j) = 0.5^(abs(i-j));
    end
end
A = zeros(sample,p,d);
y = zeros(sample,d);
hatx = zeros([p,1]);
for i=1:20
    hatx(i) = 1;
end
for i=1:d
    A(:,:,i) = mvnrnd(mu1,sigma1,sample);
    y(:,i) = A(:,:,i)*hatx + randn(sample,1);
end

%% set parameters and initialization

L = 0;
for i=1:d
     Xi = A(:,:,i);
     L = max(L,norm(Xi'*Xi));
end
gamma = 1.9/L;
[B,vd] = get_group(p,3);
lam = 1/norm(full(B*B'));

mu = 1;
x0 = zeros([p,1]);
v0 = zeros([vd,1]);
nablafx0 = zeros([p,1]);
for i=1:d
    Ai = A(:,:,i);
    yi = y(:,i);
    nablafx0 = nablafx0 + Ai'*(Ai*x0-yi)/d;
end
gt = nablafx0;
ibbt = sparse(1:vd,1:vd,1,vd,vd)- lam*(B*B');

%% Centralized PDFPO^2
T = 4000;
for k=1:T
    energy = mu*norm(B*x0,1);
    for i=1:d
        Ai = A(:,:,i);
        yi = y(:,i);
        energy = energy + 0.5*norm(Ai*x0-yi)^2/d;
    end
    if(mod(k,fix(T/10))==1)
        fprintf('k=%d,energy=%f\n', k,energy);
    end

    vt = B*(x0-gamma*gt)+ibbt*v0;
    
    thresh = mu*gamma/lam;
    vt(vt>thresh) = thresh;
    vt(vt<-thresh) = -thresh;
    
    xt = x0 - gamma*gt - lam*(B'*vt);
    
    x0 = xt;
    
    nablafxt = zeros([p,1]);
    for i=1:d
        Ai = A(:,:,i);
        yi = y(:,i);
        nablafxt = nablafxt + Ai'*(Ai*x0-yi)/d;
    end
    gt = nablafxt;
    nablafx0 = nablafxt;
    v0 = vt;
end

x_star = x0; 

%% ---------------------------------------------------
% MD-PDFP

T = 3000;
Re1 = MD_PDFP(x_star, 2.5/L, lam, T, A, y);
Re2 = MD_PDFP(x_star, 2/L, lam, T, A, y);
Re3 = MD_PDFP(x_star, 1/L, lam, T, A, y);
Re4 = MD_PDFP(x_star, 0.5/L, lam, T, A, y);

Re0 = Re1(1);

logRe1 = log(Re1/Re0)/log(10);
logRe2 = log(Re2/Re0)/log(10);
logRe3 = log(Re3/Re0)/log(10);
logRe4 = log(Re4/Re0)/log(10);

index = 200;

range = index:index:T;
plot(1:T,logRe1,'-x','Markersize',15,'LineWidth', 3,'MarkerIndices',range);
hold on;
plot(1:T,logRe2,'-o','Markersize',15,'LineWidth', 3,'MarkerIndices',range);
hold on;
plot(1:T,logRe3,'-^','Markersize',15,'LineWidth', 3,'MarkerIndices',range);
hold on;
plot(1:T,logRe4,'-s','Markersize',15,'LineWidth', 3,'MarkerIndices',range);
hold on;
xlabel('t');
ylabel('log_{10}(RE)')
ylim([-14 2]);
legend('\gammaL=2.5','\gammaL=2.0','\gammaL=1.0','\gammaL=0.5')
set(gcf,'position',[0,0,800,600])
ax = gca;
ax.FontSize = 25;




