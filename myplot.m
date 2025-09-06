%
% Weiwei Fu's code
%
load k_n.mat
figure(10);
% contourLevels = [1.1, 1.2, 1.3, 1.4, 1.5];
contourLevels = [9,10,11,12,13,14,16,20,25,30,35,40,45,50];
[C, h] = contour(X, k, Z, contourLevels, 'k'); 
clabel(C, h, 'FontSize', 12,'LabelSpacing',300); % Label contours
hold on; 
% Labels and title
xlabel('Exponent of windspeed dependence');
ylabel('Parameter k (cm/hr)');
title('GLODAPv2023');
set(gca,'fontsize',14);
set(gcf,'color','w');
grid on;

%
% My plot
%

%
% first rename Fu's Z to Y, that way we can use Z for the
% normalization constant
%
Y  = Z;
f = Y.*Y;

% assume:
%     that f = (d - m)' * W * ( d - m)
%     where d is the GLODAPV2_2023 Delta14C
%     and m(k,n) is the model computed Delta14C

%
% First we regrid Y2 onto a mesh with uniform spacing in both n AND k
%
N = 0.1:0.01:4; dn = max(N)-min(N);
K = 4:0.01:30;  dk = max(K)-min(K);



X = X(:); % the exponents
k = k(:); % the piston velocities
dsites = [X/dn,k/dk];
ctrs = dsites;


[X,Y] = meshgrid(N,K);
epoints = [X(:)/dn,Y(:)/dk];


rbf = @(ep,r) sqrt((1 + ep*r.*r)); ep =10^9;
%rbf = @(e,r) exp(-(e*r).^2); ep = 2.5;
IM = rbf(ep,DistanceMatrix(dsites,ctrs));
EM = rbf(ep,DistanceMatrix(epoints,ctrs));

Pf = X*0;
Pf(:) = EM*(IM\sqrt(f(:)));
figure(1)

[C, h] = contour(X, Y, Pf, contourLevels, 'k'); 
clabel(C, h, 'FontSize', 12,'LabelSpacing',300); % Label contours
hold on
scatter(dsites(:,1)*dn,dsites(:,2)*dk,'o',...
        'MarkerEdgeColor',[0.5,0.5,0.5],....
        'MarkerFaceColor',[0.5 0.5 0.5],...
        'MarkerFaceAlpha',0.2);
set(gca,'FontSize',16);
xlabel('$n$','Interpreter','latex')
ylabel('$\left<k_w\right>/\mbox{(cm/hr)}$','Interpreter','latex')
grid on

%
% Normalize the pdf and find the marginal posteriors
%
dn = X(1,2)-X(1,1);
dk = Y(2,1)-Y(1,1);
dA = dn*dk;
s2 = min(Pf(:).*Pf(:)); % optimal

Z = sum( exp(-0.5*(Pf(:).*Pf(:)/s2)).*dA );
P = 0*X;
P(:) = exp(-0.5*(Pf(:).*Pf(:)/s2))/Z;


% marginalize out n
Pkw = sum(P*dn,2);

% marginalize out kw
Pn = sum(P*dk,1);

figure(2)
H = plot(X(1,:),Pn,'-k');
set(H,'LineWidth',2);
xlabel('$n$','Interpreter','latex');
ylabel('probability density');
set(gca,'FontSize',16)
grid on      


figure(3)
H = plot(Y(:,1),Pkw,'-k');
set(H,'LineWidth',2);
xlabel('$\left<k_w\right>$/(cm/hr)','Interpreter','latex');
ylabel('probability density (cm/hr)$^{-1}$','Interpreter','latex');
grid on
set(gca,'FontSize',16)

