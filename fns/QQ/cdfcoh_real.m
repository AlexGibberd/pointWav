% Script to look at coherence in real setting

Delta = 1e-2;
X = Delta:Delta:1-Delta;
dof = 5;
Y=[]; 

% [Y] = emp_cdf(X,dof,0)

for i=1:length(X)
    y = dist_R2(X(i),dof,0.5);
    Y = [Y;y];
end

Y2 = trapz(X,Y');
 
figure
plot(X,Y,'Color','r','LineWidth',2);
hold on
plot(X,Y2,'Color','r','LineWidth',2);

function [y] = dist_R2(x,dof,r2)

n=dof-1;
a = n/2;
b = n/2;
c = 1/2;
ag = [a,b];
bg = [c];
A = gamma(n/2)/(gamma(1/2)*gamma((n-1)/2));
B = (1-r2)^(n/2);
y = A*(x^(-1/2))*((1-x)^((n-3)/2))...
    *B*hypergeom(ag,bg,r2*x);

end

function [y] = cdf_R2(x,dof,r2)

n = dof-1
nb = n/2 - 1;
a = n/2;
b = n/2;
c = 1/2;
ag = [a,b];
bg = [c];
v = hypergeom(ag,bg,r2*x);
a = nb;
b = nb;
c = -1/2;
ag = [a,b];
bg = [c];
vp = (1/(2*(nb^2)*r2))*hypergeom(ag,bg,r2*x);
A = gamma(n/2)/(gamma(1/2)*gamma((n-1)/2));
B = (1-r2)^(n/2);

u = (x^(-1/2))*((1-x)^((n-3)/2));
up = -(1/2)*(x^(-3/2))*((1-x)^((n-3)/2))...
    - x^(-1/2)*((n-3)/2)*((1-x)^((n-5)/2));

y = A*B*(u-up)*vp;

end

function [Y] = emp_cdf(X,dof,r2)
% X - vector of points in (0,1)
dx = diff(X);
Y(1) = 0;
% Midpoint rule
for i=2:length(X)
   Y(i) = dx(i-1)*dist_R2((X(i)+X(i-1))/2,dof,r2) + Y(i-1);
end

end
