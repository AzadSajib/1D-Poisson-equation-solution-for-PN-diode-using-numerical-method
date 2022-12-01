clc; clear all; close all

% Parameters
Na= 1e16 * 1e6;   % Acceptor concentration
Nd= 1e15 * 1e6;   % Donor Concentration
q= 1.6e-19; 
ni= 1.5e10 * 1e6; % Si intrinsic concetration

Vt= 0.0258;       % Thermal voltage
epsilon= 11.7* 8.854e-12;
tox= 1e-9;        % oxide thickness
Vg= 1.2875;       % Gate voltage

Vbi= Vt*log((Na*Nd)/(ni^2)); 


N= input ('Number of points: ');

xn = (((2*epsilon*Vbi)/q)*(Na/Nd)*(1/(Na+Nd)))^(1/2);
xp = (((2*epsilon*Vbi)/q)*(Nd/Na)*(1/(Na+Nd)))^(1/2);
W= xp + xn; 

x= linspace(-xp, xn, N);
del_x= x(2)-x(1);


%% Finding the material boundary

x1= find(x<=0);
x2= find(x>0);


%% Energy matrix
D=[];
for i=1:N
    for j=1:N
        if(i==j)
            D(i,j)=-2;
        elseif (abs(i-j)==1)
            D(i,j)=1;
        else
            D(i,j)=0;
        end
            
    end
end

D1= D(2: N-1, :);
D2= D1/del_x^2;

[a,b]= size(D2);
c=[1, zeros(1,b-1)];
d=[zeros(1,b-1), 1];

D3= [c ; D2; d];

%% Charge Density Matrix

rho= [];

rho_p(1:length(x1))= -q*Na;
rho_n(1:length(x2))= q*Nd;

rho= [rho_p, rho_n];
rho= -rho/epsilon; 

Rho= [0, rho(2:end-1), Vbi]';

subplot(211)
plot(linspace(-xp, xn, N)', Rho, 'linewidth', 2);
title('Chareg Density Profile'), xlabel('Depth'), ylabel('Charge Density');
grid on;

%% Potential Calculation 

Potential= linsolve(D3, Rho); 

subplot(212)
plot(linspace(-xp, xn, N), Potential, 'linewidth', 2);
title('Potential Profile'), xlabel('Depth'), ylabel('Potential')
grid on;





