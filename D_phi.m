phii = 40:1:50;
amoy = 160;
Npts = 1000;
fNormMin = 0.5;
fNormMax = 6;
N = 10;
kNormMin = 0.5;
kNormMax = 6;
P = 0;
Npics = 2;
alpha0 = [0 1];
filt = 0;
seuil = 8;
NormAxis = 0;
Search = 0;
CalcD = 1;
c0 = 1.48;
c1 = 0.64;
rho0 = 1;
rho1 = 1.85;

k1p = 3+5.*exp(-(phii - 16).^2./20);
k1p = cat(2,0.*k1p.',ones(length(phii),1));

k1p = [0.02 1];


for ii = 1 : length(phii)

    [~,~,~] = calc_disp_mat(amoy,phii(ii),Npts,fNormMin,fNormMax,...
        N,kNormMin,kNormMax,P,Npics,k1p,alpha0,filt,seuil,NormAxis,Search,CalcD, ...
        c0,c1,rho0,rho1);

end