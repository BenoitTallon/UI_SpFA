function [SpFct,cm,f] = Calc_disp2D(amoy,phi,Npts,fNormMin,fNormMax,...
    N,kNormMin,kNormMax,P,Npics,alpha1,alpha0,filt,seuil,NormAxis,Search,CalcD, ...
    c0,c1,rho0,rho1)

%% Calcul de c_eff avec le mod?le de P . Sheng
%
% (Ultrasonic wave transport in a system of disordered resonant scatterers:
% Propagating resonant modes and hybridization gaps, Cowan et al., 2011, PRB)
%
% Recherche des pics de la fonction spectrale calculee a partir de
% l'amplitude de l'onde diffus?e par le systeme sphere coque dans un milieu
% effectif.

% Nouvelle version (2) : utilisation de findpeaks.m (signal processing toolbox)
% pour la detection des pics. Possibilite de detecter plusieurs pics par
% frequence

%% Definitions

amoy = amoy/1000;
fNorm=linspace(fNormMin,fNormMax,Npts);                                    % Vecteur fr?quence normalisee
f = fNorm.*c0./(4*pi*amoy);                                                % Vecteur frequences (MHz)
om0=2*pi*f;                                                                % Pulsation
P = P./100;
phi = phi/100;

seuil = seuil./100;

%% calcul eta Polydisp

switch P
    
    case 0
        
        eta = phi/(4/3.*pi*amoy^3);                                        % Concentration de gouttes
        ra = amoy;                                                         % Rayon de gouttes d'huile (mm)
        rb = amoy/(phi.^(1/3));                                            % Rayon moyen sphere de matrice (mm)
        Nr = 1;
        
    otherwise
        
        sig=amoy*P;
        amin=amoy-3*sig;
        amax=amoy+3*sig;
        
        ra=amin:0.03*(amax-amin):amax;                                     % Rayon de gouttes d'huile (mm)
        rb = ra/(phi.^(1/3));                                              % Rayon moyen sphere de matrice (mm)
        
        %Distribution gaussienne
        
        gauss = exp(-0.5*((ra - amoy)/sig).^2);
        gaussn =( 1/sum(gauss)) * gauss;
        
        eta = gaussn.* phi / sum( gaussn .* ( (4/3)*pi*(ra.^3) ) );
        eta = repmat(eta,Npts,1);
        Nr = length(ra);
        
end
%% Proprietes des materiaux

k0 = om0./c0;

k0=real(k0) + 1i.*alpha0(1).*real(k0).^alpha0(2);                          % Nombre d'ondes matrice

if alpha0(2)~=1
    
    k0=real(k0) + 1i.*alpha0(1).*f.^alpha0(2);                             % Nombre d'ondes matrice
    
end

q0=k0./rho0;

k1 = om0./c1;

k1=real(k1) + 1i.*alpha1(1).*real(k1).^alpha1(2);                          % Nombre d'ondes inclusion

if alpha0(2)~=1
    
    k1=real(k1) + 1i.*alpha1(1).*f.^alpha1(2);                             % Nombre d'ondes matrice
    
end

q1=k1./rho1;


%% Milieu effectif

kNorm = linspace(kNormMin,kNormMax,Npts);                                  % Nombre d'ondes effectif normalise
keff = kNorm./(2*amoy);                                                    % Nombre d'ondes effectif (mm^-1)

rhoeff = (rho1 + rho0*(1/phi-1))*phi;                                      % Masse Vol. effective (loi des m?langes)

An=zeros(Npts,Npts,N);                                                     % Amplitude modale pression diff. par la sph?re de matrice
F0=zeros(Npts,Nr);                                                         % Fonction de diffusion vers l'avant
SpFct = zeros(Npts,Npts);                                                  % Self-energy
ReSelfEn = SpFct;
Att1 = SpFct;

ceff = zeros(Npics,Npts);


h = waitbar(0,'','Name','Calcul Amplitudes Modales ...',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

ll = 0;

%% Calcul de la fonction spectrale

%--------------------------------------------------------------------------
% Calcul de l'amplitude modale de l'onde diffus?e par la sph?re de matrice
% Calcul analytique du d?terminant pour optimisation
%
% [M] . [Amplitudes] = [C]
%
%--------------------------------------------------------------------------

tic

AAn = zeros(Npts,N,Npics);

qeff=keff./rhoeff;

for ir = 1 : Nr
    
    for k=0:N-1
        
        % Incrementation affichage tps calcul
        
        % Vecteur C ---------------------------------------------------
        
        C1 = -qeff.*besselSp(k,keff*rb(ir),1);
        C3 = -besselS(k,keff.*rb(ir),1) ;
        
        % Creation matrice M ------------------------------------------
        
        M11 = qeff.*besselSp(k,keff*rb(ir),2); M12 = -q0.*besselSp(k,k0*rb(ir),2);
        M13 = -q0.*besselSp(k,k0*rb(ir),3);  M22 = q0.*besselSp(k,k0*ra(ir),2);
        M23 = q0.*besselSp(k,k0*ra(ir),3); M24 = -q1.*besselSp(k,k1*ra(ir),1);
        M31 = besselS(k,keff*rb(ir),2); M32 = -besselS(k,k0*rb(ir),2);
        M33 = -besselS(k,k0*rb(ir),3); M42 = besselS(k,k0*ra(ir),2);
        M43 = besselS(k,k0*ra(ir),3); M44 = -besselS(k,k1*ra(ir),1);
        
        for jj = 1 : Npts
            
            ll = ll + 1;
            
            D1 = C1(jj).*(M24.*(M32.*M43-M42.*M33)+M44.*(M22.*M33-M32.*M23)) + ...
                C3(jj).*(M44.*(M12.*M23-M22.*M13)-M24.*(M12.*M43-M42.*M13));
            
            D2 = M11(jj).*(M24.*(M32.*M43-M42.*M33)+M44.*(M22.*M33-M32.*M23)) + ...
                M31(jj).*(M44.*(M12.*M23-M22.*M13)-M24.*(M12.*M43-M42.*M13));
            
            An(:,jj,k+1) = D1./D2;                                         % Cramer
            
            waitbar(ll/((N-1)*Nr*Npts),h,' ')
            
            if getappdata(h,'canceling')
                
                delete(h)
                
                break
                
            end
            
        end
        
    end
end

clear('M11','M13','M23','M31','M33','M43','C1','C3','qeff','D1','D2',...
    'M12','M22','M24','M32','M42','M44')

ll = 0;

set(h,'Name','Calcul Fct Spectrale')

for ir = 1 : Nr
    
    for jj = 1 : Npts
        
        ll = ll + 1;
        
        NN = (2.*(0:N-1) + 1);
        
        F0(:,ir) = sum(repmat(NN,Npts,1).*squeeze(An(:,jj,:)),2)./(1i.*keff(jj).');
        
        SpFct(:,jj) = -mean(imag(1./(4.*pi.*eta.*F0)),2);                  % Calul de la fonction spectrale avec l'eq. (3.10) de
        %                                                                  %"Acoustic and electromagnetic quasimodes in dispersed random media"
        %                                                                  % (X. Jing et al PRA 1992)
        %                                                                  % Moyenne de la fonction spectrale sur les rayons de sphere
        
        Att1(:,jj) = mean(imag(2.*pi.*eta.*F0)./real(k0).',2);
        
        waitbar(ll/(Npts*Nr),h,' ')
        
        if getappdata(h,'canceling')
            
            delete(h)
            
            break
            
        end
        
    end
    
end

toc

delete(h)

alphaeff = zeros(Npics,Npts);

if Search == 1
    
    switch Npics
        
        case 1
            
            for ii = 1 : Npts
                
                try
                    
                    ceff(ii) = find(SpFct(:,ii) == max(SpFct(:,ii)));
                    
                catch
                    
                    ceff(ii) = ceff(ii-1);
                    
                end
                
                AAn(ii,:) = An(ii,ceff(:,ii),:);
                alphaeff(ii) = Att1(max(ceff(:,ii)),ii);
                
            end
            
        otherwise
            
            for ii = 1 : Npts
                
                [xx,yy] = findpks(SpFct.');
                
                try
                    
                    [~,xxi] = sort(xx,'descend');
                    xxi = xxi(1 : Npics);
                    ceff(:,ii) = yy(xxi);
                    
                catch
                    
                    ceff(1:length(yy),ii) = yy;
                    ceff(length(yy)+1:end,ii) = 1;
                    
                end
                
                for jj = 1 : Npics
                    
                    AAn(ii,:,jj) = An(ii,ceff(jj,ii),:) ;
                    
                end
                
                alphaeff(:,ii) = Att1(ceff(:,ii),ii);
                
            end
            
    end
    
else
    
    switch Npics
        
        case 1
            
            for ii = 1 : Npts
                
                try
                    
                    ceff(ii) = find(SpFct(ii,:) == max(SpFct(ii,:)));
                    
                catch
                    
                    ceff(ii) = ceff(ii-1);
                    
                end
                
                AAn(ii,:) = An(ii,ceff(:,ii),:);
                alphaeff(ii) = Att1(ii,max(ceff(:,ii)));
                
            end
            
        otherwise
            
            for ii = 1 : Npts
                
                [xx,yy] = findpks(SpFct(ii,:));
                
                try
                    
                    [~,xxi] = sort(xx,'descend');
                    xxi = xxi(1 : Npics);
                    ceff(:,ii) = yy(xxi);
                    
                catch
                    
                    ceff(1:length(yy),ii) = yy;
                    ceff(length(yy)+1:end,ii) = 1;
                    
                end
                
                for jj = 1 : Npics
                    
                    AAn(ii,:,jj) = An(ii,ceff(jj,ii),:);
                    
                end
                
                alphaeff(:,ii) = Att1(ii,ceff(:,ii));
                
            end
            
    end
    
end

keff(keff == keff(1)) = NaN;
alphaeff(isnan(alphaeff)) = 0;

%%

if Search == 1
    
    cm = om0(ceff)./repmat(k0,Npics,1);
    
else
    
    cm = (repmat(om0,Npics,1)./keff(ceff));
    
end

keff = keff(ceff);

%% Filtrage courbe de dispersion

for jj = 1 : size(keff,1)
    
    keff(jj,:)  = smooth2(keff(jj,:));
    
end

%% ---------------------- Calcul vitesse energie --------------------------

if CalcD == 1
    
    for jj = 1 : Npics
        
        [cosMoy , d2] = calcul_cosMoy2(AAn(:,:,jj),f,keff(jj,:));
        
        ltr = 1./(2.*alphaeff(jj,:).*(1 - cosMoy.'));
        
        %------------------- Calcul vgr en terme de delais ----------------
        
        F0 = foncDiff(AAn(:,:,jj).',0,keff(jj,:));
        
        dgr = [2*pi*eta*c0^2./(2*pi*f(1:Npts-1)).*diff(real(F0))./diff(2.*pi.*f) 0];
        
        vgr(jj,:) = (c0^2./cm(jj,:))./(1 + dgr);
        
        %------------------------------------------------------------------% Vitesse d'energie (J. Page et P. Sheng)
        
        d1 = [2*pi*eta*cm(jj,1:Npts-1).*vgr(jj,1:Npts-1)./(2*pi*f(1:Npts-1)).*diff(real(F0))./diff(2.*pi.*f) 0];
        
        d2 = [2*pi*eta*cm(jj,1:Npts-1).*vgr(jj,1:Npts-1).*d2.' 0];
        
        vt(jj,:) = (c0^2./cm(jj,:))./(1 + d1 + d2);
        
        %------------------------------------------------------------------% Coefficient de diffusion
        
        D(jj,:) = vt(jj,:).*ltr./3;
        
    end
    
end

disp(D(:,find(f>=2.5,1)))

%--------------------------------------------------------------------------
%% Figures

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot(f,cm,'o')

if CalcD ==1
    
    hold on
    plot(f,vt,'or')
    hold on
    plot(f,vgr,'ok')
    
    legend('v_m','v_m','v_e','v_{gr}')
    
end

xlabel(axes1,'Frequency (MHz)')
ylabel(axes1,'c_{eff} (mm/\mus)')

box(axes1,'on');
set(axes1,'Layer','top',...
    'fontsize',15);

figure2 = figure;

axes2 = axes('Parent',figure2);
hold(axes2,'on');

if NormAxis == 1
    
    im1 = imagesc(real(kNorm),real(fNorm),SpFct);
    
else
    
    im1 = imagesc(real(2.*pi.*f./c0),f,SpFct);
    
end

box(axes2,'on');
set(axes2,'Layer','top',...
    'fontsize',15);

set(axes2,'YAxisLocation','left');
set(axes2,'XAxisLocation','bottom');

caxis(axes2,[mean(alphaeff(1,:))*0.1 mean(alphaeff(1,:))*30])

if NormAxis == 1
    
    xlim(axes2,[kNormMin kNormMax])
    ylim(axes2,[fNormMin fNormMax])
    xlabel(axes2,'2ka')
    ylabel(axes2,'2\omegaa/c_0')
    
    hold on
    
    plot(real(fNorm),real(fNorm),'--r','linewidth',1.5)
    
    hold on
    
    plot(real(keff.*2*amoy),real(fNorm),'.w','markersize',8)
    
else
    
    xlim(axes2,[min(real(2.*pi.*f./c0)) max(real(2.*pi.*f./c0))])
    ylim(axes2,[min(f) max(f)])
    xlabel(axes2,'k (mm^{-1})')
    ylabel(axes2,'f (MHz)')
    
    hold on
    
    plot(real(2.*pi.*f./c0),f,'--r','linewidth',1.5)
    
    hold on
    
    if Search == 1
        
        plot(real(2.*pi.*f./c0),f(ceff),'.w','markersize',15)
        
    else
        
        plot(real(keff),real(f),'.w','markersize',15)
        
    end
    
end

figure3 = figure;
axes4 = axes('Parent',figure3);
hold(axes4,'on');
plot(f,alphaeff,'.')

xlabel(axes4,'Frequency (MHz)')
ylabel(axes4,'\alpha_{eff} (mm^{-1})')

box(axes4,'on');
set(axes4,'Layer','top',...
    'fontsize',15);

if CalcD == 1
    
    figure4 = figure;
    axes5 = axes('Parent',figure4);
    hold(axes5,'on');
    plot(f,D,'o')
    
    xlabel(axes5,'Frequency (MHz)')
    ylabel(axes5,'D (mm^2/\mus)')
    
    box(axes5,'on');
    set(axes5,'Layer','top',...
        'fontsize',15);
    
end

%% --------------------- Fonctions locales --------------------------------

function [pks, loc]=findpks(x)

I1=((x(1:end-1)-x(2:end))>0);
I2=((x(2:end)-x(1:end-1))>0);

I1=[I1,1];
I2=[1,I2];

loc = find((I1.*I2));
pks = x(loc);



function [cosMoy , dm] = calcul_cosMoy2(A,f,k)

dtheta = 0.01;
theta = 0 : dtheta : pi ;
Ftheta = zeros(size(A,1),length(theta));

for ii = 1 : length(theta)
    
    Ftheta(:,ii) = foncDiff(A.',theta(ii),k);
    
end

Ftheta(isnan(Ftheta)) = 0;

om2 = repmat(2*pi*f.',[1,ii]);
theta2 = repmat(theta,[length(f),1]);

cosMoy = sum(sin(theta2).*dtheta.*cos(theta2).*...
    (Ftheta.*conj(Ftheta)),2);

cosMoy = cosMoy./sum(sin(theta2).*dtheta.*...
    (Ftheta.*conj(Ftheta)),2);


theta2 = repmat(theta,[length(f)-1,1]);

dm = sum(sin(theta2).*dtheta.*(diff(unwrap(angle(Ftheta)))./diff(om2)).*...
    (Ftheta(1:length(f) - 1,:).*conj(Ftheta(1:length(f) - 1,:))),2);

function fTheta=foncDiff(an,theta,k)

N=size(an,1);
fTheta=0;

for n=0:N-1
    
    fTheta=(2.*n+1).*an(n+1,:).*legN(n,cos(theta))+fTheta;
    
end

fTheta=fTheta./(1i.*k);


function P=legN(n,x)

N=length(x);


P=legendre(n,x);

if n>0
    
    P=reshape(P(1,:,:),1,N); %Suppression des tableaux m>0 (sauf pour n=0)
    
end


function w=besselS(n,z,u)

if u==1
    w=sqrt(pi./(2.*z)).*besselj(n+0.5,z);
elseif u==2
    w=sqrt(pi./(2.*z)).*besselh(n+0.5,1,z);
else
    w=sqrt(pi./(2.*z)).*besselh(n+0.5,2,z);
end


function w=besselSp(n,x,u) %D�riv�es

if u==1
    w=sqrt(pi./(2.*x)).*((n+0.5).*besselj(n+0.5,x)./x-besselj(n+1.5,x)-besselj(n+0.5,x)./(2.*x));
elseif u==2
    w=sqrt(pi./(2.*x)).*((n+0.5).*besselh(n+0.5,1,x)./x-besselh(n+1.5,1,x)-besselh(n+0.5,1,x)./(2.*x));
    %w=0.5.*(besselh(n-1,1,x)-(besselh(n,1,x)+x.*besselh(n+1,1,x))./x);
else
    w=sqrt(pi./(2.*x)).*((n+0.5).*besselh(n+0.5,2,x)./x-besselh(n+1.5,2,x)-besselh(n+0.5,2,x)./(2.*x));
    %w=0.5.*(besselh(n-1,2,x)-(besselh(n,2,x)+x.*besselh(n+1,2,x))./x);
end