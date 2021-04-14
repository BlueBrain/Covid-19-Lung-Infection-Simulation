%{
Covid-19 Lung Infection Simulator (C19-LIS)

MIT License

Copyright (c) [2021] Blue Brain Project/EPFL Covid-19 Lung Infection Simulator

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%BBP SARS-CoV-2 infection glucose-dependence model for Logette et al., 2021

 
close all;
clear all;


dt=0.1;



%%%%% Time sweeps (optional) %%%%%%%
t = 0:5000000:10000000000 ; % msecs
tmsecs = 0:5000000:1000000000 ; % msecs
tsecs = 0:50000:10000000      ; % secs
tmins = 0:833.3333:166666.667          ; % mins
thrs =  0:13.8888:2777.7777            ;% hrs
tdays = 0:0.5787:115.7407            ; % days


%%%%%%% Independent variable sweeps %%%%%%%

gluc_b_range = 0:0.15:30  ; % mM, range of blood glucose, includes hyperglycemic levels > 5 mM
lectin_range = 0:0.05:10    ; % mM, range of lectin concentration in ASL
Rt_range = 0:1:100       ; % ohm*cm^2, paracellular resistivity, corresponds to glucose permeability (Garnett Baines 2012)

viral_load = 0.0: 0.000001 : 0.0002 ;% for range sims only

%%%%%%% Fixed or Initial Parameters %%%%%%%

gluc_b = 5    ; % initial glucose concentration in the blood normally
gluc_b_hyp= 10; % hyperglycemic blood concentration
gluc_epi_init = 2 ;
gluc_asl_init = 0.4;
gluc_asl_h_init = 1.2 ; %
lac_b = 1     ; % mM, initial lactate concentration in the blood
lac_epi = 1   ; % initial lactate concentration in the
lac_epi_init = 1 ;
lac_asl = 1  ; % ? initial lactate concentration in the ASL
lac_asl_init = 1 ;
lectin = 5.8  ; % mM, initial C-type lectin, or D (ie SP-D?), conc. in ASL; based on Guansbaek 2013, see calculations below
Rt_n = 453    ; % tight junction permeability to glucose in normal gluc population, ohm*cm^2 paracellular resistance (Garnett Baines 2012)
Rt_h = 225    ; % tight junction permeability to glucose in hyperglycemic population, ohm*cm^2 paracellular resistance (Garnett Baines 2012)

pH_asl_init = 7.0   ; % ? initial pH in the ASL
pH_epi_init = 7.2  ; % pH in the epithelial cells
pH_b   = 7.4  ; % blood pH

KmACE2 = .000016 ;  % mM; ACE2 receptor Km for SARS-CoV-2 binding and entry into apical airway epithelium
KmGlut2 = 17  ; % glucose transporter type 2 Km for basolateral airway epithelium (Baker 2018)
KmGlut10 = 0.3; % glucose transporter type 10 Km for apical airway epithelium (Baker 2018)
KmGlut1 = 3   ; % glucose transporter type 1 Km for basolateral airway epithelium (Baker 2018)
KmSglt1 = 0.3 ; %  Na-coupled glucose transporter Km for apical distal alveolar epithelium
KmGlutX =  17   ; % glucose transporter type X Km for basolateral distal alveolar epithelium

KmLectg = 1.8 ; % mM, Km for lectin binding to glucose (con. A as lectin; Schwartz 1993; see also Kussrow 2009)
KmLectv = 0.45; % Km for lectin binding to c19 virus (mannose as a substitute for virus; Schwartz 1992; see also Kussrow 2009)
KmMCTapi = 1.7; % mM, Km for apical endothelial lactate transporter (MCT4, Contreras-Baeza 2019)
KmMCTbas = 1.7; % Km for basolateral endothelial lactate transporter ? 
MCTapi = 1    ; % (placeholder) concentration of apical MCT transporters (combine MCTs 2 and 4 for simplicity)
MCTbas = 1    ; % (placeholder) concentration of basolateral MCT transporter
DifuseK = 5.65e3; %5.65e-9; % m^2/s difusion constant for glucose 

kLDHeminus = 0.0710 ; %
kLDHeplus = 1.5900  ; %

virus_asl = .00001  ; % initial number of virions in ASL
virus_aslb = .00005  ;
virus_aslc = .00010  ;

%%%%%%%%%%%%%%%% calculations section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  for lectin, based on Guansbaek (2013)
%  mol. weight = 43kDa ~ 43000 g/mol
%  estimated 500 ng/ml in BAL = 5e-4 g/L (Emm's calculation is 500ng/ml,
%  from paper was 10x higher 5000 ng/ml

% (1 mol/ 43000 g) * (5e-4 g/L) = 1.16e-8 M
% estimated dilution in BAL = 10-100 x ; take 500x as estimate
% 500 * 1.16e-8 M = 580e-8 M = 5.8 mM

% 1mg/ml glucose = 1g/L * 1 mol/180g = 5 mM; so 4 mg/ml = 20 mM glucose incubate and sig rise in 
% viral replication as per Kohio (2013) - 2.5 fold increase in # of
% infected cells per 2 mM increase in [glucose]

% viral replication log10 increase per 10 mM glucose (Reading 2018)

% viral replication from 10^2 to 10^4 virions/ml in 72 hrs (Yu 2011)


%%%%%%%%%%%%%%%%%%%% Matrices  %%%%%%%%%%%%%%%%%%%%

c19m=[]  ; %
ACE2Am=[]; %



%%%%%% Rates %%%%%%%%

tau1 = 10   ; % secs, arrival rise-time tau of c19 bolus into the ASL from sneeze
tau1hrs= 0.0027 ; % hrs
tau2 = 150000 ; %85000; % secs, clearance-time tau of c19 bolus from sneeze, reflects mucocilliary clearance 
tau2hrs =  41.6666  ;%
 

%%%%%%%%%%% Governing equations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                               
c19 = exp(-tsecs/tau2) - exp(-tsecs/tau1)      ; % SARS-CoV-2 virion concentration in ASL
%c19mins = exp(-tmins/tau2) - exp(-tmins/tau1) ;
c19hrs = exp(-thrs/tau2) - exp(-thrs/tau1)      ; % SARS-CoV-2 virion concentration in ASL 

virus_asl_effective_n =  virus_asl .* (0.4*lectin ./(KmLectg + lectin));  % 
virus_asl_effective_h =  virus_asl .* (1.2*lectin ./(KmLectg + lectin));  % 
ACE2V = c19*(virus_asl_effective_n./(KmACE2 + virus_asl_effective_n))  ;
ACE2V_n = (c19)*(virus_asl_effective_n./(KmACE2 + virus_asl_effective_n)) ;
ACE2V_h = (c19)*(virus_asl_effective_h./(KmACE2 + virus_asl_effective_h)) ;    


virus_asl_effective_nb =  virus_aslb .* (0.4*lectin ./(KmLectg + lectin));  % 
virus_asl_effective_hb =  virus_aslb .* (1.2*lectin ./(KmLectg + lectin));  % 
ACE2Vb = c19*(virus_asl_effective_nb./(KmACE2 + virus_asl_effective_nb))  ;
ACE2V_nb = (c19)*(virus_asl_effective_nb./(KmACE2 + virus_asl_effective_nb)) ;
ACE2V_hb = (c19)*(virus_asl_effective_hb./(KmACE2 + virus_asl_effective_hb)) ;    


virus_asl_effective_nc =  virus_aslc .* (0.4*lectin ./(KmLectg + lectin));  % 
virus_asl_effective_hc =  virus_aslc .* (1.2*lectin ./(KmLectg + lectin));  % 
ACE2Vc = c19*(virus_asl_effective_nc./(KmACE2 + virus_asl_effective_nc))  ;
ACE2V_nc = (c19)*(virus_asl_effective_nc./(KmACE2 + virus_asl_effective_nc)) ;
ACE2V_hc = (c19)*(virus_asl_effective_hc./(KmACE2 + virus_asl_effective_hc)) ;    

virus_asl_effective_nl =  viral_load .* (0.4*lectin ./(KmLectg + lectin));  % 
virus_asl_effective_hl =  viral_load .* (1.2*lectin ./(KmLectg + lectin));  % 
ACE2Vl = c19.*(virus_asl_effective_nl./(KmACE2 + virus_asl_effective_nl));
ACE2V_nl = c19.*(virus_asl_effective_nl./(KmACE2 + virus_asl_effective_nl)) ;
ACE2V_hl = c19.*(virus_asl_effective_hl./(KmACE2 + virus_asl_effective_h)) ;    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Glut1 = (gluc_b ./(KmGlut1 + gluc_b))   ; % activity of basolateral airway endothelial glucose transporter 1 
GlutX = (gluc_b ./(KmGlutX + gluc_b))      ; % activity of basolateral distal alveolar endothelial glucose transporter X

MCTapifrac = MCTapi * (lac_epi ./(KmMCTapi + lac_epi))  ; % activity of apical endothelial lactate transporter
MCTbasfrac = MCTbas * (lac_b ./(KmMCTbas + lac_b))      ; % activity of apical endothelial lactate transporter,

lectin_g =  5.8 * (lectin ./(KmLectg + lectin))    ; % fraction lectin bound to glucose in the ASL
lectin_v = 5.8 *lectin ./(KmLectv + lectin)  ; % lectin bound to virus in the ASL, 3 lectin bind to 1 mannose?
lectin_free = 5.8 - lectin_v - lectin_g ;

Jgluc_nr = -(gluc_asl_init-gluc_b_range)*DifuseK*1/Rt_n   ; % Fick mod for gluc blood --> gluc ASL via Rt normal gluc_b
Jgluc_hr = -(gluc_asl_h_init-gluc_b_range)*DifuseK*1/Rt_h   ; % Fick mod for gluc blood --> gluc ASL via Rt hyper gluc_b
Jgluc_n = -(gluc_asl_init-gluc_b)*DifuseK*1./Rt_n ;
Jgluc_h = -(gluc_asl_h_init-gluc_b_hyp)*DifuseK*1./Rt_h ;
Jgluc_nt = -(gluc_asl_init-gluc_b)*DifuseK*1./Rt_range ;
Jgluc_ht = -(gluc_asl_h_init-gluc_b_hyp)*DifuseK*1./Rt_range ;
Jgluc_hnd = Jgluc_hr - Jgluc_nr ;

junct_conduct = 1./Rt_range;

gluc_asl = gluc_asl_init + Jgluc_n ;%  
gluc_asl_h= gluc_asl_h_init + Jgluc_h ;%  
gluc_asl_r = gluc_asl_init + Jgluc_nr ;%  
gluc_asl_hr= gluc_asl_h_init + Jgluc_hr ;% 


ACE2V_r = ACE2V_n .* (gluc_asl_r);
ACE2V_rb = ACE2V_nb .* (gluc_asl_r);
ACE2V_rc = ACE2V_nc .* (gluc_asl_r);
ACE2V_rl = ACE2V_nl .* (gluc_asl_r);

ACE2V_hr = ACE2V_h .* (gluc_asl_r);
ACE2V_hrb = ACE2V_hb .* (gluc_asl_r);
ACE2V_hrc = ACE2V_hc .* (gluc_asl_r);
ACE2V_hrl = ACE2V_hl .* (gluc_asl_r);


%Glut10 = (gluc_asl_r ./(KmGlut10 + gluc_asl_r)) ; % activity of apical airway endothelial glucose transporter 10
%Glut2 = (gluc_asl_r ./(KmGlut2 + gluc_asl_r))   ; % activity of apical airway endothelial glucose transporter 2
Sglt1 = (gluc_asl_r ./(KmSglt1 + gluc_asl_r))   ; % activity of apical distal alveolar epithelium Na-coupled glucose transporter

gluc_epi = gluc_epi_init + Sglt1*(gluc_asl) ;
gluc_epi_h = gluc_epi_init + Sglt1*(gluc_asl_h) ;
gluc_epi_r = gluc_epi_init + Sglt1.*(gluc_asl_r) ;
gluc_epi_hr = gluc_epi_init + Sglt1.*(gluc_asl_hr) ;


lac_epi = 2*(kLDHeplus * gluc_epi_r - kLDHeminus * lac_epi_init); %
lac_asl = MCTapifrac * lac_epi     ; % ASL lactate concentration change
lac_epi_h = 2*(kLDHeplus * gluc_epi_hr - kLDHeminus * lac_epi_init); %
lac_asl_h = MCTapifrac * lac_epi_h    ; % ASL lactate concentration change

pH_asl = pH_asl_init -log10(lac_asl)  ; % + -log...
pH_epi = pH_epi_init +log10(lac_epi)  ; % - -log...


vir_epi_pre = ACE2V_r ; % 
vireprate_epi = ACE2V_n .*gluc_asl_r  ; % units 
vinumber_epi = ACE2V_r.*(2.^(vireprate_epi.*gluc_asl_r)).*tsecs; % 
vinumber_log = log10(vinumber_epi) ; %

vir_epi_preb = ACE2V_rb ; % 
vireprate_epib = ACE2V_nb .*gluc_asl_r  ; % 
vinumber_epib = ACE2V_rb.*(2.^(vireprate_epib.*gluc_asl_r)).*tsecs; % 
vinumber_logb = log10(vinumber_epib) ; %

vir_epi_prec = ACE2V_rc ; % 
vireprate_epic = ACE2V_nc .*gluc_asl_r  ; %  
vinumber_epic = ACE2V_rc.*(2.^(vireprate_epic.*gluc_asl_r)).*tsecs; % 
vinumber_logc = log10(vinumber_epic) ; %

vir_epi_prel = ACE2V_rl ; % 
vireprate_epil = ACE2V_nl .*gluc_asl_r  ; %  
vinumber_epil = ACE2V_rl.*(2.^(vireprate_epil.*gluc_asl_r)).*tsecs; % 
vinumber_logl = log10(vinumber_epil) ; %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1)   

plot(gluc_b_range,Jgluc_nr./100,'color', [0.4 0.4 0.4],'linewidth',3)
hold on
plot(gluc_b_range,Jgluc_hr./100,'k','linewidth',3)
hold off
legend('normoglycemic paracelluar permeability', 'hyperglycemic paracelluar permeability', 'Location', 'northwest')
xlim([3 15])
xlabel('blood glucose (mM)')
ylabel('glucose flux (mM/\Omegasec)')
%title('')


figure(2)   

plot(junct_conduct,Jgluc_nt./100,'color', [0.4 0.4 0.4],'linewidth',3)
hold on
plot(junct_conduct,Jgluc_ht./100,'k','linewidth',3)
legend('normoglycemic blood glucose', 'hyperglycemic blood glucose', 'Location', 'northwest')
xlabel('1/Rt (\Omega^{-1}cm^{-2})')
ylabel('glucose flux (mM/\Omegasec)')
%ylabel('paracellular glucose flux ((mM)(\mu^{2}))/(msec)(Rt))')
hold off
%title('')



figure (3) % 

plot(gluc_b_range,(0.01.*gluc_asl_r),'color', [0.4 0.4 0.4],'linewidth',3)
hold on
plot(gluc_b_range,(0.01.*gluc_asl_hr),'k','linewidth',3)
hold off
legend('normal Rt', 'impaired Rt','Location', 'northwest')
xlim([3 15])
xlabel('blood glucose (mM)')
ylabel('ASL glucose (mM)')
%title('')




figure (5)

plot(gluc_asl_r, ACE2V_r,'color', [0.7 0.7 0.7],'linewidth',3)
hold on
plot(gluc_asl_r, ACE2V_rb,'color',[0.4 0.4 0.4],'linewidth',3)
plot(gluc_asl_r, ACE2V_rc,'Color',[0 0 0],'linewidth',3)
hold off
legend('low viral load', 'intermediate viral load','high viral load','Location', 'northwest')
xlabel('ASL glucose (mM)')
ylabel('SARS-CoV-2 endocytsed ((units)(mM)/sec)')
xlim([0.2 2])
%title('')



figure (10)

plot(tsecs/3600,c19./0.7,'k','linewidth',3)
xlabel('hours')
ylabel('normalized ASL viral content')
xlim([0 240])
ylim([0 1])
title('SARS-CoV-2 "sneeze" stimulus and clearance')
 


figure(11)    %

subplot(1,3,1)
plot(tsecs/3600,ACE2V_n,'color', [0.4 0.4 0.4],'linewidth',3)
hold on
plot(tsecs/3600,ACE2V_h,'--k','linewidth',3)
%xlim([0 1e6])
xlim([0 240])
ylim([0 0.7])
xlabel('hours')
ylabel('SARS-CoV-2-ACE2 ((units)/(sec))')
legend('norm ASL gluc', 'hyper ASL gluc')
axis square
title('low viral load')
hold off

subplot(1,3,2)
plot(tsecs/3600,ACE2V_nb,'color',[0.4 0.4 0.4],'linewidth',3)
hold on
plot(tsecs/3600,ACE2V_hb,'--k','linewidth',3)
%xlim([0 1e6])
xlim([0 240])
ylim([0 0.7])
xlabel('hours')
ylabel('SARS-CoV-2-ACE2 ((units)/(sec))')
legend('norm ASL gluc', 'hyper ASL gluc')
axis square
title('intermediate viral load')
hold off

subplot(1,3,3)
plot(tsecs/3600,ACE2V_nc,'color',[0.4 0.4 0.4],'linewidth',3)
hold on
plot(tsecs/3600,ACE2V_hc,'--k','linewidth',3)
%xlim([0 1e6])
xlim([0 240])
ylim([0 0.7])
xlabel('hours')
ylabel('SARS-CoV-2-ACE2 ((units)/(sec))')
legend('norm ASL gluc', 'hyper ASL gluc')
axis square
title('high viral load')
hold off



figure(12)  %

plot(tsecs,ACE2V_n,'color', [0.4 0.4 0.4] ,'linewidth',3)
hold on
plot(tsecs,ACE2V_h,'--k' ,'linewidth',3)

plot(tsecs,ACE2V_nb,'color', [0.4 0.4 0.4] ,'linewidth',3)
plot(tsecs,ACE2V_hb,'--k','linewidth',3)

plot(tsecs,ACE2V_nc,'color', [0.4 0.4 0.4] ,'linewidth',3)
plot(tsecs,ACE2V_hc,'--k','linewidth',3)

hold off
legend('norm ASL gluc, low viral', 'hyper ASL gluc, low viral', 'norm ASL gluc, intermediate viral','hyper ASL gluc, intermediate viral', 'norm ASL gluc, high viral','hyper ASL gluc, high viral')
xlim([0 1e6])
xlabel('secs')
ylabel('SARS-CoV-2 endocytosis ((units)/(sec))')





figure (13)

m= gluc_asl_r;
n= lac_asl ;
o= pH_asl;
yyaxis left
plot(m, n,'color',[.03 .6 .8],'linewidth',3)
ylabel('ASL lactate (mM)')
yyaxis right
set(gca,'ycolor','k') 
plot(m, o,'k','linewidth',3)
ylabel('ASL pH')
legend('ASL lactate', 'ASL pH', 'Location', 'north')
xlabel('ASL glucose (mM)') 
xlim([0.1 2])
ylim([6 7])
set(gca,'ycolor','k') 
%title('')


figure (14)

plot(gluc_asl_r, lac_asl,'color',[.03 .6 .8], 'linewidth',3)
hold on
plot(gluc_asl_r, lac_epi,'k','linewidth',3)
xlabel('ASL glucose (mM)') 
ylabel('lactate (mM)')
xlim([0.1 2])
ylim([0 12])
legend('ASL', 'Epithelium', 'Location', 'north')




figure(15)  

plot(gluc_asl_r,vireprate_epi,'color', [0.7 0.7 0.7],'linewidth',3)
hold on
plot(gluc_asl_r,vireprate_epib,'color',[0.4 0.4 0.4],'linewidth',3)
plot(gluc_asl_r,vireprate_epic,'color', [0 0 0] ,'linewidth',3)
hold off
legend('low viral load','intermediate viral load', 'high viral load', 'Location', 'northwest')
xlim([0.15 2])
ylim([0 0.4])
xlabel('ASL glucose (mM)')
ylabel('virus reproduction rate (number/sec)')
%title('')


figure(16)

plot(gluc_asl_r,vinumber_epi,'color', [0.7 0.7 0.7],'linewidth',3)
hold on
plot(gluc_asl_r,vinumber_epib,'color',[0.4 0.4 0.4],'linewidth',3)
plot(gluc_asl_r,vinumber_epic,'color',[0 0 0],'linewidth',3)
hold off
legend('low viral ASL load', 'intermediate viral ASL load', 'high viral ASL load', 'Location', 'northwest')
xlim([0.1 2])
xlabel('ASL glucose (mM)')
ylabel('epithelial virus number')
%title('')


figure(17)  
plot(gluc_asl_r,vinumber_log,'color',[0.7 0.7 0.7], 'linewidth',3)
hold on
plot(gluc_asl_r,vinumber_logb,'color',[0.4 0.4 0.4] ,'linewidth',3)
plot(gluc_asl_r,vinumber_logc,'color', [0 0 0],'linewidth',3)
legend('low viral ASL load', 'intermediate viral ASL load', 'high viral ASL load', 'Location', 'northwest')
xlim([0.1 2])
xlabel('ASL glucose (mM)')
ylabel('log10(epithelial virus number)')
%title('') 





figure(32)


value1=interp1(gluc_asl_r,vireprate_epi,0.4);
value2=interp1(gluc_asl_r,vireprate_epib,0.4);
value3=interp1(gluc_asl_r,vireprate_epic,0.4);
value4=interp1(gluc_asl_r,vireprate_epi,1.2);
value5=interp1(gluc_asl_r,vireprate_epib,1.2);
value6=interp1(gluc_asl_r,vireprate_epic,1.2);

x = categorical({'normoglycemic' , 'hyperglycemic'});
x = reordercats(x,{'normoglycemic', 'hyperglycemic'});
vals = [value1 value2 value3; value4 value5 value6];
bar(x ,vals, 0.5)
b = bar(x,vals,0.5);
b(1).FaceColor=[0.7 0.7 0.7];
b(2).FaceColor=[0.4 0.4 0.4];
b(3).FaceColor=[0 0 0];
legend('low viral load', 'intermediate viral load', 'high viral load','location', 'northwest')
xlabel('ASL glucose')
ylabel('virus replication rate')


figure(33)


value1=interp1(gluc_asl_r,vinumber_epi,0.4);
value2=interp1(gluc_asl_r,vinumber_epib,0.4);
value3=interp1(gluc_asl_r,vinumber_epic,0.4);
value4=interp1(gluc_asl_r,vinumber_epi,1.2);
value5=interp1(gluc_asl_r,vinumber_epib,1.2);
value6=interp1(gluc_asl_r,vinumber_epic,1.2);



x = categorical({'normoglycemic' , 'hyperglycemic'});
x = reordercats(x,{'normoglycemic', 'hyperglycemic'});
vals = [value1 value2 value3; value4 value5 value6];
bar(x ,vals, 0.5)
b = bar(x,vals,0.5);
b(1).FaceColor=[0.7 0.7 0.7];
b(2).FaceColor=[0.4 0.4 0.4];
b(3).FaceColor=[0 0 0];
legend('low viral load', 'intermediate viral load', 'high viral load','location', 'northwest')
xlabel('ASL glucose')
ylabel('epithelial virus number')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('end ... but wait for the figures')
disp('fin ... mais attendez les graphiques')

