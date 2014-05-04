%Chemical Engineering 184B Design Project by Matt Carroll and Kyle Wu
%desgin level 3 mwc last edited 4/29/2014 5:23pm

clc; clear;
tic
fprintf('starting . . .\n')

%gloabals
%going into reactor 1
global ca0 %mol/L

%reactor 1
global tau1
global ccat1 %range from .26-0.78 g/L
global ph21 %range from 9.87atm to 39.48atm
%leaving reactor 1
global ca1 %mol/L IBAP
global cu1 %mol/L 
global cw1 %mol/L 
global cd1 %mol/L 

%reactor2
global tau2 %s^-1
global ccat2 %range from .26-0.78 g/L
global ph22 %range from 9.87atm to 39.48atm
%leaving reactor 2
global ca2 %mol/L IBAP
global cu2 %mol/L 
global cw2 %mol/L 
global cd2 %mol/L 

%reactor3
global tau3 %s^-1
global ccat3 %range from .26-0.78 g/L
global ph23 %range from 9.87atm to 39.48atm
%leaving reactor 3
global ca3 %mol/L IBAP
global cu3 %mol/L 
global cw3 %mol/L 
global cd3 %mol/L 

%reactor4
global tau4; %s^-1
global ccat4; %range from .26-0.78 g/L
global ph24; %range from 9.87atm to 39.48atm
%leaving reactor 3
global ca4 %mol/L IBAP
global cu4 %mol/L 
global cw4 %mol/L 
global cd4 %mol/L 

%constants

%henry's law constnat for hydrogen and methane
Hh2=.00036; %atm^-1
Hch4=.0021; %amt^-1

%inital concentration
ca0=5.4; %mol/L (PURE IBAP)

%kinetic constants
k1=1.14;
k2=.095;
ka=76.4;
kh=141;
kw=529;

%open files for logging and give assumptions and inital log
file_path=[pwd '_log.csv']; fprintf('log stored:\n%s\n\n',file_path);
fid1=fopen(file_path,'w'); fprintf(fid1, 'ChE184A MC/CW/KW\n\n\nall flowrates in moles per second based on a 8400hours/year plant runntime\nreactor volume is in liters\ninitail concentration of reactant a is %f\nprofit is in $/year\ntau (space-time) is in seconds\ninlet density is in g/L\nvolumetric flowrate q is in L/s\n\n\n\n\n', ca0);
fprintf(fid1,'NPV_proj($),NPV_percent,TI($),TCI($),POT(years),X,S,ccat1(g/L),ph21(atm),tau1(s^-1),ca1(mol/L),cu1(mol/L),cw1(mol/L),cd1(mol/L),ra1(mol/Ls),rb1(mol/Ls),x1,s1,pd1(mol/s),fa1(mol/s),q01(L/s),rho01(g/L),rho1(g/L),q1(L/s),liq_mol_flow_rate_out1(mol/s),V1(L),h2_mol_frac_liq1,ph1(atm),h2_reacted1 (mol/s),fh1(mol/s),fch41,pch41,pressure_ch41,p_total1,operating_profit_per_year($),deltaHa,deltaHb,heat_of_rxn1,ROIbt,profitbt,pa1,pu1,pw1,fa_fresh(mol/s),cost_gas_comp($),cost_heat_ex($),cost_cstr($),FC($),WC($),utilties_cost($),profit_flow_stream,cost_cstr,cost_heat_ex,cost_gas_comp($),sep_cost($),FC,diameter_flash (m),height_flash (m),height_flash (m),dist1_area (m),dist1_height (m),dist2_area (m),dist2_height (m),,volume_flash(L), dist1_r (reflux ratio), dist1_s (boilup ratio), dist1_N_real(stages), dist2_r(reflux ratio), dist2_s(boilup ratio), dist2_N_real (stages)\n');
%loop control iterator varibles and number of iterations
nph2=100;
nccat=100;
ntau=100;
iter=0;
k=1;


%reactor 1
for ph21=39.47%linspace(9.87,39.47,nph2)%39.47%linspace(9.87,39.47,nph2) %ph2 range
for ccat1=.78%linspace(.26,.78,nccat) %ccat range
for tau1=1145.4545%logspace(1,7,ntau);%linspace(1000, 1600,ntau)%%%657.9332
    %reactor 2
    for ph22=ph21
    for ccat2=ccat1
    for tau2=tau1
        %reactor 3
        for ph23=ph21
        for ccat3=ccat1
        for tau3=tau1
            %reactor 4
            for ph24=ph21
            for ccat4=ccat1%ccat range
            for tau4=tau1
                iter=iter+1;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%% REACTOR %%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %reactor 1
                concentrations1=fsolve(@CSTR_bal1,[1,1,1,1],optimset('Display','off'));
                ca1=concentrations1(1);%mol/L 
                cu1=concentrations1(2);%mol/L 
                cw1=concentrations1(3);%mol/L 
                cd1=concentrations1(4);%mol/L 

                %reactor 2
                concentrations2=fsolve(@CSTR_bal2,[1,1,1,1],optimset('Display','off'));
                ca2=concentrations2(1);%mol/L 
                cu2=concentrations2(2);%mol/L 
                cw2=concentrations2(3);%mol/L 
                cd2=concentrations2(4);%mol/L 

                %reactor 3
                concentrations3=fsolve(@CSTR_bal3,[1,1,1,1],optimset('Display','off'));
                ca3=concentrations3(1);%mol/L 
                cu3=concentrations3(2);%mol/L 
                cw3=concentrations3(3);%mol/L 
                cd3=concentrations3(4);%mol/L

                %reactor 4
                concentrations4=fsolve(@CSTR_bal4,[1,1,1,1],optimset('Display','off'));
                ca4=concentrations4(1);%mol/L 
                cu4=concentrations4(2);%mol/L 
                cw4=concentrations4(3);%mol/L 
                cd4=concentrations4(4);%mol/L

                %rate expressions
                ra1=(ccat1.*k1*ca1*ph21)./(1+ka*ca1+(kh*ph21).^(1./2)+kw*cw1).^(2);
                rb1=(ccat1*k2*cd1*ph21)./(1+ka*ca1+(kh*ph21).^(1./2)+kw*cw1).^(2);
                ra2=(ccat2.*k1*ca2.*ph22)./(1+ka*ca2+(kh*ph22).^(1./2)+kw*cw2).^(2);
                rb2=(ccat2*k2*cd2*ph22)./(1+ka*ca2+(kh*ph22).^(1./2)+kw*cw2).^(2);
                ra3=(ccat3.*k1*ca3.*ph23)./(1+ka*ca3+(kh*ph23).^(1./2)+kw*cw3).^(2);
                rb3=(ccat3*k2*cd3*ph23)./(1+ka*ca3+(kh*ph23).^(1./2)+kw*cw3).^(2);
                ra4=(ccat4*k1*ca4*ph24)./(1+ka*ca4+(kh*ph24).^(1./2)+kw*cw4).^(2);
                rb4=(ccat4*k2*cd4*ph24)./(1+ka*ca4+(kh*ph24).^(1./2)+kw*cw4).^(2);

                %conversion and selectivity
                x1=(ca0-ca1)./ca0;
                s1=cd1./(ca0-ca1);
                x2=(ca1-ca2)./ca1;
                s2=(cd2-cd1)./(ca1-ca2);
                x3=(ca2-ca3)./ca2;
                s3=(cd3-cd2)./(ca2-ca3);
                x4=(ca3-ca4)./ca3;
                s4=(cd4-cd3)./(ca3-ca4);

                X=(ca0-ca4)/ca0; % total converstion of the reactor system
                S=(cd4)./(ca0-(ca4));
                
                %outlet density in g/L
                rho1=ca1*178.25484+cd1*178.27+cu1*162.28+cw1*18.02;
                rho2=ca2*178.25484+cd2*178.27+cu2*162.28+cw2*18.02;
                rho3=ca3*178.25484+cd3*178.27+cu3*162.28+cw3*18.02;
                rho4=ca4*178.25484+cd4*178.27+cu4*162.28+cw4*18.02;

                %inlet density in g/L 
                rho01=952;%(assumeing pure IBAP)
                rho02=rho1;%out of reactor 1 =inlet of reactor 2
                rho03=rho2;%out of reactor 2 =inlet of reactor 3
                rho04=rho3;%out of reactor 3 =inlet of reactor 4

                %product flowrates in moles per second based on a 8400hours/year plant runntime
                pd4=1.006.*((5000000.*1000)./178.3)/(8400*60*60); %mol/s

                %liquid molar flowrate out reactor 3
                liq_mol_flow_rate_out4=pd4./(cd4./(ca4+cu4+cw4+cd4));%mol/s

                %flowrate out reactor 3
                q4=liq_mol_flow_rate_out4./(ca4+cu4+cw4+cd4); %L/s

                %total flowrates L/s
                q04=(q4*rho3)/rho04; %L/s
                q3=q04; %L/s
                q03=(q3*rho3)/rho03; %L/s
                q2=q03; %L/s
                q02=(q2*rho2)/rho02; %L/s
                q1=q02; %L/s
                q01=(q1*rho1)/rho01; %L/s

                %componunt flowrates
                pa1=ca1*q1; %mol/s
                pu1=cu1*q1; %mol/s
                pw1=cw1*q1; %mol/s
                pd1=cd1*q1; %mol/s
                pa2=ca2*q2; %mol/s
                pu2=cu2*q2; %mol/s
                pw2=cw2*q2; %mol/s
                pd2=cd2*q2; %mol/s
                pa3=ca3*q3; %mol/s
                pu3=cu3*q3; %mol/s
                pw3=cw3*q3; %mol/s
                pd3=cd3*q3; %mol/s
                pa4=ca4*q4; %mol/s
                pu4=cu4*q4; %mol/s
                pw4=cw4*q4; %mol/s
                pd4=cd4*q4; %mol/s

                fa_fresh=pd4./S;
                fa1=(.999.*.995).*pa4+fa_fresh; %mol/s
                fa2=pa1;
                fa3=pa2;
                fa4=pa3;

                %liquid molar flowrate out
                liq_mol_flow_rate_out1=q1*(ca1+cu1+cw1+cd1);%mol/s
                liq_mol_flow_rate_out2=q2*(ca2+cu2+cw2+cd2);%mol/s
                liq_mol_flow_rate_out3=q3*(ca3+cu3+cw3+cd3);%mol/s
                liq_mol_flow_rate_out4=q4*(ca4+cu4+cw4+cd4);%mol/s

                %reactor volume in liters
                V1=tau1*q01;
                V2=tau2*q02;
                V3=tau3*q03;
                V4=tau4*q04;

                %flow hydrogen out
                h2_mol_frac_liq1=ph21*Hh2; %mol fraction
                ph1=h2_mol_frac_liq1*liq_mol_flow_rate_out1; %mol/s
                h2_mol_frac_liq2=ph22*Hh2; %mol fraction
                ph2=h2_mol_frac_liq2*liq_mol_flow_rate_out2; %mol/s
                h2_mol_frac_liq3=ph23*Hh2; %mol fraction
                ph3=h2_mol_frac_liq3*liq_mol_flow_rate_out3; %mol/s
                h2_mol_frac_liq4=ph24*Hh2; %mol fraction
                ph4=h2_mol_frac_liq4*liq_mol_flow_rate_out4; %mol/s

                %hydrogen reacted                
                h2_reacted1=V1*ra1+V1*rb1; %mol/s  
                h2_reacted2=V2*ra2+V2*rb2; %mol/s  
                h2_reacted3=V3*ra3+V3*rb3; %mol/s  
                h2_reacted4=V4*ra4+V4*rb4; %mol/s  

                %flow hydrogen in
                fh1=h2_reacted1+ph1; % mol/s in through compresor
                fh2=h2_reacted2+ph2-ph1; % mol/s in through compresor
                fh3=h2_reacted3+ph3-ph2; % mol/s in through compresor
                fh4=h2_reacted4+ph4-ph3; % mol/s in through compresor

                %flow of ch4,in=ch4,out
                fch41=fh1*.01; %mol/s in through compresor
                pch41=fch41; %mol/s
                fch42=fh2*.01; %mol/s in through compresor
                pch42=fch42+pch41; %mol/s
                fch43=fh3*.01; %mol/s in through compresor
                pch43=fch43+pch42; %mol/s
                fch44=fh4*.01; %mol/s in through compresor
                pch44=fch44+pch43; %mol/s

                %pressure of ch4 in reactor
                ch4_mol_frac_liq2=pch42./liq_mol_flow_rate_out2;
                pressure_ch42=ch4_mol_frac_liq2/Hch4;
                ch4_mol_frac_liq1=pch41./liq_mol_flow_rate_out1;
                pressure_ch41=ch4_mol_frac_liq1/Hch4;
                ch4_mol_frac_liq3=pch43./liq_mol_flow_rate_out3;
                pressure_ch43=ch4_mol_frac_liq3/Hch4;
                ch4_mol_frac_liq4=pch44./liq_mol_flow_rate_out4;
                pressure_ch44=ch4_mol_frac_liq4/Hch4;

                %total pressure
                p_total1=pressure_ch41+ph21;%atm
                p_total2=pressure_ch42+ph22; %atm
                p_total3=pressure_ch43+ph23; %atm
                p_total4=pressure_ch44+ph24; %atm

                %heat of reaction
                deltaHa=-3.7; %kJ/molIBAP
                deltaHb=7.9; %kJ/molIBAP
                heat_of_rxn1=V1*ra1*deltaHa+V1*rb1*deltaHb; %kW
                heat_of_rxn2=V2*ra2*deltaHa+V2*rb2*deltaHb; %kW
                heat_of_rxn3=V3*ra3*deltaHa+V3*rb3*deltaHb; %kW
                heat_of_rxn4=V4*ra4*deltaHa+V4*rb4*deltaHb; %kW

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%% SEPERATOR %%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %FLASH SEPERATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cad0=ca4;  %mol/L 
                cud0=cu4;  %mol/L     
                cwd0=cw4;  %mol/L 
                cdd0=cd4;  %mol/L 
                
                %mole fractions
                flash_gx1=ph4;%h2
                flash_gx2=ch4_mol_frac_liq4; %ch4
                flash_gx3=cwd0./(cad0+cud0+cwd0+cdd0); %H2O
                flash_gx4=cud0./(cad0+cud0+cwd0+cdd0); %IBEB
                flash_gx5=cdd0./(cad0+cud0+cwd0+cdd0); %IBPE
                flash_gx6=cad0./(cad0+cud0+cwd0+cdd0); %IBAP
                flash_za=flash_gx1./(flash_gx1+flash_gx2+flash_gx3+flash_gx4+flash_gx5+flash_gx6);
                flash_zb=flash_gx2./(flash_gx1+flash_gx2+flash_gx3+flash_gx4+flash_gx5+flash_gx6);
                flash_zc=flash_gx3./(flash_gx1+flash_gx2+flash_gx3+flash_gx4+flash_gx5+flash_gx6);
                flash_zd=flash_gx4./(flash_gx1+flash_gx2+flash_gx3+flash_gx4+flash_gx5+flash_gx6);
                flash_ze=flash_gx5./(flash_gx1+flash_gx2+flash_gx3+flash_gx4+flash_gx5+flash_gx6);
                flash_zf=flash_gx6./(flash_gx1+flash_gx2+flash_gx3+flash_gx4+flash_gx5+flash_gx6);
                
                flash_z=[flash_za flash_zb flash_zc flash_zd flash_ze flash_zf];
                flash_F=liq_mol_flow_rate_out4.*(1+flash_gx2).*(1+ph4);
                
                %vapor pressures @373K
                flash_avap=22973021876598800000000000000000; %estimate
                flash_bvap=22745177; %estimate
                flash_cvap=33830944.0569.*exp(-4737.3334.*(1/373));
                flash_dvap=10324955.9950.*exp(-5625.9011.*(1/373));
                flash_evap=14831374.3039.*exp(-6502.6626.*(1/373));
                flash_fvap=13764206.8239.*exp(-6661.4395.*(1/373));
                flash_Pvap=[flash_avap flash_bvap flash_cvap flash_dvap flash_evap flash_fvap];
              
                %calculate K's
                flash_y=flash_z.*flash_Pvap./101.325;
                flash_K=flash_y./flash_z;
                
                %calculate phi
                solve_phi=@(var)(sum((flash_z.*(flash_K-1))./(1+var.*(flash_K-1)))); %solve for phi
                phi=fzero(solve_phi,[0 1]);
                
                %size flash drum
                flash_x=flash_z./(1+phi.*(flash_K-1));
                flash_y=flash_K.*flash_x;
                V_flash=phi*flash_F;
                L_flash=(1-phi).*flash_F;
                rho_l=flash_gx3*1+flash_gx4*.861+flash_gx5*.955+flash_gx6*.952;
                rho_v=flash_gx1*.0000899+flash_gx2*.000667151;
                k_drum=.35;
                u_perm=k_drum.*sqrt((rho_l-rho_v)/rho_v);
                flash_a_c=(V_flash*((flash_x(1).*1.01+flash_x(2).*16.05)./(flash_x(1)+flash_x(2))))./(u_perm*rho_v)*.0001; %m^2
                diameter_flash=sqrt((4*flash_a_c)/pi); %m
                flash_diameter_ft=diameter_flash*3.281;
                height_flash=3*(diameter_flash); %m
                flash_height_ft=height_flash.*3.281;
                volume_flash=3*(diameter_flash)*pi*(diameter_flash/2).^2*1000; %L
                
                %DIST COLUMN ONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %flow to distillation system
                dist1_F=L_flash; %mol/s
                
                %mole fractions
                dist1_xa=flash_x(5); %IBPE
                dist1_xb=flash_x(6); %IBAP
                dist1_xc=flash_x(4); %IBEB
                dist1_xd=flash_x(3); %H2O
                dist1_z_feed=[dist1_xa,dist1_xb,dist1_xc,dist1_xd];
                
                %assume liquid compotition i.e. q=1
                dist1_x=dist1_z_feed;
                dist1_q1=1; %assume liquid only feed (see pg82)
                
                %parameters of exp relation btw vp and temp
                %vapor pressure functions valid for 373K to 573K
                %ibpe ibap ibeb h20
                dist1_a=[14831374.3039,13764206.8239,10324955.9950,33830944.0569];
                dist1_b=[6661.4395,6502.6626,5625.9011,4737.3334];
                dist1_P=101.3;
                
                %calculate buble point
                dist1_Pvap=@(var)(dist1_a.*exp(-dist1_b.*(1/var))); %antoine
                bp=@(var)(sum(dist1_x.*dist1_Pvap(var)./dist1_P)-1); %bubble pt calc
                dist1_T=fzero(bp,300);
                dist1_y=dist1_x.*dist1_Pvap(dist1_T)./dist1_P;
   
                % assume 100% recovery of LK in distillate
                % assume 100% recovery of HK in bottoms
                dist1_K=dist1_y./dist1_x;
                dist1_alpha=dist1_K./min(dist1_K);
                dist1_alpha_a=dist1_alpha(1);
                dist1_alpha_b=dist1_alpha(2);
                dist1_alpha_c=dist1_alpha(3);
                dist1_alpha_d=dist1_alpha(4);

                %caluclate D and B assuming 100% sep
                dist1_D=(dist1_xc+dist1_xd).*dist1_F; %mol/s
                dist1_B=dist1_F-dist1_D; %mol/s

                %calculate R min
                %table 4.1 for AB/CD split LK=B=IBEB HK=C=IBAP
                dist1_r_min=(((dist1_alpha_b*dist1_xd)/(dist1_alpha_d-dist1_alpha_b))+((dist1_alpha_b.*(dist1_xc+dist1_xb))/(dist1_alpha_c-dist1_alpha_b)))./((dist1_xd+dist1_xc).*(1+dist1_xd.*(dist1_xb+dist1_xa)))+(dist1_xd.*((dist1_xd/(dist1_alpha_d-1))+(dist1_xc/(dist1_alpha_c-1))))/((dist1_xd+dist1_xc).^2);
                
                dist1_r=1.5.*dist1_r_min;
                dist1_s=(dist1_D./dist1_B).*(dist1_r+dist1_q1)-(1-dist1_q1);
                dist1_N_min=(log((.999./.001).*(.999./.001))/log(dist1_alpha_c/dist1_alpha_b));
                N_calc=dist1_N_min;
                error=1;
                while error > 1e-3;
                    N1old=N_calc;
                    N_calc=0.75.*(1-((dist1_r-dist1_r_min)/(dist1_r+1)).^0.5688).*(N_calc+1)+dist1_N_min;
                    error=abs(N_calc-N1old);
                end
                dist1_N_thoretical=N_calc;%fzero(dist1_FUG, [1 1000]);
                dist1_N=ceil(dist1_N_thoretical./((1-.24).*exp(-sqrt(((dist1_alpha_b/dist1_alpha_c).*.3)/1))+.24));
                
                dist1_v_B=dist1_s.*dist1_B; %mol/s
                dist1_v_T=(dist1_r+1).*dist1_D; %mol/s
                
                
                %vapor mole fractions
                dist1_y_h20=dist1_xd./(dist1_xd+.999.*dist1_xc+.001.*dist1_xb);
                dist1_y_ibeb=(.999.*dist1_xc)./(dist1_xd+.999.*dist1_xc+.001.*dist1_xb);
                dist1_y_ibap=(.001.*dist1_xb)./(dist1_xd+.999.*dist1_xc+.001.*dist1_xb);
                
                %liquid mole fractions
                dist1_x_ibpe=dist1_xa./(dist1_xa+.999.*dist1_xb+.001.*dist1_xc);
                dist1_x_ibap=(.999.*dist1_xb)./(dist1_xa+.999.*dist1_xb+.001.*dist1_xc);
                dist1_x_ibeb=(.001.*dist1_xc)./(dist1_xa+.999.*dist1_xb+.001.*dist1_xc);
                
                %calculate heat loads
                dist1_heat_cap_d=(dist1_y_h20.*88.7+dist1_y_ibeb.*428+dist1_y_ibap.*416.2).*(dist1_T-273); %J/mol
                dist1_heat_cap_b=(dist1_x_ibpe.*424.2+dist1_x_ibeb.*428+dist1_x_ibap.*416.2).*(dist1_T-273); %J/mol
                dist1_qc=dist1_heat_cap_d*dist1_v_T./1000; %kJ/s
                dist1_qR=dist1_heat_cap_b*dist1_v_T./1000; %kJ/s
                
                %sizing column 1
                dist1_mv=dist1_y_h20.*18.02+dist1_y_ibeb.*162.28+dist1_y_ibap.*176.25;
                dist1_ml=dist1_x_ibpe.*178.27+dist1_x_ibeb.*162.28+dist1_x_ibap.*176.25;
                dist1_rho_l=dist1_y_h20.*1+dist1_y_ibeb.*.861+dist1_y_ibap.*.952;
                dist1_rho_v=dist1_x_ibpe.*.955+dist1_x_ibeb.*.861+dist1_x_ibap.*.952;
                dist1_A=(dist1_mv.*dist1_v_T.*3600)./(sqrt(dist1_rho_l*dist1_rho_v).*0.6.*.8.*252*1000000); %m^2
                
                dist1_F_flood=sqrt(dist1_rho_l).*1.25*252.*2278;
                dist1_f=((dist1_rho_v/dist1_rho_l).*1000).^.5*(dist1_mv/dist1_ml).^1.5;%*(dist1_B/dist1_D)
                dist1_c=252/(1+2*dist1_f);
                dist1_F_flood_new=0.8.*sqrt((dist1_rho_v-dist1_rho_l).*1000)*dist1_c;
                dist1_area1=dist1_mv/(0.6*dist1_F_flood_new*sqrt(dist1_rho_v.*1000));
                dist1_area2=dist1_ml/(0.6*dist1_F_flood_new*sqrt(dist1_rho_v.*1000));
                dist1_area=(dist1_area1+dist1_area2)./2; %m^2
                dist1_diameter_ft=(sqrt(dist1_area/pi).*2).*3.281;
                
                tray_height=.31; %m
                dist1_height=tray_height.*(dist1_N+3); %m
                dist1_height_ft=dist1_height.*3.281;
                
                %calculate delta T of column
                %ibpe ibap ibeb h20
                dist1_a=[14831374.3039,13764206.8239,10324955.9950];
                dist1_b=[6661.4395,6502.6626,5625.9011];
                dist1_P=101.3;
                % Tbottoms
                xB=[dist1_x_ibpe, dist1_x_ibap, dist1_x_ibeb];
                dist1_Pvap=@(var)(dist1_a.*exp(-dist1_b.*(1/var))); %antoine
                bp=@(var)(sum(xB.*dist1_Pvap(var)./dist1_P)-1); %bubble pt calc
                dist1_T_bottom=fzero(bp,300);
                %ibpe ibap ibeb h20
                dist1_a=[13764206.8239,10324955.9950,33830944.0569];
                dist1_b=[6502.6626,5625.9011,4737.3334];
                dist1_P=101.3;
                % Ttop
                xD=[dist1_y_ibap, dist1_y_ibeb, dist1_y_h20];
                dist1_Pvap=@(var)(dist1_a.*exp(-dist1_b.*(1/var))); %antoine
                bp=@(var)(sum(xD.*dist1_Pvap(var)./dist1_P)-1); %bubble pt calc
                dist1_T_top=fzero(bp,300);
                dist1_delta_T=dist1_T_bottom-dist1_T_top;
                               
                %DIST COLUMN TWO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %calculate flow from bottom and compisition
                dist2_F=dist1_B;
                dist2_xa=(.001*dist1_xc.*dist1_F)./(dist1_xa.*dist1_F+.999*dist1_xb.*dist1_F+.001.*dist1_xc.*dist1_F); %IBEB
                dist2_xb=(.999*dist1_xb.*dist1_F)./(dist1_xa.*dist1_F+.999*dist1_xb.*dist1_F+.001.*dist1_xc.*dist1_F); %IBAP
                dist2_xc=(dist1_xa.*dist1_F)./(dist1_xa.*dist1_F+.999*dist1_xb.*dist1_F+.001.*dist1_xc.*dist1_F); %IBPE
                dist2_z_feed=[dist2_xa,dist2_xb,dist2_xc];
                
                %assume liquid compotition i.e. q=1
                dist2_x=dist2_z_feed;
                
                %parameters of exp relation btw vp and temp
                %vapor pressure functions valid for 373K to 573K
                dist2_a=[10324955.9950,13764206.8239,14831374.3039];
                dist2_b=[5625.9011,6502.6626,6661.4395];
                dist2_P=101.325;
                
                %calculate buble point
                dist2_Pvap=@(var)(dist2_a.*exp(-dist2_b.*(1/var))); %antoine
                bp=@(var)(sum(dist2_x.*dist2_Pvap(var)./dist2_P)-1); %bubble pt calc
                dist2_T=fzero(bp,300);
   
                % assume 100% recovery of LK in distillate
                % assume 100% recovery of HK in bottoms
                dist2_y=dist2_x.*dist2_Pvap(dist2_T)./dist2_P;
                flash_K=dist2_y./dist2_x;
                dist2_alpha=flash_K/min(flash_K);
                dist2_alpha_a=dist2_alpha(1);
                dist2_alpha_b=dist2_alpha(2);
                dist2_alpha_c=dist2_alpha(3);
                dist2_q=1; %assume liquid only feed (see pg82)
                
                %caluclate D and B
                dist2_D=(dist2_xb+dist2_xa).*dist2_F;
                dist2_B=dist2_F-dist2_D;

                %reflux ration
                dist2_r_min=(((dist2_xb+dist2_xc)./(dist2_alpha_b-1))+dist2_xa./(dist2_alpha_a-1))./((dist2_xa+dist2_xb).*(1+dist2_xa.*dist2_xc));
                
                dist2_r=1.5.*dist2_r_min;
                dist2_s=(dist2_D./dist2_B).*(dist2_r+dist2_q)-(1-dist2_q);
                dist2_N_min=(log((.995./.005).*(.995./.005))/log(dist2_alpha_b/dist2_alpha_c));
                N_calc=dist2_N_min;
                error=1;
                while error > 1e-3;
                    N2old=N_calc;
                    N_calc=0.75.*(1-((dist2_r-dist2_r_min)/(dist2_r+1)).^0.5688).*(N_calc+1)+dist2_N_min;
                    error=abs(N_calc-N2old);
                end
                dist2_N_thoretical=N_calc;%fzero(dist1_FUG, [1 1000]);
                dist2_N_real=ceil(dist2_N_thoretical./((1-.24).*exp(-sqrt(((dist2_alpha_b/dist2_alpha_c).*.3)/1))+.24));
                
                dist2_v_B=dist2_s.*dist2_B; %mol/s
                dist2_v_T=(dist2_r+1).*dist2_D; %mol/s

                check1=dist2_v_B-dist2_v_T;
                check2=(dist2_q-1).*dist2_F;
                
                %vapor mole fractions
                dist2_y_ibpe=(.001.*dist2_xc)./(dist2_xa+.999.*dist2_xb+.001.*dist2_xc);
                dist2_y_ibap=(.999.*dist2_xb)./(dist2_xa+.999.*dist2_xb+.001.*dist2_xc);
                dist2_y_ibeb=(dist2_xa)./(dist2_xa+.999.*dist2_xb+.001.*dist2_xc);

                %liquid mole fractions
                dist2_x_ibpe=dist2_xc./(dist2_xc+.001.*dist2_xb);
                dist2_x_ibap=(.001.*dist2_xb)./(dist2_xc+.001.*dist2_xb);
                
                %calculate heat loads
                dist2_heat_cap_d=(dist2_y_ibpe.*424.2+dist2_y_ibeb.*428+dist2_y_ibap.*416.2).*(dist2_T-273); %J/mol
                dist2_heat_cap_b=(dist2_x_ibpe.*424.2+dist2_x_ibap.*416.2).*(dist2_T-273); %J/mol
                dist2_qc=dist2_heat_cap_d*dist2_v_T./1000; %kJ/s
                dist2_qR=dist2_heat_cap_b*dist2_v_T./1000; %kJ/s
                
                %sizing column 2
                dist2_mv=dist2_y_ibpe.*178.27+dist2_y_ibeb.*162.28+dist2_y_ibap.*176.25;
                dist2_ml=dist2_x_ibpe.*178.27+dist2_x_ibap.*176.25;
                dist2_rho_l=dist2_x_ibpe.*.955+dist2_x_ibap.*.952+dist2_y_ibeb.*.861;
                dist2_rho_v=dist2_y_ibpe.*.955+dist2_y_ibap.*.952;
                dist2_A=(dist2_mv.*dist2_v_T.*3600)./(sqrt(dist2_rho_l*dist2_rho_v).*0.6.*.8.*252*1000000); %m^2
                
                dist2_F_flood=sqrt(dist2_rho_l).*1.25*252.*2278;
                dist2_f=((dist2_rho_v/dist2_rho_l).*1000).^.5*(dist2_mv/dist2_ml).^1.5;%*(dist2_B/dist2_D)   
                dist2_c=252/(1+2*dist2_f);
                dist2_F_flood_new=0.8.*sqrt((dist2_rho_l-dist2_rho_v).*1000)*dist2_c;
                dist2_area1=dist2_mv/(0.6*dist2_F_flood_new*sqrt(dist2_rho_v.*1000));
                dist2_area2=dist2_ml/(0.6*dist2_F_flood_new*sqrt(dist2_rho_v.*1000));
                dist2_area=(dist2_area1+dist2_area2)./2;
                dist2_diameter_ft=(sqrt(dist2_area/pi).*2).*3.281;
                                
                tray_height=.31; %m
                dist2_height=tray_height.*(dist2_N_real+3); %m
                dist2_height_ft=dist2_height.*3.281;
                
                %calculate delata T of column
                %ibpe ibap
                dist2_a=[14831374.3039,13764206.8239];
                dist2_b=[6661.4395,6502.6626];
                dist2_P=101.3;
                % Tbottoms
                xB=[dist2_x_ibpe, dist2_x_ibap];
                dist2_Pvap=@(var)(dist2_a.*exp(-dist2_b.*(1/var))); %antoine
                bp=@(var)(sum(xB.*dist2_Pvap(var)./dist2_P)-1); %bubble pt calc
                dist2_T_bottom=fzero(bp,300);
                %ibpe ibap ibeb
                dist2_a=[14831374.3039,13764206.8239,10324955.9950];
                dist2_b=[6661.4395,6502.6626,5625.9011];
                % Ttop
                xD=[dist2_y_ibpe, dist2_y_ibap, dist2_y_ibeb];
                dist2_Pvap=@(var)(dist2_a.*exp(-dist2_b.*(1/var))); %antoine
                bp=@(var)(sum(xD.*dist2_Pvap(var)./dist2_P)-1); %bubble pt calc
                dist2_T_top=fzero(bp,300);
                dist2_delta_T=dist2_T_bottom-dist2_T_top;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%% ECONOMICS %%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %fixed capital cost: marshall swift index current value
                ms=1550; 

                %all cost are based intastlatioin costs

                %Cstr costs (4 reactors)
                %reactor1
                fcstr_pressure=0.00130*(p_total1*14.6959488)+0.84;
                fcstr_material=1; %carbon steel
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V1/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr1=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $
                %reactor2
                fcstr_pressure=0.00130*(p_total2*14.6959488)+0.84;
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V2/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr2=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $
                %reactor3
                fcstr_pressure=0.00130*(p_total3*14.6959488)+0.84;
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V3/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr3=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $
                %reactor4
                fcstr_pressure=0.00130*(p_total4*14.6959488)+0.84;
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V4/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr4=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $
                %total reactor cost
                cost_cstr=cost_cstr1+cost_cstr2+cost_cstr3+cost_cstr4; % $

                %Heat Exchangers Reactors(4, one for each reactor)
                %heat xchgr 1
                V_feet=1.13.*V1*0.0353147;
                area=V_feet/(pi*0.25*diameter^2); %surface area of reactor
                fheat_material=1.3; %brass
                fheat_pressure=0; %factor for steam pressure assuming liq h2o <150psi
                fheat_design=0.85; %factor
                cost_heat_ex1=(ms/280).*(101.3*((fheat_pressure+fheat_design)*fheat_material+2.29)*area^0.65); % $
                %heat xchgr 2           
                V_feet=1.13.*V2*0.0353147;
                area=V_feet/(pi*0.25*diameter^2); %surface area of reactor
                fheat_material=1.3; %brass
                fheat_pressure=0; %factor for steam pressure assuming liq h2o <150psi
                fheat_design=0.85; %factor
                cost_heat_ex2=(ms/280).*(101.3*((fheat_pressure+fheat_design)*fheat_material+2.29)*area^0.65); % $
                %heat xchgr 3
                V_feet=1.13.*V3*0.0353147;
                area=V_feet/(pi*0.25*diameter^2); %surface area of reactor
                fheat_material=1.3; %brass
                fheat_pressure=0; %factor for steam pressure assuming liq h2o <150psi
                fheat_design=0.85; %factor
                cost_heat_ex3=(ms/280).*(101.3*((fheat_pressure+fheat_design)*fheat_material+2.29)*area^0.65); % $
                %heat xchgr 4
                V_feet=1.13.*V4*0.0353147;
                area=V_feet/(pi*0.25*diameter^2); %surface area of reactor
                fheat_material=1.3; %brass
                fheat_pressure=0; %factor for steam pressure assuming liq h2o <150psi
                fheat_design=0.85; %factor
                cost_heat_ex4=(ms/280).*(101.3*((fheat_pressure+fheat_design)*fheat_material+2.29)*area^0.65); % $
                %total of heat exchanger costs
                cost_heat_ex=cost_heat_ex1+cost_heat_ex2+cost_heat_ex3+cost_heat_ex4; % $

                %Gas Compressor
                %need to figure out how to relate pressure to bhp of gas compressor
                gamma_h2=1.405; %unitless ratio
                pressure_h2_in_atm=5;%atm
                pressure_h2_in=14.6959488.*pressure_h2_in_atm; %psia
                %V=nrt/P
                volumetric_flow_rate_h2_in=((fh1+fh2+fh3+fh4)*.08205946*298)./pressure_h2_in_atm; %L/s
                volumetric_flow_rate_h2_in_gal_per_min=volumetric_flow_rate_h2_in*15.8503231; %gal/min
                hp=(.0000303/gamma_h2).*pressure_h2_in.*1000*volumetric_flow_rate_h2_in_gal_per_min.*((ph24/pressure_h2_in_atm).^(gamma_h2)-1); %guthrie's correlation for compressor
                bhp=hp./.8; %a new wild guess
                fgas=1; % factor
                cost_gas_comp=0;%(ms/280)*517.5*(fgas+2.11)*bhp^0.82; % $

                %separating cost (1000000/X=$2,336,500)
                %flash
                flash_fm=1; %carbon steel
                flash_fp=1.6; %40bar=520psi
                flash_cost_column=(ms./280).*(101.9*flash_diameter_ft.^1.066.*flash_height_ft.^0.82.*(2.18.*flash_fm.*flash_fp));
                %columns
                %column1
                dist1_fm=1;
                dist1_cost_column=(ms./280).*(101.9*dist1_diameter_ft.^1.066.*dist1_height_ft.^0.82.*(2.18.*dist1_fm.*1));
                %trays and tower internals
                fs=2.2; %12in tray spacing
                ft=0; %normal tray type
                fm=0; %material=carbon steel
                dist1_fc=fs+ft+fm;
                dist1_cost_trays=(ms./280).*(4.7.*dist1_diameter_ft.^1.55.*dist1_height_ft.*dist1_fc);
                %column2
                dist2_fm=1;
                dist2_cost_column=(ms./280).*(101.9*dist2_diameter_ft.^1.066.*dist2_height_ft.^0.82.*(2.18.*dist2_fm.*1));
                %trays and tower internals
                fs=2.2; %12in tray spacing
                ft=0; %normal tray type
                fm=0; %material=carbon steel
                dist2_fc=fs+ft+fm;
                dist2_cost_trays=(ms./280).*(4.7.*dist2_diameter_ft.^1.55.*dist2_height_ft.*dist2_fc);
                %heat xchrs
                u=720;%w/m^2K condesnsing hydrocarbon cooled by boiling water
                logmean=120;
                dist1_condenser_heatxchr_area=(dist1_qc.*1000)./(u*logmean);
                dist1_condenser_heatxchr_area_ft=dist1_condenser_heatxchr_area.*10.7639;
                u=550;%w/m^2K boiling organic heated by low viscosity organic liquid
                dist1_reboiler_heatxchr_area=(dist1_qR.*1000)./(u*logmean);
                dist1_reboiler_heatxchr_area_ft=dist1_reboiler_heatxchr_area.*10.7639;
                u=720;%w/m^2K condesnsing hydrocarbon cooled by boiling water
                dist2_condenser_heatxchr_area=(dist2_qc.*1000)./(u*logmean);
                dist2_condenser_heatxchr_area_ft=dist2_condenser_heatxchr_area.*10.7639;
                u=550;%w/m^2K boiling organic heated by low viscosity organic liquid
                dist2_reboiler_heatxchr_area=(dist2_qR.*1000)./(u*logmean); %m^2
                dist2_reboiler_heatxchr_area_ft=dist2_reboiler_heatxchr_area.*10.7639;
                heat_xchr_cost=(ms./280).*((101.3.*dist1_condenser_heatxchr_area_ft.^0.65.*(2.29+1))+(101.3.*dist1_reboiler_heatxchr_area_ft.^0.65.*(2.29+1))+(101.3.*dist2_condenser_heatxchr_area_ft.^0.65.*(2.29+1))+(101.3.*dist2_reboiler_heatxchr_area_ft.^0.65.*(2.29+1)));
                sep_cost=dist1_cost_column+dist1_cost_trays+dist2_cost_column+dist2_cost_trays+heat_xchr_cost;
                
                
                %%%%%%%%%%%%%%%%%%%RECURNING UTLITIES COST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Gas compressor running cost
                bhp2=0;%hp./.6; % Guiethei
                cost_gas_comp_anuual=0;%bhp2.*8400*0.7456*.00588; % $/yr

                %reactor cooling cost
                cooling_annual=-(heat_of_rxn1+heat_of_rxn2+heat_of_rxn3+heat_of_rxn4).*8400*.00588; %$/yr
                
                %annual profit from flows
                profit_flow_stream=(pd4*.6-fa_fresh*.3-(fh1+fh2+fh3+fh4).*.00112)*(3600*8400); %$/yr

                %annual seperator cost (100000/X=$233,680)
                heating_kW=dist1_qc+dist1_qR+dist2_qc+dist2_qR;%kW
                heating_cost=heating_kW*8400.*.06; %$
                seperation_annual=heating_cost; %$/yr
                
                %annual catalyst cost
                catalyst_annual=ccat1*V1+ccat2*V2+ccat3*V3+ccat4*V4;%$/yr
                
                %TOTAL utitlities cost ($/year)
                utilties_cost=seperation_annual+catalyst_annual+cost_gas_comp_anuual+cooling_annual; %$/year

                profitbt=profit_flow_stream-utilties_cost; %$/yr
                operating_profit_per_year=profit_flow_stream-utilties_cost; %$/yr

                FC=cost_cstr+cost_heat_ex+cost_gas_comp+sep_cost;

                %assume WC is 2 months of product and raw material costs
                WC=utilties_cost./6+(pd3.*0.6./6+fa_fresh*.3./6+(fh1+fh2+fh3).*.00112./6)*3600*8400;

                %the factor of 8 relates the base price of equipment to the total cost
                %including installation, indirect costs and contigency of 25%
                TI=2.5*FC*1.1+WC+FC;

                %assuming construciton time is 2 years with cost split
                %evenly between two years
                a_startup=0.1; %assume start up cost is 0.1 of fixed cap
                CR=0.06787; %wall street journal plus 2%
                SU=TI*a_startup;

                n_ops=10;  %years running time
                ER=.12;  % enterprise rate defined in the problem spec 
                FR=.04787;  %Finance rate from wall street journal
                SV=.03*TI; %salvage value is 3% of FC
                TR=.48; %tax rate
                TR_c=1-TR;%tax rate compliment 1-TR
                dist1_D=(TI+SU)/-10;


                %%%Cash flows%%%
                
                %year -1 and -2 cashflows
                cfyn1=(-1)*TI*0.5*(1+CR)^1;
                cfyn0=(-1)*(TI*0.5+WC+SU);
                
                %total capital investment
                TCI=(cfyn1+cfyn0)*-1;
                
                %yearly cashflows
                cfy1=((profitbt*.7-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-1;
                cfy2=((profitbt*.9-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-2;
                cfy3=((profitbt-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-3;
                cfy4=((profitbt-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-4;
                cfy5=((profitbt-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-5;
                cfy6=((profitbt-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-6;
                cfy7=((profitbt-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-7;
                cfy8=((profitbt-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-8;
                cfy9=((profitbt-TCI*FR+dist1_D)*TR_c-(dist1_D))*(1+ER)^-9;
                cfy10=((profitbt-TCI*FR+dist1_D+SV)*TR_c-(dist1_D)+WC+SV-TCI)*(1+ER)^-10;

                %KEY economic indicators
                ROIbt=profitbt/TCI;
                NPV_proj=cfy1+cfy2+cfy3+cfy4+cfy5+cfy6+cfy7+cfy8+cfy9+cfy10;
                NPV_percent=(NPV_proj./TCI)/n_ops;
                POT=TI./(TR_c.*(profitbt-dist1_D)+dist1_D);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%% END ECONOMICS %%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

                %arrays for tracking iterations and later plotting
                graphing_S(k,iter)=S;
                graphing_X(k,iter)=X;
                graphing_pa4(k,iter)=pa4;
                graphing_pd4(k,iter)=pd4;
                graphing_pu4(k,iter)=pu4;
                graphing_ph4(k,iter)=ph4;
                graphing_pw4(k,iter)=pw4;
                graphing_pch44(k,iter)=pch44;
                graphing_V1(k,iter)=V1;
                graphing_V2(k,iter)=V2;
                graphing_V3(k,iter)=V3;
                graphing_V4(k,iter)=V4;
                graphing_fa_fresh(k,iter)=fa_fresh;
                graphing_pa4(k,iter)=pa4;
                graphing_q01(k,iter)=q01;
                graphing_q4(k,iter)=q4;
                graphing_profitbt(k,iter)=profitbt;
                graphing_TCI(k,iter)=TCI;
                graphing_ROIbt(k,iter)=ROIbt;
                graphing_NPV_percent(k,iter)=NPV_percent;
                graphing_NPV_proj(k,iter)=NPV_proj;
                graphing_ph21(k,iter)=ph21;
                graphing_ccat1(k,iter)=ccat1;
                
                
                %write data to file
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(NPV_proj),num2str(NPV_percent),num2str(TI),num2str(TCI),num2str(POT),num2str(X),num2str(S),num2str(ccat1),num2str(ph21),num2str(tau1),num2str(ca1),num2str(cu1),num2str(cw1),num2str(cd1),num2str(ra1),num2str(rb1),num2str(x1),num2str(s1),num2str(pd1),num2str(fa1),num2str(q01),num2str(rho01),num2str(rho1),num2str(q1),num2str(liq_mol_flow_rate_out1),num2str(V1),num2str(h2_mol_frac_liq1),num2str(ph1),num2str(h2_reacted1),num2str(fh1),num2str(fch41),num2str(pch41),num2str(pressure_ch41),num2str(p_total1),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn1),num2str(ROIbt),num2str(profitbt),num2str(pa1),num2str(pu1),num2str(pw1),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(diameter_flash),num2str(height_flash),num2str(height_flash),num2str(dist1_area),num2str(dist1_height),num2str(dist2_area),num2str(dist2_height),num2str(volume_flash),num2str(dist1_r),num2str(dist1_s),num2str(dist1_N),num2str(dist2_r),num2str(dist2_s),num2str(dist2_N_real));
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(NPV_proj),num2str(NPV_percent),num2str(TI),num2str(TCI),num2str(POT),num2str(X),num2str(S),num2str(ccat2),num2str(ph22),num2str(tau2),num2str(ca2),num2str(cu2),num2str(cw2),num2str(cd2),num2str(ra2),num2str(rb2),num2str(x2),num2str(s2),num2str(pd2),num2str(fa2),num2str(q02),num2str(rho02),num2str(rho2),num2str(q2),num2str(liq_mol_flow_rate_out2),num2str(V2),num2str(h2_mol_frac_liq2),num2str(ph2),num2str(h2_reacted2),num2str(fh2),num2str(fch42),num2str(pch42),num2str(pressure_ch42),num2str(p_total2),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn2),num2str(ROIbt),num2str(profitbt),num2str(pa2),num2str(pu2),num2str(pw2),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(diameter_flash),num2str(height_flash),num2str(height_flash),num2str(dist1_area),num2str(dist1_height),num2str(dist2_area),num2str(dist2_height),num2str(volume_flash),num2str(dist1_r),num2str(dist1_s),num2str(dist1_N),num2str(dist2_r),num2str(dist2_s),num2str(dist2_N_real));
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(NPV_proj),num2str(NPV_percent),num2str(TI),num2str(TCI),num2str(POT),num2str(X),num2str(S),num2str(ccat3),num2str(ph23),num2str(tau3),num2str(ca3),num2str(cu3),num2str(cw3),num2str(cd3),num2str(ra3),num2str(rb3),num2str(x3),num2str(s3),num2str(pd3),num2str(fa3),num2str(q03),num2str(rho03),num2str(rho3),num2str(q3),num2str(liq_mol_flow_rate_out3),num2str(V3),num2str(h2_mol_frac_liq3),num2str(ph3),num2str(h2_reacted3),num2str(fh3),num2str(fch43),num2str(pch43),num2str(pressure_ch43),num2str(p_total3),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn3),num2str(ROIbt),num2str(profitbt),num2str(pa3),num2str(pu3),num2str(pw3),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(diameter_flash),num2str(height_flash),num2str(height_flash),num2str(dist1_area),num2str(dist1_height),num2str(dist2_area),num2str(dist2_height),num2str(volume_flash),num2str(dist1_r),num2str(dist1_s),num2str(dist1_N),num2str(dist2_r),num2str(dist2_s),num2str(dist2_N_real));
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(NPV_proj),num2str(NPV_percent),num2str(TI),num2str(TCI),num2str(POT),num2str(X),num2str(S),num2str(ccat4),num2str(ph24),num2str(tau4),num2str(ca4),num2str(cu4),num2str(cw4),num2str(cd4),num2str(ra4),num2str(rb4),num2str(x4),num2str(s4),num2str(pd4),num2str(fa4),num2str(q04),num2str(rho04),num2str(rho4),num2str(q4),num2str(liq_mol_flow_rate_out4),num2str(V4),num2str(h2_mol_frac_liq4),num2str(ph4),num2str(h2_reacted4),num2str(fh4),num2str(fch44),num2str(pch44),num2str(pressure_ch44),num2str(p_total4),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn4),num2str(ROIbt),num2str(profitbt),num2str(pa3),num2str(pu4),num2str(pw4),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(diameter_flash),num2str(height_flash),num2str(height_flash),num2str(dist1_area),num2str(dist1_height),num2str(dist2_area),num2str(dist2_height),num2str(volume_flash),num2str(dist1_r),num2str(dist1_s),num2str(dist1_N),num2str(dist2_r),num2str(dist2_s),num2str(dist2_N_real));
                fprintf(fid1,'\n');
            end
            end
            end
        end
        end
        end
    end
    end
    end
end
end
end
%close file used for logging
fclose(fid1);

%%%% GRAPHS %%%%

% %x vs s
% figure
% hold on;
% title('Selectivity vs. reactor conversion of IBAP')
% xlabel('x')
% ylabel('s')
% x_axis=graphing_X(1,:);
% y_axis=graphing_S(1,:);
% plot(x_axis,y_axis)
% hold off;
% 
% %x vs V
% figure
% hold on;
% title('Reactor Volume in Liters vs. reactor conversion of IBAP')
% xlabel('x')
% ylabel('Reactor Volume (L)')
% x_axis=graphing_X(1,1:50);
% y_axis=graphing_V1(1,1:50)+graphing_V2(1,1:50)+graphing_V3(1,1:50)+graphing_V4(1,1:50);
% plot(x_axis,y_axis)
% hold off;
% 
% %x vs fa_fresh
% figure
% hold on;
% title('Fresh Feed IBAP (mol/s) vs. reactor conversion of IBAP')
% xlabel('x')
% ylabel('fresh feed IBAP (mol/s)')
% x_axis=graphing_X(1,1:50);
% y_axis=graphing_fa_fresh(1,1:50);
% plot(x_axis,y_axis)
% hold off;
%  
% %x vs pa
% figure
% hold on;
% title('Recycle Flow Rate IBAP (mol/s) vs. reactor conversion of IBAP')
% xlabel('x')
% ylabel('Recycle Flow Rate IBAP (mol/s)')
% x_axis=graphing_X(1,10:100);
% y_axis=graphing_pa4(1,10:100);
% plot(x_axis,y_axis)
% hold off;
% 
% %x vs q0
% figure
% hold on;
% title('Total Flow Rate In (L/s) vs. reactor conversion of IBAP')
% xlabel('x')
% ylabel('Total Flow Rate In (L/s)')
% x_axis=graphing_X(1,1:55);
% y_axis=graphing_q01(1,1:55);
% plot(x_axis,y_axis)
% hold off;
% 
% %x vs q
% figure
% hold on;
% title('Total Flow Rate To Seperation System (L/s) vs. reactor conversion of IBAP')
% xlabel('x')
% ylabel('Total Flow Rate To Seperation System (L/s)')
% x_axis=graphing_X(1,1:55);
% y_axis=graphing_q4(1,1:55);
% plot(x_axis,y_axis)
% hold off;
%  
% %x vs mol fract
% figure
% hold on;
% title('mole fraction vs. reactor conversion of IBAP')
% xlabel('x')
% ylabel('mole fraction')
% x_axis=graphing_X(1,:);
% y_axis=graphing_pa4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
% plot(x_axis,y_axis,'y')
% x_axis=graphing_X(1,:);
% y_axis=graphing_pd4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
% plot(x_axis,y_axis,'g')
% x_axis=graphing_X(1,:);
% y_axis=graphing_pu4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
% plot(x_axis,y_axis,'r')
% x_axis=graphing_X(1,:);
% y_axis=graphing_ph4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
% plot(x_axis,y_axis,'c')
% x_axis=graphing_X(1,:);
% y_axis=graphing_pw4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
% plot(x_axis,y_axis,'b')
% x_axis=graphing_X(1,:);
% y_axis=graphing_pch44(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
% plot(x_axis,y_axis,'k')
% legend('IBAP','IBPE','IBEB','hydrogen','water','methane')
% hold off;
% 
% %x vs profitbt, TCI
% figure
% hold on;
% title('Profit before tax and TCI vs. total conversion')
% xlabel('X')
% ylabel('Millions of $')
% x_axis=graphing_X(1,20:65)./1000000;
% y_axis=graphing_profitbt(1,20:65)./1000000;
% plot(x_axis,y_axis,'g')
% x_axis=graphing_X(1,10:75)./1000000;
% y_axis=graphing_TCI(1,10:75)./1000000;
% plot(x_axis,y_axis,'r')
% legend('Profit before tax','TCI')
% hold off;
% 
% %x vs ROIbt
% figure
% hold on;
% title('Return on Investment before taxes adn NPV% vs. total conversion')
% xlabel('x')
% ylabel('%')
% x_axis=graphing_X(1,:);
% y_axis=graphing_ROIbt(1,:)*100;
% plot(x_axis,y_axis,'g')
% x_axis=graphing_X(1,:);
% y_axis=graphing_NPV_percent(1,:)*100;
% plot(x_axis,y_axis,'b')
% legend('ROI bt','NPV%')
% ylim([0 max(graphing_ROIbt(1,:)*100)+10])
% hold off;
% 
% %x vs nvp_proj
% figure
% hold on;
% title('NPV project vs. total conversion')
% xlabel('x')
% ylabel('Millions of $')
% x_axis=graphing_X(1,:);
% y_axis=graphing_NPV_proj(1,:)./1000000;
% plot(x_axis,y_axis,'g')
% legend('NPV project')
% ylim([0 max(graphing_NPV_proj(1,:)./1000000).*1.1])
% hold off;

% %h2 vs profitbt, TCI
% figure
% hold on;
% title('Profit before tax and TCI vs. pressure of hydrogen in bar')
% xlabel('hydrogen pressure (bar)')
% ylabel('$')
% x_axis=1.01325.*graphing_ph21(1,:);
% y_axis=graphing_profitbt(1,:);
% plot(x_axis,y_axis,'g')
% x_axis=1.01325.*graphing_ph21(1,:);
% y_axis=graphing_TCI(1,:);
% plot(x_axis,y_axis,'r')
% legend('Profit before tax','TCI')
% hold off;
% 
% %h2 vs ROIbt
% figure
% hold on;
% title('Return on Investment before taxes adn NPV% vs. pressure of hydrogen in bar')
% xlabel('hydrogen pressure (bar)')
% ylabel('ROIbt')
% x_axis=1.01325.*graphing_ph21(1,:);
% y_axis=graphing_ROIbt(1,:);
% plot(x_axis,y_axis,'g')
% x_axis2=1.01325.*graphing_ph21(1,:);
% y_axis2=graphing_NPV_percent(1,:);
% plot(x_axis2,y_axis2,'b')
% legend('ROI bt','NPV%')
% hold off;
% 
% %h2 vs nvp_proj
% figure
% hold on;
% title('NPV project vs. pressure of hydrogen in bar')
% xlabel('hydrogen pressure (bar)')
% ylabel('$')
% x_axis=1.01325*graphing_ph21(1,:);
% y_axis=graphing_NPV_proj(1,:);
% plot(x_axis,y_axis,'g')
% hold off;
 
% %ccat vs profitbt, TCI
% figure
% hold on;
% title('Profit before tax and TCI vs. catalyst concentration in g/L')
% xlabel('g/L')
% ylabel('$')
% x_axis=graphing_ccat1(1,:);
% y_axis=graphing_profitbt(1,:);
% plot(x_axis,y_axis,'g')
% x_axis=graphing_ccat1(1,:);
% y_axis=graphing_TCI(1,:);
% plot(x_axis,y_axis,'r')
% legend('Profit before tax','TCI')
% xlim([.26,.78])
% hold off;
% 
% %ccat vs ROIbt
% figure
% hold on;
% title('Return on Investment before taxes adn NPV% vs. catalyst concentration in g/L')
% xlabel('g/L')
% ylabel('ROIbt')
% x_axis=graphing_ccat1(1,:);
% y_axis=graphing_ROIbt(1,:);
% plot(x_axis,y_axis,'g')
% x_axis=graphing_ccat1(1,:);
% y_axis=graphing_NPV_percent(1,:);
% plot(x_axis,y_axis,'b')
% legend('ROI bt','NPV%')
% xlim([.26,.78])
% hold off;
% 
% 
% %ccat vs nvp_proj
% figure
% hold on;
% title('NPV project vs. catalyst concentration in g/L')
% xlabel('g/L')
% ylabel('$')
% x_axis=graphing_ccat1(1,:);
% y_axis=graphing_NPV_proj(1,:);
% plot(x_axis,y_axis,'g')
% legend('NPV project')
% xlim([.26,.78])
% hold off;

%%%% END GRAPHS %%%%

toc
fprintf('Complete.\n\n\n')

%desgin level 3 mwc last edited 4/29/2014 5:24pm