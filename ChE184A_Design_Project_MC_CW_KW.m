%Chemical Engineering 184A Design Project by Matt Carroll, Conor Waldron,
%and Kyle Wu

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
fprintf(fid1,'X,S,ccat,ph2,tau,ca,cu,cw,cd,ra,rb,x,s,pd,fa,q0,rho0,rho,q,liq_mol_flow_rate_out,V,h2_mol_frac_liq,ph,h2_reacted,fh,fch4,pch4,pressure_ch4,p_total,operating_profit_per_year,deltaHa,deltaHb,heat_of_rxn,TI,TCI,ROIbt,profitbt,pa,pu,pw,fa_fresh,cost_gas_comp,cost_heat_ex,cost_cstr,FC,WC,utilties_cost,profit_flow_stream,cost_cstr,cost_heat_ex,cost_gas_comp,sep_cost,FC,NPV_proj,NPV_percent,POT,\n');

%loop control iterator varibles and number of iterations
nph2=1000;
nccat=1000;
ntau=10000;
iter=0;
k=1;


%reactor 1
for ph21=9.87 %linspace(9.87,39.47,nph2) %ph2 range
for ccat1=.78 %linspace(.26,.78,nccat) %ccat range
for tau1=2236.036;%logspace(1,7,ntau)
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
                S=cd4./(ca0-ca4);

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
                pd4=((5000000.*1000)./178.3)/(8400*60*60); %mol/s

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
                fa1=pa4+fa_fresh; %mol/s
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
                %%%%%%%%%%%%%% ECONOMICS %%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %fixed capital cost: marshall swift index current value
                ms=1550; 

                %all cost are based intastlatioin costs

                %Cstr costs (4 reactors)
                fcstr_pressure=0.00130*(p_total1*14.6959488)+0.84;
                fcstr_material=1.3; %brass
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V1/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr1=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $

                fcstr_pressure=0.00130*(p_total2*14.6959488)+0.84;
                fcstr_material=1.3; %brass
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V2/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr2=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $

                fcstr_pressure=0.00130*(p_total3*14.6959488)+0.84;
                fcstr_material=1.3; %brass
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V3/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr3=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $

                fcstr_pressure=0.00130*(p_total4*14.6959488)+0.84;
                fcstr_material=1.3; %brass
                fcstr=fcstr_material+fcstr_pressure; %design factor (dependant on material)
                diameter_m=((2.*(1.13.*V4/1000))/pi).^(1./3); %diameter in meters
                diameter=diameter_m*3.28084; %diamteer in ft.
                cost_cstr4=(ms/280).*(101.9*(diameter^1.066)*((2*diameter)^0.82)*(2.18+fcstr)); % $
                
                %total reactor cost
                cost_cstr=cost_cstr1+cost_cstr2+cost_cstr3+cost_cstr4; % $

                %Heat Exchangers (4, one for each reactor)
                V_feet=1.13.*V1*0.0353147;
                area=V_feet/(pi*0.25*diameter^2); %surface area of reactor
                fheat_material=1.3; %brass
                fheat_pressure=0; %factor for steam pressure assuming liq h2o <150psi
                fheat_design=0.85; %factor
                cost_heat_ex1=(ms/280).*(101.3*((fheat_pressure+fheat_design)*fheat_material+2.29)*area^0.65); % $

                V_feet=1.13.*V2*0.0353147;
                area=V_feet/(pi*0.25*diameter^2); %surface area of reactor
                fheat_material=1.3; %brass
                fheat_pressure=0; %factor for steam pressure assuming liq h2o <150psi
                fheat_design=0.85; %factor
                cost_heat_ex2=(ms/280).*(101.3*((fheat_pressure+fheat_design)*fheat_material+2.29)*area^0.65); % $

                V_feet=1.13.*V3*0.0353147;
                area=V_feet/(pi*0.25*diameter^2); %surface area of reactor
                fheat_material=1.3; %brass
                fheat_pressure=0; %factor for steam pressure assuming liq h2o <150psi
                fheat_design=0.85; %factor
                cost_heat_ex3=(ms/280).*(101.3*((fheat_pressure+fheat_design)*fheat_material+2.29)*area^0.65); % $

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
                cost_gas_comp=(ms/280)*517.5*(fgas+2.11)*bhp^0.82; % $

                %separating cost
                sep_cost=1000000/X; % $
                
                %%%%%%%%%%%%%%%%%%%RECURNING UTLITIES COST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Gas compressor running cost
                bhp2=hp./.6; % Guiethei
                cost_gas_comp_anuual=bhp2.*8400*0.7456*.00588; % $/yr

                %reactor cooling cost
                cooling_annual=-(heat_of_rxn1+heat_of_rxn2+heat_of_rxn3+heat_of_rxn4).*8400*.00588; %$/yr
                
                %annual profit from flows
                profit_flow_stream=(pd4*.6-fa_fresh*.3-(fh1+fh2+fh3+fh4).*.00112)*(3600*8400); %$/yr

                %annual seperator cost (given)
                seperation_annual=100000./X; %$/yr
                
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
                TI=2.5*FC;

                %assuming construciton time is 2 years with cost split
                %evenly between two years
                a_startup=0.1; %assume start up cost is 0.1 of fixed cap
                CR=0.04787; %wall street journal plus 2%
                SU=TI*a_startup;

                n_ops=10;  %years running time
                ER=.12;  % enterprise rate defined in the problem spec 
                FR=.02787;  %Finance rate from wall street journal
                SV=.03*TI; %salvage value is 3% of FC
                TR=.48; %tax rate
                TR_c=1-TR;%tax rate compliment 1-TR
                D=(TI+SU)/-10;


                %%%Cash flows%%%
                
                %year -1 and -2 cashflows
                cfyn1=(-1)*TI*0.5*(1+CR)^1;
                cfyn0=(-1)*(TI*0.5+WC+SU);
                
                %total capital investment
                TCI=(cfyn1+cfyn0)*-1;
                
                %yearly cashflows
                cfy1=((profitbt*.7-TCI*FR+D)*TR_c-(D))*(1+ER)^-1;
                cfy2=((profitbt*.9-TCI*FR+D)*TR_c-(D))*(1+ER)^-2;
                cfy3=((profitbt-TCI*FR+D)*TR_c-(D))*(1+ER)^-3;
                cfy4=((profitbt-TCI*FR+D)*TR_c-(D))*(1+ER)^-4;
                cfy5=((profitbt-TCI*FR+D)*TR_c-(D))*(1+ER)^-5;
                cfy6=((profitbt-TCI*FR+D)*TR_c-(D))*(1+ER)^-6;
                cfy7=((profitbt-TCI*FR+D)*TR_c-(D))*(1+ER)^-7;
                cfy8=((profitbt-TCI*FR+D)*TR_c-(D))*(1+ER)^-8;
                cfy9=((profitbt-TCI*FR+D)*TR_c-(D))*(1+ER)^-9;
                cfy10=((profitbt-TCI*FR+D+SV)*TR_c-(D)+WC+SV-TCI)*(1+ER)^-10;

                %KEY economic indicators
                ROIbt=profitbt/TCI;
                NPV_proj=cfy1+cfy2+cfy3+cfy4+cfy5+cfy6+cfy7+cfy8+cfy9+cfy10;
                NPV_percent=(NPV_proj./TCI)/n_ops;
                POT=TI./(TR_c.*(profitbt-D)+D);

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
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(X),num2str(S),num2str(ccat1),num2str(ph21),num2str(tau1),num2str(ca1),num2str(cu1),num2str(cw1),num2str(cd1),num2str(ra1),num2str(rb1),num2str(x1),num2str(s1),num2str(pd1),num2str(fa1),num2str(q01),num2str(rho01),num2str(rho1),num2str(q1),num2str(liq_mol_flow_rate_out1),num2str(V1),num2str(h2_mol_frac_liq1),num2str(ph1),num2str(h2_reacted1),num2str(fh1),num2str(fch41),num2str(pch41),num2str(pressure_ch41),num2str(p_total1),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn1),num2str(TI),num2str(TCI),num2str(ROIbt),num2str(profitbt),num2str(pa1),num2str(pu1),num2str(pw1),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(NPV_proj),num2str(NPV_percent),num2str(POT));
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(X),num2str(S),num2str(ccat2),num2str(ph22),num2str(tau2),num2str(ca2),num2str(cu2),num2str(cw2),num2str(cd2),num2str(ra2),num2str(rb2),num2str(x2),num2str(s2),num2str(pd2),num2str(fa2),num2str(q02),num2str(rho02),num2str(rho2),num2str(q2),num2str(liq_mol_flow_rate_out2),num2str(V2),num2str(h2_mol_frac_liq2),num2str(ph2),num2str(h2_reacted2),num2str(fh2),num2str(fch42),num2str(pch42),num2str(pressure_ch42),num2str(p_total2),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn2),num2str(TI),num2str(TCI),num2str(ROIbt),num2str(profitbt),num2str(pa2),num2str(pu2),num2str(pw2),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(NPV_proj),num2str(NPV_percent),num2str(POT));
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(X),num2str(S),num2str(ccat3),num2str(ph23),num2str(tau3),num2str(ca3),num2str(cu3),num2str(cw3),num2str(cd3),num2str(ra3),num2str(rb3),num2str(x3),num2str(s3),num2str(pd3),num2str(fa3),num2str(q03),num2str(rho03),num2str(rho3),num2str(q3),num2str(liq_mol_flow_rate_out3),num2str(V3),num2str(h2_mol_frac_liq3),num2str(ph3),num2str(h2_reacted3),num2str(fh3),num2str(fch43),num2str(pch43),num2str(pressure_ch43),num2str(p_total3),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn3),num2str(TI),num2str(TCI),num2str(ROIbt),num2str(profitbt),num2str(pa3),num2str(pu3),num2str(pw3),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(NPV_proj),num2str(NPV_percent),num2str(POT));
                fprintf(fid1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',num2str(X),num2str(S),num2str(ccat4),num2str(ph24),num2str(tau4),num2str(ca4),num2str(cu4),num2str(cw4),num2str(cd4),num2str(ra4),num2str(rb4),num2str(x4),num2str(s4),num2str(pd4),num2str(fa4),num2str(q04),num2str(rho04),num2str(rho4),num2str(q4),num2str(liq_mol_flow_rate_out4),num2str(V4),num2str(h2_mol_frac_liq4),num2str(ph4),num2str(h2_reacted4),num2str(fh4),num2str(fch44),num2str(pch44),num2str(pressure_ch44),num2str(p_total4),num2str(operating_profit_per_year),num2str(deltaHa),num2str(deltaHb),num2str(heat_of_rxn4),num2str(TI),num2str(TCI),num2str(ROIbt),num2str(profitbt),num2str(pa3),num2str(pu4),num2str(pw4),num2str(fa_fresh),num2str(cost_gas_comp),num2str(cost_heat_ex),num2str(cost_cstr),num2str(FC),num2str(WC),num2str(utilties_cost),num2str(profit_flow_stream),num2str(cost_cstr),num2str(cost_heat_ex),num2str(cost_gas_comp),num2str(sep_cost),num2str(FC),num2str(NPV_proj),num2str(NPV_percent),num2str(POT));
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

%x vs s
figure
hold on;
title('Selectivity vs. reactor conversion of IBAP')
xlabel('x')
ylabel('s')
x_axis=graphing_X(1,:);
y_axis=graphing_S(1,:);
plot(x_axis,y_axis)
hold off;

%x vs V
figure
hold on;
title('Reactor Volume in Liters vs. reactor conversion of IBAP')
xlabel('x')
ylabel('Reactor Volume (L)')
x_axis=graphing_X(1,:);
y_axis=graphing_V1(1,:)+graphing_V2(1,:)+graphing_V3(1,:)+graphing_V4(1,:);
plot(x_axis,y_axis)
hold off;

%x vs fa_fresh
figure
hold on;
title('Fresh Feed IBAP (mol/s) vs. reactor conversion of IBAP')
xlabel('x')
ylabel('fresh feed IBAP (mol/s)')
x_axis=graphing_X(1,:);
y_axis=graphing_fa_fresh(1,:);
plot(x_axis,y_axis)
hold off;

%x vs pa
figure
hold on;
title('Recycle Flow Rate IBAP (mol/s) vs. reactor conversion of IBAP')
xlabel('x')
ylabel('Recycle Flow Rate IBAP (mol/s)')
x_axis=graphing_X(1,:);
y_axis=graphing_pa4(1,:);
plot(x_axis,y_axis)
hold off;

%x vs q0
figure
hold on;
title('Total Flow Rate In (L/s) vs. reactor conversion of IBAP')
xlabel('x')
ylabel('Total Flow Rate In (L/s)')
x_axis=graphing_X(1,:);
y_axis=graphing_q01(1,:);
plot(x_axis,y_axis)
hold off;

%x vs q
figure
hold on;
title('Total Flow Rate To Seperation System (L/s) vs. reactor conversion of IBAP')
xlabel('x')
ylabel('Total Flow Rate To Seperation System (L/s)')
x_axis=graphing_X(1,:);
y_axis=graphing_q4(1,:);
plot(x_axis,y_axis)
hold off;

%x vs mol fract
figure
hold on;
title('mole fraction vs. reactor conversion of IBAP')
xlabel('x')
ylabel('mole fraction')
x_axis=graphing_X(1,:);
y_axis=graphing_pa4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
plot(x_axis,y_axis,'y')
x_axis=graphing_X(1,:);
y_axis=graphing_pd4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
plot(x_axis,y_axis,'g')
x_axis=graphing_X(1,:);
y_axis=graphing_pu4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
plot(x_axis,y_axis,'r')
x_axis=graphing_X(1,:);
y_axis=graphing_ph4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
plot(x_axis,y_axis,'c')
x_axis=graphing_X(1,:);
y_axis=graphing_pw4(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
plot(x_axis,y_axis,'b')
x_axis=graphing_X(1,:);
y_axis=graphing_pch44(1,:)./(graphing_pa4(1,:)+graphing_pd4(1,:)+graphing_pu4(1,:)+graphing_ph4(1,:)+graphing_pw4(1,:)+graphing_pch44(1,:));
plot(x_axis,y_axis,'k')
legend('IBAP','IBPE','IBEB','hydrogen','water','methane')
hold off;


%x vs profitbt, TCI
figure
hold on;
title('Profit before tax and TCI vs. total conversion')
xlabel('X')
ylabel('$')
x_axis=graphing_X(1,:);
y_axis=graphing_profitbt(1,:);
plot(x_axis,y_axis,'g')
x_axis=graphing_X(1,:);
y_axis=graphing_TCI(1,:);
plot(x_axis,y_axis,'r')
legend('Profit before tax','TCI')
hold off;

%x vs ROIbt
figure
hold on;
title('Return on Investment before taxes adn NPV% vs. total conversion')
xlabel('x')
ylabel('ROIbt')
x_axis=graphing_X(1,:);
y_axis=graphing_ROIbt(1,:);
plot(x_axis,y_axis,'g')
x_axis=graphing_X(1,:);
y_axis=graphing_NPV_percent(1,:);
plot(x_axis,y_axis,'b')
legend('ROI bt','NPV%')
hold off;

%x vs nvp_proj
figure
hold on;
title('NPV project vs. total conversion')
xlabel('x')
ylabel('$')
x_axis=graphing_X(1,:);
y_axis=graphing_NPV_proj(1,:);
plot(x_axis,y_axis,'g')
legend('NPV project')
hold off;

%x vs profitbt, TCI
figure
hold on;
title('Profit before tax and TCI vs. pressure of hydrogen in bar')
xlabel('hydrogen pressure (bar)')
ylabel('$')
x_axis=1.01325.*graphing_ph21(1,:);
y_axis=graphing_profitbt(1,:);
plot(x_axis,y_axis,'g')
x_axis=1.01325.*graphing_ph21(1,:);
y_axis=graphing_TCI(1,:);
plot(x_axis,y_axis,'r')
legend('Profit before tax','TCI')
hold off;

%x vs ROIbt
figure
hold on;
title('Return on Investment before taxes adn NPV% vs. pressure of hydrogen in bar')
xlabel('hydrogen pressure (bar)')
ylabel('ROIbt')
x_axis=1.01325.*graphing_ph21(1,:);
y_axis=graphing_ROIbt(1,:);
plot(x_axis,y_axis,'g')
x_axis=1.01325.*graphing_ph21(1,:);
y_axis=graphing_NPV_percent(1,:);
plot(x_axis,y_axis,'b')
legend('ROI bt','NPV%')
hold off;

%x vs nvp_proj
figure
hold on;
title('NPV project vs. pressure of hydrogen in bar')
xlabel('hydrogen pressure (bar)')
ylabel('$')
x_axis=1.01325*graphing_ph21(1,:);
y_axis=graphing_NPV_proj(1,:);
plot(x_axis,y_axis,'g')
legend('NPV project')
hold off;

%x vs profitbt, TCI
figure
hold on;
title('Profit before tax and TCI vs. catalyst concentration in g/L')
xlabel('g/L')
ylabel('$')
x_axis=graphing_ccat1(1,:);
y_axis=graphing_profitbt(1,:);
plot(x_axis,y_axis,'g')
x_axis=graphing_ccat1(1,:);
y_axis=graphing_TCI(1,:);
plot(x_axis,y_axis,'r')
legend('Profit before tax','TCI')
hold off;

%x vs ROIbt
figure
hold on;
title('Return on Investment before taxes adn NPV% vs. catalyst concentration in g/L')
xlabel('g/L')
ylabel('ROIbt')
x_axis=graphing_ccat1(1,:);
y_axis=graphing_ROIbt(1,:);
plot(x_axis,y_axis,'g')
x_axis=graphing_ccat1(1,:);
y_axis=graphing_NPV_percent(1,:);
plot(x_axis,y_axis,'b')
legend('ROI bt','NPV%')
hold off;


%x vs nvp_proj
figure
hold on;
title('NPV project vs. catalyst concentration in g/L')
xlabel('g/L')
ylabel('$')
x_axis=graphing_ccat1(1,:);
y_axis=graphing_NPV_proj(1,:);
plot(x_axis,y_axis,'g')
legend('NPV project')
hold off;

%%%% END GRAPHS %%%%

toc
fprintf('Complete.\n\n\n')

%desgin level 3 mwc last edited 3/10/2014 5:18pm