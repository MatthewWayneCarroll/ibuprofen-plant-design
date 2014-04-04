%Solve system of non-linear algebraic equations



function fcns = CSTR_bal1(w)
ca = w(1);
cu = w(2);
cw = w(3);
cd = w(4);

global tau1
global ca0
global ccat1
global ph21

%kinetic constants
k1=1.14;
k2=.095;
ka=76.4;
kh=141;
kw=529;

%rate expressions
ra=(ccat1*k1*ca*ph21)./(1+ka*ca+(kh*ph21).^(1./2)+kw*cw).^(2);
rb=(ccat1*k2*cd*ph21)./(1+ka*ca+(kh*ph21).^(1./2)+kw*cw).^(2);

fcns(1) =ca0-ca-tau1.*ra;
fcns(2) =cu-tau1.*rb;
fcns(3) =cw-tau1.*rb;
fcns(4) =cd-tau1.*ra+tau1.*rb;
end

    