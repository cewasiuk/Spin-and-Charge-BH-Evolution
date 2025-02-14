(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"];


(* ::Input::Initialization:: *)
nbmodes = 30;
nbener = 200;
M = 0.5;
s = 1;
nba = 50;
a = Table[0,{i,1,nba}];
For[i=30,i<nba+1,i++,a[[50+30-i]] = 1-10^(Log10[0.0001]+(Log10[0.1]-Log10[0.0001])/(nba-30)*(i-30))];
For[i=2,i<10,i++,a[[i]] = 10^(Log10[0.001]+(Log10[0.1]-Log10[0.001])/7*(i-2))];
For[i=9,i<31,i++,a[[i]]=0.1+(0.9-0.1)/21*(i-9)];
a
(*a = Table[0.9999,{i,1,1}]*)
omega = Table[10^(Log10[0.01] + (Log10[5.] - Log10[0.01])/(nbener-1)*i),{i,0,nbener-1}]
(*omega = Table[0.17,{i,1,1}]*)
l=Table[1+i,{i,0,nbmodes-1}]
Ainh = 1;


(* ::Input::Initialization:: *)
H[l_,alpha_,beta_,s_]:=If[l==0 || l==0.5 ||l==-0.5,0,(l^2-(alpha+beta)^2/4)*(l^2-s^2)*(l^2-(alpha-beta)^2/4)/(2*(l-1/2)*l^3*(l+1/2))];
coef0[l_]:=l*(l+1);
coef1[l_,m_,s_]:=If[l==0,0,-2*s^2*m/(l*(l+1.))];
coef2[l_,m_,s_]:=Module[{alpha=Abs[m+s],beta=Abs[m-s]},H[l+1,alpha,beta,s]-H[l,alpha,beta,s]-1];
coef3[l_,m_,s_]:=Module[{alpha=Abs[m+s],beta=Abs[m-s]},If[l==0,0,If[l==1,-2*s^2*m*H[l+1,alpha,beta,s]/(l*(l+1)^2*(l+2)),2*s^2*m*(H[l,alpha,beta,s]/((l-1)*l^2*(l+1))-H[l+1,alpha,beta,s]/(l*(l+1)^2*(l+2)))]]];
coef4[l_,m_,s_]:=Module[{alpha=Abs[m+s],beta=Abs[m-s]},If[l==0,1/2*H[l+1,alpha,beta,s]^2/(l+1)
-1/4*(l+2)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/((l+1)*(l+3/2)),If[l==1,4*s^4*m^2*H[l+1,alpha,beta,s]/(l^2*(l+1)^4*(l+2)^2)
+1/2*(H[l+1,alpha,beta,s]^2/(l+1)+H[l+1,alpha,beta,s]*H[l,alpha,beta,s]/(l*(l+1))-H[l,alpha,beta,s]^2/l)
+1/4*((l-1)*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]/(l*(l-1/2))-(l+2)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/((l+1)*(l+3/2))),If[l==0.5,4*s^4*m^2*(H[l+1,alpha,beta,s]/(l^2*(l+1)^4*(l+2)^2)-H[l,alpha,beta,s]/((l-1)^2*l^4*(l+1)^2))
+1/2*(H[l+1,alpha,beta,s]^2/(l+1)+H[l+1,alpha,beta,s]*H[l,alpha,beta,s]/(l*(l+1))-H[l,alpha,beta,s]^2/l)
-1/4*(l+2)*H[l+1,alpha,beta,s]*H[l+2.,alpha,beta,s]/((l+1)*(l+3/2)),4*s^4*m^2*(H[l+1,alpha,beta,s]/(l^2*(l+1)^4*(l+2)^2)-H[l,alpha,beta,s]/((l-1)^2*l^4*(l+1)^2))
+1/2*(H[l+1,alpha,beta,s]^2/(l+1)+H[l+1,alpha,beta,s]*H[l,alpha,beta,s]/(l*(l+1))-H[l,alpha,beta,s]^2/l)
+1/4*((l-1)*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]/(l*(l-1/2))-(l+2)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/((l+1)*(l+3/2)))]]]];
coef5[l_,m_,s_]:=Module[{alpha=Abs[m+s],beta=Abs[m-s]},If[l==0,0,If[l==1,-8*s^6*m^3*H[l+1,alpha,beta,s]/(l^3*(l+1)^6*(l+2)^3)
+s^2*m*(-3*H[l+1,alpha,beta,s]^2/(l*(l+1)^3*(l+2))+1/2*(3*l+7)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(l*(l+1)^3*(l+3/2)*(l+3))),If[l==0.5 || l==2,8*s^6*m^3*(H[l,alpha,beta,s]/((l-1)^3*l^6*(l+1)^3)-H[l+1,alpha,beta,s]/(l^3*(l+1)^6*(l+2)^3))
+s^2*m*(3*H[l,alpha,beta,s]^2/((l-1)*l^3*(l+1))-(7*l^2+7*l+4)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]/((l-1)*l^3*(l+1)^3*(l+2))-3*H[l+1,alpha,beta,s]^2/(l*(l+1)^3*(l+2))+1/2*(3*l+7)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(l*(l+1)^3*(l+3/2)*(l+3))),8*s^6*m^3*(H[l,alpha,beta,s]/((l-1)^3*l^6*(l+1)^3)-H[l+1,alpha,beta,s]/(l^3*(l+1)^6*(l+2)^3))
+s^2*m*(3*H[l,alpha,beta,s]^2/((l-1)*l^3*(l+1))-(7*l^2+7*l+4)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]/((l-1)*l^3*(l+1)^3*(l+2))-3*H[l+1,alpha,beta,s]^2/(l*(l+1)^3*(l+2))+1/2*((3*l+7)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(l*(l+1)^3*(l+3/2)*(l+3))-(3*l-4)*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]/((l-2)*(l-1/2)*l^3*(l+1))))]]]];
coef6[l_,m_,s_]:=Module[{alpha=Abs[m+s],beta=Abs[m-s]},If[l==0,1/4*(2*H[l+1,alpha,beta,s]^3/(l+1)^2+(l+2)^2*H[l+2,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)-(l+2)*(7*l+10)*H[l+1,alpha,beta,s]^2*H[l+2,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)+(l+3)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]*H[l+3,alpha,beta,s]/(12*(l+1)*(l+3/2)^2)),If[l==1,16*s^8*m^4*H[l+1,alpha,beta,s]/(l^4*(l+1)^8*(l+2)^4)
+4*s^4*m^2*(3*H[l+1,alpha,beta,s]^2/(l^2*(l+1)^5*(l+2)^2)-1/2*(3*l^2+14*l+17)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(l^2*(l+1)^5*(l+3/2)*(l+2)*(l+3)^2))
+1/4*(2*H[l+1,alpha,beta,s]^3/(l+1)^2+(2*l^2+4*l+3)*H[l,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(l^2*(l+1)^2)-(2*l^2+1)*H[l+1,alpha,beta,s]^2*H[l,alpha,beta,s]/(l^2*(l+1)^2)-2*H[l,alpha,beta,s]^3/l^2+(l+2)*(3*l^2+2*l-3)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(4*l*(l+1)^2*(l+3/2)^2)-(l-1)*(3*l^2+4*l-2)*H[l+1,alpha,beta,s]*H[l,alpha,beta,s]*H[l-1,alpha,beta,s]/(4*(l-1/2)^2*l^2*(l+1))+(l+2)^2*H[l+2,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)-(l-1)^2*H[l-1,alpha,beta,s]^2*H[l,alpha,beta,s]/(4*(l-1/2)^2*l^2)+(l-1)*(7*l-3)*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]^2/(4*(l-1/2)^2*l^2)-(l+2)*(7*l+10)*H[l+1,alpha,beta,s]^2*H[l+2,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)+(l+3)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]*H[l+3,alpha,beta,s]/(12*(l+1)*(l+3/2)^2)-(l-2)*H[l-2,alpha,beta,s]*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]/(12*(l-1/2)^2*l)),If[l==0.5,16*s^8*m^4*(H[l+1,alpha,beta,s]/(l^4*(l+1)^8*(l+2)^4)-H[l,alpha,beta,s]/((l-1)^4*l^8*(l+1)^4))
+4*s^4*m^2*(3*H[l+1,alpha,beta,s]^2/(l^2*(l+1)^5*(l+2)^2)-3*H[l,alpha,beta,s]^2/((l-1)^2*l^5*(l+1)^2)+(11*l^4+22*l^3+31*l^2+20*l+6)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]/((l-1)^2*l^5*(l+1)^5*(l+2)^2)-1/2*(3*l^2+14*l+17)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(l^2*(l+1)^5*(l+3/2)*(l+2)*(l+3)^2))
+1/4*(2*H[l+1,alpha,beta,s]^3/(l+1)^2+(2*l^2+4*l+3)*H[l,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(l^2*(l+1)^2)-(2*l^2+1)*H[l+1,alpha,beta,s]^2*H[l,alpha,beta,s]/(l^2*(l+1)^2)-2*H[l,alpha,beta,s]^3/l^2+(l+2)*(3*l^2+2*l-3)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(4*l*(l+1)^2*(l+3/2)^2)+(l+2)^2*H[l+2,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)-(l+2)*(7*l+10)*H[l+1,alpha,beta,s]^2*H[l+2,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)+(l+3)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]*H[l+3,alpha,beta,s]/(12*(l+1)*(l+3/2)^2)),If[l==2,16*s^8*m^4*(H[l+1,alpha,beta,s]/(l^4*(l+1)^8*(l+2)^4)-H[l,alpha,beta,s]/((l-1)^4*l^8*(l+1)^4))
+4*s^4*m^2*(3*H[l+1,alpha,beta,s]^2/(l^2*(l+1)^5*(l+2)^2)-3*H[l,alpha,beta,s]^2/((l-1)^2*l^5*(l+1)^2)+(11*l^4+22*l^3+31*l^2+20*l+6)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]/((l-1)^2*l^5*(l+1)^5*(l+2)^2)-1/2*(3*l^2+14*l+17)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(l^2*(l+1)^5*(l+3/2)*(l+2)*(l+3)^2))
+1/4*(2*H[l+1,alpha,beta,s]^3/(l+1)^2+(2*l^2+4*l+3)*H[l,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(l^2*(l+1)^2)-(2*l^2+1)*H[l+1,alpha,beta,s]^2*H[l,alpha,beta,s]/(l^2*(l+1)^2)-2*H[l,alpha,beta,s]^3/l^2+(l+2)*(3*l^2+2*l-3)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(4*l*(l+1)^2*(l+3/2)^2)-(l-1)*(3*l^2+4*l-2)*H[l+1,alpha,beta,s]*H[l,alpha,beta,s]*H[l-1,alpha,beta,s]/(4*(l-1/2)^2*l^2*(l+1))+(l+2)^2*H[l+2,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)-(l-1)^2*H[l-1,alpha,beta,s]^2*H[l,alpha,beta,s]/(4*(l-1/2)^2*l^2)+(l-1)*(7*l-3)*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]^2/(4*(l-1/2)^2*l^2)-(l+2)*(7*l+10)*H[l+1,alpha,beta,s]^2*H[l+2,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)+(l+3)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]*H[l+3,alpha,beta,s]/(12*(l+1)*(l+3/2)^2)-(l-2)*H[l-2,alpha,beta,s]*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]/(12*(l-1/2)^2
l)),16*s^8*m^4*(H[l+1,alpha,beta,s]/(l^4*(l+1)^8*(l+2)^4)-H[l,alpha,beta,s]/((l-1)^4*l^8*(l+1)^4))
+4*s^4*m^2*(3*H[l+1,alpha,beta,s]^2/(l^2*(l+1)^5*(l+2)^2)-3*H[l,alpha,beta,s]^2/((l-1)^2*l^5*(l+1)^2)+(11*l^4+22*l^3+31*l^2+20*l+6)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]/((l-1)^2*l^5*(l+1)^5*(l+2)^2)+1/2*((3*l^2-8*l+6)*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]/((l-2)^2*(l-1)*(l-1/2)*l^5*(l+1)^2)-(3*l^2+14*l+17)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(l^2*(l+1)^5*(l+3/2)*(l+2)*(l+3)^2)))
+1/4*(2*H[l+1,alpha,beta,s]^3/(l+1)^2+(2*l^2+4*l+3)*H[l,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(l^2*(l+1)^2)-(2*l^2+1)*H[l+1,alpha,beta,s]^2*H[l,alpha,beta,s]/(l^2*(l+1)^2)-2*H[l,alpha,beta,s]^3/l^2+(l+2)*(3*l^2+2*l-3)*H[l,alpha,beta,s]*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]/(4*l*(l+1)^2*(l+3/2)^2)-(l-1)*(3*l^2+4*l-2)*H[l+1,alpha,beta,s]*H[l,alpha,beta,s]*H[l-1,alpha,beta,s]/(4*(l-1/2)^2*l^2*(l+1))+(l+2)^2*H[l+2,alpha,beta,s]^2*H[l+1,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)-(l-1)^2*H[l-1,alpha,beta,s]^2*H[l,alpha,beta,s]/(4*(l-1/2)^2*l^2)+(l-1)*(7*l-3)*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]^2/(4*(l-1/2)^2*l^2)-(l+2)*(7*l+10)*H[l+1,alpha,beta,s]^2*H[l+2,alpha,beta,s]/(4*(l+1)^2*(l+3/2)^2)+(l+3)*H[l+1,alpha,beta,s]*H[l+2,alpha,beta,s]*H[l+3,alpha,beta,s]/(12*(l+1)*(l+3/2)^2)-(l-2)*H[l-2,alpha,beta,s]*H[l-1,alpha,beta,s]*H[l,alpha,beta,s]/(12*(l-1/2)^2*l))]]]]];
lambda[l_,m_,s_,gamma_]:=coef0[l]-coef1[l,m,s]*gamma+coef2[l,m,s]*gamma^2-coef3[l,m,s]*gamma^3+coef4[l,m,s]*gamma^4-coef5[l,m,s]*gamma^5+coef6[l,m,s]*gamma^6+gamma^2+2*m*gamma-2;
Delta[r_,M_,a_]:=r^2-2*M*r+a^2*M^2;
Ka[r_,M_,a_,omega_,m_]:=(r^2+a^2*M^2)*omega+m*a*M;
rplus[M_,a_]:=M*(1+Sqrt[1-a^2]);
rminus[M_,a_]:=M*(1-Sqrt[1-a^2]);
rstar[r_,M_,a_,omega_,m_]:=If[a== 0,r+2*M*Log[r/(2*M)-1],r+(2*M*rplus[M,a]+a*M*m/omega)/(rplus[M,a]-rminus[M,a])*Log[r/rplus[M,a]-1]-(2*M*rminus[M,a]+a*M*m/omega)/(rplus[M,a]-rminus[M,a])*Log[r/rminus[M,a]-1]];
invertr[rs_,M_,a_,omega_,m_]:=If[a==0,2*M*(1+ProductLog[Exp[rs/(2*M)-1]]),Module[{start=rplus[M,a]+0.00001},Re[Evaluate[u/.FindRoot[rstar[u,M,a,omega,m]==rs,{u,start}]]]]];
alpha2[M_,a_,omega_,m_]:=a^2*M^2+a*M*m/omega;
alpha[M_,a_,omega_,m_]:=If[alpha2[M,a,omega,m]>=0,Sqrt[alpha2[M,a,omega,m]],I*Sqrt[-alpha2[M,a,omega,m]]];
rho2[r_,M_,a_,omega_,m_]:=r^2+a^2*M^2+a*M*m/omega;
invertleft[rs_,M_,a_,omega_,m_,rdiv_]:=Re[Evaluate[u/.FindRoot[rstar[u,M,a,omega,m]==rs,{u,(rplus[M,a]+rdiv)/2}]]];
invertright[rs_,M_,a_,omega_,m_]:=Re[Evaluate[u/.FindRoot[rstar[u,M,a,omega,m]==rs,{u,Sqrt[-alpha2[M,a,omega,m]]*1.01}]]];


(* ::Input::Initialization:: *)
epsilon=1;
V1[r_,M_,a_,omega_,l_,m_]:=If[alpha2[M,a,omega,m]>=0,Delta[r,M,a]/rho2[r,M,a,omega,m]^2*(lambda[l,m,s,omega*M*a]+2-alpha2[M,a,omega,m]*Delta[r,M,a]/rho2[r,M,a,omega,m]^2-epsilon*I*Sqrt[alpha2[M,a,omega,m]]*(2*(r-M)/rho2[r,M,a,omega,m]-4*r*Delta[r,M,a]/rho2[r,M,a,omega,m]^2)),Delta[r,M,a]/rho2[r,M,a,omega,m]^2*(lambda[l,m,s,omega*M*a]+2+Sqrt[-alpha2[M,a,omega,m]]*(Sqrt[-alpha2[M,a,omega,m]]-epsilon*4*r)*Delta[r,M,a]/rho2[r,M,a,omega,m]^2+epsilon*2*Sqrt[-alpha2[M,a,omega,m]]*(r-M)/rho2[r,M,a,omega,m])];


(* ::Input::Initialization:: *)
coord[rs_]:=Evaluate[r/.tortoise[[1]]][rs];
coordl[rs_]:=Evaluate[r/.tortoisel[[1]]][rs];
coordr[rs_]:=Evaluate[r/.tortoiser[[1]]][rs];
func1[rs_]:=Evaluate[Xl/.soll[[1]]][rs];
coef[M_,a_,omega_,m_]:=func1[rstar[rdiv-prox,M,a,omega,m]]/Sqrt[prox];
printsteps=0;
far = 1000.;
solver[M_,a_,omega_,l_,m_]:=If[omega>-a*m/(2*rplus[M,a]),print=0;prox=10^(-5);rminf=rstar[rplus[M,a]+prox,M,a,omega,m];rpinf=rstar[rplus[M,a]+far,M,a,omega,m];If[printsteps==1,Print[1]];tortoise=NDSolve[{r'[rs]==(r[rs]^2+a^2*M^2-2*M*r[rs])/(r[rs]^2+a^2*M^2+a*M*m/omega),{r[rminf]==rplus[M,a]+prox}},r,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->18];If[printsteps==1,Print[2]];NDSolve[{X''[rs]+omega^2*X[rs]==V1[coord[rs],M,a,omega,l,m]*X[rs],{X[rminf]== Ainh*Exp[-I*omega*rminf],X'[rminf] == -I*omega*Ainh*Exp[-I*omega*rminf]}},X,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->10],print=1;rdiv=Sqrt[-alpha2[M,a,omega,m]];If[omega<0.1,If[a<0.9,prox=10^(-3),prox=10^(-2)],prox=10^(-4)];rminf=rstar[rplus[M,a]+prox,M,a,omega,m];rpinf=rstar[rplus[M,a]+far,M,a,omega,m];If[printsteps==1,Print[1]];tortoisel=NDSolve[{r'[rs]==(r[rs]^2+a^2*M^2-2*M*r[rs])/(r[rs]^2+a^2*M^2+a*M*m/omega),{r[rminf]==rplus[M,a]+prox}},r,{rs,rstar[rdiv-prox,M,a,omega,m],rminf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->18];If[printsteps==1,Print[2]];tortoiser=NDSolve[{r'[rs]==(r[rs]^2+a^2*M^2-2*M*r[rs])/(r[rs]^2+a^2*M^2+a*M*m/omega),{r[rstar[rdiv+prox,M,a,omega,m]]==rdiv+prox}},r,{rs,rstar[rdiv+prox,M,a,omega,m],rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->18];If[printsteps==1,Print[3]];soll=NDSolve[{Xl''[rs]+omega^2*Xl[rs]==V1[coordl[rs],M,a,omega,l,m]*Xl[rs],{Xl[rminf]== Ainh*Exp[I*omega*rminf],Xl'[rminf] == I*omega*Ainh*Exp[I*omega*rminf]}},Xl,{rs,rstar[rdiv-prox,M,a,omega,m],rminf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->10];If[printsteps==1,Print[4]];NDSolve[{X''[rs]+omega^2*X[rs]==V1[coordr[rs],M,a,omega,l,m]*X[rs],{X[rstar[rdiv+prox,M,a,omega,m]]==coef[M,a,omega,m]*Sqrt[prox],X'[rstar[rdiv+prox,M,a,omega,m]] == coef[M,a,omega,m]/2*Delta[rdiv+prox,M,a]/rho2[rdiv+prox,M,a,omega,m]/Sqrt[prox]}},X,{rs,rstar[rdiv+prox,M,a,omega,m],rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->10]];


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
file1=OpenWrite["test_1_fM.txt"];
file2=OpenWrite["test_1_gM.txt"];
s0=Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nba}];
Ain0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nba}];
Aout0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nba}];
contrib = Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nba}];
Al0 = Table[Table[0,{j,1,nbener}],{k,1,nba}];
Al1 = Table[Table[0,{j,1,nbener}],{k,1,nba}];
optimize=1;
For[k=1,k<nba+1,k++,For[j=1,j<nbener+1,j++,Al0[[k,j]]=0;Al1[[k,j]]=0;For[i=1,i<nbmodes+1,i++,If[optimize==1 && i>1 && Al0[[k,j]]>10^(-100) && Abs[Max[contrib[[k,j,i-1]]]/Al0[[k,j]]]<0.00001,i=nbmodes+1,For[m=-l[[i]],m<l[[i]]+1,m++,If[a[[k]]!=0 && optimize==1 && m>-l[[i]] && Al0[[k,j]]>10^(-100)  && Abs[contrib[[k,j,i,m+l[[i]]]]/Al0[[k,j]]]<0.00001,m=l[[i]]+1,If[a[[k]]==0 && m>-l[[i]],s0[[k,j,i,m+l[[i]]+1]]=s0[[k,j,i,m+l[[i]]]],s0[[k,j,i,m+l[[i]]+1]]=solver[M,a[[k]],omega[[j]],l[[i]],m]];Module[{rpinf=rstar[rplus[M,a[[k]]]+far,M,a[[k]],omega[[j]],m]},Ain0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]+Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[I*omega[[j]]*rpinf]);Aout0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]-Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[-I*omega[[j]]*rpinf])];contrib[[k,j,i,m+l[[i]]+1]]=(-1)^(print)/Abs[Aout0[[k,j,i,m+l[[i]]+1]]]^2/(Exp[4*Pi*(1+1/Sqrt[1-a[[k]]^2])*omega[[j]]*M+2*Pi*m*a[[k]]/Sqrt[1-a[[k]]^2]]-1);Al0[[k,j]]=Al0[[k,j]]+contrib[[k,j,i,m+l[[i]]+1]];Al1[[k,j]]=Al1[[k,j]]-m*contrib[[k,j,i,m+l[[i]]+1]];If[print==1,Print[{k,j,l[[i]],m,contrib[[k,j,i,m+l[[i]]+1]]}," -> superradiance"],Print[{k,j,l[[i]],m,contrib[[k,j,i,m+l[[i]]+1]]}]]]]]]];WriteLine[file1,ExportString[{Al0[[k]]},"TSV"]];WriteLine[file2,ExportString[{Al1[[k]]},"TSV"]]];
Close[file1];
Close[file2];
Print["DONE"]
