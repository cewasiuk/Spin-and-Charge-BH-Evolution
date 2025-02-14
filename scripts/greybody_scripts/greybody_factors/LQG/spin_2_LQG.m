(* ::Package:: *)

(* ::Input::Initialization:: *)
(*export PATH=$PATH:/cygdrive/c/Program\ Files/Wolfram\ Research/Mathematica/11.1/*)
ClearAll["Global`*"];


(* ::Input::Initialization:: *)
nbmodes = 30;
nbener = 50;
M = 0.5;
s = 2;
a0 = 0.11;
nbre = 10;
epsilon = {10^(-1),10^(-0.9),10^(-0.8),10^(-0.7),10^(-0.6),10^(-0.5),10^(-0.4),10^(-0.3),10^(-0.2),10^(-0.1)}
(*omega = Table[10^(Log10[0.01] + (Log10[5] - Log10[0.01])/(nbener-1)*i),{i,0,nbener-1}]*)
(*Emin = 10^(-2);
Einter = 1;
Emax = 5;
omega = Table[
  If[i < 100, 
   N[10^(Log10[Emin] + (Log10[Einter] - Log10[Emin])/100*i)], 
   N[Einter + (Emax - Einter)/(99)*(i - 100)]], {i, 0, nbener - 1}]*)
Emin = 10^(-6);
Emax = 10^(-4);
omega = Table[N[10^(Log10[Emin] + (Log10[Emax] - Log10[Emin])/(nbener - 1)*i)],{i,0,nbener-1}]
l=Table[i+2,{i,0,nbmodes-1}]
Ainh = 1;


(* ::Input::Initialization:: *)
nu2[l_]:=l*(l+1)-2;
P[epsilon_]:=(Sqrt[1+epsilon^2]-1)/(Sqrt[1+epsilon^2]+1);
rplus[M_,epsilon_]:=2*M/(1+P[epsilon])^2;
rminus[M_,epsilon_]:=rplus[M,epsilon]*P[epsilon]^2;
(*rstar[r_,M_,epsilon_,a0_]:=(r-a0^2/(rm*rp*r)+a0^2*(rp+rm)/(rp^2*rm^2)*Log[r/(rm+rp)]+(a0^2+rp^4)/(rp^2*(rp-rm))*Log[r/rp-1]+(a0^2+rm^4)/(rm^2*(rm-rp))*Log[r/rm-1])/.{rp->rplus[M,epsilon],rm->rminus[M,epsilon]};*)
rstar[r_,M_,epsilon_,a0_]:=(r-a0^2/(r*rm*rp)+a0^2*(rp+rm)/(rm^2*rp^2)*Log[r]+(a0^2+rm^4)/(rm^2*(rm-rp))*Log[r-rm]+(a0^2+rp^4)/(rp^2*(rp-rm))*Log[r-rp])/.{rp->rplus[M,epsilon],rm->rminus[M,epsilon]};
rplus[M,0.1]
rminus[M,0.1]
prox = 0.0000001
rstar[rplus[M,0.1]+prox,M,0.1,a0]


(* ::Input::Initialization:: *)
H[r_,a0_]:=r^2+a0^2/r^2;
F[r_,M_,epsilon_,a0_]:=(r-rplus[M,epsilon])*(r-rminus[M,epsilon])*r^4/(r+Sqrt[rplus[M,epsilon]*rminus[M,epsilon]])^2/(r^4+a0^2);
G[r_,M_,epsilon_,a0_]:=(r-rplus[M,epsilon])*(r-rminus[M,epsilon])*(r+Sqrt[rplus[M,epsilon]*rminus[M,epsilon]])^2/(r^4+a0^2);


(* ::Input::Initialization:: *)
V2[r_,M_,epsilon_,a0_,l_]:=nu2[l]*G[r,M,epsilon,a0]/H[r,a0] + F[r,M,epsilon,a0]*G[r,M,epsilon,a0]/(2*H[r,a0]^2)*(D[H[u,a0],u]/.{u->r})^2 - 1/2*Sqrt[F[r,M,epsilon,a0]*G[r,M,epsilon,a0]/H[r,a0]]*D[Sqrt[F[u,M,epsilon,a0]*G[u,M,epsilon,a0]/H[u,a0]]*D[H[u,a0],u],u]/.{u->r};


(* ::Input::Initialization:: *)
printsteps=0;
prox = 0.000000001;
far = 10000000.; (* warning 350. is ok at higher energies*)
coord[rs_]:=Evaluate[r/.tortoise[[1]]][rs];
solver[M_,epsilon_,a0_,omega_,l_]:=Module[{rminf = rstar[rplus[M,epsilon]+prox,M,epsilon,a0],
rpinf = rstar[rplus[M,epsilon]+far,M,epsilon,a0]},tortoise=NDSolve[{r'[rs]==r[rs]^2*(r[rs]-rplus[M,epsilon])*(r[rs]-rminus[M,epsilon])/(r[rs]^4+a0^2),{r[rminf]==rplus[M,epsilon]+prox}},r,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->18];NDSolve[{X''[rs]+omega^2*X[rs]==V2[coord[rs],M,epsilon,a0,l]*X[rs],{X[rminf]== Ainh*Exp[-I*omega*rminf],X'[rminf] == -I*omega*Ainh*Exp[-I*omega*rminf]}},X,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->10]];


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
file1=OpenWrite["fM_LQG_2_low_a0.txt"];
s0=Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbre}];
Ain0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbre}];
Aout0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbre}];
contrib = Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbre}];
Al0 = Table[Table[0,{j,1,nbener}],{k,1,nbre}];
optimize=1;
print = 0;
mLQG[M_,epsilon_]:=M/(1+P[epsilon])^2;
T[M_,epsilon_,a0_]:=4*mLQG[M,epsilon]^3*(1-P[epsilon]^2)/(32*Pi*mLQG[M,epsilon]^4+2*Pi*a0^2);


(* ::Input::Initialization:: *)
For[k=1,k<nbre+1,k++,For[j=1,j<nbener+1,j++,Al0[[k,j]]=0;For[i=1,i<nbmodes+1,i++,If[optimize==1 && i>1 && Al0[[k,j]]>10^(-100) && Abs[Max[contrib[[k,j,i-1]]]/Al0[[k,j]]]<0.00001,i=nbmodes+1,For[m=-l[[i]],m<l[[i]]+1,m++,If[ m>-l[[i]] && optimize==1 && Al0[[k,j]]>10^{-100} && Abs[contrib[[k,j,i,m+l[[i]]]]/Al0[[k,j]]]<0.00001,m=l[[i]]+1,If[m>-l[[i]],s0[[k,j,i,m+l[[i]]+1]]=s0[[k,j,i,m+l[[i]]]],s0[[k,j,i,m+l[[i]]+1]]=solver[M,epsilon[[k]],a0,omega[[j]],l[[i]]]];Module[{rpinf=rstar[rplus[M,epsilon[[k]]]+far,M,epsilon[[k]],a0]},Ain0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]+Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[I*omega[[j]]*rpinf]);Aout0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]-Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[-I*omega[[j]]*rpinf])];contrib[[k,j,i,m+l[[i]]+1]]=(-1)^print/Abs[Aout0[[k,j,i,m+l[[i]]+1]]]^2/(Exp[omega[[j]]/T[M,epsilon[[k]],a0]]-1);Al0[[k,j]]=Al0[[k,j]]+contrib[[k,j,i,m+l[[i]]+1]];If[m ==-l[[i]],Print[{k,j,l[[i]],m,contrib[[k,j,i,m+l[[i]]+1]]}]]]]]]];WriteLine[file1,ExportString[{Al0[[k]]},"TSV"]]];
Close[file1];
Print["DONE"]



