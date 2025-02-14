(* ::Package:: *)

(* ::Input::Initialization:: *)
ClearAll["Global`*"];


(* ::Input::Initialization:: *)
nbmodes = 30;
nbener = 50;
M = 0.5;
s = 1;
nbrQ = 50;
rQmin=0.01;
rQmax=0.999;
rQ=Table[(1-10^(Log10[1-rQmin]+(Log10[1-rQmax]-Log10[1-rQmin])/(nbrQ-1)*i))*M,{i,0,nbrQ-1}]
(*Emin = 10^(-2);
Einter=1;
Emax=5
omega=Table[If[i<100,N[10^(Log10[Emin]+(Log10[Einter]-Log10[Emin])/100*i)],N[Einter+(Emax-Einter)/(99)*(i-100)]],{i,0,nbener-1}]*)
Emin = 10^(-6);
Emax = 10^(-4);
omega = Table[N[10^(Log10[Emin] + (Log10[Emax] - Log10[Emin])/(nbener - 1)*i)],{i,0,nbener-1}]
l=Table[i+1,{i,0,nbmodes-1}]
Ainh = 1;


(* ::Input::Initialization:: *)
nu1[l_]:=l*(l+1);
rplus[M_,rQ_]:=M*(1+Sqrt[1-4*rQ^2/(4*M^2)]);
rminus[M_,rQ_]:=M*(1-Sqrt[1-4*rQ^2/(4*M^2)]);
rstar[r_,M_,rQ_]:=r + rplus[M,rQ]^2/(rplus[M,rQ] - rminus[M,rQ])*Log[r/rplus[M,rQ] - 1] - rminus[M,rQ]^2/(rplus[M,rQ] - rminus[M,rQ])*Log[r/rminus[M,rQ] - 1];
invertr[rs_,M_,rQ_]:=If[rQ==0,2*M*(1+ProductLog[Exp[rs/(2*M)-1]]),Module[{start=If[rs<=5,rplus[M,rQ]+0.00001,rs]},Re[Evaluate[u/.FindRoot[rstar[u,M,rQ]==rs,{u,start}]]]]];


(* ::Input::Initialization:: *)
H[r_]:=r^2;
F[r_,M_,rQ_]:=1 - 2*M/r + rQ^2/r^2;
G[r_,M_,rQ_]:=1 - 2*M/r + rQ^2/r^2;


(* ::Input::Initialization:: *)
V1[r_,M_,rQ_,l_]:=nu1[l]*G[r,M,rQ]/H[r];


(* ::Input::Initialization:: *)
printsteps=0;
prox = 0.000001;
far = 10000000.; (* warning 350. is ok at higher energies*)
coord[rs_]:=Evaluate[r/.tortoise[[1]]][rs];
solver[M_,rQ_,omega_,l_]:=Module[{rminf = rstar[rplus[M,rQ]+prox,M,rQ],
rpinf = rstar[rplus[M,rQ]+far,M,rQ]},tortoise=NDSolve[{r'[rs]==(r[rs] - rplus[M,rQ])*(r[rs]-rminus[M,rQ])/r[rs]^2,{r[rminf]==rplus[M,rQ]+prox}},r,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->18];NDSolve[{X''[rs]+omega^2*X[rs]==V1[coord[rs],M,rQ,l]*X[rs],{X[rminf]== Ainh*Exp[-I*omega*rminf],X'[rminf] == -I*omega*Ainh*Exp[-I*omega*rminf]}},X,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->10]];


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
file1=OpenWrite["fM_charged_1_low.txt"];
s0=Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrQ}];
Ain0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrQ}];
Aout0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrQ}];
contrib = Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrQ}];
Al0 = Table[Table[0,{j,1,nbener}],{k,1,nbrQ}];
optimize=1;
print = 0;
T[M_,rQ_]:=(rplus[M,rQ] - rminus[M,rQ])/(4*Pi*rplus[M,rQ]^2);


(* ::Input::Initialization:: *)
For[k=1,k<nbrQ+1,k++,For[j=1,j<nbener+1,j++,Al0[[k,j]]=0;For[i=1,i<nbmodes+1,i++,If[optimize==1 && i>1 && Al0[[k,j]]>10^(-100) && Abs[Max[contrib[[k,j,i-1]]]/Al0[[k,j]]]<0.00001,i=nbmodes+1,For[m=-l[[i]],m<l[[i]]+1,m++,If[ m>-l[[i]] && optimize==1 && Al0[[k,j]]>10^{-100} && Abs[contrib[[k,j,i,m+l[[i]]]]/Al0[[k,j]]]<0.00001,m=l[[i]]+1,If[m>-l[[i]],s0[[k,j,i,m+l[[i]]+1]]=s0[[k,j,i,m+l[[i]]]],s0[[k,j,i,m+l[[i]]+1]]=solver[M,rQ[[k]],omega[[j]],l[[i]]]];Module[{rpinf=rstar[rplus[M,rQ[[k]]]+far,M,rQ[[k]]]},Ain0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]+Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[I*omega[[j]]*rpinf]);Aout0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]-Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[-I*omega[[j]]*rpinf])];contrib[[k,j,i,m+l[[i]]+1]]=(-1)^print/Abs[Aout0[[k,j,i,m+l[[i]]+1]]]^2/(Exp[omega[[j]]/T[M,rQ[[k]]]]-1);Al0[[k,j]]=Al0[[k,j]]+contrib[[k,j,i,m+l[[i]]+1]];If[m == -l[[i]],Print[{k,j,l[[i]],m,contrib[[k,j,i,m+l[[i]]+1]]}]]]]]]];WriteLine[file1,ExportString[{Al0[[k]]},"TSV"]]];
Close[file1];
Print["DONE"]



