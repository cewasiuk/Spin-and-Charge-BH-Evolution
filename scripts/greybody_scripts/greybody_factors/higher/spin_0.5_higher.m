(* ::Package:: *)

(* export PATH=$PATH:/cygdrive/c/Program\ Files/Wolfram\ Research/Mathematica/11.1/*)
ClearAll["Global`*"];


(* ::Input::Initialization:: *)
nbmodes = 30;
nbener = 50;
M = 0.5;
s = 0;
nbrn = 7;
n = {0,1,2,3,4,5,6}
Mstar =1;
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
l=Table[i+0.5,{i,0,nbmodes-1}]
Ainh = 1;


(* ::Input::Initialization:: *)
nu12[l_]:=l*(l+1)+1/4;
rH[M_,Mstar_,n_]:=1/(Sqrt[Pi]*Mstar)*(M/Mstar)^(1/(n+1))*(8*Gamma[(n+3)/2]/(n+2))^(1/(n+1));
phip=(1+Sqrt[5])/2;
phim=(1-Sqrt[5])/2;
rstar[r_,M_,Mstar_,n_]:=Module[{result=X[u]/.DSolve[X'[u]==1/(1-(rHH/u)^(n+1)),X[u],u][[1,1]]},If[n==0,result/.{C[1]->-rHH*Log[rHH]},If[n==1,Re[result/.{C[1]->I*Pi*rHH/2}],If[n==2,Re[result/.{C[1]->-rHH/3*Log[-rHH]+rHH/6*Log[rHH^2]+Pi/(2*Sqrt[3])*rHH}],If[n==3,Re[result/.{C[1]->-rHH/4*Log[-1]+rHH*Pi/4}],If[n==4,Re[result/.{C[1]->Sqrt[Sqrt[5]]/5*(Sqrt[phip]+Sqrt[-phim])*rHH*Pi/2-rHH/5*Log[-rHH]+rHH/10*Log[rHH^2]}],If[n==5,Re[result/.{C[1]->-1/6*rHH*Log[-1]+Pi/(2*Sqrt[3])*rHH}],If[n==6,Re[result/.{C[1]->-rHH/7*Log[-1]+rHH*Pi/7*(Cos[Pi/14]+Cos[3*Pi/14]+Sin[Pi/7])}],0]]]]]]]]/.{u->r,rHH->rH[M,Mstar,n]};
(*Simplify[X[u]/.DSolve[X'[u]\[Equal]u^2/(u-rp)/(u-rm),X[u],u][[1,1]]]
rH[M,Mstar,n]
rstar[2,M,Mstar,n]
(r+2*rHH/7*(ArcTan[Cos[Pi/14]/(r/rHH+Sin[Pi/14])]*Cos[Pi/14]+ArcTan[Cos[3*Pi/14]/(r/rHH-Sin[3*Pi/14])]*Cos[3*Pi/14]+ArcTan[Sin[Pi/7]/(r/rHH+Cos[Pi/7])]*Sin[Pi/7])+rHH/7*(Log[r/rHH-1]+Sin[3*Pi/14]*Log[(r/rHH)^2-2*r/rHH*Sin[3*Pi/14]+1]-Sin[Pi/14]*Log[(r/rHH)^2+2*r/rHH*Sin[Pi/14]+1]-Cos[Pi/7]*Log[(r/rHH)^2+2*r/rHH*Cos[Pi/7]+1]))/.{rHH->rH[M,Mstar,n],r->2}*)


(* ::Input::Initialization:: *)
H[r_]:=r^2;
F[r_,M_,Mstar_,n_]:=1-(rH[M,Mstar,n]/r)^(n+1);
G[r_,M_,Mstar_,n_]:=1-(rH[M,Mstar,n]/r)^(n+1);


(* ::Input::Initialization:: *)
epsilon = 1;
V12[r_,M_,Mstar_,n_,l_]:=nu12[l]*G[r,M,Mstar,n]/H[r] + epsilon*Sqrt[nu12[l]*F[r,M,Mstar,n]*G[r,M,Mstar,n]]*D[Sqrt[G[u,M,Mstar,n]/H[u]],u]/.{u->r};


(* ::Input::Initialization:: *)
printsteps=0;
prox = 0.0000001;
coord[rs_]:=Evaluate[r/.tortoise[[1]]][rs];
solver[M_,Mstar_,n_,omega_,l_]:=Module[{rminf = rstar[rH[M,Mstar,n]+prox,M,Mstar,n],
rpinf = rstar[rH[M,Mstar,n]+350,M,Mstar,n]},tortoise=NDSolve[{r'[rs]==1-(rH[M,Mstar,n]/r[rs])^(n+1),{r[rminf]==rH[M,Mstar,n]+prox}},r,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->18];NDSolve[{X''[rs]+omega^2*X[rs]==V12[coord[rs],M,Mstar,n,l]*X[rs],{X[rminf]== Ainh*Exp[-I*omega*rminf],X'[rminf] == -I*omega*Ainh*Exp[-I*omega*rminf]}},X,{rs,rminf,rpinf},Method->"StiffnessSwitching",AccuracyGoal->18,PrecisionGoal->10]];


(* ::Input::Initialization:: *)
SetDirectory[NotebookDirectory[]];
file1=OpenWrite["fM_higher_0.5_low.txt"];
s0=Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrn}];
Ain0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrn}];
Aout0 = Table[Table[Table[Table[0,{m,1,2*l[[i]]+1}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrn}];
contrib = Table[Table[Table[Table[0,{m,-l[[i]],l[[i]]}],{i,1,nbmodes}],{j,1,nbener}],{k,1,nbrn}];
Al0 = Table[Table[0,{j,1,nbener}],{k,1,nbrn}];
optimize=1;
print = 0;
T[M_,Mstar_,n_]:=(n+1)/(4*Pi*rH[M,Mstar,n]);


Close["fM_higher_0.5_low.txt"]


(* ::Input::Initialization:: *)
For[k=1,k<nbrn+1,k++,For[j=1,j<nbener+1,j++,Al0[[k,j]]=0;For[i=1,i<nbmodes+1,i++,If[optimize==1 && i>1 && Al0[[k,j]]>10^(-100) && Abs[Max[contrib[[k,j,i-1]]]/Al0[[k,j]]]<0.00001,i=nbmodes+1,For[m=-l[[i]],m<l[[i]]+1,m++,If[ m>-l[[i]] && optimize==1 && Al0[[k,j]]>10^{-100} && Abs[contrib[[k,j,i,m+l[[i]]]]/Al0[[k,j]]]<0.00001,m=l[[i]]+1,If[m>-l[[i]],s0[[k,j,i,m+l[[i]]+1]]=s0[[k,j,i,m+l[[i]]]],s0[[k,j,i,m+l[[i]]+1]]=solver[M,Mstar,n[[k]],omega[[j]],l[[i]]]];Module[{rpinf=rstar[rH[M,Mstar,n[[k]]]+350,M,Mstar,n[[k]]]},Ain0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]+Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[I*omega[[j]]*rpinf]);Aout0[[k,j,i,m+l[[i]]+1]]=(Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]][rpinf]-Evaluate[X/.s0[[k,j,i,m+l[[i]]+1]]][[1]]'[rpinf]/(I*omega[[j]]))/(2*Exp[-I*omega[[j]]*rpinf])];contrib[[k,j,i,m+l[[i]]+1]]=(-1)^print/Abs[Aout0[[k,j,i,m+l[[i]]+1]]]^2/(Exp[omega[[j]]/T[M,Mstar,n[[k]]]]+1);Al0[[k,j]]=Al0[[k,j]]+contrib[[k,j,i,m+l[[i]]+1]];If[m ==-l[[i]],Print[{k,j,l[[i]],m,contrib[[k,j,i,m+l[[i]]+1]]}]]]]]]];WriteLine[file1,ExportString[{Al0[[k]]},"TSV"]]];
Close[file1];
Print["DONE"]






