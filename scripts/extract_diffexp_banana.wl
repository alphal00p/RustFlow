(* Scientific matrix and asymptotic-boundary extraction only. No DiffExp solver is run. *)
$HistoryLength=0;
fail[m_]:=(Print["EXTRACTION FAILURE: ",m];Exit[1]);
requiredEnvironment[name_]:=Module[{value=Environment[name]},
 If[!StringQ[value]||StringLength[value]==0,fail["Set "<>name<>" before running"]];
 ExpandFileName[value]];
base=requiredEnvironment["DIFFEXP_OUTPUT"];
root=requiredEnvironment["DIFFEXP_ROOT"];
source=FileNameJoin[{root,"Examples","Banana_EqualMass_Matrices"}];
If[FileExistsQ[base]||DirectoryQ[base],fail["DIFFEXP_OUTPUT must be a new directory"]];
If[CreateDirectory[base,CreateIntermediateDirectories->True]===$Failed,fail["Cannot create output directory"]];
sha[f_]:=IntegerString[FileHash[f,"SHA256"],16,64];
str[v_]:=ToString[v,InputForm];

matrixHashes={"d87a9039ab70226a0afbb4c4c2c37f2fbfde6c638c30a2771087734ae982ef1b","c4c597adeda5872d532d3edb82a11bdfbd04456f39e6bd1d5fd538390f96fb71","5a1d902c63c7a08a5c26afb4004c2f3b416f79d1d6166443f6e2d513c58d62ab","5a1d902c63c7a08a5c26afb4004c2f3b416f79d1d6166443f6e2d513c58d62ab","5a1d902c63c7a08a5c26afb4004c2f3b416f79d1d6166443f6e2d513c58d62ab"};
Do[file=FileNameJoin[{source,"dt_"<>ToString[k]<>".m"}];
 If[!FileExistsQ[file]||sha[file]=!=matrixHashes[[k+1]],fail["Pinned matrix missing or hash mismatch"]],{k,0,4}];
matrices=Table[Get[FileNameJoin[{source,"dt_"<>ToString[k]<>".m"}]],{k,0,4}];
If[Dimensions[matrices]=!={5,4,4}||Flatten[matrices[[3;;5]]]=!=ConstantArray[0,48],fail["Matrix shape or higher epsilon coefficients"]];
m=matrices[[1]]+eps*matrices[[2]];
mx=Together[(m/.t->-1/x)/x^2];
residue=Map[Limit[x*#,x->0]&,mx,{2}];
(* Gamma(1+z) expansion cancels EulerGamma exactly. It is used here only
   to export exact boundary expressions, not numerical starting values. *)
c[j_,n_]:=Switch[j,0,3*(-1)^n,1,2+3*(-1)^n-2^n,
 2,3+(-1)^n+(-2)^n-3^n,3,4+(-3)^n-4^n];
a={4,-12,12,-4};
g[j_]:=Exp[j*eps*logx+Sum[c[j,n]*Zeta[n]*eps^n/n,{n,2,6}]];
terms=Table[Normal[Series[g[j],{eps,0,6}]],{j,0,3}];
b3=Expand[(1+3*eps)*(1+4*eps)*Sum[a[[j+1]]*terms[[j+1]],{j,0,3}]/eps^2];
b2=Expand[(1+3*eps)*Sum[(3-j)*a[[j+1]]*terms[[j+1]],{j,0,3}]/(4*eps)];
b1=Expand[Sum[(2-j)*(3-j)*a[[j+1]]*terms[[j+1]],{j,0,3}]/12];
b4=Normal[Series[Exp[Sum[3*(-1)^n*Zeta[n]*eps^n/n,{n,2,4}]],{eps,0,4}]];
If[!And@@Flatten[Table[Coefficient[b,eps,k]===0,{b,{b1,b2,b3}},{k,-2,0}]],fail["Expected cancellation of negative and zero epsilon coefficients"]];
leading=Table[Table[Expand[Coefficient[b,eps,k]],{k,0,4}],{b,{b1,b2,b3,b4}}];
(* Independent verification of the asymptotic partial-boundary completion
   identities at fixed generic epsilon, sector by sector. *)
Do[
 f=x^(1+j*eps);
 derivedB2=Together[(-x*D[f,x]+(1+3*eps)*f)/(4*(1+4*eps))];
 If[Together[derivedB2/f-eps*(3-j)/(4*(1+4*eps))]=!=0,fail["B2 derivative identity"]],
 {j,0,3}];
records=Table[<|"epsilon_order"->k,"matrix"->Map[str,matrices[[k+1]],{2}],
 "file"->"Examples/Banana_EqualMass_Matrices/dt_"<>ToString[k]<>".m",
 "sha256"->sha[FileNameJoin[{source,"dt_"<>ToString[k]<>".m"}]]|>,{k,0,4}];
exported=Export[FileNameJoin[{base,"matrix-boundary.json"}],<|
 "schema_version"->1,"case"->"diffexp_equal_mass_three_loop_banana",
 "upstream_commit"->"784c8229bf92369a03f011a48e161522c8c54bbd",
 "notebook_sha256"->sha[FileNameJoin[{root,"Examples","Banana.nb"}]],
 "paper_sha256"->sha[FileNameJoin[{root,"Paper","main.tex"}]],
 "dimension"->4,"epsilon_variable"->"eps","dimension_convention"->"D=2-2*eps",
 "physical_variable"->"t=p^2/m^2","matrices"->records,
 "finite_singular_points"->{"0","4","16"},
 "infinity_parameterization"->"t=-1/x, x->0+",
 "infinity_transformed_matrix"->Map[str,mx,{2}],
 "infinity_residue"->Map[str,residue,{2}],
 "infinity_eigenvalues"->(str/@Eigenvalues[residue]),
 "original_partial_boundary_mask"->{False,False,True,True},
 "original_B3_leading_boundary"->"eps*(1+3*eps)*(1+4*eps)*exp(3*EulerGamma*eps)*(4*x*Gamma(eps)^3 + 6*x^(1+eps)*eps*Gamma(-eps)^2*Gamma(eps)^3/Gamma(-2*eps) + 8*x^(1+2*eps)*eps*Gamma(-eps)^3*Gamma(eps)*Gamma(2*eps)/Gamma(-3*eps) + 3*x^(1+3*eps)*eps*Gamma(-eps)^4*Gamma(3*eps)/Gamma(-4*eps))",
 "original_B4_exact_boundary"->"exp(3*EulerGamma*eps)*eps^3*Gamma(eps)^3",
 "exact_completion_identities"-><|
 "B2"->"(t*d(B3)/dt+(1+3*eps)*B3)/(4*(1+4*eps))",
 "B1"->"(t*(t-4)*d(B2)/dt+(4+t+2*eps*(8+t))*B2-(1+3*eps)*B3)/(3*t*(1+3*eps))"|>,
 "leading_epsilon_coefficients_by_master_wolfram"->Map[str,leading,{2}],
 "leading_coefficient_multipliers"->{"x","x","x","1"},
 "log_symbol"->"logx=log(x), real on x>0",
 "boundary_scope"->"B1/B2 entries were independently deduced from exact DE identities; only B3 and B4 were supplied in the notebook. B1/B2/B3 formulas are leading x powers with logarithms, not complete regular-point boundary values. Higher x orders require the differential-equation recurrence.",
 "notebook_settings"-><|"epsilon_order"->4,"working_precision"->500,"chop_precision"->250,
 "expansion_order"->50,"division_order"->3,"mobius"->True,"pade"->True,"verbosity"->2,"parallel"->False,"prescription"->"t-16+I*delta"|>,
 "notebook_routes"->{"infinity t=-1/x, x=0+ -> t=-1", "t=-1 -> t=32; save local series; evaluate saved Pade representations at t=20"},
 "no_numerical_transport_run"->True|>,"RawJSON"];
If[exported===$Failed,fail["Failed matrix/boundary export"]];
Put[leading,FileNameJoin[{base,"leading-boundary-coefficients.wl"}]];
Print["EXTRACTION PASS residue eigenvalues ",Eigenvalues[residue]," leading B3 eps1 ",leading[[3,2]]];Exit[0];
