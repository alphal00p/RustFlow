$HistoryLength=0;
(* RustFlow mathematical-data extraction only; no upstream solver code. *)
args = If[Length[$ScriptCommandLine] > 0, Rest[$ScriptCommandLine],
  scriptPosition = FirstPosition[$CommandLine, "-script", Missing["none"], {1}];
  If[MissingQ[scriptPosition], {}, Take[Drop[$CommandLine, First[scriptPosition]+1], UpTo[3]]]];
If[Length[args] != 3, Print["Expected FFlist.txt, primary form_factors_tensor_structures directory, output JSON"]; Exit[1]];
{source, primaryDirectory, outputPath} = args;
If[FileHash[source,"SHA256","HexString"] != "9792a04749d3756d59c8603acc88e851bdc4175bac9941939eb5c5fc594a17f3", Print["Wrong pinned FF input"]; Exit[1]];
fail[x_]:=(Print["FAIL ",x];Exit[1]);
str[x_]:=StringReplace[ToString[x,InputForm]," "->""];
rust[x_Integer|x_Rational]:=str[x];
rust[x_Complex]:="("<>rust[Re[x]]<>"+("<>rust[Im[x]]<>")*imaginary_unit)";
rust[x_Symbol]:=str[x];
rust[x_Plus]:="("<>StringRiffle[rust/@(List@@x),"+"]<>")";
rust[x_Times]:="("<>StringRiffle[rust/@(List@@x),"*"]<>")";
rust[x_Power]/;IntegerQ[x[[2]]]:="("<>rust[x[[1]]]<>")^("<>rust[x[[2]]]<>")";
rust[x_]:=fail[{"nonrational serialized coefficient",x}];
ff=Get[source]/.{vv[1]->s,vv[2]->t,vv[3]->b};
radicands=DeleteDuplicates[Cases[ff,Power[x_,q_Rational]/;Denominator[q]>1:>Expand[x],Infinity]];
If[!And@@(PolynomialQ[#,ieps0]&&Exponent[#,ieps0]<=1& /@ radicands),fail["nonlinear i0"]];
rootSymbols=Table[Symbol["ffroot"<>ToString[i]],{i,Length[radicands]}];
replacement=Thread[rootSymbols->(Sqrt/@radicands)];
perms={{s,b-s-t,b},{s,t,b},{b-s-t,t,b},{t,s,b}};
rows={}; primaryHashes={};
Do[
 primaryFile=primaryDirectory<>"/form_factor_"<>ToString[j]<>".m";
 primary=Get[primaryFile]/.myRoot[x_,1/2]:>Sqrt[x]/.{s12->s,mmH->b};
 If[Expand[primary-(ff[[j]]/.ieps0->0)]=!=0,fail[{"primary mismatch",j}]];
 AppendTo[primaryHashes,FileHash[primaryFile,"SHA256","HexString"]];
 terms=List@@Expand[ff[[j]]];
 pairs=Table[masters=Cases[term,_canMIPlanar|_canMINP,{0,Infinity}];
  If[Length[masters]!=1,fail["nonlinear master"]];mi=First[masters];coefficient=term/.mi->1;
  If[Expand[term-coefficient mi]=!=0,fail["term reconstruction"]];{mi,coefficient},{term,terms}];
 grouped=Map[{#[[1,1]],Total[#[[All,2]]]}&,GatherBy[pairs,First]];
 row=Table[
  mi=pair[[1]];raw=pair[[2]];
  coefficient=raw/.Power[x_,q_Rational]/;Denominator[q]>1:>Module[{where},If[Denominator[q]!=2,fail["non-square root"]];where=First@FirstPosition[radicands,Expand[x],Missing["none"],{1},Heads->False];rootSymbols[[where]]^(2q)];
  coefficient=Together[coefficient];
  If[!FreeQ[coefficient,ieps0],fail["unremoved i0"]];
  If[Expand[(coefficient/.replacement)-raw]=!=0,If[Together[(coefficient/.replacement)-raw]=!=0,fail[{"root reconstruction",j,mi}]]];
  position=FirstPosition[perms,mi[[3]],Missing["none"],{1},Heads->False];
  If[MissingQ[position],fail[{"crossing",mi}]];
  <|"family"->If[Head[mi]===canMIPlanar,"Planar_EW1","NP_EW1"],"permutation"->First[position],"index"->mi[[1]]-1,"epsilon"->mi[[2]],"coefficient"->rust[coefficient]|>,{pair,grouped}];
 AppendTo[rows,row]; Print["row ",j," terms=",Length[terms]," grouped=",Length[row]," primary+linear+root reconstruction passed"];
,{j,4}];
result=<|"schema"->"rustflow-gg-hg-generic-form-factors-v1","source"-><|"url"->"https://bitbucket.org/aschweitzer/mg5_higgs_ew_plugin/src/89f64b93d0bdbd8ee90eb85737021189229b034a/ComputationFormFacGGHGEW/mathematicaRoutines/FFlist.txt","sha256"->FileHash[source,"SHA256","HexString"]|>,"paper"->"https://arxiv.org/abs/2112.07578v1","primary_form_factor_sha256"->primaryHashes,"coordinates"->{"s","t","b"},"normalization"->"-1/(mV2^2*(4*pi)^4)","crossings"->Map[str,perms,{2}],"roots"->MapThread[<|"name"->str[#1],"radicand"->rust[#2/.ieps0->0],"i0_slope"->rust[Coefficient[#2,ieps0]]|>&,{rootSymbols,radicands}],"rows"->rows,"evidence"-><|"exact_primary_equality_at_zero_i0"->True,"exact_linear_reconstruction"->True,"exact_root_reconstruction"->True,"original_ieps0_convention"->"ieps0 is positive imaginary; native coefficients use exact side limits"|>|>;
Export[outputPath,result,"RawJSON"];
Print["DONE roots=",Length[radicands]," rows=",Length/@rows];Exit[0];
