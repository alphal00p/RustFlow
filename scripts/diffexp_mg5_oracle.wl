(* Original-DiffExp oracle on the closed two-master NP_EW1 scientific-data block.
   This driver is independently written; it loads DiffExp unchanged and does not
   reproduce the MG5 application's solver implementation.

   Required environment:
     DIFFEXP_ROOT: pinned DiffExp checkout containing DiffExp.m.
     MG5_GRID_ROOT: extracted grid root containing NP_EW1/atilde_NP_EW1.txt
       and NP_EW1/BClist/BClist_NP_EW1_p1011194946161727751.txt.
     DIFFEXP_OUTPUT: a new directory for matrices, raw results and metadata.

   Example (set the three directory variables before running):
     OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
     timeout --signal=TERM --kill-after=10s 120s \
     bern-wolfram -noinit -noprompt -script scripts/diffexp_mg5_oracle.wl

   Exact inputs and inherited error estimates are preserved. Claims are capped
   at 24 digits from the supplied grid error, not the much longer mantissas.
   This is a closed subsystem transport check, not full amplitude acceptance. *)
$HistoryLength=0;
driverFile=ExpandFileName[$InputFileName];
fail[m_]:=(Print["ORACLE FAILURE: ",m];Exit[1]);
requiredEnvironment[name_]:=Module[{value=Environment[name]},
 If[!StringQ[value]||StringLength[value]==0,fail["Set "<>name<>" before running this script"]];
 ExpandFileName[value]];
base=requiredEnvironment["DIFFEXP_OUTPUT"];
source=FileNameJoin[{requiredEnvironment["MG5_GRID_ROOT"],"NP_EW1"}];
package=FileNameJoin[{requiredEnvironment["DIFFEXP_ROOT"],"DiffExp.m"}];
sha[f_]:=IntegerString[FileHash[f,"SHA256"],16,64];
str[v_]:=ToString[v,InputForm];
scalar[v_]:=<|"real_wolfram"->str[Re[v]],"imaginary_wolfram"->str[Im[v]],"recorded_precision"->str[Precision[v]]|>;
If[!FileExistsQ[package]||sha[package]=!="67fa23a0be32f747e292debcf3142b0ecef2cf526335b5b31f5ec7dd9cdeace7",fail["Pinned package missing or hash mismatch"]];
matrixFile=FileNameJoin[{source,"atilde_NP_EW1.txt"}];
boundaryFile=FileNameJoin[{source,"BClist","BClist_NP_EW1_p1011194946161727751.txt"}];
If[!FileExistsQ[matrixFile]||sha[matrixFile]=!="becb0221008715452488a4b73a64be95ad869dc8fb54657598fe34a1f28921d9",fail["Pinned matrix missing or hash mismatch"]];
If[!FileExistsQ[boundaryFile]||sha[boundaryFile]=!="22b140927743133fc56dc7b0e6c6ab2d3e529d829f5d942556d7357e249d53ea",fail["Pinned boundary missing or hash mismatch"]];
If[FileExistsQ[base]||DirectoryQ[base],fail["DIFFEXP_OUTPUT must be a new directory"]];
If[CreateDirectory[FileNameJoin[{base,"matrices"}],CreateIntermediateDirectories->True]===$Failed,fail["Cannot create output directory"]];
a=Get[matrixFile];
If[Dimensions[a]=!={61,61}||Flatten[a[[{1,2},3;;61]]]=!=ConstantArray[0,118],fail["Non-closed block"]];
connections=Map[Together[D[a[[{1,2},{1,2}]],#]]&,{s,t,b}];
Do[Put[If[k==1,connections[[j]],ConstantArray[0,{2,2}]],FileNameJoin[{base,"matrices","d"<>ToString[{s,t,b}[[j]]]<>"_"<>ToString[k]<>".m"}]],{j,3},{k,0,4}];
rawbc=Get[boundaryFile];
p0=Thread[{s,t,b}->({vv[1],vv[2],vv[3]}/.rawbc[[1]])];
s0=s/.p0;p1=p0/.Rule[s,_]->Rule[s,s0+1];p2=p0/.Rule[s,_]->Rule[s,s0+2];
values=Transpose[rawbc[[2,All,{1,2}]]];
If[Dimensions[values]=!={2,5}||!AllTrue[Flatten[values],NumberQ],fail["Boundary shape"]];
(* Round the supplied global source-error estimate upward; do not treat the
   200+ digit mantissas as measured accuracy. Reserve three digits and cap
   input/output comparison claims at 24 digits. *)
inheritedError=10^-27;
If[!TrueQ[rawbc[[3]]<inheritedError],fail["Source error estimate exceeds conservative envelope"]];
Get[package];records={};
CheckAbort[
Do[
 wp=cfg[[1]];order=cfg[[2]];
 DiffExp`LoadConfiguration[{DiffExp`MatrixDirectory->FileNameJoin[{base,"matrices"}],
 DiffExp`EpsilonOrder->4,WorkingPrecision->wp,DiffExp`ChopPrecision->wp-20,
 DiffExp`ExpansionOrder->order,DiffExp`UseMobius->False,DiffExp`UsePade->False,
 DiffExp`DivisionOrder->3,DiffExp`Verbosity->1,AccuracyGoal->"?",
 DiffExp`DeltaPrescriptions->{s-I*\[Delta],s-1-I*\[Delta]},"Parallel"->False}];
 bc=Append[DiffExp`PrepareBoundaryConditions[N[values,wp],p0],ConstantArray[N[inheritedError,wp],{2,5}]];
 first=AbsoluteTiming[DiffExp`TransportTo[bc,p1]];
 staged=AbsoluteTiming[DiffExp`TransportTo[first[[2]],p2]];
 direct=AbsoluteTiming[DiffExp`TransportTo[bc,p2]];
 Do[If[Dimensions[r[[2]]]=!={2,5}||!AllTrue[Flatten[r[[2]]],NumberQ],fail["Malformed transport result"]],{r,{first[[2]],staged[[2]],direct[[2]]}}];
 (* Check exact returned physical coordinates and an independent analytic
    first-epsilon increment, not only equality between two solver routes. *)
 endpointChecks=Table[
   r=spec[[1]];point=spec[[2]];step=spec[[3]];
   If[!TrueQ[({s,t,b}/.Normal[r[[1]]])==({s,t,b}/.point)],fail["Returned endpoint mismatch"]];
   analytic=values[[All,2]]+{1,2}*Log[(s0+step-1)/(s0-1)];
   error=Max[Abs[r[[2,All,2]]-analytic]];
   If[!TrueQ[error<10^-24],fail[{"First epsilon analytic increment",error}]];
   str[error],
   {spec,{{first[[2]],p1,1},{staged[[2]],p2,2},{direct[[2]],p2,2}}}];
 routeError=Max[Abs[Flatten[staged[[2,2]]-direct[[2,2]]]]];
 If[!TrueQ[routeError<10^-24],fail[{"Direct/staged mismatch",routeError}]];
 Do[Put[r[[2]],FileNameJoin[{base,r[[1]]<>"-order"<>ToString[order]<>".wl"}]],{r,{{"checkpoint",first[[2]]},{"staged",staged[[2]]},{"direct",direct[[2]]}}}];
 AppendTo[records,<|"working_precision"->wp,"expansion_order"->order,"chop_precision"->wp-20,
 "checkpoint"-><|"wall_seconds"->first[[1]],"values"->Map[scalar,first[[2,2]],{2}],"reported_errors"->Map[str,first[[2,3]],{2}]|>,
 "staged"-><|"wall_seconds"->staged[[1]],"values"->Map[scalar,staged[[2,2]],{2}],"reported_errors"->Map[str,staged[[2,3]],{2}]|>,
 "direct"-><|"wall_seconds"->direct[[1]],"values"->Map[scalar,direct[[2,2]],{2}],"reported_errors"->Map[str,direct[[2,3]],{2}]|>,
 "direct_staged_max_absolute_difference"->str[routeError],
 "first_epsilon_analytic_errors_checkpoint_staged_direct"->endpointChecks,
 "raw_result_sha256"->Association[Table[kind->sha[FileNameJoin[{base,kind<>"-order"<>ToString[order]<>".wl"}]],{kind,{"checkpoint","staged","direct"}}]]|>],
 {cfg,{{80,80},{100,112}}}];
refinementErrors=Table[Max[Abs[Flatten[Get[FileNameJoin[{base,kind<>"-order80.wl"}]][[2]]-Get[FileNameJoin[{base,kind<>"-order112.wl"}]][[2]]]]],{kind,{"checkpoint","staged","direct"}}];
If[!TrueQ[Max[refinementErrors]<10^-24],fail[{"Precision/order mismatch",refinementErrors}]];
exported=Export[FileNameJoin[{base,"raw-result.json"}],<|"status"->"passed","case"->"MG5 NP_EW1 closed masters 1,2",
 "wolfram_version"->$Version,"driver_sha256"->sha[driverFile],"package_sha256"->sha[package],"diffexp_commit"->"784c8229bf92369a03f011a48e161522c8c54bbd",
 "application_commit"->"89f64b93d0bdbd8ee90eb85737021189229b034a","matrix_sha256"->sha[matrixFile],"boundary_sha256"->sha[boundaryFile],
 "coordinate_names"->{"s","t","b"},"dependent_coordinate"-><|"name"->"u","expression"->"b-s-t"|>,
 "start"->Association[(str[First[#]]->str[Last[#]])&/@p0],"checkpoint"->Association[(str[First[#]]->str[Last[#]])&/@p1],"destination"->Association[(str[First[#]]->str[Last[#]])&/@p2],
 "boundary_values"->Map[scalar,values,{2}],"original_inherited_error"->str[rawbc[[3]]],"conservative_inherited_error"->"1e-27",
 "input_verified_digits"->24,"asserted_comparison_digits"->24,
 "connections_wolfram_input"->Map[str,connections,{3}],"epsilon_order"->4,"runs"->records,
 "shared_settings"-><|"mobius"->False,"pade"->False,"division_order"->3,"parallel"->False,"verbosity"->1,"accuracy_goal"->"?","delta_prescriptions"->{"s-I*delta","s-1-I*delta"}|>,
 "exact_returned_coordinates_checked"->True,
 "independent_first_epsilon_increment"->"[1,2]*Log[(s1-1)/(s0-1)]",
 "precision_order_max_absolute_differences"->(str/@refinementErrors),
 "scope"->"Closed two-master physical-data transport only; excludes other masters, roots, form-factor assembly and full amplitude acceptance.",
 "timing_scope"->"Single runs of TransportTo; configuration, package launch and boundary preparation excluded. No throughput or speedup claim.",
 "error_policy"->"Original grid error estimate rounded up to 1e-27 per coefficient. Input/output claims capped at 24 digits, regardless of mantissa length or smaller refinement differences."|>,"RawJSON"];
If[exported===$Failed,fail["Could not export result metadata"]];
Print["ORACLE PASS ",str[refinementErrors]];Exit[0],fail["DiffExp aborted"]];
