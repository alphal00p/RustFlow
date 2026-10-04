(* Independent scientific-data export and line inventory. No DiffExp solver. *)
$HistoryLength=0;
driverFile=ExpandFileName[$InputFileName];
fail[m_]:=(Print["FAIL: ",m];Exit[1]);
requiredEnvironment[name_]:=Module[{value=Environment[name]},
 If[!StringQ[value]||StringLength[value]==0,fail["Set "<>name]];ExpandFileName[value]];
root=requiredEnvironment["DIFFEXP_ROOT"];base=requiredEnvironment["DIFFEXP_OUTPUT"];
If[FileExistsQ[base]||DirectoryQ[base],fail["DIFFEXP_OUTPUT must be new"]];
If[CreateDirectory[base,CreateIntermediateDirectories->True]===$Failed,fail["Cannot create output"]];
fixture=Import[FileNameJoin[{DirectoryName[DirectoryName[driverFile]],"fixtures","diffexp","banana-unequal.json"}],"RawJSON"];
If[!AssociationQ[fixture],fail["Missing input provenance fixture"]];
Do[file=FileNameJoin[{root,entry["file"]}];
 If[!FileExistsQ[file]||IntegerString[FileHash[file,"SHA256"],16,64]=!=entry["sha256"],fail["Pinned matrix mismatch"]],{entry,fixture["source_files"]}];

writeJSON[args___]:=If[Export[args]===$Failed,fail["JSON export failed"]];
str[v_]:=ToString[v,InputForm];sha[f_]:=IntegerString[FileHash[f,"SHA256"],16,64];
variables={psq,mm1,mm2,mm3,mm4};records={};all={};
Do[
 matrices=Table[file=FileNameJoin[{root,"Examples","Banana_Matrices","d"<>str[var]<>"_"<>ToString[k]<>".m"}];m=Get[file];
 If[Dimensions[m]=!={15,15},fail["Shape mismatch"]];
 AppendTo[records,<|"coordinate"->str[var],"epsilon_order"->k,"file"->"Examples/Banana_Matrices/"<>FileNameTake[file],"sha256"->sha[file],"bytes"->FileByteCount[file],"nonzero_entries"->Count[Flatten[m],Except[0]]|>];m,{k,0,4}];
 AppendTo[all,matrices],{var,variables}];
writeJSON[FileNameJoin[{base,"generic-matrices.json"}],<|"schema_version"->1,"dimension"->15,"epsilon_orders"->Range[0,4],"coordinates"->(str/@variables),"coordinate_meaning"->{"p^2","m1^2","m2^2","m3^2","m4^2"},"ordering"->"matrices[coordinate][epsilon_order][row][column]","matrices"->Map[str,all,{4}],"sources"->records|>,"RawJSON"];
Print["Generic export complete"];
lines={<|"name"->"mass_deformation_half","coordinates"->{1/2,1+x,1+x/2,1+x/3,1},"start"->0,"end"->1|>,<|"name"->"physical_psq","coordinates"->{x,2,3/2,4/3,1},"start"->1/2,"end"->50|>,<|"name"->"mass_deformation_50","coordinates"->{50,1+x,1+x/2,1+x/3,1},"start"->0,"end"->1|>};
lineRecords={};
Do[
 coords=line["coordinates"];rules=Thread[variables->coords];
 seconds=AbsoluteTiming[connection=Table[Map[Together,Sum[D[coords[[j]],x]*(all[[j,k]]/.rules),{j,1,5}],{2}],{k,1,5}];][[1]];
 If[Dimensions[connection]=!={5,15,15},fail["Line shape"]];
 denoms=DeleteDuplicates[Denominator[Together[#]]&/@Flatten[connection]];
 startZeros=Select[denoms,TrueQ[(#/.x->line["start"])==0]&];
 denomProduct=Fold[PolynomialLCM,1,denoms];
 degrees=Table[Max[0,Max[Exponent[#,x]&/@Flatten[Map[If[#===0,0,Numerator[#]]&,connection[[k]],{2}]]]],{k,1,5}];
 record=<|"name"->line["name"],"coordinates"->Association[Thread[(str/@variables)->(str/@coords)]],"start"->str[line["start"]],"end"->str[line["end"]],"matrix_by_epsilon"->Map[str,connection,{3}],"line_construction_seconds"->seconds,"nonzero_entries_by_epsilon"->(Count[Flatten[#],Except[0]]&/@connection),"numerator_max_degrees_by_epsilon"->degrees,"denominator_lcm"->str[Factor[denomProduct]],"denominator_lcm_degree"->Exponent[denomProduct,x],"distinct_denominators"->Length[denoms],"denominators_vanishing_at_start"->(str/@startZeros)|>;
 writeJSON[FileNameJoin[{base,line["name"]<>".json"}],record,"RawJSON"];
 AppendTo[lineRecords,KeyDrop[record,{"matrix_by_epsilon"}]];
 Print["LINE PASS ",line["name"]," seconds ",seconds," start zero denominators ",Length[startZeros]],{line,lines}];
writeJSON[FileNameJoin[{base,"matrix-line-inventory.json"}],<|"sources"->records,"lines"->lineRecords,"generic_bytes"->Total[Lookup[records,"bytes"]],"upstream_commit"->"784c8229bf92369a03f011a48e161522c8c54bbd","notebook_sha256"->sha[FileNameJoin[{root,"Examples","Banana.nb"}]],"paper_sha256"->sha[FileNameJoin[{root,"Paper","main.tex"}]],"export_driver_sha256"->sha[driverFile],"numerical_transport_run"->False|>,"RawJSON"];
Print["EXPORT PASS"];Exit[0];
