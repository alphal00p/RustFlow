$HistoryLength=0;
root=Environment["RUSTFLOW_PRISM_SYMBOLIC_DIR"];
r=Get[FileNameJoin[{root,"reduced-combinations.wl"}]];
m=r["surviving_masters"];
rows=MapThread[Function[{label,original,reduced},<|"label"->label,"original_terms"->(Function[t,<|"indices"->Rest[List@@t],"coefficient"->ToString[Factor[Coefficient[Expand[original],t]],InputForm]|>]/@DeleteDuplicates[Cases[original,_j,Infinity]]),"candidate_coefficients"->(ToString[Factor[Coefficient[Expand[reduced],#]],InputForm]&/@m)|>],{r["labels"],r["original"],r["reduced"]}];
Export[FileNameJoin[{root,"candidate-relations.json"}],<|"schema"->1,"status"->"external_exact_IBP_candidate_pending_native_replay","origin"->"bounded optional upstream AMFlow/Kira symbolic reduction; not required feature or validation dependency","dimension"->"4-2*eps","master_indices"->(Rest[List@@#]&/@m),"relations"->rows,"complete_prism_reference"->False,"native_predictions_read"->0,"oracle_records_read"->0|>,"RawJSON"];
Print["EXPORTED_CANDIDATE_RELATIONS ",Length[rows]];Exit[0];
