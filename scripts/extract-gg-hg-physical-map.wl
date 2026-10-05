SetSystemOptions["ParallelOptions" -> "ParallelThreadNumber" -> 1];
(* RustFlow mathematical-data extraction; no upstream solver code.
   Usage: WolframKernel -script this.wl PRIMARY_SCALAR_DIR PLUGIN_SYSTEM_DIR
          EXISTING_FIXTURE OUTPUT_JSON
   The existing fixture supplies declared routing and source provenance only.
   Every matrix coefficient is regenerated from the primary replacement rules. *)
args = If[Length[$ScriptCommandLine] > 0, Rest[$ScriptCommandLine],
  scriptPosition = FirstPosition[$CommandLine, "-script", Missing["none"], {1}];
  If[MissingQ[scriptPosition], {}, Take[Drop[$CommandLine, First[scriptPosition]+1], UpTo[4]]]];
If[Length[args] != 4, Print["Expected primary directory, plugin directory, existing fixture, output JSON"]; Exit[1]];
{primary, pluginRoot, templatePath, outputPath} = args;
template = Import[templatePath, "RawJSON"];
If[template["paper"]["sha256"] != "7026b05b0d8950c052a4642c3985cf03f71b838d57e25913a4f15a81fa823eb6", Print["Wrong source provenance"]; Exit[1]];
output = {};
Do[
  {name, short, n, head, file, directory} = case;
  source = primary <> "/topology_definitions_and_mappings/" <> file;
  rules = Get[source];
  Print[name, " rules=", Length[rules], " lhshead=", DeleteDuplicates[Head /@ rules[[All,1]]]];
  If[Length[rules] != n, Print["wrong row count"]; Exit[2]];
  rhs = rules[[All,2]] /. {s12 -> s, s13 -> -t, s23 -> s+t-b, mZ -> 1};
  sourceRhs = rhs /. myRoot[x_, 1/2] :> myRoot[Expand[x], 1/2];
  radicands = DeleteDuplicates[Cases[rhs, myRoot[x_, 1/2] :> Expand[x], Infinity]];
  rs = Table[Symbol["root" <> ToString[i]], {i, Length[radicands]}];
  rhs = rhs /. myRoot[x_, 1/2] :> rs[[First@FirstPosition[radicands, Expand[x], Missing["none"], {1}, Heads -> False]]];
  backRoots = Thread[rs -> (myRoot[#, 1/2]& /@ radicands)];
  rootReconstruction = And @@ (TrueQ[Together[#] === 0]& /@ ((rhs /. backRoots)-sourceRhs));
  If[!rootReconstruction, Print["root roundtrip failed"]; Exit[11]];
  Print[name, " exact original-root roundtrip=", rootReconstruction];
  coeff = Table[Together[Coefficient[Expand[rhs[[i]]], head[j]]], {i,n}, {j,n}];
  If[!And @@ Table[Together[rhs[[i]]-Sum[coeff[[i,j]]head[j],{j,n}]]===0,{i,n}],
     Print["source linear reconstruction failed"]; Exit[3]];
  Print[name, " nonzeros=", Count[coeff,x_/;x=!=0,{2}], " roots=", Length[radicands], " starting inverse"];
  bounds = {0}; reach = 0;
  Do[reach = Max[reach, Max[Flatten[Position[coeff[[i]], x_ /; x =!= 0, {1}]]]];
    If[reach == i, AppendTo[bounds,i]],{i,n}];
  If[Last[bounds] != n, Print["bad block structure"]; Exit[8]];
  Print["diagonal block sizes=", Differences[bounds]];
  inverse = ConstantArray[0,{n,n}];
  Do[
    first=bounds[[block]]+1; last=bounds[[block+1]];
    rows=Range[first,last];
    small=TimeConstrained[Map[Together, Inverse[coeff[[rows,rows]]],{2}],120,$Failed];
    If[small===$Failed, Print["small inverse timeout", rows];Exit[9]];
    inverse[[rows,rows]]=small;
    If[first>1,
      prior=Range[first-1];
      inverse[[rows,prior]]=Map[Together,-small.coeff[[rows,prior]].inverse[[prior,prior]],{2}]
    ];
    Print["inverted block ",block," rows ",first,"..",last];
    ,{block,Length[bounds]-1}];
  If[inverse === $Failed, Print["inverse timeout"]; Exit[4]];
  Print[name, " inverse computed"];
  inverse = Map[Together, inverse, {2}];
  
  (* Use explicit zero tests: zero is not a boolean identity. *)
  identity = And @@ (TrueQ[#===0]& /@ Flatten[Map[Together, coeff.inverse-IdentityMatrix[n],{2}]]);
  If[!identity, Print["inverse verification failed"]; Exit[5]];
  paperFile = primary <> "/" <> directory <> "/" <> If[n==48,"Planar","NonPlanar"] <> "_dlog_Matrix.txt";
  paperDE = ReadList[paperFile, Expression] /. mmH -> b;
  pluginFile = pluginRoot <> "/" <> name <> "/atilde_" <> name <> ".txt";
  pluginDE = Get[pluginFile];
  deDifference = Expand[paperDE-pluginDE];
  exactDE = Dimensions[paperDE] === {n,n} && And @@ (TrueQ[#===0]& /@ Flatten[deDifference]);
  Print[name, " exact published/plugin dlog equality=", exactDE];
  If[!exactDE, Print["published/plugin differential systems differ"]; Exit[6]];
  sparse[m_] := Table[Select[Table[<|"column"->j-1,"coefficient"->StringReplace[ToString[m[[i,j]],InputForm]," "->""]|>,{j,n}], #["coefficient"]!="0"&],{i,n}];
  AppendTo[output, Join[SelectFirst[template["families"], #["family"] == name&], <|"family"->name,"dimension"->n,
    "physical_integrals"->(List@@#& /@ rules[[All,1]]),
    "routing" -> SelectFirst[template["families"], #["family"] == name&]["routing"],
    "roots"->MapThread[<|"name"->ToString[#1,InputForm],"radicand"->StringReplace[ToString[#2,InputForm]," "->""]|>&,{rs,radicands}],
    "publication_to_canonical"->sparse[inverse],"canonical_to_publication"->sparse[coeff],
    "evidence"-><|"source_sha256"->FileHash[source,"SHA256","HexString"],
      "linear_reconstruction"->True,"exact_root_reconstruction"->rootReconstruction,"right_inverse"->identity,
      "published_dlog_sha256"->FileHash[paperFile,"SHA256","HexString"],
      "plugin_dlog_sha256"->FileHash[pluginFile,"SHA256","HexString"],
      "exact_dlog_equality"->exactDE|>|>]];

  , {case, {{"Planar_EW1","PL",48,canMIPlanar,"replRulePlanarMIToCan.m","planar"},{"NP_EW1","NP",61,canMINP,"replRuleNonPlanarMIToCan.m","non_planar"}}}];
result = Join[KeyDrop[template, {"families"}], <|"schema" -> "rustflow-gg-hg-plugin-physical-map-v2", "families" -> output|>];
Export[outputPath, result, "RawJSON"];
Exit[0];
