(* ::Package:: *)

(* ::Title:: *)
(*MaxCutClosedForm_LiteRed -- reference maximal-cut IBP reduction with LiteRed2 (R. N. Lee, arXiv:1212.2685, arXiv:1310.1145)*)


(* ::Text:: *)
(*This file is loaded by MaxCutClosedForm`MCCompare into a SEPARATE kernel (LaunchKernels[1]), so that LiteRed, its vector*)
(*declarations and its global settings never touch the session that runs the closed forms.  Do not load it into your own session.*)
(*Input of MCLiteRed`MCLiteRedReduce[job] (plain numbers, built by MCCompare from the MCTopology and the kinematic point):*)
(*   "L", "E" (loops, independent legs), "NProps", "Coefficients" (momentum of every denominator in the basis {loops, legs}),*)
(*   "Masses" (m2_i), "ExternalGram" (p_i.p_j), "Targets" (index vectors), "DValues" (numbers d), "WorkDir", "TimeLimit".*)
(*LiteRed basis: fresh symbols for the momenta, denominators D_i = sp[q_i] - m2_i (the convention of MaxCutClosedForm),*)
(*CutDs -> 1 on every propagator (every sub-sector vanishes: maximal cut), SectorsPattern -> ISP indices <= 0, SolvejSector -> True,*)
(*then IBPReduce of the targets.  Output: LiteRed's masters and the exact coefficient rows at the requested values of d.*)


(* ::Section:: *)
(*Load LiteRed2 (template cell of LiteRed suppressed: NextCell is blocked while loading)*)


LiteRed`$LiteRedLog = TrueQ[MCLiteRed`$Log];
LiteRed`$LiteRedMonitor = False;
If[!MemberQ[$Packages, "LiteRed`"],
  Block[{NextCell}, Quiet[Get[If[StringQ[MCLiteRed`$File] && FileExistsQ[MCLiteRed`$File], MCLiteRed`$File, "LiteRed2`"]]]]];
If[!MemberQ[$Packages, "LiteRed`"], Block[{NextCell}, Quiet[Get["LiteRed2024`"]]]];
MCLiteRed`$Status = If[MemberQ[$Packages, "LiteRed`"] && Names["LiteRed`NewDsBasis"] =!= {}, "loaded",
   "LiteRed2 could not be loaded (<<LiteRed2` failed in the reference kernel; set MaxCutClosedForm`$MCLiteRedFile to the full path of LiteRed2024.m)"];
MCLiteRed`$Version = If[MCLiteRed`$Status === "loaded", "LiteRed " <> ToString[$LiteRedVersion] <> " (" <> ToString[$LiteRedReleaseDate] <> ")", "none"];
If[!IntegerQ[MCLiteRed`$Counter], MCLiteRed`$Counter = 0];


(* ::Section:: *)
(*One reduction job*)


MCLiteRed`MCLiteRedReduce[job_Association] := Module[{odir = Directory[], res},
  res = CheckAbort[
    TimeConstrained[MCLiteRed`Private`run[job], Lookup[job, "TimeLimit", Infinity],
      <|"Status" -> "timeout", "Message" -> "LiteRed did not finish within " <> ToString[Lookup[job, "TimeLimit", Infinity]] <> " s"|>],
    <|"Status" -> "aborted", "Message" -> "LiteRed aborted the computation"|>];
  Quiet[SetDirectory[odir]];
  res];

MCLiteRed`Private`run[job_Association] := Catch[Module[
  {k, nL, nE, nP, co, nAll, vecs, lms, ems, dsym, dens, fam, dir, t0, t1, t2, targets, ints, red, used, mis, rows, rest, vals, example, symQ},
  If[MCLiteRed`$Status =!= "loaded", Throw[<|"Status" -> "failed", "Message" -> MCLiteRed`$Status|>, MCLiteRed`Private`tag]];
  k = ++MCLiteRed`$Counter;
  nL = job["L"]; nE = job["E"]; nP = job["NProps"]; co = job["Coefficients"]; nAll = Length[co];
  (* fresh symbols for every job: no scalar product of an earlier job can leak in *)
  vecs = Join[Table[Symbol["l" <> ToString[i] <> "J" <> ToString[k]], {i, nL}], Table[Symbol["p" <> ToString[i] <> "J" <> ToString[k]], {i, nE}]];
  lms = Take[vecs, nL]; ems = Drop[vecs, nL];
  dsym = Symbol["dJ" <> ToString[k]];
  SetDim[dsym];
  With[{v = vecs}, Declare[v, Vector]];
  Do[With[{a = ems[[i1]], b = ems[[i2]], val = job["ExternalGram"][[i1, i2]]}, sp[a, b] = val], {i1, nE}, {i2, i1, nE}];
  dens = Table[sp[Total[co[[i]] vecs]] - job["Masses"][[i]], {i, nAll}];
  fam = Symbol["famJ" <> ToString[k]];
  dir = FileNameJoin[{job["WorkDir"], "job" <> ToString[k] <> "_" <> ToString[$ProcessID]}];
  If[!DirectoryQ[dir], CreateDirectory[dir]];
  SetDirectory[dir];
  t0 = AbsoluteTime[];
  symQ = !(TrueQ[MCLiteRed`$NoSymmetries] || $VersionNumber == 14.1);
  With[{nm = fam, ds = dens, l = lms, cut = Join[ConstantArray[1, nP], ConstantArray[0, nAll - nP]],
        pat = Join[ConstantArray[_, nP], ConstantArray[0, nAll - nP]], bd = FileNameJoin[{dir, "basis"}]},
    If[symQ,
      (* standard LiteRed2 chain: GenerateIBP, AnalyzeSectors, FindSymmetries, SolvejSector *)
      NewDsBasis[nm, ds, l, CutDs -> cut, SectorsPattern -> pat, Directory -> bd, SolvejSector -> True],
      (* Mathematica 14.1: its Solve bug makes FindSymmetries hang (LiteRed2 discussion #9).  Without FindSymmetries every nonzero
         sector is unique and has no symmetry (the state FindSymmetries itself leaves when it finds none), then SolvejSector *)
      NewDsBasis[nm, ds, l, CutDs -> cut, SectorsPattern -> pat, Directory -> bd, GenerateIBP -> True, AnalyzeSectors -> True];
      If[TrueQ[DsBasisQ[nm]],
        nm /: UniqueSectors[nm] = NonZeroSectors[nm];
        Scan[(nm /: jSymmetries[Sequence @@ #] = {{}}) &, NonZeroSectors[nm]];
        SolvejSector[nm]]]];
  t1 = AbsoluteTime[];
  If[!TrueQ[With[{nm = fam}, DsBasisQ[nm]]],
    Throw[<|"Status" -> "failed", "Message" -> "NewDsBasis did not accept the denominators " <> ToString[dens, InputForm]|>, MCLiteRed`Private`tag]];
  mis = With[{nm = fam}, MIs[nm]];
  targets = job["Targets"];
  ints = With[{nm = fam}, (j[nm, ##] &) @@@ targets];
  red = IBPReduce[ints];
  red = red /. With[{nm = fam}, ZerojRule[nm]];
  t2 = AbsoluteTime[];
  used = Union[Cases[red, _j, Infinity]];
  If[!SubsetQ[mis, used],
    Throw[<|"Status" -> "failed", "Message" -> "LiteRed reduction incomplete; left over: " <> ToString[Complement[used, mis], InputForm],
      "Masters" -> (Rest[List @@ #] &) /@ mis|>, MCLiteRed`Private`tag]];
  rows = Table[Coefficient[red[[i]], m], {i, Length[targets]}, {m, mis}];
  rest = Together[red /. Thread[mis -> 0]];
  If[!AllTrue[rest, # === 0 &],
    Throw[<|"Status" -> "failed", "Message" -> "LiteRed result is not a combination of MIs (remainder " <> ToString[Short[rest, 3], InputForm] <> ")"|>, MCLiteRed`Private`tag]];
  vals = Association[Table[dv -> Together[rows /. dsym -> dv], {dv, job["DValues"]}]];
  example = If[rows === {}, "", StringReplace[ToString[Factor[First[rows]], InputForm], ToString[dsym] -> "d"]];
  ResetDirectory[];
  Quiet[DeleteDirectory[dir, DeleteContents -> True]];
  <|"Status" -> "ok", "Version" -> MCLiteRed`$Version <> ", Mathematica " <> ToString[$VersionNumber] <>
      If[symQ, "", ", without FindSymmetries"], "Masters" -> ((Rest[List @@ #] &) /@ mis), "Rows" -> vals,
    "Example" -> example, "TimeBasis" -> Round[t1 - t0, 0.01], "TimeReduce" -> Round[t2 - t1, 0.01], "NTargets" -> Length[targets]|>],
  MCLiteRed`Private`tag];
