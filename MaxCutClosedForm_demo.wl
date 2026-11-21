(* ::Package:: *)

(* ::Title:: *)
(*MaxCutClosedForm -- demo*)


(* ::Text:: *)
(*Closed forms, in symbolic propagator powers, of the maximal-cut IBP reduction coefficients for topologies with up to 2 ISPs*)
(*(standard Baikov representation).*)
(*Workflow:  topo = MCTopology[loops, legs, propagators, ISPs, kinematics]  ->  res = MCClosedForm[topo, point]  ->*)
(*           MCReduce[res, n, d] (coefficients), MCFormula[res] (closed formulas), MCInfo[res], MCTeX[res, file], MCCheck[res],*)
(*           MCCompare[res, targets] (comparison with LiteRed2).*)
(*Index vectors n list the propagator powers first (>= 1; n - 1 = dots), then the ISP powers (<= 0; -n = numerator power).*)
(**)
(*Independent check of every topology: MCCompare runs LiteRed2 (R. N. Lee) on THIS computer, in a separate kernel, with CutDs on every*)
(*propagator (maximal cut), reduces the listed integrals plus a random sample, and compares LiteRed's exact coefficients with the closed*)
(*forms at two values of d.  Nothing is precomputed: change the kinematic point, the targets or the topology and LiteRed is rerun.*)
(*LiteRed results are cached in the folder MaxCutClosedForm_LiteRed/cache (delete it to rerun LiteRed); the first run takes longer.*)
(**)
(*Files (same folder): MaxCutClosedForm.m (package), MaxCutClosedForm_LiteRed.wl (LiteRed side, loaded only into the separate kernel).*)
(*Needs LiteRed2 installed so that <<LiteRed2` works (otherwise set $MCLiteRedFile = "/full/path/LiteRed2024.m" after loading).*)
(*Run:  in a terminal   cd ~/Documents/Claude/IBPrec/Mathematica && wolframscript -file MaxCutClosedForm_demo.wl*)
(*      or open this file in Mathematica and choose Evaluation > Evaluate Notebook.*)
(*Everything printed and every message also goes to MaxCutClosedForm_demo_log.txt.  Each step is timed and has a time limit*)
(*(mcStep); a step that reaches its limit is reported in the log and the demo goes on with the next one.*)


(* ::Section:: *)
(*Loading*)


ClearAll[k, l, p, p1, p2, p3, s, t, u, mm, M2];
Quiet[Remove["Global`MCTopology", "Global`MCClosedForm", "Global`MCReduce", "Global`MCMasters", "Global`MCInfo", "Global`MCFormula",
  "Global`MCTeX", "Global`MCCheck", "Global`MCSelfTest", "Global`MCLauricellaFD", "Global`MCResult", "Global`MCz", "Global`MCCompare",
  "Global`MCCompareSummary", "Global`MCStopReference", "Global`$MCLiteRedFile", "Global`$MCReferenceDirectory"]];

mcDemoDir = Which[
   StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName],
   StringQ[Quiet[Check[NotebookDirectory[], $Failed]]], NotebookDirectory[],
   DirectoryQ[FileNameJoin[{$HomeDirectory, "Documents", "Claude", "IBPrec", "Mathematica"}]],
     FileNameJoin[{$HomeDirectory, "Documents", "Claude", "IBPrec", "Mathematica"}],
   True, Directory[]];
SetDirectory[mcDemoDir];


(* ::Text:: *)
(*Log file: Print output goes to it as plain text (the stream is closed and reopened after every step, so the file is always*)
(*complete up to the running step); messages are written to it by a message handler.*)


mcLogFile = FileNameJoin[{mcDemoDir, "MaxCutClosedForm_demo_log.txt"}];
mcLogOpen[opener_] := (mcLog = opener[mcLogFile, FormatType -> OutputForm, PageWidth -> 180, CharacterEncoding -> "UTF8"];
   $Output = Append[DeleteCases[$Output, os_OutputStream /; First[os] === mcLogFile], mcLog];);
mcLogClose[] := ($Output = DeleteCases[$Output, os_OutputStream /; First[os] === mcLogFile];
   Quiet[Close /@ Select[Streams[], First[#] === mcLogFile &]];);
mcFlush[] := (mcLogClose[]; mcLogOpen[OpenAppend]);
mcLogClose[]; mcLogOpen[OpenWrite];

mcMsgHandler[Hold[Message[mn : MessageName[sym_, tag_String], args___], True]] := Quiet[Module[{templ},
    templ = MessageName[sym, tag];
    If[!StringQ[templ], templ = MessageName[General, tag]];
    If[!StringQ[templ], templ = "(no text)"];
    WriteString[mcLog, "MESSAGE  ", ToString[HoldForm[mn], OutputForm], ": ",
      ToString[StringForm[templ, Sequence @@ (Function[x, Short[HoldForm[x], 2], HoldAllComplete] /@ Unevaluated[{args}])], OutputForm], "\n"]]];
mcMsgHandler[___] := Null;
Quiet[Internal`RemoveHandler["Message", mcMsgHandler]];
Internal`AddHandler["Message", mcMsgHandler];

(* one step of the demo: label, start time, time limit, duration; the log is flushed before and after *)
mcSec[x_] := ToString[NumberForm[N[x], {8, 2}]];
SetAttributes[mcStep, HoldRest];
mcStep[label_String, expr_, limit_ : 900] := Module[{t0 = AbsoluteTime[], r},
  Print[""]; Print[">> ", label, "   [", DateString[{"Hour", ":", "Minute", ":", "Second"}], ", limit ", limit, " s]"];
  mcFlush[];
  r = TimeConstrained[expr, limit, $TimedOut];
  Print["<< ", label, ": ", If[r === $TimedOut, "TIME LIMIT REACHED (" <> ToString[limit] <> " s) -- step skipped, the demo goes on",
     "done in " <> mcSec[AbsoluteTime[] - t0] <> " s"]];
  mcFlush[];
  r];
(* LiteRed comparisons: if the limit is reached, the reference kernel is closed (a new one is launched by the next MCCompare) *)
SetAttributes[mcCompareStep, HoldRest];
mcCompareStep[label_String, expr_, limit_ : 2400] := Module[{r = mcStep[label, expr, limit]}, If[r === $TimedOut, MCStopReference[]]; r];

Print["MaxCutClosedForm demo, ", DateString[], ";  folder ", mcDemoDir, ";  Mathematica ", $Version];


mcPackageFile = Which[
   FileExistsQ[FileNameJoin[{mcDemoDir, "MaxCutClosedForm.m"}]], FileNameJoin[{mcDemoDir, "MaxCutClosedForm.m"}],
   FindFile["MaxCutClosedForm`"] =!= $Failed, FindFile["MaxCutClosedForm`"],
   True, $Failed];
If[mcPackageFile === $Failed,
  Print["MaxCutClosedForm.m not found: put it in ", mcDemoDir, " (the folder of this demo)"]; mcFlush[]; Abort[]];
Get[mcPackageFile];
Print["loaded ", mcPackageFile, ";  MCTopology is in context ", Context["MCTopology"], "  (must be MaxCutClosedForm`)"];
Print["LiteRed side: ", FileNameJoin[{DirectoryName[mcPackageFile], "MaxCutClosedForm_LiteRed.wl"}],
  If[FileExistsQ[FileNameJoin[{DirectoryName[mcPackageFile], "MaxCutClosedForm_LiteRed.wl"}]], "  (found)", "  (MISSING: MCCompare will not work)"]];
dDemo = 37/11;


(* ::Text:: *)
(*Test of the numerical Lauricella F_D (series, continued Euler integral) and of the interval formula against direct quadrature.*)


mcStep["self-test of the special functions", MCSelfTest[], 600];


(* ::Section:: *)
(*1. Quick start: massless planar double box (2 ISPs)*)


(* ::Subsection:: *)
(*Definition and analysis of the maximal cut*)


(* propagators k^2, l^2, (k-l)^2, (k+p1)^2, (k+p1+p2)^2, (l+p1+p2)^2, (l-p3)^2; ISPs (l+p1)^2, (k-p3)^2; massless legs, u = (p1+p3)^2 *)
dbox = MCTopology[{k, l}, {p1, p2, p3},
   {k, l, k - l, k + p1, k + p1 + p2, l + p1 + p2, l - p3},
   {l + p1, k - p3},
   {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p1 p2 -> s/2, p1 p3 -> u/2, p2 p3 -> -(s + u)/2},
   "Name" -> "planar double box"];

resDbox = mcStep["planar double box: MCClosedForm", MCClosedForm[dbox, {s -> 7, u -> 3}], 600];   (* P_S, number of masters, cycles, masters *)

mcStep["planar double box: MCInfo", MCInfo[resDbox]; Print["masters: ", MCMasters[resDbox]], 120];


(* ::Subsection:: *)
(*The closed formula in symbolic ISP powers*)


(* G_C(a, b; rho): a = power of the block variable x, b = power of t, rho = gamma - j; one Lauricella F_D per mu *)
mcStep["planar double box: closed formulas",
  Print[MCFormula[resDbox, "G"]];
  Print[MCFormula[resDbox, "Points"]];
  Print[MCFormula[resDbox, "F"]];
  Print[MCFormula[resDbox, "Pi"][{2, 1, 1, 1, 1, 1, 1}]],      (* residue truncation Pi_j for a dot on the first propagator *)
  120];


(* ::Subsection:: *)
(*Reduction coefficients*)


mcStep["planar double box: reduction coefficients",
  Print[MCReduce[resDbox, {1, 1, 1, 1, 1, 1, 1, -3, -2}, dDemo]];          (* numerator (l+p1)^6 (k-p3)^4 *)
  Print[MCReduce[resDbox, {2, 1, 1, 1, 1, 1, 1, -1, 0}, dDemo]];           (* a dot and a numerator *)
  (* coefficients of the second master for ISP powers 0..4 *)
  Print[TableForm[N[Table[MCReduce[resDbox, {1, 1, 1, 1, 1, 1, 1, -a8, -a9}, dDemo][[2]], {a8, 0, 4}, {a9, 0, 4}], 12],
    TableHeadings -> {Range[0, 4], Range[0, 4]}]],
  300];


(* ::Subsection:: *)
(*Checks: LiteRed2 (independent), then the internal checks (IBP identities, modular Laporta)*)


(* the listed integrals + 8 random ones (dots <= 2, rank <= 4), at d = 37/11 and 43/13 *)
mcCompareStep["planar double box vs LiteRed",
  MCCompare[resDbox, {{1, 1, 1, 1, 1, 1, 1, -1, 0}, {1, 1, 1, 1, 1, 1, 1, -2, -1}, {1, 1, 1, 1, 1, 1, 1, 0, -3}, {2, 1, 1, 1, 1, 1, 1, 0, 0},
     {1, 1, 1, 2, 1, 1, 1, -1, 0}, {1, 1, 1, 1, 1, 1, 1, -1, -2}, {1, 1, 1, 1, 1, 1, 1, -3, -2}, {2, 1, 1, 1, 1, 1, 1, -1, 0}}]];
mcStep["planar double box: MCCheck", MCCheck[resDbox, Method -> "Both", "Grid" -> {1, 2}], 600];


(* ::Subsection:: *)
(*LaTeX file with the closed formulas*)


mcStep["planar double box: MCTeX", Print["LaTeX written to ", MCTeX[resDbox, FileNameJoin[{mcDemoDir, "dbox_closed_form.tex"}]]], 120];


(* ::Subsection:: *)
(*Own choice of masters (here with a dot); compared with LiteRed in OUR basis ("Masters" -> "Ours": I(n) = c(n).L(masters))*)


resDbox2 = mcStep["planar double box, dotted master: MCClosedForm",
   MCClosedForm[dbox, {s -> 7, u -> 3}, "Masters" -> {{1, 1, 1, 1, 1, 1, 1, 0, 0}, {1, 1, 2, 1, 1, 1, 1, 0, 0}}], 600];
mcCompareStep["planar double box, dotted master vs LiteRed",
  MCCompare[resDbox2, {{1, 1, 1, 1, 1, 1, 1, 0, -1}, {1, 1, 1, 1, 1, 1, 1, -2, -1}, {1, 1, 1, 2, 1, 1, 1, -1, 0}, {2, 1, 1, 1, 1, 1, 1, 0, 0}},
    "Masters" -> "Ours"]];


(* ::Section:: *)
(*2. One-loop topologies*)


(* ::Subsection:: *)
(*Triangle with an ISP (1 ISP): masses 1, 2, 3; p1^2 = 1, p2^2 = 2, p3^2 = 5*)


tri = MCTopology[{k}, {p1, p2, p3}, {{k, 1}, {k + p1, 2}, {k + p1 + p2, 3}}, {k + p1 + p2 + p3},
   {p1^2 -> 1, p2^2 -> 2, p3^2 -> 5, p1 p2 -> (s - 3)/2, p2 p3 -> (t - 7)/2, p1 p3 -> (u - 6)/2}, "Name" -> "massive triangle"];
resTri = mcStep["massive triangle: MCClosedForm", MCClosedForm[tri, {s -> 7, t -> -2, u -> 3}], 600];
mcStep["massive triangle: closed formula", Print[MCFormula[resTri, "G"]], 120];
mcCompareStep["massive triangle vs LiteRed",
  MCCompare[resTri, {{2, 1, 1, 0}, {1, 1, 1, -1}, {1, 1, 1, -3}, {1, 2, 2, -2}, {3, 1, 1, -1}}]];
mcStep["massive triangle: MCCheck", MCCheck[resTri, Method -> "Both", "Grid" -> {2, 3}], 600];


(* ::Subsection:: *)
(*Box (0 ISPs): the reduction is a rational function of d (exact comparison)*)


box = MCTopology[{k}, {p1, p2, p3}, {{k, 1}, {k + p1, 2}, {k + p1 + p2, 3}, {k + p1 + p2 + p3, 5}}, {},
   {p1^2 -> 1, p2^2 -> 2, p3^2 -> 5, p1 p2 -> (s - 3)/2, p2 p3 -> (t - 7)/2, p1 p3 -> (u - 6)/2}, "Name" -> "massive box"];
resBox = mcStep["massive box: MCClosedForm", MCClosedForm[box, {s -> 7, t -> -2, u -> 3}], 600];
mcStep["massive box: a coefficient", Print[MCReduce[resBox, {2, 1, 1, 1}, dDemo]], 120];
mcCompareStep["massive box vs LiteRed", MCCompare[resBox, {{2, 1, 1, 1}, {1, 1, 2, 2}, {3, 1, 1, 1}, {1, 2, 1, 3}}]];


(* ::Section:: *)
(*3. Two-loop vertices (1 ISP)*)


(* ::Subsection:: *)
(*Massless crossed vertex*)


xtri = MCTopology[{k, l}, {p1, p2}, {k, k + p1, l, l + p2, k + l, k + l + p1 + p2}, {k + p2},
   {p1^2 -> 0, p2^2 -> 0, p1 p2 -> s/2}, "Name" -> "crossed vertex"];
resXtri = mcStep["crossed vertex: MCClosedForm", MCClosedForm[xtri, {s -> 7}], 600];
mcCompareStep["crossed vertex vs LiteRed",
  MCCompare[resXtri, {{1, 1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1, -3}, {2, 1, 1, 1, 1, 1, 0}, {1, 2, 1, 1, 1, 1, -2}}]];
mcStep["crossed vertex: MCCheck", MCCheck[resXtri, Method -> "Both", "Grid" -> {1, 3}], 600];


(* ::Subsection:: *)
(*Massive crossed vertex: 3 masters, two of the four singular points complex*)


xtrim = MCTopology[{k, l}, {p1, p2}, {{k, 1}, k + p1, {l, 2}, l + p2, {k + l, 3}, {k + l + p1 + p2, 1}}, {k + p2},
   {p1^2 -> 0, p2^2 -> 2, p1 p2 -> (s - 2)/2}, "Name" -> "massive crossed vertex"];
resXtrim = mcStep["massive crossed vertex: MCClosedForm", MCClosedForm[xtrim, {s -> 7}], 600];
mcStep["massive crossed vertex: MCInfo", MCInfo[resXtrim], 120];
mcCompareStep["massive crossed vertex vs LiteRed",
  MCCompare[resXtrim, {{1, 1, 1, 1, 1, 1, -3}, {2, 1, 1, 1, 1, 1, 0}, {1, 1, 1, 1, 2, 1, -1}, {1, 1, 1, 1, 1, 1, -4}}]];


(* ::Subsection:: *)
(*Massless planar ladder vertex: P_S has a single root, no master (everything reduces to sub-sectors)*)


ladder = MCTopology[{k, l}, {p1, p2}, {k, k + p1, k + p1 + p2, k - l, l, l + p1 + p2}, {l + p1},
   {p1^2 -> 0, p2^2 -> 0, p1 p2 -> s/2}, "Name" -> "planar ladder vertex"];
resLadder = mcStep["planar ladder vertex: MCClosedForm", MCClosedForm[ladder, {s -> 7}], 600];
mcStep["planar ladder vertex: a reduction", Print[MCReduce[resLadder, {1, 1, 1, 1, 1, 1, -2}, dDemo]], 60];     (* {}: zero on the maximal cut *)
mcCompareStep["planar ladder vertex vs LiteRed (must reduce everything to zero)",
  MCCompare[resLadder, {{1, 1, 1, 1, 1, 1, 0}, {2, 1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1, -2}}]];
mcStep["planar ladder vertex: MCCheck", MCCheck[resLadder, Method -> "Laporta", "Grid" -> {1, 2}], 600];


(* ::Section:: *)
(*4. Two-loop self-energies*)


(* ::Subsection:: *)
(*Equal-mass kite (0 ISPs)*)


kite = MCTopology[{k, l}, {p}, {{k, 1}, {l, 1}, {k - p, 1}, {l - p, 1}, {k - l, 1}}, {}, {p^2 -> s}, "Name" -> "kite"];
resKite = mcStep["kite: MCClosedForm", MCClosedForm[kite, {s -> 7}], 600];
mcCompareStep["kite vs LiteRed", MCCompare[resKite, {{2, 1, 1, 1, 1}, {1, 1, 2, 1, 2}, {1, 1, 1, 1, 3}}]];


(* ::Subsection:: *)
(*Sunrise with three different masses (2 ISPs, 4 masters; cycles from both fibrations)*)


sunrise = MCTopology[{k, l}, {p}, {{k, 1}, {k - l, 2}, {l - p, 3}}, {k - p, l}, {p^2 -> s}, "Name" -> "sunrise"];
resSun = mcStep["sunrise, s = 50: MCClosedForm", MCClosedForm[sunrise, {s -> 50}], 600];
mcStep["sunrise, s = 50: MCInfo", MCInfo[resSun], 120];
mcCompareStep["sunrise, s = 50 vs LiteRed",
  MCCompare[resSun, {{1, 1, 1, -1, -1}, {1, 1, 1, -2, 0}, {2, 1, 1, 0, 0}, {1, 2, 1, -1, 0}, {1, 1, 1, -2, -1}}]];
mcStep["sunrise, s = 50: MCCheck", MCCheck[resSun, Method -> "Both", "Grid" -> {1, 2}], 900];


(* ::Subsection:: *)
(*Same sunrise below threshold, s = -5*)


resSunE = mcStep["sunrise, s = -5: MCClosedForm", MCClosedForm[sunrise, {s -> -5}], 600];
mcCompareStep["sunrise, s = -5 vs LiteRed",
  MCCompare[resSunE, {{1, 1, 1, -1, -1}, {1, 1, 1, -2, 0}, {2, 1, 1, 0, 0}, {1, 1, 2, 0, -1}}]];


(* ::Section:: *)
(*5. Double boxes*)


(* ::Subsection:: *)
(*Massless crossed double box*)


xbox = MCTopology[{k, l}, {p1, p2, p3}, {k, k - p1, k + p2, l - k, l - p1, l + p2, l + p2 + p3}, {k + p3, l + p1 + p2},
   {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p1 p2 -> s/2, p1 p3 -> -(s + t)/2, p2 p3 -> t/2}, "Name" -> "crossed double box"];
resXbox = mcStep["crossed double box: MCClosedForm", MCClosedForm[xbox, {s -> 7, t -> 3}], 600];
mcCompareStep["crossed double box vs LiteRed",
  MCCompare[resXbox, {{1, 1, 1, 1, 1, 1, 1, -1, 0}, {1, 1, 1, 1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, 1, 0, -2}, {1, 1, 2, 1, 1, 1, 1, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, -2, -1}}]];
mcStep["crossed double box: MCCheck", MCCheck[resXbox, "Grid" -> {1, 2}], 600];


(* ::Subsection:: *)
(*Planar double box with a massive rung (k-l)^2 - mm (3 masters)*)


rung = MCTopology[{k, l}, {p1, p2, p3}, {k, l, {k - l, mm}, k + p1, k + p1 + p2, l + p1 + p2, l - p3}, {l + p1, k - p3},
   {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p1 p2 -> s/2, p1 p3 -> u/2, p2 p3 -> -(s + u)/2}, "Name" -> "double box, massive rung"];
resRung = mcStep["massive rung: MCClosedForm", MCClosedForm[rung, {s -> 7, u -> 3, mm -> 2}], 600];
mcStep["massive rung: closed formula", Print[MCFormula[resRung, "G"][[1]]], 120];
mcCompareStep["massive rung vs LiteRed",
  MCCompare[resRung, {{1, 1, 1, 1, 1, 1, 1, -1, 0}, {1, 1, 1, 1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, 1, 0, -3}, {1, 1, 2, 1, 1, 1, 1, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, -2, -1}}]];


(* ::Subsection:: *)
(*Planar double box with a massive k-loop (4 masters; complex roots of A).  Slower: the Lauricella functions need the continued Euler integral*)


kloop = MCTopology[{k, l}, {p1, p2, p3}, {{k, mm}, l, k - l, {k + p1, mm}, {k + p1 + p2, mm}, l + p1 + p2, l - p3}, {l + p1, k - p3},
   {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p1 p2 -> s/2, p1 p3 -> u/2, p2 p3 -> -(s + u)/2}, "Name" -> "double box, massive k-loop"];
resKloop = mcStep["massive k-loop: MCClosedForm", MCClosedForm[kloop, {s -> 7, u -> 3, mm -> 2}], 900];
mcCompareStep["massive k-loop vs LiteRed",
  MCCompare[resKloop, {{1, 1, 1, 1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, 1, -2, 0}, {2, 1, 1, 1, 1, 1, 1, 0, 0}, {1, 1, 1, 1, 1, 1, 1, -1, -2}}]];


(* ::Section:: *)
(*6. A topology that was not used while the package was written: planar double box with one massive leg, p4^2 = M2*)


(* ::Text:: *)
(*Your own topology works the same way: MCTopology[loops, legs, propagators ({q, m2} or q), ISPs, rules for all external scalar*)
(*products] -> MCClosedForm[topo, point with exact numbers] -> MCCompare[res, targets].  d stays symbolic in MCFormula and is a*)
(*number in MCReduce.  Options of MCClosedForm: "Masters" -> {n1, ...}, WorkingPrecision -> 30, "TestDimension" -> 25/7, "Verbose" -> True.*)


(* p4 = -(p1 + p2 + p3):  p4^2 = s + t + u = M2,  t = (p2 + p3)^2 = M2 - s - u *)
dbox1m = MCTopology[{k, l}, {p1, p2, p3}, {k, l, k - l, k + p1, k + p1 + p2, l + p1 + p2, l - p3}, {l + p1, k - p3},
   {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p1 p2 -> s/2, p1 p3 -> u/2, p2 p3 -> (M2 - s - u)/2}, "Name" -> "planar double box, p4^2 = M2"];
resDbox1m = mcStep["one-mass double box: MCClosedForm", MCClosedForm[dbox1m, {s -> 7, u -> 3, M2 -> 2}], 900];
mcCompareStep["one-mass double box vs LiteRed",
  MCCompare[resDbox1m, {{1, 1, 1, 1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, 1, -3, 0}, {2, 1, 1, 1, 1, 1, 1, 0, -1}}, "Sample" -> {10, {2, 5}}]];


(* ::Section:: *)
(*7. Summary of all comparisons with LiteRed*)


mcStep["summary", MCCompareSummary[]; MCStopReference[], 120];
Print["demo finished ", DateString[]];
Quiet[Internal`RemoveHandler["Message", mcMsgHandler]];
mcLogClose[];
