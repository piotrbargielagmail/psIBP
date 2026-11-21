(* ::Package:: *)

(* ::Title:: *)
(*MaxCutClosedForm -- closed forms in symbolic propagator powers of maximal-cut IBP reduction coefficients (up to 2 ISPs)*)


(* ::Text:: *)
(*Usage (see MaxCutClosedForm_demo.wl):*)
(*    topo = MCTopology[{k, l}, {p1, p2, p3}, {k, l, k - l, k + p1, k + p1 + p2, l + p1 + p2, l - p3}, {l + p1, k - p3},*)
(*                      {p1^2 -> 0, p2^2 -> 0, p3^2 -> 0, p1 p2 -> s/2, p1 p3 -> u/2, p2 p3 -> -(s + u)/2}];*)
(*    res  = MCClosedForm[topo, {s -> 7, u -> 3}];        (* masters, cycles, closed formulas *)*)
(*    MCReduce[res, {2, 1, 1, 1, 1, 1, 1, -3, -1}, 37/11]  (* coefficients of the masters *)*)
(*    MCFormula[res], MCInfo[res], MCTeX[res, "file.tex"], MCCheck[res],*)
(*    MCCompare[res, targets]   (* against LiteRed2, run in a separate kernel; needs MaxCutClosedForm_LiteRed.wl next to this file *)*)
(**)
(*Conventions.  Propagators D_i = q_i^2 - m_i^2 (an entry q means m^2 = 0, an entry {q, m2} means q^2 - m2), integrals*)
(*  I(n) = Integral prod_i D_i^(-n_i),  n_i >= 1 for the propagators (n_i - 1 = dots), n_i <= 0 for the ISPs (-n_i = numerator powers).*)
(*Maximal cut of the top sector, standard Baikov representation (P = Gram determinant of loop and external momenta in the variables*)
(*z_i = D_i, gamma = (d - L - E - 1)/2), cut kernel P_S = P(z_S = 0):*)
(*  [z_S^(n_S - 1)] P^gamma = sum_j Binomial(gamma, j) P_S^(gamma - j) Pi_j(z_Y),   Pi_j = [z_S^(n_S - 1)] (P - P_S)^j     (residue truncation),*)
(*  F_C(n) = sum_j Binomial(gamma, j) sigma_C^j sum_(alpha) [Pi_j]_alpha G_C(k + alpha; gamma - j),   k = ISP powers,*)
(*  A(n) = F(n) . M^(-1),  M_(i C) = F_C(master_i)                                                   (reduction coefficients).*)
(*G_C = cycle integral of an ISP monomial against |P_S|^rho:*)
(*  0 ISPs: no integral (F = P_S^(-j) Pi_j; one master).*)
(*  1 ISP : G_C(a; rho) = int_C t^a |P_S(t)|^rho dt.*)
(*  2 ISPs: P_S = A(t) x^2 + 2 B(t) x + C(t) (Gaussian block in x), Delta = B^2 - A C,*)
(*          G_C(a, b; rho) = sum_mu Binomial(a, 2 mu) Gamma(mu + 1/2) Gamma(rho + 1)/Gamma(rho + mu + 3/2)*)
(*                           int_C t^b (-B)^(a - 2 mu) |A|^(-rho - 1) A^(-a) |Delta|^(rho + 1/2) Delta^mu dt.*)
(*The remaining t-integral is a product of powers of linear factors (t - t_k): roots of A, Delta, B, t (and of P_S for 1 ISP);*)
(*it is ONE Lauricella F_D per mu with parameters linear in (a, b, rho, mu):  closed form in all ISP powers.  Cycles C: real chambers*)
(*between consecutive real singular points (both fibrations) and straight complex segments; they are selected by rank.  The number*)
(*of masters is the Euler characteristic of {z_Y : P_S != 0}; masters are ISP monomials of lowest degree unless given.*)


BeginPackage["MaxCutClosedForm`"];

MCTopology::usage = "MCTopology[loops, legs, propagators, isps, kinematics] defines an integral family.
loops = {k, l, ...}: loop momenta; legs = {p1, p2, ...}: independent external momenta;
propagators = {q1, {q2, m2}, ...}: an entry q means 1/q^2, {q, m2} means 1/(q^2 - m2); isps: the irreducible scalar products, same format;
kinematics = rules for ALL external scalar products, e.g. {p1^2 -> 0, p1 p2 -> s/2}.
Propagators + ISPs must number L(L+1)/2 + L E and determine all scalar products.  Option \"Name\" -> string.";
MCClosedForm::usage = "MCClosedForm[topo, point] analyses the maximal cut of the top sector of topo at the kinematic point (list of rules
with exact numbers), counts the master integrals, selects masters and integration cycles and returns an MCResult object holding the
closed forms.  Options: \"Masters\" -> Automatic | {n1, n2, ...} (index vectors), \"Dimension\" -> symbol used for d in formulas,
\"TestDimension\" -> 25/7 (generic d used for the numerical selection), WorkingPrecision -> 30, \"Verbose\" -> True.";
MCReduce::usage = "MCReduce[res, n, d] gives the maximal-cut reduction coefficients {c1, ..., cM}, I(n) = c1 I(master1) + ... + cM I(masterM)
(up to sub-sector integrals), at the number d, evaluated from the closed forms.  n: index vector (propagator powers >= 1, ISP powers <= 0)
or a list of index vectors.  d is converted to an exact number; the precision is the WorkingPrecision given to MCClosedForm.
0 ISPs: the result is exact (rational in d).";
MCMasters::usage = "MCMasters[res] gives the master integrals as index vectors.";
MCInfo::usage = "MCInfo[res] prints the cut Baikov polynomial, the fibrations, the singular points with their exponents, the cycles and the masters.";
MCFormula::usage = "MCFormula[res] gives an association of the closed formulas (keys \"gamma\", \"PS\", \"Fibrations\", \"Points\", \"G\", \"F\", \"A\",
\"Masters\", \"Cycles\").  In \"G\" the ISP powers are the formal symbols a, b, the twist is rho (= gamma - j), mu is the summation index and t_k are
the singular points listed in \"Points\".  MCFormula[res, key] gives one entry.";
MCTeX::usage = "MCTeX[res, file] writes a LaTeX document with the closed formulas of res to file (and returns the file name).";
MCCheck::usage = "MCCheck[res] verifies the closed forms on a grid of index vectors.  Method -> \"IBP\": every integration-by-parts identity
(maximal cut, sub-sectors = 0) is evaluated with the closed-form coefficients (relative residuals); Method -> \"Laporta\": exact modular
Laporta reduction on the maximal cut, rows compared with the closed forms; Method -> \"Both\".  Options \"Grid\" -> {dotsMax, rankMax},
\"Dimension\" -> 37/11, \"LaportaBox\" -> Automatic.  Returns the largest deviation.";
MCCompare::usage = "MCCompare[res, targets] reduces the index vectors targets (and a random sample, option \"Sample\") with LiteRed2
(R. N. Lee) on the maximal cut -- CutDs on every propagator, run in a separate kernel launched with LaunchKernels[1] -- and compares
the exact LiteRed coefficients with MCReduce[res, ...] at the values of d given by \"Dimensions\".  If LiteRed chose other masters
(same number) the closed forms are recomputed with LiteRed's masters; if the numbers differ (LiteRed also uses symmetry relations)
the check is I(n) = c(n).L(masters), with L the LiteRed reduction of our masters.  LiteRed results are cached in the folder
$MCReferenceDirectory (default: MaxCutClosedForm_LiteRed next to the package).  Returns the largest relative difference.
Options: \"Dimensions\" -> {37/11, 43/13}, \"Sample\" -> {8, {2, 4}} (number of random targets, {dots, rank}), \"Seed\", \"Show\" -> 4,
\"TimeLimit\" -> 1800 (seconds for LiteRed), \"Masters\" -> Automatic | \"Ours\", \"LiteRedLog\" -> False, \"Verbose\" -> True.";
MCCompareSummary::usage = "MCCompareSummary[] prints a table of all MCCompare results of the session and returns them.";
MCStopReference::usage = "MCStopReference[] closes the kernel used by MCCompare for LiteRed.";
$MCLiteRedFile::usage = "$MCLiteRedFile = Automatic (<<LiteRed2`) or the full path of the LiteRed2 main file, e.g. LiteRed2024.m.";
$MCReferenceDirectory::usage = "$MCReferenceDirectory = Automatic (folder MaxCutClosedForm_LiteRed next to the package) or a directory for the
LiteRed cache (subfolder cache) and LiteRed's working files (subfolder work).";
MCSelfTest::usage = "MCSelfTest[] tests the numerical Lauricella F_D (series, analytically continued Euler integral) against direct quadrature and built-in functions.";
MCLauricellaFD::usage = "MCLauricellaFD[a, {b1, ..., bN}, c, {x1, ..., xN}] is the Lauricella function F_D^(N)(a; b1..bN; c; x1..xN)
(N = 1: Hypergeometric2F1, N = 2: AppellF1); it evaluates numerically for inexact arguments (series or analytically continued Euler integral).";
MCResult::usage = "MCResult[...] is returned by MCClosedForm; res[\"key\"] gives a part: \"Masters\", \"NumberOfMasters\", \"Cycles\", \"Selected\",
\"Fibrations\", \"PS\", \"P\", \"Point\", \"Topology\", ...";
MCz::usage = "MCz[i] is the Baikov variable z_i = D_i of the i-th propagator / ISP.";

Begin["`Private`"];

$mcPackageDir = If[StringQ[$InputFileName] && $InputFileName =!= "", DirectoryName[$InputFileName], Directory[]];


(* ::Section:: *)
(*Helpers*)


getD[asc_Association, k_, def_] := If[KeyExistsQ[asc, k], asc[k], def];
SetAttributes[fixedPrec, HoldRest];
(* at least p digits everywhere (no upper cap: Mathematica's special functions raise their internal precision);
   N::precsm only reports requests of lower precision inside built-ins and is harmless *)
fixedPrec[p_, expr_] := Quiet[Block[{$MinPrecision = p}, expr], {N::precsm}];
binomRho[g_, jo_Integer] := Product[(g - i)/(i + 1), {i, 0, jo - 1}];
leadC[p_, t_] := Coefficient[p, t, Exponent[p, t]];
Format[MCz[i_Integer]] := Subscript["z", i];
normProp[{q_, m2_}] := {q, m2};
normProp[q_] := {q, 0};
$mcCounter = 0;
$mcMemo = <||>;          (* per result: memo tables (residue terms, cycle integrals, master matrices) *)
memoGet[id_, tab_String, key_] := Module[{m = getD[$mcMemo, id, <||>]}, getD[getD[m, tab, <||>], key, Missing[]]];
memoSet[id_, tab_String, key_, val_] := Module[{m = getD[$mcMemo, id, <||>], t},
  t = getD[m, tab, <||>]; t[key] = val; m[tab] = t; $mcMemo[id] = m; val];


(* ::Section:: *)
(*MCTopology*)


Options[MCTopology] = {"Name" -> "topology"};
MCTopology::count = "Propagators + ISPs must number L(L+1)/2 + L E = `1` (given: `2`).";
MCTopology::lin = "The momentum `1` is not a linear combination of the loop and external momenta with numerical coefficients.";
MCTopology::kin = "The kinematics does not fix the external scalar product(s) `1`.";
MCTopology::basis = "The propagators and ISPs do not determine all scalar products (they are linearly dependent).";

MCTopology[loops_List, legs_List, props_List, isps_List, kin_List, OptionsPattern[]] := Catch[Module[
  {basis, nL, nE, nB, all, nAll, nProps, coeffs, ext, missing, Gm, zv, Dexpr, unknowns, sol, i1, i2},
  basis = Join[loops, legs]; nL = Length[loops]; nE = Length[legs]; nB = nL + nE;
  all = normProp /@ Join[props, isps]; nAll = Length[all]; nProps = Length[props];
  If[nAll =!= nL (nL + 1)/2 + nL nE, Message[MCTopology::count, nL (nL + 1)/2 + nL nE, nAll]; Throw[$Failed, mcTag]];
  coeffs = Table[Coefficient[all[[i1, 1]], #] & /@ basis, {i1, nAll}];
  Do[If[!VectorQ[coeffs[[i1]], NumericQ] || Expand[all[[i1, 1]] - coeffs[[i1]] . basis] =!= 0,
       Message[MCTopology::lin, all[[i1, 1]]]; Throw[$Failed, mcTag]], {i1, nAll}];
  (* external Gram matrix from the kinematic rules *)
  ext = Table[If[i1 == i2, legs[[i1]]^2, legs[[i1]] legs[[i2]]] /. kin, {i1, nE}, {i2, nE}];
  missing = Union[Flatten[Table[If[FreeQ[ext[[i1, i2]], Alternatives @@ legs], {}, {If[i1 == i2, legs[[i1]]^2, legs[[i1]] legs[[i2]]]}], {i1, nE}, {i2, i1, nE}]]];
  If[missing =!= {}, Message[MCTopology::kin, missing]; Throw[$Failed, mcTag]];
  (* Gram matrix of (loops, legs): loop entries are unknowns spv[i, j] *)
  Gm = Table[If[i1 <= nL || i2 <= nL, spv[Min[i1, i2], Max[i1, i2]], ext[[i1 - nL, i2 - nL]]], {i1, nB}, {i2, nB}];
  zv = Array[MCz, nAll];
  Dexpr = Table[Expand[coeffs[[i1]] . Gm . coeffs[[i1]] - all[[i1, 2]]], {i1, nAll}];
  unknowns = Union[Cases[Gm, _spv, {2}]];
  sol = Quiet[Solve[Thread[zv == Dexpr], unknowns]];
  If[Length[sol] =!= 1 || Length[First[sol]] =!= Length[unknowns], Message[MCTopology::basis]; Throw[$Failed, mcTag]];
  <|"Type" -> "MCTopology", "Name" -> OptionValue["Name"], "Loops" -> loops, "Legs" -> legs, "Momenta" -> all[[All, 1]],
    "Masses" -> all[[All, 2]], "NProps" -> nProps, "N" -> nAll, "L" -> nL, "E" -> nE, "Basis" -> basis, "Coefficients" -> coeffs,
    "GramSP" -> Gm, "SPofZ" -> First[sol], "GramZ" -> Expand[Gm /. First[sol]], "Kinematics" -> kin|>], mcTag];


(* ::Section:: *)
(*Structure of the maximal cut: factors, singular points, Euler characteristic*)


(* monic irreducible factors (over Q) of the polynomials in src = {{tag, poly}, ...} in the variable t, with their multiplicity per tag
   (tags: "A", "D" (Delta), "B" (-B), "P" (P_S, 1 ISP), "T" (the monomial t)) *)
factorTable[src_List, t_] := Module[{tab = <||>, fl, f, m, fm, tag, poly},
  Do[
    tag = s[[1]]; poly = s[[2]];
    If[poly =!= 0 && Exponent[poly, t] > 0,
      fl = FactorList[poly];
      Do[f = fe[[1]]; m = fe[[2]];
        If[Exponent[f, t] > 0,
          fm = Expand[f/leadC[f, t]];
          If[!KeyExistsQ[tab, fm], tab[fm] = <|"A" -> 0, "D" -> 0, "B" -> 0, "P" -> 0, "T" -> 0|>];
          tab[fm] = Append[tab[fm], tag -> tab[fm][tag] + m]],
        {fe, fl}]],
    {s, src}];
  tab];

(* exponents of the point t_k:  |t - t_k|^Eabs (t - t_k)^N  (real chambers),  (t - t_k)^(Eabs + N)  (complex segments) *)
pointExponents[m_Association, "2isp"] := {-m["A"] (\[FormalRho] + 1) + m["D"] (\[FormalRho] + 1/2),
   -m["A"] \[FormalA] + m["D"] \[FormalMu] + m["B"] (\[FormalA] - 2 \[FormalMu]) + m["T"] \[FormalB]};
pointExponents[m_Association, "1isp"] := {m["P"] \[FormalRho], m["T"] \[FormalA]};

(* the points: every root of every factor; real ones first (ascending), then the complex ones; "Label" = position *)
makePoints[facs_Association, t_, kind_String, wp_] := Module[{out = {}, rts, m, v, ex},
  Do[
    m = facs[fm];
    rts = t /. Solve[fm == 0, t, Cubics -> False, Quartics -> False];
    ex = pointExponents[m, kind];
    Do[
      v = N[r, wp + 40];
      AppendTo[out, <|"t" -> r, "v" -> v, "Real" -> (Im[v] == 0), "Factor" -> fm, "Mult" -> m, "Eabs" -> ex[[1]], "N" -> ex[[2]],
                     "Singular" -> (m["A"] + m["D"] + m["P"] > 0)|>],
      {r, rts}],
    {fm, Keys[facs]}];
  out = Join[SortBy[Select[out, #["Real"] &], Re[#["v"]] &], Select[out, !#["Real"] &]];
  Table[Append[out[[i]], "Label" -> i], {i, Length[out]}]];

(* Euler characteristic of {(x, t): P_S != 0} from the fibration over t (deg_x P_S = 2):
   chi = |S| - 1 + sum_{s in S} chi(fibre_s),  S = roots of A Delta,  chi(fibre) = 1 - #roots of P_S(x, t_s) (0 if P_S(., t_s) = 0) *)
euler2[facs_Association, aA_, bB_, cC_, t_] := Module[{nS = 0, tot = 0, k, c, zeroQ},
  zeroQ[p_, f_] := Expand[PolynomialRemainder[p, f, t]] === 0;
  Do[
    If[facs[f]["A"] > 0 || facs[f]["D"] > 0,
      k = Exponent[f, t]; nS += k;
      c = If[facs[f]["A"] > 0 && zeroQ[bB, f], If[zeroQ[cC, f], 0, 1], 0];
      tot += k c],
    {f, Keys[facs]}];
  nS - 1 + tot];

analyse1[PS_, t_, wp_] := Module[{facs, nroots, pts},
  facs = factorTable[{{"P", PS}, {"T", t}}, t];
  nroots = Total[Exponent[#, t] & /@ Select[Keys[facs], facs[#]["P"] > 0 &]];
  If[nroots == 0, Return[<|"Kind" -> 1, "NM" -> 0, "Fibs" -> {}, "Reason" -> "P_S has no root in the ISP"|>]];
  pts = makePoints[facs, t, "1isp", wp];
  <|"Kind" -> 1, "NM" -> nroots - 1, "Fibs" -> {<|"Kind" -> "1isp", "ISP" -> 1, "x" -> None, "t" -> t, "PS" -> PS, "lc" -> leadC[PS, t],
     "Points" -> pts, "Chi" -> 1 - nroots|>}, "Reason" -> If[nroots == 1, "P_S has a single root in the ISP", ""]|>];

analyse2[PS_, zY_List, wp_] := Catch[Module[{degs, fibs = {}, x, t, aA, bB, cC, dD, facs, pts, chis},
  degs = Exponent[PS, #] & /@ zY;
  If[Min[degs] <= 1, Throw[<|"Kind" -> 2, "NM" -> 0, "Fibs" -> {}, "Reason" -> "P_S is at most linear in an ISP (scaleless)"|>, zeroTag]];
  If[Min[degs] > 2, Message[MCClosedForm::nogauss, degs]; Throw[$Failed, mcTag]];
  Do[
    If[degs[[ib]] == 2,
      x = zY[[ib]]; t = zY[[3 - ib]];
      aA = Expand[Coefficient[PS, x, 2]]; bB = Expand[Coefficient[PS, x, 1]/2]; cC = Expand[Coefficient[PS, x, 0]];
      dD = Expand[bB^2 - aA cC];
      If[dD === 0, Throw[<|"Kind" -> 2, "NM" -> 0, "Fibs" -> {}, "Reason" -> "P_S is a perfect square in an ISP (scaleless)"|>, zeroTag]];
      facs = factorTable[{{"A", aA}, {"D", dD}, {"B", -bB}, {"T", t}}, t];
      pts = makePoints[facs, t, "2isp", wp];
      AppendTo[fibs, <|"Kind" -> "2isp", "ISP" -> ib, "x" -> x, "t" -> t, "A" -> aA, "B" -> bB, "C" -> cC, "Delta" -> dD,
         "lcA" -> leadC[aA, t], "lcD" -> leadC[dD, t], "lcB" -> If[bB === 0, 0, leadC[-bB, t]], "BZero" -> (bB === 0),
         "Points" -> pts, "Chi" -> euler2[facs, aA, bB, cC, t]|>]],
    {ib, 2}];
  chis = Union[#["Chi"] & /@ fibs];
  If[Length[chis] > 1, Print["MCClosedForm: warning: the two fibrations give different Euler characteristics ", chis]];
  <|"Kind" -> 2, "NM" -> Max[chis], "Fibs" -> fibs, "Reason" -> ""|>], zeroTag];


(* ::Section:: *)
(*Cycles: real chambers, complex segments (spanning tree), chamber signs, reference points*)


segDist[o_, p_, q_] := Module[{L = q - p, u}, u = Re[(o - p) Conjugate[L]]/Abs[L]^2; u = Min[Max[u, 0], 1]; Abs[o - (p + u L)]];

cycleSigma[cyc_Association, fib_Association] := Module[{pts = fib["Points"], tm, t = fib["t"]},
  If[cyc["Kind"] === "complex", Return[If[fib["Kind"] === "2isp", -1, 1]]];
  tm = Which[cyc["p"] === -Infinity, pts[[cyc["q"]]]["v"] - 1, cyc["q"] === Infinity, pts[[cyc["p"]]]["v"] + 1,
             True, (pts[[cyc["p"]]]["v"] + pts[[cyc["q"]]]["v"])/2];
  If[fib["Kind"] === "2isp", Sign[-(fib["Delta"]/fib["A"]) /. t -> tm], Sign[fib["PS"] /. t -> tm]]];

(* reference point of the Moebius map (p -> 0, q -> 1, r -> infinity) minimising max |x_k| over the singular points *)
chooseRef[cyc_Association, fib_Association] := Module[{pts = fib["Points"], sing, cands, score, vp, vq, qq, vqq, side},
  If[cyc["Kind"] === "complex", Return[<|"Type" -> "none"|>]];
  sing = #["Label"] & /@ Select[pts, #["Singular"] &];
  If[cyc["p"] =!= -Infinity && cyc["q"] =!= Infinity,
    vp = pts[[cyc["p"]]]["v"]; vq = pts[[cyc["q"]]]["v"];
    cands = Join[{Infinity}, #["Label"] & /@ Select[pts, #["Real"] && (#["v"] < vp || #["v"] > vq) &]];
    score[r_] := Max[Join[{0},
      Table[Abs[If[r === Infinity, (vq - vp)/(pts[[k]]["v"] - vp), (vq - vp) (pts[[k]]["v"] - pts[[r]]["v"])/((pts[[k]]["v"] - vp) (vq - pts[[r]]["v"]))]],
        {k, DeleteCases[sing, cyc["p"] | cyc["q"] | r]}],
      If[r === Infinity, {}, {Abs[(vq - vp)/(vq - pts[[r]]["v"])]}]]];
    With[{best = First[SortBy[cands, score]]}, Return[If[best === Infinity, <|"Type" -> "inf"|>, <|"Type" -> "point", "Label" -> best|>]]]];
  (* unbounded chamber: finite endpoint qq, reference on the other side of qq *)
  qq = If[cyc["p"] === -Infinity, cyc["q"], cyc["p"]]; vqq = pts[[qq]]["v"];
  side = If[cyc["p"] === -Infinity, (#["v"] > vqq &), (#["v"] < vqq &)];
  cands = #["Label"] & /@ Select[pts, #["Real"] && side[#] &];
  If[cands === {}, Return[<|"Type" -> "aux", "Offset" -> If[cyc["p"] === -Infinity, 1, -1]|>]];      (* auxiliary reference t_qq +- 1 *)
  score[r_] := Max[Join[{0}, Table[Abs[(pts[[k]]["v"] - pts[[r]]["v"])/(vqq - pts[[r]]["v"])], {k, DeleteCases[sing, qq | r]}]]];
  <|"Type" -> "point", "Label" -> First[SortBy[cands, score]]|>];

makeCycles[fib_Association, fi_Integer] := Module[{pts = fib["Points"], sing, real, cyc = {}, par, find, pairs, ok, i1, la, lb},
  sing = Select[pts, #["Singular"] &]; real = Select[sing, #["Real"] &];
  If[real =!= {},
    AppendTo[cyc, <|"Kind" -> "real", "Fib" -> fi, "p" -> -Infinity, "q" -> real[[1]]["Label"]|>];
    Do[AppendTo[cyc, <|"Kind" -> "real", "Fib" -> fi, "p" -> real[[i1]]["Label"], "q" -> real[[i1 + 1]]["Label"]|>], {i1, Length[real] - 1}];
    AppendTo[cyc, <|"Kind" -> "real", "Fib" -> fi, "p" -> real[[-1]]["Label"], "q" -> Infinity|>]];
  (* complex segments: spanning tree of the singular points (Kruskal, shortest first); real neighbours are already joined *)
  par = Association[Table[s["Label"] -> s["Label"], {s, sing}]];
  find[i_] := Module[{r = i}, While[par[r] =!= r, r = par[r]]; r];
  Do[par[find[real[[i1]]["Label"]]] = find[real[[i1 + 1]]["Label"]], {i1, Length[real] - 1}];
  pairs = SortBy[Select[Subsets[sing, {2}], !(#[[1]]["Real"] && #[[2]]["Real"]) &], Abs[#[[1]]["v"] - #[[2]]["v"]] &];
  Do[
    la = pr[[1]]["Label"]; lb = pr[[2]]["Label"];
    If[find[la] =!= find[lb],
      ok = AllTrue[Select[sing, !MemberQ[{la, lb}, #["Label"]] &], segDist[#["v"], pr[[1]]["v"], pr[[2]]["v"]] > 10^-8 &];
      If[ok, par[find[la]] = find[lb]; AppendTo[cyc, <|"Kind" -> "complex", "Fib" -> fi, "p" -> la, "q" -> lb|>]]],
    {pr, pairs}];
  Table[With[{c = cyc[[i1]]}, Join[c, <|"Sigma" -> cycleSigma[c, fib], "Ref" -> chooseRef[c, fib]|>]], {i1, Length[cyc]}]];


(* ::Section:: *)
(*The closed formula*)


(* int over the cycle of prod_k |t - t_k|^Eabs_k (t - t_k)^N_k  (real chambers)  or  prod_k (t - t_k)^(Eabs_k + N_k)  (complex segments,
   principal powers of fixed numbers), as Gamma functions times ONE Lauricella F_D.  val[k] = value (or display symbol) of t_k. *)
intervalExpr[cyc_Association, pts_List, val_] := Module[
  {labels = Range[Length[pts]], ea, en, ee, einf, p, q, qq, P, Q, QQ, R, ks, xs, bs, y, ref, pref},
  ea[k_] := pts[[k]]["Eabs"]; en[k_] := pts[[k]]["N"]; ee[k_] := ea[k] + en[k];
  einf = -2 - Total[ee /@ labels];
  ref = cyc["Ref"];
  Which[
    cyc["Kind"] === "complex",
      p = cyc["p"]; q = cyc["q"]; P = val[p]; Q = val[q]; ks = DeleteCases[labels, p | q];
      pref = (Q - P)^(1 + ee[p]) (P - Q)^ee[q] Times @@ ((P - val[#])^ee[#] & /@ ks);
      pref Gamma[ee[p] + 1] Gamma[ee[q] + 1]/Gamma[ee[p] + ee[q] + 2] MCLauricellaFD[ee[p] + 1, -ee /@ ks, ee[p] + ee[q] + 2, (Q - P)/(val[#] - P) & /@ ks],
    cyc["p"] =!= -Infinity && cyc["q"] =!= Infinity,
      p = cyc["p"]; q = cyc["q"]; P = val[p]; Q = val[q];
      If[ref["Type"] === "inf",
        ks = DeleteCases[labels, p | q]; y = 0; xs = (Q - P)/(val[#] - P) & /@ ks; bs = -ee /@ ks,
        R = val[ref["Label"]]; ks = DeleteCases[labels, p | q | ref["Label"]]; y = (Q - P)/(Q - R);
        xs = Append[(Q - P) (val[#] - R)/((val[#] - P) (Q - R)) & /@ ks, y]; bs = Append[-ee /@ ks, -einf]];
      pref = ((Q - P) (1 - y))^(1 + ee[p]) Abs[Q - P]^ea[q] (P - Q)^en[q] Times @@ ((Abs[P - val[#]]^ea[#] (P - val[#])^en[#]) & /@ DeleteCases[labels, p | q]);
      pref Gamma[ee[p] + 1] Gamma[ee[q] + 1]/Gamma[ee[p] + ee[q] + 2] MCLauricellaFD[ee[p] + 1, bs, ee[p] + ee[q] + 2, xs],
    True,
      qq = If[cyc["p"] === -Infinity, cyc["q"], cyc["p"]]; QQ = val[qq];
      If[ref["Type"] === "aux", R = QQ + ref["Offset"]; ks = DeleteCases[labels, qq], R = val[ref["Label"]]; ks = DeleteCases[labels, qq | ref["Label"]]];
      pref = Abs[QQ - R]^(1 + Total[ea /@ labels]) (QQ - R)^Total[en /@ labels];
      pref Gamma[einf + 1] Gamma[ee[qq] + 1]/Gamma[einf + ee[qq] + 2] MCLauricellaFD[einf + 1, -ee /@ ks, einf + ee[qq] + 2, (val[#] - R)/(QQ - R) & /@ ks]]];

(* summand of the cycle integral G_C:  2 ISPs: the mu-th term (sum over mu = 0 .. Floor[a/2]);  1 ISP: the whole integral *)
gTerm[cyc_Association, fib_Association, val_] := Module[{lcA, lcD, lcB, blk, cst, bfac, lc},
  If[fib["Kind"] === "1isp",
    lc = fib["lc"];
    Return[If[cyc["Kind"] === "real", Abs[lc]^\[FormalRho], lc^\[FormalRho]] intervalExpr[cyc, fib["Points"], val]]];
  lcA = fib["lcA"]; lcD = fib["lcD"]; lcB = fib["lcB"];
  blk = Binomial[\[FormalA], 2 \[FormalMu]] Gamma[\[FormalMu] + 1/2] Gamma[\[FormalRho] + 1]/Gamma[\[FormalRho] + \[FormalMu] + 3/2];
  bfac = If[fib["BZero"], KroneckerDelta[\[FormalA], 2 \[FormalMu]], lcB^(\[FormalA] - 2 \[FormalMu])];
  cst = If[cyc["Kind"] === "real",
    Abs[lcA]^(-\[FormalRho] - 1) lcA^(-\[FormalA]) Abs[lcD]^(\[FormalRho] + 1/2) lcD^\[FormalMu],
    lcA^(-\[FormalRho] - 1 - \[FormalA]) lcD^(\[FormalRho] + 1/2 + \[FormalMu])];
  blk cst bfac intervalExpr[cyc, fib["Points"], val]];


(* ::Section:: *)
(*Numerical Lauricella F_D*)


(* F_D^(N)(a; b; c; x) = Gamma(c)/(Gamma(a) Gamma(c - a)) int_0^1 u^(a-1) (1-u)^(c-a-1) prod_k (1 - x_k u)^(-b_k) du, continued in a, c, x *)
MCLauricellaFD[a_, {}, c_, {}] := 1;
MCLauricellaFD[a_?InexactNumberQ, bs_List, c_?InexactNumberQ, xs_List] /;
   Length[bs] == Length[xs] && Length[bs] > 0 && VectorQ[bs, NumericQ] && VectorQ[xs, NumericQ] && FreeQ[{bs, xs}, Infinity | ComplexInfinity | Indeterminate] :=
  lauricellaN[a, N[bs, Precision[a]], c, N[xs, Precision[a]]];

terminatingQ[b_, pr_] := Im[b] == 0 && Abs[b - Round[Re[b]]] < 10^(-pr/2) && Round[Re[b]] <= 0;

lauricellaN[a_, bs_List, c_, xs_List] := Module[{pr = Precision[{a, bs, c, xs}], term, other, bs2},
  If[pr === MachinePrecision, pr = $MachinePrecision];
  If[Length[bs] == 1, Return[Hypergeometric2F1[a, bs[[1]], c, xs[[1]]]]];
  term = Select[Range[Length[bs]], terminatingQ[bs[[#]], pr] &];            (* polynomial factors (1 - x u)^n, any x *)
  other = Complement[Range[Length[bs]], term];
  bs2 = bs; Do[bs2[[k]] = Round[Re[bs[[k]]]], {k, term}];                  (* exact integer exponents for the polynomial factors *)
  If[AllTrue[xs[[other]], Abs[#] < 85/100 &] && seriesCost[bs2, xs, pr] < 3000, Return[fdSeries[a, bs2, c, xs, pr]]];
  eulerContinued[a, bs2, c, xs, pr, term] Gamma[c]/(Gamma[a] Gamma[c - a])];

(* order of the nested sums: terminating variables outermost, then increasing |x|; the innermost sum is Hypergeometric2F1 *)
seriesOrder[bs_, xs_, pr_] := SortBy[Range[Length[bs]], {If[terminatingQ[bs[[#]], pr], 0, 1], Abs[xs[[#]]]} &];

(* estimated number of nested series terms (all but the innermost sum), to choose between series and Euler integral *)
seriesCost[bs_, xs_, pr_] := Module[{nt},
  nt = Table[If[terminatingQ[bs[[k]], pr], -Round[Re[bs[[k]]]] + 1, If[Abs[xs[[k]]] < 10^-30, 1, Ceiling[(pr + 5) Log[10]/(-Log[Abs[xs[[k]]]])] + 5]], {k, Length[bs]}];
  Times @@ Most[nt[[seriesOrder[bs, xs, pr]]]]];

fdSeries[a_, {b_}, c_, {x_}, pr_] := Hypergeometric2F1[a, b, c, x];
fdSeries[a_, bs_List, c_, xs_List, pr_] /; Length[bs] >= 2 := Module[{ord, b1, x1, rb, rx, tot = 0, term, coef = 1, m = 0, mmax, eps},
  ord = seriesOrder[bs, xs, pr]; b1 = bs[[ord[[1]]]]; x1 = xs[[ord[[1]]]]; rb = bs[[Rest[ord]]]; rx = xs[[Rest[ord]]];
  mmax = If[terminatingQ[b1, pr], -Round[Re[b1]], Infinity];
  eps = 10^(-pr - 3);
  While[True,
    term = coef fdSeries[a + m, rb, c + m, rx, pr];
    tot += term;
    If[m >= mmax, Break[]];
    If[m > 10 && Abs[term] <= eps Max[Abs[tot], 10^-300], Break[]];
    coef = coef (a + m) (b1 + m) x1/((c + m) (m + 1));
    m++;
    If[m > 20000, Print["MCLauricellaFD: series did not converge (|x| = ", Abs[x1], ")"]; Break[]]];
  tot];

(* Taylor coefficients (orders 0..nt-1) of prod_k (1 - y_k v)^(-b_k) *)
taylorProd[bs_List, ys_List, nt_Integer] := Module[{co = UnitVector[nt, 1], f},
  Do[
    f = FoldList[#1 (bs[[k]] + #2 - 1)/#2 ys[[k]] &, 1, Range[nt - 1]];
    co = Table[Sum[co[[j + 1]] f[[i - j + 1]], {j, 0, i}], {i, 0, nt - 1}],
    {k, Length[bs]}];
  co];

(* {B(x; p + n, q)}_{n = 0..nt-1},  B(x; p, q) = int_0^x u^(p-1) (1-u)^(q-1) du = x^p sum_k (1-q)_k x^k/(k! (p + k))  (continued in p, 0 < x <= 0.45) *)
betaIncList[p_, q_, x_, nt_Integer, pr_] := Module[{w = {1}, wk = 1, k = 0},
  While[k < 20000,
    wk = wk (1 - q + k) x/(k + 1); k++; AppendTo[w, wk];
    If[k > 10 && Abs[wk] < 10^(-pr - 10), Break[]]];
  Table[x^(p + n) Sum[w[[j + 1]]/(p + n + j), {j, 0, Length[w] - 1}], {n, 0, nt - 1}]];

(* E(a; b; c; x) = int_0^1 u^(a-1) (1-u)^(c-a-1) prod_k (1 - x_k u)^(-b_k) du = Gamma(a) Gamma(c-a)/Gamma(c) F_D, continued in a and c:
   near each endpoint (half the disc of convergence) the product is Taylor expanded and integrated term by term with the continued
   incomplete Beta function; the middle piece is a smooth quadrature *)
eulerContinued[a_, bs_List, c_, xs_List, pr_] := eulerContinued[a, bs, c, xs, pr, {}];
eulerContinued[a_, bs_List, c_, xs_List, pr_, skip_List] := Module[{nz, R0, R1, u0, u1, nt, g0, gv, i0, i1, im, u},
  nz = Select[xs[[Complement[Range[Length[xs]], skip]]], Abs[#] > 10^-30 &];     (* polynomial factors (skip) are entire *)
  R0 = Min[Join[{1}, Abs[1/#] & /@ nz]]; R1 = Min[Join[{1}, Abs[1 - 1/#] & /@ nz]];
  u0 = Min[R0/2, 45/100]; u1 = 1 - Min[R1/2, 45/100];
  u0 = SetPrecision[u0, pr]; u1 = SetPrecision[u1, pr];
  nt = Ceiling[(pr + 5) 3.33] + 10;
  g0 = taylorProd[bs, xs, nt];                                                                    (* g(u) around u = 0 *)
  gv = (Times @@ MapThread[(1 - #1)^(-#2) &, {xs, bs}]) taylorProd[bs, -xs/(1 - xs), nt];          (* g(1 - v) around v = 0 *)
  i0 = g0 . betaIncList[a, c - a, u0, nt, pr];
  i1 = gv . betaIncList[c - a, a, 1 - u1, nt, pr];
  (* middle piece: smooth integrand; its parameters get pr + 10 digits so that NIntegrate can work at pr digits *)
  im = Block[{$MinPrecision = 0},
    With[{aa = SetPrecision[a, pr + 10], cc = SetPrecision[c, pr + 10], xx = SetPrecision[xs, pr + 10], bb = SetPrecision[bs, pr + 10],
          v0 = SetPrecision[u0, pr + 10], v1 = SetPrecision[u1, pr + 10]},
      NIntegrate[u^(aa - 1) (1 - u)^(cc - aa - 1) Times @@ MapThread[(1 - #1 u)^(-#2) &, {xx, bb}], {u, v0, (v0 + v1)/2, v1},
        WorkingPrecision -> pr, PrecisionGoal -> pr - 10, MaxRecursion -> 20]]];
  i0 + im + i1];


(* ::Section:: *)
(*Residue truncation, cycle integrals, F_C(n), master selection*)


(* Pi_j = [z_S^(n_S - 1)] (P - P_S)^j for j = 0 .. w (w = total dots): <|j -> {{ISP exponent vector, coefficient}, ...}|> *)
residueTerms[st_Association, nS_List] := Module[{m, vars, nSn, target, w, rules, delta, cur, res, mul},
  m = memoGet[st["ID"], "res", nS]; If[!MissingQ[m], Return[m]];
  vars = Join[st["zS"], st["zY"]]; nSn = st["NProps"]; target = nS - 1; w = Total[target];
  rules = CoefficientRules[st["P"] - st["PS"], vars];
  delta = Select[rules, And @@ Thread[Take[#[[1]], nSn] <= target] &];       (* terms of P - P_S that can still contribute *)
  mul[aa_, bb_] := Module[{acc = <||>, e},
    Do[e = r1[[1]] + r2[[1]];
       If[And @@ Thread[Take[e, nSn] <= target], acc[e] = getD[acc, e, 0] + r1[[2]] r2[[2]]],
       {r1, aa}, {r2, bb}];
    Select[Normal[acc], #[[2]] =!= 0 &]];
  cur = {ConstantArray[0, Length[vars]] -> 1}; res = <||>;
  Do[
    If[jo > 0, cur = mul[cur, delta]];
    With[{tt = Select[cur, Take[#[[1]], nSn] == target &]}, If[tt =!= {}, res[jo] = ({Drop[#[[1]], nSn], #[[2]]} & /@ tt)]],
    {jo, 0, w}];
  memoSet[st["ID"], "res", nS, res]];

(* numerical cycle integral G_C(k; rho) of the ISP monomial z_Y^k (k >= 0), rho exact; memoised *)
gNum[st_Association, ci_Integer, kv_List, rho_] := Module[{key = {ci, kv, rho}, m, cyc, fib, a, b, wp2, mus, val, t0},
  m = memoGet[st["ID"], "G", key]; If[!MissingQ[m], Return[m]];
  cyc = st["Cycles"][[ci]]; fib = st["Fibrations"][[cyc["Fib"]]];
  wp2 = st["WorkingPrecision"] + 15;
  t0 = AbsoluteTime[];
  val = fixedPrec[wp2,
    If[fib["Kind"] === "1isp",
      N[cyc["Template"] /. {\[FormalA] -> kv[[1]], \[FormalRho] -> SetPrecision[rho, wp2]}, wp2],
      a = kv[[fib["ISP"]]]; b = kv[[3 - fib["ISP"]]];
      mus = If[fib["BZero"], If[EvenQ[a], {a/2}, {}], Range[0, Floor[a/2]]];
      Total[Table[N[cyc["Template"] /. {\[FormalA] -> a, \[FormalB] -> b, \[FormalMu] -> mu, \[FormalRho] -> SetPrecision[rho, wp2]}, wp2], {mu, mus}]]]];
  If[AbsoluteTime[] - t0 > 10, Print["   (slow cycle integral: cycle ", ci, ", ISP powers ", kv, ", rho = ", rho, ": ", NumberForm[AbsoluteTime[] - t0, {6, 1}], " s)"]];
  memoSet[st["ID"], "G", key, val]];

(* F_C(n) = sum_j Binomial(gamma, j) sigma_C^j sum_alpha [Pi_j]_alpha G_C(k + alpha; gamma - j) *)
fNum[st_Association, ci_Integer, n_List, gam_] := Module[{nP = st["NProps"], kv, terms, sig, tot = 0},
  kv = -Drop[n, nP]; terms = residueTerms[st, Take[n, nP]];
  sig = st["Cycles"][[ci]]["Sigma"];
  Do[tot += binomRho[gam, jo] sig^jo Total[(#[[2]] gNum[st, ci, kv + #[[1]], gam - jo]) & /@ terms[jo]], {jo, Keys[terms]}];
  tot];

(* 0 ISPs: F(n) = sum_j Binomial(gamma, j) P_S^(-j) Pi_j  (exact) *)
f0[st_Association, n_List, gam_] := Module[{terms = residueTerms[st, Take[n, st["NProps"]]], ps = st["PS"]},
  Total[Table[binomRho[gam, jo] ps^(-jo) Total[#[[2]] & /@ terms[jo]], {jo, Keys[terms]}]]];

(* cycles first (with a probe set of low-degree candidates), then masters; columns computed lazily *)
selectMasters[st_Association, mopt_, dTest_, verbose_] := Module[
  {nM = st["NumberOfMasters"], nP = st["NProps"], nY = st["NISP"], wp = st["WorkingPrecision"], cyc = st["Cycles"], top, cands, gam,
   cols = <||>, bad = {}, col, rankOf, order, nprobe, probe, sel = {}, chosen = {}},
  top = ConstantArray[1, nP];
  If[nM == 0, Return[{{}, {}}]];
  If[mopt =!= Automatic && Length[mopt] =!= nM,
    Message[MCClosedForm::masters, "expected " <> ToString[nM] <> " master(s), given " <> ToString[Length[mopt]]]; Throw[$Failed, mcTag]];
  If[st["Kind"] === 0, Return[{If[mopt === Automatic, {top}, mopt], {}}]];
  gam = (dTest - st["L"] - st["E"] - 1)/2;
  cands = If[mopt === Automatic, Join[top, -#] & /@ Flatten[Table[Select[Tuples[Range[0, deg], nY], Total[#] == deg &], {deg, 0, 8}], 1], mopt];
  col[ci_, j_] := Module[{v},
    If[MemberQ[bad, ci], Return[$Failed]];
    If[KeyExistsQ[cols, {ci, j}], Return[cols[{ci, j}]]];
    v = Quiet[fNum[st, ci, cands[[j]], gam]];
    If[!NumericQ[v] || !FreeQ[v, Overflow | Underflow | Indeterminate | DirectedInfinity],
      AppendTo[bad, ci]; If[verbose, Print["   cycle ", ci, " skipped (pole of the cycle integral at d = ", dTest, ")"]]; Return[$Failed]];
    cols[{ci, j}] = v];
  rankOf[rows_List, cycs_List] := Module[{mat},
    If[rows === {} || cycs === {}, Return[0]];
    mat = Table[cols[{ci, j}], {j, rows}, {ci, cycs}];
    mat = If[Max[Abs[#]] == 0, #, #/Max[Abs[#]]] & /@ mat;
    Length[SingularValueList[mat, Tolerance -> 10^(-wp/2)]]];
  order = SortBy[Range[Length[cyc]], {cyc[[#]]["Fib"], Boole[cyc[[#]]["Kind"] === "complex"], Boole[cyc[[#]]["p"] === -Infinity || cyc[[#]]["q"] === Infinity]} &];
  nprobe = If[mopt =!= Automatic, nM, Min[Length[cands], nM + 1]];
  While[True,
    probe = Range[nprobe]; sel = {};
    Do[
      If[!MemberQ[bad, ci] && AllTrue[probe, col[ci, #] =!= $Failed &],
        If[rankOf[probe, Append[sel, ci]] > Length[sel], AppendTo[sel, ci]]];
      If[Length[sel] == nM, Break[]],
      {ci, order}];
    If[Length[sel] == nM || nprobe >= Length[cands] || mopt =!= Automatic, Break[]];
    nprobe = Min[Length[cands], nprobe + Max[2, nY + 1]]];
  If[Length[sel] < nM,
    Message[MCClosedForm::masters, "only " <> ToString[Length[sel]] <> " independent cycles found for " <> ToString[nM] <> " masters"]; Throw[$Failed, mcTag]];
  Do[
    If[AllTrue[sel, col[#, j] =!= $Failed &],
      If[rankOf[Append[chosen, j], sel] > Length[chosen], AppendTo[chosen, j],
        If[mopt =!= Automatic, Message[MCClosedForm::masters, "the given masters are not independent on the maximal cut"]; Throw[$Failed, mcTag]]]];
    If[Length[chosen] == nM, Break[]],
    {j, Length[cands]}];
  If[Length[chosen] < nM, Message[MCClosedForm::masters, "master selection failed"]; Throw[$Failed, mcTag]];
  {cands[[chosen]], sel}];


(* ::Section:: *)
(*MCClosedForm*)


Options[MCClosedForm] = {"Masters" -> Automatic, "Dimension" -> Automatic, "TestDimension" -> 25/7, WorkingPrecision -> 30, "Verbose" -> True};
MCClosedForm::topo = "The first argument is not an MCTopology.";
MCClosedForm::point = "The kinematic point does not fix the symbol(s) `1`.";
MCClosedForm::isp = "`1` irreducible scalar products: only up to 2 are supported.";
MCClosedForm::nogauss = "P_S has degrees `1` in the two ISPs: no quadratic (Gaussian) block, not supported.";
MCClosedForm::masters = "`1`";

MCClosedForm[topo_Association, point_List, OptionsPattern[]] := Catch[Module[
  {verbose, wp, dsym, kinPt, nAll, nP, zv, zS, zY, gz, free, Pfull, PS, nY, data, fibs, cycles, id, st, ms, i1, t0},
  If[getD[topo, "Type", None] =!= "MCTopology", Message[MCClosedForm::topo]; Throw[$Failed, mcTag]];
  verbose = TrueQ[OptionValue["Verbose"]]; wp = OptionValue[WorkingPrecision];
  dsym = OptionValue["Dimension"]; If[dsym === Automatic, dsym = Symbol["d"]];
  kinPt = point;
  If[!FreeQ[kinPt[[All, 2]], _Real | _Complex], kinPt = Rationalize[kinPt, 0]; If[verbose, Print["MCClosedForm: kinematic point rationalised: ", kinPt]]];
  nAll = topo["N"]; nP = topo["NProps"];
  zv = Array[MCz, nAll]; zS = Take[zv, nP]; zY = Drop[zv, nP];
  gz = topo["GramZ"] /. kinPt;
  free = Complement[Variables[gz], zv];
  If[free =!= {}, Message[MCClosedForm::point, free]; Throw[$Failed, mcTag]];
  Pfull = Expand[Det[gz]];
  PS = Expand[Pfull /. Thread[zS -> 0]];
  nY = Length[zY];
  If[nY > 2, Message[MCClosedForm::isp, nY]; Throw[$Failed, mcTag]];
  data = Which[
    PS === 0, <|"Kind" -> nY, "NM" -> 0, "Fibs" -> {}, "Reason" -> "P_S = 0"|>,
    nY == 0, <|"Kind" -> 0, "NM" -> 1, "Fibs" -> {}, "Reason" -> ""|>,
    nY == 1, analyse1[PS, zY[[1]], wp],
    True, analyse2[PS, zY, wp]];
  fibs = data["Fibs"];
  cycles = Join @@ Table[makeCycles[fibs[[fi]], fi], {fi, Length[fibs]}];
  cycles = Table[With[{c = cycles[[i1]], fb = fibs[[cycles[[i1]]["Fib"]]]},
      Join[c, <|"Template" -> gTerm[c, fb, fb["Points"][[#]]["t"] &], "Display" -> gTerm[c, fb, Subscript[\[FormalT], #] &]|>]],
    {i1, Length[cycles]}];
  id = ++$mcCounter;
  st = <|"ID" -> id, "Name" -> topo["Name"], "Topology" -> topo, "Point" -> kinPt, "Dimension" -> dsym, "WorkingPrecision" -> wp,
    "L" -> topo["L"], "E" -> topo["E"], "N" -> nAll, "NProps" -> nP, "NISP" -> nY, "zS" -> zS, "zY" -> zY, "P" -> Pfull, "PS" -> PS,
    "Kind" -> data["Kind"], "NumberOfMasters" -> data["NM"], "Reason" -> data["Reason"], "Fibrations" -> fibs, "Cycles" -> cycles|>;
  If[verbose,
    Print["MCClosedForm: ", topo["Name"], ":  L = ", topo["L"], ", E = ", topo["E"], ", ", nP, " propagators, ", nY, " ISP(s)"];
    Print["   P_S = ", Factor[PS]];
    Print["   number of master integrals (Euler characteristic) = ", data["NM"], If[data["Reason"] =!= "", "   (" <> data["Reason"] <> ")", ""]];
    Do[Print["   fibration ", fi, ": ", If[fibs[[fi]]["Kind"] === "2isp", Row[{"Gaussian block in x = ", fibs[[fi]]["x"], ", endgame variable t = ", fibs[[fi]]["t"]}],
          Row[{"variable t = ", fibs[[fi]]["t"]}]], ";  ", Count[fibs[[fi]]["Points"], p_ /; p["Singular"]], " singular points (",
          Count[fibs[[fi]]["Points"], p_ /; p["Singular"] && p["Real"]], " real), ",
          Count[cycles, c_ /; c["Fib"] == fi && c["Kind"] === "real"], " real chambers, ", Count[cycles, c_ /; c["Fib"] == fi && c["Kind"] === "complex"], " complex segments"],
      {fi, Length[fibs]}]];
  t0 = AbsoluteTime[];
  ms = selectMasters[st, OptionValue["Masters"], OptionValue["TestDimension"], verbose];
  st = Join[st, <|"Masters" -> ms[[1]], "Selected" -> ms[[2]]|>];
  If[verbose,
    Print["   masters: ", ms[[1]], "   (cycles and masters selected in ", NumberForm[AbsoluteTime[] - t0, {6, 2}], " s)"];
    If[ms[[2]] =!= {}, Print["   cycles used: ", Table[cycleLabel[st, ci], {ci, ms[[2]]}]]]];
  MCResult[st]], mcTag];

cycleLabel[st_Association, ci_Integer] := Module[{c = st["Cycles"][[ci]], lab},
  lab[x_] := If[x === -Infinity || x === Infinity, x, Subscript["t", x]];
  Row[{ci, ": ", If[c["Kind"] === "real", "[", "segment "], lab[c["p"]], ", ", lab[c["q"]], If[c["Kind"] === "real", "]", ""],
       " (fibration ", c["Fib"], ", \[Sigma] = ", c["Sigma"], ")"}]];


(* ::Section:: *)
(*MCResult, MCReduce, MCMasters*)


MCResult[st_Association][key_String] := getD[st, key, Missing["KeyAbsent", key]];
MCResult /: MakeBoxes[MCResult[st_Association], fmt_] := ToBoxes[Row[{"MCResult[", st["Name"], ", ", st["NumberOfMasters"],
   If[st["NumberOfMasters"] == 1, " master", " masters"], ", ", st["NISP"], " ISP(s)]"}], fmt];
MCMasters[MCResult[st_Association]] := st["Masters"];

MCReduce::idx = "Index vector `1`: expected `2` integers, the ISP entries <= 0.";
MCReduce::dim = "The dimension `1` is not a number.";
MCReduce[MCResult[st_Association], n_?(VectorQ[#, IntegerQ] &), dval_] := Module[{r = MCReduce[MCResult[st], {n}, dval]}, If[r === $Failed, $Failed, First[r]]];
MCReduce[MCResult[st_Association], ns_?(MatrixQ[#, IntegerQ] &), dval_] := Catch[Module[{dex},
  If[!NumericQ[dval], Message[MCReduce::dim, dval]; Throw[$Failed, mcTag]];
  dex = If[ExactNumberQ[dval], dval, Rationalize[dval, 0]];
  Do[If[Length[n] =!= st["N"] || AnyTrue[Drop[n, st["NProps"]], # > 0 &], Message[MCReduce::idx, n, st["N"]]; Throw[$Failed, mcTag]], {n, ns}];
  mcRows[st, ns, dex]], mcTag];

(* rows A(n) = F(n) . M^(-1) at the exact dimension dex *)
mcRows[st_Association, ns_List, dex_] := Module[{nM = st["NumberOfMasters"], nP = st["NProps"], gam, minv, sel = st["Selected"],
   masters = st["Masters"], wp = st["WorkingPrecision"], wp2, f0m},
  If[nM == 0, Return[Table[{}, {Length[ns]}]]];
  wp2 = wp + 15;
  gam = (dex - st["L"] - st["E"] - 1)/2;
  If[st["Kind"] === 0,
    f0m = f0[st, masters[[1]], gam];
    Return[Table[If[AnyTrue[Take[n, nP], # <= 0 &], {0}, {Together[f0[st, n, gam]/f0m]}], {n, ns}]]];
  minv = memoGet[st["ID"], "Minv", dex];
  If[MissingQ[minv],
    minv = fixedPrec[wp2, Inverse[Table[fNum[st, ci, m, gam], {m, masters}, {ci, sel}]]];
    memoSet[st["ID"], "Minv", dex, minv]];
  Table[If[AnyTrue[Take[n, nP], # <= 0 &], ConstantArray[0, nM],
      realPart[SetPrecision[fixedPrec[wp2, Table[fNum[st, ci, n, gam], {ci, sel}] . minv], wp], wp]], {n, ns}]];

(* real kinematics and real d give real coefficients: drop imaginary rounding noise (complex points / complex segments) *)
realPart[row_List, wp_] := If[Max[Abs[Im[row]]] <= 10^(-wp/2) Max[Join[{10^-300}, Abs[row]]], Re[row], row];


(* ::Section:: *)
(*MCInfo, MCFormula, MCTeX*)


expTable[fb_Association] := Grid[Prepend[
   Table[{Subscript["t", p["Label"]], p["t"], N[p["v"], 8], p["Factor"], If[p["Singular"], "singular", "numerator only"], p["Eabs"], p["N"]}, {p, fb["Points"]}],
   {"point", "value", "\[TildeTilde]", "root of", "", "|t-t_k| exponent", "(t-t_k) exponent"}], Frame -> All, Alignment -> Left];

MCInfo[MCResult[st_Association]] := Module[{fibs = st["Fibrations"], cyc = st["Cycles"]},
  Print[Style[st["Name"], Bold], ":  L = ", st["L"], ", E = ", st["E"], ", ", st["NProps"], " propagators, ", st["NISP"], " ISP(s);  point ", st["Point"]];
  Print["\[Gamma] = (d - L - E - 1)/2 = ", (st["Dimension"] - st["L"] - st["E"] - 1)/2, ";   P_S = ", Factor[st["PS"]]];
  Print["number of master integrals (Euler characteristic of {P_S \[NotEqual] 0}) = ", st["NumberOfMasters"],
        If[st["Reason"] =!= "", "  (" <> st["Reason"] <> ")", ""], ";   masters: ", st["Masters"]];
  Do[With[{fb = fibs[[fi]]},
      If[fb["Kind"] === "2isp",
        Print["fibration ", fi, ":  P_S = A x^2 + 2 B x + C,  x = ", fb["x"], ",  t = ", fb["t"], ";   A = ", Factor[fb["A"]], ",  B = ", Factor[fb["B"]],
              ",  \[CapitalDelta] = B^2 - A C = ", Factor[fb["Delta"]]],
        Print["fibration ", fi, ":  t = ", fb["t"]]];
      Print[expTable[fb]]],
    {fi, Length[fibs]}];
  If[cyc =!= {},
    Print[Grid[Prepend[Table[{ci, cyc[[ci]]["Fib"], cyc[[ci]]["Kind"], cyc[[ci]]["p"], cyc[[ci]]["q"], cyc[[ci]]["Sigma"],
         Switch[cyc[[ci]]["Ref"]["Type"], "inf", "\[Infinity]", "point", Subscript["t", cyc[[ci]]["Ref"]["Label"]], "aux", "endpoint " <> If[cyc[[ci]]["Ref"]["Offset"] > 0, "+ 1", "- 1"], _, "-"],
         If[MemberQ[st["Selected"], ci], "used", ""]}, {ci, Length[cyc]}],
       {"cycle", "fibration", "kind", "from t_p", "to t_q", "\[Sigma]", "reference", ""}], Frame -> All]]];
  ];

MCFormula[MCResult[st_Association]] := Module[{fibs = st["Fibrations"], cyc = st["Cycles"], gam, gdisp},
  gam = (st["Dimension"] - st["L"] - st["E"] - 1)/2;
  gdisp = Table[With[{c = cyc[[ci]], fb = fibs[[cyc[[ci]]["Fib"]]]},
      If[fb["Kind"] === "2isp",
        With[{lhs = Subscript["G", ci][\[FormalA], \[FormalB], \[FormalRho]], rhs = Inactive[Sum][c["Display"], {\[FormalMu], 0, Floor[\[FormalA]/2]}]}, HoldForm[lhs == rhs]],
        With[{lhs = Subscript["G", ci][\[FormalA], \[FormalRho]], rhs = c["Display"]}, HoldForm[lhs == rhs]]]],
    {ci, st["Selected"]}];
  <|"gamma" -> gam, "PS" -> Factor[st["PS"]],
    "Fibrations" -> Table[If[fb["Kind"] === "2isp", <|"x" -> fb["x"], "t" -> fb["t"], "A" -> Factor[fb["A"]], "B" -> Factor[fb["B"]], "C" -> Factor[fb["C"]],
        "Delta" -> Factor[fb["Delta"]]|>, <|"t" -> fb["t"], "PS" -> Factor[fb["PS"]]|>], {fb, fibs}],
    "Points" -> Table[Table[Subscript[\[FormalT], p["Label"]] -> p["t"], {p, fb["Points"]}], {fb, fibs}],
    "Exponents" -> Table[Table[Subscript[\[FormalT], p["Label"]] -> {p["Eabs"], p["N"]}, {p, fb["Points"]}], {fb, fibs}],
    "G" -> gdisp,
    "Cycles" -> Table[cycleLabel[st, ci], {ci, st["Selected"]}],
    "F" -> "F_C(n) = Sum_{j=0}^{w} Binomial(gamma, j) sigma_C^j Sum_alpha [Pi_j(n_S)]_alpha G_C(k + alpha; gamma - j),  w = sum(n_S - 1),  k = -n_Y (ISP powers; in G: a = power of x, b = power of t),  Pi_j(n_S) = [z_S^(n_S - 1)] (P - P_S)^j  (MCFormula[res, \"Pi\"][n_S])",
    "A" -> "A(n) = F(n) . M^(-1),  M_(i,C) = F_C(master_i),  C over the cycles used",
    "Pi" -> (residueTerms[st, #] &),
    "Masters" -> st["Masters"]|>];
MCFormula[r_MCResult, key_String] := MCFormula[r][key];

(* LaTeX *)
(* display symbols a, b, rho, mu, t, z, F, D: created once here (their context is not on $ContextPath, so no General::shdw) *)
$texSym = AssociationThread[{"a", "b", "rho", "mu", "t", "z", "F", "D"},
   Symbol["MCTeXsym`" <> #] & /@ {"a", "b", "\[Rho]", "\[Mu]", "t", "z", "F", "D"}];
texRules[] := {\[FormalA] -> $texSym["a"], \[FormalB] -> $texSym["b"], \[FormalRho] -> $texSym["rho"], \[FormalMu] -> $texSym["mu"],
   \[FormalT] -> $texSym["t"], MCz[i_] :> Subscript[$texSym["z"], i], MCLauricellaFD[x__] :> Subscript[$texSym["F"], $texSym["D"]][x]};
texOf[expr_] := Quiet[Block[{$ContextPath = Join[{"MCTeXsym`"}, $ContextPath]}, ToString[TeXForm[expr /. texRules[]]]], {General::shdw}];

MCTeX[MCResult[st_Association], file_String] := Module[{fibs = st["Fibrations"], cyc = st["Cycles"], out = {}, add, gam},
  add[s__] := AppendTo[out, StringJoin[s]];
  gam = (st["Dimension"] - st["L"] - st["E"] - 1)/2;
  add["\\documentclass[10pt]{article}\n\\usepackage{amsmath,amssymb}\n\\usepackage[margin=1.5cm,landscape]{geometry}\n\\allowdisplaybreaks\n\\begin{document}\n"];
  add["\\section*{Maximal-cut closed forms: ", StringReplace[ToString[st["Name"]], {"_" -> "\\_", "&" -> "\\&"}], "}\n"];
  add["Loops $L=", ToString[st["L"]], "$, independent external momenta $E=", ToString[st["E"]], "$, ", ToString[st["NProps"]], " propagators, ",
      ToString[st["NISP"]], " ISP(s); kinematic point $", texOf[st["Point"]], "$.\n\n"];
  add["\\[ \\gamma = \\frac{d-L-E-1}{2} = ", texOf[gam], ", \\qquad P_S = ", texOf[Factor[st["PS"]]], " \\]\n"];
  add["Number of master integrals (Euler characteristic of $\\{P_S\\neq 0\\}$): ", ToString[st["NumberOfMasters"]], ". Masters (index vectors): $",
      texOf[st["Masters"]], "$.\n\n"];
  add["Reduction: $I(n)=\\sum_i A_i(n)\\, I(\\text{master}_i)$ on the maximal cut, $A(n)=F(n)\\,M^{-1}$, $M_{iC}=F_C(\\text{master}_i)$,\n",
      "\\[ F_C(n)=\\sum_{j=0}^{w}\\binom{\\gamma}{j}\\sigma_C^{\\,j}\\sum_{\\alpha}\\big[\\Pi_j(n_S)\\big]_{\\alpha}\\,G_C(k+\\alpha;\\gamma-j),\\qquad ",
      "\\Pi_j(n_S)=\\big[z_S^{\\,n_S-1}\\big](P-P_S)^j,\\quad w=\\textstyle\\sum_{s}(n_s-1),\\quad k=-n_Y. \\]\n"];
  Do[With[{fb = fibs[[fi]]},
      add["\\subsection*{Fibration ", ToString[fi], "}\n"];
      If[fb["Kind"] === "2isp",
        add["Gaussian block in $x=", texOf[fb["x"]], "$, endgame variable $t=", texOf[fb["t"]], "$: $P_S=A x^2+2Bx+C$,\n\\[ A=", texOf[Factor[fb["A"]]],
            ",\\quad B=", texOf[Factor[fb["B"]]], ",\\quad \\Delta=B^2-AC=", texOf[Factor[fb["Delta"]]], ". \\]\n",
            "In $G_C$: $a$ = power of $x$, $b$ = power of $t$; $\\rho=\\gamma-j$.\n"],
        add["Variable $t=", texOf[fb["t"]], "$, $P_S=", texOf[Factor[fb["PS"]]], "$; in $G_C$: $a$ = power of $t$, $\\rho=\\gamma-j$.\n"]];
      add["\\begin{align*}\n"];
      Do[add["t_{", ToString[p["Label"]], "} &= ", If[Exponent[p["Factor"], fb["t"]] <= 2, texOf[p["t"]], "\\text{root of } " <> texOf[p["Factor"]]],
             " \\approx ", texOf[N[p["v"], 10]], " && |t-t_{", ToString[p["Label"]], "}|^{", texOf[p["Eabs"]], "}\\,(t-t_{", ToString[p["Label"]], "})^{",
             texOf[p["N"]], "}\\\\\n"], {p, fb["Points"]}];
      add["\\end{align*}\n"]],
    {fi, Length[fibs]}];
  If[st["Selected"] =!= {},
    add["\\subsection*{Cycle integrals (closed form in the ISP powers)}\n",
        "Real chambers carry $|t-t_k|^{e}(t-t_k)^{N}$, complex segments $(t-t_k)^{e+N}$ (principal powers); $F_D$ is the Lauricella function.\n"]];
  Do[With[{c = cyc[[ci]], fb = fibs[[cyc[[ci]]["Fib"]]]},
      add["\\paragraph{Cycle ", ToString[ci], "} ", If[c["Kind"] === "real", "real chamber", "complex segment"], " from ",
          If[MemberQ[{-Infinity, Infinity}, c["p"]], "$" <> texOf[c["p"]] <> "$", "$t_{" <> ToString[c["p"]] <> "}$"], " to ",
          If[MemberQ[{-Infinity, Infinity}, c["q"]], "$" <> texOf[c["q"]] <> "$", "$t_{" <> ToString[c["q"]] <> "}$"],
          " (fibration ", ToString[c["Fib"]], ", $\\sigma_C=", ToString[c["Sigma"]], "$):\n"];
      If[fb["Kind"] === "2isp",
        add["\\begin{align*}\nG_{", ToString[ci], "}(a,b;\\rho) &= \\sum_{\\mu=0}^{\\lfloor a/2\\rfloor} ", texOf[c["Display"]], "\n\\end{align*}\n"],
        add["\\begin{align*}\nG_{", ToString[ci], "}(a;\\rho) &= ", texOf[c["Display"]], "\n\\end{align*}\n"]]],
    {ci, st["Selected"]}];
  add["\\end{document}\n"];
  Export[file, StringJoin[out], "Text"];
  file];


(* ::Section:: *)
(*MCCheck: IBP identities with the closed forms; exact modular Laporta reduction on the maximal cut*)


(* IBP relations d/dl_a . v (v = loop or external momentum) as lists {shift -> coefficient(nIdx[i], d)} at the kinematic point *)
ibpRelations[st_Association, dval_] := Module[{topo = st["Topology"], nL, nB, nAll, coeffs, Gm, sol, zv, rels = {}, terms, vqz, c0, kap, sh, unitv},
  nL = topo["L"]; nB = topo["L"] + topo["E"]; nAll = topo["N"]; coeffs = topo["Coefficients"]; Gm = topo["GramSP"]; sol = topo["SPofZ"];
  zv = Array[MCz, nAll];
  Do[
    terms = <||>;
    If[v == a, terms[ConstantArray[0, nAll]] = dval];
    Do[
      If[coeffs[[i, a]] =!= 0,
        unitv = UnitVector[nB, v];
        vqz = Expand[(unitv . Gm . coeffs[[i]]) /. sol /. st["Point"]];
        c0 = vqz /. Thread[zv -> 0];
        If[c0 =!= 0, sh = UnitVector[nAll, i]; terms[sh] = getD[terms, sh, 0] - 2 coeffs[[i, a]] c0 nIdx[i]];
        Do[kap = Coefficient[vqz, zv[[jz]]];
           If[kap =!= 0, sh = UnitVector[nAll, i] - UnitVector[nAll, jz]; terms[sh] = getD[terms, sh, 0] - 2 coeffs[[i, a]] kap nIdx[i]],
           {jz, nAll}]],
      {i, nAll}];
    AppendTo[rels, Select[Normal[terms], Expand[#[[2]]] =!= 0 &]],
    {a, nL}, {v, nB}];
  rels];

gridSeeds[nP_, nY_, dmax_, rmax_] := Flatten[Table[Join[1 + dv, -kv],
   {dv, Select[Tuples[Range[0, dmax], nP], Total[#] <= dmax &]}, {kv, Select[Tuples[Range[0, rmax], nY], Total[#] <= rmax &]}], 1];

ibpResidualCheck[st_Association, rels_, dmax_, rmax_, dval_, verbose_] := Module[{nP = st["NProps"], nAll = st["N"], nM = st["NumberOfMasters"],
   seeds, nv, needed, rows, worst = 0, cnt = 0, vec, scale, coef, mm, t0 = AbsoluteTime[]},
  seeds = gridSeeds[nP, nAll - nP, dmax, rmax]; nv = Array[nIdx, nAll];
  If[nM == 0, If[verbose, Print["MCCheck (IBP): no master on the maximal cut -- nothing to evaluate (use Method -> \"Laporta\" to confirm)"]]; Return[0]];
  needed = Reap[Do[mm = s + r[[1]]; If[AllTrue[Take[mm, nP], # >= 1 &] && AllTrue[Drop[mm, nP], # <= 0 &], Sow[mm]], {s, seeds}, {rel, rels}, {r, rel}]][[2]];
  needed = If[needed === {}, {}, Union[First[needed]]];
  rows = AssociationThread[needed, mcRows[st, needed, dval]];
  Do[
    vec = ConstantArray[0, nM]; scale = 0;
    Do[mm = s + r[[1]];
       If[AllTrue[Take[mm, nP], # >= 1 &],
         coef = r[[2]] /. Thread[nv -> s];
         If[coef =!= 0,
           If[AnyTrue[Drop[mm, nP], # > 0 &], Print["MCCheck: unexpected positive ISP index ", mm]; Continue[]];
           vec += coef rows[mm]; scale += Abs[coef] Norm[rows[mm]]]],
       {r, rel}];
    If[scale > 0, worst = Max[worst, Norm[vec]/scale]; cnt++],
    {s, seeds}, {rel, rels}];
  If[verbose, Print["MCCheck (IBP): ", cnt, " IBP identities on ", Length[seeds], " seeds (dots <= ", dmax, ", rank <= ", rmax, "), d = ", dval,
     ", ", Length[needed], " integrals from the closed forms:  max relative residual = ", N[worst, 3], "   (", Round[AbsoluteTime[] - t0, 0.1], " s)"]];
  worst];

ratrec[a_Integer, m_Integer] := Module[{r0 = m, r1 = Mod[a, m], s0 = 0, s1 = 1, q, bound = Floor[Sqrt[m/2]]},
  While[r1 > bound, q = Quotient[r0, r1]; {r0, r1} = {r1, r0 - q r1}; {s0, s1} = {s1, s0 - q s1}];
  If[s1 == 0 || Abs[s1] > bound, $Failed, r1/s1]];
modQ[x_, p_] := Mod[Numerator[x] PowerMod[Denominator[x], -1, p], p];

laportaCheck[st_Association, rels_, dmax_, rmax_, dval_, box_, verbose_] := Module[
  {nP = st["NProps"], nAll = st["N"], nM = st["NumberOfMasters"], masters = st["Masters"], nY, dd, rr, okbox, seeds, nv, eqsQ, targets, keyF, modRows,
   primes = {}, resid = {}, p = 2^61, done = False, exact = <||>, ok, row, red, closed, worst = 0, t0 = AbsoluteTime[]},
  nY = nAll - nP;
  {dd, rr} = If[box === Automatic, {dmax + 1, rmax + 2}, box];
  okbox[m_, dD_, rR_] := Which[AnyTrue[Take[m, nP], # <= 0 &], "zero", AnyTrue[Drop[m, nP], # > 0 &], "out",
     Total[Take[m, nP] - 1] > dD || -Total[Drop[m, nP]] > rR, "out", True, "var"];
  seeds = gridSeeds[nP, nY, dd, rr]; nv = Array[nIdx, nAll];
  eqsQ = Reap[Do[
      Module[{e = <||>, good = True, c, mm, stt},
        Do[c = r[[2]] /. Thread[nv -> s];
           If[c =!= 0,
             mm = s + r[[1]]; stt = okbox[mm, dd + 1, rr + 1];
             Which[stt === "out", good = False; Break[], stt === "var", e[mm] = getD[e, mm, 0] + c]],
           {r, rel}];
        If[good, e = Select[e, # =!= 0 &]; If[Length[e] > 0, Sow[e]]]],
      {s, seeds}, {rel, rels}]][[2]];
  eqsQ = If[eqsQ === {}, {}, First[eqsQ]];
  targets = Select[gridSeeds[nP, nY, dmax, rmax], !MemberQ[masters, #] &];
  keyF[m_] := If[MemberQ[masters, m], Join[{0, 0, 0}, -m], Join[{1, Total[Take[m, nP] - 1], -Total[Drop[m, nP]]}, -m]];
  modRows[pp_] := Block[{$RecursionLimit = 100000}, Module[{rules = <||>, eqs, e, hit, c, piv, inv, full, res = <||>},
    eqs = Select[Select[Association[KeyValueMap[#1 -> modQ[#2, pp] &, #]], # != 0 &] & /@ eqsQ, Length[#] > 0 &];
    eqs = SortBy[eqs, keyF[Last[SortBy[Keys[#], keyF]]] &];
    Do[
      e = eq;
      While[True,
        hit = Select[Keys[e], KeyExistsQ[rules, #] &];
        If[hit === {}, Break[]];
        Do[If[KeyExistsQ[e, m],
             c = e[m]; e = KeyDrop[e, {m}];
             Do[With[{v = Mod[getD[e, mm, 0] + c rules[m][mm], pp]}, If[v == 0, e = KeyDrop[e, {mm}], e[mm] = v]], {mm, Keys[rules[m]]}]],
           {m, hit}]];
      If[Length[e] > 0,
        piv = Last[SortBy[Keys[e], keyF]];
        If[MemberQ[masters, piv], Print["MCCheck (Laporta): an IBP relation among the masters -- they are not independent"]; Throw[$Failed, mcTag]];
        inv = PowerMod[e[piv], -1, pp];
        rules[piv] = Association[Table[mm -> Mod[-e[mm] inv, pp], {mm, DeleteCases[Keys[e], piv]}]]],
      {eq, eqs}];
    full[m_] := full[m] = If[!KeyExistsQ[rules, m], <|m -> 1|>,
       Module[{out = <||>, sub},
         Do[sub = full[mm];
            Do[With[{v = Mod[getD[out, kk, 0] + rules[m][mm] sub[kk], pp]}, If[v == 0, out = KeyDrop[out, {kk}], out[kk] = v]], {kk, Keys[sub]}],
            {mm, Keys[rules[m]]}];
         out]];
    Do[With[{r = full[t]}, res[t] = If[AllTrue[Keys[r], MemberQ[masters, #] &], Table[getD[r, mi, 0], {mi, masters}], None]], {t, targets}];
    res]];
  While[!done && Length[primes] < 24,
    p = NextPrime[p];
    If[AnyTrue[eqsQ, AnyTrue[Values[#], Mod[Denominator[#], p] == 0 &] &], Continue[]];
    AppendTo[primes, p]; AppendTo[resid, modRows[p]];
    If[Length[primes] >= 2,
      ok = True; exact = <||>;
      Do[
        If[resid[[1]][t] === None, exact[t] = None; Continue[]];
        row = Table[ratrec[ChineseRemainder[Table[resid[[kk]][t][[i]], {kk, Length[primes] - 1}], Most[primes]], Times @@ Most[primes]], {i, nM}];
        If[MemberQ[row, $Failed], ok = False; Break[]];
        If[!AllTrue[Range[nM], Mod[modQ[row[[#]], Last[primes]] - resid[[-1]][t][[#]], Last[primes]] == 0 &], ok = False; Break[]];
        exact[t] = row,
        {t, targets}];
      If[ok, done = True]]];
  If[!done, Print["MCCheck (Laporta): rational reconstruction did not stabilise"]; Return[Infinity]];
  red = Select[targets, exact[#] =!= None &];
  closed = If[red === {}, {}, mcRows[st, red, dval]];
  Do[worst = Max[worst, Max[Join[{0}, Abs[closed[[i]] - exact[red[[i]]]]/Max[1, Max[Join[{0}, Abs[exact[red[[i]]]]]]]]]], {i, Length[red]}];
  If[verbose, Print["MCCheck (Laporta): exact maximal-cut reduction, box dots <= ", dd + 1, ", rank <= ", rr + 1, ", ", Length[eqsQ], " equations, ",
     Length[primes], " primes; ", Length[red], " of ", Length[targets], " targets reduced to the masters ", masters,
     ";  max |closed form - exact| / max(1, |exact|) = ", N[worst, 3], "   (", Round[AbsoluteTime[] - t0, 0.1], " s)"]];
  If[Length[red] < Length[targets], Print["MCCheck (Laporta): not reduced to the masters (enlarge \"LaportaBox\", or the master count is wrong): ",
     Select[targets, exact[#] === None &]]];
  worst];

Options[MCCheck] = {"Grid" -> {1, 2}, "Dimension" -> 37/11, Method -> "IBP", "LaportaBox" -> Automatic, "Verbose" -> True};
MCCheck[MCResult[st_Association], OptionsPattern[]] := Catch[Module[{dmax, rmax, dval, method, verbose, rels, worst = 0},
  {dmax, rmax} = OptionValue["Grid"]; dval = OptionValue["Dimension"]; method = OptionValue[Method]; verbose = TrueQ[OptionValue["Verbose"]];
  If[!ExactNumberQ[dval], dval = Rationalize[dval, 0]];
  rels = ibpRelations[st, dval];
  If[MemberQ[{"IBP", "Both"}, method], worst = Max[worst, ibpResidualCheck[st, rels, dmax, rmax, dval, verbose]]];
  If[MemberQ[{"Laporta", "Both"}, method], worst = Max[worst, laportaCheck[st, rels, dmax, rmax, dval, OptionValue["LaportaBox"], verbose]]];
  worst], mcTag];


(* ::Section:: *)
(*MCCompare: the closed forms against an independent IBP reduction (LiteRed2 in a separate kernel, cached)*)


If[!ValueQ[$MCLiteRedFile], $MCLiteRedFile = Automatic];
If[!ValueQ[$MCReferenceDirectory], $MCReferenceDirectory = Automatic];
If[!ValueQ[$mcRefKernel], $mcRefKernel = None];
If[!ValueQ[$mcCompareResults], $mcCompareResults = {}];

mcRefDir[] := If[StringQ[$MCReferenceDirectory], $MCReferenceDirectory, FileNameJoin[{$mcPackageDir, "MaxCutClosedForm_LiteRed"}]];

(* number as a one-line string, e.g. 1.23e-31 *)
mcSci[x_] := Which[!NumericQ[x], ToString[x], x == 0, "0",
   True, Module[{ax = N[Abs[x]], e}, e = Floor[Log10[ax]];
     If[x < 0, "-", ""] <> ToString[NumberForm[N[ax/10^e], {3, 2}]] <> "e" <> ToString[e]]];

(* the reference kernel: launched once, LiteRed2 and MaxCutClosedForm_LiteRed.wl loaded there *)
mcRefKernel[logQ_] := Module[{k, stat, drv = FileNameJoin[{$mcPackageDir, "MaxCutClosedForm_LiteRed.wl"}]},
  If[MatchQ[$mcRefKernel, _KernelObject] && MemberQ[Kernels[], $mcRefKernel], Return[$mcRefKernel]];
  If[!FileExistsQ[drv], Return["MaxCutClosedForm_LiteRed.wl not found next to MaxCutClosedForm.m (" <> drv <> ")"]];
  k = Quiet[LaunchKernels[1]];
  If[!MatchQ[k, {_KernelObject}], Return["could not launch a parallel kernel (LaunchKernels[1] failed; see Evaluation > Parallel Kernel Configuration)"]];
  k = First[k];
  stat = With[{dr = drv, lf = $MCLiteRedFile, lg = TrueQ[logQ]},
    ParallelEvaluate[MCLiteRed`$File = lf; MCLiteRed`$Log = lg; Get[dr]; {MCLiteRed`$Status, MCLiteRed`$Version}, {k}]];
  stat = If[MatchQ[stat, {{_String, _String}}], First[stat], {"no answer: " <> ToString[Short[stat, 2], InputForm], "none"}];
  If[stat[[1]] =!= "loaded", Quiet[CloseKernels[{k}]]; Return["reference kernel: " <> stat[[1]]]];
  Print["   reference kernel launched: ", stat[[2]]];
  $mcRefKernel = k; k];

MCStopReference[] := (If[MatchQ[$mcRefKernel, _KernelObject], Quiet[CloseKernels[{$mcRefKernel}]]]; $mcRefKernel = None;);

(* the job: plain numbers only (momenta as coefficient vectors in the basis {loops, legs}, masses, external Gram matrix) *)
mcRefJob[st_Association, targets_List, dvals_List, tlim_] := Module[{topo = st["Topology"], nL, pt = st["Point"], masses, gram},
  nL = topo["L"];
  masses = topo["Masses"] /. pt;
  gram = topo["GramSP"][[nL + 1 ;;, nL + 1 ;;]] /. pt;
  If[!VectorQ[masses, NumericQ] || !MatrixQ[gram, NumericQ], Message[MCCompare::kin]; Throw[$Failed, mcTag]];
  <|"Name" -> ToString[st["Name"]], "L" -> nL, "E" -> topo["E"], "NProps" -> st["NProps"], "Coefficients" -> topo["Coefficients"],
    "Masses" -> masses, "ExternalGram" -> gram, "Targets" -> targets, "DValues" -> dvals, "TimeLimit" -> tlim|>];

(* cached LiteRed result, or a new run in the reference kernel *)
mcReference[job_Association, logQ_, verbose_] := Module[{dir = mcRefDir[], key, file, lr, k, t0},
  key = IntegerString[Hash[KeyDrop[job, {"Name", "TimeLimit"}]], 36];
  file = FileNameJoin[{dir, "cache", StringReplace[job["Name"], Except[WordCharacter] -> "_"] <> "_" <> key <> ".m"}];
  If[FileExistsQ[file],
    lr = Quiet[Get[file]];
    If[AssociationQ[lr] && lr["Status"] === "ok", Return[Append[lr, "Cache" -> FileNameTake[file]]]]];
  If[verbose, Print["   LiteRed: new reduction (", Length[job["Targets"]], " integrals; the result is cached in ", FileNameJoin[{dir, "cache"}], ")"]];
  k = mcRefKernel[logQ];
  If[!MatchQ[k, _KernelObject], Return[<|"Status" -> "failed", "Message" -> k|>]];
  t0 = AbsoluteTime[];
  lr = With[{jb = Append[job, "WorkDir" -> FileNameJoin[{dir, "work"}]]}, ParallelEvaluate[MCLiteRed`MCLiteRedReduce[jb], {k}]];
  lr = If[MatchQ[lr, {_Association}], First[lr],
    <|"Status" -> "failed", "Message" -> "no answer from the reference kernel: " <> ToString[Short[lr, 2], InputForm]|>];
  If[lr["Status"] === "ok",
    If[!DirectoryQ[DirectoryName[file]], CreateDirectory[DirectoryName[file]]];
    Put[lr, file]; lr = Append[lr, "Cache" -> FileNameTake[file] <> " (new, " <> ToString[Round[AbsoluteTime[] - t0, 0.1]] <> " s)"]];
  lr];

(* random targets in the box dots <= dmax, rank <= rmax (reproducible: seeded by the topology name) *)
mcSample[st_Association, n_Integer, {dmax_Integer, rmax_Integer}, seed_] := Module[{nP = st["NProps"], cands},
  cands = Select[gridSeeds[nP, st["N"] - nP, dmax, rmax], !MemberQ[st["Masters"], #] &];
  SortBy[BlockRandom[SeedRandom[seed]; RandomSample[cands, Min[n, Length[cands]]]],
    {Total[Take[#, nP] - 1], -Total[Drop[#, nP]], #} &]];

Options[MCCompare] = {"Dimensions" -> {37/11, 43/13}, "Sample" -> {8, {2, 4}}, "Seed" -> Automatic, "Show" -> 4, "TimeLimit" -> 1800,
   "Masters" -> Automatic, "LiteRedLog" -> False, "Verbose" -> True};
MCCompare::tg = "Targets must be index vectors of length `1` (propagator powers >= 1, ISP powers <= 0).";
MCCompare::kin = "The masses or the external scalar products are not numbers at the kinematic point.";
MCCompare::ref = "No reference reduction: `1`";

MCCompare[r : MCResult[st_Association], tg_List : {}, OptionsPattern[]] := Catch[Module[
  {verbose, dvals, targets, nP = st["NProps"], nAll = st["N"], sample, seed, job, lr, lrM, nLR, ourM, nM, rr, mode, mastersRR, tested,
   rowsLR, rowOf, cm, bm, pred, devs = {}, tab = {}, show, dev, ok, tol, summary},
  verbose = TrueQ[OptionValue["Verbose"]]; show = OptionValue["Show"];
  dvals = Rationalize[#, 0] & /@ Flatten[{OptionValue["Dimensions"]}];
  If[!VectorQ[dvals, ExactNumberQ], Message[MCReduce::dim, dvals]; Throw[$Failed, mcTag]];
  targets = Which[tg === {}, {}, VectorQ[tg, IntegerQ], {tg}, MatrixQ[tg, IntegerQ], tg, True, Message[MCCompare::tg, nAll]; Throw[$Failed, mcTag]];
  Do[If[Length[n] =!= nAll || AnyTrue[Take[n, nP], # < 1 &] || AnyTrue[Drop[n, nP], # > 0 &], Message[MCCompare::tg, nAll]; Throw[$Failed, mcTag]], {n, targets}];
  sample = OptionValue["Sample"];
  If[MatchQ[sample, {_Integer?Positive, {_Integer, _Integer}}],
    seed = OptionValue["Seed"]; If[seed === Automatic, seed = Hash[ToString[st["Name"]]]];
    targets = Join[targets, mcSample[st, sample[[1]], sample[[2]], seed]]];
  targets = DeleteDuplicates[Join[targets, st["Masters"]]];
  If[verbose, Print["MCCompare: ", st["Name"], " at ", st["Point"], " against LiteRed2, ", Length[targets], " integrals, d = ", dvals]];
  job = mcRefJob[st, targets, dvals, OptionValue["TimeLimit"]];
  lr = mcReference[job, OptionValue["LiteRedLog"], verbose];
  If[lr["Status"] =!= "ok", Message[MCCompare::ref, lr["Message"]]; Throw[$Failed, mcTag]];
  lrM = lr["Masters"]; nLR = Length[lrM]; ourM = st["Masters"]; nM = Length[ourM];
  If[verbose, Print["   ", lr["Version"], ": basis, IBPs, symmetries and the top sector solved in ", lr["TimeBasis"], " s, ",
     lr["NTargets"], " integrals reduced in ", lr["TimeReduce"], " s   [", lr["Cache"], "]"]];
  (* which masters to compare in *)
  rr = r;
  mode = Which[
    nM == 0, "no master",
    Sort[lrM] === Sort[ourM], "same masters",
    nLR == nM && OptionValue["Masters"] =!= "Ours",
      rr = Quiet[MCClosedForm[st["Topology"], st["Point"], "Masters" -> lrM, "Verbose" -> False, WorkingPrecision -> st["WorkingPrecision"]]];
      If[Head[rr] === MCResult, "closed forms recomputed with LiteRed's masters", rr = r; "projection"],
    True, "projection"];
  mastersRR = MCMasters[rr];
  tested = Select[targets, !MemberQ[mastersRR, #] &];
  Do[
    rowsLR = AssociationThread[targets, lr["Rows"][dv]];
    If[nM == 0,
      (* no master on the maximal cut: LiteRed must reduce everything to zero (sub-sectors) *)
      AppendTo[devs, Max[Join[{0}, Abs[Flatten[Values[rowsLR]]]]]];
      If[dv === First[dvals], tab = Table[{n, rowsLR[n], "{}  (no master)", If[rowsLR[n] === {} || AllTrue[rowsLR[n], # == 0 &], "0", "nonzero"]}, {n, Take[tested, UpTo[show]]}]],
      rowOf[e_] := If[MemberQ[lrM, e], UnitVector[nLR, First[FirstPosition[lrM, e]]], rowsLR[e]];
      If[tested =!= {},
        cm = MCReduce[rr, tested, dv];                           (* closed-form coefficients on the masters of rr *)
        bm = rowOf /@ mastersRR;                                 (* LiteRed reduction of those masters *)
        pred = cm . bm;                                          (* closed form expressed in LiteRed's masters *)
        Do[AppendTo[devs, Max[Abs[pred[[i]] - rowsLR[tested[[i]]]]]/Max[1, Max[Abs[rowsLR[tested[[i]]]]]]], {i, Length[tested]}];
        If[dv === First[dvals],
          tab = Table[{tested[[i]], N[rowsLR[tested[[i]]], 12], N[pred[[i]], 12],
              mcSci[Max[Abs[pred[[i]] - rowsLR[tested[[i]]]]]/Max[1, Max[Abs[rowsLR[tested[[i]]]]]]]}, {i, Min[show, Length[tested]]}]]]],
    {dv, dvals}];
  dev = Max[Join[{0}, devs]];
  tol = If[st["Kind"] === 0 || nM == 0, 0, 10^(-st["WorkingPrecision"]/2)];
  ok = TrueQ[dev <= tol] && (nM > 0 || nLR == 0);
  If[verbose,
    Print["   masters: closed form ", mastersRR, ";  LiteRed ", lrM, "   (", mode, ")"];
    If[mode === "projection", Print["   (different masters: checked I(n) = c(n).L(masters), L = LiteRed's reduction of the closed-form masters)"]];
    If[tab =!= {}, Print[Grid[Prepend[tab, {"n", "LiteRed, exact (12 digits shown), d = " <> ToString[First[dvals], InputForm],
        "closed form" <> If[mode === "projection", ", in LiteRed's masters", ""], "rel. difference"}], Frame -> All, Alignment -> Left]]];
    If[StringQ[lr["Example"]] && lr["Example"] =!= "", Print["   LiteRed, exact in d, for ", First[targets], ":  ", StringTake[lr["Example"], UpTo[300]],
        If[StringLength[lr["Example"]] > 300, " ...", ""]]];
    Print["   ", Length[tested], " integrals x ", Length[dvals], " values of d:  max relative difference closed form vs LiteRed = ", mcSci[dev],
      If[ok, "   -> AGREE", "   -> DISAGREE"]]];
  summary = <|"Topology" -> st["Name"], "Masters" -> nM, "LiteRedMasters" -> nLR, "Integrals" -> Length[tested], "d" -> dvals,
     "MaxDifference" -> dev, "Agree" -> ok, "LiteRedTime" -> lr["TimeBasis"] + lr["TimeReduce"], "Mode" -> mode|>;
  $mcCompareResults = Append[DeleteCases[$mcCompareResults, s_ /; s["Topology"] === st["Name"] && s["Masters"] === nM && s["Mode"] === mode], summary];
  dev], mcTag];

MCCompareSummary[] := (Print[Grid[Prepend[
     ({#["Topology"], #["Masters"], #["LiteRedMasters"], #["Integrals"], ToString[#["d"], InputForm], mcSci[#["MaxDifference"]],
       If[#["Agree"], "AGREE", "DISAGREE"], #["LiteRedTime"]} &) /@ $mcCompareResults,
     {"topology", "masters", "LiteRed masters", "integrals", "d", "max rel. difference", "", "LiteRed time [s]"}], Frame -> All, Alignment -> Left]];
   $mcCompareResults);


(* ::Section:: *)
(*MCSelfTest: special functions and the interval formula against direct quadrature*)


MCSelfTest[] := Module[{pr = 30, wp2 = 45, res = {}, rep, quad, a, c, bs, xs, pts, cyc, val, f, direct, ok, t0},
  (* every test prints its line as soon as it is done *)
  rep[name_, x_, y_, tol_] := Module[{e = Block[{$MinPrecision = 0}, N[Abs[x - y]/Max[Abs[y], 10^-300], 6]]},
    AppendTo[res, {name, ScientificForm[e, 3], ScientificForm[N[tol], 1], e < tol}];
    Print["   MCSelfTest: ", name, ":  relative difference ", ScientificForm[e, 3], If[e < tol, "  ok", "  FAILED"], "   (", Round[AbsoluteTime[] - t0, 0.01], " s)"];
    t0 = AbsoluteTime[]];
  (* direct quadrature (double-exponential rule handles the integrable endpoint singularities) *)
  quad[g_, range_List] := Block[{$MinPrecision = 0},
    NIntegrate[g, range, WorkingPrecision -> pr + 10, PrecisionGoal -> 20, MaxRecursion -> 20, Method -> "DoubleExponential"]];
  t0 = AbsoluteTime[];
  (* 1. two-variable series vs built-in AppellF1 *)
  fixedPrec[wp2,
    a = N[3/10, wp2]; c = N[19/10, wp2]; bs = N[{7/10, -2/5}, wp2]; xs = N[{3/10, -1/2}, wp2];
    rep["F_D^(2) series vs AppellF1", MCLauricellaFD[a, bs, c, xs], AppellF1[a, bs[[1]], bs[[2]], c, xs[[1]], xs[[2]]], 10^-20]];
  (* 2. three variables: series vs continued Euler integral, continued parameters *)
  fixedPrec[wp2,
    a = N[-13/10, wp2]; c = N[2/5, wp2]; bs = N[{7/10, -3/10, 13/10}, wp2]; xs = N[{3/10, -1/2, 1/4 + I/5}, wp2];
    rep["F_D^(3) series vs continued Euler integral", fdSeries[a, bs, c, xs, wp2], eulerContinued[a, bs, c, xs, wp2] Gamma[c]/(Gamma[a] Gamma[c - a]), 10^-20]];
  (* 3. large and complex arguments: continued Euler integral vs direct quadrature (convergent parameters) *)
  fixedPrec[wp2,
    a = N[13/10, wp2]; c = N[21/10, wp2]; bs = N[{7/10, -3/10, 13/10}, wp2]; xs = N[{5/2 + 6/5 I, -3, 1/2}, wp2];
    direct = quad[u^(13/10 - 1) (1 - u)^(21/10 - 13/10 - 1) (1 - (5/2 + 6/5 I) u)^(-7/10) (1 + 3 u)^(3/10) (1 - u/2)^(-13/10), {u, 0, 1}];
    rep["F_D^(3), |x| > 1, Euler vs quadrature", Gamma[a] Gamma[c - a]/Gamma[c] MCLauricellaFD[a, bs, c, xs], direct, 10^-10]];
  (* 4. polynomial factor with a large argument *)
  fixedPrec[wp2,
    a = N[13/10, wp2]; c = N[21/10, wp2]; bs = N[{7/10, -3, 13/10}, wp2]; xs = N[{-1/2, 7, 1/3}, wp2];
    direct = quad[u^(13/10 - 1) (1 - u)^(21/10 - 13/10 - 1) (1 + u/2)^(-7/10) (1 - 7 u)^3 (1 - u/3)^(-13/10), {u, 0, 1}];
    rep["F_D^(3), polynomial factor, x = 7", Gamma[a] Gamma[c - a]/Gamma[c] MCLauricellaFD[a, bs, c, xs], direct, 10^-10]];
  (* 5. interval formula: real chamber with |.| and signed powers, reference infinity and two reference points *)
  pts = {<|"Eabs" -> -13/10, "N" -> 0|>, <|"Eabs" -> 7/10, "N" -> 2|>, <|"Eabs" -> 0, "N" -> 3|>, <|"Eabs" -> -6/10, "N" -> 0|>, <|"Eabs" -> 9/20, "N" -> 1|>};
  val = {-2, 0, 1/2, 3/2, 4}[[#]] &;
  f[tt_] := Times @@ Table[Abs[tt - val[k]]^pts[[k]]["Eabs"] (tt - val[k])^pts[[k]]["N"], {k, 5}];
  direct = quad[f[tt], {tt, 0, 1/2, 3/2}];
  Do[cyc = <|"Kind" -> "real", "p" -> 2, "q" -> 4, "Ref" -> rf|>;
     fixedPrec[wp2, rep["real chamber, reference " <> If[rf["Type"] === "inf", "infinity", "t_" <> ToString[rf["Label"]]], N[intervalExpr[cyc, pts, val], wp2], direct, 10^-10]],
     {rf, {<|"Type" -> "inf"|>, <|"Type" -> "point", "Label" -> 5|>, <|"Type" -> "point", "Label" -> 1|>}}];
  (* 6. unbounded chamber *)
  pts = {<|"Eabs" -> -13/10, "N" -> 0|>, <|"Eabs" -> 7/10, "N" -> 0|>, <|"Eabs" -> -6/10, "N" -> -1|>, <|"Eabs" -> -9/20, "N" -> 1|>};
  val = {-2, 0, 3/2, 4}[[#]] &;
  f[tt_] := Times @@ Table[Abs[tt - val[k]]^pts[[k]]["Eabs"] (tt - val[k])^pts[[k]]["N"], {k, 4}];
  direct = quad[f[tt], {tt, 4, 5, Infinity}];
  cyc = <|"Kind" -> "real", "p" -> 4, "q" -> Infinity, "Ref" -> <|"Type" -> "point", "Label" -> 3|>|>;
  fixedPrec[wp2, rep["unbounded chamber", N[intervalExpr[cyc, pts, val], wp2], direct, 10^-10]];
  (* 7. complex segment (analytic convention) vs quadrature along the segment *)
  pts = {<|"Eabs" -> -3/10, "N" -> 0|>, <|"Eabs" -> 7/10, "N" -> 1|>, <|"Eabs" -> -6/10, "N" -> 0|>};
  val = {1 + I, 3 - I/2, -1}[[#]] &;
  direct = quad[(val[2] - val[1]) (val[2] - val[1])^(-3/10) u^(-3/10) (val[1] - val[2])^(17/10) (1 - u)^(17/10)
        (val[1] - val[3])^(-6/10) (1 - (val[2] - val[1])/(val[3] - val[1]) u)^(-6/10), {u, 0, 1}];
  cyc = <|"Kind" -> "complex", "p" -> 1, "q" -> 2, "Ref" -> <|"Type" -> "none"|>|>;
  fixedPrec[wp2, rep["complex segment", N[intervalExpr[cyc, pts, val], wp2], direct, 10^-10]];
  ok = And @@ res[[All, 4]];
  Print[Grid[Prepend[res, {"test", "relative difference", "tolerance", "passed"}], Frame -> All]];
  Print["MCSelfTest: ", If[ok, "all passed", "FAILED"]];
  ok];

End[];
EndPackage[];
