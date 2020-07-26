(* ::Package:: *)

Block[{Print = Identity},
      <<FeynArts`;
      <<IGraphM`;
];

(* ----------------------------------------------- *)
(* code for working with topologies in Mathematica *)


(* 0 => no messages, 1 => summary messages, 2 => all messages *)
$FAVerbose = 1;


makeTopologies::usage = "Creates topology objects.";
makeTopologies[nScalars_, nFermions_] :=
  CreateTopologies[
    0, (* only want tree-level topologies *)
    (nScalars + nFermions) -> 0,
    ExcludeTopologies -> {SelfEnergies, Tadpoles, Triangles},
    Adjacencies -> {3, 4}
  ];


makeDiagrams::usage = "Creates diagram objects with particle and momentum information.";
makeDiagrams[nScalars_, nFermions_] :=
  InsertFields[
    makeTopologies[nScalars, nFermions],
    Join[Table[F[2], {nFermions}], Table[S[1], {nScalars}]] -> {},
    InsertionLevel -> Generic,
    (* Don't want any vectors in the tree-level topology *)
    ExcludeParticles -> {V}
  ];


(* Define diagrams *)
(* Quiet[ *)
(*   Block[{Print = Identity}, *)
(*         $D5Diagrams = makeDiagrams[2, 2]; *)
(*         $D7FourFermionDiagrams = makeDiagrams[1, 4]; *)
(*         $D9SixFermionDiagrams = makeDiagrams[0, 6]; *)
(*         $D11SixFermionDiagrams = makeDiagrams[2, 6]; *)
(*   ] *)
(* ]; *)


(* Diagrams defined above invovle many that only differ by momentum arrows.
   Strategy is to use IGraphM to filter out isomorphic graphs and keep unique
   topologies. *)

(* Some auxiliary functions to convert diagrams to graph objects *)

stylise::usage = "IGraphM isomorphic graph comparison functions need edges to be labeled by integers.";
stylise[S] := 1; (* Scalar *)
stylise[F] := 2; (* Fermion *)

(* FeynArts Vertex objects keep numerical label and adjacency. Need to map these
   to a unique integer for graph comparison. *)
intMapper[Vertex[a_, b_]] := 2^a 3^b;

(* Style keeps track of fermion or scalar *)
extractStyleMapper[Rule[left_, right_]] := left -> If[Length[right] == 0, right, Head[right]];

extractStyle::usage = "IGraphM isomorphic graph comparison functions need edge list in this form.";
extractStyle[rules__] := Map[extractStyleMapper, List[rules]];


extractEdgeReplacement::usage = "Read FeynArts vertex.";
extractEdgeReplacement[Propagator[x_][Vertex[y_][z_], Vertex[w_][q_], Field[f_]]] :=
  {UndirectedEdge[intMapper[Vertex[y, z]], intMapper[Vertex[w, q]]] -> stylise[Field[f]]};

getStyle[Rule[Topology[_][topology__], Insertions[_][FeynmanGraph[__][repl__]]]] :=
  extractStyle[repl];

getEdgeReplacements[Rule[Topology[_][topology__], Insertions[_][FeynmanGraph[__][repl__]]]] :=
  Map[extractEdgeReplacement, List[topology]];


completionTopologies::usage = "Returns list of lists of graph object and edge colours (for isomorphic graph filtering)";
completionTopologies[feynmanDiagrams_] :=
  Table[{Graph[Keys[Apply[Join, getEdgeReplacements[d]]]],
         "EdgeColors" -> Association[Apply[Join, getEdgeReplacements[d]] //. getStyle[d]]},
        {d, Apply[List, feynmanDiagrams]}
  ];

removeIsomorphic::usage = "Filters isomorphic graphs from input.";
removeIsomorphic[diags_] :=
  Block[{i, j, copyDiags},
        Print["Mathematica: Filtering isomorphic graphs..."];
        copyDiags = diags;
        i = 1;
        While[i < Length[copyDiags],
              j = i + 1;
              While[j <= Length[copyDiags],
                    If[IGVF2IsomorphicQ[copyDiags[[i]], copyDiags[[j]]],
                       copyDiags = Delete[copyDiags, j];, j++; Nothing
                    ];
              ];
              i++];
        copyDiags
  ];


particleStyle[1] := Sequence[Dashed, Gray];
particleStyle[2] := Gray;

decorateGraph::usage = "Displays graph nicely with scalar and fermion lines.";
decorateGraph[gr_, vrts_, edgs_] :=
  Block[{gv},
        gv = Fold[SetProperty[{#1, #2}, {VertexSize -> Large}]&, gr, Keys[vrts]];
        Fold[SetProperty[{#1, #2}, EdgeStyle -> {particleStyle[edgs[[Key[#2]]]], Thick}] &, gv, Keys[edgs]]
  ];


(* constructPartition *)

(* getFirst[y_List] := Map[First, y]; *)
(* getRest[x_List] := getRest[Rest[x] /. First[x]]; *)
(* getRest[x_List /; Length[x] == 1] := x; *)
(* grouper[x_] := Last[x]; *)
(* constructPartition[{g_Graph, Rule["EdgeColors", assoc_]}] := *)
(*   Block[{partition, styleReplacement}, *)
(*         partition = Map[Last, getRest[SortBy[Normal[Map[getFirst, GroupBy[EdgeList[g], grouper]]], Length[#[[2]]]&]]]; *)
(*         styleReplacement = Map[#[[1, 1]] -> If[#[[2]] == 2, F[#[[1, 1]]], S[#[[1, 1]]]] &, Normal[assoc]]; *)
(*         partition /. styleReplacement *)
(*   ]; *)

externalVertices[g_] :=
  Keys[Select[MapThread[Rule, {VertexList[g], AdjacencyList[g]}],
    Length[#[[2]]] == 1 &]];

SetAttributes[immediateConnections, HoldRest];
immediateConnections[g_, v_, seen_] :=
  Block[{ext = externalVertices[g]},
        AppendTo[seen, v];
        If[MemberQ[ext, v],
           v,
           Table[
             Which[
               i[[1]] == v, If[MemberQ[seen, i[[2]]], Nothing, i[[2]]],
               i[[2]] == v, If[MemberQ[seen, i[[1]]], Nothing, i[[1]]],
               True, Nothing
             ], {i, EdgeList[g]}]
        ]
  ]

constructPartitionList[g_] :=
  Block[{v = Max[VertexList[g]], seen = {}, new = {}, out, prev,
         counter = 1},
        out = immediateConnections[g, v, seen];
        prev = out;
        While[new =!= prev,
              prev = out;
              new = Map[immediateConnections[g, #, seen] &, out, {counter}];
              out = new;
              counter++;
        ];
        out
  ];

constructPartition[{g_Graph, Rule["EdgeColors", assoc_]}] :=
  Block[{partition, styleReplacement},
        partition = constructPartitionList[g];
        styleReplacement = Map[#[[1, 1]] -> If[#[[2]] == 2, F[#[[1, 1]]], S[#[[1, 1]]]] &, Normal[assoc]];
        partition /. styleReplacement
  ];


outFileName[path_, type_, prefix_, n_] :=
  Block[{outpath},
        outpath = FileNameJoin[{path, type, prefix <> "_" <> ToString[n]}];
        outpath
  ];

exportGraphsAndPartitions[diags_, prefix_, path_] :=
  Block[
    {cleanGraphs, decoratedGraphs, edgeLists, partitions, readableEdgeLists,
     diagramPath, graphPath, partitionPath},

    cleanGraphs = removeIsomorphic[completionTopologies[diags]];
    decoratedGraphs = Map[decorateGraph[#[[1]], <||>, #[[2, 2]]]&, cleanGraphs];
    edgeLists = Map[EdgeList[First[#]]&, cleanGraphs];
    partitions = Map[constructPartition, cleanGraphs];

    readableEdgeLists = Map[{#[[1]], #[[2]]} &, edgeLists, {2}];
    Do[
      diagramPath = outFileName[path, "diagrams", prefix, i] <> ".png";
      Export[diagramPath, decoratedGraphs[[i]], ImageResolution -> 300];
      Print["Mathematica: "<>diagramPath<>" written!"];

      graphPath = outFileName[path, "graphs", prefix, i] <> ".csv";
      Export[graphPath, readableEdgeLists[[i]]];
      Print["Mathematica: "<>graphPath<>" written!"];

      partitionPath = outFileName[path, "partitions", prefix, i] <> ".dat";
      Export[partitionPath, FortranForm[partitions[[i]]]];
      Print["Mathematica: "<>partitionPath<>" written!"];

    , {i, Length[partitions]}];
  ];


(* ------------------------------ *)
(* process command-line arguments *)

output = $ScriptCommandLine[[2]];
nScalarFields = ToExpression[$ScriptCommandLine[[3]]];
nFermionFields = ToExpression[$ScriptCommandLine[[4]]];

Print["Mathematica: FeynArts generating topologies..."];
diagrams = makeDiagrams[nScalarFields, nFermionFields];
prefix = ToString[nScalarFields] <> "s" <> ToString[nFermionFields] <> "f";
exportGraphsAndPartitions[diagrams, prefix, output];
