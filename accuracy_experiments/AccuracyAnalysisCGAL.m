(* ::Package:: *)

(* Check for necessary command line arguments *)
If[Length[$ScriptCommandLine] <= 1,
    Print["Usage: wolframscript -file this_package.m <latitudes>"];
    Return[$Failed];
];

(* Extract latitude argument from the command line *)
latitudeString = $ScriptCommandLine[[2]];
latitudePoints = ToExpression[StringSplit[latitudeString, ","]];

(* Set the base path as the parent directory of the current working directory *)
ClearAll[absolutePath];
absolutePath = DirectoryName[Directory[]]

(* Convert latitude and longitude to Cartesian coordinates on a unit sphere *)
ClearAll[latLonToXYZ];
latLonToXYZ[lat_, lon_] :=
  With[{θ = (90 - lat) Degree, ϕ = lon Degree},
    {Sin[θ] Cos[ϕ], Sin[θ] Sin[ϕ], Cos[θ]}
  ];

(* Check if the span is within a given tolerance *)
ClearAll[isSpanCorrect];
isSpanCorrect[startPoint_, endPoint_, spanRadians_, toleranceRadians_] :=
  Module[{arcLength},
    arcLength = ArcCos[Normalize[startPoint] . Normalize[endPoint]];
    arcLength <= spanRadians + toleranceRadians
  ];

(* Generate random great circle arcs with a specified span *)
ClearAll[generateGreatCircleArcs];
generateGreatCircleArcs[numArcs_Integer, latRange_List, lonRange_List, precision_: MachinePrecision, spanDegrees_: -1.0, toleranceDegrees_: 0.1] :=
  Module[{arcs = {}, startPoint, endPoint, spanRadians, toleranceRadians, startLat, startLon, endLat, endLon},
    spanRadians = SetPrecision[If[spanDegrees > 0, spanDegrees Degree, spanDegrees], precision];
    toleranceRadians = SetPrecision[toleranceDegrees Degree, precision];

    While[Length[arcs] < numArcs,
      startLat = RandomReal[latRange];
      startLon = RandomReal[lonRange];
      endLat = RandomReal[latRange];
      endLon = RandomReal[lonRange];

      startPoint = latLonToXYZ[SetPrecision[startLat, precision], SetPrecision[startLon, precision]];
      endPoint = latLonToXYZ[SetPrecision[endLat, precision], SetPrecision[endLon, precision]];

      If[spanDegrees <= 0 || isSpanCorrect[startPoint, endPoint, spanRadians, toleranceRadians],
        AppendTo[arcs, {startPoint, endPoint}]
      ]
    ];
    arcs
  ];

(* Compute intersection coordinates with a specified precision *)
ClearAll[gcaConstLatIntersectionCoordinatesNewEqn];
gcaConstLatIntersectionCoordinatesNewEqn[pointA_List, pointB_List, constZ_, precision_: MachinePrecision] :=
  Module[{n, nx, ny, nz, nxSquared, nySquared, normNSquared, sTilde, px, py},
    n = Cross[pointA, pointB];
    {nx, ny, nz} = n;
    nxSquared = SetPrecision[nx^2, precision];
    nySquared = SetPrecision[ny^2, precision];
    normNSquared = nxSquared + nySquared + SetPrecision[nz^2, precision];

    sTilde = Sqrt[nxSquared + nySquared - normNSquared SetPrecision[constZ^2, precision]];
    px = -(1/(nxSquared + nySquared)) (SetPrecision[constZ, precision] nx nz + sTilde ny);
    py = -(1/(nxSquared + nySquared)) (SetPrecision[constZ, precision] ny nz - sTilde nx);

    {SetPrecision[px, precision], SetPrecision[py, precision], SetPrecision[constZ, precision]}
  ];

(* Compute relative error with adjustable precision *)
ClearAll[relativeError];
relativeError[baseline_, comparison_, constZ_, precision_] := 
  Module[{diff, nom, relError},
    diff = (SetPrecision[baseline, precision] - SetPrecision[comparison, precision])^2;
    nom = baseline[[1]]^2 + baseline[[2]]^2 + SetPrecision[constZ, precision]^2;
    relError = Sqrt[(diff[[1]] + diff[[2]]) / nom];
    N[relError, precision]
  ];

(* Format latitude as a string *)
ClearAll[formatLatitude];
formatLatitude[lat_] := StringReplace[
  ToString[NumberForm[N[lat], {Infinity, 6}, ExponentFunction -> (Null &)]],
  {"." -> "_", " " -> ""}
];

(* Import and reconstruct arcs from CSV data *)
ClearAll[importAndReconstructArcs];
importAndReconstructArcs[arcFilePath_] := Module[{data, reconstructedData},
    data = Import[arcFilePath, "CSV", "HeaderLines" -> 1];
    reconstructedData = data /. {
        {axs_, axe_, ays_, aye_, azs_, aze_, bxs_, bxe_, bys_, bye_, bzs_, bze_, czs_, cze_} :> {
            axs * 2^axe, ays * 2^aye, azs * 2^aze,
            bxs * 2^bxe, bys * 2^bye, bzs * 2^bze,
            czs * 2^cze
        }
    };
    reconstructedData
];

(* Import and process CGAL benchmark data *)
ClearAll[importAndProcessCGALBenchmarkData, reconstruct];
reconstruct[significand_, exponent_] := significand * 2^exponent;


importAndProcessCGALBenchmarkData[cgalResultsFile_] := Module[{cgalData},
    cgalData = Rest@Import[cgalResultsFile];
    cgalData = Map[reconstruct @@@ Partition[#, 2] &, cgalData];
    <|"CGALResults" -> cgalData|>
];

(* Analyze errors in CGAL results *)
ClearAll[analyzeCGALErrors, computeBaselineResults, computeRelativeErrors];

computeBaselineResults[arcData_, analysisPrecision_] :=
  Table[
    With[{pointA = arc[[1 ;; 3]], pointB = arc[[4 ;; 6]], constZ = arc[[7]]},
      gcaConstLatIntersectionCoordinatesNewEqn[pointA, pointB, constZ, analysisPrecision]
    ],
    {arc, arcData}
  ];

computeRelativeErrors[baselineResults_, benchmarkResults_, analysisPrecision_] :=
  Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkResults[[i, 1 ;; 2]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[baselineResults]}
  ];


summarizeData[data_] := Max[data] 
analyzeCGALErrors[arcData_, benchmarkData_, analysisPrecision_] :=
  Module[{baselineResults, cgalFirstErrors, cgalSecondErrors, cgalErrors, meanCgalErrors},
    baselineResults = computeBaselineResults[arcData, analysisPrecision];
    cgalFirstErrors = computeRelativeErrors[baselineResults, benchmarkData["CGALResults"][[All, 1 ;; 2]], analysisPrecision];
    cgalSecondErrors = computeRelativeErrors[baselineResults, benchmarkData["CGALResults"][[All, 3 ;; 4]], analysisPrecision];
    cgalErrors = MapThread[Min, {cgalFirstErrors, cgalSecondErrors}];
    meanCgalErrors = summarizeData[cgalErrors];
    <|"CGALMethodErrors" -> meanCgalErrors|>
  ];

(* Define output directory and analyze errors *)
analysisPrecision = 100;
outputDirectory = absolutePath <> "mathematica_data/";

cgalErrorsList = {};
For[i = 1, i < Length[latitudePoints], i++,
  startLatStr = formatLatitude[latitudePoints[[i]]];
  endLatStr = formatLatitude[latitudePoints[[i + 1]]];
  latRangeStr = startLatStr <> "_" <> endLatStr;
  
  arcData = importAndReconstructArcs[absolutePath <> "generated_arcs/" <> latRangeStr <> "Arcs_Exponent.csv"];
  benchmarkData = importAndProcessCGALBenchmarkData[absolutePath <> "benchmark_results/" <> latRangeStr <> "_cgal_double.csv"];
  
  errorAnalysisResults = analyzeCGALErrors[arcData, benchmarkData, analysisPrecision];
  latitudeRangeMidpoint = N[(latitudePoints[[i]] + latitudePoints[[i + 1]])/2, analysisPrecision];
  
  AppendTo[cgalErrorsList, {latitudeRangeMidpoint, errorAnalysisResults["CGALMethodErrors"]}];
];

Export[outputDirectory <> "CGAL_Relative_Error.csv", cgalErrorsList, "CSV"];
Print["Done CSV Writing"];
