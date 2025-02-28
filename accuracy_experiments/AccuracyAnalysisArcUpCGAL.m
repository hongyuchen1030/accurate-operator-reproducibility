(* ::Package:: *)

(* Check for necessary command line arguments *)
If[Length[$ScriptCommandLine] <= 2,
    Print["Usage: wolframscript -file this_package.m <mpfr_precisions> <offsets>"];
    Return[$Failed];
];

(* Extract arguments from the command line *)
mpfrPrecisionString = $ScriptCommandLine[[2]];
latitudeString = $ScriptCommandLine[[3]];

(* Convert the comma-separated strings to lists *)
mpfrPrecisionList = ToExpression[StringSplit[mpfrPrecisionString, ","]];
latitudePoints = ToExpression[StringSplit[latitudeString, ","]];


ClearAll[absolutePath];
absolutePath = DirectoryName[Directory[]]

(* Convert latitude and longitude to Cartesian coordinates on a unit sphere *)
ClearAll[latLonToXYZ];
latLonToXYZ[lat_, lon_] := 
  With[{\[Theta] = (90 - lat) Degree, \[CurlyPhi] = lon Degree},
    {Sin[\[Theta]] Cos[\[CurlyPhi]], Sin[\[Theta]] Sin[\[CurlyPhi]], Cos[\[Theta]]}
  ];


(* Check if the span is within the given tolerance *)
ClearAll[isSpanCorrect];
isSpanCorrect[startPoint_, endPoint_, spanRadians_, toleranceRadians_] := 
  Module[{arcLength},
    arcLength = ArcCos[Normalize[startPoint] . Normalize[endPoint]];
    arcLength <= spanRadians + toleranceRadians
  ];


(* Function to generate random great circle arcs with a specified span *)
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

(* Function to calculate intersection coordinates with a specified precision *)
ClearAll[gcaConstLatIntersectionCoordinatesNewEqn];
gcaConstLatIntersectionCoordinatesNewEqn[pointA_List, pointB_List, constZ_, precision_: MachinePrecision] := 
  Module[{n, nx, ny, nz, nxSquared, nySquared, nzSquared, normNSquared, sTilde, px, py},
    n = Cross[pointA, pointB];
    {nx, ny, nz} = n;
    nxSquared = SetPrecision[nx^2, precision];
    nySquared = SetPrecision[ny^2, precision];
    nzSquared = SetPrecision[nz^2, precision];
    normNSquared = nxSquared + nySquared + nzSquared ;
    sTilde = Sqrt[nxSquared + nySquared - normNSquared SetPrecision[constZ^2, precision]];
    px = -(1/(nxSquared + nySquared)) (SetPrecision[constZ, precision] nx nz + sTilde ny);
    py = -(1/(nxSquared + nySquared)) (SetPrecision[constZ, precision] ny nz - sTilde nx);
    {SetPrecision[px, precision], SetPrecision[py, precision], SetPrecision[constZ, precision]}
  ];
  
(* Define a function to compute the relative error with adjustable precision *)
ClearAll[relativeError];
relativeError[baseline_, comparison_, constZ_, precision_] := Module[
  {diff, nom, relError, baselineSetPrec, comparisonSetPrec, constZSetPrec},
  (* Ensure inputs are treated with specified precision *)
  baselineSetPrec = SetPrecision[baseline, precision];
  comparisonSetPrec = SetPrecision[comparison, precision];
  constZSetPrec = SetPrecision[constZ, precision];
  
  diff = (baselineSetPrec - comparisonSetPrec)^2;  (* Calculate squared differences *)
  nom = baselineSetPrec[[1]]^2 + baselineSetPrec[[2]]^2 + constZSetPrec^2;  (* Calculate the norm squared *)
  
  (* Calculate relative error with specified precision *)
  relError = Sqrt[(diff[[1]] + diff[[2]]) / nom];
  N[relError, precision]  (* Ensure final result is returned with specified precision *)
];

(* Define a function to format the latitude similar to your C++ function *)
ClearAll[formatLatitude];
formatLatitude[lat_] := 
  StringReplace[
    ToString[NumberForm[N[lat], {Infinity, 16}, ExponentFunction -> (Null &)]], 
    {"." -> "_", " " -> ""}]; (* Removes any spaces that NumberForm might add *)
   
   
   
ClearAll[importAndReconstructArcs];
importAndReconstructArcs[arcFilePath_] := Module[
    {data, reconstructedData},

    (* Import the CSV data, skipping the header *)
    data = Import[arcFilePath, "CSV", "HeaderLines" -> 1];

    (* Reconstruct the doubles from the significand and exponent pairs *)
    reconstructedData = data /. {
        {axs_, axe_, ays_, aye_, azs_, aze_, bxs_, bxe_, bys_, bye_, bzs_, bze_, czs_, cze_} :> {
            axs * 2^axe, ays * 2^aye, azs * 2^aze,  (* Point A x, y, z *)
            bxs * 2^bxe, bys * 2^bye, bzs * 2^bze,  (* Point B x, y, z *)
            czs * 2^cze                            (* Const Z *)
        }
    };

    (* Return reconstructed data *)
    reconstructedData
];

ClearAll[importAndProcessBenchmarkData, reconstruct];

(* Helper function to reconstruct numbers from significand and exponent *)
reconstruct[significand_, exponent_] := significand * 2^exponent;

(* Main function to import and process CGAL results *)
importAndProcessCGALBenchmarkData[cgalResultsFile_] := Module[
    {cgalData},

    (* Import the data, skip the header, and reconstruct the numbers *)
    cgalData = Rest@Import[cgalResultsFile];

    (* Apply the reconstruct function to each set of values *)
    cgalData = Map[
        reconstruct @@@ Partition[#, 2] &,
        cgalData
    ];

    (* Return as an association for clarity *)
    <|"CGALResults" -> cgalData|>
]


ClearAll[analyzeErrors, computeBaselineResults, computeRelativeErrors];
 summarizeData[data_] := Max[data] 
(*summarizeData[data_] := Quantile[data, 0.99]*)

(* Helper function to compute baseline results *)
computeBaselineResults[arcData_, analysisPrecision_] := 
  Table[
    With[{pointA = arc[[1 ;; 3]], pointB = arc[[4 ;; 6]], constZ = arc[[7]]},
      gcaConstLatIntersectionCoordinatesNewEqn[pointA, pointB, constZ, analysisPrecision]
    ],
    {arc, arcData}
  ];

(* Helper function to compute relative errors *)
computeRelativeErrors[baselineResults_, benchmarkResults_, analysisPrecision_] := 
  Table[
    relativeError[
      baselineResults[[i, 1 ;; 2]], benchmarkResults[[i, 1 ;; 2]], 
      baselineResults[[i, 3]], analysisPrecision
    ],
    {i, Length[baselineResults]}
  ];

(* Main function to analyze errors *)
analyzeCGALErrors[arcData_, benchmarkData_,  analysisPrecision_] := 
  Module[{baselineResults, cgalFirstErrors, cgalSecondErrors, cgalErrors, meanCgalErrors},

    (* Compute baseline results *)
    baselineResults = computeBaselineResults[arcData, analysisPrecision];

    (* Compute first point relative errors for CGAL results *)
    cgalFirstErrors = computeRelativeErrors[baselineResults, benchmarkData["CGALResults"][[All, 1 ;; 2]], analysisPrecision];

    (* Compute second point relative errors for CGAL results *)
    cgalSecondErrors = computeRelativeErrors[baselineResults, benchmarkData["CGALResults"][[All, 3 ;; 4]], analysisPrecision];

    (* Compute the minimum relative error between first and second points *)
    cgalErrors = MapThread[Min, {cgalFirstErrors, cgalSecondErrors}];

    (* Compute the mean CGAL errors *)
    meanCgalErrors = summarizeData[cgalErrors];

    (* Combine results into an association *)
    <|"CGALMethodErrors" -> meanCgalErrors|>
  ]


analysisPrecision=100;


outputDirectory = absolutePath <> "mathematica_data/";

(* Create CSV data list holders *)
cgalErrorsList = {};


For[i = 1, i < Length[latitudePoints], i++,
  (* Format start and end latitudes *)
  startLatStr = formatLatitude[latitudePoints[[i]]];
  endLatStr = formatLatitude[latitudePoints[[i + 1]]];
  
  (* Generate filenames with precision *)
  latRangeStr = startLatStr <> "_" <> endLatStr;
  arcFilename = latRangeStr <> "ArcUp_Arcs_Exponent.csv";
  benchmarkCGALFilename = latRangeStr <> "_ArcUp_cgal_double.csv";
  (* Construct the full paths to the files and import data *)
  arcData = importAndReconstructArcs[absolutePath <> "generated_arcs/" <> arcFilename];
  
  benchmarkData = importAndProcessCGALBenchmarkData[absolutePath <> "benchmark_results/" <> benchmarkCGALFilename];
  errorAnalysisResults = analyzeCGALErrors[arcData, benchmarkData, analysisPrecision];
    (* Calculate midpoint of the latitude range *)
  latitudeRangeMidpoint = N[(latitudePoints[[i]] + latitudePoints[[i + 1]])/2, analysisPrecision];
  
  (* Append results to lists *)
  AppendTo[cgalErrorsList, {latitudeRangeMidpoint, errorAnalysisResults["CGALMethodErrors"]}];
];

(* Export results to CSV *)
Export[outputDirectory <> "CGAL_Relative_Error_ArcUp.csv", cgalErrorsList, "CSV"];


Print["Done CSV Writing"];
