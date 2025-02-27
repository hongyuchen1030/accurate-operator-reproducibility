(* ::Package:: *)

(* Check for necessary command line arguments *)
If[Length[$ScriptCommandLine] <= 2,
    Print["Usage: wolframscript -file this_package.m <mpfr_precisions> <latitudes>"];
    Return[$Failed];
];

(* Extract arguments from the command line *)
mpfrPrecisionString = $ScriptCommandLine[[2]];
latitudeString = $ScriptCommandLine[[3]];

(* Convert the comma-separated strings to lists *)
mpfrPrecisionList = ToExpression[StringSplit[mpfrPrecisionString, ","]];
latitudePoints = ToExpression[StringSplit[latitudeString, ","]];


ClearAll[absolutePath];
absolutePath = "/home/hyvchen/AccuracyBenchmarkEFT/"

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
  
 summarizeData[data_] := Mean[data] 
(*summarizeData[data_] := Quantile[data, 0.99]*)
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
    ToString[NumberForm[N[lat], {Infinity, 6}, ExponentFunction -> (Null &)]], 
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

ClearAll[importAndProcessBenchmarkData];
importAndProcessBenchmarkData[ourResultsFile_, float64ResultsFile_, quadrupleResultsFile_, mpfrFiles_] := Module[
    {
        ourData, float64Data,quadrupleData, mpfrData, reconstruct, mpfrPrecisionList, mpfrAssociation
    },
    
    (* Function to reconstruct numbers from significand and exponent *)
    reconstruct[significand_, exponent_] := significand * 2^exponent;

    (* Import and reconstruct our results data *)
    ourData = Import[ourResultsFile];
    ourData = ourData[[2 ;; All, All]]; (* Skip header *)
    ourData = Map[{reconstruct[#[[1]], #[[2]]], reconstruct[#[[3]], #[[4]]]} &, ourData];
    
    (* Import and reconstruct float64 results data *)
    float64Data = Import[float64ResultsFile];
    float64Data = float64Data[[2 ;; All, All]]; (* Skip header *)
    float64Data = Map[
        {
            reconstruct[#[[1]], #[[2]]], reconstruct[#[[3]], #[[4]]],
            reconstruct[#[[5]], #[[6]]], reconstruct[#[[7]], #[[8]]],
            reconstruct[#[[9]], #[[10]]], reconstruct[#[[11]], #[[12]]]
        } &,
        float64Data
    ];

    (* Import and reconstruct quadruple results data *)
    quadrupleData = Import[quadrupleResultsFile];
    quadrupleData = quadrupleData[[2 ;; All, All]]; (* Skip header *)
    quadrupleData = Map[
        {
            reconstruct[#[[1]], #[[2]]], reconstruct[#[[3]], #[[4]]],
            reconstruct[#[[5]], #[[6]]], reconstruct[#[[7]], #[[8]]],
            reconstruct[#[[9]], #[[10]]], reconstruct[#[[11]], #[[12]]]
        } &,
        quadrupleData
    ];
    
    (* Extract precision list from mpfrFiles names *)
    mpfrPrecisionList = StringCases[mpfrFiles, RegularExpression["mpfr_(\\d+)_"] :> ToExpression["$1"]];

    (* Import and reconstruct MPFR results for each precision level, each in its specific entry *)
    mpfrAssociation = Association @@ Table[
        "MPFR" <> IntegerString[prec] <> "Results" -> Map[
            {
                reconstruct[#[[1]], #[[2]]], reconstruct[#[[3]], #[[4]]],
                reconstruct[#[[5]], #[[6]]], reconstruct[#[[7]], #[[8]]],
                reconstruct[#[[9]], #[[10]]], reconstruct[#[[11]], #[[12]]]
            } &,
            Import[SelectFirst[mpfrFiles, StringContainsQ[#, "mpfr_" <> IntegerString[prec] <> "_"] &]][[2 ;; All]] (* Skip header and import data *)
        ],
        {prec, mpfrPrecisionList}
    ];
    
    (* Combine all data into an association *)
    Join[
        <|"OurResults" -> ourData, "Float64Results" -> float64Data, "QuadrupleResults" ->quadrupleData|>,
        mpfrAssociation
    ]
]


ClearAll[analyzeErrors];
analyzeErrors[arcData_, benchmarkData_, mpfrPrecisionList_, analysisPrecision_] := Module[
  {
    baselineResults, ourErrors, meanOurErrors,
    float64NewErrors, float64OldErrors,float64BaselineErrors, meanFloat64NewErrors, meanFloat64OldErrors,meanFloat64BaselineErrors,
    quadrupleNewErrors, quadrupleOldErrors, quadrupleBaselineErrors, meanQuadrupleNewErrors, meanQuadrupleOldErrors,  meanQuadrupleBaselineErrors,
    mpfrNewErrors, mpfrOldErrors, mpfrBaselineErrors,meanMPFRNewErrors, meanMPFROldErrors, meanMPFRBaselineErrors,results
  },

  baselineResults = {};
  baselineResults = Table[
    With[
      {pointA = arc[[1;;3]],
       pointB = arc[[4;;6]],
       constZ = arc[[7]]},
      gcaConstLatIntersectionCoordinatesNewEqn[pointA, pointB, constZ, analysisPrecision]
    ],
    {arc, arcData}
  ];
 (* Print the first 5 rows of each relevant dataset *)
(* Print["baseline result (First 5 rows): ", baselineResults[[1 ;; 5]]];
Print["MPFR16 result (First 5 rows): ", N[benchmarkData["MPFR16Results"][[1 ;; 5]]]];
Print["MPFR18 result (First 5 rows): ", N[benchmarkData["MPFR18Results"][[1 ;; 5]]]];*)
  

  (* Compute relative errors for our method *)
  ourErrors = Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["OurResults"][[i]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[arcData]}
  ];
  meanOurErrors = summarizeData[ourErrors];
 
  
  (* Compute relative errors for Float64 results *)
  float64NewErrors = Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["Float64Results"][[i, 1 ;; 2]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[arcData]}
  ];
  float64OldErrors = Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["Float64Results"][[i, 3 ;; 4]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[arcData]}
  ];
   float64BaselineErrors = Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["Float64Results"][[i, 5 ;; 6]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[arcData]}
  ];
  meanFloat64NewErrors = summarizeData[float64NewErrors];
  meanFloat64OldErrors = summarizeData[float64OldErrors];
  meanFloat64BaselineErrors = summarizeData[float64BaselineErrors];

  (* Compute relative errors for Quadruple results *)
  quadrupleNewErrors = Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["QuadrupleResults"][[i, 1 ;; 2]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[arcData]}
  ];
  quadrupleOldErrors = Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["QuadrupleResults"][[i, 3 ;; 4]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[arcData]}
  ];
  
   quadrupleBaselineErrors = Table[
    relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["QuadrupleResults"][[i, 5 ;; 6]], baselineResults[[i, 3]], analysisPrecision],
    {i, Length[arcData]}
  ];
  meanQuadrupleNewErrors = summarizeData[quadrupleNewErrors];
  meanQuadrupleOldErrors = summarizeData[quadrupleOldErrors];
  meanQuadrupleBaselineErrors = summarizeData[quadrupleBaselineErrors];
  
  (* Initialize dictionaries for MPFR results *)
  meanMPFRNewErrors = Association[];
  meanMPFROldErrors = Association[];
  meanMPFRBaselineErrors = Association[];

  (* Compute relative errors for MPFR results for each precision *)
  Do[
    mpfrNewErrors = Table[
      relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["MPFR" <> IntegerString[prec] <> "Results"][[i, 1 ;; 2]], baselineResults[[i, 3]], analysisPrecision],
      {i, Length[arcData]}
    ];
    mpfrOldErrors = Table[
      relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["MPFR" <> IntegerString[prec] <> "Results"][[i, 3 ;; 4]], baselineResults[[i, 3]], analysisPrecision],
      {i, Length[arcData]}
    ];
    mpfrBaselineErrors = Table[
      relativeError[baselineResults[[i, 1 ;; 2]], benchmarkData["MPFR" <> IntegerString[prec] <> "Results"][[i, 5 ;; 6]], baselineResults[[i, 3]], analysisPrecision],
      {i, Length[arcData]}
    ];
    meanMPFRNewErrors["MPFR" <> IntegerString[prec] <> "NewMethodErrors"] = summarizeData[mpfrNewErrors];
    meanMPFROldErrors["MPFR" <> IntegerString[prec] <> "OldMethodErrors"] = summarizeData[mpfrOldErrors];
    meanMPFRBaselineErrors["MPFR" <> IntegerString[prec] <> "BaselineMethodErrors"] = summarizeData[mpfrBaselineErrors];
    ,
    {prec, mpfrPrecisionList}
  ];
  
    (* Combine all results into a single association *)
  results = Join[
    <|
      "OurMethodErrors" -> meanOurErrors,
      "Float64NewMethodErrors" -> meanFloat64NewErrors,
      "Float64OldMethodErrors" -> meanFloat64OldErrors,
      "Float64BaselineMethodErrors" -> meanFloat64BaselineErrors,
      "QuadrupleNewMethodErrors" -> meanQuadrupleNewErrors,
      "QuadrupleOldMethodErrors" -> meanQuadrupleOldErrors,
      "QuadrupleBaselineMethodErrors" -> meanQuadrupleBaselineErrors
    |>,
    meanMPFRNewErrors,
    meanMPFROldErrors,
    meanMPFRBaselineErrors
  ];

  results
]


analysisPrecision=100;


outputDirectory = absolutePath <> "mathematica_data/";

(* Create CSV data list holders *)
ourErrorsList = {};
float64ErrorsList = {};
quadrupleErrorsList = {};
mpfrErrorsLists = AssociationMap[Function[prec, {}], mpfrPrecisionList];

For[i = 1, i < Length[latitudePoints], i++,
  (* Format start and end latitudes *)
  startLatStr = formatLatitude[latitudePoints[[i]]];
  endLatStr = formatLatitude[latitudePoints[[i + 1]]];
  
  (* Generate filenames with precision *)
  latRangeStr = startLatStr <> "_" <> endLatStr;
  arcFilename = latRangeStr <> "Arcs_Exponent.csv";
  benchmarkOurFilename = latRangeStr <> "_our_results_double.csv";
  benchmarkDoubleFilename = latRangeStr <> "_float64_results_double.csv";
  benchmarkQuadrupleFilename = latRangeStr <> "_quadruple_results_double.csv";
  benchmarkMPFRFilenames = Table[absolutePath <> "benchmark_results/" <> latRangeStr <> "_mpfr_" <> IntegerString[mpfrPrecision] <> "_results.csv", {mpfrPrecision, mpfrPrecisionList}];
  
  (* Construct the full paths to the files and import data *)
  arcData = importAndReconstructArcs[absolutePath <> "generated_arcs/" <> arcFilename];
  
  benchmarkData = importAndProcessBenchmarkData[absolutePath <> "benchmark_results/" <> benchmarkOurFilename, absolutePath <> "benchmark_results/" <> benchmarkDoubleFilename, absolutePath <> "benchmark_results/" <> benchmarkQuadrupleFilename, benchmarkMPFRFilenames];
  errorAnalysisResults = analyzeErrors[arcData, benchmarkData, mpfrPrecisionList, analysisPrecision];
    (* Calculate midpoint of the latitude range *)
  latitudeRangeMidpoint = N[(latitudePoints[[i]] + latitudePoints[[i + 1]])/2, analysisPrecision];
  
  (* Append results to lists *)
  AppendTo[ourErrorsList, {latitudeRangeMidpoint, errorAnalysisResults["OurMethodErrors"]}];
  AppendTo[float64ErrorsList, {latitudeRangeMidpoint, errorAnalysisResults["Float64NewMethodErrors"], errorAnalysisResults["Float64OldMethodErrors"], errorAnalysisResults["Float64BaselineMethodErrors"]}];
  AppendTo[quadrupleErrorsList, {latitudeRangeMidpoint, errorAnalysisResults["QuadrupleNewMethodErrors"], errorAnalysisResults["QuadrupleOldMethodErrors"], errorAnalysisResults["QuadrupleBaselineMethodErrors"]}];
  Do[
    AppendTo[mpfrErrorsLists[prec], {latitudeRangeMidpoint, errorAnalysisResults["MPFR" <> IntegerString[prec] <> "NewMethodErrors"], errorAnalysisResults["MPFR" <> IntegerString[prec] <> "OldMethodErrors"], errorAnalysisResults["MPFR" <> IntegerString[prec] <> "BaselineMethodErrors"]}],
    {prec, mpfrPrecisionList}
  ];
];

(* Export results to CSV *)
Export[outputDirectory <> "Our_Relative_Error.csv", ourErrorsList, "CSV"];
Export[outputDirectory <> "Float64_Relative_Error.csv", float64ErrorsList, "CSV"];
Export[outputDirectory <> "Quadruple_Relative_Error.csv", quadrupleErrorsList, "CSV"];
Do[
  Export[outputDirectory <> "MPFR" <> IntegerString[prec] <> "_Relative_Error.csv", mpfrErrorsLists[prec], "CSV"],
  {prec, mpfrPrecisionList}
];

Print["Done CSV Writing"];
