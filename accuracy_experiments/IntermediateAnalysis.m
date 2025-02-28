(* ::Package:: *)

(* Check for necessary command line arguments *)
If[Length[$ScriptCommandLine] <= 3,
Print["Usage: wolframscript -file this_package.m <mpfr_precisions> <latitudes><QueryRangeList>"];
Return[$Failed];
];

(* Extract arguments from the command line *)
mpfrPrecisionString = $ScriptCommandLine[[2]];
latitudeString = $ScriptCommandLine[[3]];
queryRangeString = $ScriptCommandLine[[4]];

Print[mpfrPrecisionString]
Print[latitudeString]
Print[queryRangeString]

(* Convert the comma-separated strings to lists *)
mpfrPrecisionList = ToExpression[StringSplit[mpfrPrecisionString, ","]];
latitudePoints = ToExpression[StringSplit[latitudeString, ","]];
(* Parse the QueryRange string properly *)
queryRangeString = StringReplace[queryRangeString, {"{" -> "{", "}" -> "}"}];
QueryRange = ToExpression["{" <> queryRangeString <> "}"];

Print[QueryRange]

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


ClearAll[absolutePath];
absolutePath = DirectoryName[Directory[]]

(* Define the baseline function *)
Clear[gcaConstLatIntersectionCoordinatesNewEqnIntermediate];
gcaConstLatIntersectionCoordinatesNewEqnIntermediate[pointA_List, pointB_List, constZ_, precision_: MachinePrecision] := Module[
  {n, nx, ny, nz, nxSquared, nySquared, nzSquared, normNSquared, sTilde, sSquare, znxnz, sny, znxnzsny, nxSquarenySquare, px},

  n = Cross[pointA, pointB];
  {nx, ny, nz} = n;
  nxSquared = nx^2;
  nySquared = ny^2;
  nzSquared = nz^2;
  nxSquarenySquare = nxSquared + nySquared;
  normNSquared = nxSquared + nySquared + nzSquared;
  
  sTilde = Sqrt[nxSquared + nySquared - normNSquared constZ^2];
  sSquare = nxSquared + nySquared - normNSquared constZ^2;
  
  znxnz = constZ nx nz;
  sny = sTilde ny;
  znxnzsny = znxnz + sny;
  
  px = -(1/(nxSquared + nySquared)) * znxnzsny;
  
  (* Set precision for all outputs *)
  {nx, ny, nz, nxSquared, nySquared, nzSquared, normNSquared, sTilde, sSquare, znxnz, sny, znxnzsny, nxSquarenySquare, px} = 
    SetPrecision[{nx, ny, nz, nxSquared, nySquared, nzSquared, normNSquared, sTilde, sSquare, znxnz, sny, znxnzsny, nxSquarenySquare, px}, precision];
  
  (* Return results as an association *)
  <|
    "nx" -> nx,
    "normNSquared" -> normNSquared,
    "sSquare" -> sSquare,
    "s" -> sTilde,
    "znxnz" -> znxnz,
    "sny" -> sny,
    "znxnzsny" -> znxnzsny,
    "nxSquarenySquare" -> nxSquarenySquare,
    "px" -> px
  |>
];

(* Define a function to read the CSV and reconstruct the values *)
ClearAll[ImportAndReconstructValues];
ImportAndReconstructValues[filePath_String, precision_] := Module[
    {data, headers, rows, reconstructedData},
    
    (* Import the CSV file, specifying that it has headers *)
    data = Import[filePath, {"CSV", "Data"}];
    headers = First[data];
    rows = Rest[data];

    (* Define a function to reconstruct a number from its significand and exponent *)
    reconstruct[significand_, exponent_] := significand * 2^exponent;

    (* Map the reconstruction function across all relevant columns for each row *)
    reconstructedData = rows /. {
        {a_, b_, c_, d_, e_, f_, g_, h_, i_, j_, k_, l_, m_, n_, o_, p_, q_, r_, s_, t_, u_, v_, w_, x_, y_, z_, aa_, ab_, ac_, ad_} :> <|
            "nx" -> SetPrecision[reconstruct[a, b] + reconstruct[c, d], precision],
            "normNSquared" -> SetPrecision[reconstruct[e, f] + reconstruct[g, h], precision],
            "sSquare" -> SetPrecision[reconstruct[i, j] + reconstruct[k, l], precision],
            "s" -> SetPrecision[reconstruct[m, n] + reconstruct[o, p], precision],
            "znxnz" -> SetPrecision[reconstruct[q, r] + reconstruct[s, t], precision],
            "sny" -> SetPrecision[reconstruct[u, v] + reconstruct[w, x], precision],
            "znxnzsny" -> SetPrecision[reconstruct[y, z], precision],
            "nxSquarenySquare" -> SetPrecision[reconstruct[aa, ab], precision],
            "px" -> SetPrecision[reconstruct[ac, ad], precision]
        |>
    };
    (* Return the reconstructed data *)
    reconstructedData
]

(* Define a function to import and construct MPFR values from a CSV file *)
ClearAll[importAndConstructMPFRValues];
importAndConstructMPFRValues[filePath_String, precision_] := Module[
    {data, headers, rows, records},
    
    (* Import the CSV file, specifying that it has headers *)
    data = Import[filePath, {"CSV", "Data"}];
    headers = First[data];
    rows = Rest[data];
    
    (* Define the header mapping *)
    headerMapping = <|
        "nx" -> "nx",
        "normNSquared" -> "nNorm",
        "sSquare" -> "sSquare",
        "s" -> "s",
        "znxnz" -> "znxnz",
        "sny" -> "sny",
        "znxnzsny" -> "znxnzsny",
        "nxSquarenySquare" -> "nxSquarenySquare",
        "px" -> "px"
    |>;
    
    (* Construct associations for each row with specified precision and new headers *)
    records = Map[
        AssociationThread[Keys[headerMapping], SetPrecision[#, precision] & /@ Values[headerMapping /. AssociationThread[headers, #]]] &,
        rows
    ];
   
    
    (* Return the constructed MPFR values *)
    records
]



(* Define a function to read the CSV and reconstruct the values *)
ClearAll[ImportAndReconstructDoubleValues];
ImportAndReconstructDoubleValues[filePath_String, precision_] := Module[
    {data, headers, rows, reconstructedData},
    
    (* Import the CSV file, specifying that it has headers *)
    data = Import[filePath, {"CSV", "Data"}];
    headers = First[data];
    rows = Rest[data];
    
    (* Define a function to reconstruct a number from its significand and exponent *)
    reconstruct[significand_, exponent_] := significand * 2^exponent;
    
    (* Map the reconstruction function across all relevant columns for each row *)
    reconstructedData = rows /. {
        {nxSig_, nxE_, nNormSig_, nNormE_,
         sSquareSig_, sSquareE_, sSig_, sE_, znxnzSig_, znxnzE_, 
         snySig_, snyE_, znxnzsnySig_, znxnzsnyE_, nxSquarenySquareSig_, nxSquarenySquareE_, pxSig_, pxE_} :> <|
            "nx" -> SetPrecision[reconstruct[nxSig, nxE], precision],
            "normNSquared" -> SetPrecision[reconstruct[nNormSig, nNormE], precision],
            "sSquare" -> SetPrecision[reconstruct[sSquareSig, sSquareE], precision],
            "s" -> SetPrecision[reconstruct[sSig, sE], precision],
            "znxnz" -> SetPrecision[reconstruct[znxnzSig, znxnzE], precision],
            "sny" -> SetPrecision[reconstruct[snySig, snyE], precision],
            "znxnzsny" -> SetPrecision[reconstruct[znxnzsnySig, znxnzsnyE], precision],
            "nxSquarenySquare" -> SetPrecision[reconstruct[nxSquarenySquareSig, nxSquarenySquareE], precision],
            "px" -> SetPrecision[reconstruct[pxSig, pxE], precision]
        |>
    };

    (* Return the reconstructed data *)
    reconstructedData
]



(* Function to compute relative errors *)
computeRelativeErrors[ourResult_, baselineResults_] := Module[
    {relativeErrors},
    relativeErrors = Table[
        KeyValueMap[
            (#1 -> If[KeyExistsQ[baselineResults[[i]], #1],
                      Abs[(ourResult[[i]][#1] - baselineResults[[i]][#1]) / baselineResults[[i]][#1]],
                      Missing["KeyAbsent", #1]
                     ]) &,
            ourResult[[i]]
        ],
        {i, Length[ourResult]}
    ];
    relativeErrors
]

(* Function to compute averages *)
computeAverages[data_List] := Module[
    {keys, values},
    keys = Keys[data[[1]]];
    values = Table[
        key -> Mean[DeleteCases[Lookup[data, key], _Missing]],
        {key, keys}
    ];
    Association[values]
]

(* Function to filter data based on a latitude range *)
filterDataByLatitude[data_List, latitudes_List, range_List] := Module[
    {filteredIndices},
    filteredIndices = Flatten[Position[latitudes, _?(range[[1]] <= # <= range[[2]] &)]];
    data[[filteredIndices]]
]

  (* Compute averages for each query range *)
averageResultsForQueryRange[data_List, latitudes_List, queryRange_List] := Module[
    {filteredData},
    filteredData = filterDataByLatitude[data, latitudes, queryRange];
    computeAverages[filteredData]
]


(* Convert a single value to scientific notation with 'E' *)
toScientificString[value_, precision_] := ScientificForm[value, precision, NumberFormat -> (#1 <> "e" <> #3 &)] // ToString;

createCSVRow[method_, dataList_] := Prepend[
  Map[toScientificString[#, analysisDesiredPrecision] &, dataList],
  method
];


analysisDesiredPrecision = 64
(* Define the filename for the CSV output *)
outputCSVFile = FileNameJoin[{absolutePath, "/mathematica_data/error_analysis_intermediate_results_equator.csv"}];

(* Collect latitude data to filter later *)
latitudes = Table[(latitudePoints[[j]] + latitudePoints[[j + 1]])/2, {j, Length[latitudePoints] - 1}];

csvLabels = {"range", "nx", "normNSquared", "sSquare", "s", "znxnz", "sny", "znxnzsny", "nxSquarenySquareDenom", "px"};
allRelativeErrorsAvgOur = {};
allRelativeErrorsAvgDouble = {};
allRelativeErrorsAvgQuadruple = {};

  For[j = 1, j < Length[latitudePoints], j++,
    (* Format start and end latitudes *)
    startLatStr = formatLatitude[latitudePoints[[j]]];
    endLatStr = formatLatitude[latitudePoints[[j + 1]]];
  
    (* Generate filenames with precision *)
    latRangeStr = startLatStr <> "_" <> endLatStr;
    fileNamePrefix = "intermediate_results/" <> latRangeStr;
    ourData = ImportAndReconstructValues[absolutePath <> fileNamePrefix <> "_intermediate_results_our.csv", analysisDesiredPrecision];
    doubleData = ImportAndReconstructDoubleValues[absolutePath <> fileNamePrefix <> "_intermediate_results_double.csv", analysisDesiredPrecision];
    quadrupleData =  importAndConstructMPFRValues[absolutePath <> fileNamePrefix <> "_intermediate_results_quadruple.csv", analysisDesiredPrecision];
  
    arcFilename = latRangeStr <> "Arcs_Exponent.csv";
    (* Construct the full paths to the files and import data *)
    arcData = importAndReconstructArcs[absolutePath <> "generated_arcs/" <> arcFilename];
  
    (* Generate baseline results from arcs data *)
    baselineResults = Table[
      With[
        {pointA = arc[[1 ;; 3]], pointB = arc[[4 ;; 6]], constZ = arc[[7]]},
        gcaConstLatIntersectionCoordinatesNewEqnIntermediate[pointA, pointB, constZ, analysisDesiredPrecision]
      ],
      {arc, arcData}
    ];
  
    relativeErrorsOurs = computeRelativeErrors[ourData, baselineResults];
    relativeErrorsDouble = computeRelativeErrors[doubleData, baselineResults];
    relativeErrorsQuadruple = computeRelativeErrors[quadrupleData, baselineResults];
  
    relativeErrorsAvgOur = computeAverages[relativeErrorsOurs];
    relativeErrorsAvgDouble = computeAverages[relativeErrorsDouble];
    relativeErrorsAvgQuadruple = computeAverages[relativeErrorsQuadruple];
    (* Collect the averages for each iteration *)
    AppendTo[allRelativeErrorsAvgOur, relativeErrorsAvgOur];
    AppendTo[allRelativeErrorsAvgDouble, relativeErrorsAvgDouble];
    AppendTo[allRelativeErrorsAvgQuadruple, relativeErrorsAvgQuadruple];
  ];
  

	finalAveragesOur = Table[
	    Join[<|"range" -> ToString[range]|>, averageResultsForQueryRange[allRelativeErrorsAvgOur, latitudes, range]],
	    {range, QueryRange}
	];
	finalAveragesDouble = Table[
	    Join[<|"range" -> ToString[range]|>, averageResultsForQueryRange[allRelativeErrorsAvgDouble, latitudes, range]],
	    {range, QueryRange}
	];
	finalAveragesQuadruple = Table[
	    Join[<|"range" -> ToString[range]|>, averageResultsForQueryRange[allRelativeErrorsAvgQuadruple, latitudes, range]],
	    {range, QueryRange}
	];
	
  (* Convert final averages to a list of lists suitable for CSV export *)
	finalAveragesListOur = Prepend[
	    Table[
	        Values[finalAveragesOur[[i]]],
	        {i, Length[finalAveragesOur]}
	    ],
	    csvLabels
	];
	  (* Convert final averages to a list of lists suitable for CSV export *)
	finalAveragesListDouble = Prepend[
	    Table[
	        Values[finalAveragesDouble[[i]]],
	        {i, Length[finalAveragesDouble]}
	    ],
	    csvLabels
	];
	
		  (* Convert final averages to a list of lists suitable for CSV export *)
	finalAveragesListQuadruple = Prepend[
	    Table[
	        Values[finalAveragesQuadruple[[i]]],
	        {i, Length[finalAveragesQuadruple]}
	    ],
	    csvLabels
	];

	(* Define the output file path *)
	outputFilePathOur = absolutePath <> "mathematica_data/" <>queryRangeString<> "_intermediate_analysis_our.csv";
	outputFilePathDouble = absolutePath <> "mathematica_data/"  <>queryRangeString<> "_intermediate_analysis_double.csv";
	outputFilePathQuadruple = absolutePath <> "mathematica_data/"  <>queryRangeString<> "_intermediate_analysis_quadruple.csv";
	
	(* Export the final averages to CSV *)
	Export[outputFilePathOur, finalAveragesListOur];
	Export[outputFilePathDouble, finalAveragesListDouble];
	Export[outputFilePathQuadruple, finalAveragesListQuadruple];
  
  

allRelativeErrorsAvgOneMPFR = {};  
For[i = 1, i <= Length[mpfrPrecisionList], i++,
  precision = mpfrPrecisionList[[i]];
  allRelativeErrorsAvgOneMPFR = {};

  For[j = 1, j < Length[latitudePoints], j++,
    (* Format start and end latitudes *)
    startLatStr = formatLatitude[latitudePoints[[j]]];
    endLatStr = formatLatitude[latitudePoints[[j + 1]]];
  
    (* Generate filenames with precision *)
    latRangeStr = startLatStr <> "_" <> endLatStr;
    fileNamePrefix = "intermediate_results/" <> latRangeStr;
    arcFilename = latRangeStr <> "Arcs_Exponent.csv";
    (* Construct the full paths to the files and import data *)
    arcData = importAndReconstructArcs[absolutePath <> "generated_arcs/" <> arcFilename];
  
    (* Generate baseline results from arcs data *)
    baselineResults = Table[
      With[
        {pointA = arc[[1 ;; 3]], pointB = arc[[4 ;; 6]], constZ = arc[[7]]},
        gcaConstLatIntersectionCoordinatesNewEqnIntermediate[pointA, pointB, constZ, analysisDesiredPrecision]
      ],
      {arc, arcData}
    ];
  
    (* Generate filename based on precision *)
    mpfrFileName = absolutePath <> fileNamePrefix <> "_intermediate_result_MPFR_" <> ToString[precision] <> ".csv";
  
    (* Read and process MPFR CSV data *)
    mpfrData = importAndConstructMPFRValues[mpfrFileName, analysisDesiredPrecision];
  
    relativeErrorMPFR = computeRelativeErrors[mpfrData, baselineResults];
    AppendTo[allRelativeErrorsAvgOneMPFR, computeAverages[relativeErrorMPFR]];
  ];
  
	finalAveragesOneMPFR = Table[
	    Join[<|"range" -> ToString[range]|>, averageResultsForQueryRange[allRelativeErrorsAvgOneMPFR, latitudes, range]],
	    {range, QueryRange}
	];
	
	(* Convert final averages to a list of lists suitable for CSV export *)
	finalAveragesListOneMPFR = Prepend[
	    Table[
	        Values[finalAveragesOneMPFR[[i]]],
	        {i, Length[finalAveragesOneMPFR]}
	    ],
	    csvLabels
	];
	
	outputFilePathOneMPFR  = absolutePath <> "mathematica_data/" <>queryRangeString<> "_intermediate_result_MPFR_" <> ToString[precision] <> ".csv";
	Export[outputFilePathOneMPFR, finalAveragesListOneMPFR];
]

Print["Done Intermediate results"];

