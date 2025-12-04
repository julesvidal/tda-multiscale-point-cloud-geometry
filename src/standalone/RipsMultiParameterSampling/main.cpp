/// TODO 12: Add your information
/// \author Your Name Here <Your Email Address Here>.
/// \date The Date Here.
///
/// \brief dummy program example.

// TTK Includes
#include <CommandLineParser.h>
#include <filesystem>
#include <string>
#include <ttkRipsMultiParameterSampling.h>

// VTK Includes
#include <vtkPLYReader.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLDataObjectWriter.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLGenericDataObjectReader.h>

int main(int argc, char **argv) {

  // ---------------------------------------------------------------------------
  // Program variables
  // ---------------------------------------------------------------------------
  std::vector<std::string> inputFilePaths;
  std::vector<std::string> inputArrayNames;
  std::string outputPathPrefix{"output"};

  std::string normalArrayName{"Normals"};

  std::string ratios{""};
  // int nbcc{10};
  int sizeMin{-1};

  int WriteTimings{false};
  std::string TimingFilePrefix{"timings"};

  // sampling
  // int useSampling{true};
  double sampling_min{0.01};
  double sampling_max{0.1};
  int sampling_count{50};
  double sampling_factor{0.1};
  SimplexId NbPointsMax{100000};

  // metrics option
  int FitType{1}; // 0: plane, 1: sphere, 2: quadric
  int FitDimension{4};
  int UseKnnGraph{false};
  int UseRIMLS{false};
  // int UseASO{true};
  // int UseNormalizedCoords{true};
  // int UseProjectiveNorm{false};

  // approx parameters
  double MaxEdgeLength{0.01};
  double DensitySigma{0.01};


  // epsilon sampling
  double Epsilon{1e-5};
  double EpsilonNbMaxPoint{10000};

  // output parameter
  double ShrinkingDiagMax{0.1};
  int NbPrimitives{100};
  int GrowPrimitives{1};
  int GrowPrimitiveIterations{10};
  float GrowPrimitiveThreshold{0.1};

  double ProductPersistenceMin{0.1};


  // double min_jaccard{0.7};
  // double gridMode{1};

  int nb_primitives{100};
  // int outputGrid{false};


  bool listArrays{false};

  // ---------------------------------------------------------------------------
  // Set program variables based on command line arguments
  // ---------------------------------------------------------------------------
  {
    ttk::CommandLineParser parser;

    // -------------------------------------------------------------------------
    // Standard options and arguments
    // -------------------------------------------------------------------------
    parser.setArgument(
      "i", &inputFilePaths, "Input data-sets (*.vti, *vtu, *vtp)", false);
    parser.setOption("l", &listArrays, "List available arrays");
    parser.setArgument(
      "o", &outputPathPrefix, "outputPathPrefix", true);
    // parser.setArgument(
    //   "o", &outputPathPrefix, "Output file prefix (no extension)", true);

    // -------------------------------------------------------------------------
    // TODO 13: Declare custom arguments and options
    // -------------------------------------------------------------------------
    // parser.setArgument("O", &outputArrayName, "Output array name", true);
    parser.setArgument("n", &normalArrayName, "normal array name", true);

    // parser.setArgument("S", &useSampling, "use sampling", false);
    parser.setArgument("N", &NbPointsMax, "max nb points", true);
    parser.setArgument("Smin", &sampling_min, "min sampling", true);
    parser.setArgument("Smax", &sampling_max, "max sampling", true);
    parser.setArgument("Scount", &sampling_count, "count sampling", true);
    // parser.setArgument("Sfactor", &sampling_factor, "sampling factor", true);
    parser.setArgument("knn", &UseKnnGraph, "use knn-graph", true);
    parser.setArgument("rimls", &UseRIMLS, "use RIMLS", true);


    parser.setArgument("lmax", &MaxEdgeLength, "max edge length", true);
    parser.setArgument("e", &Epsilon, "epsilon sampling parameter", true);
    parser.setArgument("m", &EpsilonNbMaxPoint, "epsilon sampling max point", true);
    parser.setArgument("S", &sizeMin, "Size min", true);

    parser.setArgument("s", &DensitySigma, "density sigma", true);

    parser.setArgument("fitDimension", &FitDimension, "Fit Dimension", true);
    parser.setArgument("fitType", &FitType, "Fit Type (0: plane, 1: sphere, 2: quadric) ", true);

    parser.setArgument("P", &ProductPersistenceMin, "Product Persistence Min Threshold (relative to highest)", true);

    parser.setArgument("tm", &WriteTimings, "write timings", true);
    parser.setArgument("tmf", &TimingFilePrefix, "TimingFilePrefix", true);
    // parser.setArgument("g", &outputGrid, "output grid", false);
    // parser.setArgument("R", &ratios, "Ratios (comma-separated list of values)", false);
    // parser.setArgument("G", &gridMode, "grid mode (0: compact, 1: full)", false);
    // parser.setArgument("J", &min_jaccard, "Jaccard constraint", true);
    // parser.setArgument("N", &nbcc, "Max number of components", false);
    // parser.setArgument("P", &nb_primitives, "output max number of primitives", true);

    parser.parse(argc, argv);
  }

  // ---------------------------------------------------------------------------
  // Command line output messages.
  // ---------------------------------------------------------------------------
  ttk::Debug msg;
  msg.setDebugMsgPrefix("RipsMultiParameterSampling");

  // ---------------------------------------------------------------------------
  // Initialize ttkRipsMultiParameterSampling module (adjust parameters)
  // ---------------------------------------------------------------------------
  auto ripsMultiParameterSampling = vtkSmartPointer<ttkRipsMultiParameterSampling>::New();

  // ---------------------------------------------------------------------------
  // TODO 14: Pass custom arguments and options to the module
  // ---------------------------------------------------------------------------
  // ripsMultiParameterSampling->SetOutputArrayName(outputArrayName);
  // ripsMultiParameterSampling->SetRatioQuery(ratios.data());
    ripsMultiParameterSampling->SetInputArrayToProcess(0, 0, 0, 0, normalArrayName.data());
    ripsMultiParameterSampling->SetScaleSamplingMinMaxCount(sampling_min, sampling_max, sampling_count);
    ripsMultiParameterSampling->SetNbPointsMax(NbPointsMax);
    ripsMultiParameterSampling->SetUseKnnGraph(UseKnnGraph);
    ripsMultiParameterSampling->SetUseRIMLS(UseRIMLS);
    ripsMultiParameterSampling->SetSizeMin(sizeMin);
    ripsMultiParameterSampling->SetWriteTimings(WriteTimings);
    ripsMultiParameterSampling->SetTimingFilePrefix(TimingFilePrefix);


    ripsMultiParameterSampling->SetDiagMax(ShrinkingDiagMax);
    ripsMultiParameterSampling->SetNbPrimitives(NbPrimitives);

    ripsMultiParameterSampling->SetMaxEdgeLength(MaxEdgeLength);
    ripsMultiParameterSampling->SetDensitySigma(DensitySigma);
    ripsMultiParameterSampling->SetEpsilon(Epsilon);
    ripsMultiParameterSampling->SetEpsilonNbMaxPoint(EpsilonNbMaxPoint);

    ripsMultiParameterSampling->SetGrowPrimitives(GrowPrimitives);
    ripsMultiParameterSampling->SetGrowPrimitiveIterations(GrowPrimitiveIterations);
    ripsMultiParameterSampling->SetGrowPrimitiveThreshold(GrowPrimitiveThreshold);
    ripsMultiParameterSampling->SetProductPersistenceMin(ProductPersistenceMin);

    ripsMultiParameterSampling->SetFitDimension(FitDimension);
    ripsMultiParameterSampling->SetFitType(FitType);



  // ---------------------------------------------------------------------------
  // Read input vtkDataObjects (optionally: print available arrays)
  // ---------------------------------------------------------------------------
  vtkDataArray *defaultArray = nullptr;
  for(size_t i = 0; i < inputFilePaths.size(); i++) {

    std::filesystem::path fsPath(inputFilePaths[i]);
    vtkDataObject*  inputDataObject{};

    vtkNew<vtkPLYReader> reader_ply;
    vtkNew<vtkXMLGenericDataObjectReader> reader;
    if(fsPath.extension()==".ply"){
      reader_ply->SetFileName(inputFilePaths[i].data());
      reader_ply->Update();
      inputDataObject = reader_ply->GetOutput();
    }else{
      reader->SetFileName(inputFilePaths[i].data());
      reader->Update();
      inputDataObject = reader->GetOutput();
    }

    // check if input vtkDataObject was successfully read
    if(!inputDataObject) {
      msg.printErr("Unable to read input file `" + inputFilePaths[i] + "' :(");
      return 1;
    }

    auto inputAsVtkDataSet = vtkDataSet::SafeDownCast(inputDataObject);

    // if requested print list of arrays, otherwise proceed with execution
    if(listArrays) {
      msg.printMsg(inputFilePaths[i] + ":");
      if(inputAsVtkDataSet) {
        // Point Data
        msg.printMsg("  PointData:");
        auto pointData = inputAsVtkDataSet->GetPointData();
        for(int j = 0; j < pointData->GetNumberOfArrays(); j++)
          msg.printMsg("    - " + std::string(pointData->GetArrayName(j)));

        // Cell Data
        msg.printMsg("  CellData:");
        auto cellData = inputAsVtkDataSet->GetCellData();
        for(int j = 0; j < cellData->GetNumberOfArrays(); j++)
          msg.printMsg("    - " + std::string(cellData->GetArrayName(j)));
      } else {
        msg.printErr("Unable to list arrays on file `" + inputFilePaths[i]
                     + "'");
        return 1;
      }
    } else {
      // feed input object to ttkRipsMultiParameterSampling filter
      ripsMultiParameterSampling->SetInputDataObject(i, inputAsVtkDataSet);

      // default arrays
      if(!defaultArray) {
        defaultArray = inputAsVtkDataSet->GetPointData()->GetArray(0);
        if(!defaultArray)
          defaultArray = inputAsVtkDataSet->GetCellData()->GetArray(0);
      }
    }
  }

  // terminate program if it was just asked to list arrays
  if(listArrays) {
    return 0;
  }

  // ---------------------------------------------------------------------------
  // Specify which arrays of the input vtkDataObjects will be processed
  // ---------------------------------------------------------------------------
  if(!inputArrayNames.size()) {
    if(defaultArray)
      inputArrayNames.emplace_back(defaultArray->GetName());
  }
  // for(size_t i = 0; i < inputArrayNames.size(); i++)
  //   ripsMultiParameterSampling->SetInputArrayToProcess(i, 0, 0, 0, inputArrayNames[i].data());

  // ---------------------------------------------------------------------------
  // Execute ttkRipsMultiParameterSampling filter
  // ---------------------------------------------------------------------------
  ripsMultiParameterSampling->Update();

  // ---------------------------------------------------------------------------
  // If output prefix is specified then write all output objects to disk
  // ---------------------------------------------------------------------------
  if(!outputPathPrefix.empty()) {
    std::string suffix[2] = {"primitives", "diagram"};
    int ports[2] = {5, 2};
    for(int i = 0; i < 2; i++) {
      auto output = ripsMultiParameterSampling->GetOutputDataObject(ports[i]);
      // auto writer = vtkSmartPointer<vtkXMLWriter>::Take(
      //   vtkXMLMultiBlockDataWriter::NewWriter(output->GetDataObjectType()));
      auto writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();

      std::string outputFileName = outputPathPrefix + "_" + suffix[i] + "."
                                   + writer->GetDefaultFileExtension();
      msg.printMsg("Writing output file `" + outputFileName + "'...");
      writer->SetInputDataObject(output);
      writer->SetFileName(outputFileName.data());
      writer->Update();
    }
  }

  return 0;
}
