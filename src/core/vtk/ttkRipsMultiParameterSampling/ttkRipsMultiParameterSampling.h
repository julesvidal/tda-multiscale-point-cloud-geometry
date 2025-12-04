/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
/// \ingroup vtk
/// \class ttkRipsMultiParameterSampling
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::RipsMultiParameterSampling module.
///
/// This VTK filter uses the ttk::RipsMultiParameterSampling module to compute an averaging of
/// the data values of an input point data array defined on the input
/// vtkDataSet.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/RipsMultiParameterSampling/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::RipsMultiParameterSampling
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include "PCAnalysis.h"
#include <fstream>
#include <sstream>
#include <ttkRipsMultiParameterSamplingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkMacros.h>
#include <vtkSetGet.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkMultiBlockDataSet.h>

/* Note on including VTK modules
 *
 * Each VTK module that you include a header from needs to be specified in this
 * module's vtk.module file, either in the DEPENDS or PRIVATE_DEPENDS (if the
 * header is included in the cpp file only) sections.
 *
 * In order to find the corresponding module, check its location within the VTK
 * source code. The VTK module name is composed of the path to the header. You
 * can also find the module name within the vtk.module file located in the same
 * directory as the header file.
 *
 * For example, vtkSphereSource.h is located in directory VTK/Filters/Sources/,
 * so its corresponding VTK module is called VTK::FiltersSources. In this case,
 * the vtk.module file would need to be extended to
 *
 * NAME
 *   ttkRipsMultiParameterSampling
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <RipsMultiParameterSampling.h>

using namespace ttk;

enum class GRIDMODE{COMPACT=0, FULL=1, SPARSE=2};

enum class DENSITYMODE{EUCLIDEAN=0, EUCLIDEAN_CONSTRAINED=1, GRAPH=2, MST=3};
using PrimitiveType = Primitive::PrimitiveType;
using STEP = Timings::STEP;

class TTKRIPSMULTIPARAMETERSAMPLING_EXPORT ttkRipsMultiParameterSampling
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::RipsMultiParameterSampling // and we inherit from the base class
{
private:
  std::string PersistenceQuery{0};
  std::string RatioQuery{0};
  std::string OutputFileName{0};
  std::string fittingStringId_{""};
  std::string epsSamplingStringId_{""};
  std::string savefile_path_{"./.cache/savefiles/"};
  std::string TimingFilePrefix{"timings"};

  float FitScale{0.1};
  float Spacing{1};

  float ScaleThresholdMin{0.0};
  float ScaleThresholdMax{2};
  float ScaleRangeMin{-1};

  double DensitySigma{0.01};
  float MinDensity{0};
  double MaxEdgeLength{-1};
  double Epsilon{1e-5};
  int GrowPrimitives{0};
  float GrowPrimitiveThreshold{0.01};
  int GrowPrimitiveIterations{1};
  double max_epsilon_{1e-5};
  SimplexId EpsilonNbMaxPoint{1000000};

  double NoiseLevel{false};

  // match one component options
  SimplexId InputCC{0};
  SimplexId InputTreeId{0};
  SimplexId InputThreshId{0};
  float ScoreMin{0.7};
  double BirthMax{1}; // all births should be negative
  double PersistenceMin{0.0};
  double ProductPersistenceMin{0.0};
  int NbPrimitives{1};
  int NbCCMaxPV{-1};

  double TauWeight{1};
  double KappaWeight{1};

  double CoverCriteria{0.8};

  bool OutputSimplicialComplex{false};

  // boolean options
  bool WriteTimings{false};
  bool UpdateSampling{false};
  bool UpdateFitting{true};
  bool UpdateEpsilonSampling{true};
  bool UpdatePers{true};
  bool UpdateDensity{true};
  bool UpdateDisplay{true};

  bool AddRandomNoise{false};

  bool UseSparseGrid{false};
  bool UseFullGrid{false};
  bool UseSampling{true};
  bool UseProjectiveNorm{false};
  bool UseProductPersistenceScore{false};
  GRIDMODE GridMode{GRIDMODE::COMPACT};
  DENSITYMODE DensityMode{DENSITYMODE::MST};
  bool UseMLSCoordinates{false};
  bool DisplayFit{false};
  bool UseKnnGraph{false};
  bool UseRIMLS{false};
  bool UseScaleDistanceDensity{false};
  bool UseEdgeDensity{false};
  bool UseMixAllScales{false};
  bool UseConstrainedDensity{true};
  bool UseGraphDensity{false};
  bool LoadExistingFits{true};
  bool LoadExistingDensity{true};
  bool WriteMatrix{false};

  bool Compute3DEmbedding{false};

  bool RatioVersion{false}; //use the prototype with hybrid metrics and ratios
  bool UsePlane{false};
  bool UseNormalizedCoords{true};
  bool EvaluationAtPoint{true};
  bool ChangeBasis{true};
  bool LocalMetrics{false};
  bool ShowOnlyLeaves{false};
  bool UseSmallScaleMatrixOrdering{false};
  int MatrixOrdering{0};
  SimplexId BaryMaxId{-1};

  std::string ExternalLabelFile{"default"};
  std::string ExportMatrixFileName{"./distanceMatrix.csv"};
;
  SimplexId NbPointsMax{10000};

  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> diagramPerRatio_{};
  vtkSmartPointer<vtkMultiBlockDataSet> outputMB_{};
  std::vector<std::vector<vtkSmartPointer<vtkUnstructuredGrid>>> outputGrid_{};

  // =================================================
  //                    Params
  // =================================================
  // - subsampling
  ScaleSpaceSampling<pointCloud::MyPoint::Scalar, pointCloud::MyPoint> sampler_subsampling_{}; // 3D sampler for the initial subsampling
  std::vector<std::vector<SimplexId>> sampling_indexes_{};
  int subsampling_scale_index_{-1};
  // - fitting
  int FitDimension{1};
  PrimitiveType FitType{PrimitiveType::SPHERE};
  std::vector<std::vector<double>> fitCoordsVec_{};
  std::vector<double> densityVec_{};
  std::vector<std::vector<int>> isStableVec_{};
  std::vector<std::vector<float>> mlsCoordsVec_{};
  std::vector<std::vector<float>> mlsNormalsVec_{};
  std::vector<float> resultBarycenter_{};
  std::vector<std::vector<double>> embeddingCoordinates_{};
  std::vector<std::vector<SimplexId>> filteredToLocalId_{};
  ScaleSpaceSampling<pointCloud::MyPoint::Scalar, pointCloud::MyPoint> sampler_fitting_{}; // 3D sampler for the fitting
  // - epsilon sampling
  std::vector<SimplexId> pointsAtScaleStartIdx_{};
  std::vector<std::vector<SimplexId>> epsilonIdToFilteredId_{};
  std::vector<std::vector<SimplexId>> filteredIdxToEpsSampling_{};
  std::vector<std::vector<std::vector<double>>> distanceMatrix_{};
  std::vector<EpsilonSampling<pointCloud::MyPointND<1>::Scalar, pointCloud::MyPointND<1>>> epsilonSamplers1_{};
  std::vector<EpsilonSampling<pointCloud::MyPointND<2>::Scalar, pointCloud::MyPointND<2>>> epsilonSamplers2_{};
  std::vector<EpsilonSampling<pointCloud::MyPointND<3>::Scalar, pointCloud::MyPointND<3>>> epsilonSamplers3_{};
  std::vector<EpsilonSampling<pointCloud::MyPointND<4>::Scalar, pointCloud::MyPointND<4>>> epsilonSamplers4_{};
  std::vector<EpsilonSampling<pointCloud::MyPointND<5>::Scalar, pointCloud::MyPointND<5>>> epsilonSamplers5_{};
  std::vector<EpsilonSampling<pointCloud::MyPointND<10>::Scalar, pointCloud::MyPointND<10>>> epsilonSamplers10_{};
  std::vector<SimplexId> nbPointsSampling_{};
	
 // - filtration
	std::vector<ttk::edgeTuple> mergingEdges_{};
	std::vector<ttk::edgeTuple> ripsMSEdges_{};
	std::vector<SimplexId> sizes_ms_{};
	std::vector<ttk::edgeTuple> mergingEdgesDensity_{};
  // =================================================
  SimplexId nbPointsMS_{-1};

  



  vtkMTimeType MTime_{};
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};

  double findDiameter(float* coordinates, SimplexId nbPoints);

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */

  void checkFitDimension(){
    switch (FitType) {
      case PrimitiveType::PLANE:
        if(FitDimension != 4){
          printWrn("Wrong dimension for plane fitting, default to 4");
          FitDimension = 4;
        }
        break;
      case PrimitiveType::QUADRIC:
        if(FitDimension != 10 and FitDimension != 2){
          printWrn("Wrong dimension for quadric fitting, default to 10");
          FitDimension = 10;
        }
        break;
      default:
        if(FitDimension != 1 and FitDimension!=2 and FitDimension!=5 and FitDimension!=3){
          printWrn("Wrong dimension for sphere fitting, default to 5");
          FitDimension = 5;
        }
    }
  }

    inline void SetExternalLabelFile(const std::string &path) {
      ExternalLabelFile = path;
      Modified();
      UpdateDisplay=true;
    }
    inline void SetExportMatrixFileName(const std::string &path) {
      ExportMatrixFileName = path;
      UpdateFitting=true;
      Modified();
    }

  vtkGetMacro(UseSmallScaleMatrixOrdering, int);
  void SetUseSmallScaleMatrixOrdering(int data){
    Modified();
    UpdateDisplay=true;
    UseSmallScaleMatrixOrdering = data;
  }
  vtkGetMacro(WriteMatrix, int);
  void SetWriteMatrix(int data){
    Modified();
    UpdateFitting=true;
    WriteMatrix = data;
  }
  vtkGetMacro(NoiseLevel, double);
  void SetNoiseLevel(double data){
    Modified();
    UpdateFitting=true;
    NoiseLevel = data;
  }
  vtkGetMacro(FitDimension, int);
  void SetFitDimension(int data){
    Modified();
    clearEpsilonSamplers();
    UpdateFitting=true;
    FitDimension = data;
  }
  vtkGetMacro(OutputSimplicialComplex, int);
  void SetOutputSimplicialComplex(int data){
    Modified();
    UpdateDisplay=true;
    OutputSimplicialComplex = data;
  }

  vtkGetMacro(AddRandomNoise, int);
  void SetAddRandomNoise(int data){
    Modified();
    UpdateFitting=true;
    AddRandomNoise = data;
  }
  vtkGetMacro(MatrixOrdering, int);
  void SetMatrixOrdering(int data){
    Modified();
    UpdateDisplay=true;
    MatrixOrdering = data;
  }

  vtkGetMacro(RatioVersion, int);
  void SetRatioVersion(int data){
    Modified();
    UpdateFitting=true;
    RatioVersion = data;
  }

  vtkGetMacro(ShowOnlyLeaves, int);
  void SetShowOnlyLeaves(int data){
    Modified();
    UpdateDisplay=true;
    ShowOnlyLeaves = data;
  }
  vtkGetMacro(LocalMetrics, int);
  void SetLocalMetrics(int data){
    Modified();
    UpdateFitting=true;
    LocalMetrics = data;
  }

  vtkGetMacro(BaryMaxId, int);
  void SetBaryMaxId(int data){
    Modified();
    UpdateFitting=true;
    BaryMaxId = data;
  }
  vtkGetMacro(ChangeBasis, int);
  void SetChangeBasis(int data){
    Modified();
    UpdateFitting=true;
    ChangeBasis = data;
  }

  vtkGetMacro(UseMLSCoordinates, int);
  void SetUseMlSCoordinates(int data){
    Modified();
    UpdateDisplay=true;
    UseMLSCoordinates = data;
  }

  vtkGetMacro(UseGraphDensity, int);
  void SetUseGraphDensity(int data){
    Modified();
    UpdateDensity=true;
    UseGraphDensity = data;
  }
  vtkGetMacro(UseMixAllScales, int);
  void SetUseMixAllScales(int data){
    Modified();
    UpdatePers=true;
    UseMixAllScales = data;
  }
  vtkGetMacro(UseEdgeDensity, int);
  void SetUseEdgeDensity(int data){
    Modified();
    UpdateDensity=true;
    UseEdgeDensity = data;
  }
  vtkGetMacro(UseScaleDistanceDensity, int);
  void SetUseScaleDistanceDensity(int data){
    Modified();
    UpdateDensity=true;
    UseScaleDistanceDensity = data;
  }
  vtkGetMacro(UseConstrainedDensity, int);
  void SetUseConstrainedDensity(int data){
    Modified();
    UpdateDensity=true;
    UseConstrainedDensity = data;
  }

  vtkGetMacro(UseRIMLS, int);
  void SetUseRIMLS(int data){
    Modified();
    UpdateFitting=true;
    UseRIMLS = data;
  }

  vtkGetMacro(LoadExistingDensity, int);
  void SetLoadExistingDensity(int data){
    Modified();
    LoadExistingDensity = data;
    UpdateDensity=true;
  }

  vtkGetMacro(LoadExistingFits, int);
  void SetLoadExistingFits(int data){
    Modified();
    LoadExistingFits = data;
    UpdateFitting = true;
  }

  vtkGetMacro(UseKnnGraph, int);
  void SetUseKnnGraph(int data){
    Modified();
    UpdateFitting=true;
    UseKnnGraph = data;
  }

  vtkGetMacro(DisplayFit, int);
  void SetDisplayFit(int data){
    Modified();
    UpdateDisplay=true;
    DisplayFit = data;
  }

  vtkGetMacro(UseNormalizedCoords, int);
  void SetUseNormalizedCoords(int data){
    Modified();
    UpdateFitting=true;
    UseNormalizedCoords = data;
  }

  vtkGetMacro(UsePlane, int);
  void SetUsePlane(int data){
    Modified();
    UpdateFitting=true;
    UsePlane = data;
  }

  vtkGetMacro(EvaluationAtPoint, int);
  void SetEvaluationAtPoint(int data){
    Modified();
    UpdateFitting=true;
    EvaluationAtPoint = data;
  }

  vtkGetEnumMacro(FitType, PrimitiveType);
  void SetFitType(int type){
    if(this->FitType != static_cast<PrimitiveType>(type)) {
      this->FitType = static_cast<PrimitiveType>(type);
      this->Modified();
      UpdateFitting = true;
    }
  }
  void SetDensityMode(int mode){
    if(this->DensityMode != static_cast<DENSITYMODE>(mode)) {
      this->DensityMode = static_cast<DENSITYMODE>(mode);
      this->Modified();
      UpdateDensity = true;
      if(DensityMode==DENSITYMODE::GRAPH && ripsMSEdges_.size()==0){
        UpdatePers = true;
      }
    }
  }
  void SetGridMode(int mode){
    if(this->GridMode != static_cast<GRIDMODE>(mode)) {
      this->GridMode = static_cast<GRIDMODE>(mode);
      this->Modified();
      UpdatePers = true;
    }
  }

  void GraphOutput(){
    outputGraph( OutputFileName, persistenceGrid_);
  }

  void SetCoverCriteria(double data){
    CoverCriteria = data;
    Modified();
    UpdateDisplay=true;
  }

  vtkGetMacro(CoverCriteria, double);

  void SetWriteTimings(int data){
    Modified();
    WriteTimings=data;
    UpdateFitting=true;
    UpdateSampling=true;
  }
  vtkGetMacro(WriteTimings, int);

  void SetUseProductPersistenceScore(int data){
    Modified();
    UseProductPersistenceScore=data;
    UpdateDisplay=true;
  }
  vtkGetMacro(UseProductPersistenceScore, int);

  void SetUseProjectiveNorm(int data){
    Modified();
    UseProjectiveNorm=data;
    UpdateDisplay=true;
    UpdateFitting=true;
    UpdatePers=true;
  }
  vtkGetMacro(UseProjectiveNorm, int);

  void SetUseSampling(int data){
    Modified();
    UseSampling=data;
    UpdateDisplay=true;
    UpdateFitting=true;
    UpdatePers=true;
  }
  vtkGetMacro(UseSampling, int);

  void SetCompute3DEmbedding(int data){
    Modified();
    Compute3DEmbedding=data;
    UpdateDisplay=true;
    UpdateFitting=true;
    UpdatePers=true;
  }
  vtkGetMacro(Compute3DEmbedding, int);

  void SetEpsilonNbMaxPoint(int data){
    Modified();
    EpsilonNbMaxPoint=data;
    UpdateEpsilonSampling=true;
  }
  vtkGetMacro(EpsilonNbMaxPoint, int);

  void SetBirthMax(double data){
    Modified();
    BirthMax=data;
    UpdateDisplay=true;
  }
  vtkGetMacro(BirthMax, double);

  void SetProductPersistenceMin(double data){
    Modified();
    ProductPersistenceMin=data;
    UpdateDisplay=true;
  }
  vtkGetMacro(ProductPersistenceMin, double);
  void SetPersistenceMin(double data){
    Modified();
    PersistenceMin=data;
    UpdateDisplay=true;
  }
  vtkGetMacro(PersistenceMin, double);

  void SetScoreMin(float data){
    Modified();
    ScoreMin=data;
    UpdateDisplay=true;
  }
  vtkGetMacro(ScoreMin, float);

  void SetNbPrimitives(int data){
    Modified();
    NbPrimitives=data;
    UpdateDisplay=true;
  }
  vtkGetMacro(NbPrimitives, int);

  void SetMetricWeights(double tauWeight, double kappaWeight){
    TauWeight = tauWeight;
    KappaWeight = kappaWeight;
    Modified();
    UpdateFitting = true;
  }
  void SetInputCCParams(int inputcc, int row, int col){
    InputCC = inputcc;
    InputTreeId=row;
    InputThreshId=col;
    Modified();
    UpdatePers=true;
  }

  vtkGetMacro(GrowPrimitiveIterations, int);
  void SetGrowPrimitiveIterations(int data){
    GrowPrimitiveIterations = data;
    Modified();
    UpdateDisplay=true;
  }
  vtkGetMacro(GrowPrimitiveThreshold, float);
  void SetGrowPrimitiveThreshold(float data){
    GrowPrimitiveThreshold = data;
    Modified();
    UpdateDisplay=true;
  }
  vtkGetMacro(GrowPrimitives, int);
  void SetGrowPrimitives(int data){
    GrowPrimitives = data;
    if(GrowPrimitives and subsampling_scale_index_<0){
      UpdateSampling = true;
    }
    Modified();
    UpdateDisplay=true;
  }

  vtkGetMacro(Epsilon, double);
  void SetEpsilon(double data){
    Epsilon = data;
    Modified();
    UpdateEpsilonSampling=true;
  }

  vtkGetMacro(DensitySigma, double);
  void SetDensitySigma(double data){
    DensitySigma = data;
    Modified();
    UpdateDensity=true;
  }

  vtkGetMacro(ScaleRangeMin, float);
  void SetScaleRangeMin(float data){
    ScaleRangeMin = data;
    Modified();
    UpdateDisplay=true;
  }
  vtkGetMacro(ScaleThresholdMax, float);
  void SetScaleThresholdMax(float data){
    ScaleThresholdMax = data;
    Modified();
    UpdateDisplay=true;
  }
  vtkGetMacro(ScaleThresholdMin, float);
  void SetScaleThresholdMin(float data){
    ScaleThresholdMin = data;
    Modified();
    UpdateDisplay=true;
  }
  vtkGetMacro(MaxEdgeLength, double);
  void SetMaxEdgeLength(double data){
    MaxEdgeLength = data;
    Modified();
    UpdatePers=true;
  }
  vtkGetMacro(MinDensity, float);
  void SetMinDensity(float data){
    MinDensity = data;
    Modified();
    UpdateFitting=true;
  }

  vtkGetMacro(FitScale, float);
  void SetFitScale(float data){
    FitScale = data;
    Modified();
    UpdateFitting=true;
  }

  vtkSetMacro(InputCC, int);
  vtkGetMacro(InputCC, int);

  vtkSetMacro(InputThreshId, int);
  vtkGetMacro(InputThreshId, int);

  vtkSetMacro(InputTreeId, int);
  vtkGetMacro(InputTreeId, int)

  // vtkSetEnumMacro(GridMode, GRIDMODE);
  vtkGetEnumMacro(GridMode, GRIDMODE);

  void SetNbCCMax(const int data){
    // UpdatePers=true;
    UpdateDisplay=true;
    NbCCMaxPV = data;
    Modified();
  }

  void SetNbPointsMax(const int data){
    UpdateSampling=true;
    UpdateFitting=true;
    NbPointsMax = data;
    Modified();
  }
  vtkGetMacro(NbPointsMax, int);

  void SetDiagMax(const double data){
    // UpdatePers=true;
    UpdateDisplay=true;
    DiagMax = data;
    Modified();
  }
  vtkGetMacro(DiagMax, double);
  void SetDiagMaxDisplay(const double data){
    UpdateDisplay=true;
    DiagMaxDisplay = data;
    Modified();
  }
  vtkGetMacro(DiagMaxDisplay, double);

  void SetSizeMin(const int data){
    // UpdatePers=true;
    UpdateDisplay=true;
    SizeMin = data;
    Modified();
  }
  vtkGetMacro(SizeMin, int);
  vtkGetMacro(NbCCMax, int);

  void SetPersistenceQuery(const std::string & data){
    UpdatePers=true;
    UpdateDisplay=true;
    PersistenceQuery = data;
    Modified();
  }
  vtkGetMacro(PersistenceQuery, std::string);

  void SetRatioQuery(const std::string & data){
    UpdateFitting=true;
    UpdatePers=true;
    UpdateDisplay=true;
    RatioQuery = data;
    Modified();
  }
  vtkGetMacro(RatioQuery, std::string);

  void SetTimingFilePrefix(const std::string& data){
    TimingFilePrefix = data;
  }
  vtkGetMacro(TimingFilePrefix, std::string);

  void SetOutputFileName(const std::string& data){
    OutputFileName = data;
  }
  vtkGetMacro(OutputFileName, std::string);

  void SetSamplingFactor(double data){
    UpdateFitting=true;
    Modified();
    SamplingFactor=data;
  }
  vtkGetMacro(SamplingFactor, double);

  vtkGetMacro(Spacing, float);
  void SetSpacing(float data){
    Spacing = data;
    Modified();
    UpdateDisplay=true;
  }
  void SetScaleSamplingMinMaxCount(double min, double max, double c){
   ScaleSamplingMin=min;
   ScaleSamplingMax=max;
   ScaleSamplingCount=c;
   Modified();
   UpdateFitting=true;
   UpdateEpsilonSampling=true;
  }
  vtkSetMacro(ScaleSamplingBase, double)
  vtkGetMacro(ScaleSamplingBase, double)

  vtkSetMacro(ScaleSamplingMin, double)
  vtkGetMacro(ScaleSamplingMin, double)

  vtkSetMacro(ScaleSamplingMax, double)
  vtkGetMacro(ScaleSamplingMax, double)

  vtkSetMacro(ScaleSamplingCount, double)
  vtkGetMacro(ScaleSamplingCount, double)

  void SetJaccardConstraint(double data){
    JaccardConstraint = data;
    UpdatePers=true;
    Modified();
  }
  vtkGetMacro(JaccardConstraint, double)

  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  void SetUseNF(int data){
    UseNF = data;
    UpdateFitting = true;
    Modified();
  }
  vtkGetMacro(UseNF, int)

  inline std::string fittingStringId()const {
    if(fittingStringId_==""){
      printErr("Run generateFittingStringId() first");
      return "{}";
    }
    return fittingStringId_;
  }

  inline std::string epsSamplingStringId()const {
    if(epsSamplingStringId_==""){
      printErr("Run generateEpsSamplingStringId() first");
      return "{}";
    }
    return epsSamplingStringId_;
  }

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkRipsMultiParameterSampling *New();
  vtkTypeMacro(ttkRipsMultiParameterSampling, ttkAlgorithm);

  void resetSelection(vtkUnstructuredGrid* filteredTree);
  void applySelection(SimplexId nbCC, vtkUnstructuredGrid* filteredTree, std::vector<SimplexId> selection);
  void readSelection(std::vector<SimplexId>& selectionVec, vtkPointSet* inputVTKSelection);
  void translate(vtkUnstructuredGrid* vtu, int I, int J, int K, float spacing, float* coordinates);
  void fillSparseGridPersistences(const std::vector<std::vector<double>>& persBigVec,std::vector<std::vector<double>>& sparseGrid, int nbRows);
  void fillFullGridPersistences(const std::vector<std::vector<double>>& persBigVec,std::vector<std::vector<double>>& grid, int nbRows);

// void computeStatistics(const RipsTree& tree, double threshold, float* meanNormals, int* parentSize, int* parent, int* parent2d);

  void computeStatistics(const RipsTree& tree,
      int iRatio, int iThreshold, const std::vector<std::vector<double>>& persistenceGrid,
      float* normals, float* meanNormals, int* parentSize, int* fullSize, int* parent, int* parent2d, int* fullId, double* pers1d, double* pers2d, int* pers2dId, int* localIdPtr, int* primitiveTypePtr, int* displayPtr);
      void setSingleComponentAttributes(vtkUnstructuredGrid* vtu, SimplexId fullId, double persistence2D, SimplexId persistence2DId,int criticalType, double diag, double scale, double autosimilarity);
void createOutputCommonComponentPoints(vtkUnstructuredGrid* vtu, SimplexId fullId, double CoverCriteria, std::vector<std::vector<double>>& persistenceGrid, std::vector<double>& ratioVec, float* coordinates, float* normals);

int ProductDiagramToVTU(vtkUnstructuredGrid *vtu,
                 const ttk::DiagramType &diagram,
                 const std::vector<SimplexId>& flags,
                 const std::vector<SimplexId>& chosenPair,
                 const std::vector<std::pair<int, int>>& scalePairs,
                 const std::vector<float>& scales,
                 const std::vector<float>& scores);
void createOutputPrimitive(const Primitive& primitive, vtkUnstructuredGrid* vtu, SimplexId primId, const float* coordinates, const float* normals,const std::vector<float>& pointCloudBarycenter);
void createOutputPrimitiveSegmentation(const PrimitiveSegmentation& seg, vtkMultiBlockDataSet* mb, const float* coordinates, const float* normals, const std::vector<float>& pointCloudBarycenter);

void clearEpsilonSamplers(){
  for(auto& e : epsilonSamplers1_){
    e.clear_state();
  }
  epsilonSamplers1_.clear();
  
  for(auto& e : epsilonSamplers2_){
    e.clear_state();
  }
  epsilonSamplers2_.clear();

  for(auto& e : epsilonSamplers4_){
    e.clear_state();
  }
  epsilonSamplers4_.clear();

  for(auto& e : epsilonSamplers5_){
    e.clear_state();
  }
  epsilonSamplers5_.clear();

  for(auto& e : epsilonSamplers10_){
    e.clear_state();
  }
  epsilonSamplers10_.clear();
}

template<class PointND>
bool createOutputMergingEdges(vtkPolyData* vtp,
    std::vector<ttk::edgeTuple>& mergingEdges,
                                              float* coordinates,
                                              float* normals,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              std::vector<std::vector<float>>& mlsCoordinates,
                                              std::vector<std::vector<float>>& mlsNormals,
                                              std::vector<double>& densityVec,
                                              std::vector<std::vector<double>>& embeddingCoordinates,
                                              const std::vector<float>& resultBarycenter,
                                              std::vector<std::vector<SimplexId>>& sampling_indexes,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalIdVec,
                                              std::vector<SimplexId>& pointsAtScaleStartIdx,
                                              std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId);
  void writeCoordinatesForTomato(std::string name,
                                              float* coordinates,
                                              float* normals,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              std::vector<std::vector<float>>& mlsCoordinates,
                                              std::vector<std::vector<float>>& mlsNormals,
                                              std::vector<double>& densityVec,
                                              std::vector<std::vector<double>>& embeddingCoordinates,
                                              std::vector<std::vector<SimplexId>>& sampling_indexes,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalIdVec,
                                              std::vector<SimplexId>& pointsAtScaleStartIdx,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId);
void writeFitCoordinates(float* coords, double* fitCoords, SimplexId nbPoints, int dim, double relative_scale, std::string name);
  void ripsArcsToMergeTree(vtkMultiBlockDataSet* outputMTMB, std::vector<RipsArc>& ripsArcs, std::vector<SimplexId>& arcIds, std::vector<SimplexId>& leafIds);

  std::string generateFittingStringId(float* coordinates, float* normals, SimplexId nbPoints);

void listAvailableFits();
  bool loadFits();
  void saveFits();

  bool loadDensity();
  void saveDensity();


protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkRipsMultiParameterSampling();
  ~ttkRipsMultiParameterSampling() override = default;

  /**
   * TODO 8: Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 9: Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};

int readIntFromString(std::string& input, int n_digit){
  int v = std::atoi(input.substr(0,n_digit).c_str());
  input.erase(0,n_digit);
  return v;
}

int readBoolFromString(std::string& input){
  return readIntFromString(input, 1);
}


double readDoubleFromString(std::string& input, int n_digit, int prec){
  double v = readIntFromString(input, n_digit);
  v /= std::pow(10, prec);
  return v;
}

inline std::stringstream intToStringStream(int n, int max_digit=-1){
  std::stringstream ss;

  std::string sign = "";
  if(n<0){
    sign = "-";
  }

  if(max_digit==-1){
    ss << sign << std::abs(n);
  }else{
      ss << sign << std::setfill('0') << std::setw(max_digit) << std::abs(n);
  }
  return ss;
}

inline std::stringstream doubleToStringStream(double x, int prec = 2){
  int n = x*std::pow(10, prec);
  return intToStringStream(n);
}

inline std::stringstream doubleToStringStream(double x, int max_digit, int prec){
  int n = x*std::pow(10, prec);
  return intToStringStream(n, max_digit);
}

inline std::stringstream doubleValuesToStringStream(std::vector<double>& input, int max_digit=4, int prec = 3){
  std::stringstream ss;
  for(auto v : input){
    ss << doubleToStringStream(v, max_digit, prec).str();
  }
  return ss;
}
inline std::stringstream valuesToStringStream(std::vector<double>& input){
  std::stringstream ss;
  for(auto v : input){
    ss << doubleToStringStream(v).str();
  }
  return ss;
}

inline std::stringstream intValuesToStringStream(std::vector<int>& input, int max_digit){
  std::stringstream ss;
  for(auto v : input){
    ss << intToStringStream(v, max_digit).str();
  }
  return ss;
}

inline std::stringstream boolValuesToStringStream(std::vector<int>& input){
  std::stringstream ss;
  for(auto v : input){
    ss << (v ? "1" : "0");
  }
  return ss;
}

template<typename T>
void vectorToOfStream(const std::vector<T>& input, std::ofstream& out){
  size_t size = input.size();
  out.write((char*)&size, sizeof(size_t));
  out.write(reinterpret_cast<const char*>(&input[0]), input.size()*sizeof(T));
}

template<typename T>
void readVectorFromStream(std::vector<T>& output, std::istream& in){
  size_t size_read;
  in.read((char*)&size_read, sizeof(size_t));
  output.resize(size_read);
  in.read((char*)&output[0], sizeof(T)*size_read);
}

template<typename T>
void writeMultiscaleVector(const std::vector<std::vector<T>>& input, std::ofstream& out){
  size_t size = input.size();
  out.write((char*)&size, sizeof(size_t));
  for(auto& v : input){
    vectorToOfStream(v, out);
  }
}

template<typename T>
void readMultiscaleVector(std::vector<std::vector<T>>& output, std::istream& in){
  size_t size_read;
  in.read((char*)&size_read, sizeof(size_t));
  output.resize(size_read);
  for(auto& v : output){
    readVectorFromStream(v, in);
  }
}
