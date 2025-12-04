#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <limits>
#include <numeric>
#include <ostream>
#include <ratio>
#include <sstream>
#include <filesystem>
#include <string>
#include <thread>
#include <ttkRipsMultiParameterSampling.h>

#include <utility>
#include <vector>
#include <vtkCellType.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkIndent.h>
#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataSet.h>
#include <vtkInformationVector.h>
#include <vtkAlgorithmOutput.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkSystemIncludes.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkPCellDataToPointData.h>
#include <vtkMaskPoints.h>
#include <vtkPoints.h>
#include <vtkThreshold.h>
#include <vtkCellData.h>
#include <vtkIndent.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <PCAnalysisUtils.h>
#include <PCAnalysis.h>
#include <Sampling.h>
#include "../../base/ripsPersistence/nanoflann.hpp"
#include "DataTypes.h"
#include "Debug.h"
#include "DeprecatedDataTypes.h"
#include "FTMDataTypes.h"
#include "OpenMP.h"
#include "RIMLS.h"
#include "RipsMultiParameterSampling.h"
#include <RipsPersistence.h>
#include <FTMTree.h>
#include <FTMTreeUtils.h>
#include <FTMTreePPUtils.h>
#include <MergeTree.h>
#include <RipsTree.h>
#include <ttkRipsPersistence.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkRipsTree.h>
#include <ttkMergeTreeVisualization.h>
#include <Eigen/Eigen>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/Map.h>

namespace fs = std::filesystem;

using namespace ttk;
// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkRipsMultiParameterSampling);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkRipsMultiParameterSampling::ttkRipsMultiParameterSampling() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(8);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkRipsMultiParameterSampling::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
    return 1;
  }else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkRipsMultiParameterSampling::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }else if(port==1){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }else if(port==2){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }else if(port==3){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }else if(port==4){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }else if(port==5){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }else if(port==6){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }else if(port==7){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkRipsMultiParameterSampling::RequestInformation(
    vtkInformation *,
    vtkInformationVector **ttkNotUsed(inputVector),
    vtkInformationVector *outputVector) {


  int extent[6];
  extent[0] = 0;
  extent[1] = 99;
  extent[2] = 0;
  extent[3] = 99;
  extent[4] = 0;
  extent[5] = 99;

  auto outInfo = outputVector->GetInformationObject(1);
  outInfo->Set(
      vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

  return 1;
}

int ttkRipsMultiParameterSampling::RequestData(vtkInformation *request,
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {

  if(WriteTimings){
    LoadExistingFits=false;
  }

  Timings timings{};
  timings.init();
  timings.clock(STEP::TOTAL);
  timings.clock(STEP::INIT);

  // vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::GetData(outputVector, 6);
  // int nblocks = mb->GetNumberOfBlocks();
  // for(int iblock=0; iblock<nblocks; iblock++){
  //   mb->GetBlock(iblock)->Delete();
  // }
// vtkMultiBlockDataSet::GetData(outputVector, 6)->SetNumberOfBlocks(0);

  Timer tm_total{};

  vtkPointSet *inputDataSet = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet *inputSelection = vtkPointSet::GetData(inputVector[1]);
  if(!inputDataSet)
    return 0;

  vtkDataArray* inputCoordinates = inputDataSet->GetPoints()->GetData();
  int nbPoints = inputDataSet->GetNumberOfPoints();

  vtkDataArray *inputNormals = this->GetInputArrayToProcess(0, inputVector);
  if(!inputNormals){
    printWrn("No Normal field was found -- Normals won't be used");
  }else{
    printMsg("Using normal field: "+std::string(inputNormals->GetName()));
  }

  float* coordinates = (float*) ttkUtils::GetVoidPointer(inputCoordinates);
  float* normals = (float*) ttkUtils::GetVoidPointer(inputNormals);


  // ==================================================
  // Lazy computation setup: we detect a change in the upstream pipeline and put UpdateFitting to true if so
  auto nbOfInputPorts = this->GetInputAlgorithm()->GetNumberOfInputPorts();
  // don't proceed if the input algorithm has zero input ports (standalone exec)
  if(nbOfInputPorts){
    // we have to check for the input of the input as paraview seems to have some weird algorithms between filters
    // (such as "vtkPVPostFilter")
    auto MTime = this->GetInputAlgorithm()->GetInputAlgorithm()->GetMTime();
    auto name = this->GetInputAlgorithm()->GetInputAlgorithm();
    std::cout<<"New time "<<MTime<<" vs old time "<<MTime_<<std::endl;
    if(MTime != MTime_){ // modified time -> input data has changed
      std::cout<<" => toggle update"<<std::endl;
      UpdateSampling = true; // recompute everything
      UpdateFitting = true;
      UpdateEpsilonSampling = true;
      UpdatePers = true;
      UpdateDisplay = true;
    }
    MTime_ = MTime;
  }

  checkFitDimension();

  // ==================================================
  generateFittingStringId(coordinates, normals, nbPoints);
  // ==================================================

  if(SizeMin==-1){
    SizeMin = UseSampling ? NbPointsMax/200 : nbPoints/200;
  }




  // double bounds[6];
  // inputDataSet->GetBounds(bounds);

  // std::string outputPointsFileName = "output_points.csv";
  // std::ofstream ofs;
  // ofs.open(outputPointsFileName, std::ofstream::out | std::ofstream::trunc);
  // ofs.close();



  double diag_aabb = findDiameter(coordinates, nbPoints);//(bounds[1]-bounds[0])*(bounds[1]-bounds[0]) + (bounds[3]-bounds[2])*(bounds[3]-bounds[2]) + (bounds[5]-bounds[4])*(bounds[5]-bounds[4]);

  // std::cout<<"aabb diagonal: "<<std::sqrt(diag_aabb)<<std::endl;

  // point selection
  std::vector<SimplexId> selection{};
  std::vector<std::vector<double>>& persistenceGrid = persistenceGrid_;
  readSelection(selection, inputSelection);



  // float* fitCoordinates; 

  std::vector<double>&  densityVec = densityVec_;
  // float* density;
  std::vector<std::vector<SimplexId>>& filteredToLocalIdVec = filteredToLocalId_;
  std::vector<EpsilonSampling<pointCloud::MyPointND<1>::Scalar, pointCloud::MyPointND<1>>>& epsilonSamplers1 = epsilonSamplers1_;
  std::vector<EpsilonSampling<pointCloud::MyPointND<2>::Scalar, pointCloud::MyPointND<2>>>& epsilonSamplers2 = epsilonSamplers2_;
  std::vector<EpsilonSampling<pointCloud::MyPointND<3>::Scalar, pointCloud::MyPointND<3>>>& epsilonSamplers3 = epsilonSamplers3_;
  std::vector<EpsilonSampling<pointCloud::MyPointND<4>::Scalar, pointCloud::MyPointND<4>>>& epsilonSamplers4 = epsilonSamplers4_;
  std::vector<EpsilonSampling<pointCloud::MyPointND<5>::Scalar, pointCloud::MyPointND<5>>>& epsilonSamplers5 = epsilonSamplers5_;
  std::vector<EpsilonSampling<pointCloud::MyPointND<10>::Scalar, pointCloud::MyPointND<10>>>& epsilonSamplers10 = epsilonSamplers10_;

  std::vector<std::vector<SimplexId>>& filteredIdxToEpsSampling = filteredIdxToEpsSampling_;
  std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId = epsilonIdToFilteredId_;

  std::vector<SimplexId>& pointsAtScaleStartIdx = pointsAtScaleStartIdx_;


  std::vector<std::vector<int>>& isStableVec = isStableVec_;
  // int* isStable{}; 

  // std::vector<std::vector<std::vector<double>>>& distanceMatrix = distanceMatrix_;
  std::vector<SimplexId>& nbPointsSampling = nbPointsSampling_;


  std::vector<std::vector<float>>& mlsCoordsVec = mlsCoordsVec_;
  std::vector<std::vector<float>>& mlsNormalsVec = mlsNormalsVec_;
  // float* mlsCoordinates; 
  std::vector<std::vector<double>>& embeddingCoordinates = embeddingCoordinates_;
  std::vector<std::vector<SimplexId>>& sampling_indexes = sampling_indexes_;
  ScaleSpaceSampling<pointCloud::MyPoint::Scalar, pointCloud::MyPoint>& sampler_fitting = sampler_fitting_;




  vtkMultiBlockDataSet *outputScaleMB = vtkMultiBlockDataSet::GetData(outputVector, 0);
  // vtkImageData *outputVTI = vtkImageData::GetData(outputVector, 1);
  // int dimX=100, dimY=100, dimZ=100;
  // outputVTI->SetDimensions(dimX,dimY,dimZ);
  // double offsetBounds[6];
  // double offsetSpacing=0.05;
  // offsetBounds[0] = bounds[0]- offsetSpacing*sqrt(diag_aabb);
  // offsetBounds[1] = bounds[1]+ offsetSpacing*sqrt(diag_aabb);
  // offsetBounds[2] = bounds[2]- offsetSpacing*sqrt(diag_aabb);
  // offsetBounds[3] = bounds[3]+ offsetSpacing*sqrt(diag_aabb);
  // offsetBounds[4] = bounds[4]- offsetSpacing*sqrt(diag_aabb);
  // offsetBounds[5] = bounds[5]+ offsetSpacing*sqrt(diag_aabb);
  // // outputVTI->SetOrigin(offsetBounds[0], offsetBounds[2], offsetBounds[4]);
  // double ox=offsetBounds[0],oy=offsetBounds[2], oz=offsetBounds[4];
  // double dx = (offsetBounds[1]-offsetBounds[0])/(dimX-1);
  // double dy = (offsetBounds[3]-offsetBounds[2])/(dimY-1);
  // double dz = dimZ > 1 ? (offsetBounds[5]-offsetBounds[4])/(dimZ-1) : 0;
  // outputVTI->SetSpacing(dx,dy,dz);

  // vtkNew<vtkDoubleArray> fitScalarField{};
  // fitScalarField->SetNumberOfComponents(1);
  // fitScalarField->SetNumberOfTuples(dimX*dimY*dimZ);
  // fitScalarField->SetName("potential");
  // fitScalarField->Fill(-1);
  // double* fitSFPtr = (double*)ttkUtils::GetVoidPointer(fitScalarField);
  // outputVTI->GetPointData()->AddArray(fitScalarField);

  std::vector<float>& resultBarycenter = resultBarycenter_;


  RipsTree& treeMS = treeMS_;
  DiagramType& diagramMS = diagramMS_;
  DiagramType& ripsDiagramMS = ripsDiagramMS_;
  RipsTree treeDensity{};
	std::vector<ttk::edgeTuple>& mergingEdges = mergingEdges_;
	std::vector<ttk::edgeTuple>& ripsMSEdges = ripsMSEdges_;
  std::cout<<" number of initial merging edges "<<mergingEdges.size()<<std::endl;
  std::vector<SimplexId>& sizes_ms = sizes_ms_;
	std::vector<ttk::edgeTuple>& mergingEdgesDensity = mergingEdgesDensity_;




  // std::vector<float> fitCoordsVec(5*nbPoints);
  // float* fitCoordinates = fitCoordsVec.data(); 
  // PCAnalysis pca{};
  // pca.computeFitCoordinates<pointCloud::myOrientedSphereFitGLS>(coordinates, normals, nbPoints, fitCoordinates, 0.0);

  fitCoordsVec_.resize(ScaleSamplingCount);
  isStableVec.resize(ScaleSamplingCount);
  // isStable = isStableVec.data();
  mlsCoordsVec.resize(ScaleSamplingCount);
  mlsNormalsVec.resize(ScaleSamplingCount);
  sampling_indexes.resize(ScaleSamplingCount);

  switch(FitDimension){
    case 1:
      epsilonSamplers1.resize(ScaleSamplingCount);
      break;
    case 2:
      epsilonSamplers2.resize(ScaleSamplingCount);
      break;
    case 3:
      epsilonSamplers3.resize(ScaleSamplingCount);
      break;
    case 4:
      epsilonSamplers4.resize(ScaleSamplingCount);
      break;
    case 5:
      epsilonSamplers5.resize(ScaleSamplingCount);
      break;
    case 10:
      epsilonSamplers10.resize(ScaleSamplingCount);
      break;
    default:
      printErr("Invalid FitDimension");
      break;
  }
  filteredIdxToEpsSampling.resize(ScaleSamplingCount);
  epsilonIdToFilteredId.resize(ScaleSamplingCount);
  filteredToLocalIdVec.resize(ScaleSamplingCount);
  pointsAtScaleStartIdx.resize(ScaleSamplingCount,0);
  // distanceMatrix.resize(ScaleSamplingCount);
  nbPointsSampling.resize(ScaleSamplingCount);


  timings.clock(STEP::INIT);

  SimplexId nbScales = ScaleSamplingCount;

  timings.clock(STEP::SUBSAMPLING);

  ScaleSpaceSampling<pointCloud::MyPoint::Scalar, pointCloud::MyPoint>& sampler_subsampling = sampler_subsampling_;
  int& subsampling_scale_index = subsampling_scale_index_;
  SimplexId& nbPointsMS = nbPointsMS_;

  if(UpdateFitting){ // try to load existing fit file
    if(LoadExistingFits and loadFits()){
      UpdateFitting = false;
      UpdateEpsilonSampling=true;
    }else{
      listAvailableFits();
      if(subsampling_scale_index<0){ // was never sampled before
        UpdateSampling=true; //needed for now because the kd-tree of the subsampling is used for the fitting
      }
    }
  }

  if(UpdateSampling){


    // TODO: Use farthest point sampling
    sampler_subsampling.setDebugLevel(debugLevel_);
    sampler_subsampling.setSamplingFactor(SamplingFactor);
    sampler_subsampling.sampleScaleLog(0.01*sqrt(diag_aabb), 1*sqrt(diag_aabb), 50);

    sampler_subsampling.poissonDiskSamplingND(nbPoints, coordinates);

    // find appropriate scale
    subsampling_scale_index = sampler_subsampling.findMinScale(NbPointsMax);
    std::cout<<"scale "<<subsampling_scale_index<<" with "<<sampler_subsampling.kdTree(subsampling_scale_index).sample_count()<<std::endl;

    const std::vector<SimplexId>& indexes = sampler_subsampling.kdTree(subsampling_scale_index).samples();
    sampling_indexes[0].clear();
    std::copy(indexes.begin(), indexes.end(), std::back_inserter(sampling_indexes[0]));
    TTK_PSORT(threadNumber_, sampling_indexes[0].begin(), sampling_indexes[0].end());

    for(int iScale=1; iScale<nbScales; iScale++){
      sampling_indexes[iScale].clear();
      sampling_indexes[iScale]=sampling_indexes[0];
    }
    // std::cout<<" subsampling_scale_index "<<subsampling_scale_index<<std::endl;

    UpdateSampling = false;

    timings.clock(STEP::SUBSAMPLING);
  }
  
  sampler_fitting.setDebugLevel(debugLevel_);
  sampler_fitting.setSamplingFactor(SamplingFactor);
  sampler_fitting.sampleScaleLog(ScaleSamplingMin*sqrt(diag_aabb), ScaleSamplingMax*sqrt(diag_aabb), ScaleSamplingCount);
if(UpdateFitting){

  timings.clock(STEP::FITTING);

  // Sampling in 3D to explore different scales for the fitting
  printMsg("Sampling scales");

  timings.clock(STEP::POISSON);

  sampler_fitting.poissonDiskSampling(nbPoints, coordinates, normals);
  sampler_fitting.setUseKnnGraph(UseKnnGraph);
  if(UseKnnGraph){
    sampler_fitting.buildKnnGraph();
  }
  timings.clock(STEP::POISSON);

  // trees_.resize(nbScales);
  // persVec.resize(nbScales);

  // outputMB_ = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  // outputMB_->SetNumberOfBlocks(ratioVec.size());



  std::vector<float> mlsAllCoordinates{};
  mlsAllCoordinates.resize(3*nbPoints);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i=0; i<mlsAllCoordinates.size(); i++){
    mlsAllCoordinates[i]=coordinates[i];
  }


  for(SimplexId iScale = 0; iScale < nbScales; iScale++){


    double* fitCoordinates = fitCoordsVec_[iScale].data();
    // float* density = densityVec[iScale].data();
    int* isStable = isStableVec[iScale].data();
    float* mlsCoordinates = mlsCoordsVec[iScale].data();
    float* mlsNormals = mlsNormalsVec[iScale].data();
    std::vector<SimplexId>& filteredToLocalId = filteredToLocalIdVec[iScale];

    Timer tm_scale{};

    double scale = sampler_fitting.scale(iScale);
    printMsg("Scale "+std::to_string(iScale)+": t = " + std::to_string(scale));

    PCAnalysis pca{};
    pca.setDebugLevel(3);
    pca.setThreadNumber(threadNumber_);

    RipsPersistence ripsPersistence{};
    ripsPersistence.setDebugLevel(debugLevel_);
    ripsPersistence.setThreadNumber(threadNumber_);
    ripsPersistence.setUseProjectiveNorm(UseProjectiveNorm);




    SimplexId nbSampledPoints = sampling_indexes[iScale].size();
    fitCoordsVec_[iScale].resize(FitDimension*nbSampledPoints);
    // densityVec[iScale].resize(nbSampledPoints);
    mlsCoordsVec[iScale].resize(3*nbSampledPoints);
    mlsNormalsVec[iScale].resize(3*nbSampledPoints);
    isStableVec[iScale].resize(nbSampledPoints);
    std::cout<<"three"<<std::endl;

    fitCoordinates = fitCoordsVec_[iScale].data();
    // density = densityVec[iScale].data();
    isStable = isStableVec[iScale].data();
    mlsCoordinates = mlsCoordsVec[iScale].data();
    mlsNormals = mlsNormalsVec[iScale].data();
    std::cout<<"four"<<std::endl;


    // std::vector<std::vector<double>>& distMat = distanceMatrix[iScale];
    switch (FitType) {
      case PrimitiveType::PLANE:
        pca.computeFitCoordinatesPlane<pointCloud::MeanPlaneFit, pointCloud::MeanPlaneFitRIMLS>(coordinates,
            normals,
            nbPoints,
            FitDimension,
            fitCoordinates,
            isStable,
            mlsCoordinates,
            mlsNormals,
            sampler_fitting,
            iScale,
            resultBarycenter,
            sampling_indexes[iScale],
            mlsAllCoordinates.data(),
            BaryMaxId,
            UseNormalizedCoords,
            UseRIMLS);
        break;
      case PrimitiveType::QUADRIC:
        pca.computeFitCoordinatesEllipsoid<pointCloud::FitEllipsoid, pointCloud::FitEllipsoidRIMLS>(coordinates,
            normals,
            nbPoints,
            FitDimension,
            fitCoordinates,
            isStable,
            mlsCoordinates,
            mlsNormals,
            sampler_fitting,
            iScale,
            resultBarycenter,
            sampling_indexes[iScale],
            mlsAllCoordinates.data(),
            BaryMaxId,
            UseNormalizedCoords,
            UseRIMLS);
        break;
      default: // Sphere fit
        pca.computeFitCoordinatesSphere<pointCloud::FitASODiff, pointCloud::FitASODiffRIMLS>(coordinates,
            normals,
            nbPoints,
            FitDimension,
            fitCoordinates,
            isStable,
            mlsCoordinates,
            mlsNormals,
            sampler_fitting,
            iScale,
            resultBarycenter,
            sampling_indexes[iScale],
            mlsAllCoordinates.data(),
            BaryMaxId,
            UseNormalizedCoords,
            UseRIMLS);
    }

    std::vector<SimplexId> densitySampling{};
    sampler_subsampling.filterUnstablePoints(subsampling_scale_index, coordinates, nbPoints, densitySampling, isStable);

    std::cout<<" nb points "<<sampling_indexes[iScale].size()<<" -> "<<densitySampling.size()<<std::endl;
    nbPointsSampling[iScale] = densitySampling.size();
    filteredToLocalId.clear();
    std::copy(densitySampling.begin(), densitySampling.end(), std::back_inserter(filteredToLocalId));





  } // end for scale (fitting)
  
  mlsAllCoordinates.clear();


  UpdateFitting=false;
  UpdateEpsilonSampling = true;

  timings.clock(STEP::FITTING);

  saveFits();

} // end fitting

  std::vector<std::vector<double>> noisyFitVec{};
  if(AddRandomNoise){


    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution d{0.0, 0.5*NoiseLevel};

    printMsg("Add random noise with level " + std::to_string(NoiseLevel));
    for(int itest=0; itest<10; itest++){
    std::cout<<" test "<<d(gen)<<std::endl;
    }
    noisyFitVec.resize(nbScales);

    std::srand(0);
    // GET BOUNDS
    std::vector<double> fitCoordsBounds(2*FitDimension);
    for(int iDim=0;iDim<FitDimension; iDim++){
      fitCoordsBounds[2*iDim]=std::numeric_limits<double>::max();
      fitCoordsBounds[2*iDim+1]=std::numeric_limits<double>::lowest();
    }
    for(int iScale=0; iScale<nbScales; iScale++){
      SimplexId nbSampledPoints = sampling_indexes[iScale].size();
      noisyFitVec[iScale].resize(fitCoordsVec_[iScale].size());
      for(int iPoint=0; iPoint<nbSampledPoints; iPoint++){
        for(int iDim=0;iDim<FitDimension; iDim++){
          if(fitCoordsBounds[2*iDim]>fitCoordsVec_[iScale][FitDimension*iPoint+iDim]){
            fitCoordsBounds[2*iDim]=fitCoordsVec_[iScale][FitDimension*iPoint+iDim];
          }
          if(fitCoordsBounds[2*iDim+1]<fitCoordsVec_[iScale][FitDimension*iPoint+iDim]){
            fitCoordsBounds[2*iDim+1]=fitCoordsVec_[iScale][FitDimension*iPoint+iDim];
          }
        }
      }
    } // end for iscale
    std::vector<double> ranges(FitDimension);
    double range=0;
    for(int iDim=0;iDim<FitDimension; iDim++){
      ranges[iDim] = fitCoordsBounds[2*iDim+1]-fitCoordsBounds[2*iDim];
      range+=std::pow(ranges[iDim], 2);
    }
    range = std::sqrt(range);
    range = 1;
    
    // generate random noise
    for(int iScale=0; iScale<nbScales; iScale++){
      SimplexId nbSampledPoints = sampling_indexes[iScale].size();
      for(int iPoint=0; iPoint<nbSampledPoints; iPoint++){
        std::vector<double> noises(FitDimension);
        double noiseRange=0;
        for(int iDim=0;iDim<FitDimension; iDim++){
          int r = std::rand() % (int)1e6; // random unsigned int < 1e6
          double rf = 2*(double)r*1e-6 - 1; // random double btw 1 and -1
          // rf = rf*std::pow((1-rf*rf), 2);
          // double n = rf*range*NoiseLevel;
          noises[iDim] = rf;
          noiseRange+=std::pow(noises[iDim], 2);
         }
        noiseRange=std::sqrt(noiseRange);
        for(int iDim=0;iDim<FitDimension; iDim++){
          // generate random noise
          double n = noises[iDim];
          // n = std::abs(n) > 1e-7 ? n/std::abs(n)*std::pow((1-noiseRange*noiseRange), 2) : 1.0;
          noisyFitVec[iScale][FitDimension*iPoint+iDim] = fitCoordsVec_[iScale][FitDimension*iPoint+iDim] + (double)d(gen);
          if(std::isnan(noisyFitVec[iScale][FitDimension*iPoint+iDim])){
            std::cout<<" error nan "<<std::endl;
          }
        } // end for iDim
      } // end for iPoint
    } // end for iscale

    UpdateEpsilonSampling=true;
  } // end if AddRandomNoise

  std::vector<std::vector<double>>& fitCoordsVec = AddRandomNoise ? noisyFitVec : fitCoordsVec_;

  if(UpdateEpsilonSampling){
  timings.clock(STEP::EPSILON_SAMPLING);

      switch(FitDimension){
        case 1:
         performEpsilonSampling(epsilonSamplers1, Epsilon, EpsilonNbMaxPoint, nbPointsMS, fitCoordsVec, pointsAtScaleStartIdx,filteredToLocalIdVec, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 2:
         performEpsilonSampling(epsilonSamplers2, Epsilon, EpsilonNbMaxPoint, nbPointsMS, fitCoordsVec, pointsAtScaleStartIdx, filteredToLocalIdVec,filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 3:
         performEpsilonSampling(epsilonSamplers3, Epsilon, EpsilonNbMaxPoint, nbPointsMS, fitCoordsVec, pointsAtScaleStartIdx, filteredToLocalIdVec,filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 4:
         performEpsilonSampling(epsilonSamplers4, Epsilon, EpsilonNbMaxPoint, nbPointsMS, fitCoordsVec, pointsAtScaleStartIdx, filteredToLocalIdVec,filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 5:
         performEpsilonSampling(epsilonSamplers5,
             Epsilon,
             EpsilonNbMaxPoint,
             nbPointsMS,
             fitCoordsVec,
             pointsAtScaleStartIdx,
             filteredToLocalIdVec,
             filteredIdxToEpsSampling,
             epsilonIdToFilteredId);
         std::cout<<"finished perform "<<std::endl;
          break;
        case 10:
         performEpsilonSampling(epsilonSamplers10, Epsilon, EpsilonNbMaxPoint, nbPointsMS, fitCoordsVec, pointsAtScaleStartIdx, filteredToLocalIdVec, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }

    UpdateEpsilonSampling=false;
    UpdatePers=true;
  timings.clock(STEP::EPSILON_SAMPLING);
  } // end epsilon sampling
    
  // DensitySigma = max_epsilon_*10;

    if(WriteMatrix){
      writeCoords(ExportMatrixFileName,
          nbPointsMS,
          epsilonIdToFilteredId,
          filteredToLocalIdVec,
          sampling_indexes,
          pointsAtScaleStartIdx,
          fitCoordsVec,
          FitDimension,
          coordinates);
    }




  RipsPersistence ripsMS{};
  ripsMS.setDebugLevel(debugLevel_);
  ripsMS.setThreadNumber(threadNumber_);
  if(UpdatePers){


    std::cout<<" Multiscale rips with "<<nbPointsMS<<" points"<<std::endl;

    if(Compute3DEmbedding){

      switch(FitDimension){
        case 1:
          ripsMS.compute3DEmbeddingMultiscale(epsilonIdToFilteredId,
              filteredToLocalIdVec,
              epsilonSamplers1,
              fitCoordsVec,
              mlsCoordsVec,
              resultBarycenter,
              &ttk::pointCloud::computeNDDistanceHelper,
              FitDimension,
              embeddingCoordinates);
          break;
        case 2:
          ripsMS.compute3DEmbeddingMultiscale(epsilonIdToFilteredId,
              filteredToLocalIdVec,
              epsilonSamplers2,
              fitCoordsVec,
              mlsCoordsVec,
              resultBarycenter,
              &ttk::pointCloud::computeNDDistanceHelper,
              FitDimension,
              embeddingCoordinates);
          break;
        case 3:
          ripsMS.compute3DEmbeddingMultiscale(epsilonIdToFilteredId,
              filteredToLocalIdVec,
              epsilonSamplers3,
              fitCoordsVec,
              mlsCoordsVec,
              resultBarycenter,
              &ttk::pointCloud::computeNDDistanceHelper,
              FitDimension,
              embeddingCoordinates);
          break;
        case 4:
          ripsMS.compute3DEmbeddingMultiscale(epsilonIdToFilteredId,
              filteredToLocalIdVec,
              epsilonSamplers4,
              fitCoordsVec,
              mlsCoordsVec,
              resultBarycenter,
              &ttk::pointCloud::computeNDDistanceHelper,
              FitDimension,
              embeddingCoordinates);
          break;
        case 5:
          ripsMS.compute3DEmbeddingMultiscale(epsilonIdToFilteredId,
              filteredToLocalIdVec,
              epsilonSamplers5,
              fitCoordsVec,
              mlsCoordsVec,
              resultBarycenter,
              &ttk::pointCloud::computeNDDistanceHelper,
              FitDimension,
              embeddingCoordinates);
          break;
        case 10:
          ripsMS.compute3DEmbeddingMultiscale(epsilonIdToFilteredId,
              filteredToLocalIdVec,
              epsilonSamplers10,
              fitCoordsVec,
              mlsCoordsVec,
              resultBarycenter,
              &ttk::pointCloud::computeNDDistanceHelper,
              FitDimension,
              embeddingCoordinates);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }
    }



    // writeCoordinatesForTomato("tomato_coords.txt", coordinates, normals, fitCoordsVec,
    //     mlsCoordsVec,
    //     mlsNormalsVec,
    //     densityVec,
    //     embeddingCoordinates,
    //     sampling_indexes,
    //     filteredToLocalIdVec,
    //     pointsAtScaleStartIdx,
    //     epsilonIdToFilteredId);


    timings.clock(STEP::EDGES);

    if(UseMixAllScales){
      switch(FitDimension){
        case 1:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers1, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 2:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers2, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 3:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers3, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 4:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers4, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 5:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers5, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 10:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers10, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }

    }else{
      switch(FitDimension){
        case 1:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlann(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers1, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 2:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlann(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers2, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 3:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlann(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers3, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 4:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlann(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers4, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 5:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlann(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers5, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        case 10:
          ripsMS.computeApproximatedEdgeVectorMultiScaleNanoFlann(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers10, fitCoordsVec, mlsCoordsVec, resultBarycenter, &ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }
    }

    timings.clock(STEP::EDGES);

    timings.clock(STEP::MST);

    { // initialize sizes and representants from sampling
      std::cout<<"npoints MS "<<nbPointsMS<<std::endl;
      sizes_ms.resize(nbPointsMS);
      std::fill(sizes_ms.begin(), sizes_ms.end(), 0);

      switch(FitDimension){
        case 1:
          retrieveSizeMSVector(epsilonSamplers1, sizes_ms, pointsAtScaleStartIdx, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 2:
        retrieveSizeMSVector(epsilonSamplers2, sizes_ms, pointsAtScaleStartIdx, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 3:
        retrieveSizeMSVector(epsilonSamplers3, sizes_ms, pointsAtScaleStartIdx, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 4:
        retrieveSizeMSVector(epsilonSamplers4, sizes_ms, pointsAtScaleStartIdx, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 5:
        retrieveSizeMSVector(epsilonSamplers5, sizes_ms, pointsAtScaleStartIdx, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        case 10:
        retrieveSizeMSVector(epsilonSamplers10, sizes_ms, pointsAtScaleStartIdx, filteredIdxToEpsSampling, epsilonIdToFilteredId);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }

      ripsMS.setInitialSizesCopy(sizes_ms);
    }

    ripsMS.computePersistenceDiagram(MaxEdgeLength, nbPointsMS, mergingEdges);
    timings.clock(STEP::MST);

    ripsMS.sortPersistenceDiagram(nbPointsMS);
    ripsDiagramMS=ripsMS.getDiagram();

    if(DensityMode==DENSITYMODE::GRAPH){
      ripsMSEdges = std::move(ripsMS.getEdges());
    }else{
      ripsMS.getEdges().clear();
    }

    // treeMS.buildRipsTree(ripsMS.getDiagram());
    
    UpdatePers=false;
    UpdateDensity=true;
  } // end Distance Rips

if(UpdateDensity){
  timings.clock(STEP::DENSITY);

std::vector<edgeTuple> edges;
// = mergingEdges; // copy

  switch(DensityMode) {
    case DENSITYMODE::EUCLIDEAN:

      switch(FitDimension){
        case 1:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers1, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, densityVec);
          break;
        case 2:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers2, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, densityVec);
          break;
        case 3:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers3, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, densityVec);
          break;
        case 4:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers4, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, densityVec);
          break;
        case 5:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers5, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, densityVec);
          break;
        case 10:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers10, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, densityVec);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }
      break;

    case DENSITYMODE::EUCLIDEAN_CONSTRAINED:

      switch(FitDimension){
        case 1:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers1, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ true, DensitySigma, densityVec);
          break;
        case 2:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers2, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ true, DensitySigma, densityVec);
          break;
        case 3:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers3, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ true, DensitySigma, densityVec);
          break;
        case 4:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers4, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ true, DensitySigma, densityVec);
          break;
        case 5:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers5, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ true, DensitySigma, densityVec);
          break;
        case 10:
          ripsMS.computeDensity(epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers10, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ true, DensitySigma, densityVec);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }
      break;
    case DENSITYMODE::GRAPH:

      switch(FitDimension){
        case 1:
          ripsMS.computeDensityWithGraph(nbPointsMS, ripsMSEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers1, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, sampler_fitting_.scales(),densityVec);
          break;
        case 2:
          ripsMS.computeDensityWithGraph(nbPointsMS, ripsMSEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers2, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, sampler_fitting_.scales(), densityVec);
          break;
        case 3:
          ripsMS.computeDensityWithGraph(nbPointsMS, ripsMSEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers3, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, sampler_fitting_.scales(), densityVec);
          break;
        case 4:
          ripsMS.computeDensityWithGraph(nbPointsMS, ripsMSEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers4, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, sampler_fitting_.scales(), densityVec);
          break;
        case 5:
          ripsMS.computeDensityWithGraph(nbPointsMS, ripsMSEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers5, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, sampler_fitting_.scales(), densityVec);
          break;
        case 10:
          ripsMS.computeDensityWithGraph(nbPointsMS, ripsMSEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers10, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, sampler_fitting_.scales(), densityVec);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }
      break;
    default:

      std::vector<float> scales = sampler_fitting.scales();
      for(auto& s : scales){
        s/=sqrt(diag_aabb);
      }

      switch(FitDimension){
        case 1:
          ripsMS.computeDensityWithTree(nbPointsMS, mergingEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers1, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, scales, densityVec);
          break;
        case 2:
          ripsMS.computeDensityWithTree(nbPointsMS, mergingEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers2, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, scales, densityVec);
          break;
        case 3:
          ripsMS.computeDensityWithTree(nbPointsMS, mergingEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers3, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, scales, densityVec);
          break;
        case 4:
          ripsMS.computeDensityWithTree(nbPointsMS, mergingEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers4, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, scales, densityVec);
          break;
        case 5:
          ripsMS.computeDensityWithTree(nbPointsMS, mergingEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers5, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, scales, densityVec);
          break;
        case 10:
          ripsMS.computeDensityWithTree(nbPointsMS, mergingEdges, edges, epsilonIdToFilteredId, filteredToLocalIdVec, epsilonSamplers10, fitCoordsVec, mlsCoordsVec, resultBarycenter, ttk::pointCloud::computeNDDistanceHelper, FitDimension, MaxEdgeLength, /*UseConstrainedDensity*/ false, DensitySigma, UseScaleDistanceDensity, UseEdgeDensity, scales, densityVec);
          break;
        default:
          printErr("Invalid FitDimension");
          break;
      }
  }
  timings.clock(STEP::DENSITY);

  // saveDensity();

  // } // end UpdateDensity
  timings.clock(STEP::PERSISTENCE);

RipsPersistence ripsDensity{};
ripsDensity.setDebugLevel(debugLevel_);
ripsDensity.setThreadNumber(threadNumber_);

std::cout<<"edges"<<std::endl;
if(!UseEdgeDensity){

  for(int iEdge=0; iEdge<edges.size(); iEdge++){
    auto& e = edges[iEdge];
    SimplexId id_ms_i = std::get<0>(e);
    SimplexId id_ms_j = std::get<1>(e);

    double edge_density = std::min(densityVec[id_ms_i], densityVec[id_ms_j]);
    std::get<2>(e) = -edge_density;
  }
}
for(auto& d : densityVec){
  d = -d;
}
std::cout<<"diagram"<<std::endl;
ripsDensity.setInitialSizesCopy(sizes_ms);
ripsDensity.setEdges(edges);
ripsDensity.sortEdgesAscending();


ripsDensity.computePersistenceDiagramWithBirth(-1, nbPointsMS, mergingEdgesDensity, densityVec);
ripsDensity.sortPersistenceDiagram(nbPointsMS);

timings.clock(STEP::PERSISTENCE);

diagramMS = ripsDensity.getDiagram();

const auto ripsEdges = ripsDensity.getEdges();

for(auto& pair : diagramMS){
  const auto& birthId = pair.birth.id;
  const auto& deathId = pair.death.id;
  int iScale = findScaleId(birthId, pointsAtScaleStartIdx);
  SimplexId id_eps = birthId - pointsAtScaleStartIdx[iScale];
  SimplexId id_filtered = epsilonIdToFilteredId[iScale][id_eps];
  SimplexId id_local = filteredToLocalIdVec[iScale][id_filtered];
  pair.birth.coords[0] = fitCoordsVec[iScale][FitDimension*id_local];
  pair.birth.coords[1] = fitCoordsVec[iScale][FitDimension*id_local+1];
  pair.birth.coords[2] = 0;

  // std::cout<<"death Id "<<deathId<<" outta "<<ripsEdges.size()<<std::endl;
  if(deathId > 0){

    auto& e = ripsEdges[deathId];
    SimplexId id_ms_i = std::get<0>(e);
    SimplexId id_ms_j = std::get<1>(e);
    SimplexId id_ms_death = densityVec[id_ms_i] > densityVec[id_ms_j] ? id_ms_i : id_ms_j;

    int scale_death = findScaleId(id_ms_death, pointsAtScaleStartIdx);
    SimplexId id_eps_death = id_ms_death - pointsAtScaleStartIdx[scale_death];
    SimplexId id_filtered_death = epsilonIdToFilteredId[scale_death][id_eps_death];
    SimplexId id_local_death = filteredToLocalIdVec[scale_death][id_filtered_death];
    pair.death.coords[0] = fitCoordsVec[scale_death][FitDimension*id_local_death];
    pair.death.coords[1] = fitCoordsVec[scale_death][FitDimension*id_local_death+1];
    pair.death.coords[2] = 0;
  }
}

// treeDensity.buildRipsTree(ripsDensity)|
// RipsTree treeD{};
// std::cout<<"tree density"<<std::endl;
// treeD.buildRipsTree(ripsDensity.getDiagram());
std::cout<<"\n\ntree ms"<<std::endl;

timings.clock(STEP::RIPS_TREE);
treeMS.buildRipsTree(ripsDensity.getDiagram());
treeMS.buildAABBDiagonals(epsilonIdToFilteredId, filteredToLocalIdVec, pointsAtScaleStartIdx, fitCoordsVec, FitDimension);
timings.clock(STEP::RIPS_TREE);


UpdateDensity = false;
UpdateDisplay=true;
} // end Update Density


  std::cout<<"\n Now prepare the outputs"<<std::endl;

  UpdateDisplay = true;
  if(UpdateDisplay){

    // printMsg("Go for most persistent components");

    // outputScaleMB ->SetNumberOfBlocks(ScaleSamplingCount);
    // for(int iScale=0; iScale<ScaleSamplingCount; iScale++){

    //   //               float* fitCoordinates = fitCoordsVec[iScale].data();
    //   //               float* density = densityVec[iScale].data();
    //   //               int* isStable = isStableVec[iScale].data();
    //   //               float* mlsCoordinates = mlsCoordsVec[iScale].data();


    //   //               // std::cout<<" scale "<<iScale<<std::endl;
    //   //               vtkNew<vtkMultiBlockDataSet> outputMB{};
    //   //               outputScaleMB->SetBlock(iScale, outputMB);
    //   //               std::vector<std::vector<double>>& distMat = distanceMatrix[iScale];
    //   //               std::vector<SimplexId>& filteredToLocalId = filteredToLocalIdVec[iScale];

    //   //               if(!filteredToLocalId.empty()){
    //   //                 trees_[iScale].computeArcAutoSimilarity(distMat,filteredToLocalId);
    //   //               }
    //   //               // std::cout<<"t0"<<std::endl;

    //   //               const auto& ripsArcs = trees_[iScale].getRipsArcs();
    //   //               SimplexId nbOutputPrimitives = NbPrimitives;
    //   //               if(nbOutputPrimitives>ripsArcs.size()){
    //   //                 nbOutputPrimitives = ripsArcs.size();
    //   //               }
    //   //               outputMB->SetNumberOfBlocks(nbOutputPrimitives);
    //   //               std::cout<<"t1"<<std::endl;

    //   //               double scale = sampler_fitting_.scale(iScale);

    //   //               for(SimplexId iArc=0; iArc<nbOutputPrimitives; iArc++){


    //   //                   // std::cout<<"    arc "<<iArc<<"/"<<nbOutputPrimitives<<std::endl;
    //   //                   auto ccId = ripsArcs[iArc].ccId;
    //   //                   auto threshold = ripsArcs[iArc].death;
    //   //                   auto persistence = ripsArcs[iArc].death-ripsArcs[iArc].birth;
    //   //                   vtkNew<vtkUnstructuredGrid> block;

    //   //                   createOutputConnectedComponent(block, trees_[iScale], ccId, threshold,  coordinates, normals, fitCoordinates, isStable, density, SizeMin, sampling_indexes[iScale]);

    //   //                   vtkNew<vtkFloatArray> mlsCoordsArray;
    //   //                   mlsCoordsArray->SetName("mlsCoords");
    //   //                   mlsCoordsArray->SetNumberOfComponents(3);
    //   //                   mlsCoordsArray->SetNumberOfTuples(block->GetNumberOfPoints());
    //   //                   block->GetPointData()->AddArray(mlsCoordsArray);

    //   //                     std::vector<SimplexId> pointIds;
    //   //                     trees_[iScale].getComponentPointIdsAtThreshold(ccId, threshold, pointIds, SizeMin);

    //   //                   for(SimplexId id=0; id<block->GetNumberOfPoints(); id++){
    //   //                     SimplexId iPoint = pointIds[id];
    //   //                     mlsCoordsArray->SetTuple(id, mlsCoordinates+3*sampling_indexes[iScale][iPoint]);
    //   //                   }

    //   //                   if(Compute3DEmbedding){
    //   //                       vtkNew<vtkFloatArray> embeddingArray;
    //   //                       embeddingArray->SetName("3dEmbeddingCoordinates");
    //   //                       embeddingArray->SetNumberOfComponents(3);
    //   //                       embeddingArray->SetNumberOfTuples(block->GetNumberOfPoints());
    //   //                       block->GetPointData()->AddArray(embeddingArray);

    //   //                     std::vector<SimplexId> pointIds;
    //   //                     trees_[iScale].getComponentPointIdsAtThreshold(ccId, threshold, pointIds, SizeMin);

    //   //                       for(SimplexId id=0; id<block->GetNumberOfPoints(); id++){
    //   //                       SimplexId iPoint = pointIds[id];
    //   //                       embeddingArray->SetTuple3(id, embeddingCoordinates[iScale][0][iPoint], embeddingCoordinates[iScale][1][iPoint], embeddingCoordinates[iScale][2][iPoint]);
    //   //                       }


    //   //                   }
    //   //                   int dummy = -1;
    //   //                   setSingleComponentAttributes(block, dummy, persistence, iArc, -1, scale, ripsArcs[iArc].autosimilarity);

    //   //                   outputMB->SetBlock(iArc, block);

    //   //               }

    // }
    UpdateDisplay=false;


    vtkUnstructuredGrid* pd;
    // PERSISTENCE DIAGRAM OUTPUT
    vtkMultiBlockDataSet *outputPDMB = vtkMultiBlockDataSet::GetData(outputVector,2);

    outputPDMB->SetNumberOfBlocks(2);
    {

      // for(int iScale = 0; iScale<nbScales; iScale++){

      //   //         float* fitCoordinates = fitCoordsVec[iScale].data();
      //   //         const auto& tree = trees_[iScale];

      //   //         if(0 /*MergeTree*/){


      //   //         }else{
      vtkNew<vtkUnstructuredGrid> outputArcPD{};
      pd = outputArcPD;
      outputPDMB->SetBlock(0,outputArcPD);


      //   //           // arcs to persistent pairs
      vtkNew<vtkDoubleArray> scalars{};
      scalars->SetNumberOfComponents(1);
      DiagramToVTU(outputArcPD, diagramMS, scalars, *this, 3, 0);
    }
    vtkUnstructuredGrid* pd2;

    // PERSISTENCE DIAGRAM OUTPUT RIPS
    {
      vtkMultiBlockDataSet *outputPDMB2 = vtkMultiBlockDataSet::GetData(outputVector,4);

      outputPDMB2->SetNumberOfBlocks(1);

      vtkNew<vtkUnstructuredGrid> outputArcPD2{};
      pd2 = outputArcPD2;
      outputPDMB2->SetBlock(0,outputArcPD2);
      vtkNew<vtkDoubleArray> scalars{};
      scalars->SetNumberOfComponents(1);
      DiagramToVTU(outputArcPD2, ripsDiagramMS, scalars, *this, 3, 0);
    }
    // }

    // std::cout<<"diagram is done"<<std::endl;

    // {
    // vtkMultiBlockDataSet *outputMTMB = vtkMultiBlockDataSet::GetData(outputVector,3);
    // outputMTMB->SetNumberOfBlocks(2);

    // std::cout<<"test1"<<std::endl;
    // const auto& ripsArcs = trees_[0].getRipsArcs();
    // SimplexId nbArcs = ripsArcs.size();
    // std::cout<<"test2"<<std::endl;

    // // nodes
    // vtkNew<vtkUnstructuredGrid> vtuNodes{};
    // vtuNodes->Allocate();
    // outputMTMB->SetBlock(0, vtuNodes);
    // std::cout<<"here0"<<std::endl;
    // vtkNew<vtkPoints> pointsNodes{};
    // pointsNodes->SetNumberOfPoints(nbArcs+1);
    // vtuNodes->SetPoints(pointsNodes);
    // std::cout<<"here1"<<std::endl;

    // vtkNew<vtkDoubleArray> nodeScalarArray{};
    // nodeScalarArray->SetNumberOfComponents(1);
    // nodeScalarArray->SetNumberOfTuples(nbArcs+1);
    // nodeScalarArray->SetName("Scalar");
    // vtuNodes->GetPointData()->AddArray(nodeScalarArray);

    // vtkNew<vtkIntArray> nodeIdArray{};
    // nodeIdArray->SetNumberOfComponents(1);
    // nodeIdArray->SetNumberOfTuples(nbArcs+1);
    // nodeIdArray->SetName("NodeId");
    // vtuNodes->GetPointData()->AddArray(nodeIdArray);

    // vtkNew<vtkIntArray> nodeCriticalTypeArray{};
    // nodeCriticalTypeArray->SetNumberOfComponents(1);
    // nodeCriticalTypeArray->SetNumberOfTuples(nbArcs+1);
    // nodeCriticalTypeArray->SetName("CriticalType");
    // vtuNodes->GetPointData()->AddArray(nodeCriticalTypeArray);

    // std::cout<<"here2"<<std::endl;
    // // edges
    // vtkNew<vtkUnstructuredGrid> vtuEdges{};
    // vtuEdges->Allocate();
    // vtkNew<vtkPoints> pointsEdges{};
    // vtuEdges->SetPoints(pointsNodes);
    // outputMTMB->SetBlock(1, vtuEdges);
    // std::cout<<"here3"<<std::endl;

    // vtkNew<vtkIntArray> edgeUpNodeIdArray{};
    // edgeUpNodeIdArray->SetNumberOfComponents(1);
    // edgeUpNodeIdArray->SetNumberOfTuples(nbArcs);
    // edgeUpNodeIdArray->SetName("upNodeId");
    // vtuEdges->GetCellData()->AddArray(edgeUpNodeIdArray);

    // vtkNew<vtkIntArray> edgeDownNodeIdArray{};
    // edgeDownNodeIdArray->SetNumberOfComponents(1);
    // edgeDownNodeIdArray->SetNumberOfTuples(nbArcs);
    // edgeDownNodeIdArray->SetName("downNodeId");
    // vtuEdges->GetCellData()->AddArray(edgeDownNodeIdArray);
    // vtuEdges->GetPointData()->AddArray(nodeScalarArray);


    // pointsEdges->SetNumberOfPoints(nbArcs+1);
    // std::cout<<"nb arcs "<<nbArcs<<std::endl;
    // // for(SimplexId iArc=nbArcs-1; iArc>=0; iArc--){
    // for(SimplexId iArc=0; iArc<nbArcs; iArc++){
    //   SimplexId i = iArc;
    //   const auto& arc = ripsArcs[iArc];
    //   std::cout<<"ARC "<<arc.birthId<<"-"<<arc.deathId<<" ("<<arc.birth<<"-"<<arc.death<<")"<<std::endl;
    //   SimplexId pointId = trees_[0].getPointId(arc.ccId);
    //   pointsNodes->SetPoint(arc.birthId, coordinates[3*pointId],coordinates[3*pointId+1],coordinates[3*pointId+2]);
    //   pointsNodes->SetPoint(arc.deathId, coordinates[3*pointId],coordinates[3*pointId+1],coordinates[3*pointId+2]);
    //   int criticalType = (arc.type == 0)  ? 0 : 1;
    //   nodeCriticalTypeArray->SetTuple1(arc.birthId,criticalType);
    //   nodeIdArray->SetTuple1(arc.birthId, arc.birthId);
    //   nodeIdArray->SetTuple1(arc.deathId, arc.deathId);
    //   nodeScalarArray->SetTuple1(arc.birthId, arc.birth);
    //   nodeScalarArray->SetTuple1(arc.deathId, arc.death);

    //   // if(iArc==0){
    //   // vtuNodes->GetPoints()->SetPoint(nbArcs, 0,0,0);
    //   // nodeCriticalTypeArray->SetTuple1(nbArcs,3);
    //   // nodeIdArray->SetTuple1(nbArcs, arc.deathId);
    //   // std::cout<<"set death id "<<arc.deathId<<std::endl;
    //   // nodeScalarArray->SetTuple1(nbArcs, arc.death);
    //   // }


    //   // vtkIdType id0 = pointsEdges->InsertNextPoint(0,0,0);
    //   // vtkIdType id1 = pointsEdges->InsertNextPoint(0,0,0);
    //   // std::cout<<"id0 "<<id0<<" vs iArc "<<iArc<<std::endl;
    //   vtkNew<vtkIdList> idList{};
    //   idList->InsertNextId(arc.birthId);
    //   idList->InsertNextId(arc.deathId);
    //   vtuEdges->InsertNextCell(VTK_LINE, idList);
    //   edgeUpNodeIdArray->SetTuple1(i, arc.deathId);
    //   edgeDownNodeIdArray->SetTuple1(i, arc.birthId);

    // }
    // // pointsNodes->SetPoint(nbArcs, 0,0,0);
    // ttkUtils::CellVertexFromPoints(vtuNodes, pointsNodes);



    // }

    //  DISTANCE MATRIX OUTPUT
    {
      // std::cout<<"GO FOR MATRIX"<<std::endl;

      // std::cout<<"Go for matrix"<<std::endl;
      // vtkMultiBlockDataSet *outputMatrixMB = vtkMultiBlockDataSet::GetData(outputVector,4);
      // outputMatrixMB->SetNumberOfBlocks(ScaleSamplingCount);
      // for(int iScale=0; iScale<ScaleSamplingCount; iScale++){

      // vtkNew<vtkPolyData> outputMatrix;
      // outputMatrixMB->SetBlock(iScale, outputMatrix);

      // std::vector<std::vector<double>>& distMat = distanceMatrix[iScale];
      // std::vector<SimplexId>& filteredToLocalId = filteredToLocalIdVec[iScale];


      // // const auto nRows = static_cast<size_t>(distMat.size());
      // const auto nRows = static_cast<size_t>(sampling_indexes[iScale].size());

      // // grid points
      // vtkNew<vtkPoints> points{};
      // points->SetNumberOfPoints((nRows + 1) * (nRows + 1));
      // // grid cells connectivity arrays
      // vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
      // offsets->SetNumberOfComponents(1);
      // offsets->SetNumberOfTuples(nRows * nRows + 1);
      // connectivity->SetNumberOfComponents(1);
      // connectivity->SetNumberOfTuples(4 * nRows * nRows);
      // // distance and proximity arrays
      // vtkNew<vtkDoubleArray> dist{}, prox{};
      // dist->SetNumberOfComponents(1);
      // dist->SetName("Distance");
      // dist->SetNumberOfTuples(nRows * nRows);
      // prox->SetNumberOfComponents(1);
      // prox->SetName("Proximity");
      // prox->SetNumberOfTuples(nRows * nRows);

      // vtkNew<vtkIntArray> idArray{};
      // idArray->SetNumberOfComponents(1);
      // idArray->SetNumberOfTuples(nRows*nRows);
      // idArray->SetName("id");

      // #ifdef TTK_ENABLE_OPENMP
      // #pragma omp parallel for num_threads(this->threadNumber_)
      // #endif // TTK_ENABLE_OPENMP
      // for(size_t i = 0; i < nRows + 1; ++i) {
      //   for(size_t j = 0; j < nRows + 1; ++j) {
      //     points->SetPoint(i * (nRows + 1) + j, i, j, 0.0);
      //   }
      // }



      // std::vector<SimplexId> matrixOrdering(nRows);

      // switch(MatrixOrdering){

      //   case 1: // Rips Merge Order
      //     {
      //       std::cout<<"Rips Matrix"<<std::endl;

      //       auto& tree = UseSmallScaleMatrixOrdering ? trees_.back() : trees_[iScale];
      //       const std::vector<SimplexId>& mergeOrder = tree.getMergeOrder();
      // #ifdef TTK_ENABLE_OPENMP
      // #pragma omp parallel for num_threads(this->threadNumber_)
      // #endif // TTK_ENABLE_OPENMP
      //       for(SimplexId i=0; i<nRows; i++){
      //         matrixOrdering[i] = tree.getComponent(mergeOrder[i]).birthId();
      //       }
      //       std::cout<<"rips ordering done"<<std::endl;
      //       break;
      //     }

      //   case 2: // external order
      //     {
      //       std::ifstream input(ExternalLabelFile);
      //       std::cout<<"Ordering matrix from external file"<<std::endl;
      //       if(!input){
      //         printWrn("External order file not found, default ordering used");
      //         std::iota(matrixOrdering.begin(), matrixOrdering.end(),0);
      //         break;
      //       }
      //       std::vector<SimplexId> labels{};
      //       double lbl=-2;
      //       while(input >> lbl){
      //         labels.push_back(lbl);
      //         std::cout<<" "<<lbl;
      //       }
      //       std::cout<<std::endl;
      //       std::cout<<"label vec size "<<labels.size()<<" vs nRows "<<nRows<<std::endl;
      //       std::iota(matrixOrdering.begin(), matrixOrdering.end(),0);
      //       std::sort(matrixOrdering.begin(), matrixOrdering.end(), 
      //        [=](const auto& i, const auto& j) { return labels[i]< labels[j]; } );
      //       break;
      //     }

      //   default: // original point order
      //     std::cout<<"original order"<<std::endl;
      //     std::iota(matrixOrdering.begin(), matrixOrdering.end(),0);
      //     break;


      // }

      // std::cout<<"matrix ordering size "<<matrixOrdering.size()<<std::endl;

      // #ifdef TTK_ENABLE_OPENMP
      // #pragma omp parallel for num_threads(this->threadNumber_)
      // #endif // TTK_ENABLE_OPENMP
      // for(size_t i = 0; i < nRows; ++i) {
      //   for(size_t j = 0; j < nRows; ++j) {
      //     // build the cell connectivity and offsets array
      //     const auto nptline{static_cast<vtkIdType>(nRows + 1)};
      //     const auto curr{static_cast<vtkIdType>(i * nptline + j)};
      //     const auto o = i * nRows + j;
      //     connectivity->SetTuple1(4 * o + 0, curr);
      //     connectivity->SetTuple1(4 * o + 1, curr + nptline);
      //     connectivity->SetTuple1(4 * o + 2, curr + nptline + 1);
      //     connectivity->SetTuple1(4 * o + 3, curr + 1);
      //     offsets->SetTuple1(o, 4 * o);

      //     // copy distance matrix to heat map cell data
      //     const auto val = distMat[filteredToLocalId[matrixOrdering[i]]][filteredToLocalId[matrixOrdering[j]]];

      //     dist->SetTuple1(i * nRows + j, val);
      //     prox->SetTuple1(i * nRows + j, std::exp(-val));
      //     idArray->SetTuple1(i*nRows + j, matrixOrdering[j]);
      //   }
      // }
      // offsets->SetTuple1(nRows * nRows, connectivity->GetNumberOfTuples());

      // // gather arrays to make the PolyData
      // vtkNew<vtkCellArray> cells{};
      // cells->SetData(offsets, connectivity);
      // outputMatrix->SetPoints(points);
      // outputMatrix->SetPolys(cells);
      // outputMatrix->GetCellData()->AddArray(dist);
      // outputMatrix->GetCellData()->AddArray(prox);
      // outputMatrix->GetCellData()->AddArray(idArray);
      // }
      // std::cout<<"end distance matrices"<<std::endl;
    }

    // TRANSITION MATRIX OUTPUT
    {
      //std::cout<<"GO FOR TRANSITION MATRICES"<<std::endl;
      //vtkMultiBlockDataSet *outputMatrixMB = vtkMultiBlockDataSet::GetData(outputVector,5);
      //outputMatrixMB->SetNumberOfBlocks(ScaleSamplingCount-1);

      //for(int iScale=0; iScale<ScaleSamplingCount-1; iScale++){
      // std::cout<<"scale "<<iScale<<" to "<<iScale+1<<std::endl;
      // int test_count=0;
      // vtkNew<vtkPolyData> outputMatrix;
      // outputMatrixMB->SetBlock(iScale, outputMatrix);

      // std::vector<SimplexId>& filteredToLocalId0 = filteredToLocalIdVec[iScale];
      // std::vector<SimplexId>& filteredToLocalId1 = filteredToLocalIdVec[iScale+1];

      // std::cout<<" test count "<<++test_count<<std::endl;

      // // std::vector<SimplexId>& localIdToGlobalId0 = sampling_indexes[iScale];
      // // std::vector<SimplexId>& localIdToGlobalId1 = sampling_indexes[iScale+1];

      // std::vector<float>& fitCoords0 = fitCoordsVec[iScale];
      // std::vector<float>& fitCoords1 = fitCoordsVec[iScale+1];
      // std::cout<<" test count "<<++test_count<<std::endl;

      // const auto nRows = static_cast<size_t>(sampling_indexes[iScale].size());
      // auto& tree0 = trees_[iScale];
      // auto& tree1 = trees_[iScale+1];
      // std::cout<<" test count "<<++test_count<<std::endl;

      // const std::vector<SimplexId>& mergeOrder0 = tree0.getMergeOrder();
      // const std::vector<SimplexId>& mergeOrder1 = tree1.getMergeOrder();
      // std::cout<<" test count "<<++test_count<<std::endl;

      // // grid points
      // vtkNew<vtkPoints> points{};
      // points->SetNumberOfPoints((nRows + 1) * (nRows + 1));
      // // grid cells connectivity arrays
      // vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
      // offsets->SetNumberOfComponents(1);
      // offsets->SetNumberOfTuples(nRows * nRows + 1);
      // connectivity->SetNumberOfComponents(1);
      // connectivity->SetNumberOfTuples(4 * nRows * nRows);
      // // distance and proximity arrays
      // vtkNew<vtkDoubleArray> dist{}, prox{};
      // dist->SetNumberOfComponents(1);
      // dist->SetName("Distance");
      // dist->SetNumberOfTuples(nRows * nRows);
      // prox->SetNumberOfComponents(1);
      // prox->SetName("Proximity");
      // prox->SetNumberOfTuples(nRows * nRows);

      // //ids array
      //vtkNew<vtkIntArray> idArray0{};
      //idArray0->SetNumberOfComponents(1);
      //idArray0->SetNumberOfTuples(nRows*nRows);
      //idArray0->SetName("id0");

      //vtkNew<vtkIntArray> idArray1{};
      //idArray1->SetNumberOfComponents(1);
      //idArray1->SetNumberOfTuples(nRows*nRows);
      //idArray1->SetName("id1");

      //#ifdef TTK_ENABLE_OPENMP
      //#pragma omp parallel for num_threads(this->threadNumber_)
      //#endif // TTK_ENABLE_OPENMP
      //for(size_t i = 0; i < nRows + 1; ++i) {
      //  for(size_t j = 0; j < nRows + 1; ++j) {
      //    points->SetPoint(i * (nRows + 1) + j, i, j, 0.0);
      //  }
      //}

      // std::vector<SimplexId> matrixOrdering0(nRows);
      // std::vector<SimplexId> matrixOrdering1(nRows);
      //#ifdef TTK_ENABLE_OPENMP
      //#pragma omp parallel for num_threads(this->threadNumber_)
      //#endif // TTK_ENABLE_OPENMP
      //     for(SimplexId i=0; i<nRows; i++){
      //       matrixOrdering0[i] = tree0.getComponent(mergeOrder0[i]).birthId();
      //       matrixOrdering1[i] = tree1.getComponent(mergeOrder1[i]).birthId();
      //     }

      // // std::vector<std::vector<double>>& transMat;
      // // transMat.resize(nRows, std::vector<double>(nRow));
      // // std::cout<<" test count "<<++test_count<<std::endl;

      // // for(int i=0; i<fitCoords0.size(); i++){
      // //   std::cout<<" "<<fitCoords0[i];
      // // }
      // // std::cout<<std::endl;

      // //      std::cout<<"CHECKING FIT COORDS 2"<<std::endl;

      // // for(auto i : sampling_indexes[iScale]){
      // //   std::cout<<" point "<<i<<": ";
      // //   for(int j=0; j<5; j++){
      // //     std::cout<<" "<<(fitCoords0.data())[5*i+j];
      // //   }
      // //   std::cout<<std::endl;
      // // }

      // // std::cout<<"Checking tree birth ids"<<std::endl;
      // // std::cout<<"merge order size "<<mergeOrder0.size()<<std::endl;
      // // for(int i=0; i<mergeOrder0.size(); i++){
      // //   std::cout<<i<<" "<<mergeOrder0[i]<<" "<<tree0.getComponent(mergeOrder0[i]).birthId()<<std::endl;
      // // }
      // // std::cout<<"checking matching indices"<<mergeOrder0.size()<<std::endl;
      // // std::vector<SimplexId> test{};
      // // for(int i=0; i<mergeOrder0.size(); i++){
      // //   test.push_back(sampling_indexes[iScale][tree0.getComponent(mergeOrder0[i]).birthId()]);
      // // }
      // // std::sort(test.begin(), test.end());
      // // for(int i=0; i<mergeOrder0.size(); i++){
      // //   std::cout<<i<<" "<<test[i]<<std::endl;
      // // }

      // // std::cout<<"sampling size "<<sampling_indexes[iScale].size()<<std::endl;
      // // for(int i=0; i<sampling_indexes[iScale].size(); i++){
      //   // std::cout<<i<<" "<<sampling_indexes[iScale][i]<<std::endl;
      // // }


      //#ifdef TTK_ENABLE_OPENMP
      //#pragma omp parallel for num_threads(this->threadNumber_)
      //#endif // TTK_ENABLE_OPENMP
      // for(size_t i = 0; i < nRows; ++i) {


      //   // transition matrix values
      //   // transMat[iRow].resize(nRows);
      //   float* p0 = fitCoords0.data()+5*sampling_indexes[iScale][tree0.getComponent(mergeOrder0[i]).birthId()];
      //   for(size_t j = 0; j < nRows; ++j) {
      //     float* p1 = fitCoords1.data()+5*sampling_indexes[iScale+1][tree1.getComponent(mergeOrder1[j]).birthId()];

      //     // build the cell connectivity and offsets array
      //     const auto nptline{static_cast<vtkIdType>(nRows + 1)};
      //     const auto curr{static_cast<vtkIdType>(i * nptline + j)};
      //     const auto o = i * nRows + j;
      //     connectivity->SetTuple1(4 * o + 0, curr);
      //     connectivity->SetTuple1(4 * o + 1, curr + nptline);
      //     connectivity->SetTuple1(4 * o + 2, curr + nptline + 1);
      //     connectivity->SetTuple1(4 * o + 3, curr + 1);
      //     offsets->SetTuple1(o, 4 * o);

      //     // transMat[iRow][iCol]=pointDistance(p0,p1);
      //     const double val =  pointDistance(p0,p1,5);
      //     // std::cout<<"val "<<val<<std::endl;
      //     dist->SetTuple1(i * nRows + j, val);
      //     prox->SetTuple1(i * nRows + j, std::exp(-val));
      //     idArray0->SetTuple1(i*nRows + j, matrixOrdering0[i]);
      //     idArray1->SetTuple1(i*nRows + j, matrixOrdering1[j]);
      //   }
      // }
      // // std::cout<<" test count "<<++test_count<<std::endl;

      // offsets->SetTuple1(nRows * nRows, connectivity->GetNumberOfTuples());

      // // gather arrays to make the PolyData
      // vtkNew<vtkCellArray> cells{};
      // cells->SetData(offsets, connectivity);
      // outputMatrix->SetPoints(points);
      // outputMatrix->SetPolys(cells);
      // outputMatrix->GetCellData()->AddArray(dist);
      // outputMatrix->GetCellData()->AddArray(prox);
      // outputMatrix->GetCellData()->AddArray(idArray0);
      // outputMatrix->GetCellData()->AddArray(idArray1);
      // // std::cout<<" test count "<<++test_count<<std::endl;

      //} // end for iscale
      //std::cout<<"end distance matrices"<<std::endl;
    }

		// RIPS COMPLEX OUTPUT
		if(OutputSimplicialComplex) {
      std::cout<<"Create Rips Complex"<<std::endl;;
      // std::cout<<mergingEdges.size()<<" merging edges"<<std::endl;
      vtkPolyData* outputRC = vtkPolyData::GetData(outputVector, 0);

        switch(FitDimension){
          case 1:
            createOutputMergingEdges(outputRC, mergingEdgesDensity, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, densityVec, embeddingCoordinates, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers1, epsilonIdToFilteredId);
            break;
          case 2:
            createOutputMergingEdges(outputRC, mergingEdgesDensity, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, densityVec, embeddingCoordinates, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers2, epsilonIdToFilteredId);
            break;
          case 3:
            createOutputMergingEdges(outputRC, mergingEdgesDensity, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, densityVec, embeddingCoordinates, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers3, epsilonIdToFilteredId);
            break;
          case 4:
            createOutputMergingEdges(outputRC, mergingEdgesDensity, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, densityVec, embeddingCoordinates, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers4, epsilonIdToFilteredId);
            break;
          case 5:
            createOutputMergingEdges(outputRC, mergingEdgesDensity, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, densityVec, embeddingCoordinates, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers5, epsilonIdToFilteredId);
            break;
          case 10:
            createOutputMergingEdges(outputRC, mergingEdgesDensity, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, densityVec, embeddingCoordinates, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers10, epsilonIdToFilteredId);
            break;
          default:
            printErr("Invalid FitDimension");
            break;
        }
    }

      PrimitiveSegmentation pseg{};

    // RIPS TREE MULTISCALE OUTPUT
    {
      printMsg("Rips Tree multiscale");
      // treeMS.buildRipsArcs(SizeMin);
      // timings.clock(STEP::RIPS_ARCS);
      // treeMS.buildRipsArcsWithBirth(SizeMin, -1, epsilonIdToFilteredId, filteredToLocalIdVec, pointsAtScaleStartIdx, fitCoordsVec, densityVec, FitDimension);
      // timings.clock(STEP::RIPS_ARCS);
      timings.clock(STEP::PROCESS_ARCS);
      // std::cout<<"built ripsarc"<<std::endl;
      // vtkNew<vtkMultiBlockDataSet> outputMB{}
      vtkMultiBlockDataSet* outputMB = vtkMultiBlockDataSet::GetData(outputVector,6);

      vtkNew<vtkFloatArray> barycenterFDArray{};
      barycenterFDArray->SetNumberOfComponents(3);
      barycenterFDArray->SetNumberOfTuples(1);
      barycenterFDArray->SetTuple(0, resultBarycenter.data());
      barycenterFDArray->SetName("basisCenter");
      outputMB->GetFieldData()->AddArray(barycenterFDArray);


      // std::vector<RipsArc> ripsArcs = std::vector<RipsArc>(treeMS.getRipsArcs());
      SimplexId nbOutputPrimitives = NbPrimitives;
      if(nbOutputPrimitives>treeMS.getNumberOfPoints()){
        nbOutputPrimitives = treeMS.getNumberOfPoints();
      }



      // find scale ids corresponding to thresholds
      int scaleIdThresholdMin=-1;
      int scaleIdThresholdMax=-1;
      {
        const auto& scales = sampler_fitting.scales();
        for(int iScale=0; iScale < nbScales; iScale++){
          if(scaleIdThresholdMin == -1 and scales[iScale]> ScaleThresholdMin*sqrt(diag_aabb)){
            scaleIdThresholdMin = iScale;
          }
          int jScale = nbScales-1 - iScale;
          if(scaleIdThresholdMax == -1 and scales[jScale]< ScaleThresholdMax*sqrt(diag_aabb)){
            scaleIdThresholdMax = jScale;
          }
        }
        if(scaleIdThresholdMin==-1) {
          scaleIdThresholdMin = 0; 
          printWrn("Wrong scalemin threshold, ignoring");
        }
        if(scaleIdThresholdMax==-1) {
          scaleIdThresholdMin = nbScales-1; 
          printWrn("Wrong scalemax threshold, ignoring");
        }
      }

      std::vector<SimplexId> ccIdToArc(treeMS.getNumberOfPoints(), -1);
      std::vector<double> thresholds(nbOutputPrimitives, -1);
      std::vector<double> productPersistence(nbOutputPrimitives, -1);

      double maxProductPersistence=-1;

      // ### GET CORRECT THRESHOLDS FOR DIAGRAMS
      for(SimplexId iCC=0; iCC < nbOutputPrimitives; iCC++){
        const auto& cc = treeMS.getComponent(iCC);
        auto birth = cc.birth();
        auto death = cc.death();

        // ccIdToArc[iCC] = arcCount;

        thresholds[iCC] = death;

        SimplexId iCCParent = cc.parent();

        // thresholds[iCCParent] = std::min(death, thresholds[iCCParent]);



        // find threshold that garantees a small diag
        auto& diags = cc.diag_aabb_5D();

        if(DiagMax>=0){ 
          bool found_good_diag=false;

          for(size_t iDiag=0; iDiag<diags.size(); iDiag++) {
            auto diag = diags[iDiag];
            if(diag<DiagMax){
              thresholds[iCC] = std::min(treeMS_.getComponent(cc.getChild(iDiag)).death(), thresholds[iCC]);
              found_good_diag=true;
              break;
            }
          }
          if(!found_good_diag){ // set threshold to birth point (diag = 0), as the introduction of the first child already induces diag>diagmax
            thresholds[iCC] = birth;
          }
        }
      }

      // find scale range per component
      std::vector<std::pair<int,int>> scalePairs(nbOutputPrimitives, {-1,-1});



      // COMPUTE HIGHEST COMPONENT FOR EACH POINT AND SCORE
      std::vector<float> scores(nbOutputPrimitives,0);
      std::vector<SimplexId> highestPersComponent(nbPoints,-1);
      {
        const auto& scales = sampler_fitting.scales();
        std::vector<SimplexId> sizesAtThreshold(nbOutputPrimitives, 0);
        std::vector<SimplexId> vote(nbOutputPrimitives,0);

        for(SimplexId iCC=0; iCC < nbOutputPrimitives; iCC++){


          auto& cc = treeMS.getComponent(iCC);

          double threshold = thresholds[iCC];
          double persistence = cc.death()-cc.birth();
          double birth = cc.birth();

          if(birth > BirthMax or persistence < PersistenceMin){
            continue;
          }

          // compute votes
          // find points
          std::vector<SimplexId> pointIds;
          treeMS.getComponentPointIdsAtThreshold(iCC, threshold, pointIds,0); // multiscale ids


          std::vector<SimplexId> globalIds;
          int minScale = std::numeric_limits<int>::max();
          int maxScale = std::numeric_limits<int>::lowest();

          for(auto id_ms : pointIds){
            int iScale = findScaleId(id_ms, pointsAtScaleStartIdx);
            SimplexId id_eps = id_ms-pointsAtScaleStartIdx[iScale]; // epsilon id
            SimplexId id_filtered = epsilonIdToFilteredId[iScale][id_eps];
            SimplexId id_local = filteredToLocalIdVec[iScale][id_filtered];
            SimplexId id_global = sampling_indexes[iScale][id_local];

            globalIds.push_back(id_global);

            if(iScale > maxScale){
              maxScale = iScale;
            }
            if(iScale < minScale){
              minScale = iScale;
            }


            // add also children

            const std::vector<SimplexId>* pchildren{};
            switch(FitDimension){
              case 1:
                pchildren = &(epsilonSamplers1[iScale].getChildren(id_eps));
                break;
              case 2:
                pchildren = &(epsilonSamplers2[iScale].getChildren(id_eps));
                break;
              case 3:
                pchildren = &(epsilonSamplers3[iScale].getChildren(id_eps));
                break;
              case 4:
                pchildren = &(epsilonSamplers4[iScale].getChildren(id_eps));
                break;
              case 5:
                pchildren = &(epsilonSamplers5[iScale].getChildren(id_eps));
                break;
              case 10:
                pchildren = &(epsilonSamplers10[iScale].getChildren(id_eps));
                break;
              default:
                printErr("Invalid FitDimension");
                break;
            }

            const auto& children = *pchildren;
            for(const auto& id_filtered_child : children){
              SimplexId id_local_child = filteredToLocalIdVec[iScale][id_filtered_child];
              SimplexId id_global_child = sampling_indexes[iScale][id_local_child];
              globalIds.push_back(id_global_child);
            }
          }

          scalePairs[iCC] = {minScale, maxScale};
          float scaleRange = std::log(scales[maxScale])-std::log(scales[minScale]);
          double productPers = (double)scaleRange*persistence;
          if(productPers>maxProductPersistence){
            maxProductPersistence = productPers;
          }
          productPersistence[iCC] = productPers;
          if(UseProductPersistenceScore){
            persistence = productPers;
          }

          if(maxScale<scaleIdThresholdMin 
              or minScale>scaleIdThresholdMax 
              or (ScaleRangeMin > 0 and scaleRange < ScaleRangeMin)){
            continue;
          }

          TTK_PSORT(this->threadNumber_, globalIds.begin(), globalIds.end());
          auto last = std::unique(globalIds.begin(), globalIds.end());
          globalIds.erase(last, globalIds.end());

          sizesAtThreshold[iCC] = globalIds.size();
          if(globalIds.size()==0){
            printErr("EMPTY COMPONENT, SHOULD NOT HAPPEN !!! SEGFAULT IN 3... 2... 1...");
          }

          for(auto id_global : globalIds){

            if(highestPersComponent[id_global]==-1){
              highestPersComponent[id_global]=iCC;
            }else{

              const double current_highest_ccId = highestPersComponent[id_global];
              const auto& current_highest_cc = treeMS.getComponent(current_highest_ccId);
              double current_highest_pers = current_highest_cc.death()-current_highest_cc.birth();

              if(current_highest_pers<persistence){
                highestPersComponent[id_global] = iCC;
              }
            }

          }
        } // end for iCC
          //
          //


        // compute score for each component

        for(int iPoint=0; iPoint< nbPoints; iPoint++){
          if(highestPersComponent[iPoint]!=-1){
            SimplexId ccId = highestPersComponent[iPoint];
            vote[ccId]++;
          }
        }

        for(int iCC=0; iCC<nbOutputPrimitives; iCC++){
          scores[iCC] = (float)vote[iCC]/((float)sizesAtThreshold[iCC]);
        }


      } // end block compute score


      pseg.setNumberOfPoints(nbPoints);
      pseg.setTolerance(GrowPrimitiveThreshold*sqrt(diag_aabb));


      // Filter pairs
      std::vector<SimplexId> chosenPairs{};
      std::vector<SimplexId> chosenFlag(nbOutputPrimitives, -1);
      const std::vector<float>& scales = sampler_fitting.scales();

      for(SimplexId iCC=0; iCC < nbOutputPrimitives; iCC++){
          auto& cc = treeMS.getComponent(iCC);

          double birth = cc.birth();
          double death = cc.death();
          double persistence = death-birth;
          float score = scores[iCC];
          int minScale = scalePairs[iCC].first;
          int maxScale = scalePairs[iCC].second;
          float scaleRange = std::log(scales[maxScale])-std::log(scales[minScale]);
          double productPers = productPersistence[iCC];

          // std::cout<<"looking at pair "<<iCC<<" score "<< score<<" scales "<<minScale<<"-"<<maxScale<<std::endl;
          if(persistence < PersistenceMin){
            chosenFlag[iCC] = -1;
            continue;
          }
          if(birth > BirthMax){
            chosenFlag[iCC] = -1;
            continue;
          }
          if(productPers < maxProductPersistence*ProductPersistenceMin){
            chosenFlag[iCC] = -2;
            continue;
          }
          if( (minScale > -1 and minScale > scaleIdThresholdMax)
           or (maxScale > -1 and maxScale < scaleIdThresholdMin) 
           or (ScaleRangeMin > 0 and scaleRange < ScaleRangeMin)){
            chosenFlag[iCC] = -3;
            continue;
          }
          if(score < ScoreMin){
            chosenFlag[iCC] = -4;
            continue;
          }

          chosenFlag[iCC] = chosenPairs.size();
          chosenPairs.push_back(iCC);

          // decrease parent threshold
          SimplexId iCCParent = cc.parent();
          // thresholds[iCCParent] = std::min(death, thresholds[iCCParent]);

          // set display diagmax
        auto& diags = cc.diag_aabb_5D();

        if(DiagMaxDisplay>=0){ 
          bool found_good_diag=false;

          for(size_t iDiag=0; iDiag<diags.size(); iDiag++) {
            auto diag = diags[iDiag];
            if(diag<DiagMaxDisplay){
              thresholds[iCC] = std::min(treeMS_.getComponent(cc.getChild(iDiag)).death(), thresholds[iCC]);
              found_good_diag=true;
              break;
            }
          }
          if(!found_good_diag){ // set threshold to birth point (diag = 0), as the introduction of the first child already induces diag>diagmax
            thresholds[iCC] = birth;
          }
        }

      }


      { // add score and chosen id to pd

        vtkNew<vtkFloatArray> pdScoreArray{};
        pdScoreArray->SetNumberOfComponents(1);
        pdScoreArray->SetNumberOfTuples(pd->GetNumberOfCells());
        pdScoreArray->SetName("score");
        pdScoreArray->Fill(-1);
        pd->GetCellData()->AddArray(pdScoreArray);

        vtkNew<vtkIntArray> chosenIdArray{};
        chosenIdArray->SetNumberOfComponents(1);
        chosenIdArray->SetNumberOfTuples(pd->GetNumberOfCells());
        chosenIdArray->SetName("chosenId");
        chosenIdArray->Fill(-1);
        pd->GetCellData()->AddArray(chosenIdArray);

        vtkNew<vtkFloatArray> scaleArray{};
        scaleArray->SetNumberOfComponents(2);
        scaleArray->SetNumberOfTuples(pd->GetNumberOfCells());
        scaleArray->SetName("scaleRange");
        scaleArray->Fill(-1);
        pd->GetCellData()->AddArray(scaleArray);


        // std::vector<SimplexId> ccIdToChosenId(nbOutputPrimitives, -1);
        // for(auto iChosen=0; iChosen<chosenPairs.size(); iChosen++){
        //   ccIdToChosenId[chosenPairs[iChosen]] = iChosen;
        // }

        const std::vector<float>& scales = sampler_fitting.scales();
        for(int iCell=0; iCell<pd->GetNumberOfCells(); iCell++){
          SimplexId iPair = pd->GetCellData()->GetArray("PairIdentifier")->GetTuple1(iCell);
          if(iPair>-1 and iPair < nbOutputPrimitives){
            SimplexId ccId = iPair;
            
            float score = scores[ccId];
            SimplexId iChosen = chosenFlag[ccId];
            float maxScale = scalePairs[ccId].second >= 0 ? scales[scalePairs[ccId].second] : -1;
            float minScale = scalePairs[ccId].first >= 0 ? scales[scalePairs[ccId].first] : -1;

            pdScoreArray->SetTuple1(iCell, score);
            chosenIdArray->SetTuple1(iCell, iChosen);
            scaleArray->SetTuple2(iCell, minScale, maxScale);
          }
        }
      } // end block add to pd

      // create output product diagram
      {
        vtkNew<vtkUnstructuredGrid> vtu;
        ProductDiagramToVTU(vtu,
                            diagramMS, 
                            chosenFlag,
                            chosenPairs,
                            scalePairs,
                            sampler_fitting.scales(),
                            scores
                            );
        outputPDMB->SetBlock(1, vtu);

      }

      const auto nbChosenComponents = chosenPairs.size();
      outputMB->SetNumberOfBlocks(nbChosenComponents);
      for(SimplexId iChosen=0; iChosen<nbChosenComponents; iChosen++){

        const auto ccId = chosenPairs[iChosen];

        std::cout<<" component "<<ccId<<" "<<iChosen<<" "<<chosenFlag[ccId]<<" "<<chosenFlag[iChosen]<<std::endl;

        const auto& cc = treeMS.getComponent(ccId);
        auto threshold = thresholds[ccId];
        auto birth = cc.birth();
        float score = scores[ccId];
        // auto persistence = ripsArcs[iArc].death-ripsArcs[iArc].birth;
        // auto type = ripsArcs[iArc].type;
        // auto diag = ripsArcs[iArc].diag;


        vtkNew<vtkUnstructuredGrid> block;

        bool discard = false;

        std::vector<float> primitiveScales = sampler_fitting.scales();
        // for(auto& s : primitiveScales){
        //   s/=sqrt(diag_aabb);
        // }

        switch(FitDimension){
          case 1:
            discard = createOutputConnectedComponentMultiScale(block, treeMS, ccId, threshold, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, embeddingCoordinates, isStableVec, densityVec, 1, highestPersComponent, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers1, epsilonIdToFilteredId);
            break;
          case 2:
            discard = createOutputConnectedComponentMultiScale(block, treeMS, ccId, threshold, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, embeddingCoordinates, isStableVec, densityVec, 1, highestPersComponent, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers2, epsilonIdToFilteredId);
            break;
          case 3:
            discard = createOutputConnectedComponentMultiScale(block, treeMS, ccId, threshold, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, embeddingCoordinates, isStableVec, densityVec, 1, highestPersComponent, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers3, epsilonIdToFilteredId);
            break;
          case 4:
            discard = createOutputConnectedComponentMultiScale(block, treeMS, ccId, threshold, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, embeddingCoordinates, isStableVec, densityVec, 1, highestPersComponent, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers4, epsilonIdToFilteredId);
            getPrimitives(pseg, treeMS, ccId, Primitive::PrimitiveType::PLANE, threshold, GrowPrimitiveThreshold*diag_aabb, coordinates, normals, fitCoordsVec, primitiveScales, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers4, epsilonIdToFilteredId);
            break;
          case 5:
            discard = createOutputConnectedComponentMultiScale(block, treeMS, ccId, threshold, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, embeddingCoordinates, isStableVec, densityVec, 1, highestPersComponent, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers5, epsilonIdToFilteredId);
            getPrimitives(pseg, treeMS, ccId, Primitive::PrimitiveType::SPHERE, threshold, GrowPrimitiveThreshold*diag_aabb, coordinates, normals, fitCoordsVec, primitiveScales, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers5, epsilonIdToFilteredId);
            break;
          case 10:
            discard = createOutputConnectedComponentMultiScale(block, treeMS, ccId, threshold, coordinates, normals, fitCoordsVec, mlsCoordsVec, mlsNormalsVec, embeddingCoordinates, isStableVec, densityVec, 1, highestPersComponent, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers10, epsilonIdToFilteredId);
            getPrimitives(pseg, treeMS, ccId, Primitive::PrimitiveType::QUADRIC, threshold, GrowPrimitiveThreshold*diag_aabb, coordinates, normals, fitCoordsVec, primitiveScales, resultBarycenter, sampling_indexes, filteredToLocalIdVec, pointsAtScaleStartIdx, epsilonSamplers10, epsilonIdToFilteredId);
            break;
          default:
            printErr("Invalid FitDimension");
            break;
        }

        SimplexId size = block->GetNumberOfPoints();
        int dummy = -1;
        // setSingleComponentAttributes(block, ccId, persistence, arcCount, type, diag,-1, -1);


        vtkNew<vtkFloatArray> scoreArray;
        scoreArray->SetNumberOfComponents(1);
        scoreArray->SetNumberOfTuples(block->GetNumberOfPoints());
        scoreArray->SetName("score");
        scoreArray->Fill(score);
        block->GetPointData()->AddArray(scoreArray);

        vtkNew<vtkIntArray> chosenIdArray;
        chosenIdArray->SetNumberOfComponents(1);
        chosenIdArray->SetNumberOfTuples(block->GetNumberOfPoints());
        chosenIdArray->SetName("chosenId");
        chosenIdArray->Fill(iChosen);
        block->GetPointData()->AddArray(chosenIdArray);

        outputMB->SetBlock(iChosen, block);
      }

      timings.clock(STEP::PROCESS_ARCS);



    }

    // SEGMENTATION OUTPUT
    {
    vtkMultiBlockDataSet *outputSegMB = vtkMultiBlockDataSet::GetData(outputVector,5);

    timings.clock(STEP::GROW_COVERAGE);
    if(GrowPrimitives){
      for(int it=0; it<GrowPrimitiveIterations; it++){
        for(int iPrim=0; iPrim<pseg.getNumberOfPrimitives(); iPrim++){
          pseg.growPrimitive(iPrim, coordinates, normals, sampler_subsampling.kdTree(subsampling_scale_index));
        }
      }
    }
    timings.clock(STEP::GROW_COVERAGE);
    createOutputPrimitiveSegmentation(pseg, outputSegMB, coordinates, normals, resultBarycenter);
    }

    // implicit Primitive
    {
      vtkMultiBlockDataSet *outputSegMB = vtkMultiBlockDataSet::GetData(outputVector,5);
      vtkMultiBlockDataSet *outputVTIMB = vtkMultiBlockDataSet::GetData(outputVector,7);

      size_t nBlocks = outputSegMB->GetNumberOfBlocks();
      outputVTIMB->SetNumberOfBlocks(nBlocks);

      DimensionReduction dr{};
      auto method = DimensionReduction::METHOD::PCA;
      dr.setIsInputDistanceMatrix(false);
      dr.setInputMethod(method);
      dr.setDebugLevel(debugLevel_);
      dr.setInputModulePath("/home/julius/ttk/core/base/dimensionReduction");
      dr.setInputNumberOfComponents(3);

      for(int iBlock = 0; iBlock<nBlocks; iBlock++){
        // std::cout<<" iblock "<<iBlock<<std::endl;
        vtkNew<vtkImageData> outputVTI;
        outputVTIMB->SetBlock(iBlock, outputVTI);
        vtkUnstructuredGrid* primVtu = static_cast<vtkUnstructuredGrid*>(outputSegMB->GetBlock(iBlock));

        double primBounds[6];
        int Resolution=20;
        primVtu->GetBounds(primBounds);
        float aabb = float(primBounds[1]-primBounds[0])*float(primBounds[1]-primBounds[0]) + float(primBounds[3]-primBounds[2])*float(primBounds[3]-primBounds[2]) + float(primBounds[5]-primBounds[4])*float(primBounds[5]-primBounds[4]);
        int dimX=Resolution, dimY=Resolution, dimZ=Resolution;
        float offsetBounds[6];
        float offsetSpacing=Spacing*(0.05*sqrt(aabb));
        offsetBounds[0] = primBounds[0]- offsetSpacing;
        offsetBounds[1] = primBounds[1]+ offsetSpacing;
        offsetBounds[2] = primBounds[2]- offsetSpacing;
        offsetBounds[3] = primBounds[3]+ offsetSpacing;
        offsetBounds[4] = primBounds[4]- offsetSpacing;
        offsetBounds[5] = primBounds[5]+ offsetSpacing;
        float ox=offsetBounds[0],oy=offsetBounds[2], oz=offsetBounds[4];
        float dx = (offsetBounds[1]-offsetBounds[0])/(dimX-1);
        float dy = (offsetBounds[3]-offsetBounds[2])/(dimY-1);
        float dz = dimZ > 1 ? (offsetBounds[5]-offsetBounds[4])/(dimZ-1) : 0;
        std::vector<float> dxyz = {dx,dy,dz};

        outputVTI->SetDimensions(dimX,dimY,dimZ);
        outputVTI->SetOrigin(offsetBounds[0], offsetBounds[2], offsetBounds[4]);
        outputVTI->SetSpacing(dx,dy,dz);

        vtkNew<vtkFloatArray> fitScalarField{};
        fitScalarField->SetNumberOfComponents(1);
        fitScalarField->SetNumberOfTuples(dimX*dimY*dimZ);
        fitScalarField->SetName("f");
        fitScalarField->Fill(-1);
        float* fitSFPtr = (float*)ttkUtils::GetVoidPointer(fitScalarField);

        outputVTI->GetPointData()->AddArray(fitScalarField);

        const auto& prim = pseg.getPrimitive(iBlock);

        for(SimplexId I=0; I<dimX; I++){
          for(SimplexId J=0; J<dimY; J++){
            for(SimplexId K=0; K<dimZ; K++){
              SimplexId ID = I+J*dimX+K*dimX*dimY;
              float p[3] = {ox + I*dx, oy + J*dy, oz + K*dz};
              // double x = ox + I*dx;
              // double y = oy + J*dy;
              // double z = oz + K*dz;

              double fx = prim.f(p);
              fitSFPtr[ID]= fx;
            }
          }
        }



        SimplexId nbPrimPoints = prim.getNumberOfPoints();
        std::vector<double> mat(nbPrimPoints*3, 0);

        // std::cout<<" values ";
        int dim = 3;
        std::vector<float> covVector(dim*dim, 0);
        std::vector<float> mean(dim, 0);
        for(int ii=0; ii<nbPrimPoints; ii++){
          SimplexId iPoint = prim.point_indices()[ii];
          for(int iDim=0; iDim<3; iDim++){
            mean[iDim]  += coordinates[3*iPoint+iDim]/nbPrimPoints;
          }
        }
        Primitive::VectorType meanVec;
        meanVec << mean[0], mean[1], mean[2];
        for(int ii=0; ii<nbPrimPoints; ii++){

          SimplexId iPoint = prim.point_indices()[ii];


          // float p[3] = { coordinates[3*iPoint+0] - ox, coordinates[3*iPoint+1]-oy, coordinates[3*iPoint+2]-oz};

          for(int i=0; i<dim; i++){
            for(int j=i; j<dim; j++){
              covVector[i*dim+j] += (coordinates[dim*iPoint+i]-mean[i])*(coordinates[dim*iPoint+j]-mean[j])/nbPrimPoints;
              covVector[j*dim+i] = covVector[i*dim+j];
            }
          }
        }

        for(int i=0; i<dim; i++){
          for(int j=i; j<dim; j++){
            covVector[j*dim+i] = covVector[i*dim+j];
          }
        }

        Eigen::Matrix<float, -1, -1, 0, 6, 6> covMatrix = Eigen::Map<Eigen::Matrix<float, -1, -1>, Eigen::RowMajor>(covVector.data(), dim, dim);
        Eigen::EigenSolver<Eigen::Matrix<float, -1, -1, 0, 6, 6>> solver(covMatrix);
        Primitive::VectorType ev = covMatrix.eigenvalues().real();
        Eigen::Matrix<float, -1, -1, 0, 6, 6> P = solver.eigenvectors().real(); // std::cout<<"eigen values "<<ev<<std::endl;
        // std::cout<<"eigen values "<<ev<<std::endl;
        // std::cout<<"eigen vectors "<<P<<std::endl;

        // find bounding ids
        std::vector<SimplexId> boundIdMin(3,-1);
        std::vector<SimplexId> boundIdMax(3,-1);
        std::vector<float> boundsMin(3,-1);
        std::vector<float> boundsMax(3,-1);

        float evs[3] ={ev.coeff(0), ev.coeff(1), ev.coeff(2)};

        if(evs[0]<evs[1]){
          P.col(0).swap(P.col(1));
          std::swap(evs[0], evs[1]);
        }
        if(evs[1]<evs[2]){
          P.col(1).swap(P.col(2));
          std::swap(evs[1], evs[2]);
        }
        if(evs[0]<evs[1]){
          P.col(0).swap(P.col(1));
          std::swap(evs[0], evs[1]);
        }
        // std::cout<<"sorted"<<std::endl;
        // std::cout<<"eigen values "<<ev<<std::endl;
        // std::cout<<"eigen vectors "<<P<<std::endl;

        //     // std::cout<<" dr size  "<<dr_output.size()<<std::endl;
        //     // std::cout<<" dr size "<<dr_output[0].size()<<std::endl;

        for(int ii=0; ii<nbPrimPoints; ii++){
          SimplexId iPoint = prim.point_indices()[ii];
          // std::cout<<" point "<<ii<<" "<<iPoint<<": ";
          Primitive::VectorType x = {coordinates[3*iPoint+0], coordinates[3*iPoint+1], coordinates[3*iPoint+2]};
          x = x-meanVec;
          auto y = P.transpose()*x;
          // std::cout<<" y "<<y.coeff(0)<<" "<<y.coeff(1)<<" "<<y.coeff(2)<<std::endl;

          for(int iDim=0; iDim<3; iDim++){
            if(boundIdMin[iDim]==-1 or y.coeff(iDim)<boundsMin[iDim]){
              boundIdMin[iDim] = iPoint;
              boundsMin[iDim] = y.coeff(iDim);
            }
            if(boundIdMax[iDim]==-1 or y.coeff(iDim)>boundsMax[iDim]){
              boundIdMax[iDim] = iPoint;
              boundsMax[iDim] = y.coeff(iDim);
            }
          }
        }
        //     std::cout<<"here"<<std::endl;

        //     // std::cout<<" boundId "<<boundId[0]<<std::endl;
        //     // std::cout<<" boundId "<<boundId[1]<<std::endl;
        //     // std::cout<<" boundId "<<boundId[2]<<std::endl;
        //     // std::cout<<" boundId "<<boundId[3]<<std::endl;
        Primitive::VectorType p0 = {boundsMin[0], boundsMin[1], boundsMin[2]};
        Primitive::VectorType p1 = {boundsMax[0], boundsMin[1], boundsMin[2]};
        Primitive::VectorType p2 = {boundsMax[0], boundsMax[1], boundsMin[2]};
        Primitive::VectorType p3 = {boundsMin[0], boundsMax[1], boundsMin[2]};

        // std::cout<<" points: "<<std::endl;
        // std::cout<<p0<<" "<<p1<<" "<<p2<<" "<<p3<<std::endl;


        p0 = P*p0+meanVec;
        p1 = P*p1+meanVec;
        p2 = P*p2+meanVec;
        p3 = P*p3+meanVec;

        Primitive::VectorType n0, n1, n2, n3, q0, q1, q2, q3;

        prim.project(p0, q0, n0);
        prim.project(p1, q1, n1);
        prim.project(p2, q2, n2);
        prim.project(p3, q3, n3);

        //     // std::cout<<" points "<<p0[0]<<","<<p0[1]<<","<<p0[2]<<std::endl;
        //     // std::cout<<" points "<<p1[0]<<","<<p1[1]<<","<<p1[2]<<std::endl;
        //     // std::cout<<" points "<<p2[0]<<","<<p2[1]<<","<<p2[2]<<std::endl;
        //     // std::cout<<" points "<<p3[0]<<","<<p3[1]<<","<<p3[2]<<std::endl;

        vtkNew<vtkUnstructuredGrid> vtu;
        vtkNew<vtkPoints> points;
        points->InsertNextPoint(q0.data());
        points->InsertNextPoint(q1.data());
        points->InsertNextPoint(q2.data());
        points->InsertNextPoint(q3.data());

        // vtu->Allocate();

        vtkNew<vtkIdList> list;
        list->InsertNextId(0);
        list->InsertNextId(1);
        list->InsertNextId(2);
        list->InsertNextId(3);

        // std::cout<<"cell"<<std::endl;
        vtu->SetPoints(points);
        vtu->InsertNextCell(VTK_QUAD, list);
        // outputVTIMB->SetBlock(iBlock, vtu);








      // try to remove parts of the volume
        // std::vector<int> here(dimX*dimY*dimZ,-1);
        // for(int ii=0; ii<nbPrimPoints; ii++){

        //   SimplexId iPoint = prim.point_indices()[ii];


        //   std::vector<float> p(3);
        //   p[0] = coordinates[3*iPoint+0]-ox;
        //   p[1] = coordinates[3*iPoint+1]-oy;
        //   p[2] = coordinates[3*iPoint+2]-oz;


        //   std::vector<int> ID(3,0);
        //   for(int iDim=0; iDim<3; iDim++){
        //     while(p[iDim]>0){
        //       p[iDim] -= dxyz[iDim];
        //       ID[iDim]++;
        //       std::cout<<"ID["<<iDim<<"]"<<ID[iDim]<<std::endl;
        //     }
        //   }
        //   SimplexId I = ID[0]+ID[1]*dimX+ID[2]*dimX*dimY;
        //   std::cout<<"I "<<I<<" outta "<<here.size()<<std::endl;
        //   here[I]++;
        // }
        // for(SimplexId I=0; I<dimX; I++){
        //   for(SimplexId J=0; J<dimY; J++){
        //     for(SimplexId K=0; K<dimZ; K++){
        //       SimplexId ID = I+J*dimX+K*dimX*dimY;
        //       if(here[ID]==-1){
        //         fitSFPtr[ID] = NAN;
        //       }
        //     }
        //   }
        // }



      }

    } // end implicit primitive

    // MERGE TREE OUTPUT
    {



//       std::cout<<"GO FOR MERGE TREE"<<std::endl;

//       vtkMultiBlockDataSet *outputScaleMTMB = vtkMultiBlockDataSet::GetData(outputVector,3);
//       const auto& ripsArcs = treeMS.getRipsArcs();

//       std::vector<std::vector<RipsArc>> arcs{};
//       std::vector<std::vector<SimplexId>> arcIds{};
//       std::vector<std::vector<SimplexId>> leafIds{};
//       SimplexId leafCount=0;
//       for(size_t iArc = 0; iArc<ripsArcs.size(); iArc++){
//         int leafId = -1;
//         auto& arc = ripsArcs[iArc];
//         if(arc.type<1){
//           leafId = leafCount;
//           leafCount++;
//         }
//         SimplexId treeId = arc.treeId;
//         if(treeId>=arcs.size()){
//           arcs.resize(treeId+1);
//           arcIds.resize(treeId+1);
//           leafIds.resize(treeId+1);
//         }
//         arcs[treeId].push_back(arc);
//         arcIds[treeId].push_back(iArc);
//         leafIds[treeId].push_back(leafId);
//       }

//       SimplexId nbTrees = arcs.size();

//       outputScaleMTMB->SetNumberOfBlocks(nbTrees);

//       for(int treeId = 0; treeId<nbTrees; treeId++){
//         vtkNew<vtkMultiBlockDataSet> outputMTMB;
//         outputScaleMTMB->SetBlock(treeId, outputMTMB);
//         // std::cout<<" \n\n MT "<<treeId<<std::endl;
//         ripsArcsToMergeTree(outputMTMB, arcs[treeId], arcIds[treeId], leafIds[treeId]);
//         // std::cout<<" done MT "<<treeId<<std::endl;
//       }
//       std::cout<<" done all MTs"<<std::endl;



    }


    //{
//if(DisplayFit and !selection.empty()){
    //          // vtkNew<vtkIntArray> neighborhoodArray{};
    //          // neighborhoodArray->SetNumberOfComponents(1);
    //          // neighborhoodArray->SetNumberOfTuples(nbPoints);
    //          // neighborhoodArray->SetName("fitNeighborhood");
    //          // neighborhoodArray->Fill(0);
    //          // int* neighborhoodPtr = static_cast<int*>(ttkUtils::GetVoidPointer(neighborhoodArray));
    //          // outputBlock->GetPointData()->AddArray(neighborhoodArray);

    //          ScaleSpaceSampling<pointCloud::MyPoint::Scalar, pointCloud::MyPoint>& sampler3D = sampler_fitting_;
    //          double scale = sampler3D.scaleiRatio);
    //          auto& tree = sampler3D.kdTree(iRatio);

    //          for(auto pointId : selection){
    //            std::cout<<"Display fit for "<<pointId<<std::endl;

    //              // std::vector<SimplexId> neighborHood{};
    //              // if(sampler3D.isUseKnnGraph()){
    //              //   pointCloud::getNeighborhoodKnnGraph<ttk::pointCloud::MyPoint>(pointId,scale, sampler3D.knnGraph(), neighborHood);
    //              // }else{
    //              //   pointCloud::getNeighborhood<ttk::pointCloud::MyPoint>(pointId,scale,coordinates, tree, neighborHood);
    //              // }

    //              // for(auto& nId : neighborHood){
    //              //   neighborhoodPtr[nId] = 1;
    //              // }

    //              //retrieve coefficients of algebraic sphere
    //              auto uc = fitCoordinates[5*pointId];
    //              auto uq = fitCoordinates[5*pointId+1];
    //              auto ulx = fitCoordinates[5*pointId+2];
    //              auto uly = fitCoordinates[5*pointId+3];
    //              auto ulz = fitCoordinates[5*pointId+4];

    //              for(SimplexId I=0; I<dimX; I++){
    //                for(SimplexId J=0; J<dimY; J++){
    //                  for(SimplexId K=0; K<dimZ; K++){
    //                    SimplexId ID = I+J*dimX+K*dimX*dimY;
    //                    double x = ox + I*dx - resultBarycenter[0];
    //                    double y = oy + J*dy - resultBarycenter[1];
    //                    double z = oz + K*dz - resultBarycenter[2];
    //                    double fx = uc + ulx*x+uly*y+ulz*z + uq*(x*x+y*y+z*z);
    //                    fitSFPtr[ID]= fx;
    //                  }
    //                }
    //              }
    //          }
    //          // for(auto pointId : selection){
    //          //   neighborhoodPtr[pointId] = 2;
    //          // }
    //        }
 
    //}

  }


  // Apply selection
  {
    // if(!selection.empty()){
    //   printMsg("Apply selection");
    // }
    // for(int iRatio = 0; iRatio< outputMB->GetNumberOfBlocks(); iRatio++){
    //   vtkMultiBlockDataSet* blockRatio = (vtkMultiBlockDataSet*)outputMB->GetBlock(iRatio);
    //   for(int iPers = 0; iPers< blockRatio->GetNumberOfBlocks(); iPers++){
    //     vtkUnstructuredGrid* blockPers = (vtkUnstructuredGrid*)blockRatio->GetBlock(iPers);
    //     if(!selection.empty()){
    //       int nbConnectedComponents =
    //       blockPers->GetFieldData()->GetArray("NumberOfComponents")->GetTuple1(0);
    //       applySelection(nbConnectedComponents, blockPers, selection);
    //     }else{
    //       resetSelection(blockPers);
    //     }
    //   }
    // }
  }



  timings.clock(STEP::TOTAL);
  if(WriteTimings){
    std::stringstream filename;
    filename << TimingFilePrefix << "_" << fittingStringId_ << ".txt";
    timings.save(filename.str());
  }

  printMsg(debug::output::INVERTED+debug::output::PINK+"the end :)"+debug::output::ENDCOLOR+debug::output::ENDCOLOR, 1, tm_total.getElapsedTime(), threadNumber_);
  return 1;
  }

  void ttkRipsMultiParameterSampling::resetSelection(vtkUnstructuredGrid* filteredTree){
    if(filteredTree->GetNumberOfPoints()==0){
      return;
    }
    auto isSelectedArray = filteredTree->GetPointData()->GetArray("isSelected");
    if(isSelectedArray==nullptr){
      printErr("No selection array to fill");
      return;
    }
    for(SimplexId i=0; i<filteredTree->GetNumberOfPoints(); i++){
      isSelectedArray->SetTuple1(i, 1);
    }
  }

  void ttkRipsMultiParameterSampling::applySelection(SimplexId nbConnectedComponents, vtkUnstructuredGrid* filteredTree, std::vector<SimplexId> selection){

    // empty selection : copy all connected components
    if(selection.empty()){
      return;
    }


    auto totalNbPair = filteredTree->GetNumberOfPoints();

    std::vector<SimplexId> pairId(totalNbPair);
    auto pd = filteredTree->GetPointData();
    auto currentParentArray = pd->GetArray("componentParentAtThreshold");

    std::vector<int> isSelected(nbConnectedComponents,0);

    for(SimplexId i=0; i<totalNbPair; i++){
      auto birthId = pd->GetArray("BirthId")->GetTuple1(i);
      pairId[birthId] = i;
    }


    for(SimplexId v : selection){
      isSelected[currentParentArray->GetTuple1(pairId[v])] = 1;
    }

    auto isSelectedArray = pd->GetArray("isSelected");

    for(SimplexId i=0; i<totalNbPair; i++){
      SimplexId currentParent = currentParentArray->GetTuple1(i);
      isSelectedArray->SetTuple1(i, isSelected[currentParent]);
    }

  }

  void ttkRipsMultiParameterSampling::readSelection(std::vector<SimplexId>& selectionVec, vtkPointSet* inputVTKSelection){

    if(inputVTKSelection == nullptr){
      return;
    }

    auto pd = inputVTKSelection->GetPointData();
    SimplexId nbPoints = inputVTKSelection->GetNumberOfPoints();
    if(nbPoints==0){
      printWrn("Selection contains no points -- skipping");
      return;
    }

    auto idArray = pd->GetArray("PointIds");
    auto birthIdArray = pd->GetArray("BirthId");

    vtkDataArray* selectionArray = birthIdArray == nullptr ? idArray : birthIdArray;

    if(selectionArray == nullptr){
      printErr("Selection must contain a valid id array -- skipping");
      return;
    }

    for(SimplexId i=0; i<inputVTKSelection->GetNumberOfPoints();i++){
      selectionVec.push_back(selectionArray->GetTuple1(i));
    }

    std::sort(selectionVec.begin(), selectionVec.end());

  }

  void ttkRipsMultiParameterSampling::translate(vtkUnstructuredGrid* vtu, int
      I, int J, int K, float spacing, float* coordinates){

    if(vtu->GetNumberOfPoints()==0){
      return;
    }
    float* coords = (float*) ttkUtils::GetVoidPointer(vtu->GetPoints()->GetData());

    // #ifdef TTK_ENABLE_OPENMP
    // #pragma omp parallel for num_threads(threadNumber_)
    // #endif // TTK_ENABLE_OPENMP
    for(SimplexId i =0; i<vtu->GetNumberOfPoints(); i++){
      // SimplexId birthId = vtu->GetPointData()->GetArray("BirthId")->GetTuple1(i);
      coords[3*i+0] = I*spacing + coordinates[3*i+0];
      coords[3*i+1] = J*spacing + coordinates[3*i+1];
      coords[3*i+2] = K*spacing + coordinates[3*i+2];
    }
  }

  void ttkRipsMultiParameterSampling::fillSparseGridPersistences(const std::vector<std::vector<double>>& persBigVec,std::vector<std::vector<double>>& sparseGrid, int nbRows){
    sparseGrid.clear();

    std::vector<double> allPersistences{};
    for(const auto& pv : persBigVec){
      for(const auto& p : pv){
        allPersistences.push_back(p);
      }
    }
    std::sort(allPersistences.begin(), allPersistences.end());

    for(int i=0; i<nbRows; i++ ){
      sparseGrid.push_back(allPersistences);
    }

    for(int iRow=0; iRow<nbRows; iRow++){
      auto iP = 0;
      for(int iCol=0; iCol<allPersistences.size(); iCol++){
        if(iP<persBigVec[iRow].size() && sparseGrid[iRow][iCol]==persBigVec[iRow][iP]){
          iP++;
        }else{
          sparseGrid[iRow][iCol]=-1;
        }
      }
    }
  }

  void ttkRipsMultiParameterSampling::fillFullGridPersistences(const std::vector<std::vector<double>>& persBigVec,std::vector<std::vector<double>>& grid, int nbRows){
    grid.clear();

    std::vector<double> allPersistences{};
    for(const auto& pv : persBigVec){
      for(const auto& p : pv){
        allPersistences.push_back(p);
      }
    }
    std::sort(allPersistences.begin(), allPersistences.end());

    for(int i=0; i<nbRows; i++ ){
      grid.push_back(allPersistences);
    }
  }

  void ttkRipsMultiParameterSampling::computeStatistics(const RipsTree& tree, int iRatio, int iThreshold, const std::vector<std::vector<double>>& persistenceGrid, float* normals, float* meanNormals, int* parentSize, int* fullsize, int* parent, int* parent2d, int* fullId, double* pers1d, double* pers2d, int* pers2dId, int* localIdPtr, int* primitiveTypePtr, int* displayPtr){

    double threshold = persistenceGrid[iRatio][iThreshold];
    SimplexId nCol = persistenceGrid[iRatio].size();

    if(tree.maxDeath()<threshold){
      threshold=tree.maxDeath();
    }
    //process connected components until the threshold is met
    SimplexId localCCId=0;
    for(int ccId=0; ccId<tree.getNumberOfPoints(); ccId++){
      if(tree.getPersistence(ccId)<threshold or localCCId>=NbCCMax){
        break;
      }
      // std::vector<SimplexId> pointIds;
      tree.computeMeanFieldAtThreshold(ccId, threshold, meanNormals, 3);
      SimplexId size = tree.augmentedSizeAtThreshold(ccId, threshold, SizeMin);
      SimplexId size2 = tree.sizeAtThreshold(ccId, threshold);
      SimplexId size3 = tree.fullSize(ccId);
      tree.dispatchScalarComponentAtThreshold(ccId, threshold, parentSize, (int)size, SizeMin);
      // tree.dispatchScalarComponentAtThreshold(ccId, threshold, fullsize, (int)tree.getSize(ccId), SizeMin);
      if(size>=SizeMin){
        tree.dispatchScalarComponentAtThreshold(ccId, threshold, parent, (int)ccId, SizeMin);
        // if(tree.sizeAtThreshold(ccId, threshold)>=SizeMin){
        // if(ccId==0){
        tree.dispatchScalarComponentAtThreshold(ccId, threshold, displayPtr, 1, 1);
        // }
        // else{
        //   tree.dispatchScalarComponentAtThreshold(ccId, threshold, displayPtr, 1, SizeMin);
        // }
        // }
        SimplexId fullGridId = getFullGridId(localCCId, iRatio, iThreshold, nCol, NbCCMax);
        tree.dispatchScalarComponentAtThreshold(ccId, threshold, fullId, fullGridId, SizeMin);
        // tree.dispatchScalarComponentAtThreshold(ccId, threshold,  pers1dId, parentList_[parents2D_[fullGridId]]);
        if(RatioVersion){

          SimplexId primitiveType=0;
          tree.computeCovarianceMatrixAtThreshold(ccId, threshold, normals, meanNormals, 3, primitiveType);
          tree.dispatchScalarComponentAtThreshold(ccId, threshold, primitiveTypePtr, primitiveType, SizeMin);
          tree.dispatchScalarComponentAtThreshold(ccId, threshold, parent2d, parents2D_[fullGridId], SizeMin);
          tree.dispatchScalarComponentAtThreshold(ccId, threshold,  pers2dId, parentList_[parents2D_[fullGridId]], SizeMin);
          tree.dispatchScalarComponentAtThreshold(ccId, threshold,  pers2d, persistence2d_[parents2D_[fullGridId]], SizeMin);
          tree.dispatchScalarComponentAtThreshold(ccId, threshold,  localIdPtr, localCCId, SizeMin);
        }
        localCCId++;
      }
    }
  }

  void ttkRipsMultiParameterSampling::setSingleComponentAttributes(vtkUnstructuredGrid* vtu, SimplexId fullId, double persistence2D, SimplexId persistence2DId, int criticalType, double diag, double scale, double autosimilarity){

    SimplexId nbPoints = vtu->GetNumberOfPoints();
    float* normalsPtr = (float*)ttkUtils::GetVoidPointer(vtu->GetPointData()->GetArray("Normals"));

    vtkNew<vtkIntArray> criticalTypeArray;
    criticalTypeArray->SetName("criticalType");
    criticalTypeArray->SetNumberOfComponents(1);
    criticalTypeArray->SetNumberOfTuples(nbPoints);
    criticalTypeArray->Fill(criticalType);
    vtu->GetPointData()->AddArray(criticalTypeArray);

    vtkNew<vtkDoubleArray> scaleArray;
    scaleArray->SetNumberOfComponents(1);
    scaleArray->SetNumberOfTuples(nbPoints);
    scaleArray->SetName("scaleRadius");
    scaleArray->Fill(0);
    vtu->GetPointData()->AddArray(scaleArray);

    vtkNew<vtkDoubleArray> diagArray;
    diagArray->SetNumberOfComponents(1);
    diagArray->SetNumberOfTuples(nbPoints);
    diagArray->SetName("diag");
    diagArray->Fill(-1);
    vtu->GetPointData()->AddArray(diagArray);

    vtkNew<vtkDoubleArray> autoSimArray;
    autoSimArray->SetNumberOfComponents(1);
    autoSimArray->SetNumberOfTuples(nbPoints);
    autoSimArray->SetName("AutoSimilarity");
    autoSimArray->Fill(-1);
    vtu->GetPointData()->AddArray(autoSimArray);

    vtkNew<vtkDoubleArray> pers2dArray;
    pers2dArray->SetNumberOfComponents(1);
    pers2dArray->SetNumberOfTuples(nbPoints);
    pers2dArray->SetName("Persistence");
    pers2dArray->Fill(-1);
    vtu->GetPointData()->AddArray(pers2dArray);
    double* pers2dPtr = (double*)ttkUtils::GetVoidPointer(pers2dArray);

    vtkNew<vtkFloatArray> meanNormArray;
    meanNormArray->DeepCopy(vtu->GetPointData()->GetArray("Normals"));
    meanNormArray->SetName("MeanNormals");
    vtu->GetPointData()->AddArray(meanNormArray);
    float* meanNormalsPtr = (float*)ttkUtils::GetVoidPointer(meanNormArray);

    vtkNew<vtkIntArray> fullSizeArray{};
    fullSizeArray->SetNumberOfComponents(1);
    fullSizeArray->SetNumberOfTuples(nbPoints);
    fullSizeArray->SetName("Size");
    vtu->GetPointData()->AddArray(fullSizeArray);
    int* fullSizePtr = (int*)ttkUtils::GetVoidPointer(fullSizeArray);

    vtkNew<vtkIntArray> fullIdArray{};
    fullIdArray->SetNumberOfComponents(1);
    fullIdArray->SetNumberOfTuples(nbPoints);
    fullIdArray->SetName("fullId");
    vtu->GetPointData()->AddArray(fullIdArray);
    int* fullIdPtr = (int*)ttkUtils::GetVoidPointer(fullIdArray);

    vtkNew<vtkIntArray> pers2dIdArray{};
    pers2dIdArray->SetNumberOfComponents(1);
    pers2dIdArray->SetNumberOfTuples(nbPoints);
    pers2dIdArray->SetName("arcId");
    vtu->GetPointData()->AddArray(pers2dIdArray);
    pers2dIdArray->Fill(-1);
    int* pers2dIdPtr = static_cast<int*>(ttkUtils::GetVoidPointer(pers2dIdArray));

    vtkNew<vtkIntArray> primitiveTypeArray{};
    primitiveTypeArray->SetNumberOfComponents(1);
    primitiveTypeArray->SetNumberOfTuples(nbPoints);
    primitiveTypeArray->SetName("primitiveType");
    vtu->GetPointData()->AddArray(primitiveTypeArray);
    primitiveTypeArray->Fill(0);
    int* primitiveTypePtr = static_cast<int*>(ttkUtils::GetVoidPointer(primitiveTypeArray));

    std::vector<float> meanNormalsVec;
    computeSingleComponentMean(normalsPtr, meanNormalsVec, 3, nbPoints);

    SimplexId primitiveType=0;

    vtkNew<vtkFloatArray> covArray{};
    covArray->SetNumberOfComponents(3);
    covArray->SetNumberOfTuples(nbPoints);
    covArray->SetName("covariance");
    vtu->GetPointData()->AddArray(covArray);
    float* covPtr = (float*)ttkUtils::GetVoidPointer(covArray);

    std::vector<float> covVec;
    computeSinglePrimitiveType(normalsPtr, meanNormalsVec, nbPoints, primitiveType, covVec);


    for(int iPoint=0; iPoint<nbPoints; iPoint++){
      primitiveTypePtr[iPoint] = primitiveType;
      pers2dIdPtr[iPoint] = persistence2DId;
      pers2dPtr[iPoint] = persistence2D;
      fullSizePtr[iPoint] = nbPoints;
      scaleArray->SetTuple1(iPoint, scale);
      autoSimArray->SetTuple1(iPoint,autosimilarity);
      diagArray->SetTuple1(iPoint,diag);
      fullIdPtr[iPoint] = fullId;
      for(int i=0; i<3; i++){
        meanNormalsPtr[3*iPoint+i] = meanNormalsVec[i];
        covPtr[3*iPoint+i] = covVec[i];
      }
    }
  }

  void ttkRipsMultiParameterSampling::createOutputCommonComponentPoints(vtkUnstructuredGrid* vtu, SimplexId fullId, double percentCover, std::vector<std::vector<double>>& persistenceGrid, std::vector<double>& ratioVec, float* coordinates, float* normals){

    SimplexId nbCC = children2D_[fullId].size();
    SimplexId nCol = persistenceGrid[0].size();

    // std::vector<SimplexId> totalPointIds(0);
    std::map<SimplexId, double> pointCount;

    double total_weight=0;
    for(auto child : children2D_[fullId]){
      SimplexId localId, iRatio, iThreshold;
      getGridTuple(child, nCol, NbCCMax, localId, iRatio, iThreshold);
      SimplexId ccId = globalToLocal_[child];

      if(ccId>-1){

        std::vector<SimplexId> pointIds;
        trees_[iRatio].getComponentPointIdsAtThreshold(ccId, persistenceGrid[iRatio][iThreshold], pointIds);
        double weight = 0;
        if(iThreshold>0){
          weight = persistenceGrid[iRatio][iThreshold]-persistenceGrid[iRatio][iThreshold-1];
        }
        total_weight += weight;

        for(auto p : pointIds){
          auto it = pointCount.find(p);
          if(it == pointCount.end()){
            pointCount.emplace(p,weight);
          }else{
            it->second += weight;
          }
        }

      }
    }

    // std::sort(totalPointIds.begin(), totalPointIds.end());

    // SimplexId nbPoints = totalPointIds.size();

    auto pd = vtu->GetPointData();
    vtkNew<vtkPoints> points;
    vtu->SetPoints(points);
    vtu->Allocate();
    // points->SetNumberOfPoints(nbPoints);

    vtkNew<vtkFloatArray> normArray;
    normArray->SetName("Normals");
    normArray->SetNumberOfComponents(3);
    // normArray->SetNumberOfTuples(nbPoints);
    pd->AddArray(normArray);

    double criteria = total_weight*percentCover;
    // create geometry
    // for(int i=0; i<nbPoints; i++){
    for(auto& it : pointCount){
      if(it.second >= criteria){

        SimplexId id = it.first;
        auto i = points->InsertNextPoint(coordinates[3*id+0], coordinates[3*id+1], coordinates[3*id+2]);

        vtkNew<vtkIdList> idList;
        idList->InsertNextId(i);
        vtu->InsertNextCell(1, idList);

        normArray->InsertNextTuple(normals+3*id);
      }
    }
  }

  double ttkRipsMultiParameterSampling::findDiameter(float* coordinates, SimplexId nbPoints){
    Timer tm{};

    double diameter =0;
    for(SimplexId i=0; i<nbPoints; i++)  {
      for(SimplexId j=i+1; j<nbPoints; j++)  {
        double d=0;
        for(int k=0; k<3; k++){
          d += (coordinates[3*i+k]-coordinates[3*j+k])*(coordinates[3*i+k]-coordinates[3*j+k]);
        }

        if(d>diameter){
          diameter=d;
        }
      }
      return diameter;

    }

    printMsg("Found diameter", 1, tm.getElapsedTime(), threadNumber_);
  }

  void ttkRipsMultiParameterSampling::writeCoordinatesForTomato(
      std::string name,
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
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId){
    std::cout<<"write tomato"<<std::endl;

    std::ofstream myfile;
    myfile.open (name);

    SimplexId nbScales = sampling_indexes.size();
    SimplexId nbPoints = 0;
    for(int iScale=0; iScale<nbScales; iScale++){
      nbPoints += epsilonIdToFilteredId[iScale].size();
    }

    myfile << "x,y,z,nx,ny,nz,mls_x,mls_y,mls_z,f0,f1,f2,f3,f4,pca_x,pca_y,pca_z,d\n";
    for(SimplexId id_ms = 0; id_ms<nbPoints; id_ms++){
      int scale = findScaleId(id_ms, pointsAtScaleStartIdx);
      SimplexId id_eps = id_ms-pointsAtScaleStartIdx[scale];
      SimplexId id_fltrd = epsilonIdToFilteredId[scale][id_eps];
      SimplexId id_local = filteredToLocalIdVec[scale][id_fltrd];
      SimplexId id_global = sampling_indexes[scale][id_local];


      for(int k=0; k<3; k++){
        myfile << coordinates[3*id_global+k]<<", ";
      }
      for(int k=0; k<3; k++){
        myfile << normals[3*id_global+k]<<", ";
      }
      for(int k=0; k<3; k++){
        myfile << mlsCoordinates[scale][3*id_local+k]<<", ";
      }
      for(int k=0; k<5; k++){
        myfile << fitCoordinates[scale][5*id_local+k]<<", ";
      }
      for(int k=0; k<3; k++){
        myfile << embeddingCoordinates[scale][3*id_local+k]<<", ";
      }
      myfile << densityVec[id_ms]<<"\n";

    }
    myfile.close();
    std::cout<<"write tomato done"<<std::endl;
  }

  // int ttkRipsMultiParameterSampling::getSkeletonNodes(ttk::ftm::MergeTree<double>& mergeTree, vtkUnstructuredGrid *outputSkeletonNodes) {

  //   vtkNew<vtkUnstructuredGrid> skeletonNodes{};
  //   vtkNew<vtkPoints> points{};

  //   ttk::ftm::NodeData nodeData;
  //   nodeData.init(mergeTree, params_);
  //   nodeData.setScalarType(inputScalars_[0]->GetDataType());

  //   for(int cc = 0; cc < nbCC_; cc++) {
  //     auto tree = mergeTree[cc].tree.getTree(GetTreeType());
  //     vtkDataArray *idMapper
  //       = connected_components_[cc]->GetPointData()->GetArray(
  //         ttk::VertexScalarFieldName);

  //     const auto numberOfNodes = tree->getNumberOfNodes();
  // #ifndef TTK_ENABLE_KAMIKAZE
  //     if(!numberOfNodes) {
  //       this->printErr("Error : tree has no nodes.");
  //       return 0;
  //     }
  // #endif

  //     for(ttk::ftm::idNode nodeId = 0; nodeId < numberOfNodes; ++nodeId) {
  //       const auto node = tree->getNode(nodeId);
  // #ifndef TTK_ENABLE_KAMIKAZE
  //       if(!node) {
  //         this->printMsg({"Error : node ", std::to_string(nodeId), " is null."},
  //                        ttk::debug::Priority::ERROR);
  //         return 0;
  //       }
  // #endif
  //       const auto local_vertId = node->getVertexId();
  //       float point[3];
  //       triangulation_[cc]->getVertexPoint(
  //         local_vertId, point[0], point[1], point[2]);
  //       const auto nextPoint = points->InsertNextPoint(point);
  //       nodeData.fillArrayPoint(
  //         nextPoint, nodeId, ftmTree_[cc], idMapper, triangulation_[cc], params_);
  //     }
  //   }

  //   ttkUtils::CellVertexFromPoints(skeletonNodes, points);

  //   vtkPointData *pointData = skeletonNodes->GetPointData();
  //   nodeData.addArray(pointData, params_);

  //   outputSkeletonNodes->ShallowCopy(skeletonNodes);

  //   return 1;
  // }
  void ttkRipsMultiParameterSampling::ripsArcsToMergeTree(vtkMultiBlockDataSet* outputMTMB, std::vector<RipsArc>& ripsArcs, std::vector<SimplexId>& arcIds, std::vector<SimplexId>& leafIds){

    // std::cout<<" enter ripsArcsToMergeTree"<<std::endl;

    SimplexId nbArcs = ripsArcs.size();// + ripsArcs2.size();
                                      
    // std::cout<<"sort"<<std::endl;
    const auto cmp_pers = [=](const RipsArc &a0, const RipsArc &a1){
      return a0.ccId<a1.ccId or (a0.ccId == a1.ccId and (a0.death > a1.death));
    };

    TTK_PSORT(this->threadNumber_, ripsArcs.begin(), ripsArcs.end(), cmp_pers);

    

    // for(SimplexId i=0; i<nbArcs; i++){
    //   // Pour créer un noeud avec identifiant i (ça peut aller de 0 à numberOfNodes non compris)
    //   auto& arc = ripsArcs[i];
    //   std::cout<<" iarc "<<i<<" ccId "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<" ("<<arc.birth<<"-"<<arc.death<<"), pers "<<arc.death-arc.birth<<" type "<<arc.type<<"   treeId "<<arc.treeId<<", arcid "<<arcIds[i]<<std::endl;
    // }
    // std::cout<<nbArcs<<" arcs"<<std::endl;
    SimplexId numberOfNodes = nbArcs+1;
    // std::cout<<numberOfNodes<<" nodes"<<std::endl;

    SimplexId nbCC = 0;
    SimplexId oldNbCC = 0;
    double min_pers = std::numeric_limits<double>::max();


    // std::cout<<"find oldNbcc"<<std::endl;
    for(auto& arc : ripsArcs){
      if(arc.type<2){ // leaf arc
        if(arc.death-arc.birth<min_pers){
          min_pers = arc.death-arc.birth;
          oldNbCC = arc.deathId-arc.birthId;
          // std::cout<<" oldNbCC "<<oldNbCC<<" min_pers "<<min_pers<<std::endl;
        }else if(arc.death-arc.birth==min_pers){
          // std::cout<<" WARNING oldNbCC "<<oldNbCC<<" min_pers "<<min_pers<<std::endl;
        }
        nbCC++;
      }
    }
    // std::cout<<"old nbcc "<<oldNbCC<<std::endl;
    // std::cout<<nbCC<<" nbCC"<<std::endl;

    // =========
    // std::cout<<"find new ids"<<std::endl;
    std::vector<SimplexId> oldIdToNewId(2*oldNbCC,-1);
    int count=0;
    for(auto& arc : ripsArcs){
      if(arc.type<2){ // leaf arc
        oldIdToNewId[arc.birthId] = count;
        // oldIdToNewId[arc.deathId] = count + nbCC;
        // std::cout<<" ID "<<arc.birthId<<" -> "<<oldIdToNewId[arc.birthId]<<std::endl;
        count++;
      }
    }
    for(auto& arc : ripsArcs){
      oldIdToNewId[arc.deathId] = oldIdToNewId[arc.deathId-oldNbCC]+nbCC;
      // std::cout<<" ID "<<arc.deathId<<" -> "<<oldIdToNewId[arc.deathId]<<std::endl;
      if(arc.type==2){ // not a leaf arc
        oldIdToNewId[arc.birthId] = oldIdToNewId[arc.birthId-oldNbCC]+nbCC;
        // std::cout<<" ID "<<arc.birthId<<" -> "<<oldIdToNewId[arc.birthId]<<std::endl;
      }
    }

    // std::cout<<"set newid"<<std::endl;

    // for(auto& arc : ripsArcs){
    //   if(arc.type<2){ //leaf arc
    //     std::cout<<" type 0/1 "<<arc.birthId<<" - "<<arc.deathId<<"  ->  "<<oldIdToNewId[arc.birthId]<<"-"<<oldIdToNewId[arc.deathId-oldNbCC]+nbCC<<std::endl;
    //     arc.birthId = oldIdToNewId[arc.birthId];
    //     arc.deathId = oldIdToNewId[arc.deathId-oldNbCC]+nbCC;
    //   }else{
    //     std::cout<<" type  2 "<<arc.birthId<<" - "<<arc.deathId<<"  ->  "<<oldIdToNewId[arc.birthId-oldNbCC]+nbCC<<"-"<<oldIdToNewId[arc.deathId-oldNbCC]+nbCC<<std::endl;
    //     arc.birthId = oldIdToNewId[arc.birthId-oldNbCC]+nbCC;
    //     arc.deathId = oldIdToNewId[arc.deathId-oldNbCC]+nbCC;
    //   }
    // }

    // =========

    // std::cout<<"scalars"<<std::endl;

    auto scalars = std::make_shared<ttk::ftm::Scalars>();

    // std::cout<<"    values"<<std::endl;
    auto scalarsValues = std::make_shared<std::vector<double>>(numberOfNodes);
    // std::cout<<" ripsArcs.size() "<<ripsArcs.size()<<std::endl;
    for(SimplexId iArc=0; iArc<ripsArcs.size(); iArc++){
      auto& arc = ripsArcs[iArc];
      if(std::isnan(arc.birth)){
        std::cout<<" ERROR NAN VALUE"<<std::endl;
        std::cout<<" iarc "<<iArc<<" "<<oldIdToNewId[arc.birthId]<<"-"<<oldIdToNewId[arc.deathId]<<" outta "<<numberOfNodes<<std::endl;
      }
      (*scalarsValues)[oldIdToNewId[arc.birthId]] = arc.birth;
      (*scalarsValues)[oldIdToNewId[arc.deathId]] = arc.death;
    }
    for(SimplexId iNode=0; iNode<numberOfNodes; iNode++){
      double v = (*scalarsValues)[iNode];
      if(std::isnan(v)){
        std::cout<<" HERE ERROR NAN VALUE"<<std::endl;
      }
    }
    // Tu rempli le vecteur des valeurs scalaires ici avec (*scalarsValues)[index] = value
    scalars->size = numberOfNodes;
    scalars->values = (void *)(scalarsValues->data());

    auto params = std::make_shared<ttk::ftm::Params>();
    params->treeType = ftm::Join_Split;

    // auto tree = std::make_shared<ftm::FTMTree_MT>(params, scalars, ftm::Join_Split);

    // std::cout<<"build mt"<<std::endl;
    // params->treeType = ttk::ftm::TreeType::Join_Split;
    ttk::ftm::MergeTree<double> mergeTree(scalars, scalarsValues, params);
    // tree->makeAlloc();

    // std::cout<<"nodes"<<std::endl;
    for(SimplexId i=0; i<nbArcs; i++){
      // Pour créer un noeud avec identifiant i (ça peut aller de 0 à numberOfNodes non compris)
      auto& arc = ripsArcs[i];
      // std::cout<<" iarc "<<i<<" "<<oldIdToNewId[arc.birthId]<<"-"<<oldIdToNewId[arc.deathId]<<" type "<<arc.type<<"   treeId "<<arc.treeId<<std::endl;
      mergeTree.tree.makeNode(oldIdToNewId[arc.birthId]);
      mergeTree.tree.makeNode(oldIdToNewId[arc.deathId]);
    }
    // std::cout<<"arcs"<<std::endl;

    for(SimplexId iArc=0; iArc<nbArcs; iArc++){
      auto& arc = ripsArcs[iArc];
      // Pour créer un arc de l'enfant downId au parent upId
      // tree->makeSuperArc(arc.birthId,arc.deathId);

      if(arc.birthId>=oldIdToNewId.size() or arc.deathId>=oldIdToNewId.size()){
        std::cout<<"ERROR"<<std::endl;
      }
      // std::cout<<" iarc "<<iArc<<" "<<oldIdToNewId[arc.birthId]<<"-"<<oldIdToNewId[arc.deathId]<<" type "<<arc.type<<std::endl;
      mergeTree.tree.makeSuperArc(oldIdToNewId[arc.birthId],oldIdToNewId[arc.deathId]);
    }
    // std::cout<<"inconsistencies"<<std::endl;
    // Manage inconsistent arcs
    ttk::ftm::manageInconsistentArcsMultiParent(&(mergeTree.tree));
    // std::cout<<"remve selflink"<<std::endl;
    // Remove self link
    ttk::ftm::removeSelfLink(&(mergeTree.tree));
    // std::cout<<".type==-2pers pairs"<<std::endl;

    // std::cout<<"pairs"<<std::endl;
    ttk::ftm::computePersistencePairs<double>(&mergeTree.tree);
    // for(SimplexId iArc=0; iArc<nbArcs; iArc++){
    //   auto& arc = ripsArcs[iArc];

    //   auto bid = arc.birthId;
    //   auto did = arc.deathId;
    //   mergeTree.tree.getNode(bid)->setOrigin(ripsArcs.back().deathId);
    //   mergeTree.tree.getNode(did)->setOrigin(ripsArcs.back().deathId);
    // }

    // mergeTree.tree.getNode(ripsArcs.back().deathId)->setOrigin(ripsArcs[0].birthId);
    // std::cout<<"origins"<<std::endl;
    for(int iNode=0; iNode<nbCC; iNode++){
      // std::cout<<"i node "<<iNode<<"-"<<iNode+nbCC<<std::endl;
      mergeTree.tree.getNode(iNode+nbCC)->setOrigin(iNode);
      mergeTree.tree.getNode(iNode)->setOrigin(iNode+nbCC);
    }

    // std::cout<<"visu"<<std::endl;
    outputMTMB->SetNumberOfBlocks(2);
    vtkNew<vtkUnstructuredGrid> vtuNodes{};
    outputMTMB->SetBlock(0, vtuNodes);
    vtkNew<vtkUnstructuredGrid> vtuEdges{};
    outputMTMB->SetBlock(1, vtuEdges);

    // std::cout<<"visu maker"<<std::endl;
    ttkMergeTreeVisualization visuMaker;
    visuMaker.setPlanarLayout(true);
    visuMaker.setOutputSegmentation(false);
    visuMaker.setBranchDecompositionPlanarLayout(false);
    visuMaker.setImportantPairsSpacing(10);
    visuMaker.setVtkOutputNode(vtuNodes);
    visuMaker.setVtkOutputArc(vtuEdges);
    visuMaker.setDebugLevel(this->debugLevel_);
    visuMaker.makeTreesOutput<double>(&mergeTree.tree);



    // std::cout<<"persistence balck magic for reordering arcs"<<std::endl;
    // Add arcId sorted by persistences:
    // std::cout<<"ordering"<<std::endl;

    std::vector<SimplexId> arc_pers_indexing(ripsArcs.size());
    std::iota(arc_pers_indexing.begin(), arc_pers_indexing.end(), 0);
    const auto cmp2 = [=](const SimplexId &i0, const SimplexId &i1){
      auto& a0 = ripsArcs[i0];
      auto& a1 = ripsArcs[i1];
      return ((a0.death-a0.birth) > (a1.death-a1.birth)) || (((a0.death-a0.birth) == (a1.death-a1.birth)) && a0.ccId<a1.ccId);
    };
    const auto cmp_death = [=](const SimplexId &i0, const SimplexId &i1){
      auto& a0 = ripsArcs[i0];
      auto& a1 = ripsArcs[i1];
      return ((a0.death) > (a1.death)) || (((a0.death) == (a1.death)) && a0.ccId<a1.ccId);
    };
    TTK_PSORT(this->threadNumber_, arc_pers_indexing.begin(), arc_pers_indexing.end(), cmp_death);
    // std::cout<<" sort done "<<std::endl;

    // const auto cmp = [=](const RipsArc &a0, const RipsArc &a1){
    //   return ((a0.death-a0.birth) > (a1.death-a1.birth)) || (((a0.death-a0.birth) == (a1.death-a1.birth)) && a0.ccId<a1.ccId);
    // };
    // TTK_PSORT(this->threadNumber_, ripsArcs.begin(), ripsArcs.end(), cmp);

    // std::cout<<"again"<<std::endl;
    // get edges from tree:
    int nEdges = vtuEdges->GetNumberOfCells();
    std::vector<SimplexId> arcIdToEdgeId(nbArcs, -1); //nEdges > nArcs
    // arcIdToEdgeId = arc_pers_indexing;
    // std::iota(arcIdToEdgeId.begin(), arcIdToEdgeId.end(), 0);
    std::vector<double> edgeToPers(nEdges, 0);
    std::vector<double> edgeToBranchPers(nEdges, 0);
    std::vector<double> edgeToDeath(nEdges, 123);
    std::vector<int> edgeBranchId(nEdges, -1);
    std::vector<int> edgeIsDummy(nEdges, 0);

    std::vector<SimplexId> downNodeToEdge(nbArcs);

    for(int iEdge=0; iEdge<nEdges; iEdge++){
      int i0 =  vtuEdges->GetCell(iEdge)->GetPointId(0);
      int i1 =  vtuEdges->GetCell(iEdge)->GetPointId(1);
      double p0[3];
      vtuEdges->GetPoint(i0, p0);
      double p1[3];
      vtuEdges->GetPoint(i1, p1);
      double pers  = std::abs(p0[1]-p1[1]);
      edgeToPers[iEdge] = pers;
      edgeIsDummy[iEdge] = vtuEdges->GetCellData()->GetArray("isDummyArc")->GetTuple1(iEdge);
      edgeToDeath[iEdge] = std::max(p0[1], p1[1]);
      edgeBranchId[iEdge] = vtuEdges->GetCellData()->GetArray("BranchID")->GetTuple1(iEdge);
      edgeToBranchPers[iEdge] = vtuEdges->GetCellData()->GetArray("Persistence")->GetTuple1(iEdge);
      SimplexId trueDownNodeId = vtuEdges->GetCellData()->GetArray("trueDownNodeId")->GetTuple1(iEdge);

      if(edgeIsDummy[iEdge]==0){
        downNodeToEdge[trueDownNodeId] = iEdge;
      }

    }
    for(int iArc=0; iArc<ripsArcs.size(); iArc++){
      auto& arc = ripsArcs[iArc];
      arcIdToEdgeId[iArc] = downNodeToEdge[oldIdToNewId[arc.birthId]];
    }
    // std::cout<<"here"<<std::endl;
    // std::sort(arcIdToEdgeId.begin(), arcIdToEdgeId.end(),
    //     [&](int A, int B) -> bool {
    //     return edgeToPers[A] > edgeToPers[B] or (edgeToPers[A] == edgeToPers[B] 
    //         and edgeBranchId[A]<edgeBranchId[B]);
    //     });
    // std::cout<<" sorting "<<std::endl;
    // std::sort(arcIdToEdgeId.begin(), arcIdToEdgeId.end(),
    //     [&](int A, int B) -> bool {
    //         return (edgeIsDummy[A]<edgeIsDummy[B] or (edgeIsDummy[A]==edgeIsDummy[B]
    //         and (edgeToBranchPers[A] > edgeToBranchPers[B] or (edgeToBranchPers[A] == edgeToBranchPers[B]
    //         and edgeToDeath[A] > edgeToDeath[B] or (edgeToDeath[A] == edgeToPers[B] 
    //         and edgeBranchId[A]<edgeBranchId[B])))));
    //     });
    // std::cout<<" done \n\n"<<std::endl;
    for(int iArc = 0; iArc<nbArcs; iArc++){
      auto& arc = ripsArcs[iArc];
      // std::cout<<" iarc "<<iArc<<" ccId "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<" ("<<arc.birth<<"-"<<arc.death<<"), pers "<<arc.death-arc.birth<<" type "<<arc.type<<"   treeId "<<arc.treeId<<", arcid "<<arcIds[iArc]<<std::endl;
      SimplexId iEdge = arcIdToEdgeId[iArc];
      int i0 =  vtuEdges->GetCell(iEdge)->GetPointId(0);
      int i1 =  vtuEdges->GetCell(iEdge)->GetPointId(1);
      double p0[3];
      vtuEdges->GetPoint(i0, p0);
      double p1[3];
      vtuEdges->GetPoint(i1, p1);
      double b  = std::min(p0[1], p1[1]);
      double d  = std::max(p0[1], p1[1]);
      double pers = edgeToPers[iEdge];
      int dummy = edgeIsDummy[iEdge];
      int branchid = edgeBranchId[iEdge];
      double branchPers = edgeToBranchPers[iEdge];
      int trueDownNodeId = vtuEdges->GetCellData()->GetArray("trueDownNodeId")->GetTuple1(iEdge);
      // std::cout<<"     with edge "<<iEdge<<" branchid "<<branchid<<" "<<" ("<<b<<"-"<<d<<"), pers "<<pers<<" branchPers "<<branchPers<< " trueDownNodeId "<<trueDownNodeId<<std::endl;
      if(trueDownNodeId != oldIdToNewId[arc.birthId]){
        std::cout<<"ERROR"<<std::endl;
      }



     }


    // }
    // for(int ie : arcIdToEdgeId){
      // std::cout<<" edge "<<ie<<" pers "<<edgeToPers[ie]<<std::endl;
    // }
    //
    //
    std::vector<SimplexId> arc_branch_indexing(ripsArcs.size());
    std::iota(arc_branch_indexing.begin(), arc_branch_indexing.end(), 0);
    const auto cmp_arc_type = [=](const SimplexId &i0, const SimplexId &i1){
      return ripsArcs[i0].type<ripsArcs[i1].type;
    };
    const auto cmp_birth = [=](const SimplexId &i0, const SimplexId &i1){
      auto& a0 = ripsArcs[i0];
      auto& a1 = ripsArcs[i1];
      return a0.birth<a1.birth;
    };
    TTK_PSORT(this->threadNumber_, arc_branch_indexing.begin(), arc_branch_indexing.end(), cmp_arc_type);
    TTK_PSORT(this->threadNumber_, arc_branch_indexing.begin()+nbCC, arc_branch_indexing.end(), cmp_birth);
    // std::vector<SimplexId> arcIdToBranchId(nbArcs);
    std::vector<SimplexId> branchIdToArcId(nbCC);

    // std::cout<<"branch id to arc id"<<std::endl;
    for(int iArc=0; iArc<nbArcs; iArc++){
      // arcIdToBranchId[iArc] = edgeBranchId[arcIdToEdgeId[iArc]];
      if(ripsArcs[iArc].type<2){
        branchIdToArcId[edgeBranchId[arcIdToEdgeId[iArc]]]=iArc;

        auto& arc = ripsArcs[iArc];
        // std::cout<<" iArc "<<iArc<<std::endl;
        // std::cout<<"    edge id "<<arcIdToEdgeId[iArc]<<std::endl;
        // std::cout<<"    branch id "<<edgeBranchId[arcIdToEdgeId[iArc]]<<std::endl;
        // std::cout<<"    iarc "<<iArc<<" ccId "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<" ("<<arc.birth<<"-"<<arc.death<<"), pers "<<arc.death-arc.birth<<" type "<<arc.type<<"   treeId "<<arc.treeId<<", arcid "<<arcIds[iArc]<<std::endl;

        // std::cout<<" BRANCH ID = "<<edgeBranchId[arcIdToEdgeId[iArc]]<<std::endl;
      }
    }
    // for(int iBranch=0; iBranch<nbCC; iBranch++){
    //   branchIdToArcId[]
    //   auto& arc = ripsArcs[arc_branch_indexing[i]];
    //   std::cout<<" iarc "<<arc_branch_indexing[i]<<" ccId "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<" ("<<arc.birth<<"-"<<arc.death<<"), pers "<<arc.death-arc.birth<<" type "<<arc.type<<"   treeId "<<arc.treeId<<", arcid "<<arcIds[i]<<std::endl;
    // }

    std::vector<int> nodeType(nbCC*2, 3);
    std::vector<int> edgeMinDownType(nEdges,3);

    // std::cout<<" propagate type "<<std::endl;
    for(int i = 0; i < nbCC; i++){
      SimplexId iArc = arc_branch_indexing[i];
      auto& arc = ripsArcs[iArc];
      // std::cout<<"    iarc "<<iArc<<" ccId "<<arc.ccId<<" "<<arc.birthId<<"-"<<arc.deathId<<" ("<<arc.birth<<"-"<<arc.death<<"), pers "<<arc.death-arc.birth<<" type "<<arc.type<<"   treeId "<<arc.treeId<<", arcid "<<arcIds[iArc]<<std::endl;

      SimplexId iNode = oldIdToNewId[arc.birthId];
      if(arc.type < nodeType[iNode]){
        nodeType[iNode]=arc.type;
        // std::cout<<" type of "<<iNode<<" is "<<arc.type<<std::endl;
      }
      if(arc.type < nodeType[iNode+nbCC]){
        nodeType[iNode+nbCC]=arc.type;
        // std::cout<<" type of "<<iNode+nbCC<<" is "<<arc.type<<std::endl;
      }
      if(arc.type < nodeType[oldIdToNewId[arc.deathId]]){
        nodeType[oldIdToNewId[arc.deathId]]=arc.type;
        // std::cout<<" type of "<<oldIdToNewId[arc.deathId]<<" is "<<arc.type<<std::endl;
      }
    }

    // propagate minnodetype up the nodes
    for(auto& i = arc_branch_indexing[nbCC]; i<arc_branch_indexing.size(); i++){
      auto& arc = ripsArcs[arc_branch_indexing[i]];
      SimplexId bid = oldIdToNewId[arc.birthId];
      SimplexId did = oldIdToNewId[arc.deathId];
      if(nodeType[bid]<nodeType[did]){
        nodeType[did]=nodeType[bid];
      }
    }
    for(int iEdge = 0; iEdge < nEdges; iEdge++){
      SimplexId trueDownNodeId = vtuEdges->GetCellData()->GetArray("trueDownNodeId")->GetTuple1(iEdge);
      SimplexId downType = nodeType[trueDownNodeId];
      if(edgeMinDownType[iEdge]>downType){
        edgeMinDownType[iEdge] = downType;
        // std::cout<<" edge "<<iEdge<<" has downNode "<<trueDownNodeId<<" of type "<< downType<<std::endl;
      }
    }

    // std::cout<<" vtk arrays"<<std::endl;
    // set correct arcId
    vtkNew<vtkIntArray> nodeIdArray{};
    nodeIdArray->SetNumberOfComponents(1);
    nodeIdArray->SetNumberOfTuples(nEdges);
    nodeIdArray->SetName("nodeId");
    nodeIdArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(nodeIdArray);

    vtkNew<vtkIntArray> localArcIdArray{};
    localArcIdArray->SetNumberOfComponents(1);
    localArcIdArray->SetNumberOfTuples(nEdges);
    localArcIdArray->SetName("localArcId");
    localArcIdArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(localArcIdArray);

    vtkNew<vtkIntArray> arcIdArray{};
    arcIdArray->SetNumberOfComponents(1);
    arcIdArray->SetNumberOfTuples(nEdges);
    arcIdArray->SetName("arcId");
    arcIdArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(arcIdArray);

    vtkNew<vtkIntArray> chosenArcIdArray{};
    chosenArcIdArray->SetNumberOfComponents(1);
    chosenArcIdArray->SetNumberOfTuples(nEdges);
    chosenArcIdArray->SetName("chosenArcId");
    chosenArcIdArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(chosenArcIdArray);

    vtkNew<vtkIntArray> typeArray{};
    typeArray->SetNumberOfComponents(1);
    typeArray->SetNumberOfTuples(nEdges);
    typeArray->SetName("type2");
    typeArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(typeArray);

    vtkNew<vtkIntArray> leafArcIdArray{};
    leafArcIdArray->SetNumberOfComponents(1);
    leafArcIdArray->SetNumberOfTuples(nEdges);
    leafArcIdArray->SetName("leafArcId");
    leafArcIdArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(leafArcIdArray);

    vtkNew<vtkIntArray> leafCCIdIdArray{};
    leafCCIdIdArray->SetNumberOfComponents(1);
    leafCCIdIdArray->SetNumberOfTuples(nEdges);
    leafCCIdIdArray->SetName("leafCCIdId");
    leafCCIdIdArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(leafCCIdIdArray);

    vtkNew<vtkIntArray> leafArcTypeArray{};
    leafArcTypeArray->SetNumberOfComponents(1);
    leafArcTypeArray->SetNumberOfTuples(nEdges);
    leafArcTypeArray->SetName("leafArcType");
    leafArcTypeArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(leafArcTypeArray);

    vtkNew<vtkIntArray> minDownTypeArray{};
    minDownTypeArray->SetNumberOfComponents(1);
    minDownTypeArray->SetNumberOfTuples(nEdges);
    minDownTypeArray->SetName("minDownType");
    minDownTypeArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(minDownTypeArray);

    vtkNew<vtkDoubleArray> diagArray{};
    diagArray->SetNumberOfComponents(1);
    diagArray->SetNumberOfTuples(nEdges);
    diagArray->SetName("diag");
    diagArray->Fill(-1);
    vtuEdges->GetCellData()->AddArray(diagArray);

    for(SimplexId iEdge=0;iEdge<nEdges; iEdge++){
      leafArcIdArray->SetTuple1(iEdge, arcIds[branchIdToArcId[edgeBranchId[iEdge]]]);
      leafArcTypeArray->SetTuple1(iEdge, ripsArcs[branchIdToArcId[edgeBranchId[iEdge]]].type);
      // leafCCIdIdArray->SetTuple1(iEdge, ripsArcs[branchIdToArcId[edgeBranchId[iEdge]]].ccId);
      minDownTypeArray->SetTuple1(iEdge, edgeMinDownType[iEdge]);
    }
    for(SimplexId iArc=0;iArc<nbArcs; iArc++){
      chosenArcIdArray->SetTuple1(arcIdToEdgeId[iArc], leafIds[iArc]);
      arcIdArray->SetTuple1(arcIdToEdgeId[iArc], arcIds[iArc]);
      localArcIdArray->SetTuple1(arcIdToEdgeId[iArc], iArc);
      typeArray->SetTuple1(arcIdToEdgeId[iArc], ripsArcs[iArc].type);
      diagArray->SetTuple1(arcIdToEdgeId[iArc], ripsArcs[iArc].diag);
      nodeIdArray->SetTuple1(arcIdToEdgeId[iArc], ripsArcs[iArc].birthId);
      leafCCIdIdArray->SetTuple1(arcIdToEdgeId[iArc], ripsArcs[iArc].ccId);

      vtuEdges->GetCellData()->GetArray("Persistence")->SetTuple1(arcIdToEdgeId[iArc], edgeToPers[arcIdToEdgeId[iArc]]);
    }

    //         vtkNew<vtkIntArray> sizeArray{};
    //         sizeArray->SetNumberOfComponents(1);
    //         sizeArray->SetNumberOfTuples(nEdges);
    //         sizeArray->SetName("size");
    //         sizeArray->Fill(-1);
    //         vtuEdges->GetCellData()->AddArray(sizeArray);

    //         for(SimplexId iArc=0;iArc<nbArcs; iArc++){
    //           SimplexId ccId = ripsArcs[iArc].ccId;
    //           double threshold = ripsArcs[iArc].death;
    //           SimplexId size = trees_[iScale].sizeAtThreshold(ccId, threshold);
    //           sizeArray->SetTuple1(arcIdToEdgeId[iArc], size);
    //         }

    // }
  }


  std::string ttkRipsMultiParameterSampling::generateFittingStringId(float* coordinates, float* normals, SimplexId nbPoints){
    std::vector<double> dataset_values;
    dataset_values.push_back(nbPoints);
    for(int i=0; i<6; i++){
      dataset_values.push_back(coordinates[i]);
    }
    if(normals){
      for(int i=0; i<6; i++){
        dataset_values.push_back(normals[i]);
      }
    }
    std::stringstream dataset_ss;
    dataset_ss << valuesToStringStream(dataset_values).str();

    std::vector<int> fitting_params_bool;
    std::vector<int> fitting_params_int;
    std::vector<double> fitting_params_double;

    fitting_params_bool.push_back((int)UseSampling);
    fitting_params_bool.push_back(UseNormalizedCoords);
    fitting_params_bool.push_back(ChangeBasis);
    fitting_params_bool.push_back(UseKnnGraph);
    fitting_params_bool.push_back(UseRIMLS);
    fitting_params_bool.push_back(UseProjectiveNorm);

    fitting_params_int.push_back(FitDimension);
    fitting_params_int.push_back(std::min(NbPointsMax, nbPoints));
    fitting_params_int.push_back(ScaleSamplingCount);

    fitting_params_double.push_back(ScaleSamplingMin);
    fitting_params_double.push_back(ScaleSamplingMax);

    // for(auto v_bool : fitting_params_bool){
    //   std::cout<<" "<<v_bool;
    // }
    // std::cout<<std::endl;

    std::stringstream fitting_params_ss;
    fitting_params_ss << boolValuesToStringStream(fitting_params_bool).str();
    fitting_params_ss << intValuesToStringStream(fitting_params_int, 9).str();
    fitting_params_ss << doubleValuesToStringStream(fitting_params_double, 4, 3).str();

    // std::cout<<" string ids: "<<dataset_ss.str()<<" "<<fitting_params_ss.str()<<std::endl;
    fittingStringId_ = dataset_ss.str()+"_"+fitting_params_ss.str();
    return fittingStringId_;
  }

  // std::string ttkRipsMultiParameterSampling::generateEpsSamplingStringId(){

  //   std::vector<double> epsSampling_params_values;
  //   epsSampling_params_values.push_back(Epsilon);
  //   epsSampling_params_values.push_back(EpsilonNbMaxPoint);

  //   std::stringstream epsSampling_params_ss;
  //   epsSampling_params_ss << valuesToStringStream(epsSampling_params_values).str();

  //   // std::cout<<" string ids: "<<dataset_ss.str()<<" "<<epsSampling_params_ss.str()<<std::endl;
  //   epsSamplingStringId_ = fittingStringId()+"_"+epsSampling_params_ss.str();
  //   return epsSamplingStringId_;
  // }

  bool ttkRipsMultiParameterSampling::loadDensity(){
    std::stringstream filename;
    filename<<savefile_path_<<"savefile_"<<fittingStringId()<<".dens";

    Timer tm_read;
    std::cout<<"Loading density results from "<<filename.str()<<std::endl;

    std::ifstream file_density(filename.str());

    if(file_density.good()){
      readVectorFromStream(densityVec_, file_density); 
      file_density.close();

      std::cout<<" done in "<<tm_read.getElapsedTime()<<" s."<<std::endl;
      return true;
    }else{
      std::cout<<" couldn't find file -- skipping"<<std::endl;
      return false;
    }
  }

  void ttkRipsMultiParameterSampling::saveDensity(){
    std::stringstream filename;
    filename<<savefile_path_<<"savefile_"<<fittingStringId()<<".dens";

    Timer tm_save;
    std::cout<<"Writing density results to "<<filename.str()<<"...";

    std::ofstream file_density(filename.str());
    vectorToOfStream(densityVec_, file_density); 
    file_density.close();

    std::cout<<" done in "<<tm_save.getElapsedTime()<<" s."<<std::endl;
    return; 
  }



  bool ttkRipsMultiParameterSampling::loadFits(){
    std::stringstream filename;
    filename<<savefile_path_<<"savefile_"<<fittingStringId()<<".fit";

    Timer tm_read_fit;
    std::cout<<"Loading fitting results from "<<filename.str()<<std::endl;

    std::ifstream file_fitting(filename.str());

    if(file_fitting.good()){
      readMultiscaleVector(sampling_indexes_, file_fitting); 
      readMultiscaleVector(fitCoordsVec_, file_fitting); 
      readMultiscaleVector(mlsCoordsVec_, file_fitting); 
      readMultiscaleVector(mlsNormalsVec_, file_fitting); 
      readMultiscaleVector(filteredToLocalId_, file_fitting); 
      readMultiscaleVector(isStableVec_, file_fitting); 
      readVectorFromStream(resultBarycenter_, file_fitting); 
      file_fitting.close();

      std::cout<<" done in "<<tm_read_fit.getElapsedTime()<<" s."<<std::endl;
      return true;
    }else{
      std::cout<<" couldn't find file -- skipping"<<std::endl;
      return false;
    }
  }

  void ttkRipsMultiParameterSampling::saveFits(){
    std::stringstream filename;
    filename<<savefile_path_<<"savefile_"<<fittingStringId()<<".fit";

    Timer tm_save_fit;
    std::cout<<"Writing fitting results to "<<filename.str()<<"...";

    std::ofstream file_fitting(filename.str());
    writeMultiscaleVector(sampling_indexes_, file_fitting); 
    writeMultiscaleVector(fitCoordsVec_, file_fitting); 
    writeMultiscaleVector(mlsCoordsVec_, file_fitting); 
    writeMultiscaleVector(mlsNormalsVec_, file_fitting); 
    writeMultiscaleVector(filteredToLocalId_, file_fitting); 
    writeMultiscaleVector(isStableVec_, file_fitting); 
    vectorToOfStream(resultBarycenter_, file_fitting); 
    file_fitting.close();

    std::cout<<" done in "<<tm_save_fit.getElapsedTime()<<" s."<<std::endl;
    return; 
  }

  // void ttkRipsMultiParameterSampling::saveEpsilonSampling(){
  //   std::stringstream filename;
  //   filename<<savefile_path_<<"savefile_"<<epsSamplingStringId()<<".epsSampling";

  //   Timer tm_save_epsSampling;
  //   std::cout<<"Writing Epsilon Sampling results to "<<filename.str()<<"...";

  //   std::ofstream file_epsSampling(filename.str());
  //   writeMultiscaleVector(sampling_indexes_, file_epsSampling); 
  //   vectorToOfStream(resultBarycenter_, file_epsSampling); 
  //   file_epsSampling.close();

  //   std::cout<<" done in "<<tm_save_epsSampling.getElapsedTime()<<" s."<<std::endl;
  //   return; 
  // }

template<class PointND>
bool ttkRipsMultiParameterSampling::createOutputMergingEdges(vtkPolyData* vtp,
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
                                              std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSamplers,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId){

	int dim = PointND::Dim;

  SimplexId nbScales = sampling_indexes.size();
  if(epsilonIdToFilteredId.size()!=nbScales or filteredToLocalIdVec.size()!=nbScales or pointsAtScaleStartIdx.size()!=nbScales){
    std::cout<<"SIZE DOES NOT MATCH SCALE NUMBER"<<std::endl;
  }
  SimplexId nbPoints = 0;
  for(int iScale=0; iScale<nbScales; iScale++){
    nbPoints += epsilonIdToFilteredId[iScale].size();
  }

  vtkNew<vtkFloatArray> barycenterFDArray{};
  barycenterFDArray->SetNumberOfComponents(3);
  barycenterFDArray->SetNumberOfTuples(1);
  barycenterFDArray->SetTuple(0, resultBarycenter.data());
  barycenterFDArray->SetName("basisCenter");
  vtp->GetFieldData()->AddArray(barycenterFDArray);

  std::vector<std::vector<SimplexId>> toAdd(nbScales);


  bool compute3DEmbedding = embeddingCoordinates.size()==nbScales;

  auto pd = vtp->GetPointData();

  vtkNew<vtkIntArray> scaleIdArray;
  scaleIdArray->SetName("scaleId");
  scaleIdArray->SetNumberOfComponents(1);
  scaleIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(scaleIdArray);

  vtkNew<vtkIntArray> epsIdArray;
  epsIdArray->SetName("epsId");
  epsIdArray->SetNumberOfComponents(1);
  epsIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(epsIdArray);

  vtkNew<vtkIntArray> localIdArray;
  localIdArray->SetName("localId");
  localIdArray->SetNumberOfComponents(1);
  localIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(localIdArray);

  vtkNew<vtkIntArray> idArray;
  idArray->SetName("originalIdx");
  idArray->SetNumberOfComponents(1);
  idArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(idArray);

  vtkNew<vtkDoubleArray> densityArray;
  densityArray->SetName("density");
  densityArray->SetNumberOfComponents(1);
  densityArray->SetNumberOfTuples(nbPoints);
  densityArray->Fill(-1);
  pd->AddArray(densityArray);

  vtkNew<vtkIntArray> filteredIdArray;
  filteredIdArray->SetName("filteredIdx");
  filteredIdArray->SetNumberOfComponents(1);
  filteredIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(filteredIdArray);


  vtkNew<vtkFloatArray> normArray;
  normArray->SetName("Normals");
  normArray->SetNumberOfComponents(3);
  normArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(normArray);


  vtkNew<vtkFloatArray> mlsCoordsArray;
  mlsCoordsArray->SetName("mlsCoords");
  mlsCoordsArray->SetNumberOfComponents(3);
  mlsCoordsArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(mlsCoordsArray);

  vtkNew<vtkFloatArray> mlsNormalsArray;
  mlsNormalsArray->SetName("Normals_mls");
  mlsNormalsArray->SetNumberOfComponents(3);
  mlsNormalsArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(mlsNormalsArray);

  vtkNew<vtkDoubleArray> fitCoordsArray;
  fitCoordsArray->SetName("FitCoords");
  fitCoordsArray->SetNumberOfComponents(dim);
  fitCoordsArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(fitCoordsArray);

  vtkNew<vtkIntArray> isEpsilonSampledArray;
  isEpsilonSampledArray->SetName("isEpsilonSampled");
  isEpsilonSampledArray->SetNumberOfComponents(1);
  isEpsilonSampledArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(isEpsilonSampledArray);



  vtkNew<vtkDoubleArray> embedCoordsArray;
  if(compute3DEmbedding){
    embedCoordsArray->SetName("3d-Embedding");
    embedCoordsArray->SetNumberOfComponents(3);
    embedCoordsArray->SetNumberOfTuples(nbPoints);
    pd->AddArray(embedCoordsArray);
  }

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nbPoints);
  vtp->SetPoints(points);
	vtp->Allocate(mergingEdges.size());

  for(SimplexId id_ms =0; id_ms<nbPoints; id_ms++){

    int scale = findScaleId(id_ms, pointsAtScaleStartIdx);
    SimplexId id_eps = id_ms-pointsAtScaleStartIdx[scale];
    SimplexId id_fltrd = epsilonIdToFilteredId[scale][id_eps];
    SimplexId id_local = filteredToLocalIdVec[scale][id_fltrd];
    SimplexId id_global = sampling_indexes[scale][id_local];

    // geometry
    points->SetPoint(id_ms, coordinates[3*id_global+0], coordinates[3*id_global+1], coordinates[3*id_global+2]);
    vtkNew<vtkIdList> idList;
		idList->InsertNextId(id_ms);
		// vtp->InsertNextCell(1, idList);
		
		// arrays
		normArray->SetTuple(id_ms, normals+3*id_global);
		fitCoordsArray->SetTuple(id_ms, fitCoordinates[scale].data()+dim*id_local);
		mlsCoordsArray->SetTuple(id_ms,mlsCoordinates[scale].data()+3*id_local);
		densityArray->SetTuple(id_ms,densityVec.data()+id_ms);
    mlsNormalsArray->SetTuple(id_ms,mlsNormals[scale].data()+3*id_local);
    if(compute3DEmbedding){
      embedCoordsArray->SetTuple(id_ms,embeddingCoordinates[scale].data()+3*id_local);
    }

    isEpsilonSampledArray->SetTuple1(id_ms,1);
    idArray->SetTuple1(id_ms,id_global);
    epsIdArray->SetTuple1(id_ms,id_eps);
    localIdArray->SetTuple1(id_ms,id_local);
    filteredIdArray->SetTuple1(id_ms,id_fltrd);
    scaleIdArray->SetTuple1(id_ms, scale);

    // find children in epsilon sampling
    auto& children = epsilonSamplers[scale].getChildren(id_eps);
    for(SimplexId iChild=1; iChild<children.size(); iChild++){
      auto child = children[iChild];
      toAdd[scale].push_back(child);
    }
  }

  auto cd = vtp->GetCellData(); // celldata
	SimplexId nbEdges = mergingEdges.size();
	
  vtkNew<vtkIntArray> scaleTypeArray;
  scaleTypeArray->SetName("scaleType");
  scaleTypeArray->SetNumberOfComponents(1);
  scaleTypeArray->SetNumberOfTuples(nbEdges);
  cd->AddArray(scaleTypeArray);

  vtkNew<vtkIntArray> mergeOrderArray;
  mergeOrderArray->SetName("mergeOrder");
  mergeOrderArray->SetNumberOfComponents(1);
  mergeOrderArray->SetNumberOfTuples(nbEdges);
  cd->AddArray(mergeOrderArray);

  vtkNew<vtkDoubleArray> lengthArray;
  lengthArray->SetName("length");
  lengthArray->SetNumberOfComponents(1);
  lengthArray->SetNumberOfTuples(nbEdges);
  cd->AddArray(lengthArray);


	

	for(SimplexId iEdge = 0; iEdge<nbEdges; iEdge++){
		auto& e = mergingEdges[iEdge];
		SimplexId id_ms_i = std::get<0>(e);
		SimplexId id_ms_j = std::get<1>(e);
		double length = std::get<2>(e);

		int scale_i = findScaleId(id_ms_i, pointsAtScaleStartIdx);
		SimplexId id_eps_i = id_ms_i-pointsAtScaleStartIdx[scale_i];
		SimplexId id_fltrd_i = epsilonIdToFilteredId[scale_i][id_eps_i];
		SimplexId id_local_i = filteredToLocalIdVec[scale_i][id_fltrd_i];
		SimplexId id_global_i = sampling_indexes[scale_i][id_local_i];

		int scale_j = findScaleId(id_ms_j, pointsAtScaleStartIdx);
		SimplexId id_eps_j = id_ms_j-pointsAtScaleStartIdx[scale_j];
		SimplexId id_fltrd_j = epsilonIdToFilteredId[scale_j][id_eps_j];
		SimplexId id_local_j = filteredToLocalIdVec[scale_j][id_fltrd_j];
		SimplexId id_global_j = sampling_indexes[scale_j][id_local_j];

		scaleTypeArray->SetTuple1(iEdge, (scale_i == scale_j ? 0 : 1));
		lengthArray->SetTuple1(iEdge, length);
		mergeOrderArray->SetTuple1(iEdge, iEdge);


		vtkNew<vtkIdList> idList;
		idList->InsertNextId(id_ms_i);
		idList->InsertNextId(id_ms_j);
		vtp->InsertNextCell(3, idList);
	}



  // add non-epsilonsampled points
  SimplexId nbExtraPoints = 0;
  for(int iScale=0; iScale<nbScales; iScale++){
    nbExtraPoints+=toAdd[iScale].size();
  }
  SimplexId nbTotalPoints = nbPoints+nbExtraPoints;


  // std::cout<<"currently "<<points->GetNumberOfPoints()<<" points"<<std::endl;

  // std::cout<<"  size  "<<nbPoints<<" -> "<<nbTotalPoints<<std::endl;
  // Reallocation for extra points
  points->Resize(nbTotalPoints);
  normArray->Resize(nbTotalPoints);
  fitCoordsArray->Resize(nbTotalPoints);
  filteredIdArray->Resize(nbTotalPoints);
  idArray->Resize(nbTotalPoints);
  scaleIdArray->Resize(nbTotalPoints);
  mlsCoordsArray->Resize(nbTotalPoints);
  mlsNormalsArray->Resize(nbTotalPoints);
  epsIdArray->Resize(nbTotalPoints);
  isEpsilonSampledArray->Resize(nbTotalPoints);
  localIdArray->Resize(nbTotalPoints);
  densityArray->Resize(nbTotalPoints);

  points->SetNumberOfPoints(nbTotalPoints);
  normArray->SetNumberOfTuples(nbTotalPoints);
  fitCoordsArray->SetNumberOfTuples(nbTotalPoints);
  filteredIdArray->SetNumberOfTuples(nbTotalPoints);
  idArray->SetNumberOfTuples(nbTotalPoints);
  scaleIdArray->SetNumberOfTuples(nbTotalPoints);
  mlsCoordsArray->SetNumberOfTuples(nbTotalPoints);
  mlsNormalsArray->SetNumberOfTuples(nbTotalPoints);
  epsIdArray->SetNumberOfTuples(nbTotalPoints);
  isEpsilonSampledArray->SetNumberOfTuples(nbTotalPoints);
  localIdArray->SetNumberOfTuples(nbTotalPoints);
  densityArray->SetNumberOfTuples(nbTotalPoints);

  if(compute3DEmbedding){
    embedCoordsArray->Resize(nbTotalPoints);
    embedCoordsArray->SetNumberOfTuples(nbTotalPoints);
  }

  SimplexId i = nbPoints;
  for(int iScale=0; iScale<nbScales; iScale++){


    for(SimplexId j = 0; j<toAdd[iScale].size(); j++){

      SimplexId filteredId = toAdd[iScale][j];
      SimplexId I = filteredToLocalIdVec[iScale][filteredId]; // localId
      SimplexId id = sampling_indexes[iScale][I]; // globalId
      points->SetPoint(i, coordinates[3*id+0], coordinates[3*id+1], coordinates[3*id+2]);


      // std::cout<<"         there"<<std::endl;
      normArray->SetTuple(i, normals+3*id);
      fitCoordsArray->SetTuple(i, fitCoordinates[iScale].data()+dim*I);
      mlsCoordsArray->SetTuple(i,mlsCoordinates[iScale].data()+3*I);
      mlsNormalsArray->SetTuple(i,mlsNormals[iScale].data()+3*I);
      if(compute3DEmbedding){
        double dum[3] = {-1,-1,-1};
        embedCoordsArray->SetTuple(i,dum);
      }
      // isStableArray->SetTuple1(i, isStable[id]);
      // densityArray->SetTuple(i, density+id);
      densityArray->SetTuple1(i, -1);
      idArray->SetTuple1(i,id);
      filteredIdArray->SetTuple1(i, filteredId);
      scaleIdArray->SetTuple1(i, iScale);
      epsIdArray->SetTuple1(i, 0);
      isEpsilonSampledArray->SetTuple1(i, 0);
      localIdArray->SetTuple1(i, I);
      i++;
    }
  }

		return 0;
}


void ttkRipsMultiParameterSampling::listAvailableFits(){
  
  std::string dataset_string = fittingStringId_.substr(0, fittingStringId_.find("_"));
  std::cout<<"Available saves: "<<std::endl;
    std::string path = savefile_path_;
    for (const auto & entry : fs::directory_iterator(path))
    {
      std::string data = entry.path().string();
      data = data.erase(0,savefile_path_.length()+std::string("savefile_").length());
      std::string dataset = data.substr(0, data.find("_"));

      if(dataset == dataset_string and data.find(".fit")!=data.npos){
        data = data.erase(0, data.find("_")+1);
        std::string fitting_params = data.substr(0, data.find(".fit"));
        std::cout<<entry.path()<<std::endl;
        
        int useSampling = readBoolFromString(fitting_params);
        int useNormalizedCoords = readBoolFromString(fitting_params);
        int changeBasis = readBoolFromString(fitting_params);
        int useKnnGraph = readBoolFromString(fitting_params);
        int useRIMLS = readBoolFromString(fitting_params);
        int useProjectiveNorm = readBoolFromString(fitting_params);

        int fitDimension = readIntFromString(fitting_params, 9);
        int nbPoints = readIntFromString(fitting_params, 9);
        int scaleSamplingCount = readIntFromString(fitting_params, 9);

        double scaleSamplingMin = readDoubleFromString(fitting_params, 4, 3);
        double scaleSamplingMax = readDoubleFromString(fitting_params, 4, 3);

        std::cout<<"    useSampling: "<<useSampling<<std::endl;
        std::cout<<"    useNormalizedCoords: "<<useNormalizedCoords<<std::endl;
        std::cout<<"    changeBasis: "<<changeBasis<<std::endl;
        std::cout<<"    useKnnGraph: "<<useKnnGraph<<std::endl;
        std::cout<<"    useRIMLS: "<<useRIMLS<<std::endl;
        std::cout<<"    useProjectiveNorm: "<<useProjectiveNorm<<std::endl;
        
        std::cout<<"    fitDimension: "<<fitDimension<<std::endl;
        std::cout<<"    nbPoints: "<<nbPoints<<std::endl;
        std::cout<<"    scaleSamplingCount: "<<scaleSamplingCount<<std::endl;

        std::cout<<"    scaleSamplingMin: "<<scaleSamplingMin<<std::endl;
        std::cout<<"    scaleSamplingMax: "<<scaleSamplingMax<<std::endl;
      }

    }
}

void ttkRipsMultiParameterSampling::createOutputPrimitiveSegmentation(const PrimitiveSegmentation& seg, vtkMultiBlockDataSet* mb, const float* coordinates, const float* normals, const std::vector<float>& pointCloudBarycenter){
  
  SimplexId nbPrimitives = seg.getNumberOfPrimitives();
  mb->SetNumberOfBlocks(nbPrimitives);

  for(int iPrim = 0; iPrim<nbPrimitives; iPrim++){
    vtkNew<vtkUnstructuredGrid> vtu;
    const Primitive& prim = seg.getPrimitive(iPrim);
    createOutputPrimitive(prim, vtu, iPrim, coordinates, normals, pointCloudBarycenter);
    mb->SetBlock(iPrim, vtu);
  }


  
}

void ttkRipsMultiParameterSampling::createOutputPrimitive(const Primitive& prim, vtkUnstructuredGrid* vtu, SimplexId primId, const float* coordinates, const float* normals,
    const std::vector<float>& pointCloudBarycenter){

  using V = Primitive::VectorType;

  SimplexId nbPoints = prim.getNumberOfPoints();

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nbPoints+1); // +1 because I add the barycenter
  vtu->SetPoints(points);

  auto pd = vtu->GetPointData();
  vtkNew<vtkFloatArray> normArray;
  normArray->SetName("Normals");
  normArray->SetNumberOfComponents(3);
  normArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(normArray);

  vtkNew<vtkIntArray> sizeArray;
  sizeArray->SetName("size");
  sizeArray->SetNumberOfComponents(1);
  sizeArray->SetNumberOfTuples(nbPoints+1);
  sizeArray->Fill(nbPoints);
  pd->AddArray(sizeArray);

  vtkNew<vtkIntArray> isBaryArray;
  isBaryArray->SetName("isBaryals");
  isBaryArray->SetNumberOfComponents(1);
  isBaryArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(isBaryArray);

  vtkNew<vtkFloatArray> meanNormArray;
  meanNormArray->SetName("meanNormals");
  meanNormArray->SetNumberOfComponents(3);
  meanNormArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(meanNormArray);

  vtkNew<vtkFloatArray> baryOffsetArray;
  baryOffsetArray->SetName("baryOffset");
  baryOffsetArray->SetNumberOfComponents(3);
  baryOffsetArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(baryOffsetArray);

  vtkNew<vtkFloatArray> primIdArray;
  primIdArray->SetName("primitiveId");
  primIdArray->SetNumberOfComponents(1);
  primIdArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(primIdArray);

  vtkNew<vtkFloatArray> chosenIdArray;
  chosenIdArray->SetName("chosenId");
  chosenIdArray->SetNumberOfComponents(1);
  chosenIdArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(chosenIdArray);

  vtkNew<vtkFloatArray> projCoordsArray;
  projCoordsArray->SetName("projCoords");
  projCoordsArray->SetNumberOfComponents(3);
  projCoordsArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(projCoordsArray);

  vtkNew<vtkFloatArray> projNormalsArray;
  projNormalsArray->SetName("projNormals");
  projNormalsArray->SetNumberOfComponents(3);
  projNormalsArray->SetNumberOfTuples(nbPoints+1);
  pd->AddArray(projNormalsArray);

  float bary[3] = {0,0,0};


  std::vector<float> meanNormal(3, 0);
  for(SimplexId i=0; i<nbPoints; i++){
    SimplexId globalId = prim.point_indices()[i];
    Primitive::VectorType p = {coordinates[3*globalId+0], coordinates[3*globalId+1], coordinates[3*globalId+2]};
    Primitive::VectorType p_proj, normal_proj;
    prim.project(p, p_proj, normal_proj);

    points->SetPoint(i, p.coeff(0), p.coeff(1), p.coeff(2));

    projCoordsArray->SetTuple3(i, p_proj.coeff(0), p_proj.coeff(1), p_proj.coeff(2));
    projNormalsArray->SetTuple3(i, normal_proj.coeff(0), normal_proj.coeff(1), normal_proj.coeff(2));

    normArray->SetTuple(i, normals+3*globalId);
    primIdArray->SetTuple1(i, primId);
    chosenIdArray->SetTuple1(i, primId);
    isBaryArray->SetTuple1(i, 0);

    for(int iDim=0; iDim<3; iDim++){
      meanNormal[iDim] += normals[3*globalId+iDim]/nbPoints;
      bary[iDim] += p[iDim]/nbPoints;
    }
  }

  std::vector<float> baryOffset(3, 0);
  for(int iDim=0; iDim<3; iDim++){
    baryOffset[iDim] = bary[iDim]-pointCloudBarycenter[iDim];
  }

  for(int iDim=0; iDim<3; iDim++){
    meanNormArray->FillComponent(iDim, meanNormal[iDim]);
    baryOffsetArray->FillComponent(iDim, baryOffset[iDim]);
  }

  // barycenter info
  points->SetPoint(nbPoints, bary);
  normArray->SetTuple(nbPoints, meanNormal.data());
  primIdArray->SetTuple1(nbPoints, primId);
  chosenIdArray->SetTuple1(nbPoints, primId);
  isBaryArray->SetTuple1(nbPoints, 1);

  vtu->Allocate(nbPoints);
  for(SimplexId iCell=0; iCell<points->GetNumberOfPoints(); iCell++){
      vtkNew<vtkIdList> idList;
      idList->InsertNextId(iCell);
      vtu->InsertNextCell(1, idList);
  }
}

int ttkRipsMultiParameterSampling::ProductDiagramToVTU(vtkUnstructuredGrid *vtu,
                 const ttk::DiagramType &diagram,
                 const std::vector<SimplexId>& flags,
                 const std::vector<SimplexId>& chosenPair,
                 const std::vector<std::pair<int, int>>& scalePairs,
                 const std::vector<float>& scales,
                 const std::vector<float>& scores) {

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();

  const auto n = flags.size();

  vtkNew<vtkFloatArray> scoreArray{};
  scoreArray->SetNumberOfComponents(1);
  scoreArray->SetNumberOfTuples(n);
  scoreArray->SetName("score");
  cd->AddArray(scoreArray);

  vtkNew<vtkFloatArray> persistenceArray{};
  persistenceArray->SetNumberOfComponents(1);
  persistenceArray->SetNumberOfTuples(n);
  persistenceArray->SetName("Persistence");
  cd->AddArray(persistenceArray);

  vtkNew<vtkFloatArray> birthArray{};
  birthArray->SetNumberOfComponents(1);
  birthArray->SetNumberOfTuples(n);
  birthArray->SetName("birth");
  cd->AddArray(birthArray);

  vtkNew<vtkFloatArray> deathArray{};
  deathArray->SetNumberOfComponents(1);
  deathArray->SetNumberOfTuples(n);
  deathArray->SetName("death");
  cd->AddArray(deathArray);

  vtkNew<vtkFloatArray> minScaleArray{};
  minScaleArray->SetNumberOfComponents(1);
  minScaleArray->SetNumberOfTuples(n);
  minScaleArray->SetName("minScale");
  cd->AddArray(minScaleArray);

  vtkNew<vtkFloatArray> maxScaleArray{};
  maxScaleArray->SetNumberOfComponents(1);
  maxScaleArray->SetNumberOfTuples(n);
  maxScaleArray->SetName("maxScale");
  cd->AddArray(maxScaleArray);

  vtkNew<vtkFloatArray> scalePersistenceArray{};
  scalePersistenceArray->SetNumberOfComponents(1);
  scalePersistenceArray->SetNumberOfTuples(n);
  scalePersistenceArray->SetName("scalePersistence");
  cd->AddArray(scalePersistenceArray);

  vtkNew<vtkFloatArray> productPersistenceArray{};
  productPersistenceArray->SetNumberOfComponents(1);
  productPersistenceArray->SetNumberOfTuples(n);
  productPersistenceArray->SetName("productPersistence");
  cd->AddArray(productPersistenceArray);

  vtkNew<vtkIntArray> chosenIdArray{};
  chosenIdArray->SetNumberOfComponents(1);
  chosenIdArray->SetNumberOfTuples(n);
  chosenIdArray->SetName("chosenId");
  cd->AddArray(chosenIdArray);

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(4*n);

  for(SimplexId iCC=0; iCC<n; iCC++){
    const auto flag = flags[iCC];

    float birth = (float)diagram[iCC].birth.sfValue;
    float death = (float)diagram[iCC].death.sfValue;
    float p0[3] = {birth,birth,0};
    float p1[3] = {birth,death,0};
    float p2[3] = {birth,death,0};
    float p3[3] = {birth,birth,0};

    float score = -1;
    float productPersistence = -1;
    float scalePersistence = -1;
    float persistence = death-birth;
    float minScale=-1, maxScale=-1;

    if(flag != -1){

      score = scores[iCC];
      int minScaleId = scalePairs[iCC].first;
      int maxScaleId = scalePairs[iCC].second;

      minScale = scales[minScaleId];
      maxScale = scales[maxScaleId];

      // set points
      p0[2] = log(minScale);
      p1[2] = log(minScale);
      p2[2] = log(maxScale);
      p3[2] = log(maxScale);

      scalePersistence = log(maxScale/minScale);
      productPersistence = persistence*scalePersistence;

    }

    points->InsertPoint(4*iCC+0, p0);
    points->InsertPoint(4*iCC+1, p1);
    points->InsertPoint(4*iCC+2, p2);
    points->InsertPoint(4*iCC+3, p3);

    vtkNew<vtkIdList> list;
    list->InsertNextId(4*iCC+0);
    list->InsertNextId(4*iCC+1);
    list->InsertNextId(4*iCC+2);
    list->InsertNextId(4*iCC+3);

    vtu->InsertNextCell(VTK_QUAD, list);


    persistenceArray->SetTuple1(iCC, persistence);
    scoreArray->SetTuple1(iCC, score);
    minScaleArray->SetTuple1(iCC, minScale);
    maxScaleArray->SetTuple1(iCC, maxScale);
    scalePersistenceArray->SetTuple1(iCC, scalePersistence);
    productPersistenceArray->SetTuple1(iCC, productPersistence);
    birthArray->SetTuple1(iCC, birth);
    deathArray->SetTuple1(iCC, death);
    chosenIdArray->SetTuple1(iCC, flag);
  }
  vtu->SetPoints(points);
  return 0;
}
