/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
/// \ingroup vtk
/// \class ttkRipsTree
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::RipsTree module.
///
/// This VTK filter uses the ttk::RipsTree module to compute an averaging of
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
///   - standalone/RipsTree/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::RipsTree
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkRipsTreeModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkUnstructuredGrid.h>

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
 *   ttkRipsTree
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <RipsTree.h>
#include <Sampling.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPCellDataToPointData.h>
#include <vtkMaskPoints.h>
#include <vtkThreshold.h>
#include <vtkPoints.h>

using namespace ttk;

class TTKRIPSTREE_EXPORT ttkRipsTree
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::RipsTree // and we inherit from the base class
{
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};
  std::string PersistenceRequest{};
  double PersistenceRequestDouble{0};
  bool UseRelativePersistence{false};
  bool UseMultiPersRequest{false};
  bool SortPersistences{true};
  bool DecreasingPersistences{false};

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */
  // vtkSetMacro(OutputArrayName, const std::string &);
  // vtkGetMacro(OutputArrayName, std::string);

  vtkSetMacro(PersistenceRequestDouble, double);
  vtkGetMacro(PersistenceRequestDouble, double);

  vtkSetMacro(SortPersistences, int);
  vtkGetMacro(SortPersistences, int);

  vtkSetMacro(DecreasingPersistences, int);
  vtkGetMacro(DecreasingPersistences, int);

  vtkSetMacro(UseMultiPersRequest, int);
  vtkGetMacro(UseMultiPersRequest, int);

  vtkSetMacro(UseRelativePersistence, int);
  vtkGetMacro(UseRelativePersistence, int);

  vtkSetMacro(PersistenceRequest, const std::string &);
  vtkGetMacro(PersistenceRequest, std::string);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkRipsTree *New();
  vtkTypeMacro(ttkRipsTree, ttkAlgorithm);


protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkRipsTree();
  ~ttkRipsTree() override = default;

  /**
   * TODO 8: Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   *         (see cpp file)
   * TODO 9: Specify the data object type of each output port
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;


};



// double AABBToDiag(std::vector<double>& aabb, int dim);
// void compute5DAABB(const std::vector<SimplexId>& local_ids, const std::vector<double>& fitCoords, const int dim, std::vector<double>& outputAABB);
// double compute5DAABBDiagonal(const std::vector<SimplexId>& local_ids, const std::vector<double>& fitCoords, const int dim);
// void compute5DAABBMultiscale(const std::vector<std::vector<SimplexId>>& local_ids, const std::vector<std::vector<double>>& fitCoordsVec, const int dim, std::vector<double>& outputAABB);
// double compute5DAABBDiagonalMultiscale(const std::vector<std::vector<SimplexId>>& local_ids, const std::vector<std::vector<double>>& fitCoordsVec, const int dim);


void createOutput(vtkUnstructuredGrid* vtu, const ttk::RipsTree& tree, float* coordinates, float* normals);
void createOutputConnectedComponent(vtkUnstructuredGrid* vtu, const RipsTree& tree, SimplexId ccId, double threshold, float* coordinates, float* normals, float* fitCoordinates, int* isStable, float* density, SimplexId sizeMin, std::vector<SimplexId>&);
double getMaxPersistence(vtkUnstructuredGrid* diag);
void getAbsoluteValues(std::vector<double>& persistences, double max_persistence);
void computeCurrentSizes(vtkUnstructuredGrid* vtu, double persThresh, const ttk::Debug &dbg);
void parseStringMultiQuery(std::string& query, std::vector<double>& output,
    bool sortOutput=false, bool descendingOrder=false, const ttk::Debug &dbg=Debug());
int VTUToRipsDiagram(ttk::DiagramType &diagram, vtkUnstructuredGrid *vtu, const ttk::Debug &dbg);

bool createOutputConnectedComponentRips(vtkUnstructuredGrid* vtu,
                                              const RipsTree& tree,
                                              SimplexId ccId,
                                              double threshold,
                                              float* coordinates,
                                              std::vector<double>& density){




  std::vector<SimplexId> pointIds;
  tree.getComponentPointIdsAtThreshold(ccId, threshold, pointIds, 0);

  std::vector<SimplexId> nodeIds;
  tree.getComponentNodeIdsAtThreshold(ccId, threshold, nodeIds, 0);


  SimplexId nbPoints = pointIds.size();
  if(nodeIds.size()!=nbPoints){
    std::cout<<"MISMATCH SIZE "<<nodeIds.size()<<" vs "<<nbPoints<<std::endl;
  }

  auto pd = vtu->GetPointData();
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nbPoints);
  vtu->SetPoints(points);


  vtkNew<vtkDoubleArray> densityArray;
  densityArray->SetName("density");
  densityArray->SetNumberOfComponents(1);
  densityArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(densityArray);




  vtkNew<vtkIntArray> idArray;
  idArray->SetName("originalIdx");
  idArray->SetNumberOfComponents(1);
  idArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(idArray);



  vtkNew<vtkIntArray> mergeIdArray;
  mergeIdArray->SetName("mergeId");
  mergeIdArray->SetNumberOfComponents(1);
  mergeIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(mergeIdArray);

  vtkNew<vtkIntArray> nodeIdArray;
  nodeIdArray->SetName("nodeId");
  nodeIdArray->SetNumberOfComponents(1);
  nodeIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(nodeIdArray);

  vtkNew<vtkIntArray> fullsizeArray;
  fullsizeArray->SetName("fullsize");
  fullsizeArray->SetNumberOfComponents(1);
  fullsizeArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(fullsizeArray);

  vtkNew<vtkIntArray> ccIdArray;
  ccIdArray->SetName("ccId");
  ccIdArray->SetNumberOfComponents(1);
  ccIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(ccIdArray);

  vtkNew<vtkFloatArray> deathArray;
  deathArray->SetName("death");
  deathArray->SetNumberOfComponents(1);
  deathArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(deathArray);


  // create geometry
  for(int i=0; i<nbPoints; i++){

    SimplexId id = pointIds[i];
    points->SetPoint(i, coordinates[3*id+0], coordinates[3*id+1], coordinates[3*id+2]);
    
    densityArray->SetTuple1(i, density[id]);
    idArray->SetTuple1(i,id);
    nodeIdArray->SetTuple1(i, nodeIds[i]);
    const auto& c = tree.getComponentSafe(nodeIds[i]);
    mergeIdArray->SetTuple1(i,c.startNodeId());
    fullsizeArray->SetTuple1(i, c.size());
    ccIdArray->SetTuple1(i, ccId);
    deathArray->SetTuple1(i, c.death());

  }


  vtu->Allocate(nbPoints);
  for(SimplexId iCell=0; iCell<points->GetNumberOfPoints(); iCell++){
      vtkNew<vtkIdList> idList;
      idList->InsertNextId(iCell);
      vtu->InsertNextCell(1, idList);
  }


  bool discard = 0;

  double b[6];
  vtu->GetBounds(b);

  return discard;
}
template<class PointND>
bool createOutputConnectedComponentMultiScale(vtkUnstructuredGrid* vtu,
                                              const RipsTree& tree,
                                              SimplexId ccId,
                                              double threshold,
                                              float* coordinates,
                                              float* normals,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              std::vector<std::vector<float>>& mlsCoordinates,
                                              std::vector<std::vector<float>>& mlsNormals,
                                              std::vector<std::vector<double>>& embeddingCoordinates,
                                              std::vector<std::vector<int>>& isStable,
                                              std::vector<double>& density,
                                              SimplexId sizeMin,
                                              std::vector<SimplexId>& highestPersComponent,
                                              std::vector<std::vector<SimplexId>>& sampling_indexes,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalIdVec,
                                              std::vector<SimplexId>& pointsAtScaleStartIdx,
                                              std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId){


  SimplexId nbScales = sampling_indexes.size();
  if(epsilonIdToFilteredId.size()!=nbScales or filteredToLocalIdVec.size()!=nbScales or epsilonSampler.size()!=nbScales or pointsAtScaleStartIdx.size()!=nbScales){
    std::cout<<"SIZE DOES NOT MATCH SCALE NUMBER"<<std::endl;
  }


  int dim = PointND::Dim;


  bool compute3DEmbedding = embeddingCoordinates.size()==nbScales;
  // for(int iScale=0; iScale<nbScales; iScale++){
  //   std::cout<<" filtered size "<<filteredToLocalIdVec[iScale].size()<<std::endl;
  // }

  std::vector<SimplexId> pointIds;
  tree.getComponentPointIdsAtThreshold(ccId, threshold, pointIds, sizeMin);

  std::vector<SimplexId> nodeIds;
  tree.getComponentNodeIdsAtThreshold(ccId, threshold, nodeIds, sizeMin);


  SimplexId nbPoints = pointIds.size();
  if(nodeIds.size()!=nbPoints){
    std::cout<<"MISMATCH SIZE "<<nodeIds.size()<<" vs "<<nbPoints<<std::endl;
  }
    // std::cout<<"node and points ids SIZE "<<nodeIds.size()<<" vs "<<nbPoints<<std::endl;

  auto pd = vtu->GetPointData();
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nbPoints);
  vtu->SetPoints(points);

  vtkNew<vtkFloatArray> normArray;
  normArray->SetName("Normals");
  normArray->SetNumberOfComponents(3);
  normArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(normArray);

  vtkNew<vtkDoubleArray> densityArray;
  densityArray->SetName("density");
  densityArray->SetNumberOfComponents(1);
  densityArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(densityArray);


  vtkNew<vtkDoubleArray> fitCoordsArray;
  fitCoordsArray->SetName("FitCoords");
  fitCoordsArray->SetNumberOfComponents(PointND::Dim);
  fitCoordsArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(fitCoordsArray);

  // vtkNew<vtkIntArray> isStableArray;
  // isStableArray->SetName("isFitStable");
  // isStableArray->SetNumberOfComponents(1);
  // isStableArray->SetNumberOfTuples(nbPoints);
  // pd->AddArray(isStableArray);

  vtkNew<vtkIntArray> filteredIdArray;
  filteredIdArray->SetName("filteredIdx");
  filteredIdArray->SetNumberOfComponents(1);
  filteredIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(filteredIdArray);

  vtkNew<vtkIntArray> idArray;
  idArray->SetName("originalIdx");
  idArray->SetNumberOfComponents(1);
  idArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(idArray);

  vtkNew<vtkIntArray> multiscaleIdArray;
  multiscaleIdArray->SetName("multiscaleId");
  multiscaleIdArray->SetNumberOfComponents(1);
  multiscaleIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(multiscaleIdArray);


  vtkNew<vtkIntArray> hghestPersCompArray;
  hghestPersCompArray->SetName("hghestPersComp");
  hghestPersCompArray->SetNumberOfComponents(1);
  hghestPersCompArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(hghestPersCompArray);

  vtkNew<vtkIntArray> mergeIdArray;
  mergeIdArray->SetName("mergeId");
  mergeIdArray->SetNumberOfComponents(1);
  mergeIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(mergeIdArray);

  vtkNew<vtkIntArray> nodeIdArray;
  nodeIdArray->SetName("nodeId");
  nodeIdArray->SetNumberOfComponents(1);
  nodeIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(nodeIdArray);

  vtkNew<vtkIntArray> fullsizeArray;
  fullsizeArray->SetName("fullsize");
  fullsizeArray->SetNumberOfComponents(1);
  fullsizeArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(fullsizeArray);

  vtkNew<vtkIntArray> ccIdArray;
  ccIdArray->SetName("ccId");
  ccIdArray->SetNumberOfComponents(1);
  ccIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(ccIdArray);

  vtkNew<vtkIntArray> scaleIdArray;
  scaleIdArray->SetName("scaleId");
  scaleIdArray->SetNumberOfComponents(1);
  scaleIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(scaleIdArray);

  vtkNew<vtkFloatArray> deathArray;
  deathArray->SetName("death");
  deathArray->SetNumberOfComponents(1);
  deathArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(deathArray);

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

  vtkNew<vtkDoubleArray> embedCoordsArray;
  if(compute3DEmbedding){
    embedCoordsArray->SetName("3d-Embedding");
    embedCoordsArray->SetNumberOfComponents(3);
    embedCoordsArray->SetNumberOfTuples(nbPoints);
    pd->AddArray(embedCoordsArray);
  }

  vtkNew<vtkIntArray> epsilonIdArray;
  epsilonIdArray->SetName("isEpsilonSampled");
  epsilonIdArray->SetNumberOfComponents(1);
  epsilonIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(epsilonIdArray);

  vtkNew<vtkIntArray> localIdArray;
  localIdArray->SetName("localId");
  localIdArray->SetNumberOfComponents(1);
  localIdArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(localIdArray);

  std::vector<std::vector<SimplexId>> toAdd(nbScales);
  std::vector<std::vector<SimplexId>> localIds(nbScales);
  // create geometry
  for(int i=0; i<nbPoints; i++){
    int iScale = findScaleId(pointIds[i], pointsAtScaleStartIdx);
    SimplexId pointId = pointIds[i]-pointsAtScaleStartIdx[iScale]; // epsilon id
    SimplexId filteredId = epsilonIdToFilteredId[iScale][pointId];


    SimplexId I = filteredToLocalIdVec[iScale][filteredId];

    localIds[iScale].push_back(I);

    SimplexId id = sampling_indexes[iScale][I];
    points->SetPoint(i, coordinates[3*id+0], coordinates[3*id+1], coordinates[3*id+2]);
    

    normArray->SetTuple(i, normals+3*id);
    fitCoordsArray->SetTuple(i, fitCoordinates[iScale].data()+dim*I);
    mlsCoordsArray->SetTuple(i,mlsCoordinates[iScale].data()+3*I);
    mlsNormalsArray->SetTuple(i,mlsNormals[iScale].data()+3*I);
    if(compute3DEmbedding){
      embedCoordsArray->SetTuple(i,embeddingCoordinates[iScale].data()+3*I);
    }

    // isStableArray->SetTuple1(i, isStable[id]);
    densityArray->SetTuple1(i, density[pointIds[i]]);
    idArray->SetTuple1(i,id);
    multiscaleIdArray->SetTuple1(i,pointIds[i]);
    filteredIdArray->SetTuple1(i,filteredId); 
    nodeIdArray->SetTuple1(i, nodeIds[i]);
    const auto& c = tree.getComponentSafe(nodeIds[i]);
    mergeIdArray->SetTuple1(i,c.startNodeId());
    hghestPersCompArray->SetTuple1(i, highestPersComponent[id]);
    fullsizeArray->SetTuple1(i, c.size());
    ccIdArray->SetTuple1(i, ccId);
    scaleIdArray->SetTuple1(i, iScale);
    deathArray->SetTuple1(i, c.death());
    epsilonIdArray->SetTuple1(i, 1);
    localIdArray->SetTuple1(i, I);

    // find children in epsilon sampling
    auto& children = epsilonSampler[iScale].getChildren(pointId);
    if(filteredId != children[0]){
      std::cout<<"ERROR rep "<<filteredId<<" of scale "<<iScale<<" is not its own child"<<std::endl;
      for(auto c : children){
        std::cout<<" "<<c;
      }
      std::cout<<std::endl;
      // std::cout<<" BUT... eps array of "<<children[0]<<" is actually "<<filteredIdxToEpsSampling[iScale][children[0]]<<std::endl;
    }
    // std::cout<<" size of filtered->local "<<filteredToLocalIdVec[iScale].size()<<std::endl;
    for(SimplexId iChild=1; iChild<children.size(); iChild++){
      auto child = children[iChild];
      toAdd[iScale].push_back(child);
      // localIds[iScale].push_back(filteredToLocalIdVec[iScale][child]);
      // std::cout<<" to add at scale "<<iScale<<": point (filtered)"<<child<<" , (local)"<<filteredToLocalIdVec[iScale][child]<<std::endl;
    }

  }


  int d=dim;
  std::vector<double> aabb5D(2*d);
  for(int k=0; k<d; k++){
    fitCoordsArray->GetFiniteRange(aabb5D.data()+2*k, k+dim-d);
  }
  std::cout<<" aabb 1: "<<std::endl;
  for(int ii=0; ii<5; ii++){
    std::cout<<" ["<<aabb5D[2*ii]<<", "<<aabb5D[2*ii+1]<<"] ";
  }
  std::cout<<std::endl;
  double diag = AABBToDiag(aabb5D, d);
  std::cout<<" diag "<<diag<<std::endl;


  double diag2{};
  // std::vector<double> aabb5D2(2*2);
  // compute5DAABBMultiscale(localIds, fitCoordinates, dim, aabb5D2);
  // std::cout<<" aabb 2: "<<std::endl;
  // for(int ii=0; ii<5; ii++){
  //   std::cout<<" ["<<aabb5D2[2*ii]<<", "<<aabb5D2[2*ii+1]<<"] ";
  // }
  // std::cout<<std::endl;
  diag2 = compute5DAABBDiagonalMultiscale(localIds, fitCoordinates, dim);

  // std::cout<<" 5d aabb diag "<< diag<<"  vs  "<< diag2<<std::endl;



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
  mergeIdArray->Resize(nbTotalPoints);
  hghestPersCompArray->Resize(nbTotalPoints);
  multiscaleIdArray->Resize(nbTotalPoints);
  nodeIdArray->Resize(nbTotalPoints);
  fullsizeArray->Resize(nbTotalPoints);
  ccIdArray->Resize(nbTotalPoints);
  scaleIdArray->Resize(nbTotalPoints);
  deathArray->Resize(nbTotalPoints);
  mlsCoordsArray->Resize(nbTotalPoints);
  mlsNormalsArray->Resize(nbTotalPoints);
  epsilonIdArray->Resize(nbTotalPoints);
  localIdArray->Resize(nbTotalPoints);
  densityArray->Resize(nbTotalPoints);

  points->SetNumberOfPoints(nbTotalPoints);
  normArray->SetNumberOfTuples(nbTotalPoints);
  fitCoordsArray->SetNumberOfTuples(nbTotalPoints);
  filteredIdArray->SetNumberOfTuples(nbTotalPoints);
  idArray->SetNumberOfTuples(nbTotalPoints);
  mergeIdArray->SetNumberOfTuples(nbTotalPoints);
  hghestPersCompArray->SetNumberOfTuples(nbTotalPoints);
  multiscaleIdArray->SetNumberOfTuples(nbTotalPoints);
  nodeIdArray->SetNumberOfTuples(nbTotalPoints);
  fullsizeArray->SetNumberOfTuples(nbTotalPoints);
  ccIdArray->SetNumberOfTuples(nbTotalPoints);
  scaleIdArray->SetNumberOfTuples(nbTotalPoints);
  deathArray->SetNumberOfTuples(nbTotalPoints);
  mlsCoordsArray->SetNumberOfTuples(nbTotalPoints);
  mlsNormalsArray->SetNumberOfTuples(nbTotalPoints);
  epsilonIdArray->SetNumberOfTuples(nbTotalPoints);
  localIdArray->SetNumberOfTuples(nbTotalPoints);
  densityArray->SetNumberOfTuples(nbTotalPoints);

  if(compute3DEmbedding){
    embedCoordsArray->Resize(nbTotalPoints);
    embedCoordsArray->SetNumberOfTuples(nbTotalPoints);
  }
  
  // std::cout<<" Add extra points "<<nbExtraPoints<<", total "<<nbTotalPoints<<std::endl;
  SimplexId i = nbPoints;
  int minScale = nbScales;
  int maxScale = 0;
  for(int iScale=0; iScale<nbScales; iScale++){

    if(localIds[iScale].size()>0){
      if(iScale>maxScale) maxScale = iScale;
      if(iScale<minScale) minScale = iScale;
    }
    // std::cout<<"     scale "<<iScale<<": adding "<<toAdd[iScale].size()<<" points"<<std::endl;

    for(SimplexId j = 0; j<toAdd[iScale].size(); j++){

      // std::cout<<"     i, j, "<<i<<" "<<j<<std::endl;
      SimplexId filteredId = toAdd[iScale][j];
      SimplexId I = filteredToLocalIdVec[iScale][filteredId]; // localId
      // std::cout<<"         I="<<I<<" outta "<<sampling_indexes[iScale].size()<<std::endl;
      SimplexId id = sampling_indexes[iScale][I]; // globalId
      // std::cout<<"         id="<<id<<std::endl;
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
      nodeIdArray->SetTuple1(i, -1);
      mergeIdArray->SetTuple1(i,-1);
      hghestPersCompArray->SetTuple1(i, highestPersComponent[id]);
      multiscaleIdArray->SetTuple1(i,-1);
      fullsizeArray->SetTuple1(i, 1);
      ccIdArray->SetTuple1(i, ccId);
      scaleIdArray->SetTuple1(i, iScale);
      deathArray->SetTuple1(i, epsilonSampler[iScale].epsilon());
      epsilonIdArray->SetTuple1(i, 0);
      localIdArray->SetTuple1(i, I);
      i++;
    }
  }
  // std::cout<<"now "<<points->GetNumberOfPoints()<<" points"<<std::endl;
  vtu->Allocate(nbTotalPoints);
  for(SimplexId iCell=0; iCell<points->GetNumberOfPoints(); iCell++){
      vtkNew<vtkIdList> idList;
      idList->InsertNextId(iCell);
      vtu->InsertNextCell(1, idList);
  }


  vtkNew<vtkDoubleArray> aabbDiagArray;
  aabbDiagArray->SetName("5DaabbDiag");
  aabbDiagArray->SetNumberOfComponents(1);
  aabbDiagArray->SetNumberOfTuples(nbTotalPoints);
  aabbDiagArray->Fill(diag);
  pd->AddArray(aabbDiagArray);

  vtkNew<vtkDoubleArray> aabbDiagArray2;
  aabbDiagArray2->SetName("5DaabbDiag2");
  aabbDiagArray2->SetNumberOfComponents(1);
  aabbDiagArray2->SetNumberOfTuples(nbTotalPoints);
  aabbDiagArray2->Fill(diag2);
  pd->AddArray(aabbDiagArray2);

  vtkNew<vtkIntArray> maxScaleArray;
  maxScaleArray->SetName("maxScale");
  maxScaleArray->SetNumberOfComponents(1);
  maxScaleArray->SetNumberOfTuples(nbTotalPoints);
  maxScaleArray->Fill(maxScale);
  pd->AddArray(maxScaleArray);

  vtkNew<vtkIntArray> minScaleArray;
  minScaleArray->SetName("minScale");
  minScaleArray->SetNumberOfComponents(1);
  minScaleArray->SetNumberOfTuples(nbTotalPoints);
  minScaleArray->Fill(minScale);
  pd->AddArray(minScaleArray);

  bool discard = 0;

  double b[6];
  vtu->GetBounds(b);
  // std::cout<<"created components with bounds "<<b[0]<<" "<<b[1]<<" "<<b[2]<<" "<<b[3]<<" "<<b[4]<<" "<<b[5]<<std::endl;

  return discard;
}

void findThresholds(const RipsTree& tree, const int nbCC, std::vector<double>& persVec, SimplexId sizeMin);
