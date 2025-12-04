#include "DataTypes.h"
#include "PersistenceDiagramUtils.h"
#include <RipsTree.h>
#include <Sampling.h>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <ttkRipsTree.h>

#include <vtkInformation.h>

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

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <ttkPersistenceDiagramUtils.h>

using namespace ttk;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkRipsTree);

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
ttkRipsTree::ttkRipsTree() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkRipsTree::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkRipsTree::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkRipsTree::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkUnstructuredGrid *inputDataSet = vtkUnstructuredGrid::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  std::vector<double> persVec{};
  if(UseMultiPersRequest){
    parseStringMultiQuery(this->PersistenceRequest, persVec, SortPersistences, DecreasingPersistences, *this);
  }else{
    persVec.push_back(PersistenceRequestDouble);
  }

  auto max_persistence = getMaxPersistence(inputDataSet);

  if(UseRelativePersistence){
    getAbsoluteValues(persVec, max_persistence);
  }

  // printMsg("TEST RIPS TREE");
  // DiagramType diagram;
  // VTUToRipsDiagram(diagram, inputDataSet, *this);
  // buildRipsTree(diagram);
  // printNodeIds();
  // printSizes(persVec[0]);
  // printMsg("TEST RIPS TREE FINISHED");

  vtkMultiBlockDataSet *outputMB = vtkMultiBlockDataSet::GetData(outputVector, 0);
  outputMB->SetNumberOfBlocks(persVec.size());

  // vtkNew<vtkUnstructuredGrid> projDiag{};
  // ProjectDiagramInsideDomain(inputDataSet, projDiag, *this);

  // vtkNew<vtkCellDataToPointData> cellDataToPointData{};
  // cellDataToPointData->SetInputDataObject(projDiag);
  // cellDataToPointData->Update();

  // vtkNew<vtkMaskPoints> maskPoints{};
  // maskPoints->SetInputConnection(cellDataToPointData->GetOutputPort());
  // maskPoints->SetMaximumNumberOfPoints(projDiag->GetNumberOfCells());
  // maskPoints->SetGenerateVertices(true);
  // maskPoints->SetSingleVertexPerCell(true);
  // maskPoints->Update();



  for(int iBlock=0; iBlock<persVec.size(); iBlock++){

    // vtkNew<vtkThreshold> threshold{};
    // threshold->SetInputConnection(maskPoints->GetOutputPort());
    // threshold->SetInputArrayToProcess(0, 0, 0,
    //                                   vtkDataObject::FIELD_ASSOCIATION_CELLS,
    //                                   "CriticalType");
// #if VTK_VERSION_NUMBER < VTK_VERSION_CHECK(9, 2, 0)

// #else
    // threshold->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    // threshold->SetLowerThreshold(0);
    // threshold->SetUpperThreshold(0);
// #endif

    // threshold->Update();

    vtkNew<vtkUnstructuredGrid> outputBlock{};
    // outputBlock->DeepCopy(maskPoints->GetOutput());

    // computeCurrentSizes(outputBlock, persVec[iBlock], *this);

    // vtkNew<vtkDoubleArray> persThreshArray{};
    // persThreshArray->SetNumberOfComponents(1);
    // persThreshArray->SetNumberOfTuples(1);
    // persThreshArray->SetName("persThresh");
    // persThreshArray->SetTuple1(0,persVec[iBlock]);
    // outputBlock->GetFieldData()->AddArray(persThreshArray);

    // createOutput(outputBlock,  coordinates, normals);
    outputMB->SetBlock(iBlock,outputBlock);

    std::cout<<"final test"<<std::endl;
  }
    // vtkNew<vtkDoubleArray> persThreshArray{};
    // persThreshArray->SetNumberOfComponents(1);
    // persThreshArray->SetNumberOfTuples(1);
    // persThreshArray->SetName("persThresh");
    // persThreshArray->SetTuple1(0,0);
    // outputMB->GetFieldData()->AddArray(persThreshArray);
  


  // return success
  return 1;
}

void findThresholds(const RipsTree& tree, const int nbCC, std::vector<double>& persVec, SimplexId sizeMin){

  int foundCC = 0;
  persVec.clear();

  SimplexId nbPoints = tree.getNumberOfPoints();

  std::cout<<"find thresholds, size min "<<sizeMin<<std::endl;
  for(size_t i=0; i<nbPoints; i++){
    SimplexId size = tree.getFullSize(i);
    auto p = tree.getPersistence(i);
    if(size >= sizeMin){
      persVec.push_back(p);
      foundCC++;
    }
    if(foundCC==nbCC){
      break;
    }
  }
}

// SimplexId getCurrentSize(const vtkUnstructuredGrid* vtu, double persThresh){

// }

void computeCurrentSizes(vtkUnstructuredGrid* vtu, double persThresh, const ttk::Debug &dbg){

  dbg.printMsg("Computing current sizes");

  auto nbPair = vtu->GetNumberOfPoints();
  auto pd = vtu->GetPointData();


  vtkDataArray* componentCurrentSize = pd->GetArray("componentSizeAtThreshold");

  auto
    * componentCurrentParent = pd->GetArray("componentParentAtThreshold");

  vtkDataArray * currentParentSize = pd->GetArray("currentParentSize");

  vtkNew<vtkDoubleArray> currentParentPersistence{};
  currentParentPersistence->SetNumberOfTuples(nbPair);
  currentParentPersistence->SetNumberOfComponents(1);
  currentParentPersistence->SetName("currentParentPersistence");

  vtkNew<vtkFloatArray> meanNormalArray{};
  meanNormalArray->SetNumberOfComponents(3);
  meanNormalArray->SetNumberOfTuples(nbPair);
  meanNormalArray->DeepCopy(pd->GetArray("Normals"));
  meanNormalArray->SetName("meanNormals");
  pd->AddArray(meanNormalArray);

  vtkNew<vtkFloatArray> meanCoordinatesArray{};
  meanCoordinatesArray->SetNumberOfComponents(3);
  meanCoordinatesArray->SetNumberOfTuples(nbPair);
  meanCoordinatesArray->SetName("meanCoordinates");
  pd->AddArray(meanCoordinatesArray);

  vtkNew<vtkIntArray> isSelectedArray{};
  isSelectedArray->SetNumberOfComponents(1);
  isSelectedArray->SetNumberOfTuples(nbPair);
  isSelectedArray->SetName("isSelected");
  pd->AddArray(isSelectedArray);

  vtkNew<vtkIntArray> ccIdArray{};
  ccIdArray->SetNumberOfComponents(1);
  ccIdArray->SetNumberOfTuples(nbPair);
  ccIdArray->SetName("ccId");
  pd->AddArray(ccIdArray);

  float* meanNormals = (float*)ttkUtils::GetVoidPointer(meanNormalArray);
  float* normals = (float*)ttkUtils::GetVoidPointer(pd->GetArray("Normals") );
  float* meanCoords = (float*)ttkUtils::GetVoidPointer(meanCoordinatesArray);

  std::vector<int> sizesAtThreshold(nbPair, 1);

  for(size_t i =nbPair-1; 0<i; i--){
    auto p = pd->GetArray("Persistence")->GetTuple1(i);
    if(p<persThresh){
      SimplexId size = pd->GetArray("componentSize")->GetTuple1(i);
      SimplexId parent = pd->GetArray("parentComponent")->GetTuple1(i);
      sizesAtThreshold[i] = 0;
      sizesAtThreshold[parent] += size;
   }
  }


  SimplexId nbConnectedComponents=1;


  {
    // global pair is always here
    componentCurrentParent->SetTuple1(0,0);
    ccIdArray->SetTuple1(0,0);
    currentParentSize->SetTuple1(0,sizesAtThreshold[0]);
    currentParentPersistence->SetTuple1(0,pd->GetArray("Persistence")->GetTuple1(0));
    meanNormals[3*0+0]=0;
    meanNormals[3*0+1]=0;
    meanNormals[3*0+2]=0;
  }

  // can probably parallelized
  for(size_t i =1; i<nbPair; i++){
    auto p = pd->GetArray("Persistence")->GetTuple1(i);
    SimplexId size = pd->GetArray("componentSize")->GetTuple1(i);
    SimplexId parent = pd->GetArray("parentComponent")->GetTuple1(i);
    if(p>=persThresh){
      componentCurrentParent->SetTuple1(i,i);
      ccIdArray->SetTuple1(i,i);
      nbConnectedComponents++;
    }else{
      componentCurrentParent->SetTuple1(i,componentCurrentParent->GetTuple1(parent));
    }
    currentParentSize->SetTuple1(i,sizesAtThreshold[componentCurrentParent->GetTuple1(i)]);
    currentParentPersistence->SetTuple1(i,pd->GetArray("Persistence")->GetTuple1(componentCurrentParent->GetTuple1(i)));

    meanNormals[3*i+0]=0;
    meanNormals[3*i+1]=0;
    meanNormals[3*i+2]=0;
  }


  for(size_t i =0; i<nbPair; i++){
    SimplexId currentParent = componentCurrentParent->GetTuple1(i);
    // meanCoords[3*currentParent+0]+=meanCoords[3*i+0]/sizesAtThreshold[currentParent];
    // meanCoords[3*currentParent+1]+=meanCoords[3*i+1]/sizesAtThreshold[currentParent];
    // meanCoords[3*currentParent+2]+=meanCoords[3*i+2]/sizesAtThreshold[currentParent];
    meanNormals[3*currentParent+0]+=normals[3*i+0]/sizesAtThreshold[currentParent];
    meanNormals[3*currentParent+1]+=normals[3*i+1]/sizesAtThreshold[currentParent];
    meanNormals[3*currentParent+2]+=normals[3*i+2]/sizesAtThreshold[currentParent];

    ccIdArray->SetTuple1(i, ccIdArray->GetTuple1(currentParent));
    isSelectedArray->SetTuple1(i,1);
  }



// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < nbPair; ++i) {
    componentCurrentSize->SetTuple1(i, sizesAtThreshold[i]); // careful, not thread-safe
    SimplexId currentParent = componentCurrentParent->GetTuple1(i);
    meanNormals[3*i+0]=meanNormals[3*currentParent+0];
    meanNormals[3*i+1]=meanNormals[3*currentParent+1];
    meanNormals[3*i+2]=meanNormals[3*currentParent+2];

  }


  // void* currentParentPtr = ttkUtils::GetVoidPointer(componentCurrentParent);

  //compute bounds
  vtkDataArray* coordinatesArray = vtu->GetPoints()->GetData();
  float* coordinates = (float*) ttkUtils::GetVoidPointer(coordinatesArray);

  // array xMin, xMax, yMin, yMax, zMin, zMax
  std::array<float, 6> initBoundArray = {std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
                                         std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
                                         std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()};
  std::vector<std::array<float, 6>> coordsBounds(nbConnectedComponents, initBoundArray);
  std::vector<std::array<float, 6>> normBounds(nbConnectedComponents, initBoundArray);



  for(SimplexId i=0; i<nbPair; i++){
    SimplexId currentParent = componentCurrentParent->GetTuple1(i);
    auto& bounds = coordsBounds[currentParent];
    auto& nbounds = normBounds[currentParent];

    if(coordinates[3*i+0]<bounds[0])      {bounds[0] = coordinates[3*i+0];}
    if(coordinates[3*i+0]>bounds[1]) {bounds[1] = coordinates[3*i+0];}
    if(coordinates[3*i+1]<bounds[2])      {bounds[2] = coordinates[3*i+1];}
    if(coordinates[3*i+1]>bounds[3]) {bounds[3] = coordinates[3*i+1];}
    if(coordinates[3*i+2]<bounds[4])      {bounds[4] = coordinates[3*i+2];}
    if(coordinates[3*i+2]>bounds[5]) {bounds[5] = coordinates[3*i+2];}

    if(normals[3*i+0]<nbounds[0])      {nbounds[0] = normals[3*i+0];}
    if(normals[3*i+0]>nbounds[1]) {nbounds[1] = normals[3*i+0];}
    if(normals[3*i+1]<nbounds[2])      {nbounds[2] = normals[3*i+1];}
    if(normals[3*i+1]>nbounds[3]) {nbounds[3] = normals[3*i+1];}
    if(normals[3*i+2]<nbounds[4])      {nbounds[4] = normals[3*i+2];}
    if(normals[3*i+2]>nbounds[5]) {nbounds[5] = normals[3*i+2];}

  }

  vtkNew<vtkFloatArray> diagonalArray{};
  diagonalArray->SetNumberOfComponents(1);
  diagonalArray->SetNumberOfTuples(nbPair);
  diagonalArray->SetName("aabbDiagonalLength");
  pd->AddArray(diagonalArray);
  float* diagonal = (float*)ttkUtils::GetVoidPointer(diagonalArray);

  vtkNew<vtkFloatArray> diagonalNormalArray{};
  diagonalNormalArray->SetNumberOfComponents(1);
  diagonalNormalArray->SetNumberOfTuples(nbPair);
  diagonalNormalArray->SetName("aabbNormalDiagonalLength");
  pd->AddArray(diagonalNormalArray);
  float* normDiagonal = (float*)ttkUtils::GetVoidPointer(diagonalNormalArray);

  for(SimplexId i=0; i<nbConnectedComponents; i++){
    auto& bounds = coordsBounds[i];
    auto& nbounds = normBounds[i];
    diagonal[i] = std::sqrt( (bounds[1]-bounds[0])*(bounds[1]-bounds[0])
                  + (bounds[3]-bounds[2])*(bounds[3]-bounds[2])
                  + (bounds[5]-bounds[4])*(bounds[5]-bounds[4]) );
    normDiagonal[i] = std::sqrt( (nbounds[1]-nbounds[0])*(nbounds[1]-nbounds[0])
                  + (nbounds[3]-nbounds[2])*(nbounds[3]-nbounds[2])
                  + (nbounds[5]-nbounds[4])*(nbounds[5]-nbounds[4]) );
  }

  for(SimplexId i=0; i<nbPair; i++){
    SimplexId currentParent = componentCurrentParent->GetTuple1(i);
    diagonal[i] = diagonal[currentParent];
    normDiagonal[i] = normDiagonal[currentParent];

  }
    
  // add total nb of components as fielddata array
  vtkNew<vtkIntArray> nbCCArray{};
  nbCCArray->SetNumberOfComponents(1);
  nbCCArray->SetNumberOfTuples(1);
  nbCCArray->SetName("NumberOfComponents");
  nbCCArray->SetTuple1(0,nbConnectedComponents);
  vtu->GetFieldData()->AddArray(nbCCArray);


  pd->AddArray(currentParentPersistence);
}

double getMaxPersistence(vtkUnstructuredGrid* diag){

  auto persArray = diag->GetCellData()->GetArray("Persistence");
  if(persArray->GetNumberOfTuples()<2){
    return 0.0;
  }

  return persArray->GetTuple1(1);

}

void getAbsoluteValues(std::vector<double>& persistences, double max_persistence){
  for(int i=0; i<persistences.size(); i++) {
    persistences[i] = persistences[i]*max_persistence;
  }
}

void parseStringMultiQuery(std::string& query, std::vector<double>& output,
    bool sortOutput, bool descendingOrder,
                 const ttk::Debug &dbg){

    output.clear();

    // parse PersistenceRequest into a vector of double
    // code copied from ttkTriangulationRequest (thanks P. Guillou)
    std::istringstream iss(query);
    double v;
    bool start_range=false;
    bool end_range=false;
    double step;

    while(iss >> v){
      if(iss.peek()==':'){
        iss.ignore();
        double step, end;
        iss >> step;
        if(iss.peek()!=':'){
          dbg.printErr("Wrong formatting in range query");
          return;
        }
        iss.ignore();
        iss >> end;
        for(double w=v; w<=end; w+=step){
          output.emplace_back(w);
        }
      }else{
        output.emplace_back(v);
      }
      if(iss.peek() == ',' or iss.peek()==':') {
        iss.ignore();
      }
    }
    if(sortOutput){
      if(descendingOrder){
        std::sort(output.begin(), output.end(), std::greater<double>());
      }else{
        std::sort(output.begin(), output.end());
      }
    }

}


int VTUToRipsDiagram(ttk::DiagramType &diagram,
                 vtkUnstructuredGrid *vtu,
                 const ttk::Debug &dbg) {

  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();
  const auto points = vtu->GetPoints();

  if(pd == nullptr) {
    dbg.printErr("VTU diagram with NULL Point Data");
    return -1;
  }
  if(cd == nullptr) {
    dbg.printErr("VTU diagram with NULL Cell Data");
    return -2;
  }
  if(points == nullptr) {
    dbg.printErr("VTU with no points");
    return -3;
  }

  // cell data
  const auto pairId = vtkIntArray::SafeDownCast(
    cd->GetArray(ttk::PersistencePairIdentifierName));
  const auto pairType
    = vtkIntArray::SafeDownCast(cd->GetArray(ttk::PersistencePairTypeName));
  const auto pairPers = cd->GetArray(ttk::PersistenceName);
  const auto birthScalars = cd->GetArray(ttk::PersistenceBirthName);
  const auto isFinite = cd->GetArray(ttk::PersistenceIsFinite);
  const auto parent = vtkIntArray::SafeDownCast(cd->GetArray("parentComponent"));
  const auto pairSize = vtkIntArray::SafeDownCast(cd->GetArray("componentSize"));

  // point data
  const auto vertexId
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::VertexScalarFieldName));
  const auto critType
    = vtkIntArray::SafeDownCast(pd->GetArray(ttk::PersistenceCriticalTypeName));
  const auto coords = vtkFloatArray::SafeDownCast(
    pd->GetArray(ttk::PersistenceCoordinatesName));

  const bool embed = coords == nullptr;

  if(pairId == nullptr) {
    dbg.printErr("Missing PairIdentifier cell data array");
    return -5;
  }
  if(pairType == nullptr) {
    dbg.printErr("Missing PairType cell data array");
    return -6;
  }
  if(pairPers == nullptr) {
    dbg.printErr("Missing Persistence cell data array");
    return -7;
  }
  if(vertexId == nullptr) {
    dbg.printErr("Missing ttkVertexScalarField point data array");
    return -8;
  }
  if(critType == nullptr) {
    dbg.printErr("Missing CriticalType point data array");
    return -9;
  }
  if(birthScalars == nullptr) {
    dbg.printErr("Missing Birth cell data array");
    return -10;
  }
  if(isFinite == nullptr) {
    dbg.printErr("Missing IsFinite cell data array");
    return -12;
  }
  if(parent == nullptr) {
    dbg.printErr("Missing Parent cell data array");
    return -13;
  }
  if(pairSize == nullptr) {
    dbg.printErr("Missing pairSize cell data array");
    return -14;
  }

  int nPairs = pairId->GetNumberOfTuples();

  // compact pairIds in [0, nPairs - 1] (diagonal excepted)
  for(int i = 0; i < nPairs; i++) {
    if(pairId->GetTuple1(i) != -1) {
      pairId->SetTuple1(i, i);
    } else {
      // detect diagram diagonal
      nPairs -= 1;
    }
  }

  if(nPairs < 1) {
    dbg.printErr("Diagram has no pairs");
    return -4;
  }

  diagram.resize(nPairs);

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < nPairs; ++i) {

    const auto pId = pairId->GetValue(i);
    // skip diagram diagonal if present
    if(pId == -1) {
      continue;
    }

    const auto i0 = 2 * i + 0;
    const auto i1 = 2 * i + 1;

    const auto v0 = vertexId->GetValue(i0);
    const auto v1 = vertexId->GetValue(i1);
    const auto ct0 = static_cast<ttk::CriticalType>(critType->GetValue(i0));
    const auto ct1 = static_cast<ttk::CriticalType>(critType->GetValue(i1));

    const auto pType = pairType->GetValue(i);
    const auto pers = pairPers->GetTuple1(i);
    const auto size = pairSize->GetTuple1(i);
    const auto par = parent->GetTuple1(i);
    const auto birth = birthScalars->GetTuple1(i);
    const auto isFin = static_cast<bool>(isFinite->GetTuple1(i));

    std::array<float, 3> coordsBirth{}, coordsDeath{};
    // no vtkPoints::GetPoint() taking a float array, have to do the
    // conversion by hand...
    std::array<double, 3> tmp{};

    if(embed) {
      points->GetPoint(i0, tmp.data());
      coordsBirth = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
      points->GetPoint(i1, tmp.data());
      coordsDeath = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
    } else {
      coords->GetTuple(i0, tmp.data());
      coordsBirth = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
      coords->GetTuple(i1, tmp.data());
      coordsDeath = {static_cast<float>(tmp[0]), static_cast<float>(tmp[1]),
                     static_cast<float>(tmp[2])};
    }

    // put pairs in diagram
    diagram[i] = ttk::PersistencePair{
      ttk::CriticalVertex{v0, ct0, birth, coordsBirth},
      ttk::CriticalVertex{v1, ct1, birth + pers, coordsDeath}, pType, isFin, size, par};
  }

  return 0;
}

void createOutput(vtkUnstructuredGrid* vtu, const RipsTree& tree, float* coordinates, float* normals){

  SimplexId nbPoints = tree.getNumberOfPoints();
  auto pd = vtu->GetPointData();
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nbPoints);
  vtu->SetPoints(points);
  vtu->Allocate(nbPoints);

  vtkNew<vtkDoubleArray> persArray;
  persArray->SetName("Persistence");
  persArray->SetNumberOfComponents(1);
  persArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(persArray);

  vtkNew<vtkFloatArray> normArray;
  normArray->SetName("Normals");
  normArray->SetNumberOfComponents(3);
  normArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(normArray);

  vtkNew<ttkSimplexIdTypeArray> ccId;
  ccId->SetNumberOfComponents(1);
  ccId->SetNumberOfTuples(nbPoints);
  ccId->SetName("CCIdentifier");
  pd->AddArray(ccId);

  vtkNew<ttkSimplexIdTypeArray> fullSizeArray;
  fullSizeArray->SetNumberOfComponents(1);
  fullSizeArray->SetNumberOfTuples(nbPoints);
  fullSizeArray->SetName("fullSize");
  pd->AddArray(fullSizeArray);


  // create geometry
  for(int i=0; i<nbPoints; i++){
    SimplexId id = tree.getPointId(i);
    points->SetPoint(id, coordinates[3*id+0], coordinates[3*id+1], coordinates[3*id+2]);

    vtkNew<vtkIdList> idList;
    idList->InsertNextId(id);
    vtu->InsertNextCell(1, idList);

    normArray->SetTuple(id, normals+3*id);
    ccId->SetTuple1(id, i);
    
    double p = tree.getPersistence(i);
    persArray->SetTuple1(id,p);

    fullSizeArray->SetTuple1(id, tree.getSize(i));

  }
}

void createOutputConnectedComponent(vtkUnstructuredGrid* vtu, const RipsTree& tree, SimplexId ccId, double threshold, float* coordinates, float* normals, float* fitCoordinates, int* isStable, float* density, SimplexId sizeMin, std::vector<SimplexId>& sampling_indexes){

  std::vector<SimplexId> pointIds;
  tree.getComponentPointIdsAtThreshold(ccId, threshold, pointIds, sizeMin);

  std::vector<SimplexId> nodeIds;
  tree.getComponentNodeIdsAtThreshold(ccId, threshold, nodeIds, sizeMin);


  SimplexId nbPoints = pointIds.size();
  if(nodeIds.size()!=nbPoints){
    std::cout<<"MISMATCH SIZE "<<nodeIds.size()<<" vs "<<nbPoints<<std::endl;
  }

  auto pd = vtu->GetPointData();
  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nbPoints);
  vtu->SetPoints(points);
  vtu->Allocate(nbPoints);

  vtkNew<vtkFloatArray> normArray;
  normArray->SetName("Normals");
  normArray->SetNumberOfComponents(3);
  normArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(normArray);

  // vtkNew<vtkFloatArray> densityArray;
  // densityArray->SetName("density");
  // densityArray->SetNumberOfComponents(1);
  // densityArray->SetNumberOfTuples(nbPoints);
  // pd->AddArray(densityArray);

  vtkNew<vtkFloatArray> fitCoordsArray;
  fitCoordsArray->SetName("FitCoords");
  fitCoordsArray->SetNumberOfComponents(5);
  fitCoordsArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(fitCoordsArray);

  vtkNew<vtkIntArray> isStableArray;
  isStableArray->SetName("isFitStable");
  isStableArray->SetNumberOfComponents(1);
  isStableArray->SetNumberOfTuples(nbPoints);
  pd->AddArray(isStableArray);

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
    SimplexId id = sampling_indexes[pointIds[i]];
    points->SetPoint(i, coordinates[3*id+0], coordinates[3*id+1], coordinates[3*id+2]);

    vtkNew<vtkIdList> idList;
    idList->InsertNextId(i);
    vtu->InsertNextCell(1, idList);

    normArray->SetTuple(i, normals+3*id);
    fitCoordsArray->SetTuple(i, fitCoordinates+5*id);
    isStableArray->SetTuple1(i, isStable[id]);
    // densityArray->SetTuple(i, density+id);
    idArray->SetTuple1(i,id);
    filteredIdArray->SetTuple1(i,pointIds[i]);
    nodeIdArray->SetTuple1(i, nodeIds[i]);
    const auto& c = tree.getComponentSafe(nodeIds[i]);
    mergeIdArray->SetTuple1(i,c.startNodeId());
    fullsizeArray->SetTuple1(i, c.size());
    ccIdArray->SetTuple1(i, ccId);
    deathArray->SetTuple1(i, c.death());

  }
  // double b[6];
  // vtu->GetBounds(b);
  // std::cout<<"created components with bounds "<<b[0]<<" "<<b[1]<<" "<<b[2]<<" "<<b[3]<<" "<<b[4]<<" "<<b[5]<<std::endl;
}



