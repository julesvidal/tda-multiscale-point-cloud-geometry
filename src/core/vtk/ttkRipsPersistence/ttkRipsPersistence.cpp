#include <ttkRipsPersistence.h>

#include <utility>
#include <vtkIdList.h>
#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <PCAnalysisUtils.h>
#include <PCAnalysis.h>
#include <Sampling.h>
#include <RipsTree.h>
#include "../../base/ripsPersistence/nanoflann.hpp"
#include "DataTypes.h"
#include "PersistenceDiagramUtils.h"
#include "RipsPersistence.h"
#include <ttkPersistenceDiagramUtils.h>
#include <ttkRipsTree.h>

using namespace ttk;

// template <typename Derived>


// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkRipsPersistence);

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
ttkRipsPersistence::ttkRipsPersistence() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkRipsPersistence::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkRipsPersistence::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  if(port == 1) {
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
int ttkRipsPersistence::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  this->clear();
  vtkPointSet *inputDataSet = vtkUnstructuredGrid::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  vtkDataArray *inputField = this->GetInputArrayToProcess(0, inputVector);
  if(!inputField){
    printWrn("No Normal field was found -- Field won't be used");
    return 0;
  }else{
    printMsg("Using normal field: "+std::string(inputField->GetName()));
  }
  vtkDataArray* inputCoordinates = inputDataSet->GetPoints()->GetData();
  float* coordinates = (float*) ttkUtils::GetVoidPointer(inputCoordinates);

  int nbPoints = inputDataSet->GetNumberOfPoints();

  std::vector<double> births(nbPoints);
  for(SimplexId iPoint=0; iPoint<nbPoints; iPoint++){
    births[iPoint] = inputField->GetTuple1(iPoint);
  }

  SimplexId nEdges = inputDataSet->GetNumberOfCells();
  std::vector<edgeTuple> edges{};

  for(SimplexId iEdge=0; iEdge<nEdges; iEdge++){
    SimplexId v0 = inputDataSet->GetCell(iEdge)->GetPointId(0);
    SimplexId v1 = inputDataSet->GetCell(iEdge)->GetPointId(1);
    if(v0>v1){
      std::swap(v0,v1);
    }
    double f0 = inputField->GetTuple1(v0);
    double f1 = inputField->GetTuple1(v1);
    edges.emplace_back(v0,v1,std::max(f0,f1));
  }

  std::cout<<edges.size()<<" edges "<<std::endl;

  this->setEdges(edges);
std:cout<<"sort edges"<<std::endl;
  this->sortEdgesAscending();

  std::vector<edgeTuple> mergingEdges{};
  std::cout<<"diagram"<<std::endl;
  this->computePersistenceDiagramWithBirth(MaxEdgeLength, nbPoints, mergingEdges, births);
std::cout<<"sort diagram"<<std::endl;
  this->sortPersistenceDiagram(nbPoints);

  std::cout<<"Rips tree"<<std::endl;
  RipsTree tree{};
  tree.buildRipsTree(this->getDiagram());
  tree.buildRipsArcsWithBirth(births);

    // PERSISTENCE DIAGRAM OUTPUT
    {
      vtkUnstructuredGrid *outputPD = vtkUnstructuredGrid::GetData(outputVector,0);

      vtkNew<vtkDoubleArray> scalars{};
      scalars->SetNumberOfComponents(1);
      DiagramToVTU(outputPD, this->getDiagram(), scalars, *this, 3, 0);
    }

    // CLUSTERING OUTPUT
    {
      printMsg("Rips Tree ");

      vtkMultiBlockDataSet* outputMB = vtkMultiBlockDataSet::GetData(outputVector,1);

      // vtkNew<vtkFloatArray> barycenterFDArray{};
      // barycenterFDArray->SetNumberOfComponents(3);
      // barycenterFDArray->SetNumberOfTuples(1);
      // barycenterFDArray->SetTuple(0, resultBarycenter.data());
      // barycenterFDArray->SetName("basisCenter");
      // outputMB->GetFieldData()->AddArray(barycenterFDArray);


      std::vector<RipsArc> ripsArcs = std::vector<RipsArc>(tree.getRipsArcs());

      SimplexId nbOutputPrimitives = NbPrimitives;
      if(nbOutputPrimitives>ripsArcs.size()){
        nbOutputPrimitives = ripsArcs.size();
      }

      outputMB->SetNumberOfBlocks(nbOutputPrimitives);


      bool branchPersistence{true};
      if(branchPersistence){
        const auto cmp_persistence = [&](const RipsArc &a0, const RipsArc &a1){
          return  a0.ccId<a1.ccId or (a0.ccId==a1.ccId and (a0.birth<a1.birth or (a0.birth==a1.birth and a0.birthId < a1.birthId))) ; 
        };

        for(int iArc=0; iArc<ripsArcs.size();iArc++){
          auto& arc = ripsArcs[iArc];
          arc.death = tree.getComponent(arc.ccId).death();
        }
        TTK_PSORT(this->threadNumber_, ripsArcs.begin(), ripsArcs.end(), cmp_persistence);

        for(SimplexId iArc=0; iArc<ripsArcs.size(); iArc++){
          const auto& arc = ripsArcs[iArc];
          auto ccId = arc.ccId;
          auto birthId = arc.birthId;
          auto deathId = arc.deathId;
          auto birth = arc.birth;
          auto death = arc.death;
          auto type = arc.type;
          auto persistence = death-birth;
          auto parent = tree.getComponent(ccId).parent();
          // std::cout<<"  iArc    "<<iArc<<std::endl;
          // std::cout<<"    | "<<"ccId "<<ccId<<std::endl;
          // std::cout<<"    | "<<"birthId "<<birthId<<std::endl;
          // std::cout<<"    | "<<"deathId "<<deathId<<std::endl;
          // std::cout<<"    | "<<"birth "<<birth<<std::endl;
          // std::cout<<"    | "<<"death "<<death<<std::endl;
          // std::cout<<"    | "<<"persistence "<<persistence<<std::endl;
          // std::cout<<"    | "<<"type "<<type<<std::endl;
          // std::cout<<"    | "<<"parent "<<parent<<std::endl;
        }
      }

      double current_threshold = ripsArcs.size()>0 ? ripsArcs[0].death : 0.0;

      std::vector<SimplexId> ccIdToArc(tree.getNumberOfPoints(), -1);
      std::vector<double> thresholds(nbOutputPrimitives, -1);
      std::vector<double> sizesAtThreshold(nbOutputPrimitives, 0);
      std::vector<SimplexId> chosenArcs(nbOutputPrimitives, -1);
      std::vector<SimplexId> highestPersComponent(nbPoints,-1);

      SimplexId arcCount=0;
      for(SimplexId iArc=0; iArc < ripsArcs.size(); iArc++){
        if(arcCount>=nbOutputPrimitives){
          break;
        }
        auto type = ripsArcs[iArc].type;
        if(branchPersistence){
          if(type!=0){
            continue;
          }
        }
        std::cout<<"arcCount "<<arcCount<<std::endl;
        auto ccId = ripsArcs[iArc].ccId;
        ccIdToArc[ccId] = arcCount;
        // std::cout<<"set ccId "<<ccId<<" to arc "<<arcCount<<std::endl;
        thresholds[arcCount] = ripsArcs[iArc].death;
        chosenArcs[arcCount] = iArc;
        // std::cout<<"     selected arc "<<iArc<<" of ccid "<<ccId<<std::endl;
        std::cout<<"     set "<<arcCount<<" thresh to "<<thresholds[arcCount]<<std::endl;
        SimplexId ccId_parent = tree.getComponent(ccId).parent();
        SimplexId arcCount_parent = ccIdToArc[ccId_parent];
        // std::cout<<"     parent "<<ccId_parent<<", arcCount "<<arcCount_parent<<std::endl;

        thresholds[arcCount_parent] = std::min(ripsArcs[iArc].death, thresholds[arcCount_parent]);

        std::cout<<"     set "<<arcCount_parent<<" thresh to "<<thresholds[arcCount_parent]<<std::endl;

        arcCount++;
      }


      arcCount = 0;
      for(SimplexId iArc=0; iArc < ripsArcs.size(); iArc++){

        if(arcCount>=nbOutputPrimitives){
          break;
        }

        auto ccId = ripsArcs[iArc].ccId;

        // auto threshold = ripsArcs[iArc].death;
        auto threshold = thresholds[arcCount];
        auto birth = ripsArcs[iArc].birth;
        auto persistence = ripsArcs[iArc].death-ripsArcs[iArc].birth;
        auto type = ripsArcs[iArc].type;
        auto diag = ripsArcs[iArc].diag;

        if(branchPersistence){
          if(type!=0){
            continue;
          }
        }
        std::cout<<" arc "<<arcCount<<"/"<<nbOutputPrimitives-1<<std::endl;
        vtkNew<vtkUnstructuredGrid> block;

        bool discard = createOutputConnectedComponentRips(block, tree,
            ccId,
            threshold,
            coordinates,
            births);


        // int dummy = -1;
        // setSingleComponentAttributes(block, iArc, persistence, arcCount, type, -1,-1, -1);
        vtkNew<vtkIntArray> arcIdArray{};
        arcIdArray->SetNumberOfComponents(1);
        arcIdArray->SetNumberOfTuples(block->GetNumberOfPoints());
        arcIdArray->SetName("arcId");
        block->GetPointData()->AddArray(arcIdArray);
        arcIdArray->Fill(arcCount);


        outputMB->SetBlock(arcCount, block);

        arcCount++;
      }

    }


std::cout<<"all done"<<std::endl;


  return 1;
}


int createOutputDiagram(vtkUnstructuredGrid *vtu,
                 const ttk::DiagramType &diagram,
                 vtkDataArray *const inputScalars,
                 float* normals,
                 const bool embedInDomain,
                 int dim,const ttk::Debug &dbg) {

  if(diagram.empty()) {
    dbg.printErr("Empty diagram");
    return -1;
  }
  dbg.printMsg("Create diagram vtu");



  const auto pd = vtu->GetPointData();
  const auto cd = vtu->GetCellData();

  if(pd == nullptr || cd == nullptr) {
    dbg.printErr("Grid has no point data or no cell data");
    return -2;
  }

  // point data arrays
  SimplexId nbPoints = diagram.size();

  vtkNew<ttkSimplexIdTypeArray> vertsId{};
  vertsId->SetName(ttk::VertexScalarFieldName);
  vertsId->SetNumberOfTuples(2 * diagram.size());
  pd->AddArray(vertsId);

  vtkNew<vtkIntArray> critType{};
  critType->SetName(ttk::PersistenceCriticalTypeName);
  critType->SetNumberOfTuples(2 * diagram.size());
  pd->AddArray(critType);

  vtkNew<vtkFloatArray> coordsScalars{};

  if(!embedInDomain) {
    coordsScalars->SetNumberOfComponents(3);
    coordsScalars->SetName(ttk::PersistenceCoordinatesName);
    coordsScalars->SetNumberOfTuples(2 * diagram.size());
    pd->AddArray(coordsScalars);
  }

  // cell data arrays

  vtkNew<ttkSimplexIdTypeArray> pairsId{};
  pairsId->SetName(ttk::PersistencePairIdentifierName);
  pairsId->SetNumberOfTuples(diagram.size());
  cd->AddArray(pairsId);

  vtkNew<ttkSimplexIdTypeArray> birthIdArray{};
  birthIdArray->SetName("BirthId");
  birthIdArray->SetNumberOfTuples(diagram.size());
  cd->AddArray(birthIdArray);

  vtkNew<vtkIntArray> pairsDim{};
  pairsDim->SetName(ttk::PersistencePairTypeName);
  pairsDim->SetNumberOfTuples(diagram.size());
  cd->AddArray(pairsDim);

  vtkSmartPointer<vtkDataArray> persistence{inputScalars->NewInstance()};
  persistence->SetName(ttk::PersistenceName);
  persistence->SetNumberOfTuples(diagram.size());
  cd->AddArray(persistence);

  vtkSmartPointer<vtkDataArray> birthScalars{inputScalars->NewInstance()};
  birthScalars->SetName(ttk::PersistenceBirthName);
  birthScalars->SetNumberOfTuples(diagram.size());
  cd->AddArray(birthScalars);

  vtkNew<vtkUnsignedCharArray> isFinite{};
  isFinite->SetName(ttk::PersistenceIsFinite);
  isFinite->SetNumberOfTuples(diagram.size());
  cd->AddArray(isFinite);

  vtkNew<vtkFloatArray> normalArray{};
  normalArray->SetNumberOfComponents(3);
  normalArray->SetNumberOfTuples( diagram.size());
  normalArray->SetName("Normals");
  cd->AddArray(normalArray);

  // grid

  vtkNew<vtkPoints> points{};
  points->SetNumberOfPoints(2 * diagram.size());
  vtkNew<vtkIdTypeArray> offsets{}, connectivity{};
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(diagram.size() + 1);
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(2 * diagram.size());

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(dbg.getThreadNumber())
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diagram.size(); ++i) {
    const auto &pair{diagram[i]};
    const auto i0{2 * i + 0}, i1{2 * i + 1};
    if(embedInDomain) {
      points->SetPoint(
        i0, pair.birth.coords[0], pair.birth.coords[1], pair.birth.coords[2]);
      points->SetPoint(
        i1, pair.death.coords[0], pair.death.coords[1], pair.death.coords[2]);
    } else {
      points->SetPoint(i0, pair.birth.sfValue, pair.birth.sfValue, 0);
      points->SetPoint(i1, pair.birth.sfValue, pair.death.sfValue, 0);
    }

    connectivity->SetTuple1(i0, i0);
    connectivity->SetTuple1(i1, i1);
    offsets->SetTuple1(i, 2 * i);

    // point data
    vertsId->SetTuple1(i0, pair.birth.id);
    vertsId->SetTuple1(i1, pair.death.id);
    critType->SetTuple1(i0, static_cast<ttk::SimplexId>(pair.birth.type));
    critType->SetTuple1(i1, static_cast<ttk::SimplexId>(pair.death.type));

    // float nx = normals[3*pair.birth.id+0];
    // float ny = normals[3*pair.birth.id+1];
    // float nz = normals[3*pair.birth.id+2];
    normalArray->SetTuple(i,normals+3*pair.birth.id);
    // normalArray->SetTuple(i1,normals+3*pair.birth.id);
    // normalArray->SetTuple(i1,normals+3*pair.birth.id);

    if(!embedInDomain) {
      coordsScalars->SetTuple3(
        i0, pair.birth.coords[0], pair.birth.coords[1], pair.birth.coords[2]);
      coordsScalars->SetTuple3(
        i1, pair.death.coords[0], pair.death.coords[1], pair.death.coords[2]);
    }

    // cell data
    pairsId->SetTuple1(i, i);
    birthIdArray->SetTuple1(i, pair.birth.id);
    persistence->SetTuple1(i, pair.persistence());
    birthScalars->SetTuple1(i, pair.birth.sfValue);
    isFinite->SetTuple1(i, pair.isFinite);
    pairsDim->SetTuple1(
      i, (pair.dim == 2 && pair.isFinite) ? dim - 1 : pair.dim);
  }
  offsets->SetTuple1(diagram.size(), connectivity->GetNumberOfTuples());

  vtkNew<vtkCellArray> cells{};
  cells->SetData(offsets, connectivity);
  vtu->SetPoints(points);
  vtu->SetCells(VTK_LINE, cells);

  if(!embedInDomain) {
    // highest birth (last pair)
    // const auto lastPair = std::max_element(diagram.begin(), diagram.end());
    // add diagonal (first point -> last birth/penultimate point)
    // std::array<vtkIdType, 2> diag{
    //   0, 2 * std::distance(diagram.begin(), lastPair)};
    // vtu->InsertNextCell(VTK_LINE, 2, diag.data());
    // pairsId->InsertTuple1(diagram.size(), -1);
    // pairsDim->InsertTuple1(diagram.size(), -1);
    // isFinite->InsertTuple1(diagram.size(), false);
    // // persistence of global min-max pair
    // const auto maxPersistence = diagram[0].persistence();
    // persistence->InsertTuple1(diagram.size(), 2 * maxPersistence);
    // // birth == death == 0
    // birthScalars->InsertTuple1(diagram.size(), 0);
    // birthIdArray->InsertTuple1(diagram.size(), pair.birth.id);
  }

  vtkNew<ttkSimplexIdTypeArray> parentComponent{};
  parentComponent->SetNumberOfTuples(diagram.size());
  parentComponent->SetNumberOfComponents(1);
  parentComponent->SetName("parentComponent");

  vtkNew<ttkSimplexIdTypeArray> componentSize{};
  componentSize->SetNumberOfTuples(diagram.size());
  componentSize->SetNumberOfComponents(1);
  componentSize->SetName("componentSize");

  vtkNew<ttkSimplexIdTypeArray> componentCurrentSize{};
  componentCurrentSize->SetNumberOfTuples(diagram.size());
  componentCurrentSize->SetNumberOfComponents(1);
  componentCurrentSize->SetName("componentSizeAtThreshold");

  vtkNew<ttkSimplexIdTypeArray> componentCurrentParent{};
  componentCurrentParent->SetNumberOfTuples(diagram.size());
  componentCurrentParent->SetNumberOfComponents(1);
  componentCurrentParent->SetName("componentParentAtThreshold");

  vtkNew<ttkSimplexIdTypeArray> currentParentSize{};
  currentParentSize->SetNumberOfTuples(diagram.size());
  currentParentSize->SetNumberOfComponents(1);
  currentParentSize->SetName("currentParentSize");

  std::vector<int> sizesAtThreshold(diagram.size(), 1);

  // compute size at current persistence threshold
  //

  // for(size_t i =diagram.size()-1; 0<i; i--){
  //  auto p = diagram[i].persistence();
  //  if(p<persThresh){
  //     sizesAtThreshold[i] = 0;
  //     sizesAtThreshold[diagram[i].cc_parent] += diagram[i].cc_size;
  //  }
  // }

  for(size_t i =0; i<diagram.size(); i++){
    auto p = diagram[i].persistence();
    componentCurrentParent->SetTuple1(i,i);
    currentParentSize->SetTuple1(i,sizesAtThreshold[componentCurrentParent->GetTuple1(i)]);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(dbg.getThreadNumber())
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diagram.size(); ++i) {
    parentComponent->SetTuple1(i, diagram[i].cc_parent);
    componentSize->SetTuple1(i, diagram[i].cc_size);
    componentCurrentSize->SetTuple1(i, sizesAtThreshold[i]);
  }


  vtu->GetCellData()->AddArray(componentSize);
  vtu->GetCellData()->AddArray(componentCurrentSize);
  vtu->GetCellData()->AddArray(componentCurrentParent);
  vtu->GetCellData()->AddArray(parentComponent);
  vtu->GetCellData()->AddArray(currentParentSize);
  dbg.printMsg("Done");

  return 0;
}

void fillDiagramValues(float* coordinates, DiagramType& diagram, const ttk::Debug &dbg){

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(dbg.getThreadNumber())
#endif // TTK_ENABLE_OPENMP
  for(int i=0; i<diagram.size(); i++){
    auto& p = diagram[i];
    auto birthId = p.birth.id;
    auto parentId = diagram[p.cc_parent].birth.id;
    p.birth.coords[0] = coordinates[3*birthId];
    p.birth.coords[1] = coordinates[3*birthId+1];
    p.birth.coords[2] = coordinates[3*birthId+2];

    p.death.coords[0] = coordinates[3*parentId];
    p.death.coords[1] = coordinates[3*parentId+1];
    p.death.coords[2] = coordinates[3*parentId+2];
  }
  // dbg.printMsg("Filled values");
}


void ttkRipsPersistence::computeCurrentSizes(vtkUnstructuredGrid* vtu, double persThresh){

  // printMsg("Computing current sizes");

  // auto &diagram = diagram_;

  NeedsUpdate = true;

}

