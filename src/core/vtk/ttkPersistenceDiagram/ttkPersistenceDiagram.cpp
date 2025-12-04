#include <cstddef>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkPersistenceDiagram.h>
#include <ttkPersistenceDiagramUtils.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistenceDiagram);

ttkPersistenceDiagram::ttkPersistenceDiagram() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

int ttkPersistenceDiagram::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkPersistenceDiagram::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <typename scalarType, typename triangulationType>
int ttkPersistenceDiagram::dispatch(
  vtkUnstructuredGrid *outputCTPersistenceDiagram,
  vtkDataArray *const inputScalarsArray,
  const scalarType *const inputScalars,
  scalarType *outputScalars,
  SimplexId *outputOffsets,
  int *outputMonotonyOffsets,
  const SimplexId *const inputOrder,
  const triangulationType *triangulation) {

  int status{};
  ttk::DiagramType CTDiagram{};

  if(BackEnd == BACKEND::APPROXIMATE_TOPOLOGY) {
    std::cout << "Chosen approx" << std::endl;
    double *range = inputScalarsArray->GetRange(0);
    this->setDeltaApproximate(range[1] - range[0]);
    this->setOutputScalars(outputScalars);
    this->setOutputOffsets(outputOffsets);
    this->setOutputMonotonyOffsets(outputMonotonyOffsets);
  }

  status = this->execute(CTDiagram, inputScalars, inputScalarsArray->GetMTime(),
                         inputOrder, triangulation);

  // something wrong in baseCode
  if(status != 0) {
    this->printErr("PersistenceDiagram::execute() error code: "
                   + std::to_string(status));
    return 0;
  }

  if(CTDiagram.empty()) {
    this->printErr("Empty diagram!");
    return 0;
  }

  vtkNew<vtkUnstructuredGrid> const vtu{};

  // convert CTDiagram to vtkUnstructuredGrid
  DiagramToVTU(vtu, CTDiagram, inputScalarsArray, *this,
               triangulation->getDimensionality(), this->ShowInsideDomain);

  // if(1){
  //   fillConnectedComponentsInfo(vtu, CTDiagram, getCCSizes(), getCCParents());
  // }

  outputCTPersistenceDiagram->ShallowCopy(vtu);

  if(this->ClearDGCache && this->BackEnd == BACKEND::DISCRETE_MORSE_SANDWICH) {
    this->printMsg("Clearing DiscreteGradient cache...");
    ttk::dcg::DiscreteGradient::clearCache(*triangulation);
  }

  return 1;
}

int ttkPersistenceDiagram::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkUnstructuredGrid *outputCTPersistenceDiagram
    = vtkUnstructuredGrid::GetData(outputVector, 0);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    this->printErr("Wrong triangulation");
    return 0;
  }
#endif

  this->preconditionTriangulation(triangulation);

  vtkDataArray *inputScalars = this->GetInputArrayToProcess(0, inputVector);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputScalars) {
    this->printErr("Wrong input scalars");
    return 0;
  }
#endif

  vtkDataArray *offsetField = this->GetOrderArray(
    input, 0, triangulation, false, 1, ForceInputOffsetScalarField);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!offsetField) {
    this->printErr("Wrong input offsets");
    return 0;
  }
  if(offsetField->GetDataType() != VTK_INT
     and offsetField->GetDataType() != VTK_ID_TYPE) {
    this->printErr("Input offset field type not supported");
    return 0;
  }
#endif

  vtkNew<ttkSimplexIdTypeArray> outputOffsets{};
  outputOffsets->SetNumberOfComponents(1);
  outputOffsets->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputOffsets->SetName("outputOffsets");

  vtkNew<vtkIntArray> outputMonotonyOffsets{};
  outputMonotonyOffsets->SetNumberOfComponents(1);
  outputMonotonyOffsets->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputMonotonyOffsets->SetName("outputMonotonyffsets");
  outputMonotonyOffsets->FillComponent(0, 0);

  vtkSmartPointer<vtkDataArray> const outputScalars
    = vtkSmartPointer<vtkDataArray>::Take(inputScalars->NewInstance());
  outputScalars->SetNumberOfComponents(1);
  outputScalars->SetNumberOfTuples(inputScalars->GetNumberOfTuples());
  outputScalars->DeepCopy(inputScalars);
  outputScalars->SetName("Cropped");

  int status{};
  ttkVtkTemplateMacro(
    inputScalars->GetDataType(), triangulation->getType(),
    status = this->dispatch(
      outputCTPersistenceDiagram, inputScalars,
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(inputScalars)),
      static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputScalars)),
      static_cast<SimplexId *>(ttkUtils::GetVoidPointer(outputOffsets)),
      static_cast<int *>(ttkUtils::GetVoidPointer(outputMonotonyOffsets)),
      static_cast<SimplexId *>(ttkUtils::GetVoidPointer(offsetField)),
      static_cast<TTK_TT *>(triangulation->getData())));

  // shallow copy input Field Data
  outputCTPersistenceDiagram->GetFieldData()->ShallowCopy(
    input->GetFieldData());

  return status;
}

void ttkPersistenceDiagram::fillConnectedComponentsInfo(vtkUnstructuredGrid* vtu,
                                                        const ttk::DiagramType &diagram,
                                                        const std::vector<SimplexId> &cc_size,
                                                        const std::vector<SimplexId> &merges){
  vtkNew<vtkIntArray> parentComponent{};
  parentComponent->SetNumberOfTuples(diagram.size()+1);
  parentComponent->SetNumberOfComponents(1);
  parentComponent->SetName("parentComponent");

  vtkNew<vtkIntArray> componentSize{};
  componentSize->SetNumberOfTuples(diagram.size()+1);
  componentSize->SetNumberOfComponents(1);
  componentSize->SetName("componentSize");

  vtkNew<vtkIntArray> componentCurrentSize{};
  componentCurrentSize->SetNumberOfTuples(diagram.size()+1);
  componentCurrentSize->SetNumberOfComponents(1);
  componentCurrentSize->SetName("componentSizeAtThreshold");

  vtkNew<vtkIntArray> componentCurrentParent{};
  componentCurrentParent->SetNumberOfTuples(diagram.size()+1);
  componentCurrentParent->SetNumberOfComponents(1);
  componentCurrentParent->SetName("componentParentAtThreshold");

  vtkNew<vtkIntArray> currentParentSize{};
  currentParentSize->SetNumberOfTuples(diagram.size()+1);
  currentParentSize->SetNumberOfComponents(1);
  currentParentSize->SetName("currentParentSize");

  std::vector<int> sizesAtThreshold(diagram.size(), 1);

  // compute size at current persistence threshold

  for(size_t i =diagram.size()-1; 0<i; i--){
   auto p = diagram[i].persistence();
   if(p<this->PersThresh){
      sizesAtThreshold[i] = 0;
      sizesAtThreshold[diagram[i].cc_parent] += diagram[i].cc_size;
   }
  }

  for(size_t i =0; i<diagram.size(); i++){
    auto p = diagram[i].persistence();
    if(p>=this->PersThresh){
      componentCurrentParent->SetTuple1(i,i);
    }else{
      componentCurrentParent->SetTuple1(i,componentCurrentParent->GetTuple1(diagram[i].cc_parent));
   }
    currentParentSize->SetTuple1(i,sizesAtThreshold[componentCurrentParent->GetTuple1(i)]);
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(size_t i = 0; i < diagram.size(); ++i) {
    parentComponent->SetTuple1(i, diagram[i].cc_parent);
    componentSize->SetTuple1(i, diagram[i].cc_size);
    componentCurrentSize->SetTuple1(i, sizesAtThreshold[i]);
  }

  //diagonal
  parentComponent->SetTuple1(diagram.size(), -1);
  componentSize->SetTuple1(diagram.size(), -1);
  componentCurrentSize->SetTuple1(diagram.size(), -1);

  vtu->GetCellData()->AddArray(componentSize);
  vtu->GetCellData()->AddArray(componentCurrentSize);
  vtu->GetCellData()->AddArray(componentCurrentParent);
  vtu->GetCellData()->AddArray(parentComponent);
  vtu->GetCellData()->AddArray(currentParentSize);

  onlyTreeUpdate_=false;
}


