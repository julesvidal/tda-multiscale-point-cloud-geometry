#include "DataTypes.h"
#include <sstream>
#include <string>
#include <ttkVisualizeAlgebraicSphere.h>

#include <vtkDoubleArray.h>
#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <ttkMacros.h>
#include <ttkUtils.h>


using namespace ttk;

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkVisualizeAlgebraicSphere);

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
ttkVisualizeAlgebraicSphere::ttkVisualizeAlgebraicSphere() {
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
int ttkVisualizeAlgebraicSphere::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
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
int ttkVisualizeAlgebraicSphere::FillOutputPortInformation(int port, vtkInformation *info) {
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
int ttkVisualizeAlgebraicSphere::RequestData(vtkInformation *ttkNotUsed(request),
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet){
    printErr("No data set");
    return 0;
  }
  
  // read_input

  auto pd = inputDataSet->GetPointData();
  SimplexId nbPoints = inputDataSet->GetNumberOfPoints();
  vtkMultiBlockDataSet *outputMB = vtkMultiBlockDataSet::GetData(outputVector,0);
  outputMB->SetNumberOfBlocks(0);

  if(nbPoints==0){
    printWrn("Selection contains no points");
    outputMB->SetNumberOfBlocks(0);
    return 0;
  }

  // auto idArray = pd->GetArray("PointIds");
  // auto birthIdArray = pd->GetArray("BirthId");
  vtkDataArray* APSSParamsArray = pd->GetArray("FitCoords");
  
  if(!APSSParamsArray){
    printErr("No sphere parameters array");
    return 0;
  }
  int fitDim = APSSParamsArray->GetNumberOfComponents();
  if(fitDim<5){
    printErr("Not enough components in sphere parameters array");
    return 0;
  }

  vtkDataArray* basisCenterArray = inputDataSet->GetFieldData()->GetArray("basisCenter");
  std::vector<double> basisCenter(3,0);
  if(basisCenterArray){
    basisCenterArray->GetTuple(0, basisCenter.data());
  }
  std::stringstream ss;
  ss << basisCenter[0] <<","<<basisCenter[1]<<","<<basisCenter[2];
  printMsg("Used basis center " + ss.str());

  // vtkDataArray* selectionArray = birthIdArray == nullptr ? idArray : birthIdArray;

  // if(selectionArray == nullptr){
  //   printErr("Selection must contain a valid id array -- skipping");
  //   return;
  // }
  
  printMsg(std::to_string(nbPoints) + " points");

  outputMB->SetNumberOfBlocks(nbPoints);

  // get bounds
  double bounds[6];
  inputDataSet->GetBounds(bounds);
  double diag_aabb = (bounds[1]-bounds[0])*(bounds[1]-bounds[0]) + (bounds[3]-bounds[2])*(bounds[3]-bounds[2]) + (bounds[5]-bounds[4])*(bounds[5]-bounds[4]);
  int dimX=Resolution, dimY=Resolution, dimZ=Resolution;
  double offsetBounds[6];
  double offsetSpacing=Spacing*(0.05*sqrt(diag_aabb)+1);
  offsetBounds[0] = bounds[0]- offsetSpacing;
  offsetBounds[1] = bounds[1]+ offsetSpacing;
  offsetBounds[2] = bounds[2]- offsetSpacing;
  offsetBounds[3] = bounds[3]+ offsetSpacing;
  offsetBounds[4] = bounds[4]- offsetSpacing;
  offsetBounds[5] = bounds[5]+ offsetSpacing;
  double ox=offsetBounds[0],oy=offsetBounds[2], oz=offsetBounds[4];
  double dx = (offsetBounds[1]-offsetBounds[0])/(dimX-1);
  double dy = (offsetBounds[3]-offsetBounds[2])/(dimY-1);
  double dz = dimZ > 1 ? (offsetBounds[5]-offsetBounds[4])/(dimZ-1) : 0;


  for(SimplexId i=0; i<nbPoints;i++){

    // std::cout<<" point "<<i<<std::endl;
    vtkNew<vtkImageData> outputVTI;

    outputMB->SetBlock(i, outputVTI);

    outputVTI->SetDimensions(dimX,dimY,dimZ);
    outputVTI->SetOrigin(offsetBounds[0], offsetBounds[2], offsetBounds[4]);
    outputVTI->SetSpacing(dx,dy,dz);

    SimplexId gridNbPoints = dimX*dimY*dimZ;

    vtkNew<vtkDoubleArray> fitScalarField{};
    fitScalarField->SetNumberOfComponents(1);
    fitScalarField->SetNumberOfTuples(gridNbPoints);
    fitScalarField->SetName("f");
    fitScalarField->Fill(-1);
    double* fitSFPtr = (double*)ttkUtils::GetVoidPointer(fitScalarField);

    outputVTI->GetPointData()->AddArray(fitScalarField);

    // selectionVec.push_back(selectionArray->GetTuple1(i));

    double params[5];
    APSSParamsArray->GetTuple(i, params);
    double uc = params[0];
    double ulx = params[1];
    double uly = params[2];
    double ulz = params[3];

    if(fitDim==5){

      double uq = params[4];

      for(SimplexId I=0; I<dimX; I++){
        for(SimplexId J=0; J<dimY; J++){
          for(SimplexId K=0; K<dimZ; K++){
            SimplexId ID = I+J*dimX+K*dimX*dimY;
            double x = ox + I*dx - basisCenter[0];
            double y = oy + J*dy - basisCenter[1];
            double z = oz + K*dz - basisCenter[2];
            double fx = uc + ulx*x+uly*y+ulz*z + uq*(x*x+y*y+z*z);
            fitSFPtr[ID]= fx;
          }
        }
      }
    }else{

      double uqxx = params[4];
      double uqxy = params[7];
      double uqyx = uqxy;
      double uqxz = params[8];
      double uqzx = uqxz;
      double uqyy = params[5];
      double uqyz = params[9];
      double uqzy = uqyz;
      double uqzz = params[6];

      for(SimplexId I=0; I<dimX; I++){
        for(SimplexId J=0; J<dimY; J++){
          for(SimplexId K=0; K<dimZ; K++){
            SimplexId ID = I+J*dimX+K*dimX*dimY;
            double x = ox + I*dx - basisCenter[0];
            double y = oy + J*dy - basisCenter[1];
            double z = oz + K*dz - basisCenter[2];

            double Uqx_x = uqxx*x + uqxy*y + uqxz*z;
            double Uqx_y = uqyx*x + uqyy*y + uqyz*z;
            double Uqx_z = uqzx*x + uqzy*y + uqzz*z;

            double xUqx = Uqx_x*x + Uqx_y*y + Uqx_z*z;
            double fx = uc + ulx*x+uly*y+ulz*z + xUqx;
            fitSFPtr[ID]= fx;
          }
        }
      }


    }


    // copy arrays
    // std::cout<<"copy array "<<std::endl;
    for(int iArray=0; iArray<pd->GetNumberOfArrays(); iArray++){
      vtkDataArray* ar = pd->GetArray(iArray);
      // std::cout<<" array "<<ar->GetName()<<std::endl;
      // std::cout<<"    add"<<std::endl;
      vtkDataArray* addedArray  = ar->NewInstance();
      outputVTI->GetPointData()->AddArray(addedArray);
      // std::cout<<"    set name "<<std::endl;
      addedArray->SetName(ar->GetName());
      int nComp = ar->GetNumberOfComponents();
      // std::cout<<"    nComp "<<nComp<<std::endl;
      // std::cout<<"    set ncomp "<<std::endl;
      addedArray->SetNumberOfComponents(nComp);
      // std::cout<<"    set ntuples "<<std::endl;
      addedArray->SetNumberOfTuples(gridNbPoints);
      

      // std::cout<<"    set tuple "<<std::endl;
      for(int iPoint = 0; iPoint<gridNbPoints; iPoint++){
        addedArray->SetTuple(iPoint, ar->GetTuple(i));
      }


      
    }


  }

  // 



  return 1;
}
