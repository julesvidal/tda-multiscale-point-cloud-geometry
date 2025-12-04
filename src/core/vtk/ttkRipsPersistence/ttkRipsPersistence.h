/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
/// \ingroup vtk
/// \class ttkRipsPersistence
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::RipsPersistence module.
///
/// This VTK filter uses the ttk::RipsPersistence module to compute an averaging of
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
///   - standalone/RipsPersistence/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::RipsPersistence
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkRipsPersistenceModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkSetGet.h>
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
 *   ttkRipsPersistence
 * DEPENDS
 *   ttkAlgorithm
 *   VTK::FiltersSources
 */

// TTK Base Includes
#include <RipsPersistence.h>

class TTKRIPSPERSISTENCE_EXPORT ttkRipsPersistence
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::RipsPersistence // and we inherit from the base class
{
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};

public:
  ttkRipsPersistence();
  ~ttkRipsPersistence() override = default;

  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  vtkSetMacro(PersistenceRequest, const std::string &);
  vtkGetMacro(PersistenceRequest, std::string);

  vtkGetMacro(MaxEdgeLength, double)
  vtkSetMacro(MaxEdgeLength, double);

  vtkGetMacro(NbPrimitives, int);
  vtkSetMacro(NbPrimitives, int);

  vtkGetMacro(PersThresh, double);
  inline void SetPersThresh(double data){
    PersThresh = data;
    Modified();
    NeedsUpdate=true;
    // if(!diagram_.empty()){
    //   NeedsUpdate = false;
    // }
  }

  // Hierarchical sampling
  vtkSetMacro(SamplingFactor, double);
  vtkGetMacro(SamplingFactor, double);

  void SetScaleSamplingMinMaxCount(double min, double max, double c){
   ScaleSamplingMin=min;
   ScaleSamplingMax=max;
   ScaleSamplingCount=c;
   Modified();
  }

  vtkSetMacro(UseNF, int)
  vtkGetMacro(UseNF, int)

  vtkSetMacro(HierarchicalSampling, int)
  vtkGetMacro(HierarchicalSampling, int)

  vtkSetMacro(OutputSampling, int)
  vtkGetMacro(OutputSampling, int)

  vtkSetMacro(ScaleSamplingBase, double)
  vtkGetMacro(ScaleSamplingBase, double)

  vtkSetMacro(ScaleSamplingMin, double)
  vtkGetMacro(ScaleSamplingMin, double)

  vtkSetMacro(ScaleSamplingMax, double)
  vtkGetMacro(ScaleSamplingMax, double)

  vtkSetMacro(ScaleSamplingCount, double)
  vtkGetMacro(ScaleSamplingCount, double)

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkRipsPersistence *New();
  vtkTypeMacro(ttkRipsPersistence, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */

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

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;



  void computeCurrentSizes(vtkUnstructuredGrid* vtu, double persThresh);

  private:
  double MaxEdgeLength{-1};
  int NbPrimitives{1};
  bool NeedsUpdate{true};
  double PersThresh{0.0};
  double Ratio{0.0};
  std::string PersistenceRequest{0};

  bool HierarchicalSampling{false};
  bool OutputSampling{false};
  bool UseNF{true};
  double ScaleSamplingMin{1};
  double ScaleSamplingMax{1};
  double ScaleSamplingBase{10};
  int ScaleSamplingCount{2};
  double SamplingFactor{0.1};
};


void fillDiagramValues(float* coordinates, DiagramType& diagram, const ttk::Debug &dbg);

int createOutputDiagram(vtkUnstructuredGrid *vtu,
                  const ttk::DiagramType &diagram,
                  vtkDataArray *const inputScalars,
                  float* normals,
                  const bool embedInDomain,int dim,const ttk::Debug &dbg);
