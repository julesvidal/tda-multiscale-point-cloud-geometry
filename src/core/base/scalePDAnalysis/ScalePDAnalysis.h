/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::ScalePDAnalysis
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %ScalePDAnalysis class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'ScalePDAnalysis'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <PersistenceDiagram.h>
#include <PersistenceDiagramClustering.h>

namespace ttk {

  /**
   * The ScalePDAnalysis class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class ScalePDAnalysis : virtual public Debug {

    public:
  enum class SCALE_METRIC {
    SUM = 0,
    MAX = 1,
  };

  protected:
    double TimeLimit{1.0};
    SCALE_METRIC Metric{SCALE_METRIC::SUM};

  public:
    ScalePDAnalysis();


void computeDistanceMatrix(std::vector<std::vector<DiagramType>>& diagrams,
                                            std::vector<std::vector<double>>& distMat);

    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    }

    /**
     * TODO 3: Implementation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computeAverages(dataType *outputData,
                        const dataType *inputData,
                        const triangulationType *triangulation) const {
    return 0;
    }

  double computeWassersteinDistanceAuction(const DiagramType& D0, const DiagramType& D1, int p);
double computeWassersteinDistance(const DiagramType& D0, const DiagramType& D1);

  }; // ScalePDAnalysis class

} // namespace ttk
