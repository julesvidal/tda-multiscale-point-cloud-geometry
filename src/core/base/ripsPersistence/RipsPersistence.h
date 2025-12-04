/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::RipsPersistence
/// \author Your Name Here <your.email@address.here>
/// \date The Date Here.
///
/// This module defines the %RipsPersistence class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'RipsPersistence'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include "DataTypes.h"
#include "PCAnalysisUtils.h"
#include <Debug.h>
#include <Triangulation.h>
#include <PersistenceDiagramUtils.h>
#include <Sampling.h>
#include <DimensionReduction.h>
#include <cinttypes>
#include <exception>
#include <limits>
#include <omp.h>
#include <queue>
#include <string>

namespace ttk {

  using edgeTuple = std::tuple<SimplexId, SimplexId, double>;
  /**
   * The RipsPersistence class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class RipsPersistence : virtual public Debug {

    public:
      RipsPersistence();



      ~RipsPersistence(){
      }

      int preconditionTriangulation(
          ttk::AbstractTriangulation *triangulation) const {
        return 0;
      }

int findScaleId(SimplexId ripsBirthId, const std::vector<SimplexId>& pointsAtScaleStartIdx) const ;

template<class PointND>
void compute3DEmbeddingMultiscale(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    std::vector<std::vector<double>>& outputCoordinates) const;
// void compute3DEmbeddingMultiscale(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
//                                                     std::vector<std::vector<SimplexId>>& filteredToLocalId,
//                                                     std::vector<EpsilonSampling<double, MyPoint5D>>& epsilonSampler,
//                                                     std::vector<std::vector<double>>& fitCoordinates,
//                                                     std::vector<std::vector<float>>& mls3DCoordinates,
//                                                     std::vector<float>& pointCloudBarycenter,
//                                                     std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
//                                                     int dim,
//                   std::vector<std::vector<double>>& outputCoordinates) const;

template<class PointND>
void computeDensityWithGraph(const SimplexId nbPoints,
                                              const std::vector<edgeTuple>& edges,
                                              std::vector<edgeTuple>& outputEdges, 
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                              std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              std::vector<std::vector<float>>& mls3DCoordinates,
                                              std::vector<float>& pointCloudBarycenter,
                                              std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                              int dim,
                                              double maxEdgeLength,
                                              bool constrained_density /* choose to compute the density over all scales (not-constrained) or only on contiguous scales (constrained)*/,
                                              double sigma,
                                              bool usescaledistance,
                                              bool useedgedensity,
                                              const std::vector<float>& scales,
                                              std::vector<double>& density);
template<class PointND>
void computeDensityWithTree(const SimplexId nbPoints,
                                              const std::vector<edgeTuple>& edges,
                                              std::vector<edgeTuple>& outputEdges,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                              std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              std::vector<std::vector<float>>& mls3DCoordinates,
                                              std::vector<float>& pointCloudBarycenter,
                                              std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                              int dim,
                                              double maxEdgeLength,
                                              bool constrained_density /* choose to compute the density over all scales (not-constrained) or only on contiguous scales (constrained)*/,
                                              double sigma,
                                              bool usescaledistance,
                                              bool useedgedensity,
                                              const std::vector<float>& scales,
                                              std::vector<double>& density);
template<class PointND>
void computeDensity(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    double maxEdgeLength,
                                                    bool constrained_density /* choose to compute the density over all scales (not-constrained) or only on contiguous scales (constrained)*/,
                                                    double sigma,
                                                    std::vector<double>& densityVec);
void sortEdgesAscending();

template<class PointND>
void computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    double maxEdgeLength);
  template<class PointND>
void computeApproximatedEdgeVectorMultiScaleNanoFlann(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    double maxEdgeLength);

      // void filterDistanceMatrixByDensity(std::vector<std::vector<double>>& distMat, float min_density, float* density, std::vector<SimplexId>& sampling_indexes);
// void computeEpsilonSamplingEdgeVectorMultiScaleND(std::vector<std::vector<SimplexId>>& sampling,
      //                                               std::vector<std::vector<double>>& fitCoordinates,
      //                                               std::vector<std::vector<float>>& mls3DCoordinates,
      //                                               std::vector<float>& pointCloudBarycenter,
      //                                               std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
      //                                               std::vector<EpsilonSampling<pointCloud::MyPoint5D::Scalar, pointCloud::MyPoint5D>>& epsilonSamplers,
      //                                               int dim,
      //                                               double maxEdgeLength);
// void computeApproximatedEdgeVectorMultiScaleND(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
      //                                           std::vector<std::vector<SimplexId>>& filteredToLocalId,
      //                                                std::vector<std::vector<double>>& fitCoordinates,
      //                                                std::vector<std::vector<float>>& mls3DCoordinates,
      //                                                std::vector<float>& pointCloudBarycenter,
      //                                               std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
      //                                               int dim,
      //                                               double maxEdgeLength);
      // void computeAcceleratedEdgeVectorMultiScaleND(const std::vector<std::vector<SimplexId>>& sampling,
      //     const std::vector<std::vector<double>>& fitCoordinates,
      //     const std::vector<std::vector<float>>& mls3DCoordinates,
      //     const std::vector<float>& pointCloudBarycenter,
      //     std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
      //     int dim,
      //     double maxEdgeLength=std::numeric_limits<double>::max());
      // void computeEdgeVectorMultiScaleND(const std::vector<std::vector<SimplexId>>& sampling,
      //     const std::vector<std::vector<double>>& fitCoordinates,
      //     const std::vector<std::vector<float>>& mls3DCoordinates,
      //     const std::vector<float>& pointCloudBarycenter,
      //     std::function<double(const double*, const double*, const float*, const float*, const float*)> distanceFunction,
      //     int dim);
      // void computeEdgeVectorND(SimplexId nbPoints, float* coordinates, int dim);
      // void computeEdgeVectorAtScaleDistanceMatrix(std::vector<std::vector<double>>& distMat, std::vector<SimplexId>& sampling_indexes, std::vector<SimplexId>& densitySampling, float* density);
      // void computeEdgeVector(SimplexId nbPoints, float* coordinates, float* normals, double ratio);
      void setInitialSizesAndRepresentants(std::vector<SimplexId>& sizes, std::vector<SimplexId>& reps);
      void setInitialSizes(std::vector<SimplexId>& sizes);
      void setInitialSizesCopy(const std::vector<SimplexId>& sizes);
      void setEdges(std::vector<edgeTuple>& edges);
      void computePersistenceDiagramWithBirth(double maxEdgeLength,
                                              SimplexId nbPoints,
                                              std::vector<edgeTuple>& mergingEdges,
                                              const std::vector<double>& point_births);
      void computePersistenceDiagram(double maxEdgeLength, SimplexId nbPoints, std::vector<edgeTuple>& mergingEdges);
      void computePersistenceDiagram(double maxEdgeLength, SimplexId nbPoints);

// template<class PointND>
// void addEpsilonPairs(std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSamplers,
//                                       std::vector<std::vector<SimplexId>>& filteredIdxToEpsSampling, 
//                                       std::vector<std::vector<SimplexId>>& filteredToLocalId,
//                                       std::vector<SimplexId>& pointsAtScaleStartIdx);
      void sortPersistenceDiagram(SimplexId nbPoints=-1);
      void replaceBirthIds(std::vector<SimplexId>& sampling_indexes);
      void computeEdgeVectorWithSamplingNF(ScaleSpaceSampling<pointCloud::MyPoint6D::Scalar, pointCloud::MyPoint6D>& sampler,
          SimplexId nbPoints, float* coordinates, float* normals, double ratio);
      void computeEdgeVectorWithSampling(ScaleSpaceSampling<pointCloud::MyPoint6D::Scalar, pointCloud::MyPoint6D>& sampler,
          SimplexId nbPoints, float* coordinates, float* normals, double ratio);
      template <typename Sampler>
        void computeEdgeVectorWithSamplingND(Sampler& sampler, SimplexId nbPoints, float* coordinates, int _dim);

      void computeEdgeVectorAtScaleND(float* coordinates, int dim, std::vector<SimplexId>& sampling_indexes, std::vector<SimplexId>& densitySampling);

      template <typename Sampler>
        void computeEdgeVectorAtScaleDistanceMatrix(Sampler& sampler, int scale_index, std::vector<std::vector<double>>& distMat);

      bool compareEdgeTuple(const edgeTuple& e0, const edgeTuple& e1);

      void compute3DEmbedding(SimplexId _nbPoints,
          int _dim,
          const std::vector<edgeTuple>& _edges,
          std::vector<std::vector<double>>& _outputCoordinates)const ;

      std::vector<edgeTuple>& getEdges(){
        return edges_;
      }
      DiagramType& getDiagram(){
        return diagram_;
      }

      inline void setUseProjectiveNorm(bool data){
        useProjectiveNorm_=data;
      }
      inline void clear(){
        edges_.clear();
        diagram_.clear();
        sizes_.clear();
        parent_.clear();
        reps_.clear();
      }

    private:
      std::vector<edgeTuple> edges_{};
    protected:
      DiagramType diagram_{};
      std::vector<SimplexId> sizes_{};
      std::vector<SimplexId> reps_{};
      std::vector<SimplexId> parent_{};
      bool useProjectiveNorm_{true};

  }; // RipsPersistence class

  // template <typename Sampler>
  //   void RipsPersistence::computeEdgeVectorWithSamplingND(Sampler& sampler,
  //       SimplexId nbPoints, float* coordinates, int dim){

  //     printMsg("Compute Rips edges with "+std::to_string(dim)+" sampling");
  //     Timer tm{};

  //     // sampler.printTreesInfosNF();
  //     // auto pt = sampler.getPTest();
  //     edges_.clear();
  //     int nScale = sampler.scale_count();
  //     for(int iScale=nScale-1; 0<iScale; iScale--){

  //       double distMax = iScale == nScale-1 ? std::numeric_limits<double>::infinity() : sampler.scale(iScale+1)*sampler.samplingFactor();
  //       auto& tree = sampler.kdTree(iScale);

  //       Timer tm_scale{};

  //       LongSimplexId nbPointsAtScale = tree.samples().size();

  //       for(LongSimplexId id=0; id<nbPointsAtScale; id++){

  //         LongSimplexId i = tree.samples()[id];

  //         double ni = 0;
  //         for(int k=0; k<dim; k++){
  //           ni += coordinates[dim*i+k]*coordinates[dim*i+k];
  //         }

  //         auto tree_query =tree.range_neighbors(i,distMax);

  //         for(LongSimplexId j : tree_query){

  //           double distCoords2 = 0;

  //           for(int idim =0; idim<dim; idim++){
  //             distCoords2 += (coordinates[dim*j+idim]-coordinates[dim*i+idim])*(coordinates[dim*j+idim]-coordinates[dim*i+idim]);
  //           }
  //           double dist = std::sqrt(distCoords2);

  //           // ============= projective case ==========
  //           double nj=0;
  //           for(int k=0; k<dim; k++){
  //             nj += coordinates[dim*j+k]*coordinates[dim*j+k];
  //           }

  //           double dist2=0;
  //           for(int k=0; k<dim; k++){
  //             dist2 += coordinates[dim*i+k]*coordinates[dim*j+k];
  //           }
  //           double dist3=dist2*dist2/(ni*nj);
  //           double dist4 = std::sqrt(dist3);
  //           if(dist4>1){
  //             dist4=1;
  //           }
  //           double dist5 = std::acos(dist4);
  //           // ========================================

  //           if(ni<=1e-5 or nj<=1e-5){
  //           }else{
  //             if(i<j){
  //               edges_.emplace_back(i, j, dist5);
  //             }
  //           }
  //         }

  //         int d_progress = (int)nbPointsAtScale*0.03+1;
  //         if(id%d_progress==0){ // write progress
  //           std::stringstream msg;
  //           msg<<"Scale "<<iScale<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
  //           printMsg(msg.str(), (double)id/nbPointsAtScale, tm_scale.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);

  //         }
  //       }
  //       std::stringstream msg;
  //       msg<<"Scale "<<iScale<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
  //       printMsg(msg.str(), 1, tm_scale.getElapsedTime(), threadNumber_);
  //     }

  //     if(1){

  //       double distMax =  nScale<1 ? std::numeric_limits<double>::infinity() : sampler.scale(1)*sampler.samplingFactor()*1.01;
  //       auto& tree = sampler.kdTree(0);

  //       Timer tm_scale{};

  //       LongSimplexId nbPointsAtScale = tree.samples().size();

  //       for(LongSimplexId id=0; id<nbPointsAtScale; id++){

  //         LongSimplexId i = tree.samples()[id];

  //         auto tree_query =tree.range_neighbors(i,distMax);

  //         double ni = 0;
  //         for(int k=0; k<dim; k++){
  //           ni += coordinates[dim*i+k]*coordinates[dim*i+k];
  //         }

  //         for(LongSimplexId j : tree_query){

  //           double distCoords2 = 0;

  //           for(int idim =0; idim<dim; idim++){
  //             distCoords2 += (coordinates[dim*j+idim]-coordinates[dim*i+idim])*(coordinates[dim*j+idim]-coordinates[dim*i+idim]);
  //           }
  //           double dist = std::sqrt(distCoords2);

  //           // ============= projective case ==========
  //           double nj=0;
  //           for(int k=0; k<dim; k++){
  //             nj += coordinates[dim*j+k]*coordinates[dim*j+k];
  //           }

  //           double dist2=0;
  //           for(int k=0; k<dim; k++){
  //             dist2 += coordinates[dim*i+k]*coordinates[dim*j+k];
  //           }
  //           double dist3=dist2*dist2/(ni*nj);
  //           double dist4 = std::sqrt(dist3);
  //           if(dist4>1){
  //             dist4=1;
  //           }
  //           double dist5 = std::acos(dist4);
  //           // ========================================

  //           if(ni<=1e-5 or nj<=1e-5){
  //           }else{
  //             if(i<j){
  //               if(std::isnan(dist5)){
  //                 std::cout<<"emplace "<<dist5<<std::endl;
  //               }
  //               if(i>nbPoints){
  //                 std::cout<<"TOO BIG I "<<i<<std::endl;
  //               }
  //               edges_.emplace_back(i, j, dist5);
  //             }
  //           }
  //         }

  //         int d_progress = (int)nbPointsAtScale*0.03+1;
  //         if(id%d_progress==0){ // write progress
  //           std::stringstream msg;
  //           msg<<"Scale "<<0<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
  //           printMsg(msg.str(), (double)id/nbPointsAtScale, tm_scale.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);

  //         }
  //       }
  //       std::stringstream msg;
  //       msg<<"Scale "<<0<<" ("<<nbPointsAtScale<<", "<<distMax<<"), size "<<tree.sample_count();
  //       printMsg(msg.str(), 1, tm_scale.getElapsedTime(), threadNumber_);
  //     }

  //     std::stringstream msg;
  //     msg<<"Hierarchical Method : "<<edges_.size()<<" edges";
  //     printMsg(msg.str(), 1, tm.getElapsedTime(), threadNumber_);

  //     const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
  //       return compareEdgeTuple(e0, e1);
  //     };

  //     TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
  //   }



  // template <typename Sampler>
  //   void RipsPersistence::computeEdgeVectorAtScaleDistanceMatrix(Sampler& sampler, int scale_index, std::vector<std::vector<double>>& distMat){

  //     Timer tm{};
  //     printMsg("Computing exact list of edges from DISTANCE MATRIX");

  //     auto& tree = sampler.kdTree(scale_index);
  //     SimplexId nbPoints = tree.sample_count();
  //     auto& sampling_indexes = tree.samples();
  //     std::cout<<"scale index "<<scale_index<<" vs sampler scale count "<<sampler.scale_count()<<std::endl;
  //     printMsg(std::to_string(nbPoints)+" points here");
  //     printMsg("size of sampling vec : "+std::to_string(sampling_indexes.size()));
  //     if(sampling_indexes.size()!=nbPoints){
  //       printErr("size != nbPoints /!\\");
  //     }

  //     LongSimplexId size = (LongSimplexId)nbPoints*((LongSimplexId)nbPoints-1)/2;
  //     edges_.clear();
  //     edges_.resize(size);


// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(1)
// #endif // TTK_ENABLE_OPENMP
  //     for(LongSimplexId I=0; I<(LongSimplexId)nbPoints-1; I++){

  //       // LongSimplexId i = sampling_indexes[I];


  //       for(LongSimplexId J=I+1; J<(LongSimplexId)nbPoints; J++){

  //         // LongSimplexId j = sampling_indexes[J];


  //         LongSimplexId id = J-1 + I*(2*nbPoints-I-3)/2;

  //         double dist = distMat[I][J];
  //         // std::cout<<" dist is  "<<dist<<std::endl;
  //         edges_[id] = std::make_tuple(I, J, dist);
  //       }
  //     }

  //     std::cout<<" Base Method: "<<edges_.size()<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;

  //     const auto cmp = [=](const edgeTuple &e0, const edgeTuple &e1){
  //       return compareEdgeTuple(e0, e1);
  //     };

  //     TTK_PSORT(this->threadNumber_, edges_.begin(), edges_.end(), cmp);
  //     printMsg("Done, found ");
  //   }

template<class PointND>
void RipsPersistence::compute3DEmbeddingMultiscale(
                                                    std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    std::vector<std::vector<double>>& outputCoordinates
                                                    ) const {

  printMsg("Computing 3d Embedding",0,0, threadNumber_);


  size_t nb = 0;

  SimplexId nbScales = epsilonIdToFilteredId.size();
  std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);

  for (int iScale=0; iScale<nbScales; iScale++){
    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();
    if(iScale+1<nbScales){
      pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;
    }
    nb+=nbPointsAtScale;
    std::cout<<" iscale "<<iScale<<" nb "<<nb<<std::endl;
    std::cout<<" filteredToLocal size "<<filteredToLocalId[iScale].size()<<std::endl;

  }

  size_t nbnb = nb*nb;
  std::cout<<"matrix size "<<nbnb<<std::endl;
  std::vector<double> distanceMatrix(nbnb,0);

  std::cout<<"compute distance matrix with "<<nbnb<<" entries"<<std::endl;

  for(SimplexId id_multiscale_i=0; id_multiscale_i<nb; id_multiscale_i++){
    int iScale = findScaleId(id_multiscale_i, pointsAtScaleStartIdx);

    SimplexId id_mono_eps_i = id_multiscale_i-pointsAtScaleStartIdx[iScale];

    const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id


    const double* coordinates_i = fitCoordinates[iScale].data();
    const double* coordsi = &coordinates_i[dim*id_mono_local_i];
    const float* mlsCoordinates_i = mls3DCoordinates[iScale].data();



    for(SimplexId id_multiscale_j=id_multiscale_i+1; id_multiscale_j<nb; id_multiscale_j++){

      int jScale = findScaleId(id_multiscale_j, pointsAtScaleStartIdx);

      SimplexId id_mono_eps_j = id_multiscale_j-pointsAtScaleStartIdx[jScale];
      const SimplexId id_mono_local_j = filteredToLocalId[jScale][epsilonIdToFilteredId[jScale][id_mono_eps_j]];


      const double* coordinates_j = fitCoordinates[jScale].data();
      const double* coordsj = &coordinates_j[dim*id_mono_local_j];
      const float* mlsCoordinates_j = mls3DCoordinates[jScale].data();

      double dist_comp = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*id_mono_local_i], &mlsCoordinates_j[3*id_mono_local_j], pointCloudBarycenter.data(), dim);

      distanceMatrix[(size_t)(id_multiscale_i + nb*id_multiscale_j)] = dist_comp;
      distanceMatrix[(size_t)(id_multiscale_j + nb*id_multiscale_i)] = dist_comp;
    }
  }
  std::cout<<"dist mat computed"<<std::endl;

  DimensionReduction dr{};
  auto method = DimensionReduction::METHOD::PCA;
  dr.setIsInputDistanceMatrix(true);
  dr.setInputMethod(method);
  dr.setDebugLevel(debugLevel_);
  if(method == DimensionReduction::METHOD::MDS){
    dr.setThreadNumber(1);
  }
  dr.setInputModulePath("/home/julius/ttk/core/base/dimensionReduction");
  dr.setInputNumberOfComponents(dim);
  std::vector<std::vector<double>> dr_output;
  dr.execute(dr_output, distanceMatrix, nb, nb);

  std::cout<<"dr output size "<<dr_output.size()<<std::endl;



  outputCoordinates.resize(nbScales);
  for (int iScale=0; iScale<nbScales; iScale++){
    outputCoordinates[iScale].resize(mls3DCoordinates[iScale].size(), -2);
  }

  for(SimplexId id_multiscale_i=0; id_multiscale_i<nb; id_multiscale_i++){
    int iScale = findScaleId(id_multiscale_i, pointsAtScaleStartIdx);

    SimplexId id_mono_eps_i = id_multiscale_i-pointsAtScaleStartIdx[iScale];

    const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id
    // std::cout<<" i_es_ms "<<id_multiscale_i<<" "<<iScale<<" :";
    for(int k=0; k<3;k++){
      outputCoordinates[iScale][3*id_mono_local_i+k] = dr_output[k][id_multiscale_i];
      // std::cout<<" "<<dr_output[k][id_multiscale_i];
    }
    // std::cout<<std::endl;

  }
}

template<class PointND>
void RipsPersistence::computeDensityWithGraph(const SimplexId nbPoints,
                                              const std::vector<edgeTuple>& edges,
                                              std::vector<edgeTuple>& outputEdges,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                              std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              std::vector<std::vector<float>>& mls3DCoordinates,
                                              std::vector<float>& pointCloudBarycenter,
                                              std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                              int dim,
                                              double maxEdgeLength,
                                              bool constrained_density /* choose to compute the density over all scales (not-constrained) or only on contiguous scales (constrained)*/,
                                              double sigma,
                                              bool useScaleDistance,
                                              bool useEdgeDensity,
                                              const std::vector<float>& scales,
                                              std::vector<double>& density){

  Timer tm{};

  // bool tomato = false;

  printMsg("Computing density from graph, sigma="+ std::to_string(sigma));

  // int K = 1000;
  SimplexId nbScales = epsilonIdToFilteredId.size();

  double sigma2 = sigma*sigma;


  // if(maxEdgeLength<0){
  //   maxEdgeLength = std::numeric_limits<double>::max();
  // }

  std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);

  SimplexId totalNbPoints=0;
  for (int iScale=0; iScale<nbScales; iScale++){
    totalNbPoints+=epsilonIdToFilteredId[iScale].size();
  }
  density.resize(nbPoints);

  // std::ofstream nFile, dFile;
  // nFile.open ("neighborGraph.txt");
  // dFile.open ("dFile.txt");
  //

  for (int iScale=0; iScale<nbScales-1; iScale++){
    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();
    pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;

  }


  // copy edges
  outputEdges = edges;

  // build neighbor list
  std::vector<std::vector<std::pair<SimplexId,double>>> directNeighbors(nbPoints);
  // for(const auto& e : edges){
  //   if(std::get<2>(e)<=sigma){
  //     SimplexId i0 = std::get<0>(e);
  //     SimplexId i1 = std::get<1>(e);
  //     double distance = std::get<2>(e);
  //     directNeighbors[i0].emplace_back(i1, distance);
  //     directNeighbors[i1].emplace_back(i0, distance);
  //   }
  // }
  for(const auto& e : edges){
    SimplexId i0 = std::get<0>(e);
    SimplexId i1 = std::get<1>(e);

    double distance = std::get<2>(e);
    if(useScaleDistance){
      double scale_i0 = (double)scales[findScaleId(i0, pointsAtScaleStartIdx)];
      double scale_i1 = (double)scales[findScaleId(i1, pointsAtScaleStartIdx)];
      distance*=(scale_i0+scale_i1)/2;
    }
    if(distance<=sigma){
      directNeighbors[i0].emplace_back(i1, distance);
      directNeighbors[i1].emplace_back(i0, distance);
    }
  }
  for(SimplexId i=0; i<nbPoints; i++){
    std::sort(directNeighbors[i].begin(), directNeighbors[i].end(),
        [](const std::pair<SimplexId,double>& a, const std::pair<SimplexId,double>&  b) {
        return a.second < b.second;
        });
  }




  std::cout<<"check time "<<tm.getElapsedTime()<<" s."<<std::endl;
  int prev_prog = -1;
  std::vector<std::vector<bool>> visited_threads(threadNumber_, std::vector<bool>(nbPoints));
  std::vector<std::vector<double>> distance_threads(threadNumber_, std::vector<double>(nbPoints));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId id_multiscale_i=0; id_multiscale_i<nbPoints; id_multiscale_i++){
    if(id_multiscale_i%1000==0){
      std::cout<<" "<<id_multiscale_i<<"/"<<nbPoints<<std::endl;
    }

    double dens=0;

    const int tid = omp_get_thread_num();

    std::vector<double>& distance = distance_threads[tid];
    std::fill(distance.begin(), distance.end(), -1);

    std::vector<bool> visited = visited_threads[tid];
    std::fill(visited.begin(), visited.end(), false);

    using tt = std::tuple<SimplexId, double>;
    auto cmp = [](const tt& a, const tt& b) {
      return std::get<1>(a) > std::get<1>(b);
    };
    std::priority_queue<tt, std::vector<tt>, decltype(cmp)> stack(cmp);
    std::vector<SimplexId> neighborhood{};
    stack.emplace(id_multiscale_i, 0);
    distance[id_multiscale_i]=0;

    while(!stack.empty()){

      const auto& top = stack.top();
      SimplexId current = std::get<0>(top);
      double dist = std::get<1>(top);

      // std::cout<<"pop "<<dist<<std::endl;

      stack.pop();

      if(!visited[current]){
        visited[current]=true;

        neighborhood.push_back(current);

        // push neighbors to stack
        int mark = 0;
        for(const auto& n : directNeighbors[current]){
          double d = n.second;
          if(dist+d>sigma){
            break;
          }
          mark++;
        } 

        for(int iNeigh = mark-1; iNeigh>=0; iNeigh--){
          const auto& n = directNeighbors[current][iNeigh];
          SimplexId neigh = n.first;
          double d = n.second;
          if(!visited[neigh]){
            // std::cout<<" neigh "<<neigh<<" dist+d "<<dist+d<<" distance "<<distance[neigh]<<std::endl;
            if((distance[neigh]<0 or distance[neigh]>dist+d)){
              stack.emplace(neigh, dist+d);
              // std::cout<<"emplace"<<std::endl;
              distance[neigh] = dist+d;
            }
          }
        }
        // std::cout<<"next"<<std::endl;
      }

    } // end while 
      
    TTK_PSORT(this->threadNumber_, neighborhood.begin(), neighborhood.end());
    auto last = std::unique(neighborhood.begin(), neighborhood.end());
    neighborhood.erase(last, neighborhood.end());


    for(auto n : neighborhood){
      int scale_n = findScaleId(n, pointsAtScaleStartIdx);
      SimplexId id_mono_eps_n = n-pointsAtScaleStartIdx[scale_n];
      SimplexId size = epsilonSampler[scale_n].size(id_mono_eps_n);
      double dist = distance[n];
      if(dist<0 or dist > sigma){
        std::cout<<"ERROR"<<std::endl;
      }
      double dist2 = dist*dist;
      dens += (dist2/sigma2 - 1)*(dist2/sigma2 - 1)*size;
    }

    density[id_multiscale_i] = std::log(dens);
    // if(id_multiscale_i == 2193 or id_multiscale_i == 298 ) {
    // std::cout<<" i "<<id_multiscale_i<<" : density "<<density[id_multiscale_i]<<" and neighbordhood of size "<<neighborhood_size<<std::endl;
    //   }
    


  } // end for id_multiscale_i

  double tm_vertices = tm.getElapsedTime();
  std::cout<<" computed density for "<<nbPoints<<" vertices in "<<tm.getElapsedTime()<<" s."<<std::endl;
  if(useEdgeDensity){

    std::cout<<" now computing density for "<<outputEdges.size()<<" edges: "<<std::endl;
    // compute density value for edges:
    for(auto& e : outputEdges){
      SimplexId i0 = std::get<0>(e);
      SimplexId i1 = std::get<1>(e);

      double dist = std::get<2>(e);
      if(useScaleDistance){
        double scale_i0 = (double)scales[findScaleId(i0, pointsAtScaleStartIdx)];
        double scale_i1 = (double)scales[findScaleId(i1, pointsAtScaleStartIdx)];
        dist*=(scale_i0+scale_i1)/2;
      }

      if(dist >= 2*sigma){
        std::get<2>(e) = -0.0;
        continue;
      }
      if(dist < 0.1*sigma){
        std::get<2>(e) = - std::min(density[i0], density[i1]);
        continue;
      }

      // compute density of the middle point
      double dens=0;
    std::vector<double> distance(nbPoints, -1);
    std::vector<bool> visited(nbPoints,false);

    using tt = std::tuple<SimplexId, double>;
    auto cmp = [](const tt& a, const tt& b) {
      return std::get<1>(a) > std::get<1>(b);
    };
    std::priority_queue<tt, std::vector<tt>, decltype(cmp)> stack(cmp);
    std::vector<SimplexId> neighborhood{};
    stack.emplace(i0, dist/2);

    distance[i0]=dist/2;
    stack.emplace(i1, dist/2);
    distance[i1]=dist/2;

    while(!stack.empty()){

      const auto& top = stack.top();
      SimplexId current = std::get<0>(top);
      double dist = std::get<1>(top);

      // std::cout<<"pop "<<dist<<std::endl;

      stack.pop();

      if(!visited[current]){
        visited[current]=true;

        neighborhood.push_back(current);

        // push neighbors to stack
        int mark = 0;
        for(const auto& n : directNeighbors[current]){
          double d = n.second;
          if(dist+d>sigma){
            break;
          }
          mark++;
        } 

        for(int iNeigh = mark-1; iNeigh>=0; iNeigh--){
          const auto& n = directNeighbors[current][iNeigh];
          SimplexId neigh = n.first;
          double d = n.second;
          if(!visited[neigh]){
            // std::cout<<" neigh "<<neigh<<" dist+d "<<dist+d<<" distance "<<distance[neigh]<<std::endl;
            if((distance[neigh]<0 or distance[neigh]>dist+d)){
              stack.emplace(neigh, dist+d);
              // std::cout<<"emplace"<<std::endl;
              distance[neigh] = dist+d;
            }
          }
        }
        // std::cout<<"next"<<std::endl;
      }

    } // end while 
      
    TTK_PSORT(this->threadNumber_, neighborhood.begin(), neighborhood.end());
    auto last = std::unique(neighborhood.begin(), neighborhood.end());
    neighborhood.erase(last, neighborhood.end());


    for(auto n : neighborhood){
      int scale_n = findScaleId(n, pointsAtScaleStartIdx);
      SimplexId id_mono_eps_n = n-pointsAtScaleStartIdx[scale_n];
      SimplexId size = epsilonSampler[scale_n].size(id_mono_eps_n);
      double dist = distance[n];
      if(dist<0 or dist > sigma){
        std::cout<<"ERROR"<<std::endl;
      }
      double dist2 = dist*dist;
      dens += (dist2/sigma2 - 1)*(dist2/sigma2 - 1)*size;
    }

      dens = std::log(dens);
      std::get<2>(e) = -std::min({dens, density[i0], density[i1]});
    }
  std::cout<<" computed density for edges in "<<tm.getElapsedTime()-tm_vertices<<" s."<<std::endl;
  }

  printMsg("density done", 1, tm.getElapsedTime(), threadNumber_);

}
template<class PointND>
void RipsPersistence::computeDensityWithTree(const SimplexId nbPoints,
                                              const std::vector<edgeTuple>& edges,
                                              std::vector<edgeTuple>& outputEdges,
                                              std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                              std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                              std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                              std::vector<std::vector<double>>& fitCoordinates,
                                              std::vector<std::vector<float>>& mls3DCoordinates,
                                              std::vector<float>& pointCloudBarycenter,
                                              std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                              int dim,
                                              double maxEdgeLength,
                                              bool constrained_density /* choose to compute the density over all scales (not-constrained) or only on contiguous scales (constrained)*/,
                                              double sigma,
                                              bool useScaleDistance,
                                              bool useEdgeDensity,
                                              const std::vector<float>& scales,
                                              std::vector<double>& density){

  Timer tm{};

  // bool tomato = false;

  printMsg("Computing density using tree version, sigma="+ std::to_string(sigma));

  // int K = 1000;
  SimplexId nbScales = epsilonIdToFilteredId.size();

  double sigma2 = sigma*sigma;


  // if(maxEdgeLength<0){
  //   maxEdgeLength = std::numeric_limits<double>::max();
  // }

  std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);

  SimplexId totalNbPoints=0;
  for (int iScale=0; iScale<nbScales; iScale++){
    totalNbPoints+=epsilonIdToFilteredId[iScale].size();
  }
  density.resize(nbPoints);

  // std::ofstream nFile, dFile;
  // nFile.open ("neighborGraph.txt");
  // dFile.open ("dFile.txt");
  //

  for (int iScale=0; iScale<nbScales-1; iScale++){
    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();
    pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;

  }

  // copy edges
  outputEdges = edges;

  // build neighbor list
  std::vector<std::vector<std::pair<SimplexId,double>>> directNeighbors(nbPoints);
  for(const auto& e : edges){
    SimplexId i0 = std::get<0>(e);
    SimplexId i1 = std::get<1>(e);

    double distance = std::get<2>(e);
    if(useScaleDistance){
      double scale_i0 = (double)scales[findScaleId(i0, pointsAtScaleStartIdx)];
      double scale_i1 = (double)scales[findScaleId(i1, pointsAtScaleStartIdx)];
      distance*=(scale_i0+scale_i1)/2;
    }
    if(distance<=sigma){
      directNeighbors[i0].emplace_back(i1, distance);
      directNeighbors[i1].emplace_back(i0, distance);
    }
  }




  int prev_prog = -1;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId id_multiscale_i=0; id_multiscale_i<nbPoints; id_multiscale_i++){


    // int progress = (double)id_multiscale_i/nbPoints*100;

    // if(progress%5==0){
    //   if(progress!=prev_prog){
    //     std::cout<<" progress "<<progress<<"% after "<<tm.getElapsedTime()<<"s."<<std::endl;
    //     prev_prog=progress;
    //   }
    // }
    //build neighborhood list
    SimplexId neighborhood_size=0;

    // int iScale = findScaleId(id_multiscale_i, pointsAtScaleStartIdx);
    // SimplexId id_mono_eps_i = id_multiscale_i-pointsAtScaleStartIdx[iScale];
    // const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id
    // const double* coordinates_i = fitCoordinates[iScale].data();
    // const double* coordsi = &coordinates_i[dim*id_mono_local_i];
    // const float* mlsCoordinates_i = mls3DCoordinates[iScale].data();

    double dens=0;

    std::stack<std::tuple<SimplexId,SimplexId, double>> stack;

    stack.emplace(id_multiscale_i, -1, 0);

    while(!stack.empty()){

      const auto& top = stack.top();
      SimplexId current = std::get<0>(top);
      SimplexId prev = std::get<1>(top);
      double dist = std::get<2>(top);

    // if(id_multiscale_i == 2193 or id_multiscale_i == 298 ) {
    //   std::cout<<"     pop "<<current<<" "<<prev<<std::endl;
    // }
      stack.pop();

    // if(id_multiscale_i == 2193 or id_multiscale_i == 298 ) {
    //   std::cout<<"        cur "<<current<<"/"<<nbPoints<<std::endl;
    // }

      int scale_current = findScaleId(current, pointsAtScaleStartIdx);
      // std::cout<<"          scale "<<scale_current<<"/"<<nbScales-1<<std::endl;
      SimplexId id_mono_eps_current = current-pointsAtScaleStartIdx[scale_current];
      // std::cout<<"          id_eps "<<id_mono_eps_current<<"/"<<epsilonIdToFilteredId[scale_current].size()<<std::endl;
      // const SimplexId id_mono_local_current = filteredToLocalId[scale_current][epsilonIdToFilteredId[scale_current][id_mono_eps_current]]; //local id
      // std::cout<<"          id_local "<<id_mono_local_current<<std::endl;
                                                                                                                                           //
      // const double* coordinates_current = fitCoordinates[scale_current].data();
      // const double* coords_current = &coordinates_current[dim*id_mono_local_current];
      // const float* mlscoordinates_current = mls3DCoordinates[scale_current].data();

      // dist = distanceFunction(coordsi, coords_current, &mlsCoordinates_i[3*id_mono_local_i], &mlscoordinates_current[3*id_mono_local_current], pointCloudBarycenter.data(), dim);

    // if(id_multiscale_i == 2193 or id_multiscale_i == 298 ) {
    //   std::cout<<"          dist  "<<dist<<std::endl;
    // }


      neighborhood_size++;

      // find size
      SimplexId size = epsilonSampler[scale_current].size(id_mono_eps_current);

      // add density
      double dist2 = dist*dist;
      dens += (dist2/sigma2 - 1)*(dist2/sigma2 - 1)*size;

      // push neighbors to stack
      for(const auto& n : directNeighbors[current]){
        SimplexId neigh = n.first;
        double d = n.second;
        if(neigh != prev and (dist+d <= sigma)){
          stack.emplace(neigh, current, dist+d);
          // if(id_multiscale_i == 2193 or id_multiscale_i == 298 ) {
          //   std::cout<<"          push "<<neigh<<", "<<current<<std::endl;

          // }
        }
      } 
      // std::cout<<"          - done"<<std::endl;
    } // end while 

    density[id_multiscale_i] = std::log(dens);
    // if(id_multiscale_i == 2193 or id_multiscale_i == 298 ) {
    // std::cout<<" i "<<id_multiscale_i<<" : density "<<density[id_multiscale_i]<<" and neighbordhood of size "<<neighborhood_size<<std::endl;
    //   }



  } // end for id_multiscale_i
    
  double tm_vertices = tm.getElapsedTime();
  std::cout<<" computed density for "<<nbPoints<<" vertices in "<<tm.getElapsedTime()<<" s."<<std::endl;
  if(useEdgeDensity){

    std::cout<<" now computing density for "<<outputEdges.size()<<" edges: "<<std::endl;
    // compute density value for edges:
    for(auto& e : outputEdges){
      SimplexId i0 = std::get<0>(e);
      SimplexId i1 = std::get<1>(e);

      double distance = std::get<2>(e);
      if(useScaleDistance){
        double scale_i0 = (double)scales[findScaleId(i0, pointsAtScaleStartIdx)];
        double scale_i1 = (double)scales[findScaleId(i1, pointsAtScaleStartIdx)];
        distance*=(scale_i0+scale_i1)/2;
      }

      if(distance >= 2*sigma){
        std::get<2>(e) = -0.0;
        continue;
      }
      if(distance < 0.1*sigma){
        std::get<2>(e) = - std::min(density[i0], density[i1]);
        continue;
      }

      // compute density of the middle point
      double dens=0;
      std::stack<std::tuple<SimplexId,SimplexId, double>> stack;
      stack.emplace(i0, i1, distance/2);
      stack.emplace(i1, i0, distance/2);

      while(!stack.empty()){

        const auto& top = stack.top();
        SimplexId current = std::get<0>(top);
        SimplexId prev = std::get<1>(top);
        double dist = std::get<2>(top);

        stack.pop();


        int scale_current = findScaleId(current, pointsAtScaleStartIdx);
        SimplexId id_mono_eps_current = current-pointsAtScaleStartIdx[scale_current];


        // find size
        SimplexId size = epsilonSampler[scale_current].size(id_mono_eps_current);

        // add density
        double dist2 = dist*dist;
        dens += (dist2/sigma2 - 1)*(dist2/sigma2 - 1)*size;

        // push neighbors to stack
        for(const auto& n : directNeighbors[current]){
          SimplexId neigh = n.first;
          double d = n.second;
          if(neigh != prev and (dist+d <= sigma)){
            stack.emplace(neigh, current, dist+d);
          }
        } 
      } // end while 

      dens = std::log(dens);
      std::get<2>(e) = -std::min({dens, density[i0], density[i1]});
    }
  std::cout<<" computed density for edges in "<<tm.getElapsedTime()-tm_vertices<<" s."<<std::endl;
  }

  printMsg("density done", 1, tm.getElapsedTime(), threadNumber_);

}

template<class PointND>
void RipsPersistence::computeDensity(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    double maxEdgeLength,
                                                    bool constrained_density /* choose to compute the density over all scales (not-constrained) or only on contiguous scales (constrained)*/,
                                                    double sigma,
                                                    std::vector<double>& density){

  Timer tm{};

  bool tomato = false;

  if(constrained_density){
    printMsg("Computing euclidean density (constrained to successive scales)");
  }else{
    printMsg("Computing euclidean density");
  }

  int K = 1000;
  SimplexId nbScales = epsilonIdToFilteredId.size();


  if(maxEdgeLength<0){
    maxEdgeLength = std::numeric_limits<double>::max();
  }

  std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);

  SimplexId totalNbPoints=0;
  for (int iScale=0; iScale<nbScales; iScale++){
    totalNbPoints+=epsilonIdToFilteredId[iScale].size();
  }

    std::ofstream nFile, dFile;
    nFile.open ("neighborGraph.txt");
    dFile.open ("dFile.txt");

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif // TTK_ENABLE_OPENMP
  for (int iScale=0; iScale<nbScales-1; iScale++){
    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();
    pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;

  }

  density.resize(totalNbPoints);


  double sigma2 = sigma*sigma;
  double radius = sigma;

  double lmin = log(0.03)/log(10);
  double lmax = log(0.2)/log(10);
  double delta = lmax-lmin;

  for (int iScale=0; iScale<nbScales; iScale++){

    std::cout<<"\nScale "<<iScale<<std::endl;
    // sigma = std::pow(10, lmin+iScale*delta/nbScales);
    std::cout<<"sigma "<<sigma<<std::endl;


    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();




    for(LongSimplexId id_mono_eps_i=0; id_mono_eps_i<nbPointsAtScale; id_mono_eps_i++){

      const SimplexId id_ms_i = pointsAtScaleStartIdx[iScale] + id_mono_eps_i;

      const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id
                                                                                                                 //

      std::vector<SimplexId> globalMatchesForTomato{};
      std::vector<SimplexId> globalDistForTomato{};

      // std::cout<<" Density for point "<<id_mono_eps_i<<" "<<id_ms_i<<std::endl;

      const double* coordinates_i = fitCoordinates[iScale].data();

        std::vector<double> queryPoint(dim);
        for(int k=0; k<dim; k++){
          queryPoint[k] = coordinates_i[dim*id_mono_local_i+k];
        }



      std::vector<unsigned int> ret_matches(nbScales*K);
      std::vector<double> ret_distances(nbScales*K);
      std::vector<SimplexId> ret_sizes(nbScales*K);
      // SimplexId n_matches_next_scale=0;
      // SimplexId n_matches_prev_scale=0;
      // SimplexId n_matches_same_scale=0;

      SimplexId n_matches=0;

      int start_scale = 0;
      int stop_scale = nbScales-1;

      if(constrained_density){
        start_scale = std::max(0, iScale-1);
        stop_scale = std::min(nbScales-1, iScale+1);
      }
      // for(int jScale = 0; jScale<nbScales; jScale++){

      double dens = 0;
      for(int jScale = start_scale; jScale<=stop_scale; jScale++){


        // sigma = std::pow(10, lmin+jScale*delta/nbScales);
        sigma2 = sigma*sigma;
        radius = sigma;

        auto tree = epsilonSampler[jScale].treeNF();
        std::vector<nanoflann::ResultItem<unsigned int, double>> tree_query_at_scale;
        // SimplexId n_matches_at_scale = tree->knnSearch(queryPoint.data(), K, ret_matches.data()+n_matches, ret_distances.data()+n_matches);
        SimplexId n_matches_at_scale = tree->radiusSearch(queryPoint.data(), radius*radius, tree_query_at_scale);

        auto& epsSampler = epsilonSampler[jScale];
        // for(int i_match = 0; i_match<n_matches_at_scale; i_match++){
        //   ret_sizes[i_match+n_matches] = epsSampler.size(ret_matches[i_match+n_matches]);
        // }

        // if(tomato){

        //   for(SimplexId i_match=n_matches; i_match<n_matches+n_matches_at_scale; i_match++){
        //     SimplexId match = ret_matches[i_match];
        //     SimplexId globalId = match + pointsAtScaleStartIdx[jScale];

        //     globalMatchesForTomato.push_back(globalId);
        //     globalDistForTomato.push_back(ret_distances[i_match]);
        //   }
        // }
        
        for(auto& m : tree_query_at_scale){

          SimplexId size = epsSampler.size(m.first);
          double dist2 = m.second;
          dens += (dist2/sigma2 - 1)*(dist2/sigma2 - 1)*size;
        }

        n_matches += n_matches_at_scale;
      }
      density[id_ms_i] = std::log(std::sqrt(dens));
      // density[id_ms_i] = dens;



      // ret_matches.resize(n_matches);
      // ret_distances.resize(n_matches);
      // ret_sizes.resize(n_matches);

      // // sort by ascending order of distances
      // std::vector<SimplexId> ascendingOrder(n_matches);
      // std::iota(ascendingOrder.begin(), ascendingOrder.end(), 0);
      // // std::cout<<"sort"<<std::endl;
      // std::sort(ascendingOrder.begin(), ascendingOrder.end(), [=](const SimplexId &a, const SimplexId &b){ return ret_distances[a]<ret_distances[b];});


      // double sum_dist = 0;
      // // double sum_size = epsilonSampler[iScale].size(id_mono_eps_i);
      // double sum_size = 0;
      // double dens=0;
      // int n_neighbors = std::min(K, n_matches);

      // for(int i=0; i<n_neighbors; i++){
      //   if(id_ms_i==59991){
      //     std::cout<<" neighbor "<<i<<", dist "<<ret_distances[ascendingOrder[i]]<<" and size "<<ret_sizes[ascendingOrder[i]]<<std::endl;
      //   }
      //   double dist2 = ret_distances[ascendingOrder[i]];

      //   if(dist2 < maxEdgeLength){
      //     dens += std::exp(-dist2/(2*sigma2))*ret_sizes[ascendingOrder[i]];
      //     sum_dist += dist2*ret_sizes[ascendingOrder[i]];
      //     sum_size+=ret_sizes[ascendingOrder[i]];
      //   }
      // }

      // if(sum_dist<1e-9){
      //   // density[id_ms_i]  = std::numeric_limits<double>::max();
      //   sum_dist = epsilonSampler[iScale].epsilon();
      //   std::cout<<"reached peak density"<<std::endl;
      // }

      // double d = std::log(dens);
      // double d = sum_size/sum_dist;
      // double d = sum_dist/sum_size;
      // density[id_ms_i] = d;
      // density[id_ms_i] = ld;

      // if(tomato){

      //   for(int i_match=0; i_match<globalMatchesForTomato.size(); i_match++){
      //     nFile<<globalMatchesForTomato[i_match];
      //     if(i_match<globalMatchesForTomato.size()-1){
      //       nFile<<", ";
      //     }else{
      //       nFile<<"\n";
      //     }
      //   }
      //   dFile<<density[id_ms_i]<<"\n";
      // }

    } // end for id_mono_eps_i

  } // end for iScale




    dFile.close();
    nFile.close();
    std::cout<<" done density in "<<tm.getElapsedTime()<<" s."<<std::endl;
}

template<class PointND>
void RipsPersistence::computeApproximatedEdgeVectorMultiScaleNanoFlannAllScalesMixed(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    double maxEdgeLength){

  Timer tm{};
  printMsg("Computing exact list of edges in "+std::to_string(dim)+"D, WITH SUBSAMPLING, all scales mixed, below "+std::to_string(maxEdgeLength));

  SimplexId nbScales = epsilonIdToFilteredId.size();

  edges_.clear();
  edges_.reserve(10000000);

  if(maxEdgeLength<0){
    maxEdgeLength = std::numeric_limits<double>::max();
  }

  long long int numberOfEdges = 0;


  std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);


  SimplexId totalNbPoints=0;
  for (int iScale=0; iScale<nbScales; iScale++){
    totalNbPoints+=epsilonIdToFilteredId[iScale].size();
  }
  std::cout<<" Total nb of points : "<<totalNbPoints<<std::endl;

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif // TTK_ENABLE_OPENMP
  for (int iScale=0; iScale<nbScales-1; iScale++){
    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();
    pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;

  }


  for (int iScale=0; iScale<nbScales; iScale++){
    std::cout<<"\nScale "<<iScale<<std::endl;


      SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();


      double distMax = std::numeric_limits<double>::max();
      distMax = std::min(distMax, (double)maxEdgeLength);

      long long unsigned int addedEdges=0;
      SimplexId addedEdges_monoscale = 0;
      SimplexId addedEdges_multiscale = 0;

      double max_added_dist=0;

      std::cout<<" nb samples: "<<nbPointsAtScale<< " (potentially "<<nbPointsAtScale*(nbPointsAtScale-1)<<" edges )"<<std::endl;


      for(LongSimplexId id_mono_eps_i=0; id_mono_eps_i<nbPointsAtScale; id_mono_eps_i++){

        const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id

        const double* coordinates_i = fitCoordinates[iScale].data();
        const float* mlsCoordinates_i = mls3DCoordinates[iScale].data();

        const double* coordsi = &coordinates_i[dim*id_mono_local_i];
        const SimplexId id_multiscale_i = id_mono_eps_i+pointsAtScaleStartIdx[iScale];

        std::vector<double> queryPoint(dim);
        for(int k=0; k<dim; k++){
          queryPoint[k] = coordinates_i[dim*id_mono_local_i+k];
        }
        // const double queryPoint[6] = {coordinates_i[5*id_mono_local_i+0],
        //                               coordinates_i[5*id_mono_local_i+1],
        //                               coordinates_i[5*id_mono_local_i+2],
        //                               coordinates_i[5*id_mono_local_i+3],
        //                               coordinates_i[5*id_mono_local_i+4]};

        for (int jScale=iScale; jScale<nbScales; jScale++){

        // Query in scale jScale >= iScale
        auto tree = epsilonSampler[jScale].treeNF();
        std::vector<nanoflann::ResultItem<unsigned int, double>> tree_query_same_scale;
        SimplexId n_matches_same_scale = tree->radiusSearch(queryPoint.data(), distMax*distMax, tree_query_same_scale);
        for(auto& match : tree_query_same_scale){

          const double dist_match = match.second;

          const SimplexId id_mono_eps_j = match.first;
          const SimplexId id_mono_local_j = filteredToLocalId[jScale][epsilonIdToFilteredId[jScale][id_mono_eps_j]];

          const double* coordinates_j = fitCoordinates[jScale].data();
          const double* coordsj = &coordinates_j[dim*id_mono_local_j];
          const float* mlsCoordinates_j = mls3DCoordinates[jScale].data();

          const SimplexId id_multiscale_j = id_mono_eps_j+pointsAtScaleStartIdx[jScale];

          if(id_multiscale_i<id_multiscale_j){ // to avoid duplicates

            double dist_comp = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*id_mono_local_i], &mlsCoordinates_j[3*id_mono_local_j], pointCloudBarycenter.data(), dim);
            // std::cout<<" dist match "<<std::sqrt(dist_match)<<"  vs  dist comp "<<dist_comp<<std::endl;

            if(dist_comp <= maxEdgeLength){
              edges_.emplace_back(id_multiscale_i, id_multiscale_j, dist_comp);

              if(dist_comp>max_added_dist){
                max_added_dist=dist_comp;
              }
              numberOfEdges++;
              addedEdges++;
              addedEdges_monoscale++;
            }
          }
        } // end for tree query same scal

        } //end for jscale
      } // end for id_mono_eps_i
      // std::cout<<"      added "<<addedEdges_monoscale<<" edges between "<<iScale<<" and "<<iScale<<std::endl;
      // std::cout<<"      added "<<addedEdges_multiscale<<" edges between "<<iScale<<" and "<<iScale+1<<std::endl;
      // std::cout<<"added "<<addedEdges<<" edges with max dist "<<max_added_dist<<std::endl;



  } // end for iScale


  std::cout<<" Found: "<<numberOfEdges<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;

  sortEdgesAscending();

  printMsg("Done, found ");
}
template<class PointND>
void RipsPersistence::computeApproximatedEdgeVectorMultiScaleNanoFlann(std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                                    std::vector<EpsilonSampling<typename PointND::Scalar, PointND>>& epsilonSampler,
                                                    std::vector<std::vector<double>>& fitCoordinates,
                                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                                    std::vector<float>& pointCloudBarycenter,
                                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                                    int dim,
                                                    double maxEdgeLength){

  Timer tm{};
  printMsg("Computing exact list of edges in "+std::to_string(dim)+"D, WITH SUBSAMPLING, below "+std::to_string(maxEdgeLength));

  SimplexId nbScales = epsilonIdToFilteredId.size();

  edges_.clear();
  edges_.reserve(10000000);

  if(maxEdgeLength<0){
    maxEdgeLength = std::numeric_limits<double>::max();
  }

  long long int numberOfEdges = 0;


  std::vector<SimplexId> pointsAtScaleStartIdx(nbScales,0);


  SimplexId totalNbPoints=0;
  for (int iScale=0; iScale<nbScales; iScale++){
    totalNbPoints+=epsilonIdToFilteredId[iScale].size();
  }
  std::cout<<" Total nb of points : "<<totalNbPoints<<std::endl;

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif // TTK_ENABLE_OPENMP
  for (int iScale=0; iScale<nbScales-1; iScale++){
    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();
    pointsAtScaleStartIdx[iScale+1] = pointsAtScaleStartIdx[iScale] + nbPointsAtScale;

  }


  for (int iScale=0; iScale<nbScales; iScale++){
    std::cout<<"\nScale "<<iScale<<std::endl;


    bool last_scale = (iScale == nbScales-1);
    SimplexId nbPointsAtScale = epsilonIdToFilteredId[iScale].size();


    double distMax = std::numeric_limits<double>::max();
    distMax = std::min(distMax, (double)maxEdgeLength);

    long long unsigned int addedEdges=0;
    SimplexId addedEdges_monoscale = 0;
    SimplexId addedEdges_multiscale = 0;

    double max_added_dist=0;

    std::cout<<" nb samples: "<<nbPointsAtScale<< " (potentially "<<nbPointsAtScale*(nbPointsAtScale-1)<<" edges )"<<std::endl;


    for(LongSimplexId id_mono_eps_i=0; id_mono_eps_i<nbPointsAtScale; id_mono_eps_i++){

      const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id

      const double* coordinates_i = fitCoordinates[iScale].data();
      const float* mlsCoordinates_i = mls3DCoordinates[iScale].data();

      const double* coordsi = &coordinates_i[dim*id_mono_local_i];
      const SimplexId id_multiscale_i = id_mono_eps_i+pointsAtScaleStartIdx[iScale];

        std::vector<double> queryPoint(dim);
        for(int k=0; k<dim; k++){
          queryPoint[k] = coordinates_i[dim*id_mono_local_i+k];
        }
      // const double queryPoint[6] = {coordinates_i[5*id_mono_local_i+0],
      //                               coordinates_i[5*id_mono_local_i+1],
      //                               coordinates_i[5*id_mono_local_i+2],
      //                               coordinates_i[5*id_mono_local_i+3],
      //                               coordinates_i[5*id_mono_local_i+4]};


      // Query in same scale iScale
      auto tree_same_scale = epsilonSampler[iScale].treeNF();
      std::vector<nanoflann::ResultItem<unsigned int, double>> tree_query_same_scale;
      SimplexId n_matches_same_scale = tree_same_scale->radiusSearch(queryPoint.data(), distMax*distMax, tree_query_same_scale);
      for(auto& match : tree_query_same_scale){

        const double dist_match = match.second;

        const SimplexId id_mono_eps_j = match.first;
        const SimplexId id_mono_local_j = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_j]];

        const double* coordinates_j = fitCoordinates[iScale].data();
        const double* coordsj = &coordinates_j[dim*id_mono_local_j];
        const float* mlsCoordinates_j = mls3DCoordinates[iScale].data();

        const SimplexId id_multiscale_j = id_mono_eps_j+pointsAtScaleStartIdx[iScale];

        if(id_multiscale_i<id_multiscale_j){ // to avoid duplicates

          double dist_comp = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*id_mono_local_i], &mlsCoordinates_j[3*id_mono_local_j], pointCloudBarycenter.data(), dim);
          // std::cout<<" dist match "<<std::sqrt(dist_match)<<"  vs  dist comp "<<dist_comp<<std::endl;

          if(dist_comp <= maxEdgeLength){
            edges_.emplace_back(id_multiscale_i, id_multiscale_j, dist_comp);

            if(dist_comp>max_added_dist){
              max_added_dist=dist_comp;
            }
            numberOfEdges++;
            addedEdges++;
            addedEdges_monoscale++;
          }
        }
      } // end for tree query same scal

      // Query in next scale iScale+1
      if(!last_scale){
        SimplexId last_added = addedEdges;
        auto tree_diff_scale = epsilonSampler[iScale+1].treeNF();
        std::vector<nanoflann::ResultItem<unsigned int, double>> tree_query_diff_scale;
        SimplexId n_matches_diff_scale = tree_diff_scale->radiusSearch(queryPoint.data(), distMax*distMax, tree_query_diff_scale);
        for(auto& match : tree_query_diff_scale){

          const double dist_match = match.second;

          const SimplexId id_mono_eps_j = match.first;
          const SimplexId id_mono_local_j = filteredToLocalId[iScale+1][epsilonIdToFilteredId[iScale+1][id_mono_eps_j]];

          const double* coordinates_j = fitCoordinates[iScale+1].data();
          const double* coordsj = &coordinates_j[dim*id_mono_local_j];
          const float* mlsCoordinates_j = mls3DCoordinates[iScale+1].data();

          const SimplexId id_multiscale_j = id_mono_eps_j+pointsAtScaleStartIdx[iScale+1];


          double dist_comp = distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*id_mono_local_i], &mlsCoordinates_j[3*id_mono_local_j], pointCloudBarycenter.data(), dim);
          // std::cout<<" dist match "<<std::sqrt(dist_match)<<"  vs  dist comp "<<dist_comp<<std::endl;
          if(dist_comp <= maxEdgeLength ){
            edges_.emplace_back(id_multiscale_i, id_multiscale_j, dist_comp);

            if(dist_comp>max_added_dist){
              max_added_dist=dist_comp;
            }
            numberOfEdges++;
            addedEdges++;
            addedEdges_multiscale++;
          }
        } // end for tree query diff scale

      }
    } // end for id_mono_eps_i
    std::cout<<"      added "<<addedEdges_monoscale<<" edges between "<<iScale<<" and "<<iScale<<std::endl;
    std::cout<<"      added "<<addedEdges_multiscale<<" edges between "<<iScale<<" and "<<iScale+1<<std::endl;
    std::cout<<"added "<<addedEdges<<" edges with max dist "<<max_added_dist<<std::endl;

  } // end for iScale


  std::cout<<" Found: "<<numberOfEdges<<" edges in "<<tm.getElapsedTime()<<" s."<<std::endl;

  sortEdgesAscending();

  printMsg("Done, found ");
}



} // namespace ttk
