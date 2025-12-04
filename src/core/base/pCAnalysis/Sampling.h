#pragma once

//ttk includes
#include <Debug.h>
#include <Shuffle.h>

// Ponca includes
#include "DataTypes.h"
#include "Ponca/src/SpatialPartitioning/KdTree/kdTree.h"
#include <Ponca/SpatialPartitioning>

//std includes
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <PCAnalysisUtils.h>
#include <numeric>
#include <random>
#include "../ripsPersistence/nanoflann.hpp"

using namespace ttk;
using namespace pointCloud;


struct PointCloudAdaptor5DEpsSampling
{
    const double* fitCoords;
    const SimplexId n;
    const std::vector<SimplexId> &epsSampling;
    const std::vector<SimplexId> &initialSampling;

    PointCloudAdaptor5DEpsSampling(const double* _fitCoords, const SimplexId _n, const std::vector<SimplexId>& _epsSampling, const std::vector<SimplexId>& _initialSampling) : fitCoords(_fitCoords), n(_n), epsSampling(_epsSampling), initialSampling(_initialSampling){}

    inline size_t kdtree_get_point_count() const
    {
        return n;
    }

    inline double kdtree_get_pt(const SimplexId i, const size_t dim) const
    {
      SimplexId id_local = initialSampling[epsSampling[i]];
        return fitCoords[5*id_local+dim];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }

};  // end of PointCloudAdaptor
struct PointCloudAdaptor5D
{
    const double* fitCoords;
    const SimplexId n;
    const std::vector<SimplexId> &initialSampling;

    PointCloudAdaptor5D(const double* _fitCoords, const SimplexId _n, const std::vector<SimplexId>& _initialSampling) : fitCoords(_fitCoords), n(_n), initialSampling(_initialSampling){}

    inline size_t kdtree_get_point_count() const
    {
        return n;
    }
    inline double kdtree_get_pt(const SimplexId i, const size_t dim) const
    {
      SimplexId id_local = initialSampling[i];
        return fitCoords[5*id_local+dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }

};  // end of PointCloudAdaptor
struct PointCloudAdaptor
{
    // using coord_t = typename Derived::coord_t;

    // const Derived& obj;  //!< A const ref to the data set origin
                         //
    const float* coords;
    const float* normals;
    const SimplexId n;
    std::vector<SimplexId> &sampling;

    /// The constructor that sets the data set source
    PointCloudAdaptor(const float* _coords, const float* _normals, const SimplexId _n, std::vector<SimplexId>& _sampling) : coords(_coords), normals(_normals), n(_n), sampling(_sampling){}

    /// CRTP helper method
    // inline const Derived& derived() const {}

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const
    {
        return n;
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate
    // value, the
    //  "if/else's" are actually solved at compile time.
    inline float kdtree_get_pt(const SimplexId i, const size_t dim) const
    {
      SimplexId idx = sampling[i];
        if (dim == 0)
            return coords[3*idx+0];
        else if (dim == 1)
            return coords[3*idx+1];
        else if (dim == 2)
            return coords[3*idx+2];
        else if (dim == 3)
            return normals[3*idx+0];
        else if (dim == 4)
            return normals[3*idx+1];
        else 
            return normals[3*idx+2];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }

};  // end of PointCloudAdaptor

template <
    class T, class DataSource, typename _DistanceType = T,
    typename IndexType = uint32_t>
struct My_Custom_Metric_Adaptor
{
    using ElementType  = T;
    using DistanceType = _DistanceType;

    const DataSource& data_source;

    double _myParam = 0.0;

    My_Custom_Metric_Adaptor(const DataSource& _data_source, double myParam)
        : data_source(_data_source), _myParam(myParam)
    {
    }

    inline DistanceType evalMetric(
        const T* a, const IndexType b_idx, size_t size) const
    {
        DistanceType result = DistanceType();
        // DistanceType r = DistanceType();
        for (size_t i = 0; i < 3; ++i)
        {
            const DistanceType diff =
                a[i] - data_source.kdtree_get_pt(b_idx, i);
            result += std::pow(diff, 2);
        }
        // r=result;
        for (size_t i = 3; i < 6; ++i)
        {
            const DistanceType diff =
                a[i] - data_source.kdtree_get_pt(b_idx, i);
            result += _myParam*std::pow(diff, 2);
        }
        // std::cout<<" distance "<<a[0]<<"-"<<b_idx<<" from inside "<<result<<"("<<r<<")"<<std::endl;
        return result;
    }

    template <typename U, typename V>
    inline DistanceType accum_dist(const U a, const V b, const size_t i) const
    {
      if(i<3){
        return std::pow((a - b), 2)/(1+_myParam);
      }else{
        return _myParam*std::pow((a - b), 2)/(1+_myParam);
      }
    }
};

template <
    class T, class DataSource, typename _DistanceType = T,
    typename IndexType = uint32_t>
struct My_Custom_Metric_Adaptor5D
{
    using ElementType  = T;
    using DistanceType = _DistanceType;

    const DataSource& data_source;


    My_Custom_Metric_Adaptor5D(const DataSource& _data_source)
        : data_source(_data_source)
    {
    }

    inline DistanceType evalMetric(
        const T* a, const IndexType b_idx, size_t size) const
    {
        DistanceType result = DistanceType();
        // DistanceType r = DistanceType();
        for (size_t i = 0; i < 5; ++i)
        {
            const DistanceType diff =
                a[i] - data_source.kdtree_get_pt(b_idx, i);
            result += std::pow(diff, 2);
        }
        return result;
    }

    template <typename U, typename V>
    inline DistanceType accum_dist(const U a, const V b, const size_t i) const
    {
        return std::pow((a - b), 2);
    }
};

template <int Dim=5>
struct PointCloudAdaptorND
{
    const double* fitCoords;
    const SimplexId n;
    const std::vector<SimplexId> &initialSampling;

    PointCloudAdaptorND(const double* _fitCoords, const SimplexId _n, const std::vector<SimplexId>& _initialSampling) : fitCoords(_fitCoords), n(_n), initialSampling(_initialSampling){}

    inline size_t kdtree_get_point_count() const
    {
        return n;
    }
    inline double kdtree_get_pt(const SimplexId i, const size_t dim) const
    {
      SimplexId id_local = initialSampling[i];
        return fitCoords[Dim*id_local+dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }

};  // end of PointCloudAdaptor
template <int Dim = 5>
struct PointCloudAdaptorNDEpsSampling
{
    const double* fitCoords;
    const SimplexId n;
    const std::vector<SimplexId> &epsSampling;
    const std::vector<SimplexId> &initialSampling;

    PointCloudAdaptorNDEpsSampling(const double* _fitCoords, const SimplexId _n, const std::vector<SimplexId>& _epsSampling, const std::vector<SimplexId>& _initialSampling) : fitCoords(_fitCoords), n(_n), epsSampling(_epsSampling), initialSampling(_initialSampling){}

    inline size_t kdtree_get_point_count() const
    {
        return n;
    }

    inline double kdtree_get_pt(const SimplexId i, const size_t iDim) const
    {
      SimplexId id_local = initialSampling[epsSampling[i]];
        return fitCoords[Dim*id_local+iDim];
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
        return false;
    }
};  // end of PointCloudAdaptor

template <
    class T, class DataSource, int Dim, typename _DistanceType = T,
    typename IndexType = uint32_t>
struct My_Custom_Metric_AdaptorND
{
    using ElementType  = T;
    using DistanceType = _DistanceType;

    const DataSource& data_source;


    My_Custom_Metric_AdaptorND(const DataSource& _data_source)
        : data_source(_data_source)
    {
    }

    inline DistanceType evalMetric(
        const T* a, const IndexType b_idx, size_t size) const
    {
        DistanceType result = DistanceType();
        // DistanceType r = DistanceType();
        for (size_t i = 0; i < Dim; ++i)
        {
            const DistanceType diff =
                a[i] - data_source.kdtree_get_pt(b_idx, i);
            result += std::pow(diff, 2);
        }
        return result;
    }

    template <typename U, typename V>
    inline DistanceType accum_dist(const U a, const V b, const size_t i) const
    {
        return std::pow((a - b), 2);
    }
};

// template <class Scalar>
// class EpsilonSamplingAbstract: public ttk::Debug{
//   private:
//   EpsilonSampling<pointCloud::MyPointND<1>::Scalar, pointCloud::MyPointND<1>> epsilonSampler1;


//   public:
//     EpsilonSamplingAbstract();
// }

template <class Scalar, typename Point>
class EpsilonSampling : public ttk::Debug{

  public:
using KDTreeNFND = nanoflann::KDTreeSingleIndexAdaptor<
  My_Custom_Metric_AdaptorND<double, PointCloudAdaptorND<Point::Dim>, Point::Dim>,
  PointCloudAdaptorND<Point::Dim>, Point::Dim /* dim */ >;

using KDTreeNFNDEpsSampling = nanoflann::KDTreeSingleIndexAdaptor<
  My_Custom_Metric_AdaptorND<double, PointCloudAdaptorNDEpsSampling<Point::Dim>, Point::Dim>,
  PointCloudAdaptorNDEpsSampling<Point::Dim>, Point::Dim /* dim */ >;

  public:
  EpsilonSampling(){
    this->setDebugMsgPrefix("EpsilonSampling");
    this->setDebugLevel(3);
  };

  
void buildSampling(const int _nbPoints, const Scalar* _coordinates, const std::vector<SimplexId>& initial_sampling, std::vector<SimplexId>& outputEpsilonSampling, int seed=0);
int buildKDTree(Ponca::KdTree<Point>& tree, Scalar* _coordinates, std::vector<SimplexId>& initial_sampling, std::vector<SimplexId>& representant, SimplexId _nbPoints);
void checkAndClean(int _nbPoints, Scalar* _coordinates, float* mlsCoords, const float* barycenter, std::vector<SimplexId>& initial_sampling);

const std::vector<std::vector<SimplexId>>& getChildren() const {
  return m_children;
}

const std::vector<SimplexId>& getChildren(SimplexId i) const {
  return m_children[i];
}
const SimplexId size(SimplexId i) const {
  return m_children[i].size();
}

inline SimplexId numberOfSamples() const {
  return m_tree.sample_count();
}

inline SimplexId sample_count() const {
  return m_tree.sample_count();
}

void setNbRepMax(SimplexId nbRepMax){
  m_nbRepMax = nbRepMax;
}

void setEpsilon(Scalar epsilon){
  m_epsilon = epsilon;
}
inline SimplexId max_representant_count() const {
  return m_nbRepMax;
}
inline Scalar epsilon() const {
  return m_epsilon;
}

const std::vector<SimplexId>& samples() const {
  return m_tree.samples();
}

const Ponca::KdTree<Point>& tree() const{
  return m_tree;
}
const KDTreeNFNDEpsSampling* treeNF() const{
  return m_treeNF;
}


std::vector<SimplexId>& representant() {
  return m_representant;
}

void clear_state(){
  m_children.clear();
  m_representant.clear();
}

  private:

  SimplexId m_nbRepMax{-1};
  Scalar m_epsilon{};
  Ponca::KdTree<Point> m_tree{};
  KDTreeNFNDEpsSampling* m_treeNF{};
  std::vector<std::vector<SimplexId>> m_children{};
  std::vector<SimplexId> m_representant{};
};

template <class Scalar, typename Point>
int EpsilonSampling<Scalar, Point>::buildKDTree(Ponca::KdTree<Point>& tree, Scalar* _coordinates, std::vector<SimplexId>& initial_sampling, std::vector<SimplexId>& representant, SimplexId _nbPoints){


  Timer tm{};
  struct VN{
      Scalar* coordinates;

      int n;
      std::vector<SimplexId>& initial_sampling;

      int size() const {
        return n;
      }
  };

  if(representant.size()!=_nbPoints or initial_sampling.size()!=_nbPoints){
    printErr("Size mismatch");
    return -1;
  }

  std::vector<int> sampling{};

  sampling.reserve(_nbPoints);
  for(int i = 0; i<_nbPoints; i++){
    if(representant[i]==i){
      sampling.push_back(i);
    }
  }

  typedef typename Ponca::KdTree<Point>::PointContainer PtContainer;
  typedef typename Ponca::KdTree<Point>::IndexContainer IdxContainer;

  tree.buildWithSampling(VN({_coordinates, _nbPoints, initial_sampling}), IdxContainer(sampling),
      [=](VN buffer, PtContainer &out){
                out.reserve(buffer.size());
                for(int i = 0; i<buffer.size(); i++){
                out.push_back(Point(buffer.coordinates, buffer.initial_sampling[i]));
                // if(buffer.initial_sampling[i]==300697){
                // std::cout<<" point 300697 from epssampling "<<out.back().pos()<<std::endl;
                // }
                }});
  printMsg("kd-tree built, point count "+std::to_string(tree.point_count()) + ", index count "+std::to_string(tree.sample_count()) , 1, tm.getElapsedTime());
  return 0;
}
template<class Scalar, typename Point>
void EpsilonSampling<Scalar, Point>::buildSampling(const int _nbPoints, const Scalar* _coordinates, const std::vector<SimplexId>& initial_sampling, std::vector<SimplexId>& outputEpsilonSampling, int seed){

  int dim = Point::Dim;
  std::cout<<"build sampling, dimension "<<dim<<std::endl;
  // printMsg("Start Epsilon Sampling, dimension "+std::to_string(dim));
  if(_nbPoints==0){
    printWrn("No points -- skipping");
    return;
  }

  // std::cout<<"epsilon sampling"<<std::endl;
  Timer tm{};

  std::mt19937 random_engine{};
  random_engine.seed(seed);
  auto& representant = m_representant;
  representant.resize(_nbPoints);
  // std::vector<SimplexId> samplingNF(_nbPoints);
  std::vector<SimplexId> shuffledIndices(_nbPoints);
  std::iota(shuffledIndices.begin(), shuffledIndices.end(), 0);
  ttk::shuffle(shuffledIndices, random_engine);

  Scalar inc = m_epsilon/10;
  SimplexId nb_inc=0;

  Scalar disk_radius = m_epsilon*0.8; 

  PointCloudAdaptorND<Point::Dim>* adaptor = new PointCloudAdaptorND<Point::Dim>(_coordinates, _nbPoints, initial_sampling);
  KDTreeNFND* ptree = new KDTreeNFND(Point::Dim, *adaptor, {10});

  m_children.clear();
  while(m_children.empty() or (m_children.size()>m_nbRepMax and m_nbRepMax>0)){

    disk_radius+=inc;
    nb_inc++;
    if(nb_inc%10==0){
      inc*=10;
    }

    std::iota(representant.begin(), representant.end(), 0);

    // Ponca::KdTree<Point> fullTree{};
    // m_tree.set_min_cell_size(1);
    // buildKDTree(m_tree, _coordinates,initial_sampling,  representant, _nbPoints);

    // std::iota(samplingNF.begin(), samplingNF.end(), 0);





    // std::cout<<"indices size "<<shuffledIndices.size()<<" nbpoints "<<_nbPoints<<std::endl;

    m_children.clear();
    for(int iPoint=0; iPoint<_nbPoints; iPoint++){
      int iShuffledPoint = shuffledIndices[iPoint];
      if(representant[iShuffledPoint]==iShuffledPoint){
        std::vector<SimplexId> children{};
        children.push_back(iShuffledPoint);


        // ponca query 
        // auto tree_query = m_tree.range_neighbors(iShuffledPoint, disk_radius);
        // SimplexId nbNeighbors=0;
        // SimplexId poncaMatches=0;
        // for(SimplexId iNeigh : tree_query){
        //   if(representant[iNeigh]==iNeigh){
        //     nbNeighbors++;
        //     representant[iNeigh] = iShuffledPoint;
        //     poncaMatches++;
        //     children.push_back(iNeigh);
        //   }
        // }
        // std::cout<<"ponca query "<<poncaMatches<<" matches"<<std::endl;

        SimplexId id_local = initial_sampling[iShuffledPoint];
        // const Scalar queryPoint[6] = {_coordinates[5*id_local+0],
        //   _coordinates[5*id_local+1],
        //   _coordinates[5*id_local+2],
        //   _coordinates[5*id_local+3],
        //   _coordinates[5*id_local+4]};
        std::vector<Scalar> queryPoint(dim);
        for(int k=0; k<dim; k++){
          queryPoint[k] = _coordinates[dim*id_local+k];
        }

        std::vector<nanoflann::ResultItem<unsigned int, Scalar>> ret_matches;
        SimplexId nMatches = ptree->radiusSearch(queryPoint.data(), disk_radius*disk_radius, ret_matches);
        SimplexId NFmatches=0;
        for(auto match : ret_matches){
          SimplexId iNeigh = match.first;
          if(iNeigh!=iShuffledPoint and representant[iNeigh]==iNeigh){
            representant[iNeigh] = iShuffledPoint;
            children.push_back(iNeigh);
          }
        }

        m_children.emplace_back();
        m_children.back() = std::move(children);
      }
    }
  std::cout<<" epsilon "<<m_epsilon<<" "<<disk_radius<<std::endl;
    std::cout<<" nb of representants: "<<m_children.size()<<std::endl;
  }
  delete adaptor;
  delete ptree;

  // new subsampling
  SimplexId nbPointsSampling=m_children.size();
  outputEpsilonSampling.resize(nbPointsSampling);
  for(SimplexId ieps=0; ieps<nbPointsSampling; ieps++){
    outputEpsilonSampling[ieps] = m_children[ieps][0];
  }

  PointCloudAdaptorNDEpsSampling<Point::Dim>* adaptorEpsSampling = new PointCloudAdaptorNDEpsSampling<Point::Dim>(_coordinates, nbPointsSampling, outputEpsilonSampling, initial_sampling);
  m_treeNF = new KDTreeNFNDEpsSampling(Point::Dim, *adaptorEpsSampling, {10});




  // printMsg("Epsilon Sampling done", 1, tm.getElapsedTime(), threadNumber_);
  // std::cout<<" sample count : "<<adaptorEpsSampling->kdtree_get_point_count()<<" vs children size "<<m_children.size()<<std::endl;
  // std::cout<<" objective was "<<m_nbRepMax<<" reps"<<std::endl;
  // std::cout<<" epsilon "<<m_epsilon<<" "<<disk_radius<<std::endl;
  this->setEpsilon(disk_radius);

  std::cout<<" end build sampling"<<std::endl;
}

template<class Scalar, typename Point>
void EpsilonSampling<Scalar, Point>::checkAndClean(int _nbPoints, Scalar* _coordinates, float* mlsCoords, const float* barycenter, std::vector<SimplexId>& initial_sampling){

  printMsg("Check Epsilon Sampling, dimension "+std::to_string(Point::Dim));
  // std::cout<<"epsilon sampling"<<std::endl;
  Timer tm{};

  Scalar disk_radius = m_epsilon;

  SimplexId nbRepresentants = m_children.size();

  std::cout<<"checking query at epsilon "<<m_epsilon<<std::endl;
  for(int iPoint=0; iPoint<nbRepresentants; iPoint++){

    int filteredId = m_children[iPoint][0];

    auto tree_query = m_tree.range_neighbors(filteredId, disk_radius);
    SimplexId nMatches = 0;
    for(int iNeigh : tree_query){
      nMatches++;
    }
    if(nMatches>0){

      std::cout<<"WARNING: filteredid "<<filteredId<<" has "<<nMatches<<" matches"<<std::endl;
      if(filteredId==207391){
        for(int iNeigh : tree_query){
          std::cout<<" "<<iNeigh;
        }
        std::cout<<std::endl;
      }

    }
  }

  // std::cout<<"checking distance"<<std::endl;
  // for(int i=0; i<nbRepresentants; i++){
  //   if(m_children[i][0]==207391){
  //     SimplexId nbNeighbors=0;
  //     SimplexId I = initial_sampling[m_children[i][0]];
  //     for(int j=0; j<nbRepresentants; j++){
  //       SimplexId J = initial_sampling[m_children[j][0]];
  //       double dist = ttk::pointCloud::compute5DDistanceHelper(_coordinates+5*I, _coordinates+5*J, mlsCoords+3*I, mlsCoords+3*J, barycenter);
  //       if(dist < m_epsilon){
  //         std::cout<<" dist "<<dist<<" btw "<<m_children[i][0]<<"("<<I<<") and "<<m_children[j][0]<<"("<<J<<")"<<std::endl;
  //         nbNeighbors++;
  //       }
  //     }
  //     std::cout<<" nb neighbors "<<nbNeighbors<<std::endl;
  //   }
  // }

printMsg("Checking Epsilon Sampling done", 1, tm.getElapsedTime(), threadNumber_);

}


template <class Scalar, typename Point>
class ScaleSpaceSampling : public ttk::Debug{


  public:
  ScaleSpaceSampling(){
    this->setDebugMsgPrefix("ScaleSpaceSampling");
  };

  using KDTreeNF5D = nanoflann::KDTreeSingleIndexAdaptor<
    My_Custom_Metric_Adaptor5D<double, PointCloudAdaptor5D>,
    PointCloudAdaptor5D, 5 /* dim */ >;

  using KDTreeNF = nanoflann::KDTreeSingleIndexAdaptor<
    My_Custom_Metric_Adaptor<float, PointCloudAdaptor>,
    PointCloudAdaptor, 6 /* dim */ >;

  // access
  public:
  void setSamplingFactor(double factor){
    m_factor = factor;
  }
  double samplingFactor(){
    return m_factor;
  }
  inline KDTreeNF* kdTreeNF(int lvl){
    return m_treesNF[lvl];
  }
  inline KDTreeNF* kdTreeNF(){
    return kdTreeNF(m_curLevel);
  }
  inline Ponca::KdTree<Point>& kdTree(int lvl){
    return lvl >= 0 ? m_trees[lvl] : m_tree_original;
  }
  inline Ponca::KdTree<Point>& kdTree(){
    return kdTree(m_curLevel);
  }
  inline Ponca::KnnGraph<Point>& knnGraph(){
    return *m_KnnGraph;
  }
  inline Scalar scale(){
    return scale(m_curLevel);
  }
  inline Scalar scale(int lvl){
    return m_scales[lvl];
  }
  inline const std::vector<Scalar>& scales(){
    return m_scales;
  }
  inline int level(int i){
    return m_levels[i];
  }
  // inline void setKDTree(Ponca::KdTree<MyPoint>* tree){
  //   m_tree = tree;
  // }
  inline int scale_count() const{
    return m_scales.size();
  }
  inline const std::vector<SimplexId>& sampling(int lvl){
    std::cout<<"size of m_sampling "<<m_sampling.size()<<std::endl;
    std::cout<<"sending scale index "<<lvl<<std::endl;
    return m_sampling[lvl];
  }
  inline bool is_sub_sampled(int iPoint, int level){
    return (level <= m_levels[iPoint]);
  };

  inline KDTreeNF* getPTest(){
  return pTest;
  }

  inline bool isUseKnnGraph() const{
    return useKnnGraph_;
  }
  inline void setUseKnnGraph(bool data){
    useKnnGraph_=data;
  }

  void printTreesInfosNF(){
    std::cout<<"======= NF TREES INFOS ======"<<std::endl;
    std::cout<<" "<<m_treesNF.size()<<" trees in memory"<<std::endl;
    for(auto p : m_treesNF){
    std::cout<<p<<" size "<<p->size_<<std::endl;
    }
    std::cout<<"============================="<<std::endl;
  }

  void print(){
    this->printMsg("ScaleSpaceSampling:");
    std::cout<<" Scales : ";
    for(auto s : m_scales){
      std::cout<<" "<<s;
    }
    std::cout<<"\n Poisson Disk levels:\n";
    for(int i=0; i< m_scales.size(); i++){
      std::cout<<" level "<<i<<": ";
      for(auto l : m_levels){
        if(l>=i) std::cout<<" "<<l;
        else std::cout<<" _";
      }
      std::cout<<std::endl;
    }
      std::cout<<std::endl;
  }

  public:

  public:
  void sampleScaleLog(Scalar _scaleMin, Scalar _scaleMax, int _nbSamples, double base = 10);
  void poissonDiskSampling(int _nbPoints, Scalar* _coordinates, Scalar* _normals, double ratio=1);
  void poissonDiskSamplingNANOFLANN(int _nbPoints, Scalar* _coordinates, Scalar* _normals, double ratio);
  void poissonDiskSamplingND(int _nbPoints, Scalar* _coordinates);
  void poissonDiskSamplingBucketIndexing(std::vector<std::vector<Scalar>>& coordinates,
                                         std::vector<std::vector<SimplexId>>& internalSampling,
                                         std::vector<std::vector<SimplexId>>& externalSampling,
                                         int startBucketId=-1,
                                         int stopBucketId=-1);
  void poissonDiskSamplingBucketIndexingNF(std::vector<std::vector<Scalar>>& coordinates,
                                         std::vector<std::vector<SimplexId>>& internalSampling,
                                         std::vector<std::vector<SimplexId>>& externalSampling,
                                         int startBucketId=-1,
                                         int stopBucketId=-1);

  void buildKDTree(int level, Scalar* _coordinates, Scalar* _normals, int _nbPoints, double ratio);
  int findMinScale(SimplexId maxPoints);
  void buildKDTreeND(int level, Scalar* _coordinates, int _nbPoints);
  void buildKDTreeBucketIndexing(int level,
                                SimplexId _nbPoints,
                                std::vector<SimplexId>& _bucketStartId,
                                std::vector<std::vector<Scalar>>& _coordinates,
                                std::vector<std::vector<SimplexId>>& _internalSampling,
                                std::vector<std::vector<SimplexId>>& _externalSampling,
                                int startBucketId=0);
  void filterScaleByDensity(int level, Scalar* _coordinates, int _nbPoints, Scalar* density, Scalar min_density, std::vector<SimplexId>& densitySampling, int* isStable);
  void filterUnstablePoints(int level, Scalar* _coordinates, int _nbPoints, std::vector<SimplexId>& densitySampling, int* isStable);
  Scalar estimateDensity(Scalar* coordinates, Scalar* normals, SimplexId nbPoints, int percent_tries, SimplexId nbNeighbors);
  void computeDensityField(int iScale, SimplexId k, Scalar* coords, int dim, SimplexId nbPoints, std::vector<SimplexId>& sampling_indexes, Scalar* output_density);
  void computeDensityField(int iScale, SimplexId k, std::vector<std::vector<double>>& distMat, Scalar radius, Scalar* output_density);

  void buildKnnGraph();

  void printDist(SimplexId i, SimplexId j, Scalar* c, Scalar* n){
      std::cout<<"       coordsI"<<c[3*i]<<" "<<c[3*i+1]<<" "<<c[3*i+2];
      std::cout<<" "<<n[3*i]<<" "<<n[3*i+1]<<" "<<n[3*i+2]<<std::endl;
      std::cout<<"       coordsJ"<<c[3*j]<<" "<<c[3*j+1]<<" "<<c[3*j+2];
      std::cout<<" "<<n[3*j]<<" "<<n[3*j+1]<<" "<<n[3*j+2]<<std::endl;
    double dist = std::sqrt( (c[3*i+0]-c[3*j+0])*(c[3*i+0]-c[3*j+0]) +
                  (c[3*i+1]-c[3*j+1])*(c[3*i+1]-c[3*j+1]) +
                  (c[3*i+2]-c[3*j+2])*(c[3*i+2]-c[3*j+2]) );
      std::cout<<"       dist "<<dist<<std::endl;

  }

  void clear(){
    for(int i=0; i<scale_count();i++){
      if(m_pointCloudAdaptors[i]){
        delete m_pointCloudAdaptors[i];
      }
      if(m_treesNF[i]){
        delete m_treesNF[i];
      }
    }
    delete m_treeNF_original;
  }

  private:

  // int m_scaleCount{0};
  // double m_scaleMin{1}, m_scaleMax{1};
  int m_curLevel{0};
  std::vector<int> m_levels{};
  std::vector<Scalar> m_scales{};
  double m_factor{0.1};
  Ponca::KdTree<Point>* m_curTree{nullptr};
  std::vector<Ponca::KdTree<Point>> m_trees{};
  Ponca::KdTree<Point> m_tree_original{};
  std::vector<KDTreeNF*> m_treesNF{};
  KDTreeNF* m_treeNF_original{};
  std::vector<std::vector<SimplexId>> m_sampling{};
  std::vector<PointCloudAdaptor*> m_pointCloudAdaptors{};
  KDTreeNF* pTest{};

  bool useKnnGraph_{false};
  Ponca::KnnGraph<Point>* m_KnnGraph{};
};

template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::sampleScaleLog(Scalar _scaleMin, Scalar _scaleMax, int _nbSamples, double base){

      std::stringstream ss;
      ss << "Sampling "<<_nbSamples<<" scales btw "<<_scaleMin<<" and "<<_scaleMax<<" (base "<<base<<" )";
      printMsg(ss.str());

      m_scales.clear();
      m_scales.resize(_nbSamples);
      m_scales[0] = _scaleMin;
      if(_nbSamples==1){
        return;
      }

      if(std::abs(base)<1e-3) {
        // uniform sampling
        Scalar delta0 = (_scaleMax-_scaleMin)/(_nbSamples-1);
        for(int i=1; i<_nbSamples; ++i) {
          m_scales[i] = _scaleMin + delta0*i;
        }

      } else {
        // log sampling
        double lMin = log(_scaleMin)/log(base);
        double lMax = log(_scaleMax)/log(base);
        Scalar delta = (lMax-lMin) / (_nbSamples-1);
        for(int i=0; i<_nbSamples; ++i) {
          m_scales[i] = std::pow(base, lMin + delta*i);
        }
      }
}


template<class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::poissonDiskSamplingNANOFLANN(int _nbPoints, Scalar* _coordinates, Scalar* _normals, double ratio){

  Timer tm{};

    printMsg("Start Poisson Disk Sampling NF");

  if(m_scales.empty()){
    printErr("No scales detected");
    return;
  }
  // m_treesNF.resize(scale_count());

  m_levels.resize(_nbPoints, 0);
  std::vector<bool> free(_nbPoints, true);

  // shuffle indices for random sampling
  std::vector<int> shuffledIndices(_nbPoints);
  std::iota(shuffledIndices.begin(), shuffledIndices.end(), 0);

  // shuffle them using the seed 0
  std::mt19937 random_engine{};
  random_engine.seed(0);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  ttk::shuffle(shuffledIndices, random_engine);

  m_pointCloudAdaptors.resize(scale_count());
  m_sampling.resize(m_scales.size());
  m_sampling[0].resize(_nbPoints);
  std::iota(m_sampling[0].begin(), m_sampling[0].end(),0);

  // KDTreeNF* ptreeTest{};

  for(int iScale = 1; iScale<m_scales.size(); iScale++){
    Scalar disk_radius = m_factor * m_scales[iScale];

    // build kd tree of previous level

    m_pointCloudAdaptors[iScale-1] = new PointCloudAdaptor(_coordinates, _normals, m_sampling[iScale-1].size(), m_sampling[iScale-1]);



    KDTreeNF* ptree = new KDTreeNF(6, *m_pointCloudAdaptors[iScale-1], {10}, ratio);
    printMsg("kd-tree (NF) built at level "+std::to_string(iScale-1) + ", point count "+std::to_string(m_sampling[iScale-1].size()) , 1, tm.getElapsedTime());
    m_treesNF.push_back(ptree);

  
    // random sampling on available points
    for(int iPoint=0; iPoint<_nbPoints; iPoint++){
      int iShuffledPoint = shuffledIndices[iPoint];
      if(free[iShuffledPoint]){
        m_levels[iShuffledPoint] = iScale; // sample this point
      const Scalar queryPoint[6] = {_coordinates[3*iShuffledPoint+0],
                                          _coordinates[3*iShuffledPoint+1],
                                          _coordinates[3*iShuffledPoint+2],
                                          _normals[3*iShuffledPoint+0],
                                          _normals[3*iShuffledPoint+1],
                                          _normals[3*iShuffledPoint+2]};
        std::vector<nanoflann::ResultItem<unsigned int, Scalar>> ret_matches;
        SimplexId nMatches = ptree->radiusSearch(queryPoint, disk_radius*disk_radius, ret_matches);
        // std::cout<<" query "<<disk_radius<<" on "<<iShuffledPoint<<" : "<<nMatches-1<<" matches"<<std::endl;
        for(auto match : ret_matches){
          SimplexId idx = m_sampling[iScale-1][match.first];
            if(idx!=iShuffledPoint){
              // std::cout<<" "<<idx;
                // printDist(iShuffledPoint, idx, _coordinates, _normals);
              free[idx] = false;
            }
        }
          // std::cout<<std::endl;
      }
    }
    // m_sampling[iScale].reserve(_nbPoints);
    for(int i = 0; i<_nbPoints; i++){
      if(is_sub_sampled(i,iScale)){
        m_sampling[iScale].push_back(i);
      }
    }

    // printTreesInfosNF();
  }

  // build tree for last level
    m_pointCloudAdaptors[scale_count()-1] = new PointCloudAdaptor(_coordinates, _normals, m_sampling[scale_count()-1].size(), m_sampling[scale_count()-1]);
    KDTreeNF* ptree = new KDTreeNF(6, *m_pointCloudAdaptors[scale_count()-1], {10}, ratio);
    printMsg("kd-tree (NF) built at level "+std::to_string(scale_count()-1) + ", point count "+std::to_string(m_sampling[scale_count()-1].size()) , 1, tm.getElapsedTime());
    m_treesNF.push_back(ptree);
    // printTreesInfosNF();

      // std::cout<<"test loop"<<std::endl;
  // SimplexId iShuffledPoint = 1264;
      // const float queryPoint[6] = {_coordinates[3*iShuffledPoint+0],
      //                                     _coordinates[3*iShuffledPoint+1],
      //                                     _coordinates[3*iShuffledPoint+2],
      //                                     _normals[3*iShuffledPoint+0],
      //                                     _normals[3*iShuffledPoint+1],
      //                                     _normals[3*iShuffledPoint+2]};
      // double dist = 0.305019;
      // for(int i =0; i<10; i++){
      //   std::cout<<"ptreeTest size "<<ptreeTest->size_<<std::endl;
      //   std::vector<nanoflann::ResultItem<unsigned int, float>> ret_matches;
      //   // ptreeTest->buildIndex();
      //   SimplexId nMatches = ptreeTest->radiusSearch(queryPoint, dist*dist, ret_matches);
      //   std::cout<<"ptreeTest size "<<ptreeTest->size_<<std::endl;
      //   std::cout<<" found "<<nMatches<<" matches"<<std::endl;

      // }

  printMsg("Poisson Disk Sampling done", 1, tm.getElapsedTime(), threadNumber_);
}
template<class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::poissonDiskSamplingND(int _nbPoints, Scalar* _coordinates){
 
  
  printMsg("Start Poisson Disk Sampling Ponca, dimension "+std::to_string(Point::Dim));

  Timer tm{};

  if(m_scales.empty()){
    printErr("No scales detected");
    return;
  }
  m_trees.resize(m_scales.size());

  // if(!m_tree){
  //   printErr("KD-Tree not built");
  //   return;
  // }

  m_levels.resize(_nbPoints, -1);
  std::vector<bool> free(_nbPoints, true);

  // shuffle indices for random sampling
  // std::vector<int> shuffledIndices(_nbPoints);
  // std::iota(shuffledIndices.begin(), shuffledIndices.end(), 0);

  // shuffle them using the seed 0
  std::mt19937 random_engine{};
  random_engine.seed(0);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  // ttk::shuffle(shuffledIndices, random_engine);



  for(int iScale = 0; iScale<m_scales.size(); iScale++){
    Scalar disk_radius = m_factor * m_scales[iScale];

    // build kd tree of previous level
    buildKDTreeND(iScale-1, _coordinates,_nbPoints);

    //get points at this level
    SimplexId localNbPoints = kdTree(iScale-1).sample_count();

    std::vector<SimplexId> indices(kdTree(iScale-1).samples());
    ttk::shuffle(indices, random_engine);

    std::cout<<"indices size "<<indices.size()<<std::endl;

    // double progress=0;
    // random sampling on available points
    for(int iPoint=0; iPoint<localNbPoints; iPoint++){
      int iShuffledPoint = indices[iPoint];
      // std::cout<<" ishuf "<<iShuffledPoint<<std::endl;
      if(free[iShuffledPoint]){
        // std::cout<<"free"<<std::endl;
        m_levels[iShuffledPoint] = iScale; // sample this point
        auto tree_query = kdTree(iScale-1).range_neighbors(iShuffledPoint, disk_radius);
        // SimplexId nmatches=0;
        // for(int iNeigh : tree_query){
        //   nmatches++;
        // }
        // std::cout<<" query "<<disk_radius<<" on "<<iShuffledPoint<<" : "<<nmatches<<" matches"<<std::endl;
        for(int iNeigh : tree_query){
          // std::cout<<" "<<iNeigh;
          // std::cout<<" neighbor "<<iNeigh<<std::endl;
            // std::cout<<"  PONCA match "<<iNeigh<<"  ("<<iShuffledPoint<<","<<disk_radius<<")"<<std::endl;
            // printDist(iShuffledPoint, iNeigh, _coordinates, _normals);
          free[iNeigh] = false;
        }
        // std::cout<<std::endl;
        // progress += 1/_nbPoints;
        // printMsg("Poisson disk sampling ", progress, 0, threadNumber_, ttk::debug::LineMode::REPLACE);
      }
    }

    // print();
  }

  // build tree for last level
  buildKDTreeND(m_scales.size()-1, _coordinates, _nbPoints);

  // {
  //   int i=0;
  //   for(auto t : m_trees){
  //     std::cout<<"tree "<<++i<<" "<<t.point_count()<<" "<<t.index_count()<<std::endl;
  //   }
  // }

  printMsg("Poisson Disk Sampling done", 1, tm.getElapsedTime(), threadNumber_);
}

template<class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::poissonDiskSamplingBucketIndexingNF(std::vector<std::vector<Scalar>>& coordinates,
                                                                          std::vector<std::vector<SimplexId>>& internalSampling,
                                                                          std::vector<std::vector<SimplexId>>& externalSampling,
                                                                          int startBucketId,
                                                                          int stopBucketId){
  
  printMsg("Start Poisson Disk Sampling Ponca Bucket indexing");


  Timer tm{};

  if(m_scales.empty()){
    printErr("No scales detected");
    return;
  }
  m_trees.resize(m_scales.size());

  SimplexId nbBuckets = coordinates.size();
  if(startBucketId < 0 or startBucketId>=nbBuckets){
    startBucketId=0;
  }
  if(stopBucketId < 0 or stopBucketId>=nbBuckets){
    stopBucketId = nbBuckets-1;
  }
  std::vector<SimplexId> bucketStartId(nbBuckets,0);
  SimplexId nbTotalPoints=0;
  for(int iBucket = 0; iBucket<nbBuckets; iBucket++){
    bucketStartId[iBucket]=nbTotalPoints;
    nbTotalPoints += internalSampling[iBucket].size();
  }
  SimplexId nbPoints=0;
  for(int iBucket = startBucketId; iBucket<=stopBucketId; iBucket++){
    nbPoints += internalSampling[iBucket].size();
  }
  std::cout<<" poisson disk for buckets btw  "<<startBucketId<<" and "<<stopBucketId<<" with "<<nbPoints<<" points "<<std::endl;

  m_levels.resize(nbPoints, 0);
  std::vector<bool> free(nbPoints, true);

  // shuffle indices for random sampling
  std::vector<SimplexId> shuffledIndices(nbPoints);
  std::iota(shuffledIndices.begin(), shuffledIndices.end(), 0);

  // shuffle them using the seed 0
  std::mt19937 random_engine{};
  random_engine.seed(0);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  ttk::shuffle(shuffledIndices, random_engine);

  for(int iScale = 1; iScale<m_scales.size(); iScale++){
    Scalar disk_radius = m_factor * m_scales[iScale-1];

    kdTree(iScale-1).set_min_cell_size(16);
    // build kd tree of previous level
    buildKDTreeBucketIndexing(iScale-1, nbPoints, bucketStartId, coordinates, internalSampling, externalSampling, startBucketId);


    // random sampling on available points
    for(int iPoint=0; iPoint<nbPoints; iPoint++){
      SimplexId iShuffledPoint = shuffledIndices[iPoint];
      if(free[iShuffledPoint]){
        m_levels[iShuffledPoint] = iScale; // sample this point
        auto tree_query = kdTree(iScale-1).range_neighbors(iShuffledPoint, disk_radius);
        for(int iNeigh : tree_query){
          free[iNeigh] = false;
        }
      }
    }
  }

  // build tree for last level
  kdTree(m_scales.size()-1).set_min_cell_size(16);
  buildKDTreeBucketIndexing(m_scales.size()-1, nbPoints, bucketStartId, coordinates, internalSampling, externalSampling, startBucketId);

  printMsg("Poisson Disk Sampling done", 1, tm.getElapsedTime(), threadNumber_);
}
template<class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::poissonDiskSamplingBucketIndexing(std::vector<std::vector<Scalar>>& coordinates,
                                                                          std::vector<std::vector<SimplexId>>& internalSampling,
                                                                          std::vector<std::vector<SimplexId>>& externalSampling,
                                                                          int startBucketId,
                                                                          int stopBucketId){
  
  printMsg("Start Poisson Disk Sampling Ponca Bucket indexing");


  Timer tm{};

  if(m_scales.empty()){
    printErr("No scales detected");
    return;
  }
  m_trees.resize(m_scales.size());

  SimplexId nbBuckets = coordinates.size();
  if(startBucketId < 0 or startBucketId>=nbBuckets){
    startBucketId=0;
  }
  if(stopBucketId < 0 or stopBucketId>=nbBuckets){
    stopBucketId = nbBuckets-1;
  }
  std::vector<SimplexId> bucketStartId(nbBuckets,0);
  SimplexId nbTotalPoints=0;
  for(int iBucket = 0; iBucket<nbBuckets; iBucket++){
    bucketStartId[iBucket]=nbTotalPoints;
    nbTotalPoints += internalSampling[iBucket].size();
  }
  SimplexId nbPoints=0;
  for(int iBucket = startBucketId; iBucket<=stopBucketId; iBucket++){
    nbPoints += internalSampling[iBucket].size();
  }
  std::cout<<" poisson disk for buckets btw  "<<startBucketId<<" and "<<stopBucketId<<" with "<<nbPoints<<" points "<<std::endl;

  m_levels.resize(nbPoints, 0);
  std::vector<bool> free(nbPoints, true);

  // shuffle indices for random sampling
  std::vector<SimplexId> shuffledIndices(nbPoints);
  std::iota(shuffledIndices.begin(), shuffledIndices.end(), 0);

  // shuffle them using the seed 0
  std::mt19937 random_engine{};
  random_engine.seed(0);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  ttk::shuffle(shuffledIndices, random_engine);

  for(int iScale = 1; iScale<m_scales.size(); iScale++){
    Scalar disk_radius = m_factor * m_scales[iScale-1];

    kdTree(iScale-1).set_min_cell_size(16);
    // build kd tree of previous level
    buildKDTreeBucketIndexing(iScale-1, nbPoints, bucketStartId, coordinates, internalSampling, externalSampling, startBucketId);


    // std::cout<<"there1"<<std::endl;
    // random sampling on available points
    for(int iPoint=0; iPoint<nbPoints; iPoint++){
      SimplexId iShuffledPoint = shuffledIndices[iPoint];
      // std::cout<<" ishuf "<<iShuffledPoint<<std::endl;
      // std::cout<<" disk radius "<<disk_radius<<std::endl;
      // if(startBucketId==0 and iShuffledPoint==0){
      //   std::cout<<" iShuffledPoint "<<iShuffledPoint<<std::endl;
      // }
      // if(startBucketId==0 and iShuffledPoint==0){
      //   std::cout<<" is free ? "<<free[iShuffledPoint]<<std::endl;
      // }
      if(free[iShuffledPoint]){
        // std::cout<<"free"<<std::endl;
        m_levels[iShuffledPoint] = iScale; // sample this point
        // std::cout<<"here1"<<std::endl;
        auto tree_query = kdTree(iScale-1).range_neighbors(iShuffledPoint, disk_radius);
        // std::cout<<" kdtree samples "<<kdTree(iScale-1).sample_count()<<" and points "<<kdTree(iScale-1).point_count()<<std::endl;
        // std::cout<<" kd tree points "<<std::endl;
        // for(LongSimplexId ii=0; ii<kdTree(iScale-1).points().size(); ii++){
        //   auto p = kdTree(iScale-1).points()[ii];
        //   std::cout<<"ii "<<ii<<" :  p "<<p.pos()<<std::endl;
        // }
        // std::cout<<" kd tree samples "<<std::endl;
        // for(LongSimplexId ii=0; ii<kdTree(iScale-1).samples().size(); ii++){
        //   auto p = kdTree(iScale-1).samples()[ii];
        //   std::cout<<"ii "<<ii<<" :  p "<<p<<std::endl;
        // }
        // std::cout<<"here2"<<std::endl;
        SimplexId nmatches=0;
        // std::cout<<"here3"<<std::endl;
        // std::cout<<" tree query "<<*tree_query.begin()<<std::endl;
        for(int iNeigh : tree_query){
      // if((startBucketId==0 and iShuffledPoint==0) or (startBucketId==0 and iNeigh==0)){
      //     std::cout<<"\n disk radius "<<disk_radius<<" iShuffled "<<iShuffledPoint<<", neighbor "<<iNeigh<<std::endl;
      // }
          nmatches++;
        }
        // std::cout<<"here4"<<std::endl;
        // std::cout<<" query "<<disk_radius<<" on "<<iShuffledPoint<<" : "<<nmatches<<" matches"<<std::endl;
        for(int iNeigh : tree_query){
          // std::cout<<" "<<iNeigh;
          // std::cout<<" neighbor "<<iNeigh<<std::endl;
          //   std::cout<<"  PONCA match "<<iNeigh<<"  ("<<iShuffledPoint<<","<<disk_radius<<")"<<std::endl;
            // printDist(iShuffledPoint, iNeigh, _coordinates, _normals);
          free[iNeigh] = false;
        }
        // std::cout<<std::endl;
      }
    }
    // std::cout<<"there2"<<std::endl;
    // random sampling on available points

    // print();
  }
  // std::cout<<"there"<<std::endl;

  // build tree for last level
  kdTree(m_scales.size()-1).set_min_cell_size(16);
  buildKDTreeBucketIndexing(m_scales.size()-1, nbPoints, bucketStartId, coordinates, internalSampling, externalSampling, startBucketId);

  // {
  //   int i=0;
  //   for(auto t : m_trees){
  //     std::cout<<"tree "<<++i<<" "<<t.point_count()<<" "<<t.index_count()<<std::endl;
  //   }
  // }

  printMsg("Poisson Disk Sampling done", 1, tm.getElapsedTime(), threadNumber_);
}
template<class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::poissonDiskSampling(int _nbPoints, Scalar* _coordinates, Scalar* _normals, double ratio){
  
  printMsg("Start Poisson Disk Sampling Ponca");

  Timer tm{};

  if(m_scales.empty()){
    printErr("No scales detected");
    return;
  }
  m_trees.resize(m_scales.size());

  // if(!m_tree){
  //   printErr("KD-Tree not built");
  //   return;
  // }

  m_levels.resize(_nbPoints, -1);
  std::vector<bool> free(_nbPoints, true);

  // shuffle indices for random sampling
  std::vector<int> shuffledIndices(_nbPoints);
  std::iota(shuffledIndices.begin(), shuffledIndices.end(), 0);

  // shuffle them using the seed 0
  std::mt19937 random_engine{};
  random_engine.seed(0);
  // use the Fisher-Yates algorithm instead of std::shuffle, whose
  // results are platform-dependent
  ttk::shuffle(shuffledIndices, random_engine);

  for(int iScale = 0; iScale<m_scales.size(); iScale++){
    Scalar disk_radius = m_factor * m_scales[iScale];

    // build kd tree of previous level
    buildKDTree(iScale-1, _coordinates, _normals, _nbPoints, ratio);


    // random sampling on available points
    for(int iPoint=0; iPoint<_nbPoints; iPoint++){
      int iShuffledPoint = shuffledIndices[iPoint];
      // std::cout<<" ishuf "<<iShuffledPoint<<std::endl;
      if(free[iShuffledPoint]){
        // std::cout<<"free"<<std::endl;
        m_levels[iShuffledPoint] = iScale; // sample this point
        auto tree_query = kdTree(iScale-1).range_neighbors(iShuffledPoint, disk_radius);
        SimplexId nmatches=0;
        for(int iNeigh : tree_query){
          nmatches++;
        }
        // std::cout<<" query "<<disk_radius<<" on "<<iShuffledPoint<<" : "<<nmatches<<" matches"<<std::endl;
        for(int iNeigh : tree_query){
          // std::cout<<" "<<iNeigh;
          // std::cout<<" neighbor "<<iNeigh<<std::endl;
            // std::cout<<"  PONCA match "<<iNeigh<<"  ("<<iShuffledPoint<<","<<disk_radius<<")"<<std::endl;
            // printDist(iShuffledPoint, iNeigh, _coordinates, _normals);
          free[iNeigh] = false;
        }
        // std::cout<<std::endl;
      }
    }

    // print();
  }

  // build tree for last level
  buildKDTree(m_scales.size()-1, _coordinates, _normals, _nbPoints, ratio);

  // {
  //   int i=0;
  //   for(auto t : m_trees){
  //     std::cout<<"tree "<<++i<<" "<<t.point_count()<<" "<<t.index_count()<<std::endl;
  //   }
  // }

  printMsg("Poisson Disk Sampling done", 1, tm.getElapsedTime(), threadNumber_);
}


template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::buildKDTree(int level, Scalar* _coordinates, Scalar* _normals, int _nbPoints, double ratio){


  Timer tm{};
  struct VN{
      Scalar* coordinates;
      Scalar* normals;
      int n;

      int size() const {
        return n;
      }
  };

  std::vector<int> sampling{};

  if(level==-1){
    sampling.resize(_nbPoints);
    std::iota(sampling.begin(), sampling.end(),0);
    }else{
    sampling.reserve(_nbPoints);
    for(int i = 0; i<_nbPoints; i++){
      if(is_sub_sampled(i,level)){
        sampling.push_back(i);
      }
    }
  }

  typedef typename Ponca::KdTree<Point>::PointContainer PtContainer;
  typedef typename Ponca::KdTree<Point>::IndexContainer IdxContainer;

  kdTree(level).buildWithSampling(VN({_coordinates, _normals, _nbPoints}), IdxContainer(sampling),
      [=](VN buffer, PtContainer &out){
                out.reserve(buffer.size());
                for(int i = 0; i<buffer.size(); i++){
                // std::cout<<" point "<<i<<" subsampled ? "<<is_sub_sampled(i,level)<<std::endl;
                out.push_back(Point(buffer.coordinates, buffer.normals, i, (float)ratio));
                }});
  printMsg("kd-tree built at level "+std::to_string(level) + ", point count "+std::to_string(kdTree(level).point_count()) + ", index count "+std::to_string(kdTree(level).sample_count()) , 1, tm.getElapsedTime());
}

template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::buildKDTreeND(int level, Scalar* _coordinates, int _nbPoints){


  Timer tm{};
  struct VN{
      Scalar* coordinates;
      int n;

      int size() const {
        return n;
      }
  };

  std::vector<int> sampling{};

  if(level==-1){
    sampling.resize(_nbPoints);
    std::iota(sampling.begin(), sampling.end(),0);
    }else{
    sampling.reserve(_nbPoints);
    for(int i = 0; i<_nbPoints; i++){
      if(is_sub_sampled(i,level)){
        sampling.push_back(i);
      }
    }
  }

  typedef typename Ponca::KdTree<Point>::PointContainer PtContainer;
  typedef typename Ponca::KdTree<Point>::IndexContainer IdxContainer;

  kdTree(level).buildWithSampling(VN({_coordinates, _nbPoints}), IdxContainer(sampling),
      [=](VN buffer, PtContainer &out){
                out.reserve(buffer.size());
                for(int i = 0; i<buffer.size(); i++){
                // std::cout<<" point "<<i<<" subsampled ? "<<is_sub_sampled(i,level)<<std::endl;
                out.push_back(Point(buffer.coordinates, i));
                }});
  printMsg("kd-tree built at level "+std::to_string(level) + ", point count "+std::to_string(kdTree(level).point_count()) + ", index count "+std::to_string(kdTree(level).sample_count()) , 1, tm.getElapsedTime());
}
template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::buildKDTreeBucketIndexing(int level,
                                                                  SimplexId _nbPoints,
                                                                  std::vector<SimplexId>& _bucketStartId,
                                                                  std::vector<std::vector<Scalar>>& _coordinates,
                                                                  std::vector<std::vector<SimplexId>>& _internalSampling, 
                                                                  std::vector<std::vector<SimplexId>>& _externalSampling, 
                                                                  int startBucketId){


  Timer tm{};
  struct VN{
      std::vector<std::vector<Scalar>>& coordinates;
      std::vector<std::vector<SimplexId>>& internalSampling;
      std::vector<std::vector<SimplexId>>& externalSampling;
      std::vector<SimplexId>& bucketStartId;
      int n;
      int startBucketId;

      int size() const {
        return n;
      }
  };

  std::vector<int> sampling{};

  if(level==0){
    sampling.resize(_nbPoints);
    std::iota(sampling.begin(), sampling.end(),0);
  }else{
    sampling.reserve(_nbPoints);
    for(int i = 0; i<_nbPoints; i++){
      if(is_sub_sampled(i,level)){
        sampling.push_back(i);
      }
    }
  }
  // std::cout<<"sampling size "<< sampling.size()<<std::endl;

  typedef typename Ponca::KdTree<Point>::PointContainer PtContainer;
  typedef typename Ponca::KdTree<Point>::IndexContainer IdxContainer;

  kdTree(level).buildWithSampling(VN({_coordinates, _internalSampling,_externalSampling, _bucketStartId, _nbPoints, startBucketId}), IdxContainer(sampling),
      [=](VN buffer, PtContainer &out){
                out.reserve(buffer.size());
                SimplexId firstId = buffer.bucketStartId[buffer.startBucketId];
                for(int i = 0; i<buffer.size(); i++){

                // find bucket
                int iBucket = buffer.bucketStartId.size()-1;
                while(iBucket>0 and (i+firstId)<buffer.bucketStartId[iBucket]){
                iBucket--;
                }
                SimplexId index_bucket_seq = i+firstId;
                SimplexId index_bucket = index_bucket_seq - buffer.bucketStartId[iBucket];
                SimplexId index_local = buffer.externalSampling[iBucket][buffer.internalSampling[iBucket][index_bucket]];

                out.push_back(Point(buffer.coordinates[iBucket].data(), index_local));
                // if(index_local==300697){
                //   std::cout<<" point 300697 from converter: "<<std::endl;
                //   out.back().print();
                //   std::cout<<" filling out "<<i<<" "<<std::endl;
                //   std::cout<<" kd tree index "<<i<<std::endl;
                //   std::cout<<"     bucket number "<<iBucket<<std::endl;
                //   std::cout<<"     index in bucket sequence "<<index_bucket_seq<<std::endl;
                //   std::cout<<"     index in bucket "<<index_bucket<<std::endl;
                //   std::cout<<"     local index "<<index_local<<std::endl;
                // }
                }
      });
  // std::cout<<" checking dist "<<std::endl;
  // SimplexId maxMatches=0;
  // for(auto id : kdTree(level).samples()){
  //   auto tree_query = kdTree(level).range_neighbors(id, 1e-4);
  //   SimplexId nmatches=0;
  //   for(auto m : tree_query){
  //     nmatches++;
  //   }
  //   if(nmatches>maxMatches){
  //     maxMatches=nmatches;
  //   }
  // }
  // if(maxMatches>1){
  //   std::cout<<" WARNING poisson disk sampling, max Matches "<<maxMatches<<std::endl;
  // }
  printMsg("kd-tree built at level "+std::to_string(level) + ", point count "+std::to_string(kdTree(level).point_count()) + ", index count "+std::to_string(kdTree(level).sample_count()) , 1, tm.getElapsedTime());
  // kdTree(level).print(std::cout, true);
}

template <class Scalar, typename Point>
Scalar ScaleSpaceSampling<Scalar, Point>::estimateDensity(Scalar* coordinates, Scalar* normals, SimplexId nbPoints, int percent_tries, SimplexId nbNeighbors){

  // build kd tree
  Ponca::KdTree<pointCloud::MyPoint> tree;

  struct VN{
      Scalar* coordinates;
      int n;

      int size() const {
        return n;
      }
  };

  tree.build(VN({coordinates, nbPoints}), [&](VN buffer, Ponca::KdTree<MyPoint>::PointContainer &out){
                out.reserve(buffer.size());
                for(int i = 0; i<buffer.size(); i++){
                  out.push_back(MyPoint(buffer.coordinates, i));
                }});

  // generate random point ids
  std::mt19937 random_engine{};
  random_engine.seed(0);
  std::uniform_int_distribution<int> dis(0, nbPoints-1);

  SimplexId nbTries = nbPoints*percent_tries/100+1;
  // for each random id, compute the average distance to its 10 nearest neighbors
  Scalar total_mean_dist=0;
  Scalar total_mean_normal_dist=0;

  for(SimplexId i=0; i<nbTries; i++){
    SimplexId id = dis(random_engine);

    auto tree_query = tree.k_nearest_neighbors(id, nbNeighbors);
    int neighborhood_size=0;
    Scalar mean_dist=0;
    Scalar mean_normal_dist=0;
    for(auto neighbor : tree_query){
     Scalar dist = (coordinates[3*id+0]-coordinates[3*neighbor+0])*(coordinates[3*id+0]-coordinates[3*neighbor+0])
            + (coordinates[3*id+1]-coordinates[3*neighbor+1])*(coordinates[3*id+1]-coordinates[3*neighbor+1])
            + (coordinates[3*id+2]-coordinates[3*neighbor+2])*(coordinates[3*id+2]-coordinates[3*neighbor+2]);
     Scalar normal_dist = (normals[3*id+0]-normals[3*neighbor+0])*(normals[3*id+0]-normals[3*neighbor+0])
            + (normals[3*id+1]-normals[3*neighbor+1])*(normals[3*id+1]-normals[3*neighbor+1])
            + (normals[3*id+2]-normals[3*neighbor+2])*(normals[3*id+2]-normals[3*neighbor+2]);
     dist = std::sqrt(dist);
     normal_dist = std::sqrt(normal_dist);

     neighborhood_size++;
     mean_dist+=dist;
     mean_normal_dist+=normal_dist;
    }
    total_mean_dist += mean_dist/neighborhood_size/nbTries;
    total_mean_normal_dist += mean_normal_dist/neighborhood_size/nbTries;
  }
  return total_mean_dist/total_mean_normal_dist;
}

template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::buildKnnGraph(){
  m_KnnGraph = new Ponca::KnnGraph<Point>(kdTree(0),10);
}

template <class Scalar, typename Point>
int ScaleSpaceSampling<Scalar, Point>::findMinScale(SimplexId maxPoints) {
  SimplexId nbPoints = -1;
  for(int iScale=-1; iScale<scale_count(); iScale++){
    nbPoints = kdTree(iScale).sample_count();
    if(nbPoints<=maxPoints){
      return iScale;
    }
  }
  printWrn("None of the scales have less than maxPoints, I chose the biggest one, with "+std::to_string(nbPoints)+" points");
  return scale_count()-1;
 }

template <class Scalar, typename Points>
void ScaleSpaceSampling<Scalar, Points>::computeDensityField(int iScale, SimplexId k, Scalar* coords, int dim, SimplexId nbPoints, std::vector<SimplexId>& sampling_indexes, Scalar* output_density){

  printMsg("computing density");
  auto& tree = kdTree(iScale);

  for(SimplexId iPoint : tree.samples()){
    auto tquery = tree.k_nearest_neighbors(iPoint, k);

    Scalar sum_dist=0;
    SimplexId nbNeighbors=0;
    // std::cout<<"point "<<iPoint;
    for(auto iNeighbor : tquery) {
      Scalar dist = 0;
      for(int i=0; i<dim; i++){
        Scalar tmp =(coords[dim*iPoint+i]-coords[dim*iNeighbor+i]);
        dist += tmp*tmp;
      }
      // dist = sqrt(dist);
      sum_dist+=dist;
      nbNeighbors++;
    }
    output_density[iPoint] = nbNeighbors/sum_dist;
  }
  printMsg("computing density done");
}

// template <class Scalar, typename Points>
// void ScaleSpaceSampling<Scalar, Points>::computeDensityFieldFromCoords(int iScale, SimplexId k, float* coords, int dim, SimplexId nbPoints, float* output_density){

//   printMsg("computing density");

//   Ponca::KdTree<pointCloud::MyPoint5D> tree;

//   Timer tm{};
//   struct VN{
//       Scalar* coordinates;
//       int n;

//       int size() const {
//         return n;
//       }
//   };

//   std::vector<int> sampling{};

//   if(level==0){
//     sampling.resize(_nbPoints);
//     std::iota(sampling.begin(), sampling.end(),0);
//     }else{
//     sampling.reserve(_nbPoints);
//     for(int i = 0; i<_nbPoints; i++){
//       if(is_sub_sampled(i,level)){
//         sampling.push_back(i);
//       }
//     }
//   }

//   typedef typename Ponca::KdTree<pointCloud::MyPoint5D>::PointContainer PtContainer;
//   typedef typename Ponca::KdTree<pointCloud::MyPoint5D>::IndexContainer IdxContainer;

//   kdTree(level).buildWithSampling(VN({coordinates, nbPoints}), IdxContainer(sampling),
//       [=](VN buffer, PtContainer &out){
//                 out.reserve(buffer.size());
//                 for(int i = 0; i<buffer.size(); i++){
//                 // std::cout<<" point "<<i<<" subsampled ? "<<is_sub_sampled(i,level)<<std::endl;
//                 out.push_back(Point(buffer.coordinates, i));
//                 }});



//   for(SimplexId iPoint : tree.index_data()){
//     auto tquery = tree.k_nearest_neighbors(iPoint, k);

//     float sum_dist=0;
//     SimplexId nbNeighbors=0;
//     // std::cout<<"point "<<iPoint;
//     for(auto iNeighbor : tquery) {
//       float dist = 0;
//       for(int i=0; i<dim; i++){
//         float tmp =(coords[dim*iPoint+i]-coords[dim*iNeighbor+i]);
//         dist += tmp*tmp;
//       }
//       // dist = sqrt(dist);
//       sum_dist+=dist;
//       nbNeighbors++;
//     }
//     output_density[iPoint] = nbNeighbors/sum_dist;
//   }
//   printMsg("computing density done");
// }

template <class Scalar, typename Points>
void ScaleSpaceSampling<Scalar, Points>::computeDensityField(int iScale, SimplexId k, std::vector<std::vector<double>>& distMat, Scalar radius, Scalar* output_density){

  printMsg("computing density with dist mat");
  
  auto& tree = kdTree(iScale);
  SimplexId nbPoints = tree.sample_count();
  if(nbPoints != distMat.size()){
    printWrn("Sizes don't match btwn kdtree and distmat");
    std::cout<<nbPoints<<" vs "<<distMat.size()<<std::endl;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId i=0; i<nbPoints; i++){
    SimplexId nbNeighbors=0;

    std::vector<double> dist;
    std::copy(distMat[i].begin(), distMat[i].end(), std::back_inserter(dist));
    std::sort(dist.begin(), dist.end());

    double sum_dist = 0;
    for(int j = 1; j<=k; j++){
      sum_dist+=dist[j];
    }


    output_density[tree.samples()[i]] = (Scalar)k/sum_dist;

  }

  // for(SimplexId i=0; i<nbPoints; i++){
  //   for(SimplexId j=i+1; j<nbPoints; j++){
  //     double d_min=std::min(output_density[tree.index_data()[i]], output_density[tree.index_data()[j]]);
  //     d_min = std::tanh(1/d_min);
  //     if(std::isnan(d_min)){
  //       std::cout<<"dmin is nan "<<output_density[tree.index_data()[i]]<<" "<<output_density[tree.index_data()[j]]<<std::endl;
  //     }
  //     std::cout<<"dmin "<<d_min<<" "<<output_density[tree.index_data()[i]]<<" "<<output_density[tree.index_data()[j]]<<std::endl;
  //     distMat[i][j] *= d_min;
  //     distMat[j][i] *= d_min;
  //   }
  // }

  printMsg("computing density done");
}

template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::filterScaleByDensity(int level, Scalar* _coordinates, int _nbPoints, Scalar* density, Scalar min_density, std::vector<SimplexId>& densitySampling, int* isStable){
  Timer tm{};
  struct VN{
      Scalar* coordinates;
      int n;

      int size() const {
        return n;
      }
  };

  const std::vector<SimplexId>& currentSampling = kdTree(level).samples();
  std::vector<int> sampling{};
  SimplexId nbPointSampling = currentSampling.size();
  densitySampling.reserve(nbPointSampling);

  // std::fill(kdTree(level).index_data().begin(), kdTree(level).index_data().end(), 0);

  sampling.reserve(nbPointSampling);
  for(SimplexId I = 0; I<nbPointSampling; I++){
    SimplexId i = currentSampling[I];
    if(is_sub_sampled(i,level) and density[i]>min_density and isStable[i]){
      sampling.push_back(i);
      densitySampling.push_back(I);
    }
  }
  // currentSampling.clear();

  typedef typename Ponca::KdTree<Point>::PointContainer PtContainer;
  typedef typename Ponca::KdTree<Point>::IndexContainer IdxContainer;

  kdTree(level).buildWithSampling(VN({_coordinates, _nbPoints}), IdxContainer(sampling),
      [=](VN buffer, PtContainer &out){
                out.reserve(buffer.size());
                for(int i = 0; i<buffer.size(); i++){
                // std::cout<<" point "<<i<<" subsampled ? "<<is_sub_sampled(i,level)<<std::endl;
                out.push_back(Point(buffer.coordinates, i));
                }});
  printMsg("kd-tree built at level "+std::to_string(level) + ", point count "+std::to_string(kdTree(level).point_count()) + ", index count "+std::to_string(kdTree(level).sample_count()) , 1, tm.getElapsedTime());
}

template <class Scalar, typename Point>
void ScaleSpaceSampling<Scalar, Point>::filterUnstablePoints(int level, Scalar* _coordinates, int _nbPoints, std::vector<SimplexId>& densitySampling, int* isStable){
  // Timer tm{};
  // struct VN{
  //     Scalar* coordinates;
  //     int n;

  //     int size() const {
  //       return n;
  //     }
  // };

  const std::vector<SimplexId>& currentSampling = kdTree(level).samples();
  std::vector<int> sampling{};
  SimplexId nbPointSampling = currentSampling.size();
  densitySampling.reserve(nbPointSampling);


  std::cout<<" Filtering unstable points"<<std::endl;
  sampling.reserve(nbPointSampling);
  for(SimplexId I = 0; I<nbPointSampling; I++){
    SimplexId i = currentSampling[I];
    if(is_sub_sampled(i,level) and isStable[I]){
      // sampling.push_back(i);
      // std::cout<<" "<<I;
      densitySampling.push_back(I);
    }
  }
  // std::cout<<std::endl;
  // currentSampling.clear();

  // typedef typename Ponca::KdTree<Point>::PointContainer PtContainer;
  // typedef typename Ponca::KdTree<Point>::IndexContainer IdxContainer;

  // kdTree(level).buildWithSampling(VN({_coordinates, _nbPoints}), IdxContainer(sampling),
  //     [=](VN buffer, PtContainer &out){
  //               out.reserve(buffer.size());
  //               for(int i = 0; i<buffer.size(); i++){
  //               // std::cout<<" point "<<i<<" subsampled ? "<<is_sub_sampled(i,level)<<std::endl;
  //               out.push_back(Point(buffer.coordinates, i));
  //               }});
  // printMsg("kd-tree built at level "+std::to_string(level) + ", point count "+std::to_string(kdTree(level).point_count()) + ", index count "+std::to_string(kdTree(level).index_count()) , 1, tm.getElapsedTime());
  std::cout<<"filtering done "<<std::endl;
}
