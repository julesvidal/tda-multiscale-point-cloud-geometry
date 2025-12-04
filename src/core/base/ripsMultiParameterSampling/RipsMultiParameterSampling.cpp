#include "DataTypes.h"
#include "Debug.h"
#include "Ponca/src/Fitting/enums.h"
#include "RipsTree.h"
#include <RipsMultiParameterSampling.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <omp.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <iostream>
#include <fstream>

ttk::RipsMultiParameterSampling::RipsMultiParameterSampling() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("RipsMultiParameterSampling");
}

double RipsMultiParameterSampling::computeJaccardIndex(const std::vector<SimplexId>& cc1, const std::vector<SimplexId>& cc2) const{
    std::vector<SimplexId> interVec;
    std::vector<SimplexId> unionVec;
    std::set_intersection(cc1.cbegin(), cc1.cend(), cc2.cbegin(), cc2.cend(), std::back_inserter(interVec));
    std::set_union(cc1.cbegin(), cc1.cend(), cc2.cbegin(), cc2.cend(), std::back_inserter(unionVec));

    return (double)interVec.size()/(double)unionVec.size();
}

void RipsMultiParameterSampling::findMatchingComponent(SimplexId ccId,
                          SimplexId inputTreeId, double inputThreshold,
                          SimplexId targetTreeId, SimplexId targetThresholdId, double targetThreshold,
                          SimplexId nCol,
                          SimplexId& matchingCC,
                          SimplexId& localCCId,
                          double& jaccard){
    // std::cout<<"\nFIND"<<std::endl;
    if(inputTreeId>=trees_.size() or targetTreeId>=trees_.size()){
    matchingCC=-1;
    jaccard=-1;
      return;
    }
    const auto& inputTree = trees_[inputTreeId];
    const auto& targetTree = trees_[targetTreeId];

    std::vector<SimplexId> inputPointIds{};
    inputTree.getComponentPointIdsAtThreshold(ccId, inputThreshold, inputPointIds);
    std::sort(inputPointIds.begin(), inputPointIds.end());
    // std::vector<SimplexId> targetCCIds{};
    // targetTree.findContainingComponents(inputPointIds, targetThreshold, targetCCIds);

    // for(SimplexId targetCC : targetCCIds){
    // std::cout<<"    find matching cc of "<<ccId<<" in "<<targetTreeId<<","<<targetThreshold<<std::endl;
    // localCCId = 0;
      if(targetThreshold>targetTree.maxDeath()){
        targetThreshold = targetTree.maxDeath();
      }

    matchingCC=-1;
    jaccard=-1;
    localCCId=-1;

// #ifdef TTK_ENABLE_OPENMP
// #pragma omp parallel for num_threads(threadNumber_)
// #endif // TTK_ENABLE_OPENMP
    for(SimplexId llCCId=0; llCCId<NbCCMax; llCCId++){

      SimplexId fullId = getFullGridId(llCCId, targetTreeId, targetThresholdId, nCol, NbCCMax);
      SimplexId targetCC = globalToLocal_[fullId];

      if(targetCC>-1){
        std::vector<SimplexId> targetPointIds{};
        targetTree.getComponentPointIdsAtThreshold(targetCC, targetThreshold, targetPointIds);
        std::sort(targetPointIds.begin(), targetPointIds.end());

        double j = computeJaccardIndex(inputPointIds, targetPointIds);

        if(j>JaccardConstraint){
          // std::cout<<"thread "<<omp_get_thread_num()<<" "<<j<<" "<<JaccardConstraint<<" in "<<ccId<<"-"<<targetCC<<" lid:"<<llCCId<<" fid:"<<fullId<<" ("<<targetTreeId<<","<<targetThresholdId<<")"<<std::endl;
          jaccard=j;
          matchingCC=targetCC;
          localCCId = llCCId;
          return;
        }

      }
    }

    // std::cout<<"done"<<std::endl;
    return;
}


SimplexId RipsMultiParameterSampling::getFullGridId(SimplexId ccId, SimplexId iRow, SimplexId iCol, SimplexId nCol, SimplexId n)const{
  SimplexId result = ccId + n*(nCol-iCol-1) + n*nCol*iRow;
  return result;
}

SimplexId RipsMultiParameterSampling::localToGlobal(SimplexId ccId, SimplexId iRow, SimplexId iCol, SimplexId nCol){
  for(int i=0; i<NbCCMax; i++){
    SimplexId fullId = getFullGridId(i, iRow, iCol, nCol, NbCCMax);
    if(globalToLocal_[fullId]==ccId){
      return fullId;
    }
  }
  return -1;
}

void RipsMultiParameterSampling::getGridTuple(SimplexId fullGridId, SimplexId nCol, SimplexId n, SimplexId& ccId, SimplexId& iRatio, SimplexId& iThreshold)const{
  iRatio = (fullGridId/n)/nCol;
  iThreshold = nCol-1-((fullGridId/n)%nCol);
  ccId = fullGridId%n;
}


void RipsMultiParameterSampling::compute2DParentsOneToOne(std::vector<std::vector<double>> persistenceGrid, const std::vector<double>& ratioVec){

  printMsg("Compute 2d pass", 0, 0, threadNumber_);

  parents2D_.clear();
  Timer tm{};
  SimplexId nbPoints = trees_[0].getNumberOfPoints();

  SimplexId nRow = persistenceGrid.size();
  SimplexId nCol = persistenceGrid[0].size();

  children2D_.clear();
  children2D_.resize(NbCCMax*nCol*nRow);
  parents2D_.clear();
  parents2D_.resize(NbCCMax*nCol*nRow);
  std::fill(parents2D_.begin(), parents2D_.end(), -1);

  // fathers_.resize(NbCCMax*ratioVec.size()*persistenceGrid[0].size());
  // sons_.resize(NbCCMax*ratioVec.size()*persistenceGrid[0].size());

  // madatory to keep track of inclusions
  buildGlobalToLocal(persistenceGrid, ratioVec);

  // persistence2d_.resize(NbCCMax*ratioVec.size()*persistenceGrid[0].size());
  // std::fill(persistence2d_.begin(), persistence2d_.end(), 0.0);

  for(SimplexId iCol=nCol-1; iCol>=0; iCol--){
    for(SimplexId iRow=0; iRow<nRow; iRow++){
    auto& tree = trees_[iRow];


      double threshold = persistenceGrid[iRow][iCol];
      // SimplexId globalCCId=0;


      // std::cout<<" Initiate propagation from "<<iRow<<" "<<iCol<<", threshold "<<threshold<<std::endl;
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
      for(SimplexId localCCId=0; localCCId<NbCCMax; localCCId++){

        SimplexId fullId = getFullGridId(localCCId, iRow, iCol, nCol,NbCCMax);
        SimplexId ccId = globalToLocal_[fullId];

        if(ccId>-1){

        // until we reach the threshold
        // if(tree.sizeAtThreshold(ccId,threshold)<SizeMin){
        //   continue;
        // }
        // if(tree.getPersistence(ccId)<threshold){
        //   break;
        // }

        // std::cout<<"trying to match "<<fullId<<" "<<globalCCId<<" "<<ccId<<","<<iRow<<","<<iCol<<std::endl;
        //   std::cout<<"   -  size  "<<tree.sizeAtThreshold(ccId, threshold)<<"/"<<SizeMin<<std::endl;
        //   std::cout<<"   -  death  "<<tree.getPersistence(ccId)<<"/"<<SizeMin<<std::endl;

        // globalCCId++;

        // const auto& it = parents2D_.find(fullId);
        // SimplexId parentId = fullId;
        // if(it!=parents2D_.end()){ // ccId has already a parent of its own
        //   continue;
        //   // parentId = it->second.first;
        // }
        // std::cout<<" from "<<fullId<<std::endl;
          if(parents2D_[fullId]==-1){
            propagate2dParent(fullId, ccId, iRow, iCol, iRow, iCol, persistenceGrid, ratioVec);
          }
        }
      }
      // printMsg("Compute 2d pass", (double)((nCol-1-iCol)*nRow+iRow)/(nRow*nCol), tm.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);
      // if(globalCCId>=NbCCMax){
      //   break;
      // }
      printMsg("Compute 2d pass", (double)((nCol-1-iCol)*nRow+iCol)/(nRow*nCol), tm.getElapsedTime(), threadNumber_, debug::LineMode::REPLACE);
    }
  }
  printMsg("Compute 2d pass", 1, tm.getElapsedTime(), threadNumber_);
}

void RipsMultiParameterSampling::compute2DParents(std::vector<std::vector<double>> persistenceGrid){

  //parents2D_.clear();
  //printMsg("Compute 2d pass", 0, 0, threadNumber_);
  //Timer tm{};
  //SimplexId nbPoints = trees_[0].getNumberOfPoints();
  //SimplexId nRow = persistenceGrid.size();

  //// for(SimplexId iRow=0; iRow<nRow; iRow++){
  //for(SimplexId iRow=nRow-1; iRow>=0; iRow--){
  //  SimplexId nCol = persistenceGrid[iRow].size();
  //  auto& tree = trees_[iRow];
  //  for(SimplexId iCol=0; iCol<nCol; iCol++){
  //    double threshold = persistenceGrid[iRow][iCol];
  //    SimplexId globalCCId=0;
  //    for(SimplexId ccId=0; ccId<tree.getNumberOfPoints(); ccId++){
  //      // until we reach the threshold
  //      if(tree.sizeAtThreshold(ccId,threshold)<SizeMin /*or parents2D_.find(getFullGridId(ccId, iRow, iCol, nCol, nbPoints))!=parents2D_.end()*/){
  //        // std::cout<<"   -  size matter "<<tree.getSize(ccId)<<std::endl;
  //        continue;
  //      }
  //      if(tree.getPersistence(ccId)<threshold){
  //        // std::cout<<"   -  thresh matter"<<std::endl;
  //        break;
  //      }
  //      SimplexId fullId = getFullGridId(globalCCId, iRow, iCol, nCol,NbCCMax);
  //      std::cout<<"trying to match "<<fullId<<" "<<globalCCId<<" "<<ccId<<","<<iRow<<","<<iCol<<std::endl;
  //        std::cout<<"   -  size  "<<tree.sizeAtThreshold(ccId, threshold)<<"/"<<SizeMin<<std::endl;
  //        std::cout<<"   -  death  "<<tree.getPersistence(ccId)<<"/"<<SizeMin<<std::endl;
  //      globalCCId++;
  //      const auto& it = parents2D_.find(fullId);
  //      SimplexId parentId = fullId;
  //      if(it!=parents2D_.end()){ // ccId has already a parent of its own
  //        // std::cout<<"        "<<fullId<<" already has parent "<<parents2D_[fullId].first<<std::endl;
  //        continue;
  //        parentId = it->second.first;
  //      }
  //      //
  //      for(SimplexId iRow2=0; iRow2<nRow; iRow2++){
  //        SimplexId nCol2 = persistenceGrid[iRow2].size();
  //        // auto& tree2 = trees_[iRow2];
  //        for(SimplexId iCol2=0; iCol2<nCol2; iCol2++){
  //          double threshold2 = persistenceGrid[iRow2][iCol2];
  //          SimplexId matchingCC=-1;
  //          double jaccard=-1;
  //          SimplexId localCCId=-1;
  //          std::cout<<"          find matching cc in "<<iRow2<<" "<<iCol2<<std::endl;
  //          findMatchingComponent(ccId, iRow, threshold, iRow2, threshold2, matchingCC, localCCId, jaccard);
  //          if(localCCId>-1 ){
  //            SimplexId fullId2 = getFullGridId(localCCId, iRow2, iCol2, nCol2, NbCCMax);
  //            std::cout<<"          matching "<<fullId2<<" to "<<fullId<<std::endl;
  //              std::cout<<"          "<<ccId<<" "<<iRow<<" "<<iCol<<std::endl;
  //              std::cout<<"          "<<matchingCC<<" "<<iRow2<<" "<<iCol2<<std::endl;
  //            if(parents2D_.find(fullId2)!=parents2D_.end()){
  //              // if matching cc already had a parent
  //              // then compare jaccards
  //              // std::cout<<"        -> "<<fullId2<<" already has parent "<<paren
  //              // if(jaccard>parents2D_[fullId2].second){
  //              //   parents2D_[fullId2] = std::make_pair(parentId, jaccard);
  //              // }
  //            }else{
  //              parents2D_.emplace(fullId2, std::make_pair(parentId, jaccard));
  //            }
  //            // tree2.setParent2D(matchingCC, getFullGridId(ccId, iRow, iCol, nCol, nbPoints));
  //          }
  //        }
  //      }
  //      // std::cout<<"         try globalCCId "<<globalCCId<<std::endl;
  //    }
  //    if(globalCCId>=NbCCMax){
  //      break;
  //    }
  //    printMsg("Compute 2d pass", (double)((nRow-1-iRow)*nCol+iCol)/(nRow*nCol), tm.getElapsedTime(), threadNumber_);
  //  }
  //}
  //printMsg("Compute 2d pass", 1, tm.getElapsedTime(), threadNumber_);
}

/*
 * check if component fullId0 is included in component fullId1
 */

bool RipsMultiParameterSampling::isPointSetIncluded(const std::vector<SimplexId>& cc1, const std::vector<SimplexId>& cc2) const{
    std::vector<SimplexId> interVec;
    std::vector<SimplexId> unionVec;
    std::set_intersection(cc1.cbegin(), cc1.cend(), cc2.cbegin(), cc2.cend(), std::back_inserter(interVec));
    std::set_union(cc1.cbegin(), cc1.cend(), cc2.cbegin(), cc2.cend(), std::back_inserter(unionVec));

    return (interVec.size()==cc1.size() or unionVec.size()==cc2.size());
}

bool RipsMultiParameterSampling::isIncluded(SimplexId ccId0, SimplexId iRatio0, SimplexId iThreshold0, SimplexId ccId1, SimplexId iRatio1, SimplexId iThreshold1) const{
  // if(iRatio0 != iRatio1){
  //   return false;
  // }
  // if(ccId0<ccId1){
  //   return false;
  // }
  // if(iThreshold0 > iThreshold1){
  //   return false;
  // }
  // if(ccId0 == ccId1){
  //   return true;
  // }

  // if(trees_[iRatio0].getComponent(ccId0).death() >= persistenceGrid_[iRatio0][iThreshold1]){
  //   return false;
  // }
  // auto parent = trees_[iRatio0].getComponent(ccId0).parent();
  // SimplexId ccId = ccId0;
  // while(parent!=ccId){
  //   if(parent == ccId1){
  //     return true;
  //   }
  //   ccId = parent;
  //   parent = trees_[iRatio0].getComponent(parent).parent(); 
  // }
  // return false;
}

void RipsMultiParameterSampling::outputGraph(std::string outputFileName, const std::vector<std::vector<double>>& persistenceGrid){

  printMsg("Computing output graph matrix");

  SimplexId nRow = persistenceGrid.size();
  SimplexId nCol = persistenceGrid[0].size();
  SimplexId totalNbIds = nRow*nCol*NbCCMax;

  std::vector<std::vector<double>> distanceMatrix(nbComponents_, std::vector<double>(nbComponents_, 0));

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(SimplexId fullId0=0; fullId0<totalNbIds; fullId0++){
    SimplexId ccId0 = globalToLocal_[fullId0];
    if(ccId0==-1) {
      continue;
    }
    SimplexId localId0, iRatio0, iThreshold0;
    getGridTuple(fullId0, nCol, NbCCMax, localId0, iRatio0, iThreshold0);

    SimplexId fullCompactId0 = globalToFullCompact_[fullId0];

    // get components point ids
    std::vector<SimplexId> inputPointIds0{};
    trees_[iRatio0].getComponentPointIdsAtThreshold(ccId0, persistenceGrid[iRatio0][iThreshold0], inputPointIds0);
    std::sort(inputPointIds0.begin(), inputPointIds0.end());

    for(SimplexId fullId1=fullId0+1; fullId1<totalNbIds; fullId1++){
    SimplexId ccId1 = globalToLocal_[fullId1];
      if(ccId1==-1) {
        continue;
      }
      SimplexId localId1, iRatio1, iThreshold1;
      getGridTuple(fullId1, nCol, NbCCMax, localId1, iRatio1, iThreshold1);

      SimplexId fullCompactId1 = globalToFullCompact_[fullId1];

      // get components point ids
      std::vector<SimplexId> inputPointIds1{};
      trees_[iRatio1].getComponentPointIdsAtThreshold(ccId1, persistenceGrid[iRatio1][iThreshold1], inputPointIds1);
      std::sort(inputPointIds1.begin(), inputPointIds1.end());

      //compute jaccard
      double jaccard = computeJaccardIndex(inputPointIds0, inputPointIds1);
      distanceMatrix[fullCompactId0][fullCompactId1] = jaccard;
      distanceMatrix[fullCompactId1][fullCompactId0] = jaccard;

      // check if there is inclusion 
      if(jaccard<JaccardConstraint){
        if(isPointSetIncluded(inputPointIds0, inputPointIds1) or isPointSetIncluded(inputPointIds1, inputPointIds0)){
          distanceMatrix[fullCompactId0][fullCompactId1] = 1;
          distanceMatrix[fullCompactId1][fullCompactId0] = 1;
        }
      }
      // if(isIncluded(ccId0, iRatio0, iThreshold0, ccId1, iRatio1, iThreshold1) or isIncluded(ccId1, iRatio1, iThreshold1, ccId0, iRatio0, iThreshold0)){
      //   distanceMatrix[fullCompactId0][fullCompactId1] = 1;
      //   distanceMatrix[fullCompactId1][fullCompactId0] = 1;
      // }
    }
  }
  // output matrix to file
  printMsg("Writing matrix to disk");
  std::ofstream outputFile;
  outputFile.open(outputFileName);
  std::stringstream header;
  header << "fullId";
  for(auto i=0; i<nbComponents_; i++){
    header<<", component"<<std::setw(3)<<std::setfill('0')<<i;
  }
  outputFile << header.str() << std::endl;

  for(auto i=0; i<nbComponents_; i++){
    std::stringstream row;
    row << fullCompactToGlobal_[i];
    for(auto j=0; j<nbComponents_; j++){
      row << ", " << distanceMatrix[i][j];
    }
    outputFile << row.str() << std::endl;
  }
  outputFile.close();
  printMsg("all done");
}


void RipsMultiParameterSampling::findPersistentComponents(const std::vector<std::vector<double>>& persistenceGrid, const std::vector<double>& ratioVec, std::vector<double>& persistence2d){

  Timer tm{};

  SimplexId nRow = persistenceGrid.size();
  SimplexId nCol = persistenceGrid[0].size();
  parentList_.clear();
  parentOrder_.clear();

  // std::vector<std::vector<SimplexId>> result(NbCCMax*nCol*nRow);
  persistence2d.resize(NbCCMax*nCol*nRow);
  std::fill(persistence2d.begin(), persistence2d.end(), 0);

  // for(const auto& it : parents2D_){
  //     result[it.second.first].push_back(it.first);
  // }
  for(int fullIdParent=0; fullIdParent<children2D_.size(); fullIdParent++){
    if(!children2D_[fullIdParent].empty()){
      parentOrder_.push_back(fullIdParent);
      // std::cout<<"persistence of "<<fullIdParent<<std::endl;
      double minRatio = std::numeric_limits<double>::max();
      double maxRatio = std::numeric_limits<double>::min();
      double minThresh = std::numeric_limits<double>::max();
      double maxThresh = std::numeric_limits<double>::min();
      for(auto fullIdChild : children2D_[fullIdParent]){
        SimplexId iRatio, iThreshold, localIdChild;
        getGridTuple(fullIdChild, nCol, NbCCMax, localIdChild, iRatio, iThreshold);
        double ratio = ratioVec[iRatio];
        double threshold = persistenceGrid[iRatio][iThreshold];
        // if(iRatio>0){

        // }
        // std::cout<<"       - child "<<fullIdChild<<" r "<<ratio<<" t "<<threshold<<std::endl;
        if(ratio>maxRatio){
          maxRatio=ratio;
        }
        if(ratio < minRatio){
          minRatio=ratio;
        }
        if(threshold>maxThresh){
          maxThresh=threshold;
        }
        if(threshold<minThresh){
          minThresh=threshold;
        }
        if(iThreshold>0){
          persistence2d[fullIdParent] += threshold-persistenceGrid[iRatio][iThreshold-1];
        }

        {
        }
        // else{
        // }
      }
      // persistence2d[fullIdParent] = (1+(maxRatio-minRatio)*(maxThresh-minThresh))*(double)children2D_[fullIdParent].size();
      // persistence2d[fullIdParent] = (maxRatio-minRatio)*(maxThresh-minThresh); //(double)children2D_[fullIdParent].size();
      // persistence2d[fullIdParent] = (double)children2D_[fullIdParent].size();
      // std::cout<<"     final persistence "<<persistence2d[fullIdParent]<<std::endl;
      // std::cout<<"      "<<minRatio<<","<<maxRatio<<"  "<<minThresh<<","<<maxThresh<<std::endl;


      { // for inclusions
        /*
        for(int ifather=0; ifather<fathers_[fullIdParent].size(); ifather++){
          SimplexId f = fathers_[fullIdParent][ifather];
          f = parents2D_[f].first;
          fathers_[fullIdParent][ifather] = f;
        }
        std::sort(fathers_[fullIdParent].begin(), fathers_[fullIdParent].end());
        auto last = std::unique(fathers_[fullIdParent].begin(), fathers_[fullIdParent].end());
        fathers_[fullIdParent].erase(last, fathers_[fullIdParent].end());

        for(auto f : fathers_[fullIdParent]){
          sons_[f].push_back(fullIdParent);
        }
        */
      }

    }
  }

  const auto cmp = [&](const SimplexId &i, const SimplexId &j){
    return persistence2d[i]>persistence2d[j];
  };

  std::sort(parentOrder_.begin(), parentOrder_.end(), cmp);
  parentList_.resize(NbCCMax*nCol*nRow);
  std::fill(parentList_.begin(), parentList_.end(), -1);
  for(int i=0; i<parentOrder_.size(); i++){
    parentList_[parentOrder_[i]]=i;
  }

  printMsg("computed 2d persistence", 1, tm.getElapsedTime(), threadNumber_);

}

void RipsMultiParameterSampling::buildGlobalToLocal(std::vector<std::vector<double>> persistenceGrid, const std::vector<double>& ratioVec){
  std::cout<<"build global to local"<<std::endl;
  globalToLocal_.resize(NbCCMax*ratioVec.size()*persistenceGrid[0].size());
  globalToFullCompact_.resize(NbCCMax*ratioVec.size()*persistenceGrid[0].size());
  fullCompactToGlobal_.resize(0);
  std::fill(globalToLocal_.begin(), globalToLocal_.end(), -1);
  std::fill(globalToFullCompact_.begin(), globalToFullCompact_.end(), -1);
  SimplexId nRow = persistenceGrid.size();
  SimplexId nCol = persistenceGrid[0].size();
  SimplexId nbPoints = trees_[0].getNumberOfPoints();

  SimplexId fullCompactCount=0;
  for(SimplexId iRow=0; iRow<nRow; iRow++){
    for(SimplexId iCol=nCol-1; iCol>=0; iCol--){
    auto& tree = trees_[iRow];

      double threshold = persistenceGrid[iRow][iCol];
      SimplexId globalCCId=0;

      for(SimplexId ccId=0; ccId<nbPoints; ccId++){
        SimplexId fullId = getFullGridId(globalCCId, iRow, iCol, nCol,NbCCMax);

        if(tree.sizeAtThreshold(ccId,threshold)<SizeMin){
          continue;
        }
        if(tree.getPersistence(ccId)<threshold){
          break;
        }
        globalToLocal_[fullId] = ccId;
        globalToFullCompact_[fullId] = fullCompactCount;
        fullCompactToGlobal_.push_back(fullId);
        fullCompactCount++;
        globalCCId++;
      }
    }
  }
  nbComponents_=fullCompactCount;
  std::cout<<"build global to local done"<<std::endl;
  printMsg("Total number of components: "+std::to_string(nbComponents_));
}

void RipsMultiParameterSampling::propagate2dParent(SimplexId parentId, SimplexId ccId, SimplexId originRow, SimplexId originCol, SimplexId targetRow, SimplexId targetCol, const std::vector<std::vector<double>>& persistenceGrid, const std::vector<double>& ratioVec){

    if(targetRow<0 or targetRow>=persistenceGrid.size() or targetCol<0 or targetCol>=persistenceGrid[originRow].size()){
      return;
    }

    double originThreshold = persistenceGrid[originRow][originCol];
    double targetThreshold = persistenceGrid[targetRow][targetCol];
    SimplexId nCol = persistenceGrid[targetRow].size();

    SimplexId matchingCC=-1;
    double jaccard=-1;
    SimplexId localCCId=-1;

    // std::cout<<"          find cc matching to "<<ccId<<","<<originRow<<","<<originCol<<" in "<<targetRow<<" "<<targetCol<<std::endl;
    findMatchingComponent(ccId, originRow, originThreshold, targetRow, targetCol, targetThreshold, nCol, matchingCC, localCCId, jaccard);

    if(localCCId>-1 ){
      SimplexId fullId = getFullGridId(localCCId, targetRow, targetCol, nCol, NbCCMax);

        // std::cout<<"          matching "<<fullId<<" to "<<parentId<<std::endl;
        // std::cout<<"          "<<matchingCC<<" "<<targetRow<<" "<<targetCol<<std::endl;
        // std::cout<<"          "<<ccId<<" "<<originRow<<" "<<originCol<<std::endl;
      if(parents2D_[fullId]==-1){
        // if matching cc doesn't already have a parent
        // if(originRow==targetRow){
        //   persistence2d_[parentId]+=std::abs(persistenceGrid[originRow][originCol]-persistenceGrid[targetRow][targetCol]);
        // }else if(originCol==targetCol){
        //   persistence2d_[parentId]+=std::abs(ratioVec[originRow]-ratioVec[targetRow]);

        // }
        parents2D_[fullId] = parentId;
        children2D_[parentId].push_back(fullId);


        if(originRow==targetRow){
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow-1, targetCol, persistenceGrid, ratioVec); //propagate below
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow+1, targetCol, persistenceGrid, ratioVec); //propagate above
        }else if(originRow==targetRow+1){
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow-1, targetCol, persistenceGrid, ratioVec); //propagate below
        }else if(originRow==targetRow-1){
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow+1, targetCol, persistenceGrid, ratioVec); //propagate above
        }

        if(originCol==targetCol){
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow, targetCol-1, persistenceGrid, ratioVec); //propagate on left
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow, targetCol+1, persistenceGrid, ratioVec); //propagate on right
        }else if(originCol==targetCol+1){
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow, targetCol-1, persistenceGrid, ratioVec); //propagate on left
        }else if(originCol==targetCol-1){
          propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow, targetCol+1, persistenceGrid, ratioVec); //propagate on right
        }

        // propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow-1, targetCol, persistenceGrid, ratioVec); //propagate below
        // propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow+1, targetCol, persistenceGrid, ratioVec); //propagate above
        // propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow-1, targetCol-1, persistenceGrid, ratioVec); //propagate diagonal
        // propagate2dParent(parentId, matchingCC, targetRow, targetCol, targetRow-1, targetCol+1, persistenceGrid, ratioVec); //propagate other diagonal

      }
    }
    // else{
    //   if(originCol==targetCol-1){ // look for inclusion
    //       SimplexId parent1DCCId = ccId;
    //         if(parentId==1423){
    //           std::cout<<" look for parent in "<<targetRow<<","<<targetCol<<std::endl;
    //         }
    //         SimplexId parent1DFullId = localToGlobal(parent1DCCId, targetRow, targetCol, nCol);
    //         if(parent1DFullId == -1){ //try with parent
    //           parent1DFullId = localToGlobal(trees_[originRow].getComponent(ccId).parent(), targetRow, targetCol, nCol);
    //         }
    //           if(parentId==1423){
    //             std::cout<<" parent of current cc "<<parentId<<"("<<ccId<<","<<originRow<<","<<originCol<<" is "<<parent1DFullId<<std::endl;
    //           }
    //         if(parent1DFullId!=-1){
    //           if(parentId==1423){
    //             // std::cout<<" parent of current cc "<<parentId<<"("<<globalToLocal_[parentId]<<","<<originRow<<","<<originCol<<" is "<<parent1DCCId<<std::endl;
    //             std::cout<<"set "<<parent1DFullId<<" ( 2d parent of "<<parent1DFullId<<")"<<" as father of 1423"<<std::endl;
    //           }
    //         // if(parents2D_.find(parent1DFullId)!=parents2D_.end()){
    //         fathers_[parentId].push_back(parent1DFullId);
    //         }
    //     }
    // }
}

void RipsMultiParameterSampling::writeCoords(std::string outputFileName,
                                    const SimplexId nb,
                                    std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                    std::vector<std::vector<SimplexId>>& sampling_indexes,
                                    std::vector<SimplexId>& pointsAtScaleStartIdx,
                                    std::vector<std::vector<double>>& fitCoordinates,
                                    int dim,
                                    const float* coordinates) const {


  // SimplexId nbScales = epsilonIdToFilteredId.size();

  SimplexId nbScales = filteredToLocalId.size();
  // std::vector<SimplexId> filteredPointsAtScaleStartIdx(nbScales,0);

  // SimplexId n;
  // for (int iScale=0; iScale<nbScales; iScale++){
  //   SimplexId nbPointsAtScale = filteredToLocalId[iScale].size();
  //   if(iScale+1<nbScales){
  //     filteredPointsAtScaleStartIdx[iScale+1] = filteredPointsAtScaleStartIdx[iScale] + nbPointsAtScale;
  //   }
  //   n+=nbPointsAtScale;
  // }



  printMsg("Writing matrix to disk");
  // SimplexId nbRows = filteredToLocalIdx.size();
  std::ofstream outputFile;
  outputFile.open(outputFileName);
  std::stringstream header;
  header << "x, y, z";
  for(auto i=0; i<dim; i++){
    header<<", param"<<std::setw(5)<<std::setfill('0')<<i;
  }
  outputFile << header.str() << std::endl;


  for (int iScale=0; iScale<nbScales; iScale++){
    for(SimplexId id_fltrd_i=0; id_fltrd_i<filteredToLocalId[iScale].size(); id_fltrd_i++){
      // int iScale = findScaleId(id_fltrd_i, pointsAtScaleStartIdx);

      // SimplexId id_mono_eps_i = id_fltrd_i-pointsAtScaleStartIdx[iScale];

      const SimplexId id_mono_local_i = filteredToLocalId[iScale][id_fltrd_i]; //local id
      const SimplexId id_mono_global_i = sampling_indexes[iScale][id_mono_local_i];

      std::stringstream row;

      row << coordinates[3*id_mono_global_i] << "," << coordinates[3*id_mono_global_i+1] << "," << coordinates[3*id_mono_global_i+2];

      for(SimplexId i=0; i<dim; i++){

        row << ", " << fitCoordinates[iScale][dim*id_mono_local_i+i];

      }

      outputFile << row.str() << std::endl;
    }
  }

  outputFile.close();

}
void RipsMultiParameterSampling::writeMatrixAndCoords(std::string outputFileName,
                                    const SimplexId nb,
                                    std::vector<std::vector<SimplexId>>& epsilonIdToFilteredId,
                                    std::vector<std::vector<SimplexId>>& filteredToLocalId,
                                    std::vector<std::vector<SimplexId>>& sampling_indexes,
                                    std::vector<SimplexId>& pointsAtScaleStartIdx,
                                    std::vector<std::vector<double>>& fitCoordinates,
                                    std::vector<std::vector<float>>& mls3DCoordinates,
                                    std::vector<float>& pointCloudBarycenter,
                                    std::function<double(const double*, const double*, const float*, const float*, const float*, int)> distanceFunction,
                                    int dim,
                                    const float* coordinates) const {


  // SimplexId nbScales = epsilonIdToFilteredId.size();


  size_t n = (size_t)nb*(size_t)(nb-1)/2;
  std::vector<float> distanceMatrix(n,0);

  std::cout<<"compute distance matrix with "<<nb*nb<<" entries"<<std::endl;

  for(LongSimplexId id_multiscale_i=0; id_multiscale_i<nb; id_multiscale_i++){
    int iScale = findScaleId(id_multiscale_i, pointsAtScaleStartIdx);

    SimplexId id_mono_eps_i = id_multiscale_i-pointsAtScaleStartIdx[iScale];

    const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id

    const double* coordinates_i = fitCoordinates[iScale].data();
    const double* coordsi = &coordinates_i[dim*id_mono_local_i];
    const float* mlsCoordinates_i = mls3DCoordinates[iScale].data();



    for(LongSimplexId id_multiscale_j=id_multiscale_i+1; id_multiscale_j<nb; id_multiscale_j++){

      int jScale = findScaleId(id_multiscale_j, pointsAtScaleStartIdx);

      SimplexId id_mono_eps_j = id_multiscale_j-pointsAtScaleStartIdx[jScale];
      const SimplexId id_mono_local_j = filteredToLocalId[jScale][epsilonIdToFilteredId[jScale][id_mono_eps_j]];


      const double* coordinates_j = fitCoordinates[jScale].data();
      const double* coordsj = &coordinates_j[dim*id_mono_local_j];
      const float* mlsCoordinates_j = mls3DCoordinates[jScale].data();

      float dist_comp = (float)distanceFunction(coordsi, coordsj, &mlsCoordinates_i[3*id_mono_local_i], &mlsCoordinates_j[3*id_mono_local_j], pointCloudBarycenter.data(), dim);

      // SimplexId id = id_multiscale_j-1 + id_multiscale_i*(2*nb-id_multiscale_i-3)/2;
      size_t id = (nb*(nb-1)/2) - (nb-id_multiscale_i)*((nb-id_multiscale_i)-1)/2 + id_multiscale_j - id_multiscale_i - 1;
      // std::cout<<" i j "<<id_multiscale_i<<" "<<id_multiscale_j<<std::endl;
      // std::cout<<" nb "<<nb<<std::endl;
      // std::cout<<id_multiscale_i+nb*id_multiscale_j<<" outta "<<distanceMatrix.size()<<std::endl;
      // std::cout<<(size_t)id<<" outta "<<distanceMatrix.size()<<std::endl;
      distanceMatrix[id] = dist_comp;
      // distanceMatrix[(size_t)(id_multiscale_j + nb*id_multiscale_i)] = dist_comp;
    }
  }

  // output matrix to file
  printMsg("Writing matrix to disk");
  // SimplexId nbRows = filteredToLocalIdx.size();
  std::ofstream outputFile;
  outputFile.open(outputFileName);
  std::stringstream header;
  header << "x, y, z";
  for(auto i=0; i<nb; i++){
    header<<", dist"<<std::setw(5)<<std::setfill('0')<<i;
  }
  outputFile << header.str() << std::endl;

  for(SimplexId id_multiscale_i=0; id_multiscale_i<nb; id_multiscale_i++){
    // std::cout<<id_multiscale_i<<"/"<<nb<<std::endl;
    int iScale = findScaleId(id_multiscale_i, pointsAtScaleStartIdx);

    SimplexId id_mono_eps_i = id_multiscale_i-pointsAtScaleStartIdx[iScale];

    const SimplexId id_mono_local_i = filteredToLocalId[iScale][epsilonIdToFilteredId[iScale][id_mono_eps_i]]; //local id
    const SimplexId id_mono_global_i = sampling_indexes[iScale][id_mono_local_i];

    std::stringstream row;

    row << coordinates[3*id_mono_global_i] << "," << coordinates[3*id_mono_global_i+1] << "," << coordinates[3*id_mono_global_i+2];

    for(SimplexId id_multiscale_j=0; id_multiscale_j<nb; id_multiscale_j++){

      int jScale = findScaleId(id_multiscale_j, pointsAtScaleStartIdx);

      SimplexId id_mono_eps_j = id_multiscale_j-pointsAtScaleStartIdx[jScale];
      const SimplexId id_mono_local_j = filteredToLocalId[jScale][epsilonIdToFilteredId[jScale][id_mono_eps_j]];

      // float dist = distanceMatrix[id_multiscale_i+nb*id_mono_local_j];
      float dist = 0;

      if(id_multiscale_i<id_multiscale_j){

        SimplexId id = id_multiscale_j-1 + id_multiscale_i*(2*nb-id_multiscale_i-3)/2;

        // std::cout<<" i "<<id_multiscale_i<<", j "<<id_multiscale_j<<": "<<id<<std::endl;
        // std::cout<<" outta "<<distanceMatrix.size()<<std::endl;

        dist = distanceMatrix[id];

      }else if(id_multiscale_i<id_multiscale_j){

        SimplexId id = id_multiscale_i-1 + id_multiscale_j*(2*nb-id_multiscale_j-3)/2;

        // std::cout<<" i "<<id_multiscale_i<<", j "<<id_multiscale_j<<": "<<id<<std::endl;
        // std::cout<<" outta "<<distanceMatrix.size()<<std::endl;

        dist = distanceMatrix[id];

      }
      row << ", " << dist;

    }

    outputFile << row.str() << std::endl;
  }

  outputFile.close();

}

void ttk::RipsMultiParameterSampling::coordsToDistanceMatrix(float* coords,
                                                int dim,
                                                std::vector<SimplexId> sampling,
                                                std::vector<std::vector<double>>& distMat){

  SimplexId nbPoints = sampling.size();
  distMat.resize(nbPoints, std::vector<double>(nbPoints,0));

  for(SimplexId I=0; I<nbPoints; I++){
    SimplexId i = I;

    for(SimplexId J=I+1; J<nbPoints; J++){
      SimplexId j = J;

      double dist = 0;
      for(int k=0; k<dim; k++){
        dist += (coords[dim*i+k]-coords[dim*j+k])*(coords[dim*i+k]-coords[dim*j+k]);
      }
      dist = sqrt(dist);
      //d
      if(std::isnan(dist)){
        std::cout<<"here is nan "<<I<<" "<<J<<std::endl;
      }

      distMat[I][J] = dist;
      distMat[J][I] = dist;
      
    }
  }
}

