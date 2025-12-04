#include <ScalePDAnalysis.h>
#include <PersistenceDiagramUtils.h>
#include <PersistenceDiagramClustering.h>
#include <string>

using namespace ttk;

ScalePDAnalysis::ScalePDAnalysis() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("ScalePDAnalysis");
}

void ScalePDAnalysis::computeDistanceMatrix(std::vector<std::vector<DiagramType>>& diagrams,
                                            std::vector<std::vector<double>>& distMat){

  if(diagrams.empty() or diagrams[0].empty()){
    printErr("No sequence or empty sequence");
    return;
  }
   size_t nbSequences = diagrams.size();
   size_t nbScales = diagrams[0].size();
   printMsg("Compute dist mat with "+std::to_string(nbSequences)+" sequences");
   distMat.resize(nbSequences);

   for(size_t iSeq=0; iSeq<nbSequences; iSeq++){
     distMat[iSeq].resize(nbSequences,0);
   }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
   for(size_t iSeq=0; iSeq<nbSequences; iSeq++){
     for(size_t jSeq=iSeq+1; jSeq<nbSequences; jSeq++){
       double dist = 0;
       double dist_max = 0;
       for(size_t iScale=0; iScale<nbScales; iScale++){
         DiagramType& diag0 = diagrams[iSeq][iScale];
         DiagramType& diag1 = diagrams[jSeq][iScale];
         double d = computeWassersteinDistance(diag0, diag1);
         if(d>dist_max){
           dist_max = d;
         }
         std::cout<<"scale "<<iScale<<": d="<<d<<std::endl;
         dist += d;
       }
       if(Metric == SCALE_METRIC::MAX){
         dist = dist_max;
       }
       distMat[iSeq][jSeq] = dist;
       distMat[jSeq][iSeq] = dist;
     }
   }
}

// computes wasserstein distance btw two diagramTypes
// using the auction algorithm w/ kdtree
  double ScalePDAnalysis::computeWassersteinDistanceAuction(const DiagramType& D0, const DiagramType& D1, int p) {

    Timer tm{};
      
    BidderDiagram bidderDiag{};
    GoodDiagram goodDiag{};

    // printMsg("Compute Wasserstein distance w/ Auction", 0, 0.0, 1);

    // fill bidder diagram
    for(size_t j = 0; j < D0.size(); j++) {
      // Add bidder to bidders
      Bidder b(D0[j], j, 1.0 /*lambda*/);
      b.setPositionInAuction(bidderDiag.size());
      bidderDiag.emplace_back(b);
    }
    // fill good diagram
    for(size_t j = 0; j < D1.size(); j++) {
      // Add bidder to bidders
      Good g(D1[j], j, 1.0 /*lambda*/);
      g.setPrice(0.0);
      goodDiag.emplace_back(g);
    }

    PersistenceDiagramAuction auction(p,/*alpha*/0.0, /*lambda*/ 1.0, /*deltaLim*/1e-5, true);
    auction.setDebugLevel(5);

    auction.BuildAuctionDiagrams(bidderDiag, goodDiag);

    double d=auction.run();
    printMsg("Compute Wasserstein distance w/ Auction", 1, tm.getElapsedTime(),1);

    return sqrt(d);
  }



double ScalePDAnalysis::computeWassersteinDistance(const DiagramType& D0, const DiagramType& D1) {

  std::vector<DiagramType> diags{};
  std::vector<DiagramType> c{};
  std::vector<std::vector<std::vector<MatchingType>>> m{};

  diags.push_back(D0);
  diags.push_back(D1);

  PersistenceDiagramClustering pdc{};
  pdc.setDebugLevel(debugLevel_);
  pdc.setTimeLimit(TimeLimit);
  pdc.setUseAccelerated(false);
  pdc.setUseKmeansppInit(false);
  pdc.setUseAdditionalPrecision(true);
  pdc.setDeltaLim(0.001);
  // pdc.setForceUseOfAlgorithm(false);
  pdc.execute(diags, c, m);
  double dist = pdc.getWassersteinDistance();
  std::cout<<"DIST "<<dist<<std::endl;
  return dist;

}
