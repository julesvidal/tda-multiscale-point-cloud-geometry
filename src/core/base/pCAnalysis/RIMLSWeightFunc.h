#pragma once

#include <Ponca/Fitting>
// #include <PCAnalysisUtils.h>


namespace ttk::pointCloud {

template<class DataPoint, class WeightKernel>
class RIMLSWeightFunc
{
public:
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;
    using MatrixType = typename DataPoint::MatrixType;

    /*! \brief Return type of the method #w() */
    using WeightReturnType = std::pair<Scalar, VectorType>;

public:
    inline RIMLSWeightFunc(Scalar t = 1, Scalar sigma = 1);

public:
    /*!
     * \brief Initialization method, called by the fitting procedure
     * @param _evalPos Basis center
     */
    PONCA_MULTIARCH inline void init( const VectorType& _evalPos )
    {
        return m_weightFunc.init(_evalPos);
    }



    PONCA_MULTIARCH inline const VectorType& basisCenter() const { return m_weightFunc.basisCenter(); }

    PONCA_MULTIARCH inline VectorType convertToLocalBasis(const VectorType& _q) const;
    PONCA_MULTIARCH inline VectorType convertToGlobalBasis(const VectorType& _q) const;

    inline WeightReturnType w(const VectorType& _q, const DataPoint& attributes) const;

    inline VectorType spacedw(const VectorType& _q, const DataPoint& attributes) const;
    inline MatrixType spaced2w(const VectorType& _q, const DataPoint& attributes) const;

    inline Scalar scaledw(const VectorType& _q, const DataPoint& attributes) const;
    inline Scalar scaled2w(const VectorType& _q, const DataPoint& attributes) const;

    inline VectorType scaleSpaced2w(const VectorType& _q, const DataPoint& attributes) const;

    inline Scalar evalScale() const;

    /*! \brief Access to the evaluation position set during the initialization */
    PONCA_MULTIARCH inline const VectorType & evalPos() const { return m_weightFunc.evalPos(); }

private:
    inline Scalar reweight(const VectorType& _q, const DataPoint& attributes) const;

protected:
    Ponca::DistWeightFunc<DataPoint, WeightKernel> m_weightFunc;
    Scalar m_invSigma2;
    Scalar m_invSigmaScale2;
};

} // namespace pdpc

#include <RIMLSWeightFunc.hpp>
