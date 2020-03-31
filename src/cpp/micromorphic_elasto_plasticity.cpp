/*!
 * micromorphic_elasto_plasticity.cpp
 *
 * An implementation of a elasto-plastic micromorphic constitutive model 
 * following the derivations of Farhad Shahabi in his dissertation.
 */

#include<micromorphic_elasto_plasticity.h>

namespace micromorphicElastoPlasticity{


    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue ){
        /*!
         * Compute the second-order Drucker Prager Yield equation
         *
         * F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
         * 
         * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
         *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
         *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration
         * :param const variableType &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableType &yieldValue: The yield value.
         */

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;
        
        errorOut error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = vectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        return NULL;
    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG ){
        /*!
         * Compute the second-order Drucker Prager Yield equation
         *
         * F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
         * 
         * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
         *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
         *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
         *
         *  Also compute the Jacobians
         * \frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
         * \frac{ \partial F }{ \partial \bar{c} } = -A^{\phi}
         * \frac{ \partial F }{ \partial C_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial C_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial C_{IJ} }
         *
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration
         * :param const variableType &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableType &yieldValue: The yield value.
         * :param variableVector &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
         * :param variableType &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
         * :param variableVector &dFdElasticRCG: The Jacobian of the yield surface w.r.t. the elastic 
         */

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableVector dPressuredStress, dPressuredRCG;
        
        errorOut error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = vectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Evaluate the jacobians
        variableVector devStressDirection = deviatoricReferenceStress / normDevStress;

        dFdStress = vectorTools::Tdot( dDevStressdStress, devStressDirection )
                  + BAngle * dPressuredStress;

        dFdc = - AAngle;

        dFdElasticRCG = vectorTools::Tdot( dDevStressdRCG, devStressDirection )
                      + BAngle * dPressuredRCG;

        return NULL;
    }

    errorOut computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdElasticRCG, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdElasticRCG ){
        /*!
         * Compute the second-order Drucker Prager Yield equation
         *
         * F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
         * 
         * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
         *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
         *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         *  Also compute the Jacobians
         * \frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
         * \frac{ \partial F }{ \partial \bar{c} } = -A^{\phi}
         * \frac{ \partial F }{ \partial C_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure )_{AB} }{ \partial C_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial C_{IJ} }
         *
         *  The second deriatives of \frac{ \partial F }{ \partial \stressMeasure_{IJ}  } are
         *  \frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial stressMeasure_{KL} } = \frac{ \partial^2 || dev( stressMeasure ) || }{ \partial dev( stressMeasure )_{AB} \partial dev( stressMeasure )_{CD} } \frac{ \partial dev( stressMeasure )_{AB} } { \partial stressMeasure_{IJ} } \frac{ \partial dev( stressMeasure )_{CD} } { \partial stressMeasure_{KL} } + \frac{ dev ( stressMeasure )_{AB} }{ || dev( stressMeasure ) || } \frac{ \partial^2 dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} \partial stressMeasure_{KL} }
         *  \frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial RCG_{KL} } = \frac{ \partial^2 || dev( stressMeasure ) || }{ \partial dev( stressMeasure )_{AB} \partial dev( stressMeasure )_{CD} } \frac{ \partial dev( stressMeasure )_{AB} } { \partial stressMeasure_{IJ} } \frac{ \partial dev( stressMeasure )_{CD} } { \partial C_{KL} } + \frac{ dev ( stressMeasure )_{AB} }{ || dev( stressMeasure ) || } \frac{ \partial^2 dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} \partial C_{KL} } + B^{\phi} \frac{ \partial^2 \bar{p} }{ \partial stressMeasure_{IJ} \partial C_{KL} }
         *  
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration
         * :param const variableType &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableType &yieldValue: The yield value.
         * :param variableVector &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
         * :param variableType &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
         * :param variableVector &dFdElasticRCG: The Jacobian of the yield surface w.r.t. the elastic 
         *     right Cauchy-Green deformation tensor.
         * :param variableMatrix &d2FdStress2: The second derivative of the flow direction w.r.t. the stress. This 
         *     is useful if one is using this expression as the flow potential and wants the jacobian of the flow direction \frac{ \partial G }{\partial Stress_{IJ} }
         * :param variableMatrix &d2FdStressdElasticRCG: The second derivative of the flow direction w.r.t. the stress.
         */
        //Assume 3D
        unsigned int dim = 3;

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableType pressure;
        variableVector deviatoricReferenceStress;
        
        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableVector dPressuredStress, dPressuredRCG;

        variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;

        errorOut error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeSecondOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of second-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableType normDevStress = vectorTools::l2norm( deviatoricReferenceStress );

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Evaluate the jacobians
        variableVector devStressDirection = deviatoricReferenceStress / normDevStress;

        dFdStress = vectorTools::Tdot( dDevStressdStress, devStressDirection )
                  + BAngle * dPressuredStress;

        dFdc = - AAngle;

        dFdElasticRCG = vectorTools::Tdot( dDevStressdRCG, devStressDirection )
                      + BAngle * dPressuredRCG;

        //Evaluate the second-order jacobians
        constantMatrix EYE = vectorTools::eye< constantType >( dim * dim );
        variableMatrix dDevStressDirectiondDevStress = ( EYE - vectorTools::dyadic( devStressDirection, devStressDirection ) ) / normDevStress;

        d2FdStress2 = vectorTools::Tdot( dDevStressdStress, vectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ) );

        d2FdStressdElasticRCG = vectorTools::Tdot( vectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ), dDevStressdRCG )
                              + BAngle * d2PressuredStressdRCG;

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int A = 0; A < dim; A++ ){
                            for ( unsigned int B = 0; B < dim; B++ ){
                                d2FdStressdElasticRCG[ dim * I + J ][ dim * K + L ] += devStressDirection[ dim * A + B ] * d2DevStressdStressdRCG[ dim * A + B ][ dim * dim * dim * I + dim * dim * J + dim * K + L ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress, 
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue ){
        /*!
         * Compute the higher-order Drucker Prager Yield equation
         *
         * F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
         * 
         * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }
         * where the K's aren't summed.
         *  dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         * :param const variableVector &referenceHigherOrderStress: The higher-order stress in the reference configuration
         * :param const variableVector &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableVector &yieldValue: The yield value.
         */

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;
        
        errorOut error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation",
                                             "Error in computation of higher-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableVector normDevStress;
        error = micromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation",
                                             "Error in computation of the deviatoric higher-order stress norm" );
            result->addNext( error );
            return result;
        }

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        return NULL;
    }

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress, 
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdElasticRCG ){
        /*!
         * Compute the higher-order Drucker Prager Yield equation
         *
         * F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
         * 
         * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }
         * where the K's aren't summed.
         *  dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         *  Also computes the Jacobians
         *
         * :param const variableVector &referenceHigherOrderStress: The higher-order stress in the reference configuration
         * :param const variableVector &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableVector &yieldValue: The yield value.
         * :param variableMatrix &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
         * :param variableMatrix &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
         * :param variableMatrix &dFdElasticRCG: The Jacobian of the yield function w.r.t. the elastic right Cauchy-Green
         *     deformation tensor.
         */

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableMatrix dPressuredStress, dPressuredRCG;
        
        errorOut error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of higher-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableVector normDevStress;
        variableMatrix dNormDevStressdDevStress;
        error = micromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress, dNormDevStressdDevStress );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (jacobian)",
                                             "Error in computation of the deviatoric higher-order stress norm" );
            result->addNext( error );
            return result;
        }

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Construct the Jacobians
        dFdStress = vectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                  + BAngle * dPressuredStress;

        dFdc = -AAngle * vectorTools::eye< constantType >( cohesion.size() );

        dFdElasticRCG = vectorTools::dot( dNormDevStressdDevStress, dDevStressdRCG )
                      + BAngle * dPressuredRCG;

        return NULL;
    }

    errorOut computeHigherOrderDruckerPragerYieldEquation( const variableVector &referenceHigherOrderStress, 
                                                           const variableVector &cohesion,
                                                           const variableVector &elasticRightCauchyGreen,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdElasticRCG, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdElasticRCG ){
        /*!
         * Compute the higher-order Drucker Prager Yield equation
         *
         * F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
         * 
         * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }
         * where the K's aren't summed.
         *  dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
         *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
         *  A^{angle} = \beta^{angle} \cos( frictionAngle )
         *  B^{angle} = \beta^{angle} \sin( frictionAngle )
         *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
         *
         *  Also computes the Jacobians
         *
         * :param const variableVector &referenceHigherOrderStress: The higher-order stress in the reference configuration
         * :param const variableVector &cohesion: The cohesion measure.
         * :param const variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor.
         * :param const parameterType &frictionAngle: The friction angle
         * :param const parameterType &beta: The beta parameter
         * :param variableVector &yieldValue: The yield value.
         * :param variableMatrix &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
         * :param variableMatrix &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
         * :param variableMatrix &dFdElasticRCG: The Jacobian of the yield function w.r.t. the elastic right Cauchy-Green
         *     deformation tensor.
         * :param variableMatrix &d2FdStress2: The second order Jacobian of the yield function w.r.t. the reference 
         *     higher order stress.
         * :param variableMatrix &d2FdStressdElasticRCG: The second order Jacobian of the yield function w.r.t. the 
         *     reference higher order stress and the elastic right Cauchy-Green Deformation metric.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the parameters
        parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );

        parameterType AAngle = betaAngle * std::cos( frictionAngle );

        parameterType BAngle = betaAngle * std::sin( frictionAngle );

        //Compute the decomposition of the stress
        variableVector pressure;
        variableVector deviatoricReferenceStress;

        variableMatrix dDevStressdStress, dDevStressdRCG;
        variableMatrix dPressuredStress, dPressuredRCG;

        variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;
        
        errorOut error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( referenceHigherOrderStress,
                             elasticRightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                             dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of higher-order stress decomposition" );
            result->addNext( error );
            return result;
        }

        //Compute the l2norm of the deviatoric stress
        variableVector normDevStress;
        variableMatrix dNormDevStressdDevStress;
        variableMatrix d2NormDevStressdDevStress2;
        error = micromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress,
                                                                 dNormDevStressdDevStress,
                                                                 d2NormDevStressdDevStress2 );

        if ( error ){
            errorOut result = new errorNode( "computeHigherOrderDruckerPragerYieldEquation (second order jacobian)",
                                             "Error in computation of the deviatoric higher-order stress norm" );
            result->addNext( error );
            return result;
        }

        //Evaluate the yield equation
        yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        //Construct the Jacobians
        dFdStress = vectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                  + BAngle * dPressuredStress;

        dFdc = -AAngle * vectorTools::eye< constantType >( cohesion.size() );

        dFdElasticRCG = vectorTools::dot( dNormDevStressdDevStress, dDevStressdRCG )
                      + BAngle * dPressuredRCG;

        //Construct the second-order jacobians
        d2FdStress2 = variableMatrix( dim, variableVector( dim * dim * dim * dim * dim * dim, 0 ) );
        d2FdStressdElasticRCG = vectorTools::dot( dNormDevStressdDevStress, d2DevStressdStressdRCG )
                              + BAngle * d2PressuredStressdRCG;

        for ( unsigned int K = 0; K < 3; K++ ){
            for ( unsigned int L = 0; L < 3; L++ ){
                for ( unsigned int M = 0; M < 3; M++ ){
                    for ( unsigned int N = 0; N < 3; N++ ){
                        for ( unsigned int O = 0; O < 3; O++ ){
                            for ( unsigned int P = 0; P < 3; P++ ){
                                for ( unsigned int Q = 0; Q < 3; Q++ ){
                                    for ( unsigned int A = 0; A < 3; A++ ){
                                        for ( unsigned int B = 0; B < 3; B++ ){
                                            for ( unsigned int C = 0; C < 3; C++ ){
                                                for ( unsigned int D = 0; D < 3; D++ ){
                                                    for ( unsigned int E = 0; E < 3; E++ ){
                                                        d2FdStressdElasticRCG[ K ][ dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ]
                                                            += d2NormDevStressdDevStress2[ K ][ dim * dim * dim * dim * dim * Q + dim * dim * dim * dim * A + dim * dim * dim * B + dim * dim * C + dim * D + E ]
                                                             * dDevStressdStress[ dim * dim * Q + dim * A + B ][ dim * dim * L + dim * M + N ]
                                                             * dDevStressdRCG[ dim * dim * C + dim * D + E ][ dim * O + P ];
                                                        for ( unsigned int F = 0; F < 3; F++ ){
                                                            d2FdStress2[ K ][ dim * dim * dim * dim * dim * L + dim * dim * dim * dim * M + dim * dim * dim * N + dim * dim * O + dim * P + Q ]
                                                                += d2NormDevStressdDevStress2[ K ][ dim * dim * dim * dim * dim * A + dim * dim * dim * dim * B + dim * dim * dim * C + dim * dim * D + dim * E + F ]
                                                                 * dDevStressdStress[ dim * dim * A + dim * B + C ][ dim * dim * L + dim * M + N ]
                                                                 * dDevStressdStress[ dim * dim * D + dim * E + F ][ dim * dim * O + dim * P + Q ];
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &inversePlasticDeformationGradient,
                                              variableVector &inversePlasticMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation, 
                                              variableVector &elasticGradientMicroDeformation ){
        /*!
         * Compute the elastic parts of the various deformation measures.
         *
         * F_{i\bar{I}}^e = F_{iI} F_{I \bar{I}}^{p, -1}
         * \chi_{i\bar{I}}^e = \chi_{iI} \chi_{I \bar{I}}^{p, -1}
         * \chi_{ k\bar{ K },\bar{ L } } = \left( \chi_{ kK, L } F_{ L \bar{ L } }^{p, -1} - \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \right) \chi_{ K \bar{ K } }^{p,-1}
         *
         * :param const variableVector &deformationGradient: The macroscale deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation tensor $\chi$
         * :param const variableVector &gradientMicroDeformation: The gradient of the micro-deformation tensor
         *     $\chi$ with respect to the reference configuration.
         * :param const variableVector &plasticDeformationGradient: The plastic part of the macroscale 
         *     deformation gradient.
         * :param const variableVector &plasticMicroDeformation: The plastic part of the micro-deformation tensor
         *     $\chi$.
         * :param const variableVector &plasticGradientMicroDeformation: The plastic part of the gradient of the 
         *     micro-deformation tensor $\chi$ with respect to the intermediate configuration.
         * :param variableVector &inversePlasticDeformationGradient: The inverse of the plastic part of the macro deformation gradient.
         * :param variableVector &inversePlasticMicroDeformation: The inverse of the plastic part of the micro deformation.
         * :param variableVector &elasticDeformationGradient: The elastic part of the macroscale deformation gradient.
         * :param variableVector &elasticMicroDeformation: The elastic part of the micro-deformation tensor $\chi$
         * :param variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the micro-deformation
         *     tensor $\chi$ w.r.t. the reference configuration.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The deformation gradient is not the correct size (3D)" );
        }

        if ( microDeformation.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The micro-deformation is not the correct size (3D)" );
        }

        if ( gradientMicroDeformation.size() != dim * dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The gradient of the micro-deformation is not the correct size (3D)" );
        }

        if ( plasticDeformationGradient.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The plastic deformation gradient is not the correct size (3D)" );
        }

        if ( plasticMicroDeformation.size() != dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The plastic micro-deformation is not the correct size (3D)" );
        }

        if ( plasticGradientMicroDeformation.size() != dim * dim * dim ){
            return new errorNode( "computeElasticPartOfDeformation",
                                  "The plastic gradient of the micro-deformation is not the correct size (3D)" );
        }

        //Compute the inverses of the plastic deformation gradient and micro-deformation
        inversePlasticDeformationGradient = vectorTools::inverse( plasticDeformationGradient, dim, dim );

        inversePlasticMicroDeformation = vectorTools::inverse( plasticMicroDeformation, dim, dim );

        //Assemble the elastic parts of the deformation measures
        elasticDeformationGradient = vectorTools::matrixMultiply( deformationGradient, inversePlasticDeformationGradient,
                                                                  dim, dim, dim, dim );

        elasticMicroDeformation = vectorTools::matrixMultiply( microDeformation, inversePlasticMicroDeformation,
                                                               dim, dim, dim, dim );

        elasticGradientMicroDeformation = variableVector( dim * dim * dim, 0 );
        for ( unsigned int k = 0; k < dim; k++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int Ab = 0; Ab < dim; Ab++ ){
                            elasticGradientMicroDeformation[ dim * dim * k + dim * Kb + Lb ]
                                += ( gradientMicroDeformation[ dim * dim * k + dim * K + Ab ]
                                 *   inversePlasticDeformationGradient[ dim * Ab + Lb ]
                                 -   elasticMicroDeformation[ dim * k + Ab ]
                                 *   plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ] )
                                 * inversePlasticMicroDeformation[ dim * K + Kb ];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation, 
                                              variableVector &elasticGradientMicroDeformation ){
        /*!
         * Compute the elastic parts of the various deformation measures.
         * F_{i\bar{I}}^e = F_{iI} F_{I \bar{I}}^{p, -1}
         * \chi_{i\bar{I}}^e = \chi_{iI} \chi_{I \bar{I}}^{p, -1}
         * \chi_{ k\bar{ K },\bar{ L } } = \left( \chi_{ kK, L } F_{ L \bar{ L } }^{p, -1} - \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \right) \chi_{ K \bar{ K } }^{p,-1}
         *
         * :param const variableVector &deformationGradient: The macroscale deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation tensor $\chi$
         * :param const variableVector &gradientMicroDeformation: The gradient of the micro-deformation tensor
         *     $\chi$ with respect to the reference configuration.
         * :param const variableVector &plasticDeformationGradient: The plastic part of the macroscale 
         *     deformation gradient.
         * :param const variableVector &plasticMicroDeformation: The plastic part of the micro-deformation tensor
         *     $\chi$.
         * :param const variableVector &plasticGradientMicroDeformation: The plastic part of the gradient of the 
         *     micro-deformation tensor $\chi$ with respect to the intermediate configuration.
         * :param variableVector &elasticDeformationGradient: The elastic part of the macroscale deformation gradient.
         * :param variableVector &elasticMicroDeformation: The elastic part of the micro-deformation tensor $\chi$
         * :param variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the micro-deformation
         *     tensor $\chi$ w.r.t. the reference configuration.
         */

        variableVector inversePlasticDeformationGradient, inversePlasticMicroDeformation;

        return computeElasticPartOfDeformation( deformationGradient, microDeformation, gradientMicroDeformation,
                                                plasticDeformationGradient, plasticMicroDeformation, plasticGradientMicroDeformation,
                                                inversePlasticDeformationGradient, inversePlasticMicroDeformation,
                                                elasticDeformationGradient, elasticMicroDeformation, elasticGradientMicroDeformation );
    }

    errorOut computeElasticPartOfDeformation( const variableVector &deformationGradient, const variableVector &microDeformation,
                                              const variableVector &gradientMicroDeformation,
                                              const variableVector &plasticDeformationGradient,
                                              const variableVector &plasticMicroDeformation,
                                              const variableVector &plasticGradientMicroDeformation,
                                              variableVector &elasticDeformationGradient, variableVector &elasticMicroDeformation,
                                              variableVector &elasticGradientMicroDeformation,
                                              variableMatrix &dElasticFdF, variableMatrix &dElasticFdPlasticF,
                                              variableMatrix &dElasticChidChi, variableMatrix &dElasticChidPlasticChi,
                                              variableMatrix &dElasticGradChidGradChi, variableMatrix &dElasticGradChidPlasticGradChi,
                                              variableMatrix &dElasticGradChidPlasticF, variableMatrix &dElasticGradChidChi,
                                              variableMatrix &dElasticGradChidPlasticChi ){
        /*!
         * Compute the elastic parts of the various deformation measures.
         * F_{i\bar{I}}^e = F_{iI} F_{I \bar{I}}^{p, -1}
         * \chi_{i\bar{I}}^e = \chi_{iI} \chi_{I \bar{I}}^{p, -1}
         * \chi_{ k\bar{ K },\bar{ L } } = \left( \chi_{ kK, L } F_{ L \bar{ L } }^{ p, -1}  - \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \right) \chi_{ K \bar{ K } }^{p,-1}
         *
         * Also compute the Jacobians
         *
         * \frac{ \partial F_{i\bar{I}}^e }{ \partial F_{nN} } = \delta_{in} F_{N \bar{I} }^{p, -1}
         * \frac{ \partial F_{i\bar{I}}^e }{ \partial F_{\bar{N} N}^p } = -F_{iI} F_{I \bar{N}}^{p, -1} F_{N \bar{I} }^{p, -1}
         * \frac{ \partial \chi_{i\bar{I}}^e }{ \partial \chi_{nN} } = \delta_{in} \chi_{N \bar{I} }^{p, -1}
         * \frac{ \partial \chi_{i\bar{I}}^e }{ \partial \chi_{\bar{N} N}^p } = -\chi_{iI} \chi_{I \bar{N}}^{p, -1} \chi_{N \bar{I} }^{p, -1}
         * \frac{ \partial \chi_{k \bar{ K }, \bar{ L } } }{ \partial F_{ \bar{N} N }^{ p } } = \chi_{ kK, L } F_{ L \bar{ N } }^{p, -1} F_{ N \bar{ L } }^{p, -1}
         * \frac{ \partial \chi_{k\bar{ K },\bar{ L } } }{ \partial \chi_{ nN } } = - \frac{ \partial \chi_{k \bar{ A } }^e }{ \partial \chi_{ nN } } \chi_{ \bar{ A } K,\bar{ L } }^{ p } \chi_{ K \bar{ K } }^{p,-1}
         * \frac{ \partial \chi_{k\bar{ K },\bar{ L } } }{ \partial \chi_{ \bar{ N } N } } = \chi_{ k \bar{ A } }^e \chi_{ \bar{ A } K,\bar{ L } }^{ p } \chi_{ K \bar{ N } }^{p,-1} \chi_{ N \bar{ K } }^{p,-1}
         *
         * :param const variableVector &deformationGradient: The macroscale deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation tensor $\chi$
         * :param const variableVector &gradientMicroDeformation: The gradient of the micro-deformation tensor
         *     $\chi$ with respect to the reference configuration.
         * :param const variableVector &plasticDeformationGradient: The plastic part of the macroscale 
         *     deformation gradient.
         * :param const variableVector &plasticMicroDeformation: The plastic part of the micro-deformation tensor
         *     $\chi$.
         * :param const variableVector &plasticGradientMicroDeformation: The plastic part of the gradient of the 
         *     micro-deformation tensor $\chi$ with respect to the intermediate configuration.
         * :param variableVector &elasticDeformationGradient: The elastic part of the macroscale deformation gradient.
         * :param variableVector &elasticMicroDeformation: The elastic part of the micro-deformation tensor $\chi$
         * :param variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the micro-deformation
         *     tensor $\chi$ w.r.t. the intermediate configuration.
         * :param variableMatrix &dElasticFdF: The Jacobian of the elastic part of the deformation gradient w.r.t. the deformation 
         *     gradient.
         * :param variableMatrix &dElasticFdPlasticF: The Jacobian of the elastic part of the deformation gradient w.r.t. the 
         *     plastic part of the deformation gradient.
         * :param variableMatrix &dElasticChidChi: The Jacobian of the elastic part of the micro-deformation w.r.t. the
         *     micro deformation.
         * :param variableMatrix &dElasticChidPlasticChi: The Jacobian of the elastic part of the micro-deformation w.r.t.
         *     the plastic part of the micro deformation.
         * :param variableMatrix &dElasticGradChidGradChi: The Jacobian of the elastic part of the gradient of the micro-deformation
         *     w.r.t. the gradient of the micro-deformation.
         * :param variableMatrix &dElasticGradChidPlasticGradChi: The Jacobian of the elastic part of the gradient of the 
         *     micro-deformation w.r.t. the plastic part of the gradient of the micro-deformation.
         * :param variableMatrix &dElasticGradChidPlasticF: The Jacobian of the elastic part of the gradient of the micro-deformation
         *     w.r.t. the plastic deformation gradient.
         * :param variableMatrix &dElasticGradChidChi: The Jacobian of the elastic part of the gradient of the micro-deformation 
         *     w.r.t. the micro deformation.
         * :param variableMatrix &dElasticGradChidPlasticChi: The Jacobian of the elastic part of the gradient of the micro-deformation
         *     w.r.t. the plastic part of the micro-deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the required deformation measures
        variableVector inversePlasticDeformationGradient, inversePlasticMicroDeformation;

        errorOut error = computeElasticPartOfDeformation( deformationGradient, microDeformation, gradientMicroDeformation,
                                                          plasticDeformationGradient, plasticMicroDeformation, 
                                                          plasticGradientMicroDeformation, inversePlasticDeformationGradient,
                                                          inversePlasticMicroDeformation, elasticDeformationGradient,
                                                          elasticMicroDeformation, elasticGradientMicroDeformation );

        if ( error ){
            errorOut result = new errorNode( "computeElasticPartOfDeformation (jacobian)",
                                             "Error in computation of the elastic part of the deformation measures" );
            result->addNext( error );
            return result;
        }

        //Assemble the Jacobians
        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dElasticFdF = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dElasticFdPlasticF = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dElasticChidChi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dElasticChidPlasticChi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int Ib = 0; Ib < dim; Ib++ ){
                for ( unsigned int n = 0; n < dim; n++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        dElasticFdF[ dim * i + Ib ][ dim * n + N ] = eye[ dim * i + n ] * inversePlasticDeformationGradient[ dim * N + Ib ];
                        dElasticChidChi[ dim * i + Ib ][ dim * n + N ] = eye[ dim * i + n ] * inversePlasticMicroDeformation[ dim * N + Ib ];

                        for ( unsigned int I = 0; I < dim; I++ ){
                            dElasticFdPlasticF[ dim * i + Ib ][ dim * n + N ] -= deformationGradient[ dim * i + I ]
                                                                               * inversePlasticDeformationGradient[ dim * I + n ]
                                                                               * inversePlasticDeformationGradient[ dim * N + Ib ];
                            dElasticChidPlasticChi[ dim * i + Ib ][ dim * n + N ] -= microDeformation[ dim * i + I ]
                                                                                   * inversePlasticMicroDeformation[ dim * I + n ]
                                                                                   * inversePlasticMicroDeformation[ dim * N + Ib ];
                        }
                    }
                }
            }
        }

        dElasticGradChidPlasticF = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dElasticGradChidChi = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dElasticGradChidPlasticChi = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        dElasticGradChidGradChi = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dElasticGradChidPlasticGradChi = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int k = 0; k < dim; k++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int n = 0; n < dim; n++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            for ( unsigned int K = 0; K < dim; K++ ){

                                dElasticGradChidGradChi[ dim * dim * k + dim * Kb + Lb ][ dim * dim * n + dim * N + K ]
                                    = eye[ dim * k + n ] * inversePlasticDeformationGradient[ dim * K + Lb ]
                                    * inversePlasticMicroDeformation[ dim * N + Kb ];
                                
                                dElasticGradChidPlasticGradChi[ dim * dim * k + dim * Kb + Lb ][ dim * dim * n + dim * N + K ]
                                    = -elasticMicroDeformation[ dim * k + n ] * eye[ dim * Lb + K ]
                                    * inversePlasticMicroDeformation[ dim * N + Kb ];

                                for ( unsigned int Ab = 0; Ab < dim; Ab++ ){
                                    dElasticGradChidPlasticF[ dim * dim * k + dim * Kb + Lb ][ dim * n + N ]
                                        -= gradientMicroDeformation[ dim * dim * k + dim * K + Ab ]
                                         * inversePlasticDeformationGradient[ dim * Ab + n ]
                                         * inversePlasticDeformationGradient[ dim * N + Lb ]
                                         * inversePlasticMicroDeformation[ dim * K + Kb ];

                                    dElasticGradChidChi[ dim * dim * k + dim * Kb + Lb ][ dim * n + N ]
                                        -= dElasticChidChi[ dim * k + Ab ][ dim * n + N ]
                                        *  plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ]
                                        *  inversePlasticMicroDeformation[ dim * K + Kb ];

                                    dElasticGradChidPlasticChi[ dim * dim * k + dim * Kb + Lb ][ dim * n + N ]
                                        -= dElasticChidPlasticChi[ dim * k + Ab ][ dim * n + N ]
                                        *  plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ]
                                        *  inversePlasticMicroDeformation[ dim * K + Kb]
                                        +  ( gradientMicroDeformation[ dim * dim * k + dim * K + Ab ]
                                        *    inversePlasticDeformationGradient[ dim * Ab + Lb ]
                                        -    elasticMicroDeformation[ dim * k + Ab ]
                                        *    plasticGradientMicroDeformation[ dim * dim * Ab + dim * K + Lb ] )
                                        *    inversePlasticMicroDeformation[ dim * K + n ]
                                        *    inversePlasticMicroDeformation[ dim * N + Kb ] ;
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeElasticDeformationMeasures( const variableVector &elasticDeformationGradient,
                                                const variableVector &elasticMicroDeformation,
                                                const variableVector &elasticGradientMicroDeformation,
                                                variableVector &elasticRightCauchyGreen,
                                                variableVector &elasticMicroRightCauchyGreen,
                                                variableVector &elasticPsi, variableVector &elasticGamma ){
        /*!
         * Compute the elastic deformation measures required for the evolution of plasticity.
         *
         * :param const variableVector &elasticDeformationGradient: The elastic part of the deformation gradient.
         * :param const variableVector &elasticMicroDeformation: The elastic part of the micro-deformation.
         * :param const variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the
         *     micro-deformation.
         * :param variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation measure.
         * :param variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation 
         *     measure.
         * :param variableVector &elasticPsi: The elastic part of the micro-deformation measure Psi.
         * :param variableVector &elasticGamma: The elastic part of the higher order deformation measure.
         */

        errorOut error = constitutiveTools::computeRightCauchyGreen( elasticDeformationGradient, elasticRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::computeRightCauchyGreen( elasticMicroDeformation, elasticMicroRightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic micro right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computePsi( elasticDeformationGradient, elasticMicroDeformation, elasticPsi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic micro-deformation metric Psi" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computeGamma( elasticDeformationGradient, elasticGradientMicroDeformation, elasticGamma );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures",
                                             "Error in computation of the elastic higher order deformation metric Gamma" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeElasticDeformationMeasures( const variableVector &elasticDeformationGradient,
                                                const variableVector &elasticMicroDeformation,
                                                const variableVector &elasticGradientMicroDeformation,
                                                variableVector &elasticRightCauchyGreen,
                                                variableVector &elasticMicroRightCauchyGreen,
                                                variableVector &elasticPsi, variableVector &elasticGamma,
                                                variableMatrix &dElasticRCGdElasticF, variableMatrix &dElasticMicroRCGdElasticChi,
                                                variableMatrix &dElasticPsidElasticF, variableMatrix &dElasticPsidElasticChi,
                                                variableMatrix &dElasticGammadElasticF, variableMatrix &dElasticGammadElasticGradChi ){
        /*!
         * Compute the elastic deformation measures required for the evolution of plasticity.
         *
         * :param const variableVector &elasticDeformationGradient: The elastic part of the deformation gradient.
         * :param const variableVector &elasticMicroDeformation: The elastic part of the micro-deformation.
         * :param const variableVector &elasticGradientMicroDeformation: The elastic part of the gradient of the
         *     micro-deformation.
         * :param variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation measure.
         * :param variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation 
         *     measure.
         * :param variableVector &elasticPsi: The elastic part of the micro-deformation measure Psi.
         * :param variableVector &elasticGamma: The elastic part of the higher order deformation measure.
         * :param variableMatrix &dElasticRCGdElasticF: The Jacobian of the right Cauchy-Green deformation measure
         *     w.r.t. the elastic deformation gradient.
         * :param variableMatrix &dElasticMicroRCGdElasticChi: The Jacobian of the micro right Cauchy-Green deformation 
         *     measure w.r.t. the elastic micro deformation
         * :param variableMatrix &dElasticPsidElasticF: The Jacobian of the elastic micro-deformation measure Psi
         *     w.r.t. the elastic deformation gradient.
         * :param variableMatrix &dElasticPsidElasticF: The Jacobian of the elastic micro-deformation measure Psi
         *     w.r.t. the elastic micro-deformation.
         * :param variableMatrix &dElasticGammadElasticF: The Jacobian of the higher order deformation measure Gamma w.r.t.
         *     the elastic deformation gradient.
         * :param variableMatrix &dElasticGammadElasticGradChi: The Jacobian of the higher order deformation measure Gamma
         *     w.r.t. the elastic part of the gradient of the micro-deformation.
         */

        errorOut error = constitutiveTools::computeRightCauchyGreen( elasticDeformationGradient, elasticRightCauchyGreen,
                                                                     dElasticRCGdElasticF );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::computeRightCauchyGreen( elasticMicroDeformation, elasticMicroRightCauchyGreen,
                                                            dElasticMicroRCGdElasticChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic micro right Cauchy-Green deformation tensor" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computePsi( elasticDeformationGradient, elasticMicroDeformation, elasticPsi,
                                               dElasticPsidElasticF, dElasticPsidElasticChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic micro-deformation metric Psi" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computeGamma( elasticDeformationGradient, elasticGradientMicroDeformation, elasticGamma,
                                                 dElasticGammadElasticF, dElasticGammadElasticGradChi );

        if ( error ){
            errorOut result = new errorNode( "computeElasticDeformationMeasures (jacobian)",
                                             "Error in computation of the elastic higher order deformation metric Gamma" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient ){
        /*!
         * Compute the plastic macro velocity gradient in the intermediate configuration.
         *
         * \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( inverseElasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticMacroVelocityGradient",
                                  "The inverse elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( macroFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMacroVelocityGradient",
                                  "The macro flow direction tensor must be 3D" );
        }

        if ( microFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMacroVelocityGradient",
                                  "The micro flow direction tensor must be 3D" );
        }

        //Compute the macro-scale velocity gradient
        plasticMacroVelocityGradient = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    plasticMacroVelocityGradient[ dim * Bb + Kb ]
                        += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                         * ( macroGamma * macroFlowDirection[ dim * Kb + Lb ]
                         +   microGamma * microFlowDirection[ dim * Kb + Lb ] );
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma ){
        /*!
         * Compute the plastic macro velocity gradient in the intermediate configuration.
         *
         * \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMacroVelocityGradient (plastic multiplier jacobian)",
                                             "Error in computation of the plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        dPlasticMacroLdMacroGamma = variableVector( dim * dim, 0 );
        dPlasticMacroLdMicroGamma = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    dPlasticMacroLdMacroGamma[ dim * Bb + Kb ] += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                                                                * macroFlowDirection[ dim * Kb + Lb ];

                    dPlasticMacroLdMicroGamma[ dim * Bb + Kb ] += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                                                                * microFlowDirection[ dim * Kb + Lb ];
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma, 
                                                  variableMatrix &dPlasticMacroLdElasticRCG, 
                                                  variableMatrix &dPlasticMacroLdMacroFlowDirection, 
                                                  variableMatrix &dPlasticMacroLdMicroFlowDirection ){
        /*!
         * Compute the plastic macro velocity gradient in the intermediate configuration.
         *
         * \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableMatrix &dPlasticMacroLdElasticRCG: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     elastic right Cauchy-Green deformation tensor.
         * :param variableMatrix &dPlasticMacroLdMacroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro flow direction tensor.
         * :param variableMatrix &dPlasticMacroLdMicroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro flow direction tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                              dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMacroVelocityGradient (full jacobian)",
                                             "Error in computation of the plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dPlasticMacroLdElasticRCG = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMacroLdMacroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMacroLdMicroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                    for ( unsigned int Pb = 0; Pb < dim; Pb++ ){
                        dPlasticMacroLdElasticRCG[ dim * Bb + Kb ][ dim * Ob + Pb ]
                            -= inverseElasticRightCauchyGreen[ dim * Bb + Ob ]
                             * plasticMacroVelocityGradient[ dim * Pb + Kb ];
                        
                        dPlasticMacroLdMacroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                            += macroGamma * inverseElasticRightCauchyGreen[ dim * Bb + Pb ] * eye[ dim * Kb + Ob ];

                        dPlasticMacroLdMicroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                            += microGamma * inverseElasticRightCauchyGreen[ dim * Bb + Pb ] * eye[ dim * Kb + Ob ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient ){
        /*!
         * Compute the plastic micro velocity gradient
         *
         *  \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
         *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
         *
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
         * :param const variableVector &microFlowDirection: The micro plastic flow direction.
         * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticMicroRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The elastic micro right Cauchy-Green deformation tensor is not 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The elastic micro deformation tensor Psi is not 3D" );
        }

        if ( inverseElasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The inverse of the elastic micro deformation tensor Psi is not 3D" );
        }

        if ( microFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient",
                                  "The micro flow direction of the elastic micro plastic flow direction is not 3D" );
        }

        plasticMicroVelocityGradient = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int Nb = 0; Nb < dim; Nb++ ){
                        for ( unsigned int Eb = 0; Eb < dim; Eb++ ){
                            plasticMicroVelocityGradient[ dim * Bb + Kb ]
                                += microGamma
                                 * inverseElasticPsi[ dim * Bb + Lb ]
                                 * microFlowDirection[ dim * Eb + Lb ]
                                 * elasticPsi[ dim * Nb + Eb ]
                                 * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma ){
        /*!
         * Compute the plastic micro velocity gradient
         *
         *  \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
         *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
         *
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
         * :param const variableVector &microFlowDirection: The micro plastic flow direction.
         * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro plastic multiplier.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticMicroRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The elastic micro right Cauchy-Green deformation tensor is not 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The elastic micro deformation tensor Psi is not 3D" );
        }

        if ( inverseElasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The inverse of the elastic micro deformation tensor Psi is not 3D" );
        }

        if ( microFlowDirection.size() != dim * dim ){
            return new errorNode( "computePlasticMicroVelocityGradient (jacobian)",
                                  "The micro flow direction of the elastic micro plastic flow direction is not 3D" );
        }

        plasticMicroVelocityGradient = variableVector( dim * dim, 0 );
        dPlasticMicroLdMicroGamma = variableVector( dim * dim, 0 );

        for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                    for ( unsigned int Nb = 0; Nb < dim; Nb++ ){
                        for ( unsigned int Eb = 0; Eb < dim; Eb++ ){
                            dPlasticMicroLdMicroGamma[ dim * Bb + Kb ]
                                += inverseElasticPsi[ dim * Bb + Lb ]
                                 * microFlowDirection[ dim * Eb + Lb ]
                                 * elasticPsi[ dim * Nb + Eb ]
                                 * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];
                        }
                    }
                }
                plasticMicroVelocityGradient[ dim * Bb + Kb ] = microGamma * dPlasticMicroLdMicroGamma[ dim * Bb + Kb ];
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma,
                                                  variableMatrix &dPlasticMicroLdElasticMicroRCG,
                                                  variableMatrix &dPlasticMicroLdElasticPsi,
                                                  variableMatrix &dPlasticMicroLdMicroFlowDirection ){
        /*!
         * Compute the plastic micro velocity gradient
         *
         *  \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]
         *
         *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
         *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
         *
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
         * :param const variableVector &microFlowDirection: The micro plastic flow direction.
         * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro plastic multiplier.
         * :param variableMatrix &dPlasticMicroLdElasticMicroRCG: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro right Cauchy-Green deformation tensor.
         * :param variableMatrix &dPlasticMicroLdElasticPsi: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro deformation measure Psi.
         * :param variableMatrix &dPlasticMicroLdMicroFlowDirection: The Jacobian of the plastic micro velocity gradient
         *     w.r.t. the micro flow direction.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen,
                                                              elasticPsi, inverseElasticPsi, microFlowDirection,
                                                              plasticMicroVelocityGradient, dPlasticMicroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMicroVelocityGradient (full jacobian)",
                                             "Error in computation of the plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        //Assemble the Jacobians
        dPlasticMicroLdElasticMicroRCG = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroLdElasticPsi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroLdMicroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int Eb = 0; Eb < dim; Eb++ ){
            for ( unsigned int Fb = 0; Fb < dim; Fb++ ){
                for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                    for ( unsigned int Pb = 0; Pb < dim; Pb++ ){
                        dPlasticMicroLdElasticPsi[ dim * Eb + Fb ][ dim * Ob + Pb ]
                            -= inverseElasticPsi[ dim * Eb + Ob ] * plasticMicroVelocityGradient[ dim * Pb + Fb ];

                        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                            dPlasticMicroLdElasticPsi[ dim * Eb + Fb ][ dim * Ob + Pb ]
                                += microGamma * inverseElasticPsi[ dim * Eb + Lb ]
                                 * microFlowDirection[ dim * Pb + Lb ]
                                 * elasticMicroRightCauchyGreen[ dim * Ob + Fb ];

                            dPlasticMicroLdMicroFlowDirection[ dim * Eb + Fb ][ dim * Ob + Pb ]
                                += microGamma * inverseElasticPsi[ dim * Eb + Pb ]
                                 * elasticPsi[ dim * Lb + Ob ]
                                 * elasticMicroRightCauchyGreen[ dim * Lb + Fb ];

                            for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                                dPlasticMicroLdElasticMicroRCG[ dim * Eb + Fb ][ dim * Ob + Pb ]
                                    += microGamma * inverseElasticPsi[ dim * Eb + Lb ]
                                     * microFlowDirection[ dim * Kb + Lb ] * elasticPsi[ dim * Ob + Kb ]
                                     * eye[ dim * Fb + Pb ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableMatrix &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         */

        variableVector skewTerm;
        return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                            elasticGamma, microGradientFlowDirection,
                                                            plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                            skewTerm );

    }
    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableMatrix &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableVector &skewTerm: The skew term ( times 2 ) from the higher order computation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( microGradientGamma.size() != dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The micro gradient plastic multiplier must have a length of 3" );
        }

        if ( elasticPsi.size() != dim  * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The elastic micro deformation measure Psi must be 3D" );
        }

        if ( inverseElasticPsi.size() != dim  * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The inverse elastic micro deformation measure Psi must be 3D" );
        }

        if ( elasticGamma.size() != dim * dim * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The elastic higher order deformation measure Gamma must be 3D" );
        }

        if ( microGradientFlowDirection.size() != dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The micro gradient flow direction must be 3D" );
        }

        for ( unsigned int i = 0; i < dim; i++ ){
            if ( microGradientFlowDirection[ i ].size() != dim * dim * dim ){
                return new errorNode( "computePlasticVelocityGradients",
                                      "The rows of the micro gradient flow direction must be of length 27" );
            }
        }

        if ( plasticMicroVelocityGradient.size() != dim * dim ){
            return new errorNode( "computePlasticMicroGradientVelocityGradient",
                                  "The plastic micro velocity gradient must be 3D" );
        }

        //Assemble the 'skew' term
        skewTerm = variableVector( dim * dim * dim, 0 );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Cb = 0; Cb < dim; Cb++ ){
                        for ( unsigned int Fb = 0; Fb < dim; Fb++ ){
                            skewTerm[ dim * dim * Db + dim * Mb + Kb ]
                                += plasticMicroVelocityGradient[ dim * Db + Cb ]
                                 * inverseElasticPsi[ dim * Cb + Fb ]
                                 * elasticGamma[ dim * dim * Fb + dim * Mb + Kb ]
                                 - plasticMicroVelocityGradient[ dim * Cb + Mb ]
                                 * inverseElasticPsi[ dim * Db + Fb ]
                                 * elasticGamma[ dim * dim * Fb + dim * Cb + Kb ];
                        }
                    }
                }
            }
        }

        plasticMicroGradientVelocityGradient = variableVector( dim * dim * dim, 0 );

        for ( unsigned int Nb = 0; Nb < dim; Nb++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                        for ( unsigned int Ib = 0; Ib < dim; Ib++ ){
                            plasticMicroGradientVelocityGradient[ dim * dim * Nb + dim * Mb + Kb ]
                                += inverseElasticPsi[ dim * Nb + Lb ]
                                 * ( microGradientGamma[ Ib ] * microGradientFlowDirection[ Ib ][ dim * dim * Kb + dim * Lb + Mb ]
                                 +   elasticPsi[ dim * Lb + Ib ] * skewTerm[ dim * dim * Ib + dim * Mb + Kb ] );
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableMatrix &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient gamma.
         * :param variableMatrix &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the platic micro velocity gradient.
         */

        variableVector skewTerm;
        return computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi, elasticGamma,
                                                            microGradientFlowDirection, plasticMicroVelocityGradient,
                                                            plasticMicroGradientVelocityGradient, skewTerm,
                                                            dPlasticMicroGradientLdMicroGradientGamma,
                                                            dPlasticMicroGradientLdPlasticMicroL );
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableMatrix &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableVector &skewTerm,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableVector &skewTerm: Two times the skew term.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient gamma.
         * :param variableMatrix &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the platic micro velocity gradient.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                                      elasticGamma, microGradientFlowDirection,
                                                                      plasticMicroVelocityGradient, 
                                                                      plasticMicroGradientVelocityGradient, skewTerm );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMicroGradientVelocityGradient (jacobian)",
                                             "Error in computation of the plasticMicroVelocityGradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dPlasticMicroGradientLdPlasticMicroL = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroGradientLdMicroGradientGamma = variableMatrix( dim * dim * dim, variableVector( dim, 0 ) );

        for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                            dPlasticMicroGradientLdMicroGradientGamma[ dim * dim * Lb + dim * Mb + Kb ][ Ob ]
                                += inverseElasticPsi[ dim * Lb + Pb ]
                                 * microGradientFlowDirection[ Ob ][ dim * dim * Kb + dim * Pb + Mb ];

                            for ( unsigned int Qb = 0; Qb < dim; Qb++ ){
                                dPlasticMicroGradientLdPlasticMicroL[ dim * dim * Lb + dim * Mb + Kb ][ dim * Ob + Pb ]
                                    += eye[ dim * Lb + Ob ] * inverseElasticPsi[ dim * Pb + Qb ]
                                     * elasticGamma[ dim * dim * Qb + dim * Mb + Kb ]
                                     - eye[ dim * Mb + Pb ] * inverseElasticPsi[ dim * Lb + Qb ]
                                     * elasticGamma[ dim * dim * Qb + dim * Ob + Kb ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticMicroGradientVelocityGradient( const variableVector &microGradientGamma, const variableVector &elasticPsi,
                                                          const variableVector &inverseElasticPsi, const variableVector &elasticGamma,
                                                          const variableMatrix &microGradientFlowDirection,
                                                          const variableVector &plasticMicroVelocityGradient,
                                                          variableVector &plasticMicroGradientVelocityGradient,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                                          variableMatrix &dPlasticMicroGradientLdPlasticMicroL,
                                                          variableMatrix &dPlasticMicroGradientLdElasticPsi,
                                                          variableMatrix &dPlasticMicroGradientLdElasticGamma,
                                                          variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection ){
        /*!
         * Compute the plastic micro gradient velocity gradient.
         *
         * \bar{L}_{ \bar{N} \bar{M}, \bar{K} }^{\chi, p} = \bar{ \Psi }_{ \bar{N} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } }_{\bar{I} } \frac{ \partial \bar{ G }_{ \bar{I} }^{ \nabla \chi } }{ \partial \bar{ M }_{ \bar{K} \bar{L} \bar{M} } } + 2 \bar{ \Psi }_{ \bar{L} \bar{D} }^{e} \text{ skw } \left[ \bar{L}_{ \bar{D} \bar{C} }^{ \chi, p } \bar{ \Psi }_{ \bar{C} \bar{F} }^{e, -1} \Gamma_{ \bar{F} \bar{M} \bar{K} }^{e} \right]
         *
         * Note: The user must ensure that elasticPsi and inverseElasticPsi are inverses of each other. This is not checked in the code.
         * 
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
         * :param const variableVector &inverseElasticPsi: The inverse elastic micro deformation measure Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation measure Gamma.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction for the micro gradient plasticity.
         * :param const variableVector &plasticMicroVelocityGradient: The velocity gradient for micro plasticity.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient gamma.
         * :param variableMatrix &dPlasticMicroGradientLdPlasticMicroL: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the platic micro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientLdElasticPsi: The Jacobian of the plastic micro gradient
         *     velocity gradient w.r.t. the elastic micro deformation tensor Psi.
         * :param variableMatrix &dPlasticMicroGradientLdElasticGamma: The Jacobian of the plastic micro gradient
         *     velocity gradient w.r.t. the elastic higher ordrer deformation tensor Gamma.
         * :param variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection: The Jacobian of the plastic micro gradient
         *     velocity gradient w.r.t. the micro gradient flow direction.
         */

        //Assume 3D
        unsigned int dim = 3;

        variableVector skewTerm;
        errorOut error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                                      elasticGamma, microGradientFlowDirection,
                                                                      plasticMicroVelocityGradient, 
                                                                      plasticMicroGradientVelocityGradient, skewTerm,
                                                                      dPlasticMicroGradientLdMicroGradientGamma,
                                                                      dPlasticMicroGradientLdPlasticMicroL );

        if ( error ){
            errorOut result = new errorNode( "computePlasticMicroGradientVelocityGradient (jacobian)",
                                             "Error in computation of the plasticMicroVelocityGradient" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dPlasticMicroGradientLdElasticPsi = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dPlasticMicroGradientLdElasticGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dPlasticMicroGradientLdMicroGradientFlowDirection = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim * dim, 0 ) );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int Mb = 0; Mb < dim; Mb++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){
                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){
                            dPlasticMicroGradientLdElasticPsi[ dim * dim * Db + dim * Mb + Kb ][ dim * Ob + Pb ]
                                += -inverseElasticPsi[ dim * Db + Ob ]
                                 * plasticMicroGradientVelocityGradient[ dim * dim * Pb + dim * Mb + Kb ]
                                 + inverseElasticPsi[ dim * Db + Ob ] * skewTerm[ dim * dim * Pb + dim * Mb + Kb ];

                            for ( unsigned int Qb = 0; Qb < dim; Qb++ ){
                                dPlasticMicroGradientLdElasticGamma[ dim * dim * Db + dim * Mb + Kb ][ dim * dim * Ob + dim * Pb + Qb ]
                                    -= plasticMicroVelocityGradient[ dim * Pb + Mb ]
                                     * inverseElasticPsi[ dim * Db + Ob ] * eye[ dim * Kb + Qb ]; 

                                for ( unsigned int Rb = 0; Rb < dim; Rb++ ){
                                    dPlasticMicroGradientLdElasticPsi[ dim * dim * Db + dim * Mb + Kb ][ dim * Ob + Pb ]
                                        -= plasticMicroVelocityGradient[ dim * Db + Qb ]
                                         * inverseElasticPsi[ dim * Qb + Ob ] * inverseElasticPsi[ dim * Pb + Rb ]
                                         * elasticGamma[ dim * dim * Rb + dim * Mb + Kb ]
                                         - plasticMicroVelocityGradient[ dim * Qb + Mb ]
                                         * inverseElasticPsi[ dim * Db + Ob ] * inverseElasticPsi[ dim * Pb + Rb ]
                                         * elasticGamma[ dim * dim * Rb + dim * Qb + Kb ];
                                    
                                    dPlasticMicroGradientLdElasticGamma[ dim * dim * Db + dim * Mb + Kb ][ dim * dim * Ob + dim * Pb + Qb ]
                                        += plasticMicroVelocityGradient[ dim * Db + Rb ] * inverseElasticPsi[ dim * Rb + Ob ]
                                         * eye[ dim * Mb + Pb ] * eye[ dim * Kb + Qb ];

                                    dPlasticMicroGradientLdMicroGradientFlowDirection[ dim * dim * Db + dim * Mb + Kb ][ dim * dim * dim * Ob + dim * dim * Pb + dim * Qb + Rb ]
                                        += inverseElasticPsi[ dim * Db + Qb ] * microGradientGamma[ Ob ]
                                         * eye[ dim * Kb + Pb ] * eye[ dim * Mb + Rb ];
                                }
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma,
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableMatrix &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient ){
        /*!
         * Compute the plastic velocity gradients in the intermediate configuration.
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation metric Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation metric Gamma.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction of the micro gradient plasticity.
         *     Note: This is a matrix because it is computed as the gradient of the flow potential which is a vector.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &plasticMicroVelocityGradient: The plastic velocity gradient for the micro plastic deformation.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic velocity gradient for the micro gradient 
         *     plastic deformation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients",
                                  "The elastic micro deformation metric Psi must be 3D" );
        }

        //Compute the required inverses of the deformation metrics
        variableVector inverseElasticRightCauchyGreen = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = vectorTools::inverse( elasticPsi, dim, dim );

        //Compute the macro velocity gradient
        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro velocity gradient
        error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen, elasticPsi, inverseElasticPsi,
                                                     microFlowDirection, plasticMicroVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro gradient velocity gradient
        error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection, 
                                                             plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients",
                                             "Error in computation of plastic micro gradient velocity gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma,
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableMatrix &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient,
                                              variableVector &dPlasticMacroLdMacroGamma,
                                              variableVector &dPlasticMacroLdMicroGamma, variableVector &dPlasticMicroLdMicroGamma,
                                              variableVector &dPlasticMicroGradientLdMicroGamma,
                                              variableMatrix &dPlasticMicroGradientLdMicroGradientGamma ){
        /*!
         * Compute the plastic velocity gradients in the intermediate configuration.
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation metric Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation metric Gamma.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction of the micro gradient plasticity.
         *     Note: This is a matrix because it is computed as the gradient of the flow potential which is a vector.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &plasticMicroVelocityGradient: The plastic velocity gradient for the micro plastic deformation.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic velocity gradient for the micro gradient 
         *     plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGamma: The Jacobian of the plastic micro gradient velocity
         *     gradient w.r.t. the micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient plastic multiplier.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (jacobian)",
                                  "The elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (jacobian)",
                                  "The elastic micro deformation metric Psi must be 3D" );
        }

        //Compute the required inverses of the deformation metrics
        variableVector inverseElasticRightCauchyGreen = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = vectorTools::inverse( elasticPsi, dim, dim );

        //Compute the plastic macro velocity gradient
        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                              dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (jacobian)",
                                             "Error in computation of plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro velocity gradient
        error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen, elasticPsi, inverseElasticPsi,
                                                     microFlowDirection, plasticMicroVelocityGradient,
                                                     dPlasticMicroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (jacobian)",
                                             "Error in computation of plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro gradient velocity gradient
        variableMatrix dPlasticMicroGradientLdPlasticMicroL;
        error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection, 
                                                             plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                             dPlasticMicroGradientLdMicroGradientGamma,
                                                             dPlasticMicroGradientLdPlasticMicroL );

        dPlasticMicroGradientLdMicroGamma = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                              dPlasticMicroLdMicroGamma );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (jacobian)",
                                             "Error in computation of plastic micro gradient velocity gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computePlasticVelocityGradients( const variableType &macroGamma, const variableType &microGamma,
                                              const variableVector &microGradientGamma, const variableVector &elasticRightCauchyGreen,
                                              const variableVector &elasticMicroRightCauchyGreen, const variableVector &elasticPsi,
                                              const variableVector &elasticGamma, const variableVector &macroFlowDirection,
                                              const variableVector &microFlowDirection, const variableMatrix &microGradientFlowDirection,
                                              variableVector &plasticMacroVelocityGradient, variableVector &plasticMicroVelocityGradient,
                                              variableVector &plasticMicroGradientVelocityGradient,
                                              variableVector &dPlasticMacroLdMacroGamma,
                                              variableVector &dPlasticMacroLdMicroGamma, variableVector &dPlasticMicroLdMicroGamma,
                                              variableVector &dPlasticMicroGradientLdMicroGamma,
                                              variableMatrix &dPlasticMicroGradientLdMicroGradientGamma,
                                              variableMatrix &dPlasticMacroLdElasticRCG,
                                              variableMatrix &dPlasticMacroLdMacroFlowDirection,
                                              variableMatrix &dPlasticMacroLdMicroFlowDirection,
                                              variableMatrix &dPlasticMicroLdElasticMicroRCG,
                                              variableMatrix &dPlasticMicroLdElasticPsi,
                                              variableMatrix &dPlasticMicroLdMicroFlowDirection,
                                              variableMatrix &dPlasticMicroGradientLdElasticMicroRCG,
                                              variableMatrix &dPlasticMicroGradientLdElasticPsi,
                                              variableMatrix &dPlasticMicroGradientLdElasticGamma,
                                              variableMatrix &dPlasticMicroGradientLdMicroFlowDirection,
                                              variableMatrix &dPlasticMicroGradientLdMicroGradientFlowDirection ){
        /*!
         * Compute the plastic velocity gradients in the intermediate configuration.
         *
         * :param const variableType &macroGamma: The macro plastic multiplier.
         * :param const variableType &microGamma: The micro plastic multiplier.
         * :param const variableVector &microGradientGamma: The micro gradient plastic multiplier.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
         * :param const variableVector &elasticPsi: The elastic micro deformation metric Psi.
         * :param const variableVector &elasticGamma: The elastic higher order deformation metric Gamma.
         * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
         * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
         * :param const variableMatrix &microGradientFlowDirection: The flow direction of the micro gradient plasticity.
         *     Note: This is a matrix because it is computed as the gradient of the flow potential which is a vector.
         * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
         * :param variableVector &plasticMicroVelocityGradient: The plastic velocity gradient for the micro plastic deformation.
         * :param variableVector &plasticMicroGradientVelocityGradient: The plastic velocity gradient for the micro gradient 
         *     plastic deformation.
         * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     macro plastic multiplier.
         * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient w.r.t. the 
         *     micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGamma: The Jacobian of the plastic micro gradient velocity
         *     gradient w.r.t. the micro plastic multiplier.
         * :param variableVector &dPlasticMicroGradientLdMicroGradientGamma: The Jacobian of the plastic micro gradient 
         *     velocity gradient w.r.t. the micro gradient plastic multiplier.
         * :param variableVector &dPlasticMacroLdElasticRCG: The Jacobian of the plastic macro velocity gradient w.r.t. 
         *     the elastic right Cauchy-Green deformation tensor.
         * :param variableVector &dPlasticMacroLdMacroFlowDirection: The Jacobian of the plastic macro velocity gradient
         *     w.r.t. the macro flow direction.
         * :param variableVector &dPlasticMacroLdMicroFlowDirection: The Jacobian of the plastic macro velocity gradient
         *     w.r.t. the micro flow direction.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( elasticRightCauchyGreen.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                  "The elastic right Cauchy-Green deformation tensor must be 3D" );
        }

        if ( elasticPsi.size() != dim * dim ){
            return new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                  "The elastic micro deformation metric Psi must be 3D" );
        }

        //Compute the required inverses of the deformation metrics
        variableVector inverseElasticRightCauchyGreen = vectorTools::inverse( elasticRightCauchyGreen, dim, dim );
        variableVector inverseElasticPsi = vectorTools::inverse( elasticPsi, dim, dim );

        //Compute the plastic macro velocity gradient
        errorOut error = computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                              macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                              dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma,
                                                              dPlasticMacroLdElasticRCG, dPlasticMacroLdMacroFlowDirection,
                                                              dPlasticMacroLdMicroFlowDirection );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                             "Error in computation of plastic macro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro velocity gradient
        error = computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen, elasticPsi, inverseElasticPsi,
                                                     microFlowDirection, plasticMicroVelocityGradient,
                                                     dPlasticMicroLdMicroGamma, dPlasticMicroLdElasticMicroRCG,
                                                     dPlasticMicroLdElasticPsi, dPlasticMicroLdMicroFlowDirection );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                             "Error in computation of plastic micro velocity gradient" );
            result->addNext( error );
            return result;
        }

        //Compute the plastic micro gradient velocity gradient
        variableMatrix dPlasticMicroGradientLdPlasticMicroL;
        error = computePlasticMicroGradientVelocityGradient( microGradientGamma, elasticPsi, inverseElasticPsi,
                                                             elasticGamma, microGradientFlowDirection, 
                                                             plasticMicroVelocityGradient, plasticMicroGradientVelocityGradient,
                                                             dPlasticMicroGradientLdMicroGradientGamma,
                                                             dPlasticMicroGradientLdPlasticMicroL,
                                                             dPlasticMicroGradientLdElasticPsi,
                                                             dPlasticMicroGradientLdElasticGamma,
                                                             dPlasticMicroGradientLdMicroGradientFlowDirection );

        dPlasticMicroGradientLdMicroGamma = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                              dPlasticMicroLdMicroGamma );

        dPlasticMicroGradientLdElasticMicroRCG = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                                   dPlasticMicroLdElasticMicroRCG );

        dPlasticMicroGradientLdElasticPsi += vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                               dPlasticMicroLdElasticPsi );

        dPlasticMicroGradientLdMicroFlowDirection = vectorTools::dot( dPlasticMicroGradientLdPlasticMicroL,
                                                                      dPlasticMicroLdMicroFlowDirection );

        if ( error ){
            errorOut result = new errorNode( "computePlasticVelocityGradients (full jacobian)",
                                             "Error in computation of plastic micro gradient velocity gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        const parameterType alpha ){
        /*!
         * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
         *
         * :param const variableType &Dt: The change in time.
         * :param const variableVector &currentPlasticMicroDeformation: The inverse of the current micro deformation.
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
         *     velocity gradient.
         * :param const variableVector &previousPlasticMicroDeformation: The plastic micro deformation 
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradient: The micro gradient deformation in the 
         *     intermediate configuation from the last converged increment.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
         *     velocity gradient from the last converged increment.
         * :param variableVector &currentPlasticMicroGradient: The current plastic micro gradient 
         *    deformation in the intermediate configuration.
         * :param parameterType alpha: The integration parameter.
         */

        variableMatrix LHS;
        return evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                          currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                          previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                          previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                          previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient, LHS,
                                          alpha );
    }

    errorOut evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableMatrix &LHS,
                                        const parameterType alpha ){
        /*!
         * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
         *
         * :param const variableType &Dt: The change in time.
         * :param const variableVector &currentPlasticMicroDeformation: The inverse of the current micro deformation.
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
         *     velocity gradient.
         * :param const variableVector &previousPlasticMicroDeformation: The the plastic micro deformation 
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradient: The micro gradient deformation in the 
         *     intermediate configuation from the last converged increment.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
         *     velocity gradient from the last converged increment.
         * :param variableVector &currentPlasticMicroGradient: The current plastic micro gradient 
         *    deformation in the intermediate configuration.
         * :param variableMatrix &LHS: The left-hand-side matrix.
         * :param parameterType alpha: The integration parameter.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( currentPlasticMicroDeformation.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic micro-deformation must be 3D" );
        }

        if ( currentPlasticMacroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic macro velocity gradient must be 3D" );
        }

        if ( currentPlasticMicroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic micro velocity gradient must be 3D" );
        }

        if ( currentPlasticMicroGradientVelocityGradient.size() != dim * dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The plastic micro gradient velocity gradient must be 3D" );
        }

        if ( previousPlasticMicroDeformation.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro-deformation must be 3D" );
        }

        if ( previousPlasticMicroGradient.size() != dim * dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro gradient must be 3D" );
        }

        if ( previousPlasticMacroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic macro velocity gradient must be 3D" );
        }

        if ( previousPlasticMicroVelocityGradient.size() != dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro velocity gradient must be 3D" );
        }

        if ( previousPlasticMicroGradientVelocityGradient.size() != dim * dim * dim ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The previous plastic micro gradient velocity gradient must be 3D" );
        }

        //Compute the required identity terms
        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        //Assemble the A term ( forcing term ) and the fourth order A term
        variableVector DtAtilde( dim * dim * dim, 0 );
        variableVector previousFourthA( dim * dim * dim * dim, 0 );
        variableVector currentFourthA( dim * dim * dim * dim, 0 );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int B = 0; B < dim; B++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                        DtAtilde[ dim * dim * Db + dim * B + Kb ] += Dt 
                            * ( alpha * previousPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Lb + Kb ]
                                      *  previousPlasticMicroDeformation[ dim * Lb + B ]
                            + ( 1. - alpha ) * currentPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Lb + Kb ]
                                             * currentPlasticMicroDeformation[ dim * Lb + B ] );

                        previousFourthA[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Lb ]
                            = ( previousPlasticMicroVelocityGradient[ dim * Db + B ] * eye[ dim * Kb + Lb ]
                            -   previousPlasticMacroVelocityGradient[ dim * Lb + Kb ] * eye[ dim * Db + B ] );

                        currentFourthA[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Lb ]
                            = ( currentPlasticMicroVelocityGradient[ dim * Db + B ] * eye[ dim * Kb + Lb ]
                            -   currentPlasticMacroVelocityGradient[ dim * Lb + Kb ] * eye[ dim * Db + B ] );
                    }
                }
            }
        }

        //Assemble the right-hand side and left-hand side term
        variableVector RHS = DtAtilde;
        LHS = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int B = 0; B < dim; B++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){
                       for ( unsigned int Bb = 0; Bb < dim; Bb++ ){
                          RHS[ dim * dim * Db + dim * B + Kb ]
                             += ( eye[ dim * Db + Bb ] * eye[ dim * Kb + Lb ] + Dt * alpha * previousFourthA[ dim * dim * dim * Db + dim * dim * Bb + dim * Kb + Lb ] )
                              * previousPlasticMicroGradient[ dim * dim * Bb + dim * B + Lb ];
                          for ( unsigned int Sb = 0; Sb < dim; Sb++ ){
                              LHS[ dim * dim * Db + dim * B + Kb ][ dim * dim * Lb + dim * Bb + Sb ]
                                  = ( eye[ dim * Db + Lb ] * eye[ dim * Kb + Sb ] - Dt * ( 1. - alpha ) * currentFourthA[ dim * dim * dim * Db + dim * dim * Lb + dim * Kb + Sb ] ) * eye[ dim * B + Bb ];
                          }
                       } 
                    }
                }
            }
        }

        //Solve for the current plastic micro gradient
        unsigned int rank;
        currentPlasticMicroGradient = vectorTools::solveLinearSystem( LHS, RHS, rank );

        if ( rank != LHS.size() ){
            return new errorNode( "evolvePlasticMicroGradChi",
                                  "The left hand side matrix is not full rank" );
        }

        return NULL;
    }

    errorOut evolvePlasticMicroGradChi( const variableType &Dt,
                                        const variableVector &currentPlasticMicroDeformation,
                                        const variableVector &currentPlasticMacroVelocityGradient,
                                        const variableVector &currentPlasticMicroVelocityGradient,
                                        const variableVector &currentPlasticMicroGradientVelocityGradient,
                                        const variableVector &previousPlasticMicroDeformation,
                                        const variableVector &previousPlasticMicroGradient,
                                        const variableVector &previousPlasticMacroVelocityGradient,
                                        const variableVector &previousPlasticMicroVelocityGradient,
                                        const variableVector &previousPlasticMicroGradientVelocityGradient,
                                        variableVector &currentPlasticMicroGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                        variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                        const parameterType alpha ){
        /*!
         * Evolve the plastic micro gradient of the micro-deformation measure in the intermediate configuration.
         *
         * :param const variableType &Dt: The change in time.
         * :param const variableVector &currentPlasticMicroDeformation: The inverse of the current micro deformation.
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient
         *     velocity gradient.
         * :param const variableVector &previousPlasticMicroDeformation: The the plastic micro deformation 
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradient: The micro gradient deformation in the 
         *     intermediate configuation from the last converged increment.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient
         *     from the last converged increment.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient 
         *     velocity gradient from the last converged increment.
         * :param variableVector &currentPlasticMicroGradient: The current plastic micro gradient 
         *    deformation in the intermediate configuration.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroDeformation: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic micro deformation.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic macro velocity gradient.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic micro velocity gradient.
         * :param variableMatrix &dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient: The jacobian of the plastic 
         *     micro deformation w.r.t. the plastic micro gradient velocity gradient.
         * :param parameterType alpha: The integration parameter.
         */

        //Assume 3D
        unsigned int dim = 3;

        //Compute the required identity terms
        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        //Compute the new currentPlasticMicroGradient
        variableMatrix LHS;
        errorOut error = evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                                    currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                                    previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                                    previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                                    previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                                    LHS, alpha );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticMicroGradChi",
                                             "Error in the evolution of the plastic micro gradient of chi" );
            result->addNext( error );
            return result;
        }

        //Compute the negative partial derivatives w.r.t. currentFourthA and the current part of DtAtilde
        //We do this in vector form so that we can interface with Eigen easier
        variableVector negdRdCurrentDtAtilde( dim * dim * dim * dim * dim * dim, 0 );
        variableVector negdRdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim, 0 );

        //Also assemble jacobians of the A terms
        variableMatrix dCurrentDTAtildedPlasticMicroDeformation( dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dCurrentDTAtildedPlasticMicroGradientVelocityGradient( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        variableMatrix dCurrentFourthAdMacroVelocityGradient( dim * dim * dim * dim, variableVector( dim * dim, 0 ) );
        variableMatrix dCurrentFourthAdMicroVelocityGradient( dim * dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int Db = 0; Db < dim; Db++ ){
            for ( unsigned int B = 0; B < dim; B++ ){
                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){
                    for ( unsigned int Rb = 0; Rb < dim; Rb++ ){
                        for ( unsigned int S = 0; S < dim; S++ ){
                            dCurrentDTAtildedPlasticMicroDeformation[ dim * dim * Db + dim * B + Kb ][ dim * Rb + S ]
                                += Dt * ( 1. - alpha ) * currentPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * Rb + Kb ]
                                 * eye[ dim * B + S ];

                            for ( unsigned int Tb = 0; Tb < dim; Tb++ ){
                                negdRdCurrentDtAtilde[ dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * B + dim * dim * dim * Kb + dim * dim * Rb + dim * S + Tb ]
                                    += eye[ dim * Db + Rb ] * eye[ dim * B + S ] * eye[ dim * Kb + Tb ];

                                dCurrentDTAtildedPlasticMicroGradientVelocityGradient[ dim * dim * Db + dim * B + Kb ][ dim * dim * Rb + dim * S + Tb ]
                                    += Dt * ( 1. - alpha ) * eye[ dim * Db + Rb ] * eye[ dim * Kb + Tb ] * currentPlasticMicroDeformation[ dim * S + B ];

                                dCurrentFourthAdMacroVelocityGradient[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Rb ][ dim * S + Tb ]
                                    -= eye[ dim * Rb + S ] * eye[ dim * Kb + Tb ] * eye[ dim * Db + B ];

                                dCurrentFourthAdMicroVelocityGradient[ dim * dim * dim * Db + dim * dim * B + dim * Kb + Rb ][ dim * S + Tb ]
                                    += eye[ dim * Db + S ] * eye[ dim * B + Tb ] * eye[ dim * Kb + Rb ];

                                for ( unsigned int Ub = 0; Ub < dim; Ub++ ){
                                    negdRdCurrentFourthA[ dim * dim * dim * dim * dim * dim * Db + dim * dim * dim * dim * dim * B + dim * dim * dim * dim * Kb + dim * dim * dim * Rb + dim * dim * S + dim * Tb + Ub ]
                                        += Dt * ( 1. - alpha ) * eye[ dim * Db + Rb ] * eye[ dim * Kb + Tb ]
                                         * currentPlasticMicroGradient[ dim * dim * S + dim * B + Ub ];

                                }
                            }
                        }
                    }
                }
            }
        }

        //Solve for the Jacobians
        variableVector vecdCurrentPlasticMicroGradientdCurrentDTAtilde( dim * dim * dim * dim * dim * dim );
        variableVector vecdCurrentPlasticMicroGradientdCurrentFourthA( dim * dim * dim * dim * dim * dim * dim );

        variableVector floatLHS = vectorTools::appendVectors( LHS );

        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > LHSMat( floatLHS.data(), LHS.size(), LHS.size() );
        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCDA( negdRdCurrentDtAtilde.data(), LHS.size(), dim * dim * dim );
        Eigen::Map< const Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > nDRDCFA( negdRdCurrentFourthA.data(), LHS.size(), dim * dim * dim * dim );

        Eigen::ColPivHouseholderQR< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > qrSolver( LHSMat );

        Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X1( vecdCurrentPlasticMicroGradientdCurrentDTAtilde.data(), LHS.size(), dim * dim * dim );
        Eigen::Map< Eigen::Matrix< variableType, -1, -1, Eigen::RowMajor > > X2( vecdCurrentPlasticMicroGradientdCurrentFourthA.data(), LHS.size(), dim * dim * dim * dim );

        X1 = qrSolver.solve( nDRDCDA );
        X2 = qrSolver.solve( nDRDCFA );

        variableMatrix dCurrentPlasticMicroGradientdCurrentDTAtilde = vectorTools::inflate( vecdCurrentPlasticMicroGradientdCurrentDTAtilde, dim * dim * dim, dim * dim * dim );
        variableMatrix dCurrentPlasticMicroGradientdCurrentFourthA = vectorTools::inflate( vecdCurrentPlasticMicroGradientdCurrentFourthA, dim * dim * dim, dim * dim * dim * dim );

        //Assemble the final terms of the deformation
        dCurrentPlasticMicroGradientdPlasticMicroDeformation = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                 dCurrentDTAtildedPlasticMicroDeformation );

        dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                      dCurrentFourthAdMacroVelocityGradient );

        dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentFourthA,
                                                                                      dCurrentFourthAdMicroVelocityGradient );

        dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient = vectorTools::dot( dCurrentPlasticMicroGradientdCurrentDTAtilde,
                                                                                 dCurrentDTAtildedPlasticMicroGradientVelocityGradient );
        return NULL;
    }

    errorOut evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       const parameterType alphaMacro,
                                       const parameterType alphaMicro,
                                       const parameterType alphaMicroGradient ){
        /*!
         * Evolve the plastic deformation
         *
         * :param const variableType &Dt: The timestep
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient 
         *     velocity gradient.
         * :param const variableVector &previousPlasticDeformationGradient: The plastic deformation gradient at the end of the last 
         *     converged timestep.
         * :param const variableVector &previousPlasticMicroDeformation: The plastic micro deformation at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMicroGradient: The plastic micro gradient at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient 
         *     at the end of the last converged timestep.
         * :param variableVector &currentPlasticDeformationGradient: The current value of the plastic deformation gradient.
         * :param variableVector &currentPlasticMicroDeformation: The current value of the plastic micro deformation.
         * :param variableVector &currentPlasticMicroGradient: The current value of the plastic micro gradient.
         * :param parameterType alphaMacro: The integration parameter for the macro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicro: The integration parameter for the micro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicroGradient: The integration parameter for the micro gradient plasticity. Defaults to 0.5.
         */

        errorOut error = constitutiveTools::evolveF( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                                     currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                                     alphaMacro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic macro deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            alphaMicro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic micro deformation" );
            result->addNext( error );
            return result;
        }

        error = evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           alphaMicroGradient );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic micro gradient" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut evolvePlasticDeformation( const variableType &Dt,
                                       const variableVector &currentPlasticMacroVelocityGradient,
                                       const variableVector &currentPlasticMicroVelocityGradient,
                                       const variableVector &currentPlasticMicroGradientVelocityGradient,
                                       const variableVector &previousPlasticDeformationGradient,
                                       const variableVector &previousPlasticMicroDeformation,
                                       const variableVector &previousPlasticMicroGradient,
                                       const variableVector &previousPlasticMacroVelocityGradient,
                                       const variableVector &previousPlasticMicroVelocityGradient,
                                       const variableVector &previousPlasticMicroGradientVelocityGradient,
                                       variableVector &currentPlasticDeformationGradient,
                                       variableVector &currentPlasticMicroDeformation,
                                       variableVector &currentPlasticMicroGradient,
                                       variableMatrix &dPlasticFdPlasticMacroL,
                                       variableMatrix &dPlasticMicroDeformationdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMacroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroL,
                                       variableMatrix &dPlasticMicroGradientdPlasticMicroGradientL,
                                       const parameterType alphaMacro,
                                       const parameterType alphaMicro,
                                       const parameterType alphaMicroGradient ){
        /*!
         * Evolve the plastic deformation
         *
         * :param const variableType &Dt: The timestep
         * :param const variableVector &currentPlasticMacroVelocityGradient: The current plastic macro velocity gradient.
         * :param const variableVector &currentPlasticMicroVelocityGradient: The current plastic micro velocity gradient.
         * :param const variableVector &currentPlasticMicroGradientVelocityGradient: The current plastic micro gradient 
         *     velocity gradient.
         * :param const variableVector &previousPlasticDeformationGradient: The plastic deformation gradient at the end of the last 
         *     converged timestep.
         * :param const variableVector &previousPlasticMicroDeformation: The plastic micro deformation at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMicroGradient: The plastic micro gradient at the end of the last converged 
         *     timestep.
         * :param const variableVector &previousPlasticMacroVelocityGradient: The plastic macro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroVelocityGradient: The plastic micro velocity gradient at the end of the 
         *     last converged timestep.
         * :param const variableVector &previousPlasticMicroGradientVelocityGradient: The plastic micro gradient velocity gradient 
         *     at the end of the last converged timestep.
         * :param variableVector &currentPlasticDeformationGradient: The current value of the plastic deformation gradient.
         * :param variableVector &currentPlasticMicroDeformation: The current value of the plastic micro deformation.
         * :param variableVector &currentPlasticMicroGradient: The current value of the plastic micro gradient.
         * :param variableMatrix &dPlasticFdPlasticMacroL: The Jacobian of the plastic deformation gradient w.r.t. the plastic 
         *     macro velocity gradient.
         * :param variableMatrix &dPlasticMicroDeformationdPlasticMicroL: The Jacobian of the plastic micro-deformation w.r.t. 
         *     the plastic micro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientdPlasticMacroL: The Jacobian of the plastic micro gradient deformation 
         *     w.r.t. the plastic macro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientdPlasticMicroL: The Jacobian of the plastic micro gradient deformation
         *     w.r.t. the plastic micro velocity gradient.
         * :param variableMatrix &dPlasticMicroGradientdPlasticMicroGradientL: The Jacobian of the plastic micro gradient deformation
         *     w.r.t. the plastic micro gradient velocity gradient.
         * :param parameterType alphaMacro: The integration parameter for the macro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicro: The integration parameter for the micro plasticity. Defaults to 0.5.
         * :param parameterType alphaMicroGradient: The integration parameter for the micro gradient plasticity. Defaults to 0.5.
         */

        errorOut error = constitutiveTools::evolveF( Dt, previousPlasticDeformationGradient, previousPlasticMacroVelocityGradient,
                                                     currentPlasticMacroVelocityGradient, currentPlasticDeformationGradient,
                                                     dPlasticFdPlasticMacroL, alphaMacro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation (jacobian)",
                                             "Error in computation of the plastic macro deformation gradient" );
            result->addNext( error );
            return result;
        }

        error = constitutiveTools::evolveF( Dt, previousPlasticMicroDeformation, previousPlasticMicroVelocityGradient,
                                            currentPlasticMicroVelocityGradient, currentPlasticMicroDeformation,
                                            dPlasticMicroDeformationdPlasticMicroL, alphaMicro, 1 );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation (jacobian)",
                                             "Error in computation of the plastic micro deformation" );
            result->addNext( error );
            return result;
        }

        variableMatrix dPlasticMicroGradientdPlasticMicroDeformation;
        error = evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation, currentPlasticMacroVelocityGradient,
                                           currentPlasticMicroVelocityGradient, currentPlasticMicroGradientVelocityGradient,
                                           previousPlasticMicroDeformation, previousPlasticMicroGradient,
                                           previousPlasticMacroVelocityGradient, previousPlasticMicroVelocityGradient,
                                           previousPlasticMicroGradientVelocityGradient, currentPlasticMicroGradient,
                                           dPlasticMicroGradientdPlasticMicroDeformation,
                                           dPlasticMicroGradientdPlasticMacroL, dPlasticMicroGradientdPlasticMicroL,
                                           dPlasticMicroGradientdPlasticMicroGradientL, alphaMicroGradient );

        if ( error ){
            errorOut result = new errorNode( "evolvePlasticDeformation",
                                             "Error in computation of the plastic micro gradient" );
            result->addNext( error );
            return result;
        }

        dPlasticMicroGradientdPlasticMicroL += vectorTools::dot( dPlasticMicroGradientdPlasticMicroDeformation,
                                                                 dPlasticMicroDeformationdPlasticMicroL );

        return NULL;
    }

    errorOut evolveStrainStateVariables( const constantType &Dt, const variableType &currentMacroGamma,
                                         const variableType &currentMicroGamma, const variableVector &currentMicroGradientGamma,
                                         const variableType &currentdMacroGdMacroC, const variableType &currentdMicroGdMicroC,
                                         const variableMatrix &currentdMicroGradientGdMicroGradientC,
                                         const variableType &previousMacroStrainISV, const variableType &previousMicroStrainISV,
                                         const variableVector &previousMicroGradientStrainISV,
                                         const variableType &previousMacroGamma, const variableType &previousMicroGamma,
                                         const variableVector &previousMicroGradientGamma, const variableType &previousdMacroGdMacroC,
                                         const variableType &previousdMicroGdMicroC,
                                         const variableMatrix &previousdMicroGradientGdMicroGradientC,
                                         variableType &currentMacroStrainISV, variableType &currentMicroStrainISV,
                                         variableVector &currentMicroGradientStrainISV,
                                         const parameterType alphaMacro,
                                         const parameterType alphaMicro,
                                         const parameterType alphaMicroGradient ){
        /*!
         * Evolve the strain-like state variables.
         *
         * :param const constantType &Dt: The timestep
         * :param const variableType &currentMacroGamma: The current macro-scale plastic multiplier.
         * :param const variableType &currentMicroGamma: The current micro-scale plastic multiplier.
         * :param const variableVector &currentMicroGradientGamma: The current micro-scale gradient plastic multiplier.
         * :param const variableType &currentdMacroGdMacroC: The current Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &currentdMicroGdMicroC: The current Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &currentdMicroGradientGdMicroGradientC: The current Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &previousMacroStrainISV: The previous value of the macro strain like ISV.
         * :param variableType &previousMicroStrainISV: The previous value of the micro strain like ISV.
         * :param variableVector &previousMicroGradientStrainISV: The previous value of the micro gradient strain like ISV.
         * :param const variableType &previousMacroGamma: The previous macro-scale plastic multiplier.
         * :param const variableType &previousMicroGamma: The previous micro-scale plastic multiplier.
         * :param const variableVector &previousMicroGradientGamma: The previous micro-scale gradient plastic multiplier.
         * :param const variableType &previousdMacroGdMacroC: The previous Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &previousdMicroGdMicroC: The previous Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &previousdMicroGradientGdMicroGradientC: The previous Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &currentMacroStrainISV: The new value of the macro strain ISV
         * :param variableType &currentMicroStrainISV: The new value of the micro strain ISV
         * :param variableVector &currentMicroGradientStrainISV: The new value of the micro gradient strain ISV.
         * :param const variableType alphaMacro: The integration parameter for the macro strain-like ISV
         * :param const variableType alphaMicro: The integration parameter for the micro strain-like ISV
         * :param const variableType alphaMicroGradient: The integration parameter for the micro Gradient strain-like ISV
         */

        if ( currentMicroGradientGamma.size() != 3 ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The current micro gradient gamma must be 3 dimensional" );
        }

        if ( previousMicroGradientGamma.size() != 3 ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The previous micro gradient gamma must be 3 dimensional" );
        }

        if ( currentMicroGradientGamma.size() != currentdMicroGradientGdMicroGradientC.size() ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The current micro gradient gamma and the current derivative of the plastic potential function w.r.t. the micro gradient cohesion are not consistent" );
        }

        if ( previousMicroGradientGamma.size() != previousdMicroGradientGdMicroGradientC.size() ){
            return new errorNode( "evolveStrainStateVariables",
                                  "The previous micro gradient gamma and the previous derivative of the plastic potential function w.r.t. the micro gradient cohesion are not consistent" );
        }

        for ( unsigned int i = 0; i < currentdMicroGradientGdMicroGradientC.size(); i++ ){
            if ( currentdMicroGradientGdMicroGradientC[ i ].size() != previousMicroGradientStrainISV.size() ){
                return new errorNode( "evolveStrainStateVariables",
                                      "The current derivative of the plastic potential function w.r.t. the micro gradient cohesion is not a square matrix" ); 
            }
        }

        for ( unsigned int i = 0; i < previousdMicroGradientGdMicroGradientC.size(); i++ ){
            if ( previousdMicroGradientGdMicroGradientC[ i ].size() != previousMicroGradientStrainISV.size() ){
                return new errorNode( "evolveStrainStateVariables",
                                      "The previous derivative of the plastic potential function w.r.t. the micro gradient cohesion is not a square matrix" ); 
            }
        }

        //Evolve the macro-scale internal strain-like state variable
        currentMacroStrainISV = previousMacroStrainISV + Dt * (         - alphaMacro * previousMacroGamma * previousdMacroGdMacroC
                                                                - ( 1 - alphaMacro ) *  currentMacroGamma * currentdMacroGdMacroC );

        //Evolve the micro-scale internal strain-like state variable
        currentMicroStrainISV = previousMicroStrainISV + Dt * (         - alphaMicro * previousMicroGamma * previousdMicroGdMicroC
                                                                - ( 1 - alphaMicro ) *  currentMicroGamma * currentdMicroGdMicroC );

        //Evolve the micro gradient internal strain-like state variables
        currentMicroGradientStrainISV = previousMicroGradientStrainISV
        + Dt * ( 
            - alphaMicroGradient * vectorTools::Tdot( previousdMicroGradientGdMicroGradientC, previousMicroGradientGamma )
            - ( 1 - alphaMicroGradient ) * vectorTools::Tdot( currentdMicroGradientGdMicroGradientC, currentMicroGradientGamma )
        );

        return NULL;
    }

    errorOut evolveStrainStateVariables( const constantType &Dt, const variableType &currentMacroGamma,
                                         const variableType &currentMicroGamma, const variableVector &currentMicroGradientGamma,
                                         const variableType &currentdMacroGdMacroC, const variableType &currentdMicroGdMicroC,
                                         const variableMatrix &currentdMicroGradientGdMicroGradientC,
                                         const variableType &previousMacroStrainISV, const variableType &previousMicroStrainISV,
                                         const variableVector &previousMicroGradientStrainISV,
                                         const variableType &previousMacroGamma, const variableType &previousMicroGamma,
                                         const variableVector &previousMicroGradientGamma, const variableType &previousdMacroGdMacroC,
                                         const variableType &previousdMicroGdMicroC,
                                         const variableMatrix &previousdMicroGradientGdMicroGradientC,
                                         variableType &currentMacroStrainISV, variableType &currentMicroStrainISV,
                                         variableVector &currentMicroGradientStrainISV,
                                         variableType &dCurrentMacroISVdCurrentMacroGamma, variableType &dCurrentMacroISVddMacroGdMacroC,
                                         variableType &dCurrentMicroISVdCurrentMicroGamma, variableType &dCurrentMicroISVddMicroGdMicroC,
                                         variableMatrix &dCurrentMicroGradISVdCurrentMicroGradGamma,
                                         variableMatrix &dCurrentMicroGradISVddMicroGradGdMicroGradC,
                                         const parameterType alphaMacro,
                                         const parameterType alphaMicro,
                                         const parameterType alphaMicroGradient){
        /*!
         * Evolve the strain-like state variables.
         *
         * :param const constantType &Dt: The timestep
         * :param const variableType &currentMacroGamma: The current macro-scale plastic multiplier.
         * :param const variableType &currentMicroGamma: The current micro-scale plastic multiplier.
         * :param const variableVector &currentMicroGradientGamma: The current micro-scale gradient plastic multiplier.
         * :param const variableType &currentdMacroGdMacroC: The current Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &currentdMicroGdMicroC: The current Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &currentdMicroGradientGdMicroGradientC: The current Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &previousMacroStrainISV: The previous value of the macro strain like ISV.
         * :param variableType &previousMicroStrainISV: The previous value of the micro strain like ISV.
         * :param variableVector &previousMicroGradientStrainISV: The previous value of the micro gradient strain like ISV.
         * :param const variableType &previousMacroGamma: The previous macro-scale plastic multiplier.
         * :param const variableType &previousMicroGamma: The previous micro-scale plastic multiplier.
         * :param const variableVector &previousMicroGradientGamma: The previous micro-scale gradient plastic multiplier.
         * :param const variableType &previousdMacroGdMacroC: The previous Jacobian of the macro plastic potential 
         *     w.r.t. the macro cohesion.
         * :param const variableType &previousdMicroGdMicroC: The previous Jacobian of the micro plastic potential 
         *     w.r.t. the micro cohesion.
         * :param const variableVector &previousdMicroGradientGdMicroGradientC: The previous Jacobian of the micro gradient 
         *     plastic potential w.r.t. the micro gradient cohesion.
         * :param variableType &dCurrentMacroISVdCurrentMacroGamma: The Jacobian of the 
         * :param variableType &currentMacroStrainISV: The new value of the macro strain ISV
         * :param variableType &currentMicroStrainISV: The new value of the micro strain ISV
         * :param variableVector &currentMicroGradientStrainISV: The new value of the micro gradient strain ISV.
         * :param variableType &dCurrentMacroISVdCurrentMacroGamma: The Jacobian of the current macro strain ISV w.r.t. the 
         *     macro plastic multiplier.
         * :param variableType &dCurrentMacroISVddMacroGdMacroC: The Jacobian of the current macro strain ISV w.r.t. the 
         *     derivative of the macro plastic potential w.r.t. the macro cohesion.
         * :param variableType &dCurrentMicroISVdCurrentMicroGamma: The Jacobian of the current micro strain ISV w.r.t. the 
         *     micro plastic multiplier.
         * :param variableType &dCurrentMicroISVddMicroGdMicroC: The Jacobian of the current micro strain ISV w.r.t. the 
         *     derivative of the micro plastic potential w.r.t. the micro cohesion.
         * :param variableMatrix &dCurrentMicroGradientISVdCurrentMicroGradientGamma: The Jacobian of the current micro 
         *     strain ISV w.r.t. the micro plastic multiplier.
         * :param variableMatrix &dCurrentMicroGradientISVddMicroGradientGdMicroGradientC: The Jacobian of the current 
         *     micro gradient strain ISV w.r.t. the derivative of the micro gradient plastic potential w.r.t. the micro 
         *     gradient cohesion.
         * :param const variableType alphaMacro: The integration parameter for the macro strain-like ISV
         * :param const variableType alphaMicro: The integration parameter for the micro strain-like ISV
         * :param const variableType alphaMicroGradient: The integration parameter for the micro Gradient strain-like ISV
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = evolveStrainStateVariables( Dt, currentMacroGamma, currentMicroGamma, currentMicroGradientGamma,
                                                     currentdMacroGdMacroC, currentdMicroGdMicroC,
                                                     currentdMicroGradientGdMicroGradientC,
                                                     previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
                                                     previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
                                                     previousdMacroGdMacroC, previousdMicroGdMicroC,
                                                     previousdMicroGradientGdMicroGradientC,
                                                     currentMacroStrainISV, currentMicroStrainISV, currentMicroGradientStrainISV,
                                                     alphaMacro, alphaMicro, alphaMicroGradient );

        if ( error ){
            errorOut result = new errorNode( "evolveStrainStateVariables (jacobian)",
                                             "Error in the evolution of the strain-like ISVs" );
            result->addNext( error );
            return result;
        }

        //Compute the Jacobians of the macro ISV evolution w.r.t. the current macro gamma and the derivative of the plastic 
        //potential function w.r.t. the macro cohesion
        dCurrentMacroISVdCurrentMacroGamma = - Dt * ( 1 - alphaMacro ) * currentdMacroGdMacroC;
        dCurrentMacroISVddMacroGdMacroC = - Dt * ( 1 - alphaMacro ) * currentMacroGamma;

        //Compute the Jacobians of the micro ISV evolution w.r.t. the current micro gamma and the derivative of the plastic 
        //potential function w.r.t. the micro cohesion
        dCurrentMicroISVdCurrentMicroGamma = - Dt * ( 1 - alphaMicro ) * currentdMicroGdMicroC;
        dCurrentMicroISVddMicroGdMicroC = - Dt * ( 1 - alphaMicro ) * currentMicroGamma;

        //Compute the Jacobians of the micro gradient ISV evolution w.r.t. the current micro gradient gamma and the derivative
        //of the plastic potential function w.r.t. the micro gradient cohesion.
        dCurrentMicroGradISVdCurrentMicroGradGamma = variableMatrix( dim, variableVector( dim, 0 ) );
        dCurrentMicroGradISVddMicroGradGdMicroGradC = variableMatrix( dim, variableVector( dim * dim, 0 ) );

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                dCurrentMicroGradISVdCurrentMicroGradGamma[ i ][ j ] -= Dt * ( 1 - alphaMicroGradient )
                                                                      * currentdMicroGradientGdMicroGradientC[ j ][ i ];
                for ( unsigned int k = 0; k < dim; k++ ){
                    dCurrentMicroGradISVddMicroGradGdMicroGradC[ i ][ dim * j + k ]
                        -= Dt * ( 1 - alphaMicroGradient ) * currentMicroGradientGamma[ j ] * eye[ dim * i + k ];
                }
            }
        }

        return NULL;
    }

    errorOut computeFlowDirections( const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                    const variableVector &referenceHigherOrderStress, const variableType &macroCohesion,
                                    const variableType &microCohesion, const variableVector &microGradientCohesion,
                                    const variableVector &elasticRightCauchyGreen, const parameterVector &macroParameters,
                                    const parameterVector &microParameters, const parameterVector &microGradientParameters,
                                    variableVector &macroFlowDirection, variableVector &microFlowDirection,
                                    variableMatrix &microGradientFlowDirection, variableType &dGdMacroCohesion, 
                                    variableType &dGdMicroCohesion, variableMatrix &dGdMicroGradientCohesion ){
        /*!
         * Compute all of the flow directions.
         *
         * :param const variableVector &PK2Stress: The second Piola-Kirchhoff stress.
         * :param const variableVector &referenceMicroStress: The reference symmetric micro stress.
         * :param const variableVector &referenceHigherOrderStress: The reference higher order stress.
         * :param const variableType &macroCohesion: The macro cohesion value.
         * :param const variableType &microCohesion: The micro cohesion value.
         * :param const variableVector &microGradientCohesion: The micro gradient cohesion value.
         * :param const variableVector &elasticRightCauchyGreen: The elastic right cauchy green deformation.
         * :param const parameterVector &macroParameters: The macro plastic parameters.
         *     [ friction angle, beta, yield value, hardening curve values ]
         * :param const parameterVector &microParameters: The micro plastic parameters.
         *     [ friction angle, beta, yield value, hardening curve values ]
         * :param const parameterVector &microGradientParameters: The micro gradient plastic parameters.
         *     [ friction angle, beta, yield value, hardening curve values ]
         * :param variableVector &macroFlowDirection: The flow direction for the macro scale plasticity.
         * :param variableVector &microFlowDirection: The flow direction for the micro scale plasticity.
         * :param variableMatrix &microGradientFlowDirection: The flow direction for the micro gradient 
         *     plasticity.
         * :param variableType &dGdMacroCohesion: The Jacobian of the macro plastic potential w.r.t. the 
         *     macro cohesion.
         * :param variableType &dGdMicroCohesion: The Jacobian of the micro plastic potential w.r.t. the 
         *     micro cohesion.
         * :param variableVector &dGdMicroGradientCohesion: The Jacobian of the micro gradient plastic 
         *     potential w.r.t. the micro gradient cohesion.
         */

        variableVector tmp;

        parameterType macroFrictionAngle = macroParameters[ 0 ];
        parameterType macroBeta          = macroParameters[ 1 ];
        parameterType macroYieldValue    = macroParameters[ 2 ];
        errorOut error = computeSecondOrderDruckerPragerYieldEquation( PK2Stress, macroCohesion, elasticRightCauchyGreen,
                                                                       macroFrictionAngle, macroBeta, macroYieldValue, 
                                                                       macroFlowDirection, dGdMacroCohesion, tmpVec );

        if ( error ){
            errorOut result = new errorNode( "computeFlowDirections",
                                             "Error in the computation of the macro flow direction" );
            result->addNext( error );
            return result;
        }

        parameterType microFrictionAngle = microParameters[ 0 ];
        parameterType microBeta          = microParameters[ 1 ];
        parameterType microYieldValue    = microParameters[ 2 ];
        error = computeSecondOrderDruckerPragerYieldEquation( referenceMicroStress, microCohesion, elasticRightCauchyGreen,
                                                              microFrictionAngle, microBeta, microYieldValue, 
                                                              microFlowDirection, dGdMicroCohesion, tmpVec );

        if ( error ){
            errorOut result = new errorNode( "computeFlowDirections",
                                             "Error in the computation of the macro flow direction" );
            result->addNext( error );
            return result;
        }

        parameterType microGradientFrictionAngle = microGradientParameters[ 0 ];
        parameterType microGradientBeta          = microGradientParameters[ 1 ];
        parameterType microGradientYieldValue    = microGradientParameters[ 2 ];
        error = computeHigherOrderDruckerPragerYieldEquation( referenceHigherOrderStress, microGradientCohesion, elasticRightCauchyGreen,
                                                              microGradientFrictionAngle, microGradientBeta, microGradientYieldValue,
                                                              microGradientFlowDirection, dGdMicroGradientCohesion, tmpMat ); 

        return NULL;
    }
}
