/*!
 * tardigrade_micromorphic_elasto_plasticity.h
 *
 * An implementation of a elasto-plastic micromorphic constitutive model 
 * following the derivations of Farhad Shahabi in his dissertation.
 */

#ifndef TARDIGRADE_MICROMORPHIC_ELASTO_PLASTICITY_H
#define TARDIGRADE_MICROMORPHIC_ELASTO_PLASTICITY_H

#include<tardigrade_error_tools.h>
#define USE_EIGEN
#include<tardigrade_vector_tools.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
#include<tardigrade_micromorphic_linear_elasticity.h>
#include<tardigrade_hydraMicromorphicLinearElasticity.h>
#include<tardigrade_hydraMicromorphicDruckerPragerPlasticity.h>
#include<tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization.h>
//#include<micromorphic_material_library.h>

namespace tardigradeMicromorphicElastoPlasticity{

    typedef tardigradeMicromorphicTools::variableType variableType; //!< Definition of the variable type. Should be used for values which vary over the course of the solution
    typedef tardigradeMicromorphicTools::variableVector variableVector; //!< Definition of a vector of variables
    typedef tardigradeMicromorphicTools::variableMatrix variableMatrix; //!< Definition of a matrix of variables

    typedef tardigradeMicromorphicTools::parameterType parameterType; //!< Definition of the parameter type. Should be used for user-changed parameters.
    typedef tardigradeMicromorphicTools::parameterVector parameterVector; //!< Definition of a vector of parameters
    typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix; //!< Definition of a matrix of parameters

    typedef tardigradeMicromorphicTools::constantType constantType; //!< Definition of the constant type. Should be expected to never change.
    typedef tardigradeMicromorphicTools::constantVector constantVector; //!< Definition of a vector of constants
    typedef tardigradeMicromorphicTools::constantMatrix constantMatrix; //!< Definition of a matrix of constants

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation );

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation,
                                                     variableMatrix &dFdGradU, variableMatrix &dChidPhi,
                                                     variableMatrix &dGradChidGradPhi );

    //! Define the hydra version of the micromorphic elasto plasticity model
    class hydraMicromorphicElastoPlasticity : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            //! The elasticity residual class
            tardigradeHydra::micromorphicLinearElasticity::residual elasticity; //!< The elasticity configuration

            tardigradeHydra::micromorphicDruckerPragerPlasticity::residual plasticity; //!< The plasticity configuration

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4,
                                                                 5, 6, 7, 8, 9 }; //!< The indices of the state variables

            static constexpr unsigned int numPlasticParameterCollections = 9;

            static constexpr unsigned int numElasticParameters = 24;

            const unsigned int getNumPlasticParameters( ){

                unsigned int numPlasticParameters = 0;

                for ( unsigned int i = 0; i < numPlasticParameterCollections; i++ ){

                    numPlasticParameters += ( 1 + ( *getParameters( ) )[ numPlasticParameters ] );

                }

                return numPlasticParameters;

            }

            variableVector getPlasticParameters( ){
                /*!
                 * Get the plastic parameters from the parameter vector
                 */

                const unsigned int numPlasticParameters = getNumPlasticParameters( );

                return variableVector( getParameters( )->begin( ), getParameters( )->begin( ) + numPlasticParameters );

            }

            variableVector getElasticParameters( ){
                /*!
                 * Get the elastic parameters from the parameter vector
                 */

                const unsigned int numPlasticParameters = getNumPlasticParameters( );

                return variableVector( getParameters( )->begin( ) + numPlasticParameters,
                                       getParameters( )->begin( ) + numPlasticParameters + numElasticParameters );

            }

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{
                /*!
                 * Set the vector of residual classes (in this case, only elasticity)
                 */

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                TARDIGRADE_ERROR_TOOLS_CATCH( elasticity = tardigradeHydra::micromorphicLinearElasticity::residual( this, *getConfigurationUnknownCount( ), getElasticParameters( ) ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticity = tardigradeHydra::micromorphicDruckerPragerPlasticity::residual( this, *getConfigurationUnknownCount( ) + 10, 1, stateVariableIndices, getPlasticParameters( ), 1.0 ) )

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    class hydraMicromorphicElastoPlasticityOptimization : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            //! The elasticity residual class
            tardigradeHydra::micromorphicLinearElasticity::residual elasticity; //!< The elasticity configuration

            tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual plasticity; //!< The plasticity configuration

            std::vector< unsigned int > stateVariableIndices = { 0,  1,  2,  3,  4,
                                                                 5,  6,  7,  8,  9,
                                                                10, 11, 12, 13, 14 }; //!< The indices of the state variables

            static constexpr unsigned int numPlasticParameterCollections = 9;

            static constexpr unsigned int numElasticParameters = 24;

            const unsigned int getNumPlasticParameters( ){

                unsigned int numPlasticParameters = 0;

                for ( unsigned int i = 0; i < numPlasticParameterCollections; i++ ){

                    numPlasticParameters += ( 1 + ( *getParameters( ) )[ numPlasticParameters ] );

                }

                return numPlasticParameters;

            }

            variableVector getPlasticParameters( ){
                /*!
                 * Get the plastic parameters from the parameter vector
                 */

                const unsigned int numPlasticParameters = getNumPlasticParameters( );

                return variableVector( getParameters( )->begin( ), getParameters( )->begin( ) + numPlasticParameters );

            }

            variableVector getElasticParameters( ){
                /*!
                 * Get the elastic parameters from the parameter vector
                 */

                const unsigned int numPlasticParameters = getNumPlasticParameters( );

                return variableVector( getParameters( )->begin( ) + numPlasticParameters,
                                       getParameters( )->begin( ) + numPlasticParameters + numElasticParameters );

            }

            void public_setUseSQPSolver( const bool value ){

                setUseSQPSolver( value );

            }

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{
                /*!
                 * Set the vector of residual classes (in this case, only elasticity)
                 */

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                TARDIGRADE_ERROR_TOOLS_CATCH( elasticity = tardigradeHydra::micromorphicLinearElasticity::residual( this, *getConfigurationUnknownCount( ), getElasticParameters( ) ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticity = tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual( this, *getConfigurationUnknownCount( ) + 15, 1, stateVariableIndices, getPlasticParameters( ), 1.0 ) )

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };
 
    void generate_input_variable_string( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                         const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                         const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                         const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                         std::vector< double > &SDVS,
                                         const std::vector< double > &current_ADD_DOF,
                                         const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                         const std::vector< double > &previous_ADD_DOF,
                                         const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                         std::string &input_variables );

    int evaluate_hydra_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                              const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                              const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                              const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                              std::vector< double > &SDVS,
                              const std::vector< double > &current_ADD_DOF,
                              const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                              const std::vector< double > &previous_ADD_DOF,
                              const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                              std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                              std::vector< std::vector< double > > &ADD_TERMS,
                              std::string &output_message
                            );

    int evaluate_hydra_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                              const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                              const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                              const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                              std::vector< double > &SDVS,
                              const std::vector< double > &current_ADD_DOF,
                              const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                              const std::vector< double > &previous_ADD_DOF,
                              const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                              std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                              std::vector< std::vector< double > > &DPK2Dgrad_u,   std::vector< std::vector< double > > &DPK2Dphi,
                              std::vector< std::vector< double > > &DPK2Dgrad_phi,
                              std::vector< std::vector< double > > &DSIGMADgrad_u, std::vector< std::vector< double > > &DSIGMADphi,
                              std::vector< std::vector< double > > &DSIGMADgrad_phi,
                              std::vector< std::vector< double > > &DMDgrad_u,     std::vector< std::vector< double > > &DMDphi,
                              std::vector< std::vector< double > > &DMDgrad_phi,
                              std::vector< std::vector< double > > &ADD_TERMS,
                              std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                              std::string &output_message
                            );


}

#endif
