/*!
 * tardigrade_micromorphic_elasto_plasticity.cpp
 *
 * An implementation of a elasto-plastic micromorphic constitutive model 
 * following the derivations of Farhad Shahabi in his dissertation.
 */

#include<tardigrade_micromorphic_elasto_plasticity.h>

namespace tardigradeMicromorphicElastoPlasticity{

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation ){
        /*!
         * Assemble the fundamental deformation measures from the degrees of freedom.
         *
         * \param &grad_u: The macro displacement gradient w.r.t. the reference configuration.
         * \param &phi: The micro displacement.
         * \param &grad_phi: The gradient of the micro displacement w.r.t. the reference configuration.
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro deformation
         * \param &gradientMicroDeformation: The gradient of the micro deformation.
         */


        //Extract the degrees of freedom
        variableVector displacementGradient = { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ],
                                                grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ],
                                                grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] };

        variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                             phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                             phi[ 6 ], phi[ 7 ], phi[ 8 ] };

        variableVector gradientMicroDisplacement = { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ],
                                                     grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ],
                                                     grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ],
                                                     grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ],
                                                     grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ],
                                                     grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ],
                                                     grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ],
                                                     grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ],
                                                     grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] };

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation ) );

        return;
    }

    void assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                     const double ( &grad_phi )[ 9 ][ 3 ],
                                                     variableVector &deformationGradient, variableVector &microDeformation,
                                                     variableVector &gradientMicroDeformation, variableMatrix &dFdGradU,
                                                     variableMatrix &dChidPhi, variableMatrix &dGradChidGradPhi ){
        /*!
         * Assemble the fundamental deformation measures from the degrees of freedom.
         *
         * \param &grad_u: The macro displacement gradient w.r.t. the reference configuration.
         * \param &phi: The micro displacement.
         * \param &grad_phi: The gradient of the micro displacement w.r.t. the reference configuration.
         * \param &deformationGradient: The deformation gradient
         * \param &microDeformation: The micro deformation
         * \param &gradientMicroDeformation: The gradient of the micro deformation.
         * \param &dFdGradU: The Jacobian of the deformation gradient w.r.t. the gradient of the displacement
         * \param &dChidPhi: The Jacobian of the micro deformation w.r.t. the micro displacement
         * \param &dGradChidGradPhi: The Jacobian of the gradient of the micro deformation w.r.t.
         *      the gradient of the micro displacement
         */

        const unsigned int dim = 3;
        const unsigned int sot_dim = dim * dim;
        const unsigned int tot_dim = sot_dim * dim;

        //Extract the degrees of freedom
        variableVector displacementGradient = { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ],
                                                grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ],
                                                grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] };

        variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                             phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                             phi[ 6 ], phi[ 7 ], phi[ 8 ] };

        variableVector gradientMicroDisplacement = { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ],
                                                     grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ],
                                                     grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ],
                                                     grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ],
                                                     grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ],
                                                     grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ],
                                                     grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ],
                                                     grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ],
                                                     grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] };

        variableVector _dFdGradU, _dChidPhi, _dGradChidGradPhi;

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient, _dFdGradU ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation, _dChidPhi ) );

        TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation,
                                                                     _dGradChidGradPhi ) );

        dFdGradU = tardigradeVectorTools::inflate( _dFdGradU, sot_dim, sot_dim );

        dChidPhi = tardigradeVectorTools::inflate( _dChidPhi, sot_dim, sot_dim );

        dGradChidGradPhi = tardigradeVectorTools::inflate( _dGradChidGradPhi, tot_dim, tot_dim );

        return;
    }

    void generate_input_variable_string( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                         const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                         const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                         const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                         std::vector< double > &SDVS,
                                         const std::vector< double > &current_ADD_DOF,
                                         const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                         const std::vector< double > &previous_ADD_DOF,
                                         const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                         std::string &input_variables ){
        /*
         * Summarize the input variables in string form for debugging
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacment values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &input_variables: The input variables in string form
         */

        input_variables = "";
        std::stringstream s( input_variables );
        s << std::scientific;

        s << "time:\n";
        for ( auto t = time.begin( ); t != time.end( ); t++ ){ s << " " << *t << ","; }
        s << "\n\nfparams:\n";
        for ( auto f = fparams.begin( ); f != fparams.end( ); f++ ){ s << " " << *f << ","; }
        s << "\n\ncurrent_grad_u:\n";
        for ( unsigned int i = 0; i < 3; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ s << " " << current_grad_u[ i ][ j ] << ","; } s << "\n"; }
        s << "\ncurrent_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ s << " " << current_phi[ i ] << ","; }
        s << "\n\ncurrent_grad_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ s << " " << current_grad_phi[ i ][ j ] << ","; } s << "\n"; }
        s << "\nprevious_grad_u:\n";
        for ( unsigned int i = 0; i < 3; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ s << " "  << previous_grad_u[ i ][ j ] << ","; } s << "\n"; }
        s << "\nprevious_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ s << " " << previous_phi[ i ] << ","; }
        s << "\n\nprevious_grad_phi:\n";
        for ( unsigned int i = 0; i < 9; i++ ){ for ( unsigned int j = 0; j < 3; j++ ){ s << " " << previous_grad_phi[ i ][ j ] << ","; } s << "\n"; }
        s << "\nSDVS:\n";
        for ( auto _s = SDVS.begin( ); _s != SDVS.end( ); _s++ ){ s << " " << *_s << ","; }
        s << "\n\ncurrent_ADD_DOF:\n";
        for ( auto a = current_ADD_DOF.begin( ); a != current_ADD_DOF.end( ); a++ ){ s << " " << *a << ","; }
        s << "\ncurrent_ADD_grad_DOF:\n";
        for ( auto a = current_ADD_grad_DOF.begin( ); a != current_ADD_grad_DOF.end( ); a++ ){ for ( auto g = a->begin( ); g != a->end( ); g++ ){ s << " " << *g << ","; } s << "\n"; }
        s << "\nprevious_ADD_DOF:\n";
        for ( auto a = previous_ADD_DOF.begin( ); a != previous_ADD_DOF.end( ); a++ ){ s << " " << *a << ","; }
        s << "\nprevious_ADD_grad_DOF:\n";
        for ( auto a = previous_ADD_grad_DOF.begin( ); a != previous_ADD_grad_DOF.end( ); a++ ){ for ( auto g = a->begin( ); g != a->end( ); g++ ){ s << " " << *g << ","; } s << "\n"; }

        input_variables = s.str( );

        return;

    }

    int evaluate_hydra_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                              const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                              const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                              const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                              std::vector< double > &SDVS,
                              const std::vector< double > &current_ADD_DOF,
                              const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                              const std::vector< double > &previous_ADD_DOF,
                              const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                              std::vector< double > &current_PK2, std::vector< double > &current_SIGMA, std::vector< double > &current_M,
                              std::vector< std::vector< double > > &ADD_TERMS,
                              std::string &output_message
                              ){
        /*!
         * Evaluate the elasto-plastic constitutive model. Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacment values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &current_PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &current_SIGMA: The current value of the reference micro stress. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &current_M: The current value of the reference higher order stress. The format is
         *     [ \f$M_{111}\f$, \f$M_{112}\f$, \f$M_{113}\f$, \f$M_{121}\f$, \f$M_{122}\f$, \f$M_{123}\f$, \f$M_{131}\f$, \f$M_{132}\f$, \f$M_{133}\f$,
         *       \f$M_{211}\f$, \f$M_{212}\f$, \f$M_{213}\f$, \f$M_{221}\f$, \f$M_{222}\f$, \f$M_{223}\f$, \f$M_{231}\f$, \f$M_{232}\f$, \f$M_{233}\f$,
         *       \f$M_{311}\f$, \f$M_{312}\f$, \f$M_{313}\f$, \f$M_{321}\f$, \f$M_{322}\f$, \f$M_{323}\f$, \f$M_{331}\f$, \f$M_{332}\f$, \f$M_{333}\f$ ]
         * \param &ADD_TERMS: Additional terms ( unused )
         * \param &output_message: The output message string.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        variableType temperature         = 293.15; // Tardigrade doesn't have temperature for micromorphic currently so we're hardcoding these
        variableType previousTemperature = 293.15;

        variableVector currentDeformationGradient,  currentMicroDeformation,  currentGradientMicroDeformation;
        variableVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

        std::string failure_string;
        try{

            /*===============================================
            | Assemble the fundamental deformation measures |
            ================================================*/

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                        currentDeformationGradient, currentMicroDeformation,
                                                        currentGradientMicroDeformation )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                        previousDeformationGradient, previousMicroDeformation,
                                                        previousGradientMicroDeformation )
            )

            // Compute the stress
            try{

                variableVector SDVS_extend( SDVS.size( ) + 5, 0 );
                std::copy( SDVS.begin( ), SDVS.end( ), SDVS_extend.begin( ) );

                hydraMicromorphicElastoPlasticityOptimization hydra( time[ 0 ], time[ 1 ],
                                                                     temperature,                     previousTemperature,
                                                                     currentDeformationGradient,      previousDeformationGradient,
                                                                     currentMicroDeformation,         previousMicroDeformation,
                                                                     currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                                     { }, { },
                                                                     SDVS_extend, fparams, 2, 15, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );

                // Turn on projection
                for ( auto residual_ptr = hydra.getResidualClasses( )->begin( ); residual_ptr != hydra.getResidualClasses( )->end( ); residual_ptr++ ){
                    ( *residual_ptr )->setUseProjection( true );
                }

                hydra.setUseLevenbergMarquardt(false);

                hydra.setGradientBeta( 0.1 );

                hydra.setMaxGradientIterations( 30 );

                hydra.setMaxRelaxedIterations( 10 );

                hydra.setFailureVerbosityLevel( 0 );
                hydra.setFailureOutputScientific( );

                try{
                    hydra.evaluate( true );
                    failure_string += "NON-OPTIMIZE RESULTS:\n\n";
                }catch( std::exception &e ){
                    failure_string += "NON-OPTIMIZE RESULTS:\n\n";
                    failure_string += hydra.getFailureOutput( ) + "\n";
                    tardigradeErrorTools::captureNestedExceptions( e, failure_string );
                    throw;
                }

                current_PK2   = variableVector( hydra.getUnknownVector( )->begin( ) +  0,
                                                hydra.getUnknownVector( )->begin( ) +  9 );

                current_SIGMA = variableVector( hydra.getUnknownVector( )->begin( ) +  9,
                                                hydra.getUnknownVector( )->begin( ) + 18 );

                current_M     = variableVector( hydra.getUnknownVector( )->begin( ) + 18,
                                                hydra.getUnknownVector( )->begin( ) + 45 );

                SDVS          = variableVector( hydra.getUnknownVector( )->begin( ) + 45,
                                                hydra.getUnknownVector( )->begin( ) + 100 );

            }
            catch( std::exception &e ){

                variableVector SDVS_extend( SDVS.size( ) + 5, 0 );
                std::copy( SDVS.begin( ), SDVS.end( ), SDVS_extend.begin( ) );

                hydraMicromorphicElastoPlasticityOptimization hydra( time[ 0 ], time[ 1 ],
                                                                     temperature,                     previousTemperature,
                                                                     currentDeformationGradient,      previousDeformationGradient,
                                                                     currentMicroDeformation,         previousMicroDeformation,
                                                                     currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                                     { }, { },
                                                                     SDVS_extend, fparams, 2, 15, 3, 45, 1e-9, 1e-9, 40, 10, 1e-4, true, 0 );

                hydra.public_setUseSQPSolver( true );

                hydra.setMaxRelaxedIterations( 10 );

                hydra.setFailureVerbosityLevel( 0 );
                hydra.setFailureOutputScientific( );

                try{
                    hydra.evaluate( true );
                }catch( std::exception &e ){
                    failure_string += "OPTIMIZE RESULTS:\n\n";
                    failure_string += hydra.getFailureOutput( ) + "\n";
                    tardigradeErrorTools::captureNestedExceptions( e, failure_string );
                    throw;
                }

                current_PK2   = variableVector( hydra.getUnknownVector( )->begin( ) +  0,
                                                hydra.getUnknownVector( )->begin( ) +  9 );

                current_SIGMA = variableVector( hydra.getUnknownVector( )->begin( ) +  9,
                                                hydra.getUnknownVector( )->begin( ) + 18 );

                current_M     = variableVector( hydra.getUnknownVector( )->begin( ) + 18,
                                                hydra.getUnknownVector( )->begin( ) + 45 );

                SDVS          = variableVector( hydra.getUnknownVector( )->begin( ) + 45,
                                                hydra.getUnknownVector( )->begin( ) + 100 );

            }

            for ( unsigned int i = 0; i < 3; i++ ){

                SDVS[ 3 * i + i + 0 ] -= 1;

                SDVS[ 3 * i + i + 9 ] -= 1;

            }

        }
        catch( tardigradeHydra::convergence_error &e ){

            //Convergence error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            output_message += "INPUT PARAMETERS FOLLOW:\n" + input_variables + "\n";

            output_message += "ADDITIONAL MESSAGES:\n\n" + failure_string + "\n\n";

            return 1;

        }
        catch( std::exception &e ){

            //Fatal error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            output_message += "INPUT PARAMETERS FOLLOW:\n" + input_variables + "\n";

            output_message += "ADDITIONAL MESSAGES:\n\n" + failure_string + "\n\n";

#ifdef TARDIGRADE_FATAL_AS_CONVERGENCE
            return 1;
#else
            return 2;
#endif

        }

        //No errors in calculation.
        return 0;

    }

    void assembleJacobians( const variableVector *dXdD,      const unsigned int numConfigurationUnknowns,
                            const variableMatrix &dFdGradU,  const variableMatrix &dChidPhi, const variableMatrix &dGradChidGradPhi,
                            variableMatrix &DPK2Dgrad_u,     variableMatrix &DPK2Dphi,       variableMatrix &DPK2Dgrad_phi,
                            variableMatrix &DSIGMADgrad_u,   variableMatrix &DSIGMADphi,     variableMatrix &DSIGMADgrad_phi,
                            variableMatrix &DMDgrad_u,       variableMatrix &DMDphi,         variableMatrix &DMDgrad_phi,
                            std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS ){
        /*!
         * Assemble the Jacobians of the stress measures w.r.t. the deformation
         *
         * \param *dXdD: The derivative of the unknown vector w.r.t. the deformation
         * \param numConfigurationUnknowns: The number of unknowns in each configuration
         * \param &dFdGradU: The derivative of the deformation gradient w.r.t. the displacement gradient
         * \param &dChidPhi: The derivative of the micro deformation w.r.t. the micro displacement
         * \param &dGradChidGradPhi: The derivative of the gradient of the micro deformation w.r.t. the gradient of the micro displacement
         * \param &DPK2Dgrad_u: The Jacobian of the PK2 stress w.r.t. the 
         *     gradient of macro displacement.
         * \param &DPK2Dphi: The Jacobian of the PK2 stress w.r.t. the
         *     micro displacement.
         * \param &DPK2Dgrad_phi: The Jacobian of the PK2 stress w.r.t.
         *     the gradient of the micro displacement.
         * \param &DSIGMADgrad_u: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the macro displacement.
         * \param &DSIGMADphi: The Jacobian of the reference symmetric micro
         *     stress w.r.t. the micro displacement.
         * \param &DSIGMADgrad_phi: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the micro displacement.
         * \param &DMDgrad_u: The Jacobian of the reference higher order
         *     stress w.r.t. the gradient of the macro displacement.
         * \param &DMDphi: The Jacobian of the reference higher order stress
         *     w.r.t. the micro displacement.
         * \param &DMDgrad_phi: The Jacobian of the reference higher order stress
         *     w.r.t. the gradient of the micro displacement.
         * \param &ADD_JACOBIANS: The jacobians of the additional
         *     terms w.r.t. the deformation. This is currently being used to support the gradient enhanced damage work
         *     by returning the Jacobians of the plastic deformation gradients w.r.t. the deformation measures. The
         *     ordering is: DFpDgrad_u, DFpDphi, DFpDgrad_phi, DchipDgrad_u, DchipDphi, DchipDgrad_phi, Dgrad_chipDgrad_u, Dgrad_chipDchi, Dgrad_chipDgrad_chi
         */

        // Compute the consistent tangents
        DPK2Dgrad_u     = variableMatrix(  9, variableVector(  9, 0 ) );

        DSIGMADgrad_u   = variableMatrix(  9, variableVector(  9, 0 ) );

        DMDgrad_u       = variableMatrix( 27, variableVector(  9, 0 ) );

        DPK2Dphi        = variableMatrix(  9, variableVector(  9, 0 ) );

        DSIGMADphi      = variableMatrix(  9, variableVector(  9, 0 ) );

        DMDphi          = variableMatrix( 27, variableVector(  9, 0 ) );

        DPK2Dgrad_phi   = variableMatrix(  9, variableVector( 27, 0 ) );

        DSIGMADgrad_phi = variableMatrix(  9, variableVector( 27, 0 ) );

        DMDgrad_phi     = variableMatrix( 27, variableVector( 27, 0 ) );

        ADD_JACOBIANS   = std::vector< variableMatrix >( 9 );

        ADD_JACOBIANS[ 0 ] = variableMatrix(  9, variableVector(  9, 0 ) );

        ADD_JACOBIANS[ 1 ] = variableMatrix(  9, variableVector(  9, 0 ) );

        ADD_JACOBIANS[ 2 ] = variableMatrix(  9, variableVector( 27, 0 ) );

        ADD_JACOBIANS[ 3 ] = variableMatrix(  9, variableVector(  9, 0 ) );

        ADD_JACOBIANS[ 4 ] = variableMatrix(  9, variableVector(  9, 0 ) );

        ADD_JACOBIANS[ 5 ] = variableMatrix(  9, variableVector( 27, 0 ) );

        ADD_JACOBIANS[ 6 ] = variableMatrix( 27, variableVector(  9, 0 ) );

        ADD_JACOBIANS[ 7 ] = variableMatrix( 27, variableVector(  9, 0 ) );

        ADD_JACOBIANS[ 8 ] = variableMatrix( 27, variableVector( 27, 0 ) );

        for ( unsigned int i = 0; i < 9; i++ ){

            for ( unsigned int j = 0; j < 9; j++ ){

                for ( unsigned int k = 0; k < 9; k++ ){

                    DPK2Dgrad_u[ i ][ j ]   += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 ) + 0 + k ] * dFdGradU[ k ][ j ];

                    DPK2Dphi[ i ][ j ]      += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 ) + 9 + k ] * dChidPhi[ k ][ j ];

                    DSIGMADgrad_u[ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 ) + 0 + k ] * dFdGradU[ k ][ j ];

                    DSIGMADphi[ i ][ j ]    += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 ) + 9 + k ] * dChidPhi[ k ][ j ];

                    ADD_JACOBIANS[ 0 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 + numConfigurationUnknowns ) + 0 + k ] * dFdGradU[ k ][ j ];

                    ADD_JACOBIANS[ 1 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 + numConfigurationUnknowns ) + 9 + k ] * dChidPhi[ k ][ j ];

                    ADD_JACOBIANS[ 3 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 + numConfigurationUnknowns ) + 0 + k ] * dFdGradU[ k ][ j ];

                    ADD_JACOBIANS[ 4 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 + numConfigurationUnknowns ) + 9 + k ] * dChidPhi[ k ][ j ];

                }

            }

            for ( unsigned int j = 0; j < 27; j++ ){

                for ( unsigned int k = 0; k < 27; k++ ){

                    DPK2Dgrad_phi[ i ][ j ]   += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                    DSIGMADgrad_phi[ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                    ADD_JACOBIANS[ 2 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 0 + numConfigurationUnknowns ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                    ADD_JACOBIANS[ 5 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 9 + numConfigurationUnknowns ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                }

            }

        }

        for ( unsigned int i = 0; i < 27; i++ ){

            for ( unsigned int j = 0; j < 9; j++ ){

                for ( unsigned int k = 0; k < 9; k++ ){

                    DMDgrad_u[ i ][ j ]   += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 ) + 0 + k ] * dFdGradU[ k ][ j ];

                    DMDphi[ i ][ j ]      += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 ) + 9 + k ] * dChidPhi[ k ][ j ];

                    ADD_JACOBIANS[ 6 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 + numConfigurationUnknowns ) + 0 + k ] * dFdGradU[ k ][ j ];

                    ADD_JACOBIANS[ 7 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 + numConfigurationUnknowns ) + 9 + k ] * dChidPhi[ k ][ j ];

                }

            }

            for ( unsigned int j = 0; j < 27; j++ ){

                for ( unsigned int k = 0; k < 27; k++ ){

                    DMDgrad_phi[ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                    ADD_JACOBIANS[ 8 ][ i ][ j ] += ( *dXdD )[ numConfigurationUnknowns * ( i + 18 + numConfigurationUnknowns ) + 18 + k ] * dGradChidGradPhi[ k ][ j ];

                }

            }

        }

    }

    int evaluate_hydra_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                              const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                              const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                              const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                              std::vector< double > &SDVS,
                              const std::vector< double > &current_ADD_DOF,
                              const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                              const std::vector< double > &previous_ADD_DOF,
                              const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                              std::vector< double > &current_PK2, std::vector< double > &current_SIGMA, std::vector< double > &current_M,
                              std::vector< std::vector< double > > &DPK2Dgrad_u,   std::vector< std::vector< double > > &DPK2Dphi,
                              std::vector< std::vector< double > > &DPK2Dgrad_phi,
                              std::vector< std::vector< double > > &DSIGMADgrad_u, std::vector< std::vector< double > > &DSIGMADphi,
                              std::vector< std::vector< double > > &DSIGMADgrad_phi,
                              std::vector< std::vector< double > > &DMDgrad_u,     std::vector< std::vector< double > > &DMDphi,
                              std::vector< std::vector< double > > &DMDgrad_phi,
                              std::vector< std::vector< double > > &ADD_TERMS,
                              std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                              std::string &output_message
                            ){
        /*!
         * Evaluate the elasto-plastic constitutive model. Note the format of the header changed to provide a 
         * consistant interface with the material model library.
         *
         * \param &time: The current time and the timestep
         *     [ current_t, dt ]
         * \param &fparams: The parameters for the constitutive model
         *     [ num_Amatrix_parameters, Amatrix_parameters, num_Bmatrix_parameters, Bmatrix_parameters,
         *       num_Cmatrix_parameters, Cmatrix_parameters, num_Dmatrix_parameters, Dmatrix_parameters,
         *       num_macroHardeningParameters, macroHardeningParameters,
         *       num_microHardeningParameters, microHardeningParameters,
         *       num_microGradientHardeningParameters, microGradientHardeningParameters,
         *       num_macroFlowParameters, macroFlowParameters,
         *       num_microFlowParameters, microFlowParameters,
         *       num_microGradientFlowParameters, microGradientFlowParameters,
         *       num_macroYieldParameters, macroYieldParameters,
         *       num_microYieldParameters, microYieldParameters,
         *       num_microGradientYieldParameters, microGradientYieldParameters,
         *       alphaMacro, alphaMicro, alphaMicroGradient,
         *       relativeTolerance, absoluteTolerance ]
         *
         * \param &current_grad_u: The current displacement gradient
         *     Assumed to be of the form [ [ \f$u_{1,1}\f$, \f$u_{1,2}\f$, \f$u_{1,3}\f$ ],
         *                                 [ \f$u_{2,1}\f$, \f$u_{2,2}\f$, \f$u_{2,3}\f$ ],
         *                                 [ \f$u_{3,1}\f$, \f$u_{3,2}\f$, \f$u_{3,3}\f$ ] ]
         * \param &current_phi: The current micro displacement values.
         *     Assumed to be of the form [ \f$\phi_{11}\f$, \f$\phi_{12}\f$, \f$\phi_{13}\f$, \f$\phi_{21}\f$, \f$\phi_{22}\f$, \f$\phi_{23}\f$, \f$\phi_{31}\f$, \f$\phi_{32}\f$, \f$\phi_{33}\f$ ]
         * \param &current_grad_phi: The current micro displacement gradient
         *     Assumed to be of the form [ [ \f$\phi_{11,1}\f$, \f$\phi_{11,2}\f$, \f$\phi_{11,3}\f$ ],
         *                                 [ \f$\phi_{12,1}\f$, \f$\phi_{12,2}\f$, \f$\phi_{12,3}\f$ ],
         *                                 [ \f$\phi_{13,1}\f$, \f$\phi_{13,2}\f$, \f$\phi_{13,3}\f$ ],
         *                                 [ \f$\phi_{21,1}\f$, \f$\phi_{21,2}\f$, \f$\phi_{21,3}\f$ ],
         *                                 [ \f$\phi_{22,1}\f$, \f$\phi_{22,2}\f$, \f$\phi_{22,3}\f$ ],
         *                                 [ \f$\phi_{23,1}\f$, \f$\phi_{23,2}\f$, \f$\phi_{23,3}\f$ ],
         *                                 [ \f$\phi_{31,1}\f$, \f$\phi_{31,2}\f$, \f$\phi_{31,3}\f$ ],
         *                                 [ \f$\phi_{32,1}\f$, \f$\phi_{32,2}\f$, \f$\phi_{32,3}\f$ ],
         *                                 [ \f$\phi_{33,1}\f$, \f$\phi_{33,2}\f$, \f$\phi_{33,3}\f$ ] ]
         * \param &previous_grad_u: The previous displacement gradient.
         * \param &previous_phi: The previous micro displacement.
         * \param &previous_grad_phi: The previous micro displacement gradient.
         * \param &SDVS: The previously converged values of the state variables
         *     [ previousMacroStrainISV, previousMicroStrainISV, previousMicroGradientStrainISV,
         *       previousMacroGamma, previousMicroGamma, previousMicroGradientGamma,
         *       previousPlasticDeformationGradient - eye, previousPlasticMicroDeformation - eye,
         *       previousPlasticMicroGradient ]
         * \param &current_ADD_DOF: The current values of the additional degrees of freedom ( unused )
         * \param &current_ADD_grad_DOF: The current values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &previous_ADD_DOF: The previous values of the additional degrees of freedom ( unused )
         * \param &previous_ADD_grad_DOF: The previous values of the gradients of the 
         *     additional degrees of freedom ( unused )
         * \param &current_PK2: The current value of the second Piola Kirchhoff stress tensor. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &current_SIGMA: The current value of the reference micro stress. The format is
         *     [ \f$S_{11}\f$, \f$S_{12}\f$, \f$S_{13}\f$, \f$S_{21}\f$, \f$S_{22}\f$, \f$S_{23}\f$, \f$S_{31}\f$, \f$S_{32}\f$, \f$S_{33}\f$ ]
         * \param &current_M: The current value of the reference higher order stress. The format is
         *     [ \f$M_{111}\f$, \f$M_{112}\f$, \f$M_{113}\f$, \f$M_{121}\f$, \f$M_{122}\f$, \f$M_{123}\f$, \f$M_{131}\f$, \f$M_{132}\f$, \f$M_{133}\f$,
         *       \f$M_{211}\f$, \f$M_{212}\f$, \f$M_{213}\f$, \f$M_{221}\f$, \f$M_{222}\f$, \f$M_{223}\f$, \f$M_{231}\f$, \f$M_{232}\f$, \f$M_{233}\f$,
         *       \f$M_{311}\f$, \f$M_{312}\f$, \f$M_{313}\f$, \f$M_{321}\f$, \f$M_{322}\f$, \f$M_{323}\f$, \f$M_{331}\f$, \f$M_{332}\f$, \f$M_{333}\f$ ]
         * \param &DPK2Dgrad_u: The Jacobian of the PK2 stress w.r.t. the 
         *     gradient of macro displacement.
         * \param &DPK2Dphi: The Jacobian of the PK2 stress w.r.t. the
         *     micro displacement.
         * \param &DPK2Dgrad_phi: The Jacobian of the PK2 stress w.r.t.
         *     the gradient of the micro displacement.
         * \param &DSIGMADgrad_u: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the macro displacement.
         * \param &DSIGMADphi: The Jacobian of the reference symmetric micro
         *     stress w.r.t. the micro displacement.
         * \param &DSIGMADgrad_phi: The Jacobian of the reference symmetric
         *     micro stress w.r.t. the gradient of the micro displacement.
         * \param &DMDgrad_u: The Jacobian of the reference higher order
         *     stress w.r.t. the gradient of the macro displacement.
         * \param &DMDphi: The Jacobian of the reference higher order stress
         *     w.r.t. the micro displacement.
         * \param &DMDgrad_phi: The Jacobian of the reference higher order stress
         *     w.r.t. the gradient of the micro displacement.
         * \param &ADD_TERMS: Additional terms ( unused )
         * \param &ADD_JACOBIANS: The jacobians of the additional
         *     terms w.r.t. the deformation. This is currently being used to support the gradient enhanced damage work
         *     by returning the Jacobians of the plastic deformation gradients w.r.t. the deformation measures. The
         *     ordering is: DFpDgrad_u, DFpDphi, DFpDgrad_phi, DchipDgrad_u, DchipDphi, DchipDgrad_phi, Dgrad_chipDgrad_u, Dgrad_chipDchi, Dgrad_chipDgrad_chi
         * \param &output_message: The output message string.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        variableType temperature         = 293.15; // Tardigrade doesn't have temperature for micromorphic currently so we're hardcoding these
        variableType previousTemperature = 293.15;

        variableVector currentDeformationGradient,  currentMicroDeformation,  currentGradientMicroDeformation;
        variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

        variableVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;
        variableMatrix previousdFdGradU, previousdChidPhi, previousdGradChidGradPhi;

        std::string failure_string;
        try{

            /*===============================================
            | Assemble the fundamental deformation measures |
            ================================================*/

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                        currentDeformationGradient, currentMicroDeformation,
                                                        currentGradientMicroDeformation,
                                                        dFdGradU, dChidPhi, dGradChidGradPhi )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                        previousDeformationGradient, previousMicroDeformation,
                                                        previousGradientMicroDeformation,
                                                        previousdFdGradU, previousdChidPhi, previousdGradChidGradPhi )
            )

            // Compute the stress
            try{
                variableVector SDVS_extend( SDVS.size( ) + 5, 0 );
                std::copy( SDVS.begin( ), SDVS.end( ), SDVS_extend.begin( ) );

                hydraMicromorphicElastoPlasticityOptimization hydra( time[ 0 ], time[ 1 ],
                                                                     temperature,                     previousTemperature,
                                                                     currentDeformationGradient,      previousDeformationGradient,
                                                                     currentMicroDeformation,         previousMicroDeformation,
                                                                     currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                                     { }, { },
                                                                     SDVS_extend, fparams, 2, 15, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );

                // Turn on projection
                for ( auto residual_ptr = hydra.getResidualClasses( )->begin( ); residual_ptr != hydra.getResidualClasses( )->end( ); residual_ptr++ ){
                    ( *residual_ptr )->setUseProjection( true );
                }

                hydra.setUseLevenbergMarquardt(false);

                hydra.setGradientBeta( 0.1 );

                hydra.setMaxGradientIterations( 30 );

                hydra.setMaxRelaxedIterations( 10 );

                hydra.setFailureVerbosityLevel( 0 );
                hydra.setFailureOutputScientific( );

                try{
                    hydra.evaluate( true );

                    current_PK2   = variableVector( hydra.getUnknownVector( )->begin( ) +  0,
                                                    hydra.getUnknownVector( )->begin( ) +  9 );

                    current_SIGMA = variableVector( hydra.getUnknownVector( )->begin( ) +  9,
                                                    hydra.getUnknownVector( )->begin( ) + 18 );

                    current_M     = variableVector( hydra.getUnknownVector( )->begin( ) + 18,
                                                    hydra.getUnknownVector( )->begin( ) + 45 );

                    SDVS          = variableVector( hydra.getUnknownVector( )->begin( ) + 45,
                                                    hydra.getUnknownVector( )->begin( ) + 100 );

                    hydra.computeTangents( );

                    assembleJacobians( hydra.getFlatdXdD( ), *hydra.getConfigurationUnknownCount( ),
                                       dFdGradU,      dChidPhi,   dGradChidGradPhi,
                                       DPK2Dgrad_u,   DPK2Dphi,   DPK2Dgrad_phi,
                                       DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                       DMDgrad_u,     DMDphi,     DMDgrad_phi,
                                       ADD_JACOBIANS );

                }catch( std::exception &e ){
                    failure_string += "NON-OPTIMIZE RESULTS:\n\n";
                    failure_string += hydra.getFailureOutput( ) + "\n";
                    tardigradeErrorTools::captureNestedExceptions( e, failure_string );
                    throw;
                }

            }
            catch( std::exception &e ){

                variableVector SDVS_extend( SDVS.size( ) + 5, 0 );
                std::copy( SDVS.begin( ), SDVS.end( ), SDVS_extend.begin( ) );

                hydraMicromorphicElastoPlasticityOptimization hydra( time[ 0 ], time[ 1 ],
                                                                     temperature,                     previousTemperature,
                                                                     currentDeformationGradient,      previousDeformationGradient,
                                                                     currentMicroDeformation,         previousMicroDeformation,
                                                                     currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                                     { }, { },
                                                                     SDVS_extend, fparams, 2, 15, 3, 45, 1e-9, 1e-9, 40, 10, 1e-4, true, 0 );

                hydra.public_setUseSQPSolver( true );

                hydra.setMaxRelaxedIterations( 10 );

                hydra.setFailureVerbosityLevel( 0 );
                hydra.setFailureOutputScientific( );

                try{
                    hydra.evaluate( true );

                    current_PK2   = variableVector( hydra.getUnknownVector( )->begin( ) +  0,
                                                    hydra.getUnknownVector( )->begin( ) +  9 );

                    current_SIGMA = variableVector( hydra.getUnknownVector( )->begin( ) +  9,
                                                    hydra.getUnknownVector( )->begin( ) + 18 );

                    current_M     = variableVector( hydra.getUnknownVector( )->begin( ) + 18,
                                                    hydra.getUnknownVector( )->begin( ) + 45 );

                    SDVS          = variableVector( hydra.getUnknownVector( )->begin( ) + 45,
                                                    hydra.getUnknownVector( )->begin( ) + 100 );

                    hydra.computeTangents( );

                    assembleJacobians( hydra.getFlatdXdD( ), *hydra.getConfigurationUnknownCount( ),
                                       dFdGradU,      dChidPhi,   dGradChidGradPhi,
                                       DPK2Dgrad_u,   DPK2Dphi,   DPK2Dgrad_phi,
                                       DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                       DMDgrad_u,     DMDphi,     DMDgrad_phi,
                                       ADD_JACOBIANS );

                }catch( std::exception &e ){
                    failure_string += "OPTIMIZE RESULTS:\n\n";
                    failure_string += hydra.getFailureOutput( ) + "\n";
                    tardigradeErrorTools::captureNestedExceptions( e, failure_string );
                    throw;
                }

            }

            for ( unsigned int i = 0; i < 3; i++ ){

                SDVS[ 3 * i + i + 0 ] -= 1;

                SDVS[ 3 * i + i + 9 ] -= 1;

            }

        }
        catch( tardigradeHydra::convergence_error &e ){

            //Convergence error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            output_message += "INPUT PARAMETERS FOLLOW:\n" + input_variables + "\n";

            output_message += "ADDITIONAL MESSAGES:\n\n" + failure_string + "\n\n";

            return 1;

        }
        catch( std::exception &e ){

            //Fatal error
            std::string input_variables;
            generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                            previous_grad_u, previous_phi, previous_grad_phi,
                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF,
                                            input_variables );

            tardigradeErrorTools::captureNestedExceptions( e, output_message );

            output_message += "INPUT PARAMETERS FOLLOW:\n" + input_variables + "\n";

            output_message += "ADDITIONAL MESSAGES:\n\n" + failure_string + "\n\n";

#ifdef TARDIGRADE_FATAL_AS_CONVERGENCE
            return 1;
#else
            return 2;
#endif

        }

        //No errors in calculation.
        return 0;

    }

}
