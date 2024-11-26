//Tests for tardigrade_constitutive_tools

#include<tardigrade_micromorphic_elasto_plasticity.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>

#define BOOST_TEST_MODULE test_tardigrade_micromorphic_elasto_plasticity
#include <boost/test/included/unit_test.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeMicromorphicTools::constantType constantType;
typedef tardigradeMicromorphicTools::constantVector constantVector;
typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

typedef tardigradeMicromorphicTools::parameterType parameterType;
typedef tardigradeMicromorphicTools::parameterVector parameterVector;
typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

typedef tardigradeMicromorphicTools::variableType variableType;
typedef tardigradeMicromorphicTools::variableVector variableVector;
typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

BOOST_AUTO_TEST_CASE( testAssembleFundamentalDeformationMeasures ){
    /*!
     * Assemble the fundamental deformation measures from the degrees of freedom.
     *
     */

    double grad_u[ 3 ][ 3 ] = { { 1, 2, 3 },
                                { 4, 5, 6 },
                                { 7, 8, 9 } };

    double phi[ 9 ] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    double grad_phi[ 9 ][ 3 ] = { {  1,  2,  3 },
                                  {  4,  5,  6 },
                                  {  7,  8,  9 },
                                  { 10, 11, 12 },
                                  { 13, 14, 15 },
                                  { 16, 17, 18 },
                                  { 19, 20, 21 },
                                  { 22, 23, 24 },
                                  { 25, 26, 27 } };

    variableVector answerDeformationGradient = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerMicroDeformation = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector answerGradientMicroDeformation = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                     10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                     19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector resultF, resultChi, resultGradChi;

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                          resultF, resultChi, resultGradChi );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultF, answerDeformationGradient ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultChi, answerMicroDeformation ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultGradChi, answerGradientMicroDeformation ) );

    //Test the Jacobians
    variableVector resultFJ, resultChiJ, resultGradChiJ;
    variableMatrix dFdGradU, dChidPhi, dGradChidGradPhi;

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, grad_phi,
                                                                          resultFJ, resultChiJ, resultGradChiJ,
                                                                          dFdGradU, dChidPhi, dGradChidGradPhi );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultFJ, answerDeformationGradient ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultChiJ, answerMicroDeformation ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultGradChiJ, answerGradientMicroDeformation ) );

    //Test the jacobians w.r.t. the gradient of the displacement
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta( 3, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_u[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 3 ][ 3 ] =
        { 
                { grad_u[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] + delta[ 2 ][ 2 ] }
        };

        double negative_perturb[ 3 ][ 3 ] =
        { 
                { grad_u[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_u[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_u[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_u[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_u[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_u[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] - delta[ 2 ][ 2 ] }
        };

        tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( positive_perturb, phi, grad_phi,
                                                                              FP, chiP, gradChiP );

        tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( negative_perturb, phi, grad_phi,
                                                                                      FM, chiM, gradChiM );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dFdGradU[ j ][ i ] ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 9; i++ ){
        constantVector delta( 9, 0 );

        delta[ i ] = eps * fabs( phi[ i ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ] = { phi[ 0 ] + delta[ 0 ], phi[ 1 ] + delta[ 1 ], phi[ 2 ] + delta[ 2 ],
                                         phi[ 3 ] + delta[ 3 ], phi[ 4 ] + delta[ 4 ], phi[ 5 ] + delta[ 5 ],
                                         phi[ 6 ] + delta[ 6 ], phi[ 7 ] + delta[ 7 ], phi[ 8 ] + delta[ 8 ] };

        double negative_perturb[ 9 ] = { phi[ 0 ] - delta[ 0 ], phi[ 1 ] - delta[ 1 ], phi[ 2 ] - delta[ 2 ],
                                         phi[ 3 ] - delta[ 3 ], phi[ 4 ] - delta[ 4 ], phi[ 5 ] - delta[ 5 ],
                                         phi[ 6 ] - delta[ 6 ], phi[ 7 ] - delta[ 7 ], phi[ 8 ] - delta[ 8 ] };

        tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, positive_perturb, grad_phi,
                                                                                      FP, chiP, gradChiP );

        tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, negative_perturb, grad_phi,
                                                                                      FM, chiM, gradChiM );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dChidPhi[ j ][ i ] ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }
    }

    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        unsigned int ii, ij;
        ii = ( int )( i / 3 );
        ij = i % 3;
        delta[ ii ][ ij ] = eps * fabs( grad_phi[ ii ][ ij ] ) + eps;

        variableVector FP, chiP, gradChiP;
        variableVector FM, chiM, gradChiM;

        double positive_perturb[ 9 ][ 3 ] =
        { 
                { grad_phi[ 0 ][ 0 ] + delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] + delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] + delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] + delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] + delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] + delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] + delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] + delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] + delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] + delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] + delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] + delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] + delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] + delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] + delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] + delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] + delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] + delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] + delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] + delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] + delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] + delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] + delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] + delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] + delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] + delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] + delta[ 8 ][ 2 ] }
        };

        double negative_perturb[ 9 ][ 3 ] =
        { 
                { grad_phi[ 0 ][ 0 ] - delta[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ] - delta[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] - delta[ 0 ][ 2 ] },
                { grad_phi[ 1 ][ 0 ] - delta[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ] - delta[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] - delta[ 1 ][ 2 ] },
                { grad_phi[ 2 ][ 0 ] - delta[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ] - delta[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] - delta[ 2 ][ 2 ] },
                { grad_phi[ 3 ][ 0 ] - delta[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ] - delta[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] - delta[ 3 ][ 2 ] },
                { grad_phi[ 4 ][ 0 ] - delta[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ] - delta[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] - delta[ 4 ][ 2 ] },
                { grad_phi[ 5 ][ 0 ] - delta[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ] - delta[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] - delta[ 5 ][ 2 ] },
                { grad_phi[ 6 ][ 0 ] - delta[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ] - delta[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] - delta[ 6 ][ 2 ] },
                { grad_phi[ 7 ][ 0 ] - delta[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ] - delta[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] - delta[ 7 ][ 2 ] },
                { grad_phi[ 8 ][ 0 ] - delta[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ] - delta[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] - delta[ 8 ][ 2 ] }
        };


        tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, positive_perturb,
                                                                                      FP, chiP, gradChiP );

        tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( grad_u, phi, negative_perturb,
                                                                                      FM, chiM, gradChiM );

        variableVector gradCol = ( FP - FM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( chiP - chiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], 0. ) );
        }

        gradCol = ( gradChiP - gradChiM ) / ( 2 * delta[ ii ][ ij ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dGradChidGradPhi[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( test_generate_input_variable_string ){

    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                        {0.100, 0.001, 0.000 },
                                        {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                          { -0.18572739,  0.06847269,  0.22931628 },
                                          { -0.01829735, -0.48731265, -0.25277529 },
                                          {  0.26626212,  0.4844646 , -0.31965177 },
                                          {  0.49197846,  0.19051656, -0.0365349  },
                                          { -0.06607774, -0.33526875, -0.15803078 },
                                          {  0.09738707, -0.49482218, -0.39584868 },
                                          { -0.45599864,  0.08585038, -0.09432794 },
                                          {  0.23055539,  0.07564162,  0.24051469 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };


    //Initialize the state variable vector
    std::vector< double > SDVS( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    std::string result;

    std::string answer = "time:\n 1.000000e+01, 2.500000e+00,\n\nfparams:\n 2.000000e+00, 2.400000e+02, 1.500000e+01, 2.000000e+00, 1.400000e+02, 2.000000e+01, 2.000000e+00, 2.000000e+00, 2.700000e+01, 2.000000e+00, 5.600000e-01, 2.000000e-01, 2.000000e+00, 1.500000e-01, -2.000000e-01, 2.000000e+00, 8.200000e-01, 1.000000e-01, 2.000000e+00, 7.000000e-01, 3.000000e-01, 2.000000e+00, 4.000000e-01, -3.000000e-01, 2.000000e+00, 5.200000e-01, 4.000000e-01, 2.000000e+00, 6.964700e+02, 6.584000e+01, 5.000000e+00, -7.690000e+00, -5.192000e+01, 3.861000e+01, -2.731000e+01, 5.130000e+00, 1.100000e+01, 1.850000e+00, -1.900000e-01, -1.080000e+00, -1.570000e+00, 2.290000e+00, -6.100000e-01, 5.970000e+00, -2.020000e+00, 2.380000e+00, -3.200000e-01, -3.250000e+00, 2.000000e+00, -5.192000e+01, 5.130000e+00, 4.000000e-01, 3.000000e-01, 3.500000e-01, 1.000000e-08, 1.000000e-08,\n\ncurrent_grad_u:\n 2.000000e-01, 1.000000e-01, 0.000000e+00,\n 1.000000e-01, 1.000000e-03, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n\ncurrent_phi:\n 1.000000e-01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,\n\ncurrent_grad_phi:\n 1.389002e-01, -3.598602e-01, -8.048856e-02,\n -1.857274e-01, 6.847269e-02, 2.293163e-01,\n -1.829735e-02, -4.873127e-01, -2.527753e-01,\n 2.662621e-01, 4.844646e-01, -3.196518e-01,\n 4.919785e-01, 1.905166e-01, -3.653490e-02,\n -6.607774e-02, -3.352688e-01, -1.580308e-01,\n 9.738707e-02, -4.948222e-01, -3.958487e-01,\n -4.559986e-01, 8.585038e-02, -9.432794e-02,\n 2.305554e-01, 7.564162e-02, 2.405147e-01,\n\nprevious_grad_u:\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n\nprevious_phi:\n 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,\n\nprevious_grad_phi:\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n 0.000000e+00, 0.000000e+00, 0.000000e+00,\n\nSDVS:\n 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,";
    answer += " 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,\n\ncurrent_ADD_DOF:\n\ncurrent_ADD_grad_DOF:\n\nprevious_ADD_DOF:\n\nprevious_ADD_grad_DOF:\n";
//    std::string answer = "time:\n 10.000000, 2.500000,\nfparams:\n 2.000000, 240.000000, 15.000000, 2.000000, 140.000000, 20.000000, 2.000000, 2.000000, 27.000000, 2.000000, 0.560000, 0.200000, 2.000000, 0.150000, -0.200000, 2.000000, 0.820000, 0.100000, 2.000000, 0.700000, 0.300000, 2.000000, 0.400000, -0.300000, 2.000000, 0.520000, 0.400000, 2.000000, 696.470000, 65.840000, 5.000000, -7.690000, -51.920000, 38.610000, -27.310000, 5.130000, 11.000000, 1.850000, -0.190000, -1.080000, -1.570000, 2.290000, -0.610000, 5.970000, -2.020000, 2.380000, -0.320000, -3.250000, 2.000000, -51.920000, 5.130000, 0.400000, 0.300000, 0.350000, 0.000000, 0.000000,\ncurrent_grad_u:\n 0.200000, 0.100000, 0.000000,\n 0.100000, 0.001000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n\ncurrent_phi:\n 0.100000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,\ncurrent_grad_phi:\n 0.138900, -0.359860, -0.080489,\n -0.185727, 0.068473, 0.229316,\n -0.018297, -0.487313, -0.252775,\n 0.266262, 0.484465, -0.319652,\n 0.491978, 0.190517, -0.036535,\n -0.066078, -0.335269, -0.158031,\n 0.097387, -0.494822, -0.395849,\n -0.455999, 0.085850, -0.094328,\n 0.230555, 0.075642, 0.240515,\n\nprevious_grad_u:\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n\nprevious_phi:\n 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,\nprevious_grad_phi:\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n 0.000000, 0.000000, 0.000000,\n\nSDVS:\n 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,\ncurrent_ADD_DOF:\n\ncurrent_ADD_grad_DOF:\n\nprevious_ADD_DOF:\n\nprevious_ADD_grad_DOF:\n";

    tardigradeMicromorphicElastoPlasticity::generate_input_variable_string( time, fparams, current_grad_u, current_phi, current_grad_phi,
                                                                            previous_grad_u, previous_phi, previous_grad_phi,
                                                                            SDVS, current_ADD_DOF, current_ADD_grad_DOF,
                                                                            previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                            result );

    BOOST_CHECK( answer.compare( result ) == 0 );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                        {0.100, 0.001, 0.000 },
                                        {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
//                                0.44195587,  0.34121966, -0.79098944, 
//                               -0.43965428,  0.88466225,  0.1684519 };
//
//    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
//                                 -0.11111872, -0.07416114, -1.01048108,
//                                  0.1804018 , -1.01116291,  0.03248007 };

    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                          { -0.18572739,  0.06847269,  0.22931628 },
                                          { -0.01829735, -0.48731265, -0.25277529 },
                                          {  0.26626212,  0.4844646 , -0.31965177 },
                                          {  0.49197846,  0.19051656, -0.0365349  },
                                          { -0.06607774, -0.33526875, -0.15803078 },
                                          {  0.09738707, -0.49482218, -0.39584868 },
                                          { -0.45599864,  0.08585038, -0.09432794 },
                                          {  0.23055539,  0.07564162,  0.24051469 } };

//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > PK2_result( 9, 0 );

    std::vector< double > SIGMA_result( 9, 0 );

    std::vector< double > M_result( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    tardigradeSolverTools::homotopyMap homotopyDEBUG;
#endif

    tardigradeSolverTools::floatVector PK2_answer = { 1.72374955e+02,  1.53586446e+01, -9.18594711e-01,  1.34613772e+01,
        1.42758448e+02, -2.14353350e-02, -1.76155689e+00,  1.77580288e+00,
        1.41003910e+02 };

    tardigradeSolverTools::floatVector SIGMA_answer = { 176.85221932,  15.84248153,  -2.83823499,  15.84248153,
       144.5183467 ,   1.85778345,  -2.83823499,   1.85778345,
       141.99307014 };

    tardigradeSolverTools::floatVector M_answer = { 0.60016978, -0.51045338,  0.61980289,  3.23507063,  1.16925208,
        1.20665651,  0.56034359, -2.52042105,  1.62648849, -2.61882314,
       -0.61143851, -1.02227749,  0.67046998,  0.49701267, -0.23999964,
       -2.77670511,  0.75636495,  1.71897722, -0.49808019,  2.62569695,
       -0.76002078,  1.23462488, -0.00650376, -2.25591243, -0.73016414,
        0.74380723,  0.90861263 };

    tardigradeSolverTools::floatVector SDVS_answer = {  7.52790780e-03,  5.06936872e-03, -3.00562092e-04,  5.05251313e-03,
                                                       -2.63763651e-03,  2.87928507e-04, -3.72441969e-04,  2.25230220e-04,
                                                       -1.58168041e-03,  7.52714199e-03,  4.60154810e-03, -3.01710599e-04,
                                                        5.56787320e-03, -2.63687070e-03,  2.58160226e-04, -3.65069825e-04,
                                                        2.58160226e-04, -1.58168041e-03,  3.83196420e-02, -2.45387739e-02,
                                                       -8.94156844e-03,  1.04125027e-02, -5.47145710e-04,  2.21268303e-02,
                                                        1.85810005e-02, -1.14241984e-02, -7.62221572e-03,  4.54399478e-02,
                                                        3.21567189e-02,  2.10892064e-02,  1.18343147e-02,  1.63154826e-02,
                                                        1.01776270e-03,  1.47855827e-02,  7.10955065e-03, -2.40560608e-02,
                                                        9.30082927e-03, -3.52878321e-02,  2.18504370e-02, -2.49577455e-02,
                                                        9.18791192e-03,  2.13410715e-02,  1.69866062e-02,  1.77613753e-02,
                                                        2.03937074e-02, -5.42245900e-23,  5.17942700e-03,  3.02931048e-02,
                                                        1.58261418e-02,  1.60309023e-02, -1.79683900e-22,  2.11178609e-02,
                                                        8.23638353e-02,  4.30296627e-02,  4.35864383e-02 };

//    tardigradeSolverTools::floatVector PK2_answer = { 1.72381583e+02,  1.53490623e+01, -9.17869152e-01,  1.34514475e+01,
//        1.42784323e+02, -2.26855750e-02, -1.76082294e+00,  1.77464785e+00,
//        1.41021846e+02 };
//
//    tardigradeSolverTools::floatVector SIGMA_answer = { 176.85984866,  15.83223893,  -2.83753359,  15.83223893,
//       144.5458737 ,   1.85655805,  -2.83753359,   1.85655805,
//       142.01240974 };
//
//    tardigradeSolverTools::floatVector M_answer = { 0.60016027, -0.51048754,  0.61984844,  3.23507227,  1.16942361,
//        1.20662324,  0.56038461, -2.52062122,  1.6263198 , -2.61891654,
//       -0.61168873, -1.02203605,  0.67065387,  0.4970073 , -0.23999741,
//       -2.7769288 ,  0.75667767,  1.71898866, -0.49803026,  2.62584805,
//       -0.75966908,  1.23459271, -0.00667791, -2.25598899, -0.73030698,
//        0.7438434 ,  0.90872187 };
//
//    tardigradeSolverTools::floatVector SDVS_answer = {  7.53391712e-03,  5.08826292e-03, -3.02187662e-04,  5.07135003e-03,
//                                                       -2.66948991e-03,  2.90841417e-04, -3.74157919e-04,  2.27820300e-04,
//                                                       -1.59391881e-03,  7.53314870e-03,  4.61870100e-03, -3.03219139e-04,
//                                                        5.58862820e-03, -2.66872150e-03,  2.60913496e-04, -3.66895165e-04,
//                                                        2.60913496e-04, -1.59391881e-03,  3.83381848e-02, -2.45330931e-02,
//                                                       -8.95538204e-03,  1.04165665e-02, -5.32012910e-04,  2.21443208e-02,
//                                                        1.85774208e-02, -1.14594874e-02, -7.63911022e-03,  4.54127235e-02,
//                                                        3.21129124e-02,  2.11070183e-02,  1.18241526e-02,  1.63088860e-02,
//                                                        1.03455560e-03,  1.47830982e-02,  7.09499615e-03, -2.40626243e-02,
//                                                        9.31063821e-03, -3.52663045e-02,  2.18795269e-02, -2.49636009e-02,
//                                                        9.20397912e-03,  2.13550248e-02,  1.69839898e-02,  1.77600107e-02,
//                                                        2.03926694e-02,  3.09063946e-16,  5.21174400e-03,  3.02958646e-02,
//                                                        1.58295896e-02,  1.60303582e-02,  1.03245982e-15,  2.12496399e-02,
//                                                        8.23712094e-02,  4.30389579e-02,  4.35848261e-02 };

    std::vector< double > SDVS = SDVSDefault;

    int errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
                                                                                  current_grad_u,  current_phi,  current_grad_phi,
                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                                  SDVS,
                                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                                  PK2_result, SIGMA_result, M_result,
                                                                                  ADD_TERMS,
                                                                                  output_message
                                                                                  );

    BOOST_CHECK( errorCode == 0 );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SDVS_answer,  SDVS ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_answer,   PK2_result ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_answer, SIGMA_result ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_answer,     M_result ) );

    //Test the Jacobians
    PK2_result.clear();
    SIGMA_result.clear();
    M_result.clear();
    ADD_TERMS.clear();

    variableMatrix result_dPK2dGradU(      9, variableVector(  9, 0 ) );

    variableMatrix result_dPK2dPhi(        9, variableVector(  9, 0 ) );

    variableMatrix result_dPK2dGradPhi(    9, variableVector( 27, 0 ) );

    variableMatrix result_dSIGMAdGradU(    9, variableVector(  9, 0 ) );

    variableMatrix result_dSIGMAdPhi(      9, variableVector(  9, 0 ) );

    variableMatrix result_dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );

    variableMatrix result_dMdGradU(       27, variableVector(  9, 0 ) );

    variableMatrix result_dMdPhi(         27, variableVector(  9, 0 ) );

    variableMatrix result_dMdGradPhi(     27, variableVector( 27, 0 ) );

    std::vector< variableMatrix > ADD_JACOBIANS;

    SDVS = SDVSDefault;

    errorCode  = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
                                                                               current_grad_u,  current_phi,  current_grad_phi,
                                                                               previous_grad_u, previous_phi, previous_grad_phi,
                                                                               SDVS,
                                                                               current_ADD_DOF,  current_ADD_grad_DOF,
                                                                               previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                               PK2_result, SIGMA_result, M_result,
                                                                               result_dPK2dGradU, result_dPK2dPhi, result_dPK2dGradPhi,
                                                                               result_dSIGMAdGradU, result_dSIGMAdPhi, result_dSIGMAdGradPhi,
                                                                               result_dMdGradU, result_dMdPhi, result_dMdGradPhi,
                                                                               ADD_TERMS, ADD_JACOBIANS, output_message
                                                                             );

    BOOST_CHECK( errorCode <= 0 );

    if ( errorCode != 0 ){
        std::cout << "output_message:\n" << output_message << "\n";
    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_result, PK2_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_result, M_answer ) );

    variableMatrix dPK2dGradU(      9, variableVector(  9, 0 ) );

    variableMatrix dPK2dPhi(        9, variableVector(  9, 0 ) );

    variableMatrix dPK2dGradPhi(    9, variableVector( 27, 0 ) );

    variableMatrix dSIGMAdGradU(    9, variableVector(  9, 0 ) );

    variableMatrix dSIGMAdPhi(      9, variableVector(  9, 0 ) );

    variableMatrix dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );

    variableMatrix dMdGradU(       27, variableVector(  9, 0 ) );

    variableMatrix dMdPhi(         27, variableVector(  9, 0 ) );

    variableMatrix dMdGradPhi(     27, variableVector( 27, 0 ) );

    variableMatrix dFpdGradU(          9, variableVector(  9, 0 ) );

    variableMatrix dFpdPhi(            9, variableVector(  9, 0 ) );

    variableMatrix dFpdGradPhi(        9, variableVector( 27, 0 ) );

    variableMatrix dChipdGradU(        9, variableVector(  9, 0 ) );

    variableMatrix dChipdPhi(          9, variableVector(  9, 0 ) );

    variableMatrix dChipdGradPhi(      9, variableVector( 27, 0 ) );

    variableMatrix dGradChipdGradU(   27, variableVector(  9, 0 ) );

    variableMatrix dGradChipdPhi(     27, variableVector(  9, 0 ) );

    variableMatrix dGradChipdGradPhi( 27, variableVector( 27, 0 ) );

    variableType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        variableVector delta( 9, 0 );

        unsigned int row = i / 3;

        unsigned int col = i % 3;

        delta[ i ] = eps * std::fabs( current_grad_u[ row ][ col ] ) + eps;

        variableType current_grad_u_p[ 3 ][ 3 ];
        variableType current_grad_u_m[ 3 ][ 3 ];

        for ( unsigned int _i = 0; _i < 3; _i++ ){
            for ( unsigned int _j = 0; _j < 3; _j++ ){
                current_grad_u_p[ _i ][ _j ] = current_grad_u[ _i ][ _j ] + delta[ 3 * _i + _j ];
                current_grad_u_m[ _i ][ _j ] = current_grad_u[ _i ][ _j ] - delta[ 3 * _i + _j ];
            }
        }

        variableVector PK2_p,   PK2_m;
        variableVector SIGMA_p, SIGMA_m;
        variableVector M_p,     M_m;
        variableVector SDVS_p = SDVSDefault;
        variableVector SDVS_m = SDVSDefault;

        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_p, current_phi, current_grad_phi,
                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                                  PK2_p, SIGMA_p, M_p,
                                                                                  ADD_TERMS, output_message
                                                                                );

        BOOST_CHECK( errorCode <= 0 );

        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_m, current_phi, current_grad_phi,
                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                                  PK2_m, SIGMA_m, M_m,
                                                                                  ADD_TERMS, output_message
                                                                                );

        BOOST_CHECK( errorCode <= 0 );

        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){

            dPK2dGradU[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){

            dSIGMAdGradU[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < M_p.size( ); j++ ){

            dMdGradU[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            dFpdGradU[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            dChipdGradU[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            dGradChipdGradU[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradU,   result_dPK2dGradU ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradU, result_dSIGMAdGradU ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradU,     result_dMdGradU ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradU,       ADD_JACOBIANS[ 0 ] ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradU,     ADD_JACOBIANS[ 3 ] ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradU, ADD_JACOBIANS[ 6 ] ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        variableVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( current_phi[ i ] ) + eps;

        variableType current_phi_p[ 9 ];
        variableType current_phi_m[ 9 ];

        for ( unsigned int _i = 0; _i < 3; _i++ ){
            for ( unsigned int _j = 0; _j < 3; _j++ ){
                current_phi_p[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] + delta[ 3 * _i + _j ];
                current_phi_m[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] - delta[ 3 * _i + _j ];
            }
        }

        variableVector PK2_p,   PK2_m;
        variableVector SIGMA_p, SIGMA_m;
        variableVector M_p,     M_m;
        variableVector SDVS_p = SDVSDefault;
        variableVector SDVS_m = SDVSDefault;

        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_p, current_grad_phi,
                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                                  PK2_p, SIGMA_p, M_p,
                                                                                  ADD_TERMS, output_message
                                                                                );

        BOOST_CHECK( errorCode <= 0 );

        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_m, current_grad_phi,
                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                                  PK2_m, SIGMA_m, M_m,
                                                                                  ADD_TERMS, output_message
                                                                                );

        BOOST_CHECK( errorCode <= 0 );

        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){

            dPK2dPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){

            dSIGMAdPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < M_p.size( ); j++ ){

            dMdPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            dFpdPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            dChipdPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            dGradChipdPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dPhi,   result_dPK2dPhi ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdPhi, result_dSIGMAdPhi ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdPhi,     result_dMdPhi ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdPhi,       ADD_JACOBIANS[ 1 ] ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdPhi,     ADD_JACOBIANS[ 4 ] ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdPhi, ADD_JACOBIANS[ 7 ] ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        variableVector delta( 27, 0 );

        unsigned int row = i / 3;

        unsigned int col = i - 3 * row;

        delta[ i ] = eps * std::fabs( current_grad_phi[ row ][ col ] ) + eps;

        variableType current_grad_phi_p[ 9 ][ 3 ];
        variableType current_grad_phi_m[ 9 ][ 3 ];

        for ( unsigned int _i = 0; _i < 9; _i++ ){
            for ( unsigned int _j = 0; _j < 3; _j++ ){
                current_grad_phi_p[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] + delta[ 3 * _i + _j ];
                current_grad_phi_m[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] - delta[ 3 * _i + _j ];
            }
        }

        variableVector PK2_p,   PK2_m;
        variableVector SIGMA_p, SIGMA_m;
        variableVector M_p,     M_m;
        variableVector SDVS_p = SDVSDefault;
        variableVector SDVS_m = SDVSDefault;

        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_p,
                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                                  PK2_p, SIGMA_p, M_p,
                                                                                  ADD_TERMS, output_message
                                                                                );

        BOOST_CHECK( errorCode <= 0 );

        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_m,
                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                                                  PK2_m, SIGMA_m, M_m,
                                                                                  ADD_TERMS, output_message
                                                                                );

        BOOST_CHECK( errorCode <= 0 );

        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){

            dPK2dGradPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){

            dSIGMAdGradPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < M_p.size( ); j++ ){

            dMdGradPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            dFpdGradPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            dChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            dGradChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradPhi,   result_dPK2dGradPhi ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradPhi, result_dSIGMAdGradPhi ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradPhi,     result_dMdGradPhi ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradPhi,       ADD_JACOBIANS[ 2 ] ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradPhi,     ADD_JACOBIANS[ 5 ] ) );
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradPhi, ADD_JACOBIANS[ 8 ] ) );

}

BOOST_AUTO_TEST_CASE( testParameterExtraction){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    double temperature = 293.15;

    double previousTemperature = 293.15;

    //Initialize the material parameters
    std::vector< double > fparams = { 4, 2.4e2, 1.5e1, 1e-2, 0.12, //Macro hardening parameters
                                      4, 1.4e2, 2.0e1, 2e-2, 0.13, //Micro hardening parameters
                                      4, 2.0e0, 2.7e1, 3e-2, 0.14, //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                        {0.100, 0.001, 0.000 },
                                        {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
//                                0.44195587,  0.34121966, -0.79098944, 
//                               -0.43965428,  0.88466225,  0.1684519 };
//
//    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
//                                 -0.11111872, -0.07416114, -1.01048108,
//                                  0.1804018 , -1.01116291,  0.03248007 };

    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                          { -0.18572739,  0.06847269,  0.22931628 },
                                          { -0.01829735, -0.48731265, -0.25277529 },
                                          {  0.26626212,  0.4844646 , -0.31965177 },
                                          {  0.49197846,  0.19051656, -0.0365349  },
                                          { -0.06607774, -0.33526875, -0.15803078 },
                                          {  0.09738707, -0.49482218, -0.39584868 },
                                          { -0.45599864,  0.08585038, -0.09432794 },
                                          {  0.23055539,  0.07564162,  0.24051469 } };

//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > PK2_result( 9, 0 );

    std::vector< double > SIGMA_result( 9, 0 );

    std::vector< double > M_result( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

#ifdef DEBUG_MODE
    tardigradeSolverTools::homotopyMap homotopyDEBUG;
#endif

    tardigradeSolverTools::floatVector PK2_answer = { 1.72374955e+02,  1.53586446e+01, -9.18594711e-01,  1.34613772e+01,
        1.42758448e+02, -2.14353350e-02, -1.76155689e+00,  1.77580288e+00,
        1.41003910e+02 };

    tardigradeSolverTools::floatVector SIGMA_answer = { 176.85221932,  15.84248153,  -2.83823499,  15.84248153,
       144.5183467 ,   1.85778345,  -2.83823499,   1.85778345,
       141.99307014 };

    tardigradeSolverTools::floatVector M_answer = { 0.60016978, -0.51045338,  0.61980289,  3.23507063,  1.16925208,
        1.20665651,  0.56034359, -2.52042105,  1.62648849, -2.61882314,
       -0.61143851, -1.02227749,  0.67046998,  0.49701267, -0.23999964,
       -2.77670511,  0.75636495,  1.71897722, -0.49808019,  2.62569695,
       -0.76002078,  1.23462488, -0.00650376, -2.25591243, -0.73016414,
        0.74380723,  0.90861263 };

    tardigradeSolverTools::floatVector SDVS_answer = {  7.52790780e-03,  5.06936872e-03, -3.00562092e-04,  5.05251313e-03,
                                                       -2.63763651e-03,  2.87928507e-04, -3.72441969e-04,  2.25230220e-04,
                                                       -1.58168041e-03,  7.52714199e-03,  4.60154810e-03, -3.01710599e-04,
                                                        5.56787320e-03, -2.63687070e-03,  2.58160226e-04, -3.65069825e-04,
                                                        2.58160226e-04, -1.58168041e-03,  3.83196420e-02, -2.45387739e-02,
                                                       -8.94156844e-03,  1.04125027e-02, -5.47145710e-04,  2.21268303e-02,
                                                        1.85810005e-02, -1.14241984e-02, -7.62221572e-03,  4.54399478e-02,
                                                        3.21567189e-02,  2.10892064e-02,  1.18343147e-02,  1.63154826e-02,
                                                        1.01776270e-03,  1.47855827e-02,  7.10955065e-03, -2.40560608e-02,
                                                        9.30082927e-03, -3.52878321e-02,  2.18504370e-02, -2.49577455e-02,
                                                        9.18791192e-03,  2.13410715e-02,  1.69866062e-02,  1.77613753e-02,
                                                        2.03937074e-02, -5.42245900e-23,  5.17942700e-03,  3.02931048e-02,
                                                        1.58261418e-02,  1.60309023e-02, -1.79683900e-22,  2.11178609e-02,
                                                        8.23638353e-02,  4.30296627e-02,  4.35864383e-02 };

//    tardigradeSolverTools::floatVector PK2_answer = { 1.72381583e+02,  1.53490623e+01, -9.17869152e-01,  1.34514475e+01,
//        1.42784323e+02, -2.26855750e-02, -1.76082294e+00,  1.77464785e+00,
//        1.41021846e+02 };
//
//    tardigradeSolverTools::floatVector SIGMA_answer = { 176.85984866,  15.83223893,  -2.83753359,  15.83223893,
//       144.5458737 ,   1.85655805,  -2.83753359,   1.85655805,
//       142.01240974 };
//
//    tardigradeSolverTools::floatVector M_answer = { 0.60016027, -0.51048754,  0.61984844,  3.23507227,  1.16942361,
//        1.20662324,  0.56038461, -2.52062122,  1.6263198 , -2.61891654,
//       -0.61168873, -1.02203605,  0.67065387,  0.4970073 , -0.23999741,
//       -2.7769288 ,  0.75667767,  1.71898866, -0.49803026,  2.62584805,
//       -0.75966908,  1.23459271, -0.00667791, -2.25598899, -0.73030698,
//        0.7438434 ,  0.90872187 };
//
//    tardigradeSolverTools::floatVector SDVS_answer = {  7.53391712e-03,  5.08826292e-03, -3.02187662e-04,  5.07135003e-03,
//                                                       -2.66948991e-03,  2.90841417e-04, -3.74157919e-04,  2.27820300e-04,
//                                                       -1.59391881e-03,  7.53314870e-03,  4.61870100e-03, -3.03219139e-04,
//                                                        5.58862820e-03, -2.66872150e-03,  2.60913496e-04, -3.66895165e-04,
//                                                        2.60913496e-04, -1.59391881e-03,  3.83381848e-02, -2.45330931e-02,
//                                                       -8.95538204e-03,  1.04165665e-02, -5.32012910e-04,  2.21443208e-02,
//                                                        1.85774208e-02, -1.14594874e-02, -7.63911022e-03,  4.54127235e-02,
//                                                        3.21129124e-02,  2.11070183e-02,  1.18241526e-02,  1.63088860e-02,
//                                                        1.03455560e-03,  1.47830982e-02,  7.09499615e-03, -2.40626243e-02,
//                                                        9.31063821e-03, -3.52663045e-02,  2.18795269e-02, -2.49636009e-02,
//                                                        9.20397912e-03,  2.13550248e-02,  1.69839898e-02,  1.77600107e-02,
//                                                        2.03926694e-02,  3.09063946e-16,  5.21174400e-03,  3.02958646e-02,
//                                                        1.58295896e-02,  1.60303582e-02,  1.03245982e-15,  2.12496399e-02,
//                                                        8.23712094e-02,  4.30389579e-02,  4.35848261e-02 };

    std::vector< double > SDVS = SDVSDefault;

    tardigradeSolverTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeSolverTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                                                    currentDeformationGradient, currentMicroDeformation,
                                                                                    currentGradientMicroDeformation );

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                                                    previousDeformationGradient, previousMicroDeformation,
                                                                                    previousGradientMicroDeformation );

    tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticity hydra( time[ 0 ], time[ 1 ],
                                                                                     temperature,                     previousTemperature,
                                                                                     currentDeformationGradient,      previousDeformationGradient,
                                                                                     currentMicroDeformation,         previousMicroDeformation,
                                                                                     currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                                                     { }, { },
                                                                                     SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );

    tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticityOptimization hydraOpt( time[ 0 ], time[ 1 ],
                                                                                                    temperature,                     previousTemperature,
                                                                                                    currentDeformationGradient,      previousDeformationGradient,
                                                                                                    currentMicroDeformation,         previousMicroDeformation,
                                                                                                    currentGradientMicroDeformation, previousGradientMicroDeformation,
                                                                                                    { }, { },
                                                                                                    SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );

    tardigradeSolverTools::floatVector plasticParameters( fparams.begin( ), fparams.begin( ) + 33 );

    tardigradeSolverTools::floatVector elasticParameters( fparams.begin( ) + 33, fparams.begin( ) + 57 );

    BOOST_TEST( hydra.getNumPlasticParameters( ) == 33 );

    BOOST_TEST( hydra.getPlasticParameters( ) == plasticParameters, CHECK_PER_ELEMENT );

    BOOST_TEST( hydra.getElasticParameters( ) == elasticParameters, CHECK_PER_ELEMENT );

    BOOST_TEST( hydraOpt.getNumPlasticParameters( ) == 33 );

    BOOST_TEST( hydraOpt.getPlasticParameters( ) == plasticParameters, CHECK_PER_ELEMENT );

    BOOST_TEST( hydraOpt.getElasticParameters( ) == elasticParameters, CHECK_PER_ELEMENT );

}
