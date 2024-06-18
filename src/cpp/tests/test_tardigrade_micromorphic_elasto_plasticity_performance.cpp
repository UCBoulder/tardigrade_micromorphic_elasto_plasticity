//Tests for tardigrade_constitutive_tools

#include<tardigrade_micromorphic_elasto_plasticity.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>

#define BOOST_TEST_MODULE test_tardigrade_micromorphic_elasto_plasticity
#include <boost/test/included/unit_test.hpp>

typedef tardigradeMicromorphicTools::constantType constantType;
typedef tardigradeMicromorphicTools::constantVector constantVector;
typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

typedef tardigradeMicromorphicTools::parameterType parameterType;
typedef tardigradeMicromorphicTools::parameterVector parameterVector;
typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

typedef tardigradeMicromorphicTools::variableType variableType;
typedef tardigradeMicromorphicTools::variableVector variableVector;
typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

typedef tardigradeMicromorphicTools::errorNode errorNode;
typedef tardigradeMicromorphicTools::errorOut errorOut;

//BOOST_AUTO_TEST_CASE( testEvaluateHydraModel){
//    /*!
//     * Test the evaluation of the constitutive model.
//     *
//     */
//
//    //Initialize the time
//    std::vector< double > time = { 10., 2.5 };
//
//    //Initialize the material parameters
//    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
//                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
//                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
//                                      2, 0.56, 0.2,                //Macro flow parameters
//                                      2, 0.15,-0.2,                //Micro flow parameters
//                                      2, 0.82, 0.1,                //Micro gradient flow parameters
//                                      2, 0.70, 0.3,                //Macro yield parameters
//                                      2, 0.40,-0.3,                //Micro yield parameters
//                                      2, 0.52, 0.4,                //Micro gradient yield parameters
//                                      2, 696.47, 65.84,            //A stiffness tensor parameters
//                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
//                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
//                                      2, -51.92, 5.13,             //D stiffness tensor parameters
//                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
//                                    };
//
//    //Initialize the gradient of the macro displacement
////    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
////                                        { -0.56632666, -0.21399259,  0.16367238 },
////                                        { -0.29129789, -0.22367825, -2.0632945  } };
////
////    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
////                                         {  0.31303067, -1.23910631, -0.93837662 },
////                                         { -0.32571524, -0.95306342, -0.93025257 } };
//
//    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
//                                        {0.100, 0.001, 0.000 },
//                                        {0.000, 0.000, 0.000 } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
//                                         {0, 0, 0},
//                                         {0, 0, 0} };
//    //Initialize the micro displacement
////    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
////                                0.44195587,  0.34121966, -0.79098944, 
////                               -0.43965428,  0.88466225,  0.1684519 };
////
////    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
////                                 -0.11111872, -0.07416114, -1.01048108,
////                                  0.1804018 , -1.01116291,  0.03248007 };
//
//    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
//                                0.000, 0.000, 0.000,
//                                0.000, 0.000, 0.000 };
//
//    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//
//    //Initialize the gradient of the micro displacement
//    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
//                                          { -0.18572739,  0.06847269,  0.22931628 },
//                                          { -0.01829735, -0.48731265, -0.25277529 },
//                                          {  0.26626212,  0.4844646 , -0.31965177 },
//                                          {  0.49197846,  0.19051656, -0.0365349  },
//                                          { -0.06607774, -0.33526875, -0.15803078 },
//                                          {  0.09738707, -0.49482218, -0.39584868 },
//                                          { -0.45599864,  0.08585038, -0.09432794 },
//                                          {  0.23055539,  0.07564162,  0.24051469 } };
//
////    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
////                                           {  0.00294417,  0.34480654, -0.34450988 },
////                                           {  0.21056511, -0.28113967, -0.45726839 },
////                                           { -0.26431286, -0.09985721,  0.47322301 },
////                                           { -0.18156887, -0.32226199, -0.37295847 },
////                                           {  0.15062371,  0.09439471,  0.09167948 },
////                                           { -0.46869859,  0.018301  ,  0.45013866 },
////                                           { -0.15455446,  0.40552715, -0.4216042  },
////                                           { -0.38930237,  0.10974753, -0.31188239 } };
//
////    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0},
////                                          {0, 0, 0} };
//
//    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0},
//                                           {0, 0, 0} };
//                                           
//
//    //Initialize the state variable vector
//    std::vector< double > SDVSDefault( 55, 0 );
//
//    //Initialize the additional degree of freedom vectors
//    std::vector< double > current_ADD_DOF;
//    std::vector< std::vector< double > > current_ADD_grad_DOF;
//
//    std::vector< double > previous_ADD_DOF;
//    std::vector< std::vector< double > > previous_ADD_grad_DOF;
//
//    //Initialize the stress measures
//    std::vector< double > PK2_result( 9, 0 );
//
//    std::vector< double > SIGMA_result( 9, 0 );
//
//    std::vector< double > M_result( 27, 0 );
//
//    //Initialize the additional terms vector
//    std::vector< std::vector< double > > ADD_TERMS;
//
//    //Initialize the output message string
//    std::string output_message;
//
//#ifdef DEBUG_MODE
//    tardigradeSolverTools::homotopyMap homotopyDEBUG;
//#endif
//
//    tardigradeSolverTools::floatVector PK2_answer = { 1.72374955e+02,  1.53586446e+01, -9.18594711e-01,  1.34613772e+01,
//        1.42758448e+02, -2.14353350e-02, -1.76155689e+00,  1.77580288e+00,
//        1.41003910e+02 };
//
//    tardigradeSolverTools::floatVector SIGMA_answer = { 176.85221932,  15.84248153,  -2.83823499,  15.84248153,
//       144.5183467 ,   1.85778345,  -2.83823499,   1.85778345,
//       141.99307014 };
//
//    tardigradeSolverTools::floatVector M_answer = { 0.60016978, -0.51045338,  0.61980289,  3.23507063,  1.16925208,
//        1.20665651,  0.56034359, -2.52042105,  1.62648849, -2.61882314,
//       -0.61143851, -1.02227749,  0.67046998,  0.49701267, -0.23999964,
//       -2.77670511,  0.75636495,  1.71897722, -0.49808019,  2.62569695,
//       -0.76002078,  1.23462488, -0.00650376, -2.25591243, -0.73016414,
//        0.74380723,  0.90861263 };
//
//    tardigradeSolverTools::floatVector SDVS_answer = {  7.52790780e-03,  5.06936872e-03, -3.00562092e-04,  5.05251313e-03,
//                                                       -2.63763651e-03,  2.87928507e-04, -3.72441969e-04,  2.25230220e-04,
//                                                       -1.58168041e-03,  7.52714199e-03,  4.60154810e-03, -3.01710599e-04,
//                                                        5.56787320e-03, -2.63687070e-03,  2.58160226e-04, -3.65069825e-04,
//                                                        2.58160226e-04, -1.58168041e-03,  3.83196420e-02, -2.45387739e-02,
//                                                       -8.94156844e-03,  1.04125027e-02, -5.47145710e-04,  2.21268303e-02,
//                                                        1.85810005e-02, -1.14241984e-02, -7.62221572e-03,  4.54399478e-02,
//                                                        3.21567189e-02,  2.10892064e-02,  1.18343147e-02,  1.63154826e-02,
//                                                        1.01776270e-03,  1.47855827e-02,  7.10955065e-03, -2.40560608e-02,
//                                                        9.30082927e-03, -3.52878321e-02,  2.18504370e-02, -2.49577455e-02,
//                                                        9.18791192e-03,  2.13410715e-02,  1.69866062e-02,  1.77613753e-02,
//                                                        2.03937074e-02, -5.42245900e-23,  5.17942700e-03,  3.02931048e-02,
//                                                        1.58261418e-02,  1.60309023e-02, -1.79683900e-22,  2.11178609e-02,
//                                                        8.23638353e-02,  4.30296627e-02,  4.35864383e-02 };
//
//    std::vector< double > SDVS = SDVSDefault;
//
//    int errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
//                                                                                  current_grad_u,  current_phi,  current_grad_phi,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS,
//                                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_result, SIGMA_result, M_result,
//                                                                                  ADD_TERMS,
//                                                                                  output_message
//                                                                                  );
//
//    BOOST_CHECK( errorCode == 0 );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SDVS_answer,  SDVS ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_answer,   PK2_result ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_answer, SIGMA_result ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_answer,     M_result ) );
//
//    //Test the Jacobians
//    PK2_result.clear();
//    SIGMA_result.clear();
//    M_result.clear();
//    ADD_TERMS.clear();
//
//    variableMatrix result_dPK2dGradU(      9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dPK2dPhi(        9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dPK2dGradPhi(    9, variableVector( 27, 0 ) );
//
//    variableMatrix result_dSIGMAdGradU(    9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dSIGMAdPhi(      9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );
//
//    variableMatrix result_dMdGradU(       27, variableVector(  9, 0 ) );
//
//    variableMatrix result_dMdPhi(         27, variableVector(  9, 0 ) );
//
//    variableMatrix result_dMdGradPhi(     27, variableVector( 27, 0 ) );
//
//    std::vector< variableMatrix > ADD_JACOBIANS;
//
//    SDVS = SDVSDefault;
//
//    errorCode  = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
//                                                                               current_grad_u,  current_phi,  current_grad_phi,
//                                                                               previous_grad_u, previous_phi, previous_grad_phi,
//                                                                               SDVS,
//                                                                               current_ADD_DOF,  current_ADD_grad_DOF,
//                                                                               previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                               PK2_result, SIGMA_result, M_result,
//                                                                               result_dPK2dGradU, result_dPK2dPhi, result_dPK2dGradPhi,
//                                                                               result_dSIGMAdGradU, result_dSIGMAdPhi, result_dSIGMAdGradPhi,
//                                                                               result_dMdGradU, result_dMdPhi, result_dMdGradPhi,
//                                                                               ADD_TERMS, ADD_JACOBIANS, output_message
//                                                                             );
//
//    BOOST_CHECK( errorCode <= 0 );
//
//    if ( errorCode != 0 ){
//        std::cout << "output_message:\n" << output_message << "\n";
//    }
//
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_result, PK2_answer ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_result, M_answer ) );
////
////    variableMatrix dPK2dGradU(      9, variableVector(  9, 0 ) );
////
////    variableMatrix dPK2dPhi(        9, variableVector(  9, 0 ) );
////
////    variableMatrix dPK2dGradPhi(    9, variableVector( 27, 0 ) );
////
////    variableMatrix dSIGMAdGradU(    9, variableVector(  9, 0 ) );
////
////    variableMatrix dSIGMAdPhi(      9, variableVector(  9, 0 ) );
////
////    variableMatrix dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );
////
////    variableMatrix dMdGradU(       27, variableVector(  9, 0 ) );
////
////    variableMatrix dMdPhi(         27, variableVector(  9, 0 ) );
////
////    variableMatrix dMdGradPhi(     27, variableVector( 27, 0 ) );
////
////    variableMatrix dFpdGradU(          9, variableVector(  9, 0 ) );
////
////    variableMatrix dFpdPhi(            9, variableVector(  9, 0 ) );
////
////    variableMatrix dFpdGradPhi(        9, variableVector( 27, 0 ) );
////
////    variableMatrix dChipdGradU(        9, variableVector(  9, 0 ) );
////
////    variableMatrix dChipdPhi(          9, variableVector(  9, 0 ) );
////
////    variableMatrix dChipdGradPhi(      9, variableVector( 27, 0 ) );
////
////    variableMatrix dGradChipdGradU(   27, variableVector(  9, 0 ) );
////
////    variableMatrix dGradChipdPhi(     27, variableVector(  9, 0 ) );
////
////    variableMatrix dGradChipdGradPhi( 27, variableVector( 27, 0 ) );
////
////    variableType eps = 1e-6;
////
////    for ( unsigned int i = 0; i < 9; i++ ){
////
////        variableVector delta( 9, 0 );
////
////        unsigned int row = i / 3;
////
////        unsigned int col = i % 3;
////
////        delta[ i ] = eps * std::fabs( current_grad_u[ row ][ col ] ) + eps;
////
////        variableType current_grad_u_p[ 3 ][ 3 ];
////        variableType current_grad_u_m[ 3 ][ 3 ];
////
////        for ( unsigned int _i = 0; _i < 3; _i++ ){
////            for ( unsigned int _j = 0; _j < 3; _j++ ){
////                current_grad_u_p[ _i ][ _j ] = current_grad_u[ _i ][ _j ] + delta[ 3 * _i + _j ];
////                current_grad_u_m[ _i ][ _j ] = current_grad_u[ _i ][ _j ] - delta[ 3 * _i + _j ];
////            }
////        }
////
////        variableVector PK2_p,   PK2_m;
////        variableVector SIGMA_p, SIGMA_m;
////        variableVector M_p,     M_m;
////        variableVector SDVS_p = SDVSDefault;
////        variableVector SDVS_m = SDVSDefault;
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_p, current_phi, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_p, SIGMA_p, M_p,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_m, current_phi, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_m, SIGMA_m, M_m,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
////
////            dPK2dGradU[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
////
////            dSIGMAdGradU[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
////
////            dMdGradU[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dFpdGradU[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dChipdGradU[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 27; j++ ){
////
////            dGradChipdGradU[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////    }
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradU,   result_dPK2dGradU ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradU, result_dSIGMAdGradU ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradU,     result_dMdGradU ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradU,       ADD_JACOBIANS[ 0 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradU,     ADD_JACOBIANS[ 3 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradU, ADD_JACOBIANS[ 6 ] ) );
////
////    for ( unsigned int i = 0; i < 9; i++ ){
////
////        variableVector delta( 9, 0 );
////
////        delta[ i ] = eps * std::fabs( current_phi[ i ] ) + eps;
////
////        variableType current_phi_p[ 9 ];
////        variableType current_phi_m[ 9 ];
////
////        for ( unsigned int _i = 0; _i < 3; _i++ ){
////            for ( unsigned int _j = 0; _j < 3; _j++ ){
////                current_phi_p[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] + delta[ 3 * _i + _j ];
////                current_phi_m[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] - delta[ 3 * _i + _j ];
////            }
////        }
////
////        variableVector PK2_p,   PK2_m;
////        variableVector SIGMA_p, SIGMA_m;
////        variableVector M_p,     M_m;
////        variableVector SDVS_p = SDVSDefault;
////        variableVector SDVS_m = SDVSDefault;
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_p, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_p, SIGMA_p, M_p,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_m, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_m, SIGMA_m, M_m,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
////
////            dPK2dPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
////
////            dSIGMAdPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
////
////            dMdPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dFpdPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dChipdPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 27; j++ ){
////
////            dGradChipdPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////    }
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dPhi,   result_dPK2dPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdPhi, result_dSIGMAdPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdPhi,     result_dMdPhi ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdPhi,       ADD_JACOBIANS[ 1 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdPhi,     ADD_JACOBIANS[ 4 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdPhi, ADD_JACOBIANS[ 7 ] ) );
////
////    for ( unsigned int i = 0; i < 27; i++ ){
////
////        variableVector delta( 27, 0 );
////
////        unsigned int row = i / 9;
////
////        unsigned int col = i % 9;
////
////        delta[ i ] = eps * std::fabs( current_grad_phi[ row ][ col ] ) + eps;
////
////        variableType current_grad_phi_p[ 9 ][ 3 ];
////        variableType current_grad_phi_m[ 9 ][ 3 ];
////
////        for ( unsigned int _i = 0; _i < 9; _i++ ){
////            for ( unsigned int _j = 0; _j < 3; _j++ ){
////                current_grad_phi_p[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] + delta[ 3 * _i + _j ];
////                current_grad_phi_m[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] - delta[ 3 * _i + _j ];
////            }
////        }
////
////        variableVector PK2_p,   PK2_m;
////        variableVector SIGMA_p, SIGMA_m;
////        variableVector M_p,     M_m;
////        variableVector SDVS_p = SDVSDefault;
////        variableVector SDVS_m = SDVSDefault;
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_p,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_p, SIGMA_p, M_p,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_m,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_m, SIGMA_m, M_m,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
////
////            dPK2dGradPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
////
////            dSIGMAdGradPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
////
////            dMdGradPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dFpdGradPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 27; j++ ){
////
////            dGradChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////    }
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradPhi,   result_dPK2dGradPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradPhi, result_dSIGMAdGradPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradPhi,     result_dMdGradPhi ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradPhi,       ADD_JACOBIANS[ 2 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradPhi,     ADD_JACOBIANS[ 5 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradPhi, ADD_JACOBIANS[ 8 ] ) );
//
//}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_1){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time
    double t0 = 0.321504 - 0.061917;
    double dt = 0.0001 * 0.061917;

//    std::vector< double > time = { 0.321504, 0.061917 };
    std::vector< double > time = { t0 + dt, dt };

    //Initialize the material parameters
    std::vector< double > fparams = { 2.000000, 3.192203, -145.182846, 2.000000, 100000000.000000, 0.000000, 2.000000, 100000000.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 372.714287, 725.914228, 5.000000, 81.392866, 148.127522, -133.511734, -159.767592, -296.621192, 11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.683567, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 148.127522, -296.621192, 0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u[ 3 ][ 3 ] = { {  0.008816, -0.000085,  0.001147 },
                                        {  0.000318,  0.008293,  0.002643 },
                                        { -0.000092, -0.000849, -0.019404 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  0.006126, -0.000097,  0.000961 },
                                         {  0.000078,  0.005958,  0.001780 },
                                         { -0.000056, -0.000716, -0.014142 } };

    //Initialize the micro displacement
    double current_phi[ 9 ] = { 0.007489, 0.000061, 0.000875, 0.000053, 0.007479, 0.001177, -0.000030, -0.000234, -0.016632 };

    double previous_phi[ 9 ] = { 0.006138, 0.000050, 0.000737, 0.000041, 0.006098, 0.000943, -0.000044, -0.000172, -0.013600 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {  0.000125, -0.000027, -0.000208 },
                                          { -0.000006,  0.000001,  0.000029 },
                                          { -0.000814,  0.000157,  0.000438 },
                                          { -0.000003, -0.000003,  0.000024 },
                                          {  0.000138, -0.000027, -0.00015  },
                                          {  0.000062,  0.000013,  0.000536 },
                                          {  0.000429, -0.000080, -0.000206 },
                                          { -0.000038, -0.000007, -0.000232 },
                                          { -0.000099,  0.000003,  0.000103} };

    double previous_grad_phi[ 9 ][ 3 ] = { {  0.000045, -0.000025, -0.000157 },
                                           { -0.000008,  0.000002,  0.000024 },
                                           { -0.000632,  0.000119,  0.000362 },
                                           { -0.000008,  0.000001,  0.000020 },
                                           {  0.000058, -0.000029, -0.000118 },
                                           {  0.000056,  0.000006,  0.000425 },
                                           {  0.000318, -0.000061, -0.000188 },
                                           { -0.000044, -0.000003, -0.000181 },
                                           {  0.000059,  0.000008,  0.000045 } };

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = { 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 };

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

    std::cout << "SDVS:\n";
    for ( auto s = SDVS.begin( ); s != SDVS.end( ); s++ ){ std::cout << " " << *s << ","; } std::cout << "\n";

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

//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_result, PK2_answer ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_result, M_answer ) );
//
//    variableMatrix dPK2dGradU(      9, variableVector(  9, 0 ) );
//
//    variableMatrix dPK2dPhi(        9, variableVector(  9, 0 ) );
//
//    variableMatrix dPK2dGradPhi(    9, variableVector( 27, 0 ) );
//
//    variableMatrix dSIGMAdGradU(    9, variableVector(  9, 0 ) );
//
//    variableMatrix dSIGMAdPhi(      9, variableVector(  9, 0 ) );
//
//    variableMatrix dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );
//
//    variableMatrix dMdGradU(       27, variableVector(  9, 0 ) );
//
//    variableMatrix dMdPhi(         27, variableVector(  9, 0 ) );
//
//    variableMatrix dMdGradPhi(     27, variableVector( 27, 0 ) );
//
//    variableMatrix dFpdGradU(          9, variableVector(  9, 0 ) );
//
//    variableMatrix dFpdPhi(            9, variableVector(  9, 0 ) );
//
//    variableMatrix dFpdGradPhi(        9, variableVector( 27, 0 ) );
//
//    variableMatrix dChipdGradU(        9, variableVector(  9, 0 ) );
//
//    variableMatrix dChipdPhi(          9, variableVector(  9, 0 ) );
//
//    variableMatrix dChipdGradPhi(      9, variableVector( 27, 0 ) );
//
//    variableMatrix dGradChipdGradU(   27, variableVector(  9, 0 ) );
//
//    variableMatrix dGradChipdPhi(     27, variableVector(  9, 0 ) );
//
//    variableMatrix dGradChipdGradPhi( 27, variableVector( 27, 0 ) );
//
//    variableType eps = 1e-6;
//
//    for ( unsigned int i = 0; i < 9; i++ ){
//
//        variableVector delta( 9, 0 );
//
//        unsigned int row = i / 3;
//
//        unsigned int col = i % 3;
//
//        delta[ i ] = eps * std::fabs( current_grad_u[ row ][ col ] ) + eps;
//
//        variableType current_grad_u_p[ 3 ][ 3 ];
//        variableType current_grad_u_m[ 3 ][ 3 ];
//
//        for ( unsigned int _i = 0; _i < 3; _i++ ){
//            for ( unsigned int _j = 0; _j < 3; _j++ ){
//                current_grad_u_p[ _i ][ _j ] = current_grad_u[ _i ][ _j ] + delta[ 3 * _i + _j ];
//                current_grad_u_m[ _i ][ _j ] = current_grad_u[ _i ][ _j ] - delta[ 3 * _i + _j ];
//            }
//        }
//
//        variableVector PK2_p,   PK2_m;
//        variableVector SIGMA_p, SIGMA_m;
//        variableVector M_p,     M_m;
//        variableVector SDVS_p = SDVSDefault;
//        variableVector SDVS_m = SDVSDefault;
//
//        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_p, current_phi, current_grad_phi,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_p, SIGMA_p, M_p,
//                                                                                  ADD_TERMS, output_message
//                                                                                );
//
//        BOOST_CHECK( errorCode <= 0 );
//
//        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_m, current_phi, current_grad_phi,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_m, SIGMA_m, M_m,
//                                                                                  ADD_TERMS, output_message
//                                                                                );
//
//        BOOST_CHECK( errorCode <= 0 );
//
//        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
//
//            dPK2dGradU[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
//
//            dSIGMAdGradU[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
//
//            dMdGradU[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 9; j++ ){
//
//            dFpdGradU[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 9; j++ ){
//
//            dChipdGradU[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 27; j++ ){
//
//            dGradChipdGradU[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
//
//        }
//
//    }
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradU,   result_dPK2dGradU ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradU, result_dSIGMAdGradU ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradU,     result_dMdGradU ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradU,       ADD_JACOBIANS[ 0 ] ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradU,     ADD_JACOBIANS[ 3 ] ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradU, ADD_JACOBIANS[ 6 ] ) );
//
//    for ( unsigned int i = 0; i < 9; i++ ){
//
//        variableVector delta( 9, 0 );
//
//        delta[ i ] = eps * std::fabs( current_phi[ i ] ) + eps;
//
//        variableType current_phi_p[ 9 ];
//        variableType current_phi_m[ 9 ];
//
//        for ( unsigned int _i = 0; _i < 3; _i++ ){
//            for ( unsigned int _j = 0; _j < 3; _j++ ){
//                current_phi_p[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] + delta[ 3 * _i + _j ];
//                current_phi_m[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] - delta[ 3 * _i + _j ];
//            }
//        }
//
//        variableVector PK2_p,   PK2_m;
//        variableVector SIGMA_p, SIGMA_m;
//        variableVector M_p,     M_m;
//        variableVector SDVS_p = SDVSDefault;
//        variableVector SDVS_m = SDVSDefault;
//
//        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_p, current_grad_phi,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_p, SIGMA_p, M_p,
//                                                                                  ADD_TERMS, output_message
//                                                                                );
//
//        BOOST_CHECK( errorCode <= 0 );
//
//        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_m, current_grad_phi,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_m, SIGMA_m, M_m,
//                                                                                  ADD_TERMS, output_message
//                                                                                );
//
//        BOOST_CHECK( errorCode <= 0 );
//
//        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
//
//            dPK2dPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
//
//            dSIGMAdPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
//
//            dMdPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 9; j++ ){
//
//            dFpdPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 9; j++ ){
//
//            dChipdPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 27; j++ ){
//
//            dGradChipdPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
//
//        }
//
//    }
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dPhi,   result_dPK2dPhi ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdPhi, result_dSIGMAdPhi ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdPhi,     result_dMdPhi ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdPhi,       ADD_JACOBIANS[ 1 ] ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdPhi,     ADD_JACOBIANS[ 4 ] ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdPhi, ADD_JACOBIANS[ 7 ] ) );
//
//    for ( unsigned int i = 0; i < 27; i++ ){
//
//        variableVector delta( 27, 0 );
//
//        unsigned int row = i / 9;
//
//        unsigned int col = i % 9;
//
//        delta[ i ] = eps * std::fabs( current_grad_phi[ row ][ col ] ) + eps;
//
//        variableType current_grad_phi_p[ 9 ][ 3 ];
//        variableType current_grad_phi_m[ 9 ][ 3 ];
//
//        for ( unsigned int _i = 0; _i < 9; _i++ ){
//            for ( unsigned int _j = 0; _j < 3; _j++ ){
//                current_grad_phi_p[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] + delta[ 3 * _i + _j ];
//                current_grad_phi_m[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] - delta[ 3 * _i + _j ];
//            }
//        }
//
//        variableVector PK2_p,   PK2_m;
//        variableVector SIGMA_p, SIGMA_m;
//        variableVector M_p,     M_m;
//        variableVector SDVS_p = SDVSDefault;
//        variableVector SDVS_m = SDVSDefault;
//
//        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_p,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_p, SIGMA_p, M_p,
//                                                                                  ADD_TERMS, output_message
//                                                                                );
//
//        BOOST_CHECK( errorCode <= 0 );
//
//        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_m,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_m, SIGMA_m, M_m,
//                                                                                  ADD_TERMS, output_message
//                                                                                );
//
//        BOOST_CHECK( errorCode <= 0 );
//
//        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
//
//            dPK2dGradPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
//
//            dSIGMAdGradPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
//
//            dMdGradPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 9; j++ ){
//
//            dFpdGradPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 9; j++ ){
//
//            dChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
//
//        }
//
//        for ( unsigned int j = 0; j < 27; j++ ){
//
//            dGradChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
//
//        }
//
//    }
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradPhi,   result_dPK2dGradPhi ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradPhi, result_dSIGMAdGradPhi ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradPhi,     result_dMdGradPhi ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradPhi,       ADD_JACOBIANS[ 2 ] ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradPhi,     ADD_JACOBIANS[ 5 ] ) );
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradPhi, ADD_JACOBIANS[ 8 ] ) );

}

//BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_3){
//    /*!
//     * Test the evaluation of the constitutive model.
//     *
//     */
//
//    //Initialize the time
//    std::vector< double > time = { 0.311879, 1e-12 };
//
//    //Initialize the material parameters
//    std::vector< double > fparams = { 2.000000, 3.192203, -24.544075, 2.000000, 100000000.000000, 0.000000, 2.000000, 100000000.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 5813.030431, 122.720374, 5.000000, -1173.239766, -2034.910254, 431.485800, -434.264640, -1.531541, 11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 8.864684, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, -2034.910254, -1.531541, 0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };
//
//    //Initialize the gradient of the macro displacement
////    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
////                                        { -0.56632666, -0.21399259,  0.16367238 },
////                                        { -0.29129789, -0.22367825, -2.0632945  } };
////
////    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
////                                         {  0.31303067, -1.23910631, -0.93837662 },
////                                         { -0.32571524, -0.95306342, -0.93025257 } };
//
//    double current_grad_u[ 3 ][ 3 ] = { {  0.005779, -0.000507, -0.002797 },
//                                        { -0.000590,  0.007092,  0.006669 },
//                                        {  0.000617, -0.002192, -0.014933 } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { {  0.005779, -0.000507, -0.002797 },
//                                         { -0.000590,  0.007092,  0.006669 },
//                                         {  0.000617, -0.002192, -0.014933 } };
//    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.007084, -0.000013, -0.000612, -0.000057, 0.007095, 0.002610, 0.000783, -0.001177, -0.016065 };
//
//    double previous_phi[ 9 ] = { 0.007084, -0.000013, -0.000612, -0.000057, 0.007095, 0.002610, 0.000783, -0.001177, -0.016065 };
//
//    //Initialize the gradient of the micro displacement
//    double current_grad_phi[ 9 ][ 3 ] = { { -0.000239,  0.000428,  0.000628 },
//                                          {  0.000122,  0.000003,  0.000152 },
//                                          {  0.000800, -0.000382,  0.001604 },
//                                          {  0.000067, -0.000012,  0.000115 },
//                                          { -0.000147,  0.000452,  0.000388 },
//                                          {  0.000662,  0.001270, -0.004211 },
//                                          { -0.000516,  0.000408, -0.000781 },
//                                          { -0.000125, -0.001332,  0.002348 },
//                                          {  0.000174, -0.000539, -0.000433 } };
//
//    double previous_grad_phi[ 9 ][ 3 ] = { {  -0.000239,  0.000428,  0.000628 },
//                                           {   0.000122,  0.000003,  0.000152 },
//                                           {   0.000800, -0.000382,  0.001604 },
//                                           {   0.000067, -0.000012,  0.000115 },
//                                           {  -0.000147,  0.000452,  0.000388 },
//                                           {   0.000662,  0.001270, -0.004211 },
//                                           {  -0.000516,  0.000408, -0.000781 },
//                                           {  -0.000125, -0.001332,  0.002348 },
//                                           {   0.000174, -0.000539, -0.000433 } };
//
//    //Initialize the state variable vector
//    std::vector< double > SDVSDefault = { 0.000056, -0.000004, -0.000042, -0.000006, 0.000068, 0.000099, 0.000023, -0.000063, -0.000124, 0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.022458, 0.000000, -0.000000, 0.000000, -0.000000, 0.000320, 0.000000, 0.000000, -0.000000, 0.000000 };
//
//    //Initialize the additional degree of freedom vectors
//    std::vector< double > current_ADD_DOF;
//    std::vector< std::vector< double > > current_ADD_grad_DOF;
//
//    std::vector< double > previous_ADD_DOF;
//    std::vector< std::vector< double > > previous_ADD_grad_DOF;
//
//    //Initialize the stress measures
//    std::vector< double > PK2_result( 9, 0 );
//
//    std::vector< double > SIGMA_result( 9, 0 );
//
//    std::vector< double > M_result( 27, 0 );
//
//    //Initialize the additional terms vector
//    std::vector< std::vector< double > > ADD_TERMS;
//
//    //Initialize the output message string
//    std::string output_message;
//
//#ifdef DEBUG_MODE
//    tardigradeSolverTools::homotopyMap homotopyDEBUG;
//#endif
//
//    tardigradeSolverTools::floatVector PK2_answer = { 1.72374955e+02,  1.53586446e+01, -9.18594711e-01,  1.34613772e+01,
//        1.42758448e+02, -2.14353350e-02, -1.76155689e+00,  1.77580288e+00,
//        1.41003910e+02 };
//
//    tardigradeSolverTools::floatVector SIGMA_answer = { 176.85221932,  15.84248153,  -2.83823499,  15.84248153,
//       144.5183467 ,   1.85778345,  -2.83823499,   1.85778345,
//       141.99307014 };
//
//    tardigradeSolverTools::floatVector M_answer = { 0.60016978, -0.51045338,  0.61980289,  3.23507063,  1.16925208,
//        1.20665651,  0.56034359, -2.52042105,  1.62648849, -2.61882314,
//       -0.61143851, -1.02227749,  0.67046998,  0.49701267, -0.23999964,
//       -2.77670511,  0.75636495,  1.71897722, -0.49808019,  2.62569695,
//       -0.76002078,  1.23462488, -0.00650376, -2.25591243, -0.73016414,
//        0.74380723,  0.90861263 };
//
//    tardigradeSolverTools::floatVector SDVS_answer = {  7.52790780e-03,  5.06936872e-03, -3.00562092e-04,  5.05251313e-03,
//                                                       -2.63763651e-03,  2.87928507e-04, -3.72441969e-04,  2.25230220e-04,
//                                                       -1.58168041e-03,  7.52714199e-03,  4.60154810e-03, -3.01710599e-04,
//                                                        5.56787320e-03, -2.63687070e-03,  2.58160226e-04, -3.65069825e-04,
//                                                        2.58160226e-04, -1.58168041e-03,  3.83196420e-02, -2.45387739e-02,
//                                                       -8.94156844e-03,  1.04125027e-02, -5.47145710e-04,  2.21268303e-02,
//                                                        1.85810005e-02, -1.14241984e-02, -7.62221572e-03,  4.54399478e-02,
//                                                        3.21567189e-02,  2.10892064e-02,  1.18343147e-02,  1.63154826e-02,
//                                                        1.01776270e-03,  1.47855827e-02,  7.10955065e-03, -2.40560608e-02,
//                                                        9.30082927e-03, -3.52878321e-02,  2.18504370e-02, -2.49577455e-02,
//                                                        9.18791192e-03,  2.13410715e-02,  1.69866062e-02,  1.77613753e-02,
//                                                        2.03937074e-02, -5.42245900e-23,  5.17942700e-03,  3.02931048e-02,
//                                                        1.58261418e-02,  1.60309023e-02, -1.79683900e-22,  2.11178609e-02,
//                                                        8.23638353e-02,  4.30296627e-02,  4.35864383e-02 };
//
//    std::vector< double > SDVS = SDVSDefault;
//
//    int errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
//                                                                                  current_grad_u,  current_phi,  current_grad_phi,
//                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//                                                                                  SDVS,
//                                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
//                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                                  PK2_result, SIGMA_result, M_result,
//                                                                                  ADD_TERMS,
//                                                                                  output_message
//                                                                                  );
//
//    std::cout << "SDVS:\n";
//    for ( auto s = SDVS.begin( ); s != SDVS.end( ); s++ ){ std::cout << " " << *s << ","; } std::cout << "\n";
//
//    BOOST_CHECK( errorCode == 0 );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SDVS_answer,  SDVS ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_answer,   PK2_result ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_answer, SIGMA_result ) );
//
//    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_answer,     M_result ) );
//
//    //Test the Jacobians
//    PK2_result.clear();
//    SIGMA_result.clear();
//    M_result.clear();
//    ADD_TERMS.clear();
//
//    variableMatrix result_dPK2dGradU(      9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dPK2dPhi(        9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dPK2dGradPhi(    9, variableVector( 27, 0 ) );
//
//    variableMatrix result_dSIGMAdGradU(    9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dSIGMAdPhi(      9, variableVector(  9, 0 ) );
//
//    variableMatrix result_dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );
//
//    variableMatrix result_dMdGradU(       27, variableVector(  9, 0 ) );
//
//    variableMatrix result_dMdPhi(         27, variableVector(  9, 0 ) );
//
//    variableMatrix result_dMdGradPhi(     27, variableVector( 27, 0 ) );
//
//    std::vector< variableMatrix > ADD_JACOBIANS;
//
//    SDVS = SDVSDefault;
//
//    errorCode  = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
//                                                                               current_grad_u,  current_phi,  current_grad_phi,
//                                                                               previous_grad_u, previous_phi, previous_grad_phi,
//                                                                               SDVS,
//                                                                               current_ADD_DOF,  current_ADD_grad_DOF,
//                                                                               previous_ADD_DOF, previous_ADD_grad_DOF,
//                                                                               PK2_result, SIGMA_result, M_result,
//                                                                               result_dPK2dGradU, result_dPK2dPhi, result_dPK2dGradPhi,
//                                                                               result_dSIGMAdGradU, result_dSIGMAdPhi, result_dSIGMAdGradPhi,
//                                                                               result_dMdGradU, result_dMdPhi, result_dMdGradPhi,
//                                                                               ADD_TERMS, ADD_JACOBIANS, output_message
//                                                                             );
//
//    BOOST_CHECK( errorCode <= 0 );
//
//    if ( errorCode != 0 ){
//        std::cout << "output_message:\n" << output_message << "\n";
//    }
//
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_result, PK2_answer ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_result, M_answer ) );
////
////    variableMatrix dPK2dGradU(      9, variableVector(  9, 0 ) );
////
////    variableMatrix dPK2dPhi(        9, variableVector(  9, 0 ) );
////
////    variableMatrix dPK2dGradPhi(    9, variableVector( 27, 0 ) );
////
////    variableMatrix dSIGMAdGradU(    9, variableVector(  9, 0 ) );
////
////    variableMatrix dSIGMAdPhi(      9, variableVector(  9, 0 ) );
////
////    variableMatrix dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );
////
////    variableMatrix dMdGradU(       27, variableVector(  9, 0 ) );
////
////    variableMatrix dMdPhi(         27, variableVector(  9, 0 ) );
////
////    variableMatrix dMdGradPhi(     27, variableVector( 27, 0 ) );
////
////    variableMatrix dFpdGradU(          9, variableVector(  9, 0 ) );
////
////    variableMatrix dFpdPhi(            9, variableVector(  9, 0 ) );
////
////    variableMatrix dFpdGradPhi(        9, variableVector( 27, 0 ) );
////
////    variableMatrix dChipdGradU(        9, variableVector(  9, 0 ) );
////
////    variableMatrix dChipdPhi(          9, variableVector(  9, 0 ) );
////
////    variableMatrix dChipdGradPhi(      9, variableVector( 27, 0 ) );
////
////    variableMatrix dGradChipdGradU(   27, variableVector(  9, 0 ) );
////
////    variableMatrix dGradChipdPhi(     27, variableVector(  9, 0 ) );
////
////    variableMatrix dGradChipdGradPhi( 27, variableVector( 27, 0 ) );
////
////    variableType eps = 1e-6;
////
////    for ( unsigned int i = 0; i < 9; i++ ){
////
////        variableVector delta( 9, 0 );
////
////        unsigned int row = i / 3;
////
////        unsigned int col = i % 3;
////
////        delta[ i ] = eps * std::fabs( current_grad_u[ row ][ col ] ) + eps;
////
////        variableType current_grad_u_p[ 3 ][ 3 ];
////        variableType current_grad_u_m[ 3 ][ 3 ];
////
////        for ( unsigned int _i = 0; _i < 3; _i++ ){
////            for ( unsigned int _j = 0; _j < 3; _j++ ){
////                current_grad_u_p[ _i ][ _j ] = current_grad_u[ _i ][ _j ] + delta[ 3 * _i + _j ];
////                current_grad_u_m[ _i ][ _j ] = current_grad_u[ _i ][ _j ] - delta[ 3 * _i + _j ];
////            }
////        }
////
////        variableVector PK2_p,   PK2_m;
////        variableVector SIGMA_p, SIGMA_m;
////        variableVector M_p,     M_m;
////        variableVector SDVS_p = SDVSDefault;
////        variableVector SDVS_m = SDVSDefault;
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_p, current_phi, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_p, SIGMA_p, M_p,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_m, current_phi, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_m, SIGMA_m, M_m,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
////
////            dPK2dGradU[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
////
////            dSIGMAdGradU[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
////
////            dMdGradU[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dFpdGradU[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dChipdGradU[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 27; j++ ){
////
////            dGradChipdGradU[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////    }
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradU,   result_dPK2dGradU ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradU, result_dSIGMAdGradU ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradU,     result_dMdGradU ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradU,       ADD_JACOBIANS[ 0 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradU,     ADD_JACOBIANS[ 3 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradU, ADD_JACOBIANS[ 6 ] ) );
////
////    for ( unsigned int i = 0; i < 9; i++ ){
////
////        variableVector delta( 9, 0 );
////
////        delta[ i ] = eps * std::fabs( current_phi[ i ] ) + eps;
////
////        variableType current_phi_p[ 9 ];
////        variableType current_phi_m[ 9 ];
////
////        for ( unsigned int _i = 0; _i < 3; _i++ ){
////            for ( unsigned int _j = 0; _j < 3; _j++ ){
////                current_phi_p[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] + delta[ 3 * _i + _j ];
////                current_phi_m[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] - delta[ 3 * _i + _j ];
////            }
////        }
////
////        variableVector PK2_p,   PK2_m;
////        variableVector SIGMA_p, SIGMA_m;
////        variableVector M_p,     M_m;
////        variableVector SDVS_p = SDVSDefault;
////        variableVector SDVS_m = SDVSDefault;
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_p, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_p, SIGMA_p, M_p,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_m, current_grad_phi,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_m, SIGMA_m, M_m,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
////
////            dPK2dPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
////
////            dSIGMAdPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
////
////            dMdPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dFpdPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dChipdPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 27; j++ ){
////
////            dGradChipdPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////    }
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dPhi,   result_dPK2dPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdPhi, result_dSIGMAdPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdPhi,     result_dMdPhi ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdPhi,       ADD_JACOBIANS[ 1 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdPhi,     ADD_JACOBIANS[ 4 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdPhi, ADD_JACOBIANS[ 7 ] ) );
////
////    for ( unsigned int i = 0; i < 27; i++ ){
////
////        variableVector delta( 27, 0 );
////
////        unsigned int row = i / 9;
////
////        unsigned int col = i % 9;
////
////        delta[ i ] = eps * std::fabs( current_grad_phi[ row ][ col ] ) + eps;
////
////        variableType current_grad_phi_p[ 9 ][ 3 ];
////        variableType current_grad_phi_m[ 9 ][ 3 ];
////
////        for ( unsigned int _i = 0; _i < 9; _i++ ){
////            for ( unsigned int _j = 0; _j < 3; _j++ ){
////                current_grad_phi_p[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] + delta[ 3 * _i + _j ];
////                current_grad_phi_m[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] - delta[ 3 * _i + _j ];
////            }
////        }
////
////        variableVector PK2_p,   PK2_m;
////        variableVector SIGMA_p, SIGMA_m;
////        variableVector M_p,     M_m;
////        variableVector SDVS_p = SDVSDefault;
////        variableVector SDVS_m = SDVSDefault;
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_p,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_p, SIGMA_p, M_p,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_m,
////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
////                                                                                  PK2_m, SIGMA_m, M_m,
////                                                                                  ADD_TERMS, output_message
////                                                                                );
////
////        BOOST_CHECK( errorCode <= 0 );
////
////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
////
////            dPK2dGradPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
////
////            dSIGMAdGradPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
////
////            dMdGradPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dFpdGradPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 9; j++ ){
////
////            dChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////        for ( unsigned int j = 0; j < 27; j++ ){
////
////            dGradChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
////
////        }
////
////    }
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPK2dGradPhi,   result_dPK2dGradPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dSIGMAdGradPhi, result_dSIGMAdGradPhi ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMdGradPhi,     result_dMdGradPhi ) );
////
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFpdGradPhi,       ADD_JACOBIANS[ 2 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dChipdGradPhi,     ADD_JACOBIANS[ 5 ] ) );
////    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dGradChipdGradPhi, ADD_JACOBIANS[ 8 ] ) );
//
//}
