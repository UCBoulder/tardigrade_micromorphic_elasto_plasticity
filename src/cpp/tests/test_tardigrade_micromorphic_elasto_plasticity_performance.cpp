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

typedef tardigradeMicromorphicTools::errorNode errorNode;
typedef tardigradeMicromorphicTools::errorOut errorOut;

bool tolerantCheck( const std::vector< double > &v1, const std::vector< double > &v2, double eps = 1e-6, double tol = 1e-9 ){

    if ( v1.size( ) != v2.size( ) ){

        return false;

    }

    BOOST_CHECK( v1.size( ) == v2.size( ) );

    const unsigned int len = v1.size( );

    for ( unsigned int i = 0; i < len; i++ ){

        if ( ( std::fabs( v1[ i ] ) < tol ) || ( std::fabs( v2[ i ] ) < tol ) ){

            if ( std::fabs( v1[ i ] - v2[ i ] ) > eps ){

                return false;

            }

        }
        else{

            if ( ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) > eps ) ||
                 ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) > eps ) ){

                std::cout << "v1: " << v1[ i ] << "\n";
                std::cout << "v2: " << v2[ i ] << "\n";
                std::cout << "r1: " << std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) << "\n";
                std::cout << "r2: " << std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) << "\n";

                return false;

            }

        }

    }

    return true;

}

bool tolerantCheck( const double &v1, const double &v2, double eps = 1e-6, double tol = 1e-9 ){

    std::vector< double > _v1 = { v1 };

    std::vector< double > _v2 = { v2 };

    return tolerantCheck( _v1, _v2, eps, tol );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_1, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.0;//0.96875;

    std::vector< double > _time = { 0.321504, 0.061917 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = { 2.000000, 3.192203, -145.182846,
                                      2.000000, 100000000.000000, 0.000000,
                                      2.000000, 100000000.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 372.714287, 725.914228,
                                      5.000000, 81.392866, 148.127522, -133.511734, -159.767592, -296.621192,
                                      11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.683567, 0.000000, 0.000000, 0.000000, 0.000000,
                                      2.000000, 148.127522, -296.621192,
                                      0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double _current_grad_u[ 3 ][ 3 ] = { {  0.008816, -0.000085,  0.001147 },
                                         {  0.000318,  0.008293,  0.002643 },
                                         { -0.000092, -0.000849, -0.019404 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  0.006126, -0.000097,  0.000961 },
                                         {  0.000078,  0.005958,  0.001780 },
                                         { -0.000056, -0.000716, -0.014142 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.007489, 0.000061, 0.000875, 0.000053, 0.007479, 0.001177, -0.000030, -0.000234, -0.016632 };

    double previous_phi[ 9 ] = { 0.006138, 0.000050, 0.000737, 0.000041, 0.006098, 0.000943, -0.000044, -0.000172, -0.013600 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  0.000125, -0.000027, -0.000208 },
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

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

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

    tardigradeSolverTools::floatVector PK2_answer = { -1.40279, 0.0121485, 0.0266627, 0.00381226, -1.45011, 0.0418199, 0.0338264, 0.0864571, -2.91531 };

    tardigradeSolverTools::floatVector SIGMA_answer = { -0.493708, 0.0121523, 0.0789793, 0.0121523, -0.517776, 0.106717, 0.0789793, 0.106717, -4.81709 };

    tardigradeSolverTools::floatVector M_answer = { 8.45418e-05, -4.04181e-06, -0.000551314, -2.62816e-06, 9.35849e-05, 4.20583e-05, 0.000290068, -2.54552e-05, -6.72292e-05, -1.81165e-05, 6.55509e-07, 0.000106142, -1.94223e-06, -1.8193e-05, 8.36339e-06, -5.40107e-05, -4.57942e-06, 2.03912e-06, -0.000144554, 2.02458e-05, 0.00030473, 1.70095e-05, -0.000104047, 0.00037256, -0.000143114, -0.000161128, 7.22794e-05 };

    tardigradeSolverTools::floatVector SDVS_answer = { 0.00458246, 3.26705e-05, 0.000297186, 0.000112011, 0.00414933, 0.000768014, 0.000230621, 0.000356534, -0.00861812, 1.43254e-11, 1.68409e-13, 8.44243e-13, 2.5468e-13, 1.31286e-11, 1.89224e-12, 7.70904e-13, 1.41952e-12, -2.64072e-11, 2.79233e-17, 5.56405e-20, -7.87013e-18, -3.37087e-18, 1.51828e-17, -1.01737e-17, 2.41855e-17, -1.34546e-18, -5.24753e-17, -2.32027e-18, 1.54445e-17, 6.40194e-18, -3.58006e-19, -4.96139e-18, -4.77126e-18, -3.94101e-18, -1.629e-17, -4.58044e-17, -1.27547e-16, 1.15812e-17, 6.96373e-17, 2.14292e-17, -2.77694e-17, 1.05305e-16, 6.73349e-18, 1.54864e-17, 3.42911e-17, 0.170641, 4.09748e-26, 2.09674e-24, 0, 0, 0.0172535, 0, 5.60883e-25, 0, 1.61587e-23 };

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

    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );

    BOOST_TEST( PK2_answer ==  PK2_result, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_answer == SIGMA_result, CHECK_PER_ELEMENT );

    BOOST_TEST( M_answer == M_result, CHECK_PER_ELEMENT );

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

    BOOST_TEST( PK2_result == PK2_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_result == SIGMA_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( M_result == M_answer, CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dGradU   ) == tardigradeVectorTools::appendVectors( result_dPK2dGradU ), CHECK_PER_ELEMENT );
    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdGradU ) == tardigradeVectorTools::appendVectors( result_dSIGMAdGradU ), CHECK_PER_ELEMENT );
    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdGradU     ) == tardigradeVectorTools::appendVectors( result_dMdGradU ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dFpdGradU       ) == tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 0 ] ), CHECK_PER_ELEMENT );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-6, 1e-6 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-6, 1e-6 ) );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dPhi   ) == tardigradeVectorTools::appendVectors( result_dPK2dPhi ), CHECK_PER_ELEMENT );
    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdPhi ) == tardigradeVectorTools::appendVectors( result_dSIGMAdPhi ), CHECK_PER_ELEMENT );
    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdPhi     ) == tardigradeVectorTools::appendVectors( result_dMdPhi ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdPhi )      , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 1 ] ), 5e-5, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-6, 1e-6 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 7 ] ), 1e-6, 1e-6 ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        variableVector delta( 27, 0 );

        unsigned int row = i / 9;

        unsigned int col = i % 9;

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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   1e-5, 1e-3 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 1e-5, 1e-3 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     1e-5, 1e-3 ) );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 1e-5, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 1e-5, 1e-6 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 1e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_2, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    std::cout << "\n\nDIFFICULT 2\n\n";

    //Initialize the time

    double s = 1.0;//0.96875;

    std::vector< double > _time = { 0.259587, 0.051598 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = { 2.000000, 100000000.000000, 0.000000,
                                      2.000000, 3.192203, -31.678450,
                                      2.000000, 100000000.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 696.441593, 126.713800,
                                      5.000000, -18.674980, -37.817315, 15.177654, -24.071197, -5.861821,
                                      11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 792.523471, 0.000000, 0.000000, 0.000000, 0.000000,
                                      2.000000, -37.817315, -5.861821,
                                      0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double _current_grad_u[ 3 ][ 3 ] = { { -0.005286,  0.007085, 0.002737 },
                                         {  0.007764, -0.025503, 0.005214 },
                                         {  0.022879, -0.031735, 0.043807 } };

    double previous_grad_u[ 3 ][ 3 ] = { { 0.005248, 0.000254,  -0.000553 },
                                         { 0.000184, 0.004900,   0.001188 },
                                         { 0.000570, -0.001337, -0.011783 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.012031, -0.003652, -0.010154, -0.003904, 0.000124, -0.005378, -0.012541, -0.013293, -0.004787 };

    double previous_phi[ 9 ] = { 0.015309, 0.000000, -0.000073, -0.000000, 0.015309, 0.000153, 0.000058, -0.000121, -0.033897 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { { -0.000016, -0.000057,  0.000361 },
                                           {  0.000004,  0.000018, -0.000148 },
                                           { -0.000387, -0.000076, -0.000091 },
                                           { -0.000019, -0.000025, -0.000004 },
                                           {  0.000010, -0.000046,  0.000397 },
                                           {  0.000235, -0.000011,  0.000739 },
                                           {  0.000841,  0.000230,  0.000407 },
                                           { -0.000195,  0.000055, -0.001201 },
                                           {  0.000065, -0.000164, -0.000076 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {  0.000001, -0.000003, -0.000018 },
                                           {  0.000000, -0.000000, -0.000002 },
                                           { -0.000050, -0.000019, -0.000089 },
                                           { -0.000000, -0.000000, -0.000002 },
                                           {  0.000001, -0.000003, -0.000015 },
                                           { -0.000017, -0.000020,  0.000187 },
                                           {  0.000040,  0.000015,  0.000070 },
                                           {  0.000013,  0.000016, -0.000147 },
                                           {  0.000002, -0.000005,  0.000003 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

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

    tardigradeSolverTools::floatVector PK2_answer = { -1.40279, 0.0121485, 0.0266627, 0.00381226, -1.45011, 0.0418199, 0.0338264, 0.0864571, -2.91531 };

    tardigradeSolverTools::floatVector SIGMA_answer = { -0.493708, 0.0121523, 0.0789793, 0.0121523, -0.517776, 0.106717, 0.0789793, 0.106717, -4.81709 };

    tardigradeSolverTools::floatVector M_answer = { 8.45418e-05, -4.04181e-06, -0.000551314, -2.62816e-06, 9.35849e-05, 4.20583e-05, 0.000290068, -2.54552e-05, -6.72292e-05, -1.81165e-05, 6.55509e-07, 0.000106142, -1.94223e-06, -1.8193e-05, 8.36339e-06, -5.40107e-05, -4.57942e-06, 2.03912e-06, -0.000144554, 2.02458e-05, 0.00030473, 1.70095e-05, -0.000104047, 0.00037256, -0.000143114, -0.000161128, 7.22794e-05 };

    tardigradeSolverTools::floatVector SDVS_answer = { 0.00458246, 3.26705e-05, 0.000297186, 0.000112011, 0.00414933, 0.000768014, 0.000230621, 0.000356534, -0.00861812, 1.43254e-11, 1.68409e-13, 8.44243e-13, 2.5468e-13, 1.31286e-11, 1.89224e-12, 7.70904e-13, 1.41952e-12, -2.64072e-11, 2.79233e-17, 5.56405e-20, -7.87013e-18, -3.37087e-18, 1.51828e-17, -1.01737e-17, 2.41855e-17, -1.34546e-18, -5.24753e-17, -2.32027e-18, 1.54445e-17, 6.40194e-18, -3.58006e-19, -4.96139e-18, -4.77126e-18, -3.94101e-18, -1.629e-17, -4.58044e-17, -1.27547e-16, 1.15812e-17, 6.96373e-17, 2.14292e-17, -2.77694e-17, 1.05305e-16, 6.73349e-18, 1.54864e-17, 3.42911e-17, 0.170641, 4.09748e-26, 2.09674e-24, 0, 0, 0.0172535, 0, 5.60883e-25, 0, 1.61587e-23 };

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

    if ( errorCode != 0 ){
        std::cout << "output_message:\n" << output_message << "\n";
    }

    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );

    BOOST_TEST( PK2_answer ==  PK2_result, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_answer == SIGMA_result, CHECK_PER_ELEMENT );

    BOOST_TEST( M_answer == M_result, CHECK_PER_ELEMENT );

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
//    BOOST_TEST( PK2_result == PK2_answer, CHECK_PER_ELEMENT );
//
//    BOOST_TEST( SIGMA_result == SIGMA_answer, CHECK_PER_ELEMENT );
//
//    BOOST_TEST( M_result == M_answer, CHECK_PER_ELEMENT );
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
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dGradU   ) == tardigradeVectorTools::appendVectors( result_dPK2dGradU ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdGradU ) == tardigradeVectorTools::appendVectors( result_dSIGMAdGradU ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdGradU     ) == tardigradeVectorTools::appendVectors( result_dMdGradU ), CHECK_PER_ELEMENT );
//
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dFpdGradU       ) == tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 0 ] ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-6, 1e-6 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-6, 1e-6 ) );
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
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dPhi   ) == tardigradeVectorTools::appendVectors( result_dPK2dPhi ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdPhi ) == tardigradeVectorTools::appendVectors( result_dSIGMAdPhi ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdPhi     ) == tardigradeVectorTools::appendVectors( result_dMdPhi ), CHECK_PER_ELEMENT );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdPhi )      , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 1 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-6, 1e-6 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 7 ] ), 1e-6, 1e-6 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   1e-5, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 1e-5, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     1e-5, 1e-3 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 1e-5, 1e-6 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 1e-5, 1e-6 ) );
//
}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_4, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    std::cout << "\n\nDIFFICULT 4\n\n";

    //Initialize the time

    double s = 1.0;//0.96875;

    std::vector< double > _time = { 0.430000, 0.002500 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = { 2.000000, 100000000.000000, 0.000000,
                                      2.000000, 3.192203, -31.678450,
                                      2.000000, 100000000.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 0.000000, 0.000000,
                                      2.000000, 696.441593, 126.713800,
                                      5.000000, -18.674980, -37.817315, 15.177654, -24.071197, -5.861821,
                                      11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 792.523471, 0.000000, 0.000000, 0.000000, 0.000000,
                                      2.000000, -37.817315, -5.861821,
                                      0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] = { { 0.009449,  0.006308,  0.001352 },
                                         { 0.006418,  0.009968, -0.003246 },
                                         { 0.001656, -0.001550, -0.023535 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  0.020475, -0.000683, -0.005150 },
                                         {  0.000439,  0.020423,  0.002642 },
                                         { -0.000052,  0.000475, -0.042953 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.026201, -0.000003, -0.000039, -0.000012, 0.026207, 0.000028, 0.000033, -0.000027, -0.056975 };

    double previous_phi[ 9 ] = { 0.024204, -0.000000, -0.000036, 0.000008, 0.024206, 0.000025, 0.000029, -0.000020, -0.052344 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  0.000003,  0.000004,  0.000003 },
                                           { -0.000021,  0.000009,  0.000009 },
                                           {  0.000010,  0.000005,  0.000581 },
                                           {  0.000012,  0.000003, -0.000003 },
                                           { -0.000008, -0.000007,  0.000004 },
                                           {  0.000005,  0.000015, -0.000457 },
                                           { -0.000010, -0.000002, -0.000520 },
                                           { -0.000002, -0.000015,  0.000421 },
                                           {  0.000017, -0.000027,  0.000003 } };

    double previous_grad_phi[ 9 ][ 3 ] = { { -0.000010,  0.000012,  0.000005 },
                                           {  0.000010,  0.000005, -0.000004 },
                                           {  0.000019,  0.000003,  0.000566 },
                                           { -0.000012,  0.000000,  0.000004 },
                                           { -0.000016,  0.000010,  0.000005 },
                                           {  0.000004,  0.000020, -0.000389 },
                                           { -0.000017, -0.000003, -0.000459 },
                                           { -0.000003, -0.000017,  0.000306 },
                                           { -0.000010,  0.000014,  0.000003 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = { 0.018733, 0.000239, -0.002176, 0.000239, 0.018860, 0.001321, -0.002215, 0.001340, -0.036451, 0.018742, 0.000234, -0.001957, 0.000234, 0.018862, 0.001186, -0.002266, 0.001380, -0.036461, 0.000000, 0.000000, 0.000002, -0.000000, 0.000000, -0.000001, 0.000001, 0.000000, 0.000028, -0.000000, 0.000000, -0.000001, -0.000000, -0.000000, 0.000001, 0.000000, 0.000001, -0.000020, 0.000001, 0.000000, 0.000025, 0.000000, 0.000001, -0.000018, -0.000000, 0.000000, -0.000003, -0.000000, 2.211167, 0.000000, 0.000000, -0.000000, 0.000000, 0.073665, 0.000000, 0.000000, -0.000000 };

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

    tardigradeSolverTools::floatVector PK2_answer = { -1.40279, 0.0121485, 0.0266627, 0.00381226, -1.45011, 0.0418199, 0.0338264, 0.0864571, -2.91531 };

    tardigradeSolverTools::floatVector SIGMA_answer = { -0.493708, 0.0121523, 0.0789793, 0.0121523, -0.517776, 0.106717, 0.0789793, 0.106717, -4.81709 };

    tardigradeSolverTools::floatVector M_answer = { 8.45418e-05, -4.04181e-06, -0.000551314, -2.62816e-06, 9.35849e-05, 4.20583e-05, 0.000290068, -2.54552e-05, -6.72292e-05, -1.81165e-05, 6.55509e-07, 0.000106142, -1.94223e-06, -1.8193e-05, 8.36339e-06, -5.40107e-05, -4.57942e-06, 2.03912e-06, -0.000144554, 2.02458e-05, 0.00030473, 1.70095e-05, -0.000104047, 0.00037256, -0.000143114, -0.000161128, 7.22794e-05 };

    tardigradeSolverTools::floatVector SDVS_answer = { 0.00458246, 3.26705e-05, 0.000297186, 0.000112011, 0.00414933, 0.000768014, 0.000230621, 0.000356534, -0.00861812, 1.43254e-11, 1.68409e-13, 8.44243e-13, 2.5468e-13, 1.31286e-11, 1.89224e-12, 7.70904e-13, 1.41952e-12, -2.64072e-11, 2.79233e-17, 5.56405e-20, -7.87013e-18, -3.37087e-18, 1.51828e-17, -1.01737e-17, 2.41855e-17, -1.34546e-18, -5.24753e-17, -2.32027e-18, 1.54445e-17, 6.40194e-18, -3.58006e-19, -4.96139e-18, -4.77126e-18, -3.94101e-18, -1.629e-17, -4.58044e-17, -1.27547e-16, 1.15812e-17, 6.96373e-17, 2.14292e-17, -2.77694e-17, 1.05305e-16, 6.73349e-18, 1.54864e-17, 3.42911e-17, 0.170641, 4.09748e-26, 2.09674e-24, 0, 0, 0.0172535, 0, 5.60883e-25, 0, 1.61587e-23 };

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

    if ( errorCode != 0 ){
        std::cout << "output_message:\n" << output_message << "\n";
    }

    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );

    BOOST_TEST( PK2_answer ==  PK2_result, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_answer == SIGMA_result, CHECK_PER_ELEMENT );

    BOOST_TEST( M_answer == M_result, CHECK_PER_ELEMENT );

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
//    BOOST_TEST( PK2_result == PK2_answer, CHECK_PER_ELEMENT );
//
//    BOOST_TEST( SIGMA_result == SIGMA_answer, CHECK_PER_ELEMENT );
//
//    BOOST_TEST( M_result == M_answer, CHECK_PER_ELEMENT );
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
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dGradU   ) == tardigradeVectorTools::appendVectors( result_dPK2dGradU ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdGradU ) == tardigradeVectorTools::appendVectors( result_dSIGMAdGradU ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdGradU     ) == tardigradeVectorTools::appendVectors( result_dMdGradU ), CHECK_PER_ELEMENT );
//
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dFpdGradU       ) == tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 0 ] ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-6, 1e-6 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-6, 1e-6 ) );
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
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dPhi   ) == tardigradeVectorTools::appendVectors( result_dPK2dPhi ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdPhi ) == tardigradeVectorTools::appendVectors( result_dSIGMAdPhi ), CHECK_PER_ELEMENT );
//    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdPhi     ) == tardigradeVectorTools::appendVectors( result_dMdPhi ), CHECK_PER_ELEMENT );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdPhi )      , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 1 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-6, 1e-6 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 7 ] ), 1e-6, 1e-6 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   1e-5, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 1e-5, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     1e-5, 1e-3 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 1e-5, 1e-6 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 1e-5, 1e-6 ) );
//
}
