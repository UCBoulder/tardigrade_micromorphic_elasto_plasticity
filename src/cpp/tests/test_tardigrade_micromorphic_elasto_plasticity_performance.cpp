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

void cleanAnswer( std::vector< double > &a, double tol = 1e-9 ){

    for ( auto v = a.begin( ); v != a.end( ); v++ ){

        if ( std::fabs( *v ) < tol ){

            *v = 0.;

        }

    }

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

    cleanAnswer( SDVS_answer );

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

    tardigradeSolverTools::floatVector PK2_answer = { 7.20108, -0.010419, 0.40178, -0.0312379, 7.31865, -0.501715, -0.443789, 0.641519, 6.95971 };

    tardigradeSolverTools::floatVector SIGMA_answer = { 7.1785, -0.0545498, -0.0987932, -0.0545498, 7.33568, 0.109286, -0.0987932, 0.109286, 6.81404 };

    tardigradeSolverTools::floatVector M_answer = { -0.0263527, 0.00471758, -0.308025, -0.0137948, 0.00234765, 0.183845, 0.625938, -0.133865, 0.0583556, -0.043483, 0.0126769, -0.0571493, -0.020039, -0.0335559, 0.00404204, 0.180415, 0.0190195, -0.124223, 0.270528, -0.109253, -0.0708277, -0.0162357, 0.321918, 0.550908, 0.267257, -0.827036, -0.0448194 };

    tardigradeSolverTools::floatVector SDVS_answer = { -0.00884271, 0.00720015, 0.0131995, 0.00722664, -0.0283276, -0.0133543, 0.0132149, -0.0134292, 0.0404522, -0.00913627, 0.00712621, 0.0127927, 0.00796593, -0.0286991, -0.0116354, 0.0123124, -0.0129925, 0.0411172, 1.54245e-05, 3.73484e-06, 6.78738e-06, -6.96842e-06, -2.49434e-07, -1.78553e-05, 2.20578e-05, 2.42514e-06, 1.63091e-06, -1.23956e-05, -2.43171e-06, -1.27993e-05, 4.98752e-06, -1.02922e-08, 2.05693e-05, -1.89819e-05, 1.07569e-06, -4.26493e-05, 4.306e-05, 1.37901e-05, 3.24546e-05, -1.70277e-05, -1.20242e-06, -8.77688e-05, -2.0412e-05, -3.72455e-06, -2.73567e-05, 1.75165e-24, 1.11041, 3.23707e-24, -4.92829e-24, 1.85212e-24, 4.6433e-25, 0.093562, 5.34969e-24, 1.98121e-23, -4.48109e-24 };

    cleanAnswer( SDVS_answer );

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

//    std::cout << "PK2  : "; for ( auto v = PK2_result.begin( );   v != PK2_result.end( );   v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SIGMA: "; for ( auto v = SIGMA_result.begin( ); v != SIGMA_result.end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "M    : "; for ( auto v = M_result.begin( );     v != M_result.end( );     v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SDVS : "; for ( auto v = SDVS.begin( );         v != SDVS.end( );         v++ ){ std::cout << *v << ", "; } std::cout << "\n";
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 7 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     5e-4, 1e-3 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 1e-4, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 1e-4, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 1e-4, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_3, * boost::unit_test::tolerance( 1e-2 ) ){ //TODO: Maybe there is an error in the gradients w.r.t. the displacement?
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

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

    tardigradeSolverTools::floatVector PK2_answer = { -3.34995, 0.0663803, 0.0199674, 0.0616158, -3.34511, 0.0108698, 0.0150697, -0.0607778, -3.4235 };

    tardigradeSolverTools::floatVector SIGMA_answer = { -3.32432, 0.0463883, 0.0132939, 0.0463883, -3.32068, -0.0195832, 0.0132939, -0.0195832, -3.53427 };

    tardigradeSolverTools::floatVector M_answer = { 0.00222828, -0.0157869, 0.005696, 0.00896147, -0.0058176, 0.00513337, -0.00688722, -0.00230328, 0.0143165, 0.00294149, 0.00708115, 0.00636042, 0.00208749, -0.00531097, 0.00927057, -0.00347852, -0.0104175, -0.0224464, 0.00335021, 0.00674239, 0.488404, -0.00505116, 0.00444042, -0.387893, -0.440723, 0.355776, 0.000391993 };

    tardigradeSolverTools::floatVector SDVS_answer = { 0.0106226, 0.00652946, 0.00147774, 0.00653205, 0.0111543, -0.00240772, 0.00159416, -0.00256405, -0.0210245, 0.0106306, 0.00652352, 0.00154492, 0.00652073, 0.0111647, -0.00241912, 0.00171859, -0.00270385, -0.0210408, 1.44984e-07, -8.65097e-08, -2.32573e-06, -4.4939e-08, -8.31985e-08, 2.9239e-06, 7.8105e-07, -1.11731e-07, 1.13201e-05, 7.2454e-08, 3.90015e-08, 2.80645e-06, -1.60316e-07, 1.3586e-07, -2.40339e-06, -1.95085e-07, 7.90054e-07, -5.66129e-06, 6.51471e-07, 1.59658e-07, 8.84056e-06, 3.58077e-08, 6.19628e-07, -3.91159e-06, 2.11321e-08, -5.32943e-08, 4.66545e-06, -1.45165e-27, 8.98758, -4.45744e-28, 4.86672e-28, -4.81152e-28, 5.13943e-27, 0.110357, 1.73479e-27, -2.21576e-27, -4.70888e-27 };

    cleanAnswer( SDVS_answer );

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

//    std::cout << "PK2  : "; for ( auto v = PK2_result.begin( );   v != PK2_result.end( );   v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SIGMA: "; for ( auto v = SIGMA_result.begin( ); v != SIGMA_result.end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "M    : "; for ( auto v = M_result.begin( );     v != M_result.end( );     v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SDVS : "; for ( auto v = SDVS.begin( );         v != SDVS.end( );         v++ ){ std::cout << *v << ", "; } std::cout << "\n";

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

    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );

    BOOST_TEST( PK2_result == PK2_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_result == SIGMA_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( M_result == M_answer, CHECK_PER_ELEMENT );

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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-2, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-2, 1e-3 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     5e-4, 1e-3 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 5e-5, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_4, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.0;//0.96875;

    std::vector< double > _time = { 0.420000, 0.01 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = { 4.000000, 100000000.000000, 0.000000, 1e-2, 0.1,
                                      4.000000, 3.192203, -31.678450, 1e-2, 0.1,
                                      4.000000, 100000000.000000, 0.000000, 1e-2, 0.1,
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
    double _current_grad_u[ 3 ][ 3 ] = { {  0.026125,  0.026303, -0.011998 },
                                         {  0.028529, -0.008631, -0.014277 },
                                         { -0.032270, -0.003117, -0.015042 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  0.011685, -0.000193, -0.003084 },
                                         { -0.000166,  0.011788,  0.005098 },
                                         {  0.000382, -0.000711, -0.025964 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.025131, 0.000168, -0.015482, 0.000148, 0.020700, 0.001777, -0.025419, 0.002269, -0.049488 };

    double previous_phi[ 9 ] = { 0.025320, -0.000000, -0.000123, -0.000000, 0.025319, 0.000204, 0.000098, -0.000163, -0.055321 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { { -0.000023, -0.000037,  0.000070 },
                                           { -0.000009, -0.000070,  0.000070 },
                                           {  0.000520, -0.000009,  0.001549 },
                                           { -0.000069,  0.000324, -0.000014 },
                                           {  0.000045,  0.000075,  0.000046 },
                                           { -0.000118,  0.000380,  0.000111 },
                                           { -0.000367,  0.000078, -0.001478 },
                                           {  0.000176, -0.000432, -0.000175 },
                                           {  0.000058,  0.000047, -0.000009 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {  -0.000004,  0.000007,  0.000031 },
                                           {   0.000000,  0.000000,  0.000001 },
                                           {   0.000265,  0.000005,  0.000091 },
                                           {  -0.000000,  0.000000,  0.000001 },
                                           {  -0.000004,  0.000006,  0.000030 },
                                           {   0.000004,  0.000262, -0.000148 },
                                           {  -0.000210, -0.000002, -0.000070 },
                                           {  -0.000002, -0.000211,  0.000113 },
                                           {  -0.000001,  0.000001, -0.000019 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = {  0.003776, -0.000034, -0.000368, -0.000034, 0.003809, 0.000597, -0.000370, 0.000601, -0.007534, 0.003777, -0.000036, -0.000324, -0.000036, 0.003812, 0.000526, -0.000382, 0.000620, -0.007538, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000003, 0.000000, 0.000001, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.000003, -0.000002, 0.000002, 0.000000, 0.000001, 0.000000, 0.000002, -0.000001, -0.000000, 0.000000, -0.000000, 0.000000, 0.176535, 0.000000, -0.000000, 0.000000, 0.000000, 0.014894, -0.000000, 0.000000, -0.000000 };

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

    tardigradeSolverTools::floatVector PK2_answer = { -0.42083, 0.188719, -0.284488, 0.113652, -0.606198, 0.198919, 0.154624, -0.272384, -0.412977 };

    tardigradeSolverTools::floatVector SIGMA_answer = { -0.49984, 0.07922, -0.0612982, 0.07922, -0.600274, -0.024636, -0.0612982, -0.024636, -0.592951 };

    tardigradeSolverTools::floatVector M_answer = { -0.0112254, 0.00535196, 0.439085, -0.0777306, 0.0376894, -0.132622, -0.334233, 0.158964, 0.0606504, -0.0371522, -0.067423, -0.0265724, 0.27262, 0.0507418, 0.335067, 0.0709937, -0.374721, 0.0389545, 0.0377269, 0.0561764, 1.24871, -0.0416176, 0.0361598, 0.023124, -1.20441, -0.108903, 0.0370146 };

    tardigradeSolverTools::floatVector SDVS_answer = { 0.026465, 0.0278464, -0.0216844, 0.0278898, -0.00914431, -0.00885806, -0.0219855, -0.00894281, -0.0138897, 0.0261713, 0.0276253, -0.0208302, 0.0283241, -0.00897724, -0.00911191, -0.0228368, -0.00914724, -0.0137652, 1.91504e-05, 9.4537e-06, 6.63661e-05, 3.22352e-06, 1.01035e-05, 2.17514e-05, 1.32479e-05, 8.0438e-06, 4.98021e-05, 2.4575e-06, -6.58701e-06, 1.81752e-05, -8.00738e-07, -2.80077e-06, 5.03653e-06, 1.42835e-05, 7.20726e-06, 4.16155e-05, 1.20598e-05, 9.93543e-06, 5.18231e-05, 1.35541e-05, -1.7557e-06, 4.0168e-05, -1.8236e-05, -6.62803e-06, -7.06047e-05, 9.05459e-23, 5.67313, -9.38528e-24, -5.87456e-24, -6.34399e-24, -5.13721e-23, 0.107536, 9.57509e-23, 4.31927e-23, 9.86484e-23 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS = SDVSDefault;

    // Explore continuation approach

    tardigradeSolverTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeSolverTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                                                    currentDeformationGradient, currentMicroDeformation,
                                                                                    currentGradientMicroDeformation );

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                                                    previousDeformationGradient, previousMicroDeformation,
                                                                                    previousGradientMicroDeformation );

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

//    std::cout << "PK2  : "; for ( auto v = PK2_result.begin( );   v != PK2_result.end( );   v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SIGMA: "; for ( auto v = SIGMA_result.begin( ); v != SIGMA_result.end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "M    : "; for ( auto v = M_result.begin( );     v != M_result.end( );     v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SDVS : "; for ( auto v = SDVS.begin( );         v != SDVS.end( );         v++ ){ std::cout << *v << ", "; } std::cout << "\n";

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

    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );

    BOOST_TEST( PK2_result == PK2_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_result == SIGMA_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( M_result == M_answer, CHECK_PER_ELEMENT );

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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     5e-4, 1e-3 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 5e-5, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_5, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.0;//0.96875;

    std::vector< double > _time = { 0.428125, 0.000625 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = {  2.000000, 100000000.000000, 0.000000, 2.000000, 3.192203, -31.678450, 2.000000, 100000000.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 696.441593, 126.713800, 5.000000, -18.674980, -37.817315, 15.177654, -24.071197, -5.861821, 11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 792.523471, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, -37.817315, -5.861821, 0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] =
{ { -0.008841,  0.021031,  0.001568 },
  {  0.018150, -0.010246, -0.000261 },
  {  0.011430, -0.016923,  0.009132 } };

    double previous_grad_u[ 3 ][ 3 ] = 
{ {  0.021098, -0.001426,  0.003623 },
  { -0.000580,  0.021356, -0.004474 },
  { -0.002557,  0.004294, -0.044437 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.027871, -0.000358, 0.000102, -0.000380, 0.027874, -0.000076, -0.000093, 0.000069, -0.060979 };

    double previous_phi[ 9 ] = { 0.023978, -0.000002, 0.000136, 0.000003, 0.023979, -0.000129, -0.000113, 0.000110, -0.051816 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  0.000017, -0.000008, -0.000004 },
                                           { -0.000066, -0.000023,  0.000003 },
                                           { -0.000099, -0.000025,  0.000528 },
                                           {  0.000049,  0.000041, -0.000008 },
                                           {  0.000006, -0.000020, -0.000003 },
                                           { -0.000028, -0.000082, -0.000506 },
                                           {  0.000098,  0.000027, -0.000533 },
                                           {  0.000027,  0.000076,  0.000522 },
                                           {  0.000050, -0.000063, -0.000010 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {  -0.000009,  0.000012, -0.000003 },
                                           {   0.000010,  0.000010, -0.000001 },
                                           {  -0.000039,  0.000002,  0.000540 },
                                           {  -0.000014, -0.000005,  0.000005 },
                                           {  -0.000014,  0.000009, -0.000003 },
                                           {   0.000005, -0.000043, -0.000466 },
                                           {   0.000028, -0.000012, -0.000435 },
                                           {  -0.000013,  0.000032,  0.000372 },
                                           {  -0.000011,  0.000012,  0.000003 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = {  0.019995, -0.000409, 0.000719, -0.000409, 0.020256, -0.000417, 0.000715, -0.000402, -0.038958, 0.019998, -0.000410, 0.000545, -0.000411, 0.020256, -0.000245, 0.000651, -0.000312, -0.038961, 0.000000, 0.000000, -0.000001, -0.000000, 0.000000, 0.000000, -0.000002, -0.000000, 0.000029, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, -0.000000, -0.000000, -0.000003, -0.000025, -0.000002, -0.000000, 0.000026, -0.000000, -0.000002, -0.000023, -0.000000, 0.000000, 0.000001, 0.000000, 4.055659, 0.000000, 0.000000, -0.000000, -0.000000, 0.078750, -0.000000, -0.000000, -0.000000 };

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

    tardigradeSolverTools::floatVector PK2_answer = { -9.62479, 0.38928, 0.361314, 0.513209, -9.65339, -0.531695, -0.0320208, 0.105003, -8.98122 };

    tardigradeSolverTools::floatVector SIGMA_answer = { -9.34451, 0.375142, 0.127635, 0.375142, -9.37008, -0.165708, 0.127635, -0.165708, -8.99023 };

    tardigradeSolverTools::floatVector M_answer = { 0.0134268, -0.0543831, -0.0879468, 0.038986, 0.00690052, -0.0112154, 0.0777729, 0.0124799, 0.0411076, -0.00806668, -0.0164, -0.00616632, 0.0329602, -0.0166184, -0.0747552, 0.0103869, 0.062287, -0.0513045, -0.0024238, 0.00359776, 0.421367, -0.0065458, -0.00293546, -0.407405, -0.40262, 0.393485, -0.0191425 };

    tardigradeSolverTools::floatVector SDVS_answer = { -0.00540907, 0.0201288, 0.00644233, 0.0201198, -0.00676143, -0.00849111, 0.00690441, -0.00908192, 0.0157183, -0.00540473, 0.0202057, 0.00563266, 0.0202134, -0.00695702, -0.00741807, 0.00634141, -0.00837407, 0.0158932, 3.54253e-06, 1.37056e-06, -6.61657e-06, -1.12583e-06, 8.89558e-08, 6.56662e-06, 5.86212e-06, -1.08602e-06, -2.05309e-05, -5.75537e-07, 5.87735e-07, 6.47761e-06, -2.62242e-06, -2.63729e-06, -7.29034e-06, -1.04981e-06, 4.21174e-06, 2.41893e-05, 5.82389e-06, -2.89086e-07, -2.83413e-05, -4.85834e-07, 4.38003e-06, 3.15114e-05, -8.90182e-07, 1.24124e-06, 1.37233e-05, 5.15022e-24, 116.919, 4.49973e-23, 3.5662e-22, 0, 0, 0.19808, 1.18479e-21, 0, 1.68029e-21 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeSolverTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeSolverTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                                                    currentDeformationGradient, currentMicroDeformation,
                                                                                    currentGradientMicroDeformation );

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                                                    previousDeformationGradient, previousMicroDeformation,
                                                                                    previousGradientMicroDeformation );

    class hydraMock : public tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticityOptimization{

        public:

            using tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticityOptimization::hydraMicromorphicElastoPlasticityOptimization;

            void public_setInitializeUnknownVector( const bool value ){ setInitializeUnknownVector( value ); }

            void public_setUnknownVector( const tardigradeSolverTools::floatVector &value ){ updateUnknownVector( value ); }

    };

//    std::cout << "hydra optimize evaluate\n";
//    double temperature = 293.15;
//    double previousTemperature = 293.15;
//    hydraMock hydra( time[ 0 ], time[ 1 ],
//                         temperature,                     previousTemperature,
//                         currentDeformationGradient,      previousDeformationGradient,
//                         currentMicroDeformation,         previousMicroDeformation,
//                         currentGradientMicroDeformation, previousGradientMicroDeformation,
//                         { }, { },
//                         SDVS, fparams, 2, 15, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
//    try{
//    hydra.evaluate( );
//    }catch(std::exception &e){tardigradeErrorTools::printNestedExceptions(e); throw;}
//
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

    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );

    BOOST_TEST( PK2_result == PK2_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_result == SIGMA_answer, CHECK_PER_ELEMENT );

    BOOST_TEST( M_result == M_answer, CHECK_PER_ELEMENT );

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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     5e-4, 1e-3 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 5e-5, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_6, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.0;//0.96875;

    std::vector< double > _time = { 0.420361, 0.000005 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = { 2.000000, 10000.000000, 0.000000, 2.000000, 3.192203, -31.678450, 2.000000, 10000.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 696.441593, 126.713800, 5.000000, -18.674980, -37.817315, 15.177654, -24.071197, -5.861821, 11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 792.523471, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, -37.817315, -5.861821, 0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] = { {  0.018563, 0.001159, -0.007310 },
{ 0.002328, 0.017604, -0.004804 },
{ 0.006537, 0.003996, -0.041366 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  0.018304, 0.001129, -0.007037 },
{ 0.002247, 0.017494, -0.004474 },
{ 0.006155, 0.003807, -0.040882 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.023783, 0.000002, -0.000115, 0.000005, 0.023785, -0.000059, 0.000141, 0.000079, -0.051489 };

    double previous_phi[ 9 ] = { 0.023937, 0.000002, -0.000111, 0.000006, 0.023939, -0.000055, 0.000135, 0.000074, -0.051848 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  -0.000011,  0.000009,  0.000005 },
                                           {   0.000008,  0.000004, -0.000035 },
                                           {  -0.000026, -0.000019,  0.000534 },
                                           {  -0.000009,  0.000002,  0.000036 },
                                           {  -0.000018,  0.000009,  0.000000 },
                                           {   0.000010, -0.000039, -0.000321 },
                                           {   0.000019,  0.000012, -0.000428 },
                                           {  -0.000001,  0.000056,  0.000251 },
                                           {  -0.000012,  0.000011, -0.000006 } };

    double previous_grad_phi[ 9 ][ 3 ] = { { -0.000011,  0.000010,  0.000005 },
                                           {  0.000008,  0.000004, -0.000033 },
                                           { -0.000024, -0.000018,  0.000532 },
                                           { -0.000009,  0.000002,  0.000034 },
                                           { -0.000017,  0.000009,  0.000000 },
                                           {  0.000010, -0.000036, -0.000322 },
                                           {  0.000017,  0.000011, -0.000428 },
                                           { -0.000001,  0.000052,  0.000252 },
                                           { -0.000011,  0.000010, -0.000006 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = {  0.015840, 0.001312, -0.000720, 0.001311, 0.015914, -0.000008, -0.000719, -0.000024, -0.030964, 0.015842, 0.001313, -0.000542, 0.001314, 0.015916, 0.000016, -0.000648, 0.000039, -0.030968, 0.000000, 0.000000, 0.000001, -0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.000021, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, -0.000013, 0.000000, -0.000000, 0.000019, 0.000000, -0.000000, -0.000012, -0.000000, 0.000000, -0.000000, -0.000000, 21.525070, 0.000000, -0.000000, -0.000000, -0.000000, 0.062018, -0.000000, -0.000000, 0.000000 };

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

    tardigradeSolverTools::floatVector PK2_answer = { -3.346, 0.131246, 0.280253, 0.0873976, -3.5572, 0.0595398, -0.256392, -0.286154, -6.42817 };

    tardigradeSolverTools::floatVector SIGMA_answer = { -3.33216, 0.103392, 0.010602, 0.103392, -3.53612, -0.110036, 0.010602, -0.110036, -6.38187 };

    tardigradeSolverTools::floatVector M_answer = { -0.00800278, 0.00583708, -0.0203863, -0.00644166, -0.0131887, 0.00775079, 0.0139905, -0.000581907, -0.00996674, 0.00665793, 0.00321173, -0.0151118, 0.00145756, 0.00678371, -0.0314114, 0.00888833, 0.0425504, 0.00949894, 0.000955772, -0.0265837, 0.455456, 0.0272351, 0.00082456, -0.273899, -0.375001, 0.221355, -0.00758327 };

    tardigradeSolverTools::floatVector SDVS_answer = { 0.0163016, 0.00136103, -0.000712161, 0.00135994, 0.0162714, -5.99447e-05, -0.000709997, -7.84994e-05, -0.0317442, 0.0163033, 0.00136228, -0.000525351, 0.00136339, 0.0162739, -2.93465e-05, -0.000630233, -1.14368e-05, -0.0317484, -8.95455e-12, 5.15936e-10, 9.87075e-07, -8.49994e-10, 4.14015e-10, 2.80928e-08, -3.08557e-08, -2.37692e-08, 2.16708e-05, 1.63627e-10, -5.0658e-11, 2.35437e-08, 1.33409e-09, -4.52712e-09, -3.23167e-08, 9.89621e-09, -4.59031e-08, -1.33538e-05, -2.32633e-08, -1.69075e-08, 1.9579e-05, 6.14424e-10, -6.61986e-08, -1.22996e-05, -1.24056e-09, 3.90059e-09, 4.5542e-08, 2.11858e-20, 196.591, 3.4893e-22, 1.97791e-21, 0, 0, 0.0636232, 0, 0, 0 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeSolverTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeSolverTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
                                                                                    currentDeformationGradient, currentMicroDeformation,
                                                                                    currentGradientMicroDeformation );

    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
                                                                                    previousDeformationGradient, previousMicroDeformation,
                                                                                    previousGradientMicroDeformation );

    class hydraMock : public tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticityOptimization{

        public:

            using tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticityOptimization::hydraMicromorphicElastoPlasticityOptimization;

            void public_setInitializeUnknownVector( const bool value ){ setInitializeUnknownVector( value ); }

            void public_setUnknownVector( const tardigradeSolverTools::floatVector &value ){ updateUnknownVector( value ); }

    };

//    std::cout << "hydra optimize evaluate\n";
//    double temperature = 293.15;
//    double previousTemperature = 293.15;
//    hydraMock hydra( time[ 0 ], time[ 1 ],
//                         temperature,                     previousTemperature,
//                         currentDeformationGradient,      previousDeformationGradient,
//                         currentMicroDeformation,         previousMicroDeformation,
//                         currentGradientMicroDeformation, previousGradientMicroDeformation,
//                         { }, { },
//                         SDVS, fparams, 2, 15, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
//    try{
//    hydra.evaluate( );
//    }catch(std::exception &e){tardigradeErrorTools::printNestedExceptions(e); throw;}
//
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

//    std::cout << "PK2  : "; for ( auto v = PK2_result.begin( );   v != PK2_result.end( );   v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SIGMA: "; for ( auto v = SIGMA_result.begin( ); v != SIGMA_result.end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "M    : "; for ( auto v = M_result.begin( );     v != M_result.end( );     v++ ){ std::cout << *v << ", "; } std::cout << "\n";
//    std::cout << "SDVS : "; for ( auto v = SDVS.begin( );         v != SDVS.end( );         v++ ){ std::cout << *v << ", "; } std::cout << "\n";

    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );

    BOOST_TEST( PK2_answer ==  PK2_result, CHECK_PER_ELEMENT );

    BOOST_TEST( SIGMA_answer == SIGMA_result, CHECK_PER_ELEMENT );

    BOOST_TEST( M_answer == M_result, CHECK_PER_ELEMENT );

    //Test the Jacobians
    PK2_result.clear();
    SIGMA_result.clear();
    M_result.clear();
    ADD_TERMS.clear();

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
//    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 5e-4, 1e-3 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     5e-4, 1e-3 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 5e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 5e-5, 1e-5 ) );

}
////BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_7, * boost::unit_test::tolerance( 5e-4 ) ){
////    /*!
////     * Test the evaluation of the constitutive model.
////     *
////     */
////
////    //Initialize the time
////
////    double s = 1.0;//0.96875;
////
////    std::vector< double > _time = { 0.420000, 0.01 };
////
////    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };
////
////    double temperature = 293.15;
////
////    double previousTemperature = 293.15;
////
////    //Initialize the material parameters
////    std::vector< double > fparams = { 4.000000, 100000000.000000, 0.000000, 1e-2, 0.1,
////                                      4.000000, 3.192203, -316.78450, 1e-2, 0.1,
////                                      4.000000, 100000000.000000, 0.000000, 1e-2, 0.1,
////                                      2.000000, 0.000000, 0.000000,
////                                      2.000000, 0.000000, 0.000000,
////                                      2.000000, 0.000000, 0.000000,
////                                      2.000000, 0.000000, 0.000000,
////                                      2.000000, 0.000000, 0.000000,
////                                      2.000000, 0.000000, 0.000000,
////                                      2.000000, 696.441593, 126.713800,
////                                      5.000000, -18.674980, -37.817315, 15.177654, -24.071197, -5.861821,
////                                      11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 792.523471, 0.000000, 0.000000, 0.000000, 0.000000,
////                                      2.000000, -37.817315, -5.861821,
////                                      0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };
////
////    //Initialize the gradient of the macro displacement
////    double _current_grad_u[ 3 ][ 3 ] = { {  0.026125,  0.026303, -0.011998 },
////                                         {  0.028529, -0.008631, -0.014277 },
////                                         { -0.032270, -0.003117, -0.015042 } };
////
////    double previous_grad_u[ 3 ][ 3 ] = { {  0.011685, -0.000193, -0.003084 },
////                                         { -0.000166,  0.011788,  0.005098 },
////                                         {  0.000382, -0.000711, -0.025964 } };
////
////    double current_grad_u[ 3 ][ 3 ];
////    for ( unsigned int i = 0; i < 3; i++ ){
////        for ( unsigned int j = 0; j < 3; j++ ){
////            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
////        }
////    }
////
////    //Initialize the micro displacement
////    double _current_phi[ 9 ] = { 0.025131, 0.000168, -0.015482, 0.000148, 0.020700, 0.001777, -0.025419, 0.002269, -0.049488 };
////
////    double previous_phi[ 9 ] = { 0.025320, -0.000000, -0.000123, -0.000000, 0.025319, 0.000204, 0.000098, -0.000163, -0.055321 };
////
////    double current_phi[ 9 ];
////    for ( unsigned int i = 0; i < 9; i++ ){
////        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
////    }
////
////    //Initialize the gradient of the micro displacement
////    double _current_grad_phi[ 9 ][ 3 ] = { { -0.000023, -0.000037,  0.000070 },
////                                           { -0.000009, -0.000070,  0.000070 },
////                                           {  0.000520, -0.000009,  0.001549 },
////                                           { -0.000069,  0.000324, -0.000014 },
////                                           {  0.000045,  0.000075,  0.000046 },
////                                           { -0.000118,  0.000380,  0.000111 },
////                                           { -0.000367,  0.000078, -0.001478 },
////                                           {  0.000176, -0.000432, -0.000175 },
////                                           {  0.000058,  0.000047, -0.000009 } };
////
////    double previous_grad_phi[ 9 ][ 3 ] = { {  -0.000004,  0.000007,  0.000031 },
////                                           {   0.000000,  0.000000,  0.000001 },
////                                           {   0.000265,  0.000005,  0.000091 },
////                                           {  -0.000000,  0.000000,  0.000001 },
////                                           {  -0.000004,  0.000006,  0.000030 },
////                                           {   0.000004,  0.000262, -0.000148 },
////                                           {  -0.000210, -0.000002, -0.000070 },
////                                           {  -0.000002, -0.000211,  0.000113 },
////                                           {  -0.000001,  0.000001, -0.000019 } };
////
////    double current_grad_phi[ 9 ][ 3 ];
////    for ( unsigned int i = 0; i < 9; i++ ){
////        for ( unsigned int j = 0; j < 3; j++ ){
////            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
////        }
////    }
////
////    //Initialize the state variable vector
////    std::vector< double > SDVSDefault = {  0.003776, -0.000034, -0.000368, -0.000034, 0.003809, 0.000597, -0.000370, 0.000601, -0.007534, 0.003777, -0.000036, -0.000324, -0.000036, 0.003812, 0.000526, -0.000382, 0.000620, -0.007538, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000003, 0.000000, 0.000001, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.000003, -0.000002, 0.000002, 0.000000, 0.000001, 0.000000, 0.000002, -0.000001, -0.000000, 0.000000, -0.000000, 0.000000, 0.176535, 0.000000, -0.000000, 0.000000, 0.000000, 0.0014894, -0.000000, 0.000000, -0.000000 };
////
////    //Initialize the additional degree of freedom vectors
////    std::vector< double > current_ADD_DOF;
////    std::vector< std::vector< double > > current_ADD_grad_DOF;
////
////    std::vector< double > previous_ADD_DOF;
////    std::vector< std::vector< double > > previous_ADD_grad_DOF;
////
////    //Initialize the stress measures
////    std::vector< double > PK2_result( 9, 0 );
////
////    std::vector< double > SIGMA_result( 9, 0 );
////
////    std::vector< double > M_result( 27, 0 );
////
////    //Initialize the additional terms vector
////    std::vector< std::vector< double > > ADD_TERMS;
////
////    //Initialize the output message string
////    std::string output_message;
////
////    tardigradeSolverTools::floatVector PK2_answer = { -3.37057, 0.0866543, 0.0323353, 0.0818984, -3.36447, -0.00167251, 0.0276005, -0.0735036, -3.36448 };
////
////    tardigradeSolverTools::floatVector SIGMA_answer = { -3.34476, 0.066115, 0.0253977, 0.066115, -3.33988, -0.031869, 0.0253977, -0.031869, -3.47732 };
////
////    tardigradeSolverTools::floatVector M_answer = { 0.00222842, -0.0157762, 0.00575456, 0.00895602, -0.00581722, 0.00508999, -0.00693765, -0.00227177, 0.0143212, 0.00293983, 0.00706979, 0.00630251, 0.00209025, -0.00530594, 0.00931891, -0.00343587, -0.0104565, -0.0224589, 0.00330374, 0.00677226, 0.488466, -0.00499722, 0.00440454, -0.3879, -0.441158, 0.356116, 0.000501602 };
////
////    tardigradeSolverTools::floatVector SDVS_answer = { 0.0107416, 0.00643479, 0.00142293, 0.00643734, 0.0112673, -0.00235179, 0.00153702, -0.00250549, -0.0212581, 0.0107496, 0.00642886, 0.00149239, 0.00642612, 0.0112775, -0.00236504, 0.00165881, -0.00264259, -0.0212743, 1.42267e-07, -8.46723e-08, -2.26279e-06, -4.38303e-08, -8.23814e-08, 2.86685e-06, 7.82471e-07, -1.08144e-07, 1.15626e-05, 7.17795e-08, 3.79839e-08, 2.75112e-06, -1.58217e-07, 1.34147e-07, -2.35395e-06, -1.90528e-07, 7.9154e-07, -5.86987e-06, 6.55041e-07, 1.58951e-07, 9.07531e-06, 3.67226e-08, 6.23863e-07, -4.11633e-06, 2.16544e-08, -5.33505e-08, 4.554e-06, -1.7134e-23, 8.8538, -8.51183e-24, -1.71876e-23, -5.08624e-24, -3.9568e-24, 0.109811, 1.84609e-23, 7.51148e-24, -5.81825e-25 };
////
////    cleanAnswer( SDVS_answer );
////
////    std::vector< double > SDVS = SDVSDefault;
////
////    // Explore continuation approach
////
////    std::cout << "\n\nDIFFICULT 5\n\n";
////    fparams[ 7 ] = -3.1678450;
////    fparams[ 8 ] = 1e-2;
////    fparams[ 9 ] = 0.8521964817713661;
////
////    tardigradeSolverTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;
////
////    tardigradeSolverTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;
////
////    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( current_grad_u, current_phi, current_grad_phi,
////                                                                                    currentDeformationGradient, currentMicroDeformation,
////                                                                                    currentGradientMicroDeformation );
////
////    tardigradeMicromorphicElastoPlasticity::assembleFundamentalDeformationMeasures( previous_grad_u, previous_phi, previous_grad_phi,
////                                                                                    previousDeformationGradient, previousMicroDeformation,
////                                                                                    previousGradientMicroDeformation );
////
////    class hydraMock : public tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticity{
////
////        public:
////
////            using tardigradeMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticity::hydraMicromorphicElastoPlasticity;
////
////            void public_setInitializeUnknownVector( const bool value ){ setInitializeUnknownVector( value ); }
////
////            void public_setUnknownVector( const tardigradeSolverTools::floatVector &value ){ updateUnknownVector( value ); }
////
////    };
////
////    try{
////        hydraMock hydra( time[ 0 ], time[ 1 ],
////                         temperature,                     previousTemperature,
////                         currentDeformationGradient,      previousDeformationGradient,
////                         currentMicroDeformation,         previousMicroDeformation,
////                         currentGradientMicroDeformation, previousGradientMicroDeformation,
////                         { }, { },
////                         SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
////
////        hydra.evaluate( );
////
////        std::cout << "X converged:\n"; for ( auto v = hydra.getUnknownVector( )->begin( ); v != hydra.getUnknownVector( )->end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
////
////        fparams[ 7 ] = -31.678450;
////        fparams[ 8 ] = 1e-2;
////        fparams[ 9 ] = 0.8521964817713661;
////
////        hydraMock hydra2( time[ 0 ], time[ 1 ],
////                          temperature,                     previousTemperature,
////                          currentDeformationGradient,      previousDeformationGradient,
////                          currentMicroDeformation,         previousMicroDeformation,
////                          currentGradientMicroDeformation, previousGradientMicroDeformation,
////                          { }, { },
////                          SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
////
////        hydra2.public_setInitializeUnknownVector( false );
////        hydra2.public_setUnknownVector( *hydra.getUnknownVector( ) );
////
////        hydra2.evaluate( );
////
////        std::cout << "X converged:\n"; for ( auto v = hydra2.getUnknownVector( )->begin( ); v != hydra2.getUnknownVector( )->end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
////
////        fparams[ 7 ] = -63.356899999999996;
////        fparams[ 8 ] = 1e-2;
////        fparams[ 9 ] = 0.8521964817713661;
////
////        hydraMock hydra3( time[ 0 ], time[ 1 ],
////                          temperature,                     previousTemperature,
////                          currentDeformationGradient,      previousDeformationGradient,
////                          currentMicroDeformation,         previousMicroDeformation,
////                          currentGradientMicroDeformation, previousGradientMicroDeformation,
////                          { }, { },
////                          SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
////
////        hydra3.public_setInitializeUnknownVector( false );
////        hydra3.public_setUnknownVector( *hydra2.getUnknownVector( ) );
////
////        hydra3.evaluate( );
////
////        std::cout << "X converged:\n"; for ( auto v = hydra3.getUnknownVector( )->begin( ); v != hydra3.getUnknownVector( )->end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
////
////        fparams[ 7 ] = -126.71379999999999;
////        fparams[ 8 ] = 1e-2;
////        fparams[ 9 ] = 0.8521964817713661;
////
////        hydraMock hydra4( time[ 0 ], time[ 1 ],
////                          temperature,                     previousTemperature,
////                          currentDeformationGradient,      previousDeformationGradient,
////                          currentMicroDeformation,         previousMicroDeformation,
////                          currentGradientMicroDeformation, previousGradientMicroDeformation,
////                          { }, { },
////                          SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
////
////        hydra4.public_setInitializeUnknownVector( false );
////        hydra4.public_setUnknownVector( *hydra3.getUnknownVector( ) );
////
////        hydra4.evaluate( );
////
////        std::cout << "X converged:\n"; for ( auto v = hydra4.getUnknownVector( )->begin( ); v != hydra4.getUnknownVector( )->end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
////
////        fparams[ 7 ] = -253.42759999999998;
////        fparams[ 8 ] = 1e-2;
////        fparams[ 9 ] = 0.8521964817713661;
////
////        hydraMock hydra5( time[ 0 ], time[ 1 ],
////                          temperature,                     previousTemperature,
////                          currentDeformationGradient,      previousDeformationGradient,
////                          currentMicroDeformation,         previousMicroDeformation,
////                          currentGradientMicroDeformation, previousGradientMicroDeformation,
////                          { }, { },
////                          SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
////
////        hydra5.public_setInitializeUnknownVector( false );
////        hydra5.public_setUnknownVector( *hydra4.getUnknownVector( ) );
////
////        hydra5.evaluate( );
////
////        std::cout << "X converged:\n"; for ( auto v = hydra5.getUnknownVector( )->begin( ); v != hydra5.getUnknownVector( )->end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
////
////        fparams[ 7 ] = -316.7845;
////        fparams[ 8 ] = 1e-2;
////        fparams[ 9 ] = 0.8521964817713661;
////
////        hydraMock hydra6( time[ 0 ], time[ 1 ],
////                          temperature,                     previousTemperature,
////                          currentDeformationGradient,      previousDeformationGradient,
////                          currentMicroDeformation,         previousMicroDeformation,
////                          currentGradientMicroDeformation, previousGradientMicroDeformation,
////                          { }, { },
////                          SDVS, fparams, 2, 10, 3, 45, 1e-9, 1e-9, 20, 10, 1e-4, true, 0 );
////
////        hydra6.public_setInitializeUnknownVector( false );
////        hydra6.public_setUnknownVector( *hydra5.getUnknownVector( ) );
////
////        hydra5.evaluate( );
////
////        std::cout << "X converged:\n"; for ( auto v = hydra6.getUnknownVector( )->begin( ); v != hydra6.getUnknownVector( )->end( ); v++ ){ std::cout << *v << ", "; } std::cout << "\n";
////
////    }catch(std::exception &e){tardigradeErrorTools::printNestedExceptions(e);throw;}
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////    int errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
//////                                                                                  current_grad_u,  current_phi,  current_grad_phi,
//////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                                  SDVS,
//////                                                                                  current_ADD_DOF,  current_ADD_grad_DOF,
//////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                                  PK2_result, SIGMA_result, M_result,
//////                                                                                  ADD_TERMS,
//////                                                                                  output_message
//////                                                                                  );
//////
//////    BOOST_CHECK( errorCode == 0 );
//////
//////    if ( errorCode != 0 ){
//////        std::cout << "output_message:\n" << output_message << "\n";
//////    }
//////
//////    BOOST_TEST( SDVS_answer == SDVS, CHECK_PER_ELEMENT );
//////
//////    BOOST_TEST( PK2_answer ==  PK2_result, CHECK_PER_ELEMENT );
//////
//////    BOOST_TEST( SIGMA_answer == SIGMA_result, CHECK_PER_ELEMENT );
//////
//////    BOOST_TEST( M_answer == M_result, CHECK_PER_ELEMENT );
////
//////    //Test the Jacobians
//////    PK2_result.clear();
//////    SIGMA_result.clear();
//////    M_result.clear();
//////    ADD_TERMS.clear();
//////
//////    variableMatrix result_dPK2dGradU(      9, variableVector(  9, 0 ) );
//////
//////    variableMatrix result_dPK2dPhi(        9, variableVector(  9, 0 ) );
//////
//////    variableMatrix result_dPK2dGradPhi(    9, variableVector( 27, 0 ) );
//////
//////    variableMatrix result_dSIGMAdGradU(    9, variableVector(  9, 0 ) );
//////
//////    variableMatrix result_dSIGMAdPhi(      9, variableVector(  9, 0 ) );
//////
//////    variableMatrix result_dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );
//////
//////    variableMatrix result_dMdGradU(       27, variableVector(  9, 0 ) );
//////
//////    variableMatrix result_dMdPhi(         27, variableVector(  9, 0 ) );
//////
//////    variableMatrix result_dMdGradPhi(     27, variableVector( 27, 0 ) );
//////
//////    std::vector< variableMatrix > ADD_JACOBIANS;
//////
//////    SDVS = SDVSDefault;
//////
//////    errorCode  = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams,
//////                                                                               current_grad_u,  current_phi,  current_grad_phi,
//////                                                                               previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                               SDVS,
//////                                                                               current_ADD_DOF,  current_ADD_grad_DOF,
//////                                                                               previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                               PK2_result, SIGMA_result, M_result,
//////                                                                               result_dPK2dGradU, result_dPK2dPhi, result_dPK2dGradPhi,
//////                                                                               result_dSIGMAdGradU, result_dSIGMAdPhi, result_dSIGMAdGradPhi,
//////                                                                               result_dMdGradU, result_dMdPhi, result_dMdGradPhi,
//////                                                                               ADD_TERMS, ADD_JACOBIANS, output_message
//////                                                                             );
//////
//////    BOOST_CHECK( errorCode <= 0 );
//////
//////    if ( errorCode != 0 ){
//////        std::cout << "output_message:\n" << output_message << "\n";
//////    }
//////
//////    BOOST_TEST( PK2_result == PK2_answer, CHECK_PER_ELEMENT );
//////
//////    BOOST_TEST( SIGMA_result == SIGMA_answer, CHECK_PER_ELEMENT );
//////
//////    BOOST_TEST( M_result == M_answer, CHECK_PER_ELEMENT );
//////
//////    variableMatrix dPK2dGradU(      9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dPK2dPhi(        9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dPK2dGradPhi(    9, variableVector( 27, 0 ) );
//////
//////    variableMatrix dSIGMAdGradU(    9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dSIGMAdPhi(      9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dSIGMAdGradPhi(  9, variableVector( 27, 0 ) );
//////
//////    variableMatrix dMdGradU(       27, variableVector(  9, 0 ) );
//////
//////    variableMatrix dMdPhi(         27, variableVector(  9, 0 ) );
//////
//////    variableMatrix dMdGradPhi(     27, variableVector( 27, 0 ) );
//////
//////    variableMatrix dFpdGradU(          9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dFpdPhi(            9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dFpdGradPhi(        9, variableVector( 27, 0 ) );
//////
//////    variableMatrix dChipdGradU(        9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dChipdPhi(          9, variableVector(  9, 0 ) );
//////
//////    variableMatrix dChipdGradPhi(      9, variableVector( 27, 0 ) );
//////
//////    variableMatrix dGradChipdGradU(   27, variableVector(  9, 0 ) );
//////
//////    variableMatrix dGradChipdPhi(     27, variableVector(  9, 0 ) );
//////
//////    variableMatrix dGradChipdGradPhi( 27, variableVector( 27, 0 ) );
//////
//////    variableType eps = 1e-6;
//////
//////    for ( unsigned int i = 0; i < 9; i++ ){
//////
//////        variableVector delta( 9, 0 );
//////
//////        unsigned int row = i / 3;
//////
//////        unsigned int col = i % 3;
//////
//////        delta[ i ] = eps * std::fabs( current_grad_u[ row ][ col ] ) + eps;
//////
//////        variableType current_grad_u_p[ 3 ][ 3 ];
//////        variableType current_grad_u_m[ 3 ][ 3 ];
//////
//////        for ( unsigned int _i = 0; _i < 3; _i++ ){
//////            for ( unsigned int _j = 0; _j < 3; _j++ ){
//////                current_grad_u_p[ _i ][ _j ] = current_grad_u[ _i ][ _j ] + delta[ 3 * _i + _j ];
//////                current_grad_u_m[ _i ][ _j ] = current_grad_u[ _i ][ _j ] - delta[ 3 * _i + _j ];
//////            }
//////        }
//////
//////        variableVector PK2_p,   PK2_m;
//////        variableVector SIGMA_p, SIGMA_m;
//////        variableVector M_p,     M_m;
//////        variableVector SDVS_p = SDVSDefault;
//////        variableVector SDVS_m = SDVSDefault;
//////
//////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_p, current_phi, current_grad_phi,
//////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
//////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                                  PK2_p, SIGMA_p, M_p,
//////                                                                                  ADD_TERMS, output_message
//////                                                                                );
//////
//////        BOOST_CHECK( errorCode <= 0 );
//////
//////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u_m, current_phi, current_grad_phi,
//////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
//////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                                  PK2_m, SIGMA_m, M_m,
//////                                                                                  ADD_TERMS, output_message
//////                                                                                );
//////
//////        BOOST_CHECK( errorCode <= 0 );
//////
//////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
//////
//////            dPK2dGradU[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
//////
//////            dSIGMAdGradU[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
//////
//////            dMdGradU[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 9; j++ ){
//////
//////            dFpdGradU[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 9; j++ ){
//////
//////            dChipdGradU[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 27; j++ ){
//////
//////            dGradChipdGradU[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////    }
//////
//////    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dGradU   ) == tardigradeVectorTools::appendVectors( result_dPK2dGradU ), CHECK_PER_ELEMENT );
//////    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdGradU ) == tardigradeVectorTools::appendVectors( result_dSIGMAdGradU ), CHECK_PER_ELEMENT );
//////    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdGradU     ) == tardigradeVectorTools::appendVectors( result_dMdGradU ), CHECK_PER_ELEMENT );
//////
//////    BOOST_TEST( tardigradeVectorTools::appendVectors( dFpdGradU       ) == tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 0 ] ), CHECK_PER_ELEMENT );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-5, 1e-5 ) );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-5, 1e-5 ) );
//////
//////    for ( unsigned int i = 0; i < 9; i++ ){
//////
//////        variableVector delta( 9, 0 );
//////
//////        delta[ i ] = eps * std::fabs( current_phi[ i ] ) + eps;
//////
//////        variableType current_phi_p[ 9 ];
//////        variableType current_phi_m[ 9 ];
//////
//////        for ( unsigned int _i = 0; _i < 3; _i++ ){
//////            for ( unsigned int _j = 0; _j < 3; _j++ ){
//////                current_phi_p[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] + delta[ 3 * _i + _j ];
//////                current_phi_m[ 3 * _i + _j ] = current_phi[ 3 * _i + _j ] - delta[ 3 * _i + _j ];
//////            }
//////        }
//////
//////        variableVector PK2_p,   PK2_m;
//////        variableVector SIGMA_p, SIGMA_m;
//////        variableVector M_p,     M_m;
//////        variableVector SDVS_p = SDVSDefault;
//////        variableVector SDVS_m = SDVSDefault;
//////
//////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_p, current_grad_phi,
//////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
//////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                                  PK2_p, SIGMA_p, M_p,
//////                                                                                  ADD_TERMS, output_message
//////                                                                                );
//////
//////        BOOST_CHECK( errorCode <= 0 );
//////
//////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi_m, current_grad_phi,
//////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
//////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                                  PK2_m, SIGMA_m, M_m,
//////                                                                                  ADD_TERMS, output_message
//////                                                                                );
//////
//////        BOOST_CHECK( errorCode <= 0 );
//////
//////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
//////
//////            dPK2dPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
//////
//////            dSIGMAdPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
//////
//////            dMdPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 9; j++ ){
//////
//////            dFpdPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 9; j++ ){
//////
//////            dChipdPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 27; j++ ){
//////
//////            dGradChipdPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////    }
//////
//////    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2dPhi   ) == tardigradeVectorTools::appendVectors( result_dPK2dPhi ), CHECK_PER_ELEMENT );
//////    BOOST_TEST( tardigradeVectorTools::appendVectors( dSIGMAdPhi ) == tardigradeVectorTools::appendVectors( result_dSIGMAdPhi ), CHECK_PER_ELEMENT );
//////    BOOST_TEST( tardigradeVectorTools::appendVectors( dMdPhi     ) == tardigradeVectorTools::appendVectors( result_dMdPhi ), CHECK_PER_ELEMENT );
//////
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdPhi )      , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 1 ] ), 5e-5, 1e-5 ) );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-5, 1e-5 ) );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 7 ] ), 1e-6, 1e-6 ) );
//////
//////    for ( unsigned int i = 0; i < 27; i++ ){
//////
//////        variableVector delta( 27, 0 );
//////
//////        unsigned int row = i / 9;
//////
//////        unsigned int col = i % 9;
//////
//////        delta[ i ] = eps * std::fabs( current_grad_phi[ row ][ col ] ) + eps;
//////
//////        variableType current_grad_phi_p[ 9 ][ 3 ];
//////        variableType current_grad_phi_m[ 9 ][ 3 ];
//////
//////        for ( unsigned int _i = 0; _i < 9; _i++ ){
//////            for ( unsigned int _j = 0; _j < 3; _j++ ){
//////                current_grad_phi_p[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] + delta[ 3 * _i + _j ];
//////                current_grad_phi_m[ _i ][ _j ] = current_grad_phi[ _i ][ _j ] - delta[ 3 * _i + _j ];
//////            }
//////        }
//////
//////        variableVector PK2_p,   PK2_m;
//////        variableVector SIGMA_p, SIGMA_m;
//////        variableVector M_p,     M_m;
//////        variableVector SDVS_p = SDVSDefault;
//////        variableVector SDVS_m = SDVSDefault;
//////
//////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_p,
//////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                                  SDVS_p, current_ADD_DOF, current_ADD_grad_DOF,
//////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                                  PK2_p, SIGMA_p, M_p,
//////                                                                                  ADD_TERMS, output_message
//////                                                                                );
//////
//////        BOOST_CHECK( errorCode <= 0 );
//////
//////        errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model( time, fparams, current_grad_u, current_phi, current_grad_phi_m,
//////                                                                                  previous_grad_u, previous_phi, previous_grad_phi,
//////                                                                                  SDVS_m, current_ADD_DOF, current_ADD_grad_DOF,
//////                                                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
//////                                                                                  PK2_m, SIGMA_m, M_m,
//////                                                                                  ADD_TERMS, output_message
//////                                                                                );
//////
//////        BOOST_CHECK( errorCode <= 0 );
//////
//////        for ( unsigned int j = 0; j < PK2_p.size( ); j++ ){
//////
//////            dPK2dGradPhi[ j ][ i ] = ( PK2_p[ j ] - PK2_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < SIGMA_p.size( ); j++ ){
//////
//////            dSIGMAdGradPhi[ j ][ i ] = ( SIGMA_p[ j ] - SIGMA_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < M_p.size( ); j++ ){
//////
//////            dMdGradPhi[ j ][ i ] = ( M_p[ j ] - M_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 9; j++ ){
//////
//////            dFpdGradPhi[ j ][ i ] = ( SDVS_p[ j ] - SDVS_m[ j ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 9; j++ ){
//////
//////            dChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 9 ] - SDVS_m[ j + 9 ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////        for ( unsigned int j = 0; j < 27; j++ ){
//////
//////            dGradChipdGradPhi[ j ][ i ] = ( SDVS_p[ j + 18 ] - SDVS_m[ j + 18 ] ) / ( 2 * delta[ i ] );
//////
//////        }
//////
//////    }
//////
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   5e-4, 1e-3 ) );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 5e-4, 1e-3 ) );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     5e-4, 1e-3 ) );
//////
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 5e-5, 1e-5 ) );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 5e-5, 1e-5 ) );
//////    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 5e-5, 1e-5 ) );
////
////}
