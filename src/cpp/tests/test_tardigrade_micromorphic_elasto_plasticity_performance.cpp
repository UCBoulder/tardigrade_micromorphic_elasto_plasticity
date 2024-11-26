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

    tardigradeMicromorphicTools::floatVector PK2_answer = { -1.40279, 0.0121485, 0.0266627, 0.00381226, -1.45011, 0.0418199, 0.0338264, 0.0864571, -2.91531 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -0.493708, 0.0121523, 0.0789793, 0.0121523, -0.517776, 0.106717, 0.0789793, 0.106717, -4.81709 };

    tardigradeMicromorphicTools::floatVector M_answer = { 8.45418e-05, -4.04181e-06, -0.000551314, -2.62816e-06, 9.35849e-05, 4.20583e-05, 0.000290068, -2.54552e-05, -6.72292e-05, -1.81165e-05, 6.55509e-07, 0.000106142, -1.94223e-06, -1.8193e-05, 8.36339e-06, -5.40107e-05, -4.57942e-06, 2.03912e-06, -0.000144554, 2.02458e-05, 0.00030473, 1.70095e-05, -0.000104047, 0.00037256, -0.000143114, -0.000161128, 7.22794e-05 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.00458246, 3.26705e-05, 0.000297186, 0.000112011, 0.00414933, 0.000768014, 0.000230621, 0.000356534, -0.00861812, 1.43254e-11, 1.68409e-13, 8.44243e-13, 2.5468e-13, 1.31286e-11, 1.89224e-12, 7.70904e-13, 1.41952e-12, -2.64072e-11, 2.79233e-17, 5.56405e-20, -7.87013e-18, -3.37087e-18, 1.51828e-17, -1.01737e-17, 2.41855e-17, -1.34546e-18, -5.24753e-17, -2.32027e-18, 1.54445e-17, 6.40194e-18, -3.58006e-19, -4.96139e-18, -4.77126e-18, -3.94101e-18, -1.629e-17, -4.58044e-17, -1.27547e-16, 1.15812e-17, 6.96373e-17, 2.14292e-17, -2.77694e-17, 1.05305e-16, 6.73349e-18, 1.54864e-17, 3.42911e-17, 0.170641, 4.09748e-26, 2.09674e-24, 0, 0, 0.0172535, 0, 5.60883e-25, 0, 1.61587e-23 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { 7.20108, -0.010419, 0.40178, -0.0312379, 7.31865, -0.501715, -0.443789, 0.641519, 6.95971 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { 7.1785, -0.0545498, -0.0987932, -0.0545498, 7.33568, 0.109286, -0.0987932, 0.109286, 6.81404 };

    tardigradeMicromorphicTools::floatVector M_answer = { -0.0263527, 0.00471758, -0.308025, -0.0137948, 0.00234765, 0.183845, 0.625938, -0.133865, 0.0583556, -0.043483, 0.0126769, -0.0571493, -0.020039, -0.0335559, 0.00404204, 0.180415, 0.0190195, -0.124223, 0.270528, -0.109253, -0.0708277, -0.0162357, 0.321918, 0.550908, 0.267257, -0.827036, -0.0448194 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { -0.00884271, 0.00720015, 0.0131995, 0.00722664, -0.0283276, -0.0133543, 0.0132149, -0.0134292, 0.0404522, -0.00913627, 0.00712621, 0.0127927, 0.00796593, -0.0286991, -0.0116354, 0.0123124, -0.0129925, 0.0411172, 1.54245e-05, 3.73484e-06, 6.78738e-06, -6.96842e-06, -2.49434e-07, -1.78553e-05, 2.20578e-05, 2.42514e-06, 1.63091e-06, -1.23956e-05, -2.43171e-06, -1.27993e-05, 4.98752e-06, -1.02922e-08, 2.05693e-05, -1.89819e-05, 1.07569e-06, -4.26493e-05, 4.306e-05, 1.37901e-05, 3.24546e-05, -1.70277e-05, -1.20242e-06, -8.77688e-05, -2.0412e-05, -3.72455e-06, -2.73567e-05, 1.75165e-24, 1.11041, 3.23707e-24, -4.92829e-24, 1.85212e-24, 4.6433e-25, 0.093562, 5.34969e-24, 1.98121e-23, -4.48109e-24 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { -3.46298, 0.177874, 0.0884305, 0.173157, -3.45115, -0.058573, 0.0844343, -0.131236, -3.09472 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -3.43632, 0.15487, 0.0802941, 0.15487, -3.42587, -0.0876043, 0.0802941, -0.0876043, -3.21697 };

    tardigradeMicromorphicTools::floatVector M_answer = { 0.00222873, -0.0157277, 0.00601891, 0.00893157, -0.00581551, 0.00489429, -0.00716736, -0.00212801, 0.0143433, 0.00293253, 0.00701827, 0.00604097, 0.0021024, -0.00528301, 0.00953699, -0.00324134, -0.0106345, -0.0225163, 0.00309185, 0.00690868, 0.488745, -0.00475145, 0.00424082, -0.387931, -0.44314, 0.357667, 0.00100115 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.0112824, 0.00600669, 0.00117447, 0.00600907, 0.0117805, -0.00209826, 0.00127799, -0.00224004, -0.0223169, 0.0112905, 0.00600083, 0.00125425, 0.00599828, 0.0117902, -0.00211992, 0.00138786, -0.00236492, -0.0223327, 1.30188e-07, -7.65556e-08, -1.97681e-06, -3.89793e-08, -7.85229e-08, 2.60755e-06, 7.89654e-07, -9.26123e-08, 1.26652e-05, 6.8574e-08, 3.35458e-08, 2.49962e-06, -1.48618e-07, 1.26265e-07, -2.12921e-06, -1.70498e-07, 7.98944e-07, -6.81775e-06, 6.71851e-07, 1.55138e-07, 1.01423e-05, 4.0342e-08, 6.43625e-07, -5.04676e-06, 2.37016e-08, -5.32782e-08, 4.04748e-06, 0, 8.24694, 1.20612e-25, 5.86763e-27, 0, 9.42791e-26, 0.107333, 0, 0, 0 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { -0.42083, 0.188719, -0.284488, 0.113652, -0.606198, 0.198919, 0.154624, -0.272384, -0.412977 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -0.49984, 0.07922, -0.0612982, 0.07922, -0.600274, -0.024636, -0.0612982, -0.024636, -0.592951 };

    tardigradeMicromorphicTools::floatVector M_answer = { -0.0112254, 0.00535196, 0.439085, -0.0777306, 0.0376894, -0.132622, -0.334233, 0.158964, 0.0606504, -0.0371522, -0.067423, -0.0265724, 0.27262, 0.0507418, 0.335067, 0.0709937, -0.374721, 0.0389545, 0.0377269, 0.0561764, 1.24871, -0.0416176, 0.0361598, 0.023124, -1.20441, -0.108903, 0.0370146 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.026465, 0.0278464, -0.0216844, 0.0278898, -0.00914431, -0.00885806, -0.0219855, -0.00894281, -0.0138897, 0.0261713, 0.0276253, -0.0208302, 0.0283241, -0.00897724, -0.00911191, -0.0228368, -0.00914724, -0.0137652, 1.91504e-05, 9.4537e-06, 6.63661e-05, 3.22352e-06, 1.01035e-05, 2.17514e-05, 1.32479e-05, 8.0438e-06, 4.98021e-05, 2.4575e-06, -6.58701e-06, 1.81752e-05, -8.00738e-07, -2.80077e-06, 5.03653e-06, 1.42835e-05, 7.20726e-06, 4.16155e-05, 1.20598e-05, 9.93543e-06, 5.18231e-05, 1.35541e-05, -1.7557e-06, 4.0168e-05, -1.8236e-05, -6.62803e-06, -7.06047e-05, 9.05459e-23, 5.67313, -9.38528e-24, -5.87456e-24, -6.34399e-24, -5.13721e-23, 0.107536, 9.57509e-23, 4.31927e-23, 9.86484e-23 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS = SDVSDefault;

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { -9.62479, 0.38928, 0.361314, 0.513209, -9.65339, -0.531695, -0.0320208, 0.105003, -8.98122 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -9.34451, 0.375142, 0.127635, 0.375142, -9.37008, -0.165708, 0.127635, -0.165708, -8.99023 };

    tardigradeMicromorphicTools::floatVector M_answer = { 0.0134268, -0.0543831, -0.0879468, 0.038986, 0.00690052, -0.0112154, 0.0777729, 0.0124799, 0.0411076, -0.00806668, -0.0164, -0.00616632, 0.0329602, -0.0166184, -0.0747552, 0.0103869, 0.062287, -0.0513045, -0.0024238, 0.00359776, 0.421367, -0.0065458, -0.00293546, -0.407405, -0.40262, 0.393485, -0.0191425 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { -0.00540907, 0.0201288, 0.00644233, 0.0201198, -0.00676143, -0.00849111, 0.00690441, -0.00908192, 0.0157183, -0.00540473, 0.0202057, 0.00563266, 0.0202134, -0.00695702, -0.00741807, 0.00634141, -0.00837407, 0.0158932, 3.54253e-06, 1.37056e-06, -6.61657e-06, -1.12583e-06, 8.89558e-08, 6.56662e-06, 5.86212e-06, -1.08602e-06, -2.05309e-05, -5.75537e-07, 5.87735e-07, 6.47761e-06, -2.62242e-06, -2.63729e-06, -7.29034e-06, -1.04981e-06, 4.21174e-06, 2.41893e-05, 5.82389e-06, -2.89086e-07, -2.83413e-05, -4.85834e-07, 4.38003e-06, 3.15114e-05, -8.90182e-07, 1.24124e-06, 1.37233e-05, 5.15022e-24, 116.919, 4.49973e-23, 3.5662e-22, 0, 0, 0.19808, 1.18479e-21, 0, 1.68029e-21 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { -3.346, 0.131246, 0.280253, 0.0873976, -3.5572, 0.0595398, -0.256392, -0.286154, -6.42817 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -3.33216, 0.103392, 0.010602, 0.103392, -3.53612, -0.110036, 0.010602, -0.110036, -6.38187 };

    tardigradeMicromorphicTools::floatVector M_answer = { -0.00800278, 0.00583708, -0.0203863, -0.00644166, -0.0131887, 0.00775079, 0.0139905, -0.000581907, -0.00996674, 0.00665793, 0.00321173, -0.0151118, 0.00145756, 0.00678371, -0.0314114, 0.00888833, 0.0425504, 0.00949894, 0.000955772, -0.0265837, 0.455456, 0.0272351, 0.00082456, -0.273899, -0.375001, 0.221355, -0.00758327 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.0163016, 0.00136103, -0.000712161, 0.00135994, 0.0162714, -5.99447e-05, -0.000709997, -7.84994e-05, -0.0317442, 0.0163033, 0.00136228, -0.000525351, 0.00136339, 0.0162739, -2.93465e-05, -0.000630233, -1.14368e-05, -0.0317484, -8.95455e-12, 5.15936e-10, 9.87075e-07, -8.49994e-10, 4.14015e-10, 2.80928e-08, -3.08557e-08, -2.37692e-08, 2.16708e-05, 1.63627e-10, -5.0658e-11, 2.35437e-08, 1.33409e-09, -4.52712e-09, -3.23167e-08, 9.89621e-09, -4.59031e-08, -1.33538e-05, -2.32633e-08, -1.69075e-08, 1.9579e-05, 6.14424e-10, -6.61986e-08, -1.22996e-05, -1.24056e-09, 3.90059e-09, 4.5542e-08, 2.11858e-20, 196.591, 3.4893e-22, 1.97791e-21, 0, 0, 0.0636232, 0, 0, 0 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }

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

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_7, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.0;//0.96875;

    std::vector< double > _time = { 0.42783, 2.2538e-11 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = {  2.000000, 10000.000000, 0.000000, 2.000000, 3.192203, -31.678450, 2.000000, 10000.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 696.441593, 126.713800, 5.000000, -18.674980, -37.817315, 15.177654, -24.071197, -5.861821, 11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 792.523471, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, -37.817315, -5.861821, 0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] = { {  0.012122, -0.000425, -0.003453 },
                                         { -0.000383,  0.012948,  0.004447 },
                                         {  0.000007,  0.000076, -0.027092 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  0.012122, -0.000425, -0.003453 },
                                         { -0.000383,  0.012948,  0.004447 },
                                         {  0.000007,  0.000076, -0.027092 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.023859, -0.000013, -0.000136, -0.000013, 0.023861, 0.000171, 0.000105, -0.000132, -0.051643 };

    double previous_phi[ 9 ] = { 0.023859, -0.000013, -0.000136, -0.000013, 0.023861, 0.000171, 0.000105, -0.000132, -0.051643 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  -0.000004, 0.000007,  0.000019 },
                                           { -0.000000, -0.000002,  0.000001 },
                                           {  0.000282,  0.000006,  0.000097 },
                                           { -0.000001, -0.000001,  0.000000 },
                                           { -0.000006,  0.000010,  0.000019 },
                                           {  0.000005,  0.000281, -0.000123 },
                                           { -0.000217, -0.000005, -0.000072 },
                                           { -0.000003, -0.000218,  0.000092 },
                                           { -0.000005,  0.000008, -0.000025 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {  -0.000004,  0.000007,  0.000019 },
                                           {  -0.000000, -0.000002,  0.000001 },
                                           {   0.000282,  0.000006,  0.000097 },
                                           {  -0.000001, -0.000001,  0.000000 },
                                           {  -0.000006,  0.000010,  0.000019 },
                                           {   0.000005,  0.000281, -0.000123 },
                                           {  -0.000217, -0.000005, -0.000072 },
                                           {  -0.000003, -0.000218,  0.000092 },
                                           {  -0.000005,  0.000008, -0.000025 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = {  0.004985, -0.000075, -0.000566, -0.000075, 0.005112, 0.000686, -0.000570, 0.000691, -0.010012, 0.004987, -0.000077, -0.000502, -0.000077, 0.005115, 0.000606, -0.000590, 0.000712, -0.010017, 0.000000, 0.000000, 0.000000, -0.000000, 0.000000, -0.000000, 0.000004, 0.000000, 0.000001, -0.000000, 0.000000, -0.000000, -0.000000, -0.000000, 0.000000, 0.000000, 0.000004, -0.000002, 0.000003, 0.000000, 0.000001, 0.000000, 0.000003, -0.000001, -0.000000, 0.000000, -0.000000, 0.000000, 0.000000, -0.000000, -0.000000, 0.000000, -0.000000, 0.019824, -0.000000, 0.000000, -0.000000 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { 0.188943, -0.0749734, -0.223576, -0.0766538, 0.350907, 0.308501, -0.349782, 0.468328, -5.0104 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { 0.145838, -0.0730737, -0.275922, -0.0730737, 0.30145, 0.373879, -0.275922, 0.373879, -5.00914 };

    tardigradeMicromorphicTools::floatVector M_answer = { -0.00304811, -0.000291436, 0.223069, -0.000555977, -0.00465674, 0.00388556, -0.169641, -0.00228998, -0.00469748, 0.005428, -0.00175762, 0.00470323, -0.000534611, 0.00769539, 0.222371, -0.00383504, -0.170343, 0.00729416, 0.0152528, 0.00078359, 0.0796262, -3.46937e-05, 0.0152673, -0.100393, -0.0582237, 0.0741628, -0.0211573 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.004985, -7.5e-05, -0.000566, -7.5e-05, 0.005112, 0.000686, -0.00057, 0.000691, -0.010012, 0.004987, -7.7e-05, -0.000502, -7.7e-05, 0.005115, 0.000606, -0.00059, 0.000712, -0.010017, 0, 0, 0, -0, 0, -0, 4e-06, 0, 1e-06, -0, 0, -0, -0, -0, 0, 0, 4e-06, -2e-06, 3e-06, 0, 1e-06, 0, 3e-06, -1e-06, -0, 0, -0, 0, 0, -0, -0, 0, -0, 0.019824, -0, 0, -0, };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }

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

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_8, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.00;//0.96875;

    std::vector< double > _time = { 0.42, 0.01 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = { 2.000000, 10000.000000, 0.000000, 2.000000, 3.192203, -3.1678450, 2.000000, 10000.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 2.000000, 696.441593, 126.713800, 5.000000, -18.674980, -37.817315, 15.177654, -24.071197, -5.861821, 11.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 792.523471, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, -37.817315, -5.861821, 0.500000, 0.500000, 0.500000, 0.000000, 0.000000 };

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] = { {  -0.015570,  0.015130, -0.001649 },
                                         {  -0.021229, -0.049745,  0.011857 },
                                         {   0.045788,  0.008404,  0.075197 } }; 

    double previous_grad_u[ 3 ][ 3 ] = { {  0.008494, -0.000070, -0.001709 },
                                         { -0.000054,  0.008905,  0.013283 },
                                         {  0.001531, -0.010064, -0.020614 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 0.026510, -0.000096, 0.000139, -0.000056, 0.022986, 0.000753, 0.000799, 0.000209, -0.053775 };

    double previous_phi[ 9 ] = { 0.025306, -0.000001, -0.000097, -0.000000, 0.025309, 0.000802, 0.000080, -0.000664, -0.055328 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  -0.000023, -0.000024,  0.000140 },
                                           {   0.000068,  0.000035, -0.000080 },
                                           {   0.000714,  0.000013, -0.000284 },
                                           {  -0.000042, -0.000041,  0.000061 },
                                           {  -0.000014, -0.000034,  0.000128 },
                                           {   0.000017,  0.000188, -0.000014 },
                                           {  -0.000727, -0.000001,  0.000297 },
                                           {   0.000041, -0.000171,  0.000065 },
                                           {  -0.000027, -0.000050,  0.000142 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {  -0.000000,  0.000001,  0.000064 },
                                           {   0.000002,  0.000000,  0.000002 },
                                           {   0.000273,  0.000028,  0.000034 },
                                           {   0.000001,  0.000000,  0.000003 },
                                           {   0.000000,  0.000001,  0.000044 },
                                           {   0.000034,  0.000051, -0.000277 },
                                           {  -0.000226, -0.000022, -0.000026 },
                                           {  -0.000027, -0.000048,  0.000207 },
                                           {   0.000001, -0.000004,  0.000001 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = {  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { 1.49261, -1.04731, 1.8083, 0.736858, -0.342405, 0.384324, 0.144605, 0.532215, 5.33617 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { 1.52344, -0.159191, 0.833773, -0.159191, -0.148749, 0.402651, 0.833773, 0.402651, 4.95921 };

    tardigradeMicromorphicTools::floatVector M_answer = { -0.0272444, 0.0579057, 0.573848, -0.0311418, -0.016227, 0.0310642, -0.544652, 0.0248917, -0.0423684, -0.0229199, 0.0353516, 0.0162531, -0.0429729, -0.0329841, 0.169744, -0.00473707, -0.139161, -0.0400953, 0.105146, -0.0627765, -0.214109, 0.0485745, 0.100853, -0.0215891, 0.20224, 0.046243, 0.0925418 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { -0.0134855, -0.00231109, 0.0181625, -0.00231944, -0.0420583, 0.00799272, 0.0183602, 0.00808636, 0.0619513, -0.0143809, -0.00298612, 0.0153429, -0.00301211, -0.041932, 0.00718093, 0.0181364, 0.00850675, 0.0627204, -2.5116e-05, -3.23837e-07, 9.19472e-06, -3.61913e-06, -2.11512e-06, 1.2987e-06, -5.67872e-05, -3.17903e-06, 2.15475e-05, -4.85412e-06, -2.55032e-06, 9.40791e-07, -2.08593e-07, -3.38863e-06, 9.69639e-07, -3.62395e-06, -2.07583e-05, 1.80779e-06, -6.17375e-05, -1.14634e-06, 2.39786e-05, 3.38009e-06, -1.97403e-05, 5.99492e-06, 2.53246e-05, 3.71247e-06, -1.01644e-05, 0, 7.98744, 0, 7.08167e-27, 3.33834e-27, 3.12771e-26, 0.130434, 0, 0, 3.85281e-27 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }

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

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_9, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.00;//0.96875;

    std::vector< double > _time = { 4.200000e-01, 1.000000e-02 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = {  2.000000e+00, 1.000000e+04, 1.000000e-08, 2.000000e+00, 3.192203e+00, -3.167845e+01, 2.000000e+00, 1.000000e+04, 1.000000e-08, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 6.964416e+02, 1.267138e+02, 5.000000e+00, -1.867498e+01, -3.781732e+01, 1.517765e+01, -2.407120e+01, -5.861821e+00, 1.100000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 7.925235e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, -3.781732e+01, -5.861821e+00, 5.000000e-01, 5.000000e-01, 5.000000e-01, 1.000000e-09, 1.000000e-09};

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] = { {  2.220177e-02,  2.070963e-03, -9.793878e-03 },
                                         { -3.984602e-04,  2.689540e-02,  7.465260e-04 },
                                         {  1.404226e-03, -4.049922e-04, -5.122192e-02 } }; 


    double previous_grad_u[ 3 ][ 3 ] = { { 9.635051e-03,  2.970750e-04, -3.179945e-03 },
                                         { 7.665442e-05,  1.449431e-02,  3.077166e-05 },
                                         { 7.028453e-04, -1.825419e-05, -2.755846e-02 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 2.497601e-02, 2.659689e-06, -4.770946e-05, -2.336211e-06, 2.499005e-02, 1.403223e-06, 3.510193e-05, -3.024517e-06, -5.423114e-02 };

    double previous_phi[ 9 ] = { 2.534622e-02, 3.685745e-08, -4.182825e-05, 2.084049e-07, 2.535998e-02, 7.350632e-07, 3.573346e-05, -5.753835e-07, -5.531934e-02 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  -1.217526e-05, -1.402528e-06,  1.005304e-05 },
                                           {  -1.179525e-06,  4.219067e-05,  2.351351e-06 },
                                           {   7.840170e-06,  2.398708e-06,  7.078417e-04 },
                                           {   6.467414e-07, -3.199255e-05, -2.096751e-06 },
                                           {  -1.496570e-05, -3.435084e-06,  8.911796e-06 },
                                           {   3.894465e-07,  2.790824e-05, -3.200066e-05 },
                                           {  -7.308141e-06, -1.751272e-06, -5.573903e-04 },
                                           {  -3.937788e-07, -2.461814e-05,  2.983384e-05 },
                                           {  -1.461687e-05, -6.245316e-06,  1.117437e-05 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {  -3.000087e-06,  6.812332e-07,  4.652051e-06 },
                                           {   3.191695e-09,  3.269879e-06, -1.156308e-07 },
                                           {   2.081071e-06,  6.435661e-07,  6.520715e-04 },
                                           {   9.524228e-08,  6.297548e-06,  1.813489e-07 },
                                           {  -6.615119e-06, -4.153298e-07,  6.714738e-06 },
                                           {   2.647349e-07,  1.345687e-05, -1.308597e-05 },
                                           {  -2.034080e-06, -5.633298e-07, -5.487315e-04 },
                                           {  -2.024323e-07, -1.091292e-05,  1.090765e-05 },
                                           {  -2.728086e-06,  3.352548e-07,  1.441408e-06 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = {  3.825127e-03, 7.937215e-05, -4.448362e-04, 7.937410e-05, 5.765822e-03, 1.779597e-06, -4.476565e-04, 1.790712e-06, -9.500946e-03, 3.826619e-03, 7.958918e-05, -3.906541e-04, 7.958451e-05, 5.765803e-03, 1.174487e-06, -4.608731e-04, 1.388454e-06, -9.502418e-03, 2.160358e-09, 6.361383e-10, 4.822937e-07, -1.980623e-10, -5.400300e-09, -5.080108e-09, 3.241278e-08, 1.110220e-08, 7.754323e-06, 5.718159e-10, 1.725566e-08, -6.359821e-09, -7.984061e-12, -5.438705e-11, 2.866011e-11, 3.969666e-09, 2.072214e-07, -1.392581e-07, 3.351035e-08, 9.891286e-09, 7.118100e-06, 3.520839e-09, 1.869695e-07, -1.280107e-07, -2.151697e-09, -5.816324e-10, -4.822508e-07, 9.929899e-24, 2.984149e-01, -1.159397e-23, -1.322973e-24, -3.047394e-25, -2.730157e-22, 1.888395e-02, -6.233582e-23, 1.669943e-23, 1.297502e-23 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { -3.5424, -0.036538, 0.175678, 0.052668, -3.49538, -0.0209274, -0.267199, 0.0249727, -4.11429 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -3.56824, 0.00601439, -0.0358833, 0.00601439, -3.53457, 0.0015479, -0.0358833, 0.0015479, -4.16335 };

    tardigradeMicromorphicTools::floatVector M_answer = { -0.00863623, -0.000820761, 0.0107715, 0.000498551, -0.0103222, 4.41678e-05, -0.0107814, -2.53282e-05, -0.0127258, -0.000940058, 0.0289788, 0.0014624, -0.0220611, -0.00232308, 0.020889, -0.00119751, -0.0195224, -0.0054188, 0.00288657, 0.00216696, 0.620341, -0.00132468, 0.00705265, -0.0283331, -0.530809, 0.0284523, 0.0113299 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.0239176, 0.000857987, -0.00450675, 0.000856695, 0.02874, 0.000187427, -0.00455364, 0.000190724, -0.0493481, 0.0239654, 0.000857013, -0.00385733, 0.00085832, 0.0287383, 0.000149621, -0.0045049, 0.00017326, -0.049394, 1.04609e-07, -3.97044e-08, 5.80952e-06, -1.12227e-09, -3.08793e-08, -2.63169e-07, 8.76244e-07, 1.59099e-07, 5.33961e-05, 2.59472e-09, 4.1081e-08, -2.47874e-07, -1.40545e-09, 4.52291e-08, 1.47874e-08, 1.90613e-08, 1.74068e-06, -1.72341e-06, 8.30942e-07, 9.12709e-08, 4.64411e-05, 2.13365e-08, 1.57542e-06, -1.81687e-06, -1.02906e-07, -5.43073e-09, -5.819e-06, 0, 5.1556, 1.01106e-24, 0, 0, 0, 0.103075, 0, 1.68732e-23, 2.0688e-23 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }

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

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_10, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.00;//0.96875;

    std::vector< double > _time = { 5.400000e-01, 1.000000e-02 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = {  2.000000e+00, 1.000000e+04, 1.000000e-08, 2.000000e+00, 3.192203e+00, -1.900707e+01, 2.000000e+00, 1.000000e+04, 1.000000e-08, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 6.964416e+02, 1.267138e+02, 5.000000e+00, -1.867498e+01, -3.781732e+01, 1.517765e+01, -2.407120e+01, -5.861821e+00, 1.100000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 7.925235e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, -3.781732e+01, -5.861821e+00, 5.000000e-01, 5.000000e-01, 5.000000e-01, 1.000000e-09, 1.000000e-09 };

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] = { {  2.016181e-02,  3.408375e-05,  1.019530e-03 },
                                         {  3.408287e-05,  2.016187e-02,  1.019557e-03 },
                                         { -4.766168e-05, -4.766289e-05, -4.363897e-02 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  2.016181e-02,  3.408375e-05,  1.019530e-03 },
                                         {  3.408287e-05,  2.016187e-02,  1.019557e-03 },
                                         { -4.766168e-05, -4.766289e-05, -4.363897e-02 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 2.818588e-02, -3.770135e-08, 1.317916e-05, -3.773425e-08, 2.818592e-02, 1.317944e-05, -1.009997e-05, -1.010019e-05, -6.148240e-02 };

    double previous_phi[ 9 ] = { 2.818588e-02, -3.770135e-08, 1.317916e-05, -3.773425e-08, 2.818592e-02, 1.317944e-05, -1.009997e-05, -1.010019e-05, -6.148240e-02 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  6.504275e-06,  8.148054e-06,  3.063838e-05 },
                                           { -1.843628e-07, -1.843628e-07, -5.058133e-08 },
                                           {  6.338505e-05, -2.557240e-07, -6.822002e-05 },
                                           { -1.845240e-07, -1.845240e-07, -5.060166e-08 },
                                           {  8.147958e-06,  6.504348e-06,  3.063842e-05 },
                                           { -2.557651e-07,  6.338521e-05, -6.822144e-05 },
                                           { -4.857690e-05,  1.223168e-07,  5.228101e-05 },
                                           {  1.223387e-07, -4.857707e-05,  5.228216e-05 },
                                           {  7.998816e-06,  7.998866e-06,  3.576544e-06 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {   6.504275e-06,  8.148054e-06,  3.063838e-05 },
                                           {  -1.843628e-07, -1.843628e-07, -5.058133e-08 },
                                           {   6.338505e-05, -2.557240e-07, -6.822002e-05 },
                                           {  -1.845240e-07, -1.845240e-07, -5.060166e-08 },
                                           {   8.147958e-06,  6.504348e-06,  3.063842e-05 },
                                           {  -2.557651e-07,  6.338521e-05, -6.822144e-05 },
                                           {  -4.857690e-05,  1.223168e-07,  5.228101e-05 },
                                           {   1.223387e-07, -4.857707e-05,  5.228216e-05 },
                                           {   7.998816e-06,  7.998866e-06,  3.576544e-06 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = { 1.469221e-02, 2.086219e-05, 3.224448e-04, 2.086219e-05, 1.469225e-02, 3.224535e-04, 3.260838e-04, 3.260926e-04, -2.871940e-02, 1.469245e-02, 2.109984e-05, 2.839340e-04, 2.109983e-05, 1.469249e-02, 2.839416e-04, 3.367983e-04, 3.368074e-04, -2.871987e-02, -2.760823e-08, 9.251805e-11, 3.110196e-08, -1.624211e-08, -1.127405e-08, 3.110273e-08, 2.000602e-06, -7.160067e-09, -2.261297e-06, -1.127431e-08, -1.624173e-08, 3.110267e-08, 9.254277e-11, -2.760903e-08, 3.110345e-08, -7.161508e-09, 2.000608e-06, -2.261346e-06, 1.645207e-06, -3.302415e-09, -1.846872e-06, -3.303227e-09, 1.645214e-06, -1.846914e-06, 2.745892e-08, 2.745975e-08, -6.207796e-08, 3.706097e-24, 1.742749e-01, 3.377384e-26, -1.281992e-25, 1.914138e-25, -1.612093e-22, 5.710195e-02, -8.967331e-24, -3.701739e-23, 9.688256e-24 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { -1.52862, 0.00313534, 0.0277891, 0.00313538, -1.52861, 0.0277898, 0.0687268, 0.0687286, -6.32312 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -1.58169, 0.00299224, 0.046074, 0.00299224, -1.58169, 0.0460752, 0.046074, 0.0460752, -6.32598 };

    tardigradeMicromorphicTools::floatVector M_answer = { 0.00481945, -0.00015429, 0.0496354, -0.000102489, 0.00602015, -0.000175201, -0.0380521, 6.8747e-05, 0.00663771, 0.00602022, -0.000102371, -0.000175168, -0.000154409, 0.00481951, 0.0496355, 6.87297e-05, -0.0380522, 0.00663775, 0.0247419, -6.30897e-05, -0.0582246, -6.31063e-05, 0.0247419, -0.0582258, 0.0447542, 0.0447552, 0.00316672 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.0146922, 2.08622e-05, 0.000322445, 2.08622e-05, 0.0146922, 0.000322454, 0.000326084, 0.000326093, -0.0287194, 0.0146925, 2.10998e-05, 0.000283934, 2.10998e-05, 0.0146925, 0.000283942, 0.000336798, 0.000336807, -0.0287199, -2.76082e-08, 9.25181e-11, 3.1102e-08, -1.62421e-08, -1.12741e-08, 3.11027e-08, 2.0006e-06, -7.16007e-09, -2.2613e-06, -1.12743e-08, -1.62417e-08, 3.11027e-08, 9.25428e-11, -2.7609e-08, 3.11034e-08, -7.16151e-09, 2.00061e-06, -2.26135e-06, 1.64521e-06, -3.30242e-09, -1.84687e-06, -3.30323e-09, 1.64521e-06, -1.84691e-06, 2.74589e-08, 2.74597e-08, -6.2078e-08, 1.50748e-21, 1.40183e-21, 0, 0, 1.38917e-20, 0, 0.0571019, 0, 1.78466e-23, 0 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }

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

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_11, * boost::unit_test::tolerance( 5e-4 ) ){
    /*!
     * Test the evaluation of the constitutive model.
     *
     */

    //Initialize the time

    double s = 1.00;//0.96875;

    std::vector< double > _time = { 6.600000e-01, 5.000000e-03 };

    std::vector< double > time = { _time[ 0 ] - _time[ 1 ] * ( 1 - s ), s * _time[ 1 ] };

    //Initialize the material parameters
    std::vector< double > fparams = { 2.000000e+00, 1.000000e+04, 1.000000e-08, 2.000000e+00, 3.192203e+00, -1.900707e+01, 2.000000e+00, 1.000000e+04, 1.000000e-08, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, 6.964416e+02, 1.267138e+02, 5.000000e+00, -1.867498e+01, -3.781732e+01, 1.517765e+01, -2.407120e+01, -5.861821e+00, 1.100000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 7.925235e+02, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 2.000000e+00, -3.781732e+01, -5.861821e+00, 5.000000e-01, 5.000000e-01, 5.000000e-01, 1.000000e-09, 1.000000e-09 };

    //Initialize the gradient of the macro displacement
    double _current_grad_u[ 3 ][ 3 ] = { {  2.386834e-02,  2.676586e-03, -8.002149e-03 },
                                         {  2.325735e-03,  2.312483e-02,  9.585999e-03 },
                                         {  1.068024e-03, -1.331543e-03, -4.598874e-02 } };

    double previous_grad_u[ 3 ][ 3 ] = { {  2.386834e-02,  2.676586e-03, -8.002149e-03 },
                                         {  2.325735e-03,  2.312483e-02,  9.585999e-03 },
                                         {  1.068024e-03, -1.331543e-03, -4.598874e-02 } };

    double current_grad_u[ 3 ][ 3 ];
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_u[ i ][ j ] = ( 1 - s ) * previous_grad_u[ i ][ j ] + s * _current_grad_u[ i ][ j ];
        }
    }

    //Initialize the micro displacement
    double _current_phi[ 9 ] = { 2.918502e-02, 8.031427e-06, -1.089626e-04, 6.145852e-06, 2.918315e-02, 1.277158e-04, 8.762238e-05, -1.025436e-04, -6.228604e-02 };

    double previous_phi[ 9 ] = { 2.918502e-02, 8.031427e-06, -1.089626e-04, 6.145852e-06, 2.918315e-02, 1.277158e-04, 8.762238e-05, -1.025436e-04, -6.228604e-02 };

    double current_phi[ 9 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        current_phi[ i ] = ( 1 - s ) * previous_phi[ i ] + s * _current_phi[ i ];
    }

    //Initialize the gradient of the micro displacement
    double _current_grad_phi[ 9 ][ 3 ] = { {  4.519668e-06, -1.772929e-05, -3.536672e-05 },
                                           { -2.169231e-07, -3.225982e-06, -4.003954e-06 },
                                           { -4.507592e-05, -2.217632e-05, -5.640293e-04 },
                                           {  6.109742e-06,  2.429224e-06, -3.329048e-06 },
                                           {  1.070898e-05, -1.192545e-05, -3.466505e-05 },
                                           { -1.763535e-05, -3.421614e-05,  6.611021e-04 },
                                           {  3.636262e-05,  1.737356e-05,  4.535645e-04 },
                                           {  1.351496e-05,  2.752710e-05, -5.308021e-04 },
                                           {  8.501887e-06, -1.669006e-05, -1.248842e-05 } };

    double previous_grad_phi[ 9 ][ 3 ] = { {  4.519668e-06, -1.772929e-05, -3.536672e-05 },
                                           { -2.169231e-07, -3.225982e-06, -4.003954e-06 },
                                           { -4.507592e-05, -2.217632e-05, -5.640293e-04 },
                                           {  6.109742e-06,  2.429224e-06, -3.329048e-06 },
                                           {  1.070898e-05, -1.192545e-05, -3.466505e-05 },
                                           { -1.763535e-05, -3.421614e-05,  6.611021e-04 },
                                           {  3.636262e-05,  1.737356e-05,  4.535645e-04 },
                                           {  1.351496e-05,  2.752710e-05, -5.308021e-04 },
                                           {  8.501887e-06, -1.669006e-05, -1.248842e-05 } };

    double current_grad_phi[ 9 ][ 3 ];
    for ( unsigned int i = 0; i < 9; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            current_grad_phi[ i ][ j ] = ( 1 - s ) * previous_grad_phi[ i ][ j ] + s * _current_grad_phi[ i ][ j ];
        }
    }

    //Initialize the state variable vector
    std::vector< double > SDVSDefault = {  1.804953e-02, 1.923236e-03, -2.580802e-03, 1.923330e-03, 1.739654e-02, 3.054039e-03, -2.613184e-03, 3.093371e-03, -3.447637e-02, 1.806610e-02, 1.904158e-03, -2.245439e-03, 1.904055e-03, 1.741872e-02, 2.655300e-03, -2.673778e-03, 3.160835e-03, -3.451498e-02, -1.541664e-07, -6.114041e-08, -2.072806e-06, 9.728821e-08, -2.531493e-11, 2.443971e-06, -1.779007e-06, -8.023710e-07, -2.145868e-05, 4.015326e-08, -4.756031e-08, 2.443417e-06, 5.405566e-08, 1.426229e-07, -2.880821e-06, -6.392950e-07, -1.384139e-06, 2.514540e-05, -1.521880e-06, -6.684446e-07, -1.854388e-05, -5.186341e-07, -1.179013e-06, 2.168409e-05, 9.985578e-08, -8.124686e-08, 4.940638e-06, 4.205015e-21, 1.627078e-01, -2.428035e-22, -3.981793e-25, 3.541683e-24, 1.101612e-21, 6.977922e-02, 1.369254e-21, -3.057028e-22, -2.150737e-21 };

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

    tardigradeMicromorphicTools::floatVector PK2_answer = { 1.02934, 0.134176, -0.100326, 0.147396, 1.00518, 0.121048, -0.451746, 0.54424, -2.85345 };

    tardigradeMicromorphicTools::floatVector SIGMA_answer = { 0.909227, 0.129612, -0.262828, 0.129612, 0.887974, 0.31678, -0.262828, 0.31678, -2.9394 };

    tardigradeMicromorphicTools::floatVector M_answer = { 0.00332625, -1.36795e-05, -0.0373681, 0.00409052, 0.00767053, -0.0108531, 0.030506, 0.00807972, 0.00721803, -0.0127424, -0.00203612, -0.0144148, 0.00164591, -0.00871738, -0.0297955, 0.0110717, 0.0244136, -0.0140619, -0.0267161, -0.00520809, -0.489739, -0.00475141, -0.0254647, 0.57392, 0.397736, -0.465501, -0.00717349 };

    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.0180495, 0.00192324, -0.0025808, 0.00192333, 0.0173965, 0.00305404, -0.00261318, 0.00309337, -0.0344764, 0.0180661, 0.00190416, -0.00224544, 0.00190406, 0.0174187, 0.0026553, -0.00267378, 0.00316084, -0.034515, -1.54166e-07, -6.11404e-08, -2.07281e-06, 9.72882e-08, -2.53153e-11, 2.44397e-06, -1.77901e-06, -8.02371e-07, -2.14587e-05, 4.01533e-08, -4.75603e-08, 2.44342e-06, 5.40557e-08, 1.42623e-07, -2.88082e-06, -6.39295e-07, -1.38414e-06, 2.51454e-05, -1.52188e-06, -6.68445e-07, -1.85439e-05, -5.18634e-07, -1.17901e-06, 2.16841e-05, 9.98558e-08, -8.12469e-08, 4.94064e-06, 4.32112e-26, 3.74833e-07, 2.51819e-25, 0, 0, 1.1015e-21, 0.0697792, 1.36829e-21, 0, 0 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 55, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }

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
////    tardigradeMicromorphicTools::floatVector PK2_answer = { -3.37057, 0.0866543, 0.0323353, 0.0818984, -3.36447, -0.00167251, 0.0276005, -0.0735036, -3.36448 };
////
////    tardigradeMicromorphicTools::floatVector SIGMA_answer = { -3.34476, 0.066115, 0.0253977, 0.066115, -3.33988, -0.031869, 0.0253977, -0.031869, -3.47732 };
////
////    tardigradeMicromorphicTools::floatVector M_answer = { 0.00222842, -0.0157762, 0.00575456, 0.00895602, -0.00581722, 0.00508999, -0.00693765, -0.00227177, 0.0143212, 0.00293983, 0.00706979, 0.00630251, 0.00209025, -0.00530594, 0.00931891, -0.00343587, -0.0104565, -0.0224589, 0.00330374, 0.00677226, 0.488466, -0.00499722, 0.00440454, -0.3879, -0.441158, 0.356116, 0.000501602 };
////
////    tardigradeMicromorphicTools::floatVector SDVS_answer = { 0.0107416, 0.00643479, 0.00142293, 0.00643734, 0.0112673, -0.00235179, 0.00153702, -0.00250549, -0.0212581, 0.0107496, 0.00642886, 0.00149239, 0.00642612, 0.0112775, -0.00236504, 0.00165881, -0.00264259, -0.0212743, 1.42267e-07, -8.46723e-08, -2.26279e-06, -4.38303e-08, -8.23814e-08, 2.86685e-06, 7.82471e-07, -1.08144e-07, 1.15626e-05, 7.17795e-08, 3.79839e-08, 2.75112e-06, -1.58217e-07, 1.34147e-07, -2.35395e-06, -1.90528e-07, 7.9154e-07, -5.86987e-06, 6.55041e-07, 1.58951e-07, 9.07531e-06, 3.67226e-08, 6.23863e-07, -4.11633e-06, 2.16544e-08, -5.33505e-08, 4.554e-06, -1.7134e-23, 8.8538, -8.51183e-24, -1.71876e-23, -5.08624e-24, -3.9568e-24, 0.109811, 1.84609e-23, 7.51148e-24, -5.81825e-25 };
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
////    tardigradeMicromorphicTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;
////
////    tardigradeMicromorphicTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;
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
////            void public_setUnknownVector( const tardigradeMicromorphicTools::floatVector &value ){ updateUnknownVector( value ); }
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
