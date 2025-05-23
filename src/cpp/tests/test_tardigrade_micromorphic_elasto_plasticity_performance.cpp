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

    tardigradeConstitutiveTools::floatVector PK2_answer = { -1.40279, 0.0121485, 0.0266627, 0.00381226, -1.45011, 0.0418199, 0.0338264, 0.0864571, -2.91531 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { -0.493708, 0.0121523, 0.0789793, 0.0121523, -0.517776, 0.106717, 0.0789793, 0.106717, -4.81709 };

    tardigradeConstitutiveTools::floatVector M_answer = { 8.45418e-05, -4.04181e-06, -0.000551314, -2.62816e-06, 9.35849e-05, 4.20583e-05, 0.000290068, -2.54552e-05, -6.72292e-05, -1.81165e-05, 6.55509e-07, 0.000106142, -1.94223e-06, -1.8193e-05, 8.36339e-06, -5.40107e-05, -4.57942e-06, 2.03912e-06, -0.000144554, 2.02458e-05, 0.00030473, 1.70095e-05, -0.000104047, 0.00037256, -0.000143114, -0.000161128, 7.22794e-05 };

    tardigradeConstitutiveTools::floatVector SDVS_answer = { 0.00458246, 3.26705e-05, 0.000297186, 0.000112011, 0.00414933, 0.000768014, 0.000230621, 0.000356534, -0.00861812, 1.43254e-11, 1.68409e-13, 8.44243e-13, 2.5468e-13, 1.31286e-11, 1.89224e-12, 7.70904e-13, 1.41952e-12, -2.64072e-11, 2.79233e-17, 5.56405e-20, -7.87013e-18, -3.37087e-18, 1.51828e-17, -1.01737e-17, 2.41855e-17, -1.34546e-18, -5.24753e-17, -2.32027e-18, 1.54445e-17, 6.40194e-18, -3.58006e-19, -4.96139e-18, -4.77126e-18, -3.94101e-18, -1.629e-17, -4.58044e-17, -1.27547e-16, 1.15812e-17, 6.96373e-17, 2.14292e-17, -2.77694e-17, 1.05305e-16, 6.73349e-18, 1.54864e-17, 3.42911e-17, 0.170641, 4.09748e-26, 2.09674e-24, 0, 0, 0.0172535, 0, 5.60883e-25, 0, 1.61587e-23 };

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

    tardigradeConstitutiveTools::floatVector PK2_answer = {
        7.20108271, -0.01041904,  0.40177972, -0.0312383 ,  7.31864868,
       -0.50171677, -0.44379043,  0.64151818,  6.95971366
    };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = {
        7.17850397, -0.05454965, -0.09879296, -0.05454965,  7.33567775,
        0.10928604, -0.09879296,  0.10928604,  6.81403522
    };

    tardigradeConstitutiveTools::floatVector M_answer = {
       -0.02630482,  0.00495418, -0.30777728, -0.01369079,  0.00198282,
        0.18441233,  0.6257582 , -0.13247079,  0.05863422, -0.04336741,
        0.01239835, -0.05725744, -0.02028124, -0.03315325,  0.00341915,
        0.18114804,  0.01745458, -0.12469868,  0.26992647, -0.1081408 ,
       -0.07141717, -0.01523851,  0.32051147,  0.55389985,  0.26501641,
       -0.82205036, -0.04292842
    };

    tardigradeConstitutiveTools::floatVector SDVS_answer = {
       -8.84270338e-03,  7.20015515e-03,  1.31995560e-02,  7.22663814e-03,
       -2.83276134e-02, -1.33543175e-02,  1.32149447e-02, -1.34292212e-02,
        4.04521652e-02, -9.13625734e-03,  7.12621145e-03,  1.27927132e-02,
        7.96593061e-03, -2.86990824e-02, -1.16353910e-02,  1.23124216e-02,
       -1.29924944e-02,  4.11171881e-02,  1.53774583e-05,  3.57374000e-06,
        7.61976200e-06, -7.29781900e-06,  1.20094000e-07, -1.94019500e-05,
        2.17407640e-05,  2.55382000e-06,  2.37251600e-06, -1.25548820e-05,
       -2.12581000e-06, -1.41241900e-05,  5.43542100e-06, -4.81407200e-07,
        2.23475100e-05, -1.98101860e-05,  1.94556800e-06, -4.69112900e-05,
        4.33254350e-05,  1.28145420e-05,  3.56474400e-05, -1.89541300e-05,
        8.52890000e-07, -9.48495700e-05, -2.08128790e-05, -3.09233300e-06,
       -2.99672700e-05,  2.22256700e-24,  1.11040652e+00,  3.24545021e-24,
       -5.00929980e-24,  1.85134758e-24,  5.07725800e-25,  9.35619438e-02,
        5.13671200e-24,  1.96936940e-23, -4.48641106e-24
    };

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
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradU     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 3 ] ), 1e-5, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradU ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 6 ] ), 1e-5, 1e-5 ) );

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
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdPhi )    , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 4 ] ), 1e-5, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 7 ] ), 1e-5, 1e-5 ) );

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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   1e-3, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 2e-3, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     1e-3, 1e-5 ) );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 1e-4, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 1e-4, 1e-5 ) );
    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 1e-4, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( testEvaluateHydraModel_difficult_3, * boost::unit_test::tolerance( 1e-5 ) ){ //TODO: Maybe there is an error in the gradients w.r.t. the displacement?
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

    tardigradeConstitutiveTools::floatVector PK2_answer = { -3.34995, 0.0663803, 0.0199674, 0.0616158, -3.34511, 0.0108698, 0.0150697, -0.0607778, -3.4235 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { -3.32432, 0.0463883, 0.0132939, 0.0463883, -3.32068, -0.0195832, 0.0132939, -0.0195832, -3.53427 };

    tardigradeConstitutiveTools::floatVector M_answer = {
        2.24148310e-03, -1.57996933e-02,  5.69204628e-03,  8.94765320e-03,
       -5.80709750e-03,  5.12849644e-03, -6.88654108e-03, -2.31106138e-03,
        1.42923482e-02,  2.92427180e-03,  7.09867640e-03,  6.35672884e-03,
        2.10304130e-03, -5.32400020e-03,  9.29336870e-03, -3.48739774e-03,
       -1.03941608e-02, -2.24154028e-02,  3.30376960e-03,  6.77260140e-03,
        4.87936082e-01, -5.02161610e-03,  4.41448610e-03, -3.87543447e-01,
       -4.41102527e-01,  3.56060741e-01,  4.69478500e-04
    };

    tardigradeConstitutiveTools::floatVector SDVS_answer = {
        1.06225659e-02,  6.52946087e-03,  1.47773829e-03,  6.53204655e-03,
        1.11543282e-02, -2.40771513e-03,  1.59416072e-03, -2.56404505e-03,
       -2.10244612e-02,  1.06305960e-02,  6.52351526e-03,  1.54492499e-03,
        6.52073143e-03,  1.11647019e-02, -2.41911503e-03,  1.71859021e-03,
       -2.70384560e-03, -2.10407978e-02,  1.28607800e-07, -6.49381000e-08,
       -2.26648120e-06, -2.89144000e-08, -1.05412500e-07,  2.88411470e-06,
        7.87694150e-07, -1.10177520e-07,  1.18550570e-05,  8.98311000e-08,
        1.93577000e-08,  2.76730870e-06, -1.73430300e-07,  1.52410200e-07,
       -2.36834840e-06, -1.90684000e-07,  7.64207400e-07, -6.05801300e-06,
        6.52302608e-07,  1.68751250e-07,  9.31367300e-06,  4.46800800e-08,
        5.90448100e-07, -4.26414300e-06,  5.06897000e-08, -9.11532000e-08,
        4.58120610e-06,  2.46040200e-27,  8.98758235e+00,  1.46244000e-27,
        2.79192900e-28,  8.03125000e-28, -2.19390000e-28,  1.10356651e-01,
        7.03070000e-27, -3.82897000e-27, -5.28030000e-28
    };

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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdPhi )      , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 1 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     1e-5, 1e-5 ) );
//
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdGradPhi       ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 2 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dChipdGradPhi     ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 5 ] ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChipdGradPhi ), tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 8 ] ), 1e-5, 1e-5 ) );

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

    tardigradeConstitutiveTools::floatVector PK2_answer = { 0.417473, 0.952059, -0.859843, 0.876019, -0.734345, -0.0504722, -0.429526, -0.516148, -0.32398 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { 0.312962, 0.821419, -0.62471, 0.821419, -0.7273, -0.264341, -0.62471, -0.264341, -0.507436 };

    tardigradeConstitutiveTools::floatVector M_answer = {
       -0.01008223,  0.00333048,  0.43691041, -0.07350997,  0.03750151,
       -0.12501433, -0.32751909,  0.15661292,  0.05710477, -0.03517283,
       -0.0639151 , -0.02273301,  0.26921522,  0.05161277,  0.32968169,
        0.06954138, -0.36992541,  0.03800122,  0.04644254,  0.05712133,
        1.25486296, -0.03663805,  0.03635676,  0.03749855, -1.20025915,
       -0.11213527,  0.02765728
    };

    tardigradeConstitutiveTools::floatVector SDVS_answer = {
        2.33657559e-02,  2.41927512e-02, -1.88951748e-02,  2.42304618e-02,
       -7.56719281e-03, -7.61511508e-03, -1.91573097e-02, -7.68842791e-03,
       -1.31780710e-02,  2.31106316e-02,  2.40005089e-02, -1.81468785e-02,
        2.46075355e-02, -7.42166948e-03, -7.84532280e-03, -1.98988044e-02,
       -7.86321690e-03, -1.30702444e-02,  1.54103900e-05,  7.88425700e-06,
        5.56283700e-05,  2.57696500e-06,  8.39921800e-06,  1.81916640e-05,
        1.09872010e-05,  6.63132600e-06,  4.20932600e-05,  1.22185200e-06,
       -5.76551160e-06,  1.55282650e-05, -1.07303500e-06, -2.49779020e-06,
        4.37094460e-06,  1.17253240e-05,  6.55734280e-06,  3.45585800e-05,
        9.63706500e-06,  8.13320500e-06,  4.35385400e-05,  1.06452520e-05,
       -1.23097900e-06,  3.35317100e-05, -1.42624800e-05, -5.37251300e-06,
       -5.93247200e-05, -1.96818000e-24,  4.97209221e+00, -9.74748000e-24,
       -5.35841000e-25,  5.70249000e-24,  3.44233000e-23,  9.60879258e-02,
       -4.84558000e-24, -1.32241600e-23,  1.81604000e-24
    };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS = SDVSDefault;

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFpdPhi )      , tardigradeVectorTools::appendVectors( ADD_JACOBIANS[ 1 ] ), 1e-5, 1e-5 ) );
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
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPK2dGradPhi   ), tardigradeVectorTools::appendVectors( result_dPK2dGradPhi ),   1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dSIGMAdGradPhi ), tardigradeVectorTools::appendVectors( result_dSIGMAdGradPhi ), 1e-5, 1e-5 ) );
//    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMdGradPhi     ), tardigradeVectorTools::appendVectors( result_dMdGradPhi ),     1e-5, 1e-5 ) );
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

    tardigradeConstitutiveTools::floatVector PK2_answer = { -9.62479, 0.38928, 0.361314, 0.513209, -9.65339, -0.531695, -0.0320208, 0.105003, -8.98122 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { -9.34451, 0.375142, 0.127635, 0.375142, -9.37008, -0.165708, 0.127635, -0.165708, -8.99023 };

    tardigradeConstitutiveTools::floatVector M_answer = {
        0.0135374 , -0.05448412, -0.08729431,  0.03888332,  0.00699627,
       -0.0119041 ,  0.07831853,  0.01190304,  0.04088716, -0.00823423,
       -0.01626964, -0.00701686,  0.03308791, -0.01670416, -0.0739713 ,
        0.0096626 ,  0.06296182, -0.05103458, -0.00208876,  0.003249  ,
        0.4228994 , -0.0068872 , -0.00254739, -0.40908723, -0.40104114,
        0.39179197, -0.01987857
    };

    tardigradeConstitutiveTools::floatVector SDVS_answer = {
       -5.40907329e-03,  2.01288291e-02,  6.44232033e-03,  2.01197997e-02,
       -6.76143024e-03, -8.49110082e-03,  6.90440383e-03, -9.08190813e-03,
        1.57182672e-02, -5.40473375e-03,  2.02056941e-02,  5.63264890e-03,
        2.02134438e-02, -6.95702466e-03, -7.41806276e-03,  6.34140280e-03,
       -8.37405885e-03,  1.58932011e-02,  3.40998000e-06,  1.57666200e-06,
       -7.03941300e-06, -9.99940000e-07, -7.69582000e-08,  7.01062900e-06,
        5.10383700e-06, -7.11000000e-08, -2.24541300e-05, -4.49461000e-07,
        4.28560000e-07,  6.90921300e-06, -2.74621200e-06, -2.52397300e-06,
       -7.78800200e-06, -2.39100000e-07,  3.28287000e-06,  2.63139000e-05,
        5.12027500e-06,  6.70066000e-07, -3.05212700e-05,  2.63208000e-07,
        3.49620600e-06,  3.38502200e-05, -6.38227000e-07,  9.26667000e-07,
        1.46590850e-05,  8.12132000e-24,  1.16918732e+02,  4.81497100e-23,
        3.58328400e-22,  8.88797000e-24, -1.80183000e-24,  1.98079681e-01,
        1.20606160e-21, -5.12651000e-23,  1.69688840e-21
    };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }

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

    tardigradeConstitutiveTools::floatVector PK2_answer = { -3.346, 0.131246, 0.280253, 0.0873976, -3.5572, 0.0595398, -0.256392, -0.286154, -6.42817 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { -3.33216, 0.103392, 0.010602, 0.103392, -3.53612, -0.110036, 0.010602, -0.110036, -6.38187 };

    tardigradeConstitutiveTools::floatVector M_answer = {
       -0.00800279,  0.0058371 , -0.02038557, -0.00644164, -0.01318874,
        0.00775058,  0.01399107, -0.00058195, -0.0099667 ,  0.00665797,
        0.00321174, -0.01511046,  0.00145756,  0.00678377, -0.03141141,
        0.00888939,  0.0425507 ,  0.00949888,  0.0009563 , -0.02658308,
        0.4554901 ,  0.02723584,  0.00082363, -0.27391843, -0.37497247,
        0.22133887, -0.00758202
    };

    tardigradeConstitutiveTools::floatVector SDVS_answer = {
        1.63015810e-02,  1.36103037e-03, -7.12160984e-04,  1.35993633e-03,
        1.62714178e-02, -5.99446620e-05, -7.09996794e-04, -7.84993852e-05,
       -3.17441716e-02,  1.63033437e-02,  1.36228278e-03, -5.25351224e-04,
        1.36339063e-03,  1.62739276e-02, -2.93465286e-05, -6.30232867e-04,
       -1.14367574e-05, -3.17484381e-02,  1.29115000e-12,  9.87671000e-10,
        9.86739096e-07, -1.73345200e-09,  8.30069000e-10,  2.71748500e-08,
       -3.16801080e-08, -2.54393300e-08,  2.16329415e-05,  3.12790000e-10,
       -1.05146000e-10,  2.27506250e-08,  1.37312030e-09, -4.59623180e-09,
       -3.12488800e-08,  1.01255830e-08, -4.58969107e-08, -1.33322646e-05,
       -2.39237090e-08, -1.83361000e-08,  1.95448928e-05,  1.25478900e-09,
       -6.66049760e-08, -1.22801266e-05, -1.28705030e-09,  3.95768920e-09,
        4.40925800e-08,  2.11770622e-20,  1.96590697e+02,  3.48453949e-22,
        1.99449390e-21, -4.43231000e-23, -9.53437000e-24,  6.36231563e-02,
       -1.47976000e-24, -2.11506000e-23, -5.41101000e-23
    };

    cleanAnswer( SDVS_answer, 1e-8 );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }

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

    tardigradeConstitutiveTools::floatVector PK2_answer = { 0.188943, -0.0749734, -0.223576, -0.0766538, 0.350907, 0.308501, -0.349782, 0.468328, -5.0104 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { 0.145838, -0.0730737, -0.275922, -0.0730737, 0.30145, 0.373879, -0.275922, 0.373879, -5.00914 };

    tardigradeConstitutiveTools::floatVector M_answer = { -0.00304811, -0.000291436, 0.223069, -0.000555977, -0.00465674, 0.00388556, -0.169641, -0.00228998, -0.00469748, 0.005428, -0.00175762, 0.00470323, -0.000534611, 0.00769539, 0.222371, -0.00383504, -0.170343, 0.00729416, 0.0152528, 0.00078359, 0.0796262, -3.46937e-05, 0.0152673, -0.100393, -0.0582237, 0.0741628, -0.0211573 };

    tardigradeConstitutiveTools::floatVector SDVS_answer = { 0.004985, -7.5e-05, -0.000566, -7.5e-05, 0.005112, 0.000686, -0.00057, 0.000691, -0.010012, 0.004987, -7.7e-05, -0.000502, -7.7e-05, 0.005115, 0.000606, -0.00059, 0.000712, -0.010017, 0, 0, 0, -0, 0, -0, 4e-06, 0, 1e-06, -0, 0, -0, -0, -0, 0, 0, 4e-06, -2e-06, 3e-06, 0, 1e-06, 0, 3e-06, -1e-06, -0, 0, -0, 0, 0, -0, -0, 0, -0, 0.019824, -0, 0, -0, };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }

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

    tardigradeConstitutiveTools::floatVector PK2_answer = { 1.49261, -1.04731, 1.8083, 0.736858, -0.342405, 0.384324, 0.144605, 0.532215, 5.33617 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { 1.52344, -0.159191, 0.833773, -0.159191, -0.148749, 0.402651, 0.833773, 0.402651, 4.95921 };

    tardigradeConstitutiveTools::floatVector M_answer = {
       -0.02807644,  0.0577785 ,  0.57195511, -0.03129482, -0.01628047,
        0.03076009, -0.54634468,  0.02474005, -0.04155135, -0.02317139,
        0.03517286,  0.01550789, -0.04319591, -0.03325799,  0.16806505,
       -0.00528952, -0.14044434, -0.03961   ,  0.105021  , -0.06275794,
       -0.21442937,  0.0486657 ,  0.10081636, -0.02139026,  0.20188963,
        0.04593065,  0.09269014
    };

    tardigradeConstitutiveTools::floatVector SDVS_answer = {
       -1.34854945e-02, -2.31109245e-03,  1.81625022e-02, -2.31944462e-03,
       -4.20582772e-02,  7.99271727e-03,  1.83601866e-02,  8.08635744e-03,
        6.19512871e-02, -1.43809004e-02, -2.98612595e-03,  1.53428794e-02,
       -3.01210996e-03, -4.19319768e-02,  7.18092889e-03,  1.81363462e-02,
        8.50674572e-03,  6.27203926e-02, -2.40823600e-05, -2.34540000e-08,
        9.41980100e-06, -3.45377200e-06, -1.91201400e-06,  1.29483789e-06,
       -5.42591300e-05, -2.22495400e-06,  2.21479810e-05, -4.67801600e-06,
       -2.26616500e-06,  8.32838000e-07, -1.51527700e-07, -3.08399600e-06,
        1.02318140e-06, -3.32987300e-06, -1.88021200e-05,  1.54856100e-06,
       -5.93264500e-05, -4.29443000e-07,  2.46497860e-05,  3.56724900e-06,
       -1.81307700e-05,  6.50180200e-06,  2.42338900e-05,  3.10745000e-06,
       -1.04429820e-05,  2.09089000e-27,  7.98743560e+00, -8.71034000e-27,
        5.76008000e-27,  2.08167000e-27,  3.25577900e-26,  1.30434277e-01,
       -6.03618000e-28, -4.66345000e-27,  5.27667000e-27
    };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }

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

    tardigradeConstitutiveTools::floatVector PK2_answer = { -3.5424, -0.036538, 0.175678, 0.052668, -3.49538, -0.0209274, -0.267199, 0.0249727, -4.11429 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { -3.56824, 0.00601439, -0.0358833, 0.00601439, -3.53457, 0.0015479, -0.0358833, 0.0015479, -4.16335 };

    tardigradeConstitutiveTools::floatVector M_answer = {
       -8.59579680e-03, -8.22779960e-04,  1.11066630e-02,  4.96429340e-04,
       -1.03220469e-02,  2.92879000e-05, -1.04886860e-02, -4.11130000e-05,
       -1.27640607e-02, -9.40792549e-04,  2.89790296e-02,  1.43909310e-03,
       -2.20623745e-02, -2.32463735e-03,  2.08256870e-02, -1.21633910e-03,
       -1.95809085e-02, -5.41662971e-03,  3.36576500e-03,  2.14470950e-03,
        6.24285610e-01, -1.34675770e-03,  7.05388901e-03, -2.84676470e-02,
       -5.27329650e-01,  2.83071850e-02,  1.08769300e-02
    };

    tardigradeConstitutiveTools::floatVector SDVS_answer = {
        2.39176147e-02,  8.57986512e-04, -4.50675201e-03,  8.56694760e-04,
        2.87399521e-02,  1.87427088e-04, -4.55364160e-03,  1.90723618e-04,
       -4.93481386e-02,  2.39653987e-02,  8.57013297e-04, -3.85732803e-03,
        8.58320251e-04,  2.87382696e-02,  1.49621241e-04, -4.50489998e-03,
        1.73260064e-04, -4.93940272e-02,  6.14600000e-08, -3.95604590e-08,
        5.31937200e-06,  1.00579000e-09, -3.22698900e-08, -2.40488700e-07,
        5.07287000e-07,  1.85741900e-07,  4.92379000e-05,  4.85936000e-09,
        4.22898300e-08, -2.25947100e-07, -1.51995300e-09,  4.76057700e-08,
        1.35041500e-08,  3.59858000e-08,  1.81969650e-06, -1.58209600e-06,
        4.81938000e-07,  1.14346400e-07,  4.24655400e-05,  4.03600000e-08,
        1.65455540e-06, -1.65619100e-06, -5.96831000e-08, -7.94436000e-09,
       -5.32824500e-06, -2.41229000e-24,  5.15560462e+00,  3.08875000e-24,
        7.15158000e-25,  6.25792000e-24, -2.43321000e-22,  1.03074621e-01,
       -1.03198000e-22,  3.80914000e-23,  2.44781600e-23
    };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }

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

    tardigradeConstitutiveTools::floatVector PK2_answer = { -1.52862, 0.00313534, 0.0277891, 0.00313538, -1.52861, 0.0277898, 0.0687268, 0.0687286, -6.32312 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { -1.58169, 0.00299224, 0.046074, 0.00299224, -1.58169, 0.0460752, 0.046074, 0.0460752, -6.32598 };

    tardigradeConstitutiveTools::floatVector M_answer = { 0.00481945, -0.00015429, 0.0496354, -0.000102489, 0.00602015, -0.000175201, -0.0380521, 6.8747e-05, 0.00663771, 0.00602022, -0.000102371, -0.000175168, -0.000154409, 0.00481951, 0.0496355, 6.87297e-05, -0.0380522, 0.00663775, 0.0247419, -6.30897e-05, -0.0582246, -6.31063e-05, 0.0247419, -0.0582258, 0.0447542, 0.0447552, 0.00316672 };

    tardigradeConstitutiveTools::floatVector SDVS_answer = { 0.0146922, 2.08622e-05, 0.000322445, 2.08622e-05, 0.0146922, 0.000322454, 0.000326084, 0.000326093, -0.0287194, 0.0146925, 2.10998e-05, 0.000283934, 2.10998e-05, 0.0146925, 0.000283942, 0.000336798, 0.000336807, -0.0287199, -2.76082e-08, 9.25181e-11, 3.1102e-08, -1.62421e-08, -1.12741e-08, 3.11027e-08, 2.0006e-06, -7.16007e-09, -2.2613e-06, -1.12743e-08, -1.62417e-08, 3.11027e-08, 9.25428e-11, -2.7609e-08, 3.11034e-08, -7.16151e-09, 2.00061e-06, -2.26135e-06, 1.64521e-06, -3.30242e-09, -1.84687e-06, -3.30323e-09, 1.64521e-06, -1.84691e-06, 2.74589e-08, 2.74597e-08, -6.2078e-08, 1.50748e-21, 1.40183e-21, 0, 0, 1.38917e-20, 0, 0.0571019, 0, 1.78466e-23, 0 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 60, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }

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

    tardigradeConstitutiveTools::floatVector PK2_answer = { 1.02934, 0.134176, -0.100326, 0.147396, 1.00518, 0.121048, -0.451746, 0.54424, -2.85345 };

    tardigradeConstitutiveTools::floatVector SIGMA_answer = { 0.909227, 0.129612, -0.262828, 0.129612, 0.887974, 0.31678, -0.262828, 0.31678, -2.9394 };

    tardigradeConstitutiveTools::floatVector M_answer = { 0.00332625, -1.36795e-05, -0.0373681, 0.00409052, 0.00767053, -0.0108531, 0.030506, 0.00807972, 0.00721803, -0.0127424, -0.00203612, -0.0144148, 0.00164591, -0.00871738, -0.0297955, 0.0110717, 0.0244136, -0.0140619, -0.0267161, -0.00520809, -0.489739, -0.00475141, -0.0254647, 0.57392, 0.397736, -0.465501, -0.00717349 };

    tardigradeConstitutiveTools::floatVector SDVS_answer = { 0.0180495, 0.00192324, -0.0025808, 0.00192333, 0.0173965, 0.00305404, -0.00261318, 0.00309337, -0.0344764, 0.0180661, 0.00190416, -0.00224544, 0.00190406, 0.0174187, 0.0026553, -0.00267378, 0.00316084, -0.034515, -1.54166e-07, -6.11404e-08, -2.07281e-06, 9.72882e-08, -2.53153e-11, 2.44397e-06, -1.77901e-06, -8.02371e-07, -2.14587e-05, 4.01533e-08, -4.75603e-08, 2.44342e-06, 5.40557e-08, 1.42623e-07, -2.88082e-06, -6.39295e-07, -1.38414e-06, 2.51454e-05, -1.52188e-06, -6.68445e-07, -1.85439e-05, -5.18634e-07, -1.17901e-06, 2.16841e-05, 9.98558e-08, -8.12469e-08, 4.94064e-06, 4.32112e-26, 3.74833e-07, 2.51819e-25, 0, 0, 1.1015e-21, 0.0697792, 1.36829e-21, 0, 0 };

    cleanAnswer( SDVS_answer );

    std::vector< double > SDVS( 55, 0 );
    std::copy( SDVSDefault.begin( ), SDVSDefault.end( ), SDVS.begin( ) );

    // Explore continuation approach

    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;

    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;

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

            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }

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
////    tardigradeConstitutiveTools::floatVector PK2_answer = { -3.37057, 0.0866543, 0.0323353, 0.0818984, -3.36447, -0.00167251, 0.0276005, -0.0735036, -3.36448 };
////
////    tardigradeConstitutiveTools::floatVector SIGMA_answer = { -3.34476, 0.066115, 0.0253977, 0.066115, -3.33988, -0.031869, 0.0253977, -0.031869, -3.47732 };
////
////    tardigradeConstitutiveTools::floatVector M_answer = { 0.00222842, -0.0157762, 0.00575456, 0.00895602, -0.00581722, 0.00508999, -0.00693765, -0.00227177, 0.0143212, 0.00293983, 0.00706979, 0.00630251, 0.00209025, -0.00530594, 0.00931891, -0.00343587, -0.0104565, -0.0224589, 0.00330374, 0.00677226, 0.488466, -0.00499722, 0.00440454, -0.3879, -0.441158, 0.356116, 0.000501602 };
////
////    tardigradeConstitutiveTools::floatVector SDVS_answer = { 0.0107416, 0.00643479, 0.00142293, 0.00643734, 0.0112673, -0.00235179, 0.00153702, -0.00250549, -0.0212581, 0.0107496, 0.00642886, 0.00149239, 0.00642612, 0.0112775, -0.00236504, 0.00165881, -0.00264259, -0.0212743, 1.42267e-07, -8.46723e-08, -2.26279e-06, -4.38303e-08, -8.23814e-08, 2.86685e-06, 7.82471e-07, -1.08144e-07, 1.15626e-05, 7.17795e-08, 3.79839e-08, 2.75112e-06, -1.58217e-07, 1.34147e-07, -2.35395e-06, -1.90528e-07, 7.9154e-07, -5.86987e-06, 6.55041e-07, 1.58951e-07, 9.07531e-06, 3.67226e-08, 6.23863e-07, -4.11633e-06, 2.16544e-08, -5.33505e-08, 4.554e-06, -1.7134e-23, 8.8538, -8.51183e-24, -1.71876e-23, -5.08624e-24, -3.9568e-24, 0.109811, 1.84609e-23, 7.51148e-24, -5.81825e-25 };
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
////    tardigradeConstitutiveTools::floatVector currentDeformationGradient, currentMicroDeformation, currentGradientMicroDeformation;
////
////    tardigradeConstitutiveTools::floatVector previousDeformationGradient, previousMicroDeformation, previousGradientMicroDeformation;
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
////            void public_setUnknownVector( const tardigradeConstitutiveTools::floatVector &value ){ updateUnknownVector( value ); }
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
