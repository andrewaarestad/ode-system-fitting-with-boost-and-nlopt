#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#define WITHOUT_NUMPY

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

double a = 250;
double b = 125;
double c = 0.0001;
double d = 10;
double e = 0.01;
double g = 0;

//struct lorenz
//{
//    template< class State , class Deriv >
//    void operator()( const State &x_ , Deriv &dxdt_ , double t ) const
//    {
//        typename boost::range_iterator< const State >::type x = boost::begin( x_ );
//        typename boost::range_iterator< Deriv >::type dxdt = boost::begin( dxdt_ );
//
//        dxdt[0] = -c * x[0] * x[1] + d;
//        dxdt[1] = c * x[0] * x[1] - e * x[1] + g;
////        dxdt[0] = sigma * ( x[1] - x[0] );
////        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
////        dxdt[2] = -b * x[2] + x[0] * x[1];
//    }
//};


typedef std::vector< double > state_type;
//typedef boost::array< double , 2 > state_type;

/* The rhs of x' = f(x) */
void two_equation_model( const state_type &x , state_type &dxdt , const double /* t */ )
{
    dxdt[0] = -c * x[0] * x[1] + d;
    dxdt[1] = c * x[0] * x[1] - e * x[1] + g;
}
//]


int N = 100;
std::vector< double > t(N);
std::vector< double > x_gen(N);
std::vector< double > y_gen(N);

void generate_simulated_data()
{
    using namespace boost::numeric::odeint;


    for (int ii=0; ii<N; ii++) {
        t[ii] = ii;

        state_type x(2);
        x[0] = a;
        x[1] = b;

        size_t steps = integrate( two_equation_model , x , 0.0 , t[ii] , 0.1 );

        std::cout << "ii: " << x[0] << " " << x[1] << std::endl;

        x_gen[ii] = x[0];
        y_gen[ii] = x[1];
    }


}


int main() {
//    std::cout << "Hello, World!" << std::endl;
//
//    state_type x(2);
//    x[0] = a; // start at x=1.0, p=0.0
//    x[1] = b;
//
//    return 0;


    using namespace std;
    using namespace boost::numeric::odeint;


    //[ state_initialization
    state_type x(2);
    x[0] = a;
    x[1] = b;
    //]

    //[ integration
//    size_t steps = integrate( two_equation_model , x , 0.0 , 10.0 , 0.1 );
    //]

    generate_simulated_data();

    plt::plot(t, x_gen);
    plt::plot(t, y_gen);
    plt::show();

}
