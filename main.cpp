#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include <nlopt.h>

#define WITHOUT_NUMPY

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;


int N = 100;

typedef struct {
    int length;
    std::vector<double> t;
    std::vector<double> x;
    std::vector<double> y;
} generated_data_t;

typedef struct {
    double a;
    double b;
    double c;
    double d;
    double e;
    double g;
} model_params_t;

typedef struct {
    generated_data_t *measurements;
    int step_count;
//    plt::Plot plot;
} objective_function_extra_t;

model_params_t truth_params = {
        .a = 250,
        .b = 125,
        .c = 0.0001,
        .d = 10,
        .e = 0.01,
        .g = 0
};




typedef std::vector< double > state_type;
//typedef boost::array< double , 2 > state_type;

/* The rhs of x' = f(x) */
//void two_equation_model( const state_type &x , state_type &dxdt , const double t, model_params_t current_params )
//{
//    dxdt[0] = -current_params.c * x[0] * x[1] + current_params.d;
//    dxdt[1] = current_params.c * x[0] * x[1] - current_params.e * x[1] + current_params.g;
//}
//]

class two_eqn_model {

    model_params_t m_params;

public:
    explicit two_eqn_model( model_params_t params ) : m_params(params) { }

    void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
    {
        dxdt[0] = -m_params.c * x[0] * x[1] + m_params.d;
        dxdt[1] = m_params.c * x[0] * x[1] - m_params.e * x[1] + m_params.g;
    }
};

void generate_simulated_data(generated_data_t *data, model_params_t params)
{
    using namespace boost::numeric::odeint;

//    printf("generating data from params: a:%f, b:%f, c:%f, d:%f, e:%f, g:%f\n", params.a, params.b, params.c, params.d, params.e, params.g);


    for (int ii=0; ii<data->length; ii++) {
        data->t[ii] = ii;

        state_type x(2);
        x[0] = params.a;
        x[1] = params.b;

        two_eqn_model model(params);

        // TODO: Make this faster by only doing a single call to `integrate` with an observer to capture intermediate time points. Also fixed step size should be used.
        size_t steps = integrate( model , x , 0.0 , data->t[ii] , 0.1 );

//        std::cout << "ii: " << x[0] << " " << x[1] << std::endl;

        data->x[ii] = x[0];
        data->y[ii] = x[1];
    }


}

void plot_solution(generated_data_t *measurements, generated_data_t *fitted_solution, model_params_t *params)
{
    char title[256];
    sprintf(title, "Fitted Numeric Model to Simulated Data with Random Noise\na=%1.5f, b=%1.5f, c=%1.5f,\n d=%1.5f, e=%1.5f, g=%1.5f", params->a, params->b, params->c, params->d, params->e, params->g);
    plt::title(title);
    plt::plot(fitted_solution->t, fitted_solution->x);
    plt::plot(fitted_solution->t, fitted_solution->y);
    plt::scatter(measurements->t, measurements->x);
    plt::scatter(measurements->t, measurements->y);
    plt::show();
}

double objective_function(unsigned n, const double *x, double *grad, void *extra_data)
{
    using namespace std;
    auto *extra = (objective_function_extra_t *) extra_data;

    model_params_t params = {
            .a = x[0],
            .b = x[1],
            .c = x[2],
            .d = x[3],
            .e = x[4],
            .g = x[5]
    };

    generated_data_t current_data;
    current_data.length = N;
    current_data.t = std::vector<double>(N);
    current_data.x = std::vector<double>(N);
    current_data.y = std::vector<double>(N);
    generate_simulated_data(&current_data, params);

    double sse = 0;

    double residual_x;
    double residual_y;
    for (int ii = 0; ii < current_data.length; ii++) {
        residual_x = current_data.x[ii] - extra->measurements->x[ii];
        residual_y = current_data.y[ii] - extra->measurements->y[ii];
        sse = sse + residual_x * residual_x;
        sse = sse + residual_y * residual_y;
//        cout << "residuals: " << residual_x << ", " << residual_y << ", " << truth_data->x[ii] << ", " << truth_data->y[ii] << endl;
    }

    extra->step_count++;

    printf("nelmin step %i: a:%f, b:%f, c:%f, d:%f, e:%f, g:%f\n", extra->step_count, params.a, params.b, params.c, params.d, params.e, params.g);

//    printf("nelmin step " << extra->step_count << " - a:" << params.a << ", b:" << params.b << ", c:" << params.c << ", d:" << params.d << ", e:" << params.e << ", g:" << params.g << " sse = " << sse << endl;

//    throw 0;

//    if (extra->step_count % 20 == 0) {
//        plot_data(extra->truth_data, &current_data);
//
//    }
//    update_plot(extra->plot, &current_data);

    return sse;

}

// TODO: Make this white gaussian instead of uniform
double get_random_number()
{
    return ((float)rand()/RAND_MAX - 0.5) * 50;
}

void inject_noise(generated_data_t *data)
{
    for (int ii=0; ii<data->length; ii++) {
        data->x[ii] = data->x[ii] + get_random_number();
        data->y[ii] = data->y[ii] + get_random_number();
    }
}

// assumes dest is properly allocated
void copy_generated_data(generated_data_t *source, generated_data_t *dest)
{
    dest->length = source->length;
    for (int ii=0; ii<source->length; ii++) {
        dest->t[ii] = source->t[ii];
        dest->x[ii] = source->x[ii];
        dest->y[ii] = source->y[ii];
    }
}




int main() {

    using namespace std;
    using namespace boost::numeric::odeint;

    generated_data_t truth_data;
    generated_data_t truth_data_with_noise;
    truth_data.length = N;
    truth_data.t = std::vector<double>(N);
    truth_data.x = std::vector<double>(N);
    truth_data.y = std::vector<double>(N);
    truth_data_with_noise.t = std::vector<double>(N);
    truth_data_with_noise.x = std::vector<double>(N);
    truth_data_with_noise.y = std::vector<double>(N);
    generate_simulated_data(&truth_data, truth_params);
    copy_generated_data(&truth_data, &truth_data_with_noise);
    inject_noise(&truth_data_with_noise);


    nlopt_opt opt;

    opt = nlopt_create(NLOPT_LN_NELDERMEAD, 6); /* algorithm and dimensionality */

    objective_function_extra_t extra = {
            .measurements = &truth_data_with_noise,
            .step_count = 0
//            .plot = plot
    };
    nlopt_set_min_objective(opt, objective_function, &extra);

    double step_sizes[6];
    step_sizes[0] = 1;
    step_sizes[1] = 1;
    step_sizes[2] = 0.00001;
    step_sizes[3] = 0.1;
    step_sizes[4] = 0.001;
    step_sizes[5] = 0.00001;
    nlopt_set_initial_step(opt, step_sizes);

    nlopt_set_xtol_rel(opt, 1e-6);

    double guesses[6] = { truth_params.a, truth_params.b, truth_params.c, truth_params.d, truth_params.e, truth_params.g };  /* `*`some` `initial` `guess`*` */
    double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */
    cout << "starting nlopt" << std::endl;
    if (nlopt_optimize(opt, guesses, &minf) < 0) {
        printf("nlopt failed!\n");
        nlopt_destroy(opt);
        return -1;
    }

    model_params_t fitted_params = {
            .a = guesses[0],
            .b = guesses[1],
            .c = guesses[2],
            .d = guesses[3],
            .e = guesses[4],
            .g = guesses[5]
    };

    printf("found minimum at: a:%f, b:%f, c:%f, d:%f, e:%f, g:%f = %0.10g\n", fitted_params.a, fitted_params.b, fitted_params.c, fitted_params.d, fitted_params.e, fitted_params.g, minf);
    printf(" - actual params: a:%f, b:%f, c:%f, d:%f, e:%f, g:%f\n", truth_params.a, truth_params.b, truth_params.c, truth_params.d, truth_params.e, truth_params.g);




    generated_data_t fitted_solution;
    fitted_solution.length = truth_data.length;
    fitted_solution.t = std::vector<double>(truth_data.length);
    fitted_solution.x = std::vector<double>(truth_data.length);
    fitted_solution.y = std::vector<double>(truth_data.length);

    generate_simulated_data(&fitted_solution, fitted_params);


//    plot_data();

//    plt::show();

    plot_solution(&truth_data_with_noise, &fitted_solution, &fitted_params);


    nlopt_destroy(opt);


}
