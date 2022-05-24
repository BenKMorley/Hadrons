#ifndef Hadrons_MScalarSUN_FourierLaplaceQmax_hpp_
#define Hadrons_MScalarSUN_FourierLaplaceQmax_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

#define PI           3.14159265358979323846

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         FourierLaplaceQmax                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class FourierLaplaceQmaxPar: Serializable
{
public:
    typedef std::pair<std::string, std::string> OpPair;
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourierLaplaceQmaxPar,
                                    std::vector<OpPair>, op,
                                    std::string, save_type,
                                    bool, debug,
                                    std::string, base_filename,
                                    double, qmax);
};

class FourierLaplaceQmaxResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourierLaplaceQmaxResult,
                                    std::string, sink,
                                    std::string, source,
                                    std::vector<std::vector<Complex>>, correlator_p,
                                    std::vector<std::vector<Complex>>, correlator_x,
                                    std::vector<std::vector<Complex>>, Laplace_p,
                                    std::vector<std::vector<std::vector<Complex>>>, P_FULL_3D);
};

template <typename SImpl>
class TFourierLaplaceQmax: public Module<FourierLaplaceQmaxPar>
{
public:
    typedef typename SImpl::Field         Field;
    typedef typename SImpl::ComplexField  ComplexField;
public:
    // constructor
    TFourierLaplaceQmax(const std::string name);
    // destructor
    virtual ~TFourierLaplaceQmax(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<std::vector<int>> mom_;
};

MODULE_REGISTER_TMP(FourierLaplaceQmaxSU2, TFourierLaplaceQmax<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceQmaxSU3, TFourierLaplaceQmax<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceQmaxSU4, TFourierLaplaceQmax<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceQmaxSU5, TFourierLaplaceQmax<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceQmaxSU6, TFourierLaplaceQmax<ScalarNxNAdjImplR<6>>, MScalarSUN);


/******************************************************************************
 *                 TFourierLaplaceQmax implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TFourierLaplaceQmax<SImpl>::TFourierLaplaceQmax(const std::string name)
: Module<FourierLaplaceQmaxPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TFourierLaplaceQmax<SImpl>::getInput(void)
{
    std::vector<std::string> in;
    std::set<std::string>    ops;

    for (auto &p: par().op)
    {
        ops.insert(p.first);
        ops.insert(p.second);
    }
    for (auto &o: ops)
    {
        in.push_back(o);
    }

    return in;
}

template <typename SImpl>
std::vector<std::string> TFourierLaplaceQmax<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TFourierLaplaceQmax<SImpl>::setup(void)
{
    const unsigned int nd = env().getDim().size();

    envTmpLat(ComplexField, "ftBuf");
    envTmpLat(ComplexField, "correlator_p");
    envTmpLat(ComplexField, "Laplace_p");
    envTmpLat(ComplexField, "Laplace_temp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TFourierLaplaceQmax<SImpl>::execute(void)
{
    const unsigned int                           nd      = env().getNd();
    const unsigned int                           nt      = env().getDim().back();
    const unsigned int                           nop     = par().op.size();
    double                                       partVol = 1.;
    const int                                    offset = 1;
    std::set<std::string>                        ops;
    std::vector<FourierLaplaceQmaxResult>        result;
    FFT                                          fft(envGetGrid(Field));
    std::vector<int>                             site(nd, 0.);
    Complex                         read_buf;
    std::vector<int>                             qt(nd,0);
    int                                          L = nt;
    FourierLaplaceQmaxResult                     r;
    std::string                                  traj = std::to_string(vm().getTrajectory());
    std::ofstream                                myfile;
    std::string                                  Laplace_file;
    std::string                                  Correlator_x_file;
    std::string                                  Correlator_p_file;
    const int                                    pn_max = par().qmax * L / (2 * PI);
    std::vector<std::vector<Complex>>  Output_vector(L);
    std::vector<std::vector<std::vector<Complex>>>  Output_3D(L);
    std::vector<std::vector<Complex>>  Laplace_Output(pn_max);
    Complex                                      bufZeroMode(0., 0.);

    envGetTmp(ComplexField, ftBuf);
    envGetTmp(ComplexField, correlator_p);
    envGetTmp(ComplexField, Laplace_p);
    envGetTmp(ComplexField, Laplace_temp);

    // This vector of vectors will be used to save the results into .h5 files
    for (int i = 0; i < L; i++) {
        std::vector<Complex> e(L, 0);
        Output_vector[i] = e;
    }

    for (int i = 0; i < L; i++) {
        std::vector<std::vector<Complex>> e(L);
        Output_3D[i] = e;

        for (int j = 0; j < L; j++) {
            std::vector<Complex> f(L, 0);
            Output_3D[i][j] = f;
        }
    }

    for (int i = 0; i < pn_max; i++) {
        std::vector<Complex> e(pn_max, 0);
        Laplace_Output[i] = e;
    }

    for (auto &p : par().op){
        if (par().save_type == "h5"){
            r.sink    = p.first;
            r.source  = p.second;
        }

        Correlator_p_file = par().base_filename + "_correlator_p_" + p.first + p.second + "." + traj + ".csv";
        Correlator_x_file = par().base_filename + "_correlator_x_" + p.first + p.second + "." + traj + ".csv";
        Laplace_file = par().base_filename + "_Laplace_p_" + p.first + p.second + "." + traj + ".csv";

        LOG(Message) << "  <" << p.first << " " << p.second << ">" << std::endl;

        auto &op1 = envGet(ComplexField, p.first);
        auto &op2 = envGet(ComplexField, p.second);

        // // Debugging against Joseph's Code
        // const int                                    window_min = 0;
        // const int                                    window_max = 0;
        // fft.FFT_all_dim(correlator_p, op1, FFT::forward);
        // fft.FFT_all_dim(ftBuf , op2, FFT::forward);
        // correlator_p *= adj(ftBuf);

        // qt = {0, 0, 0};
        // peekSite(bufZeroMode, correlator_p, qt);
        // LOG(Message) << "Before Windowing:" << bufZeroMode << std::endl;

        // fft.FFT_all_dim(ftBuf, correlator_p, FFT::backward);

        // MakeWindowField(Laplace_p, -1, -1, 1000);

        // ftBuf *= Laplace_p;

        // fft.FFT_all_dim(correlator_p, ftBuf, FFT::forward);

        // peekSite(bufZeroMode, correlator_p, qt);
        // LOG(Message) << "After Windowing:" << bufZeroMode << std::endl;

        fft.FFT_all_dim(correlator_p, op1, FFT::forward);
        fft.FFT_all_dim(ftBuf , op2, FFT::forward);

        // Calculate the momentum space value of the two-point function
        correlator_p *= adj(ftBuf);

        if (par().debug == true){
            LOG(Message) << "DEBUG CORRELATOR_P" << std::endl;

            for (int j = 0; j < L; j++){
                qt = {0, 0, j};
                peekSite(bufZeroMode, correlator_p, qt);
                LOG(Message) << "p = (0, 0, " << j << ") = " <<  bufZeroMode << std::endl;
            }
        }

        // Save the full P correlator
        if (par().save_type == "h5"){
            // Save two dimensions of the position space correlator
            for (int i = 0; i < L; i++){
                for (int j = 0; j < L; j++){
                    for (int k = 0; k < L; k++){
                        qt = {i, j, k};
                        peekSite(read_buf, correlator_p, qt);
                        Output_3D[i][j][k] = read_buf;
                    }
                }
            }

            r.P_FULL_3D = Output_3D;
        }
        

        // if (par().save_type == "h5"){
        //     // Save two dimensions of this correlator
        //     for (int i = 0; i < L; i++){
        //         for (int j = 0; j < L; j++){
        //             site = {0, i, j};
        //             peekSite(read_buf, correlator_p, site);
        //             Output_vector[i][j] = read_buf;
        //         }
        //     }

        //     r.correlator_p = Output_vector;
        // }

        // if (par().save_type == "csv"){
        //     myfile.open(Correlator_p_file, std::ofstream::out | std::ofstream::trunc);
        //     for (int i = 0; i < L; i++){
        //         for (int j = 0; j < L; j++){
        //             site = {0, i, j};
        //             peekSite(read_buf, correlator_p, site);
        //             myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
        //             myfile << ",";
        //             myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
        //             myfile << "\n";
        //         }
        //     }
        //     myfile.close();
        // }

        // fft.FFT_all_dim(Laplace_p, correlator_p, FFT::backward);

        // if (par().save_type == "h5"){
        //     // Save two dimensions of the position space correlator
        //     for (int i = 0; i < L; i++){
        //         for (int j = 0; j < L; j++){
        //             site = {0, i, j};
        //             peekSite(read_buf, Laplace_p, site);
        //             Output_vector[i][j] = read_buf;
        //         }
        //     }
        //     r.correlator_x = Output_vector;
        // }

        // if (par().save_type == "csv"){
        //     myfile.open(Correlator_x_file , std::ofstream::out | std::ofstream::trunc);
        //     for (int i = 0; i < L; i++){
        //         for (int j = 0; j < L; j++){
        //             site = {0, i, j};
        //             peekSite(read_buf, Laplace_p, site);
        //             myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
        //             myfile << ",";
        //             myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
        //             myfile << "\n";
        //         }
        //     }
        //     myfile.close();
        // }

        // Check against the old method - WARNING : VERY SLOW
        // if (par().debug == true){
        //     LOG(Message) << "Old method" << std::endl;
        //     LOG(Message) << "============================================" << std::endl;

        //     LaplaceTransform3D(Laplace_p, Laplace_temp, offset, L, par().debug);

        //     // Reset the Laplace correlator
        //     fft.FFT_all_dim(Laplace_p, correlator_p, FFT::backward);
        //     LOG(Message) << "New method" << std::endl;
        //     LOG(Message) << "============================================" << std::endl;
        // }

        // LaplaceTransform3D_qmax(Laplace_p, Laplace_temp, offset, L, par().debug, par().qmax);

        // if (par().save_type == "h5"){
        //     // Save two dimensions of the position space correlator
        //     for (int i = 0; i < pn_max; i++){
        //         for (int j = 0; j < pn_max; j++){
        //             site = {0, i, j};
        //             peekSite(read_buf, Laplace_p, site);
        //             Laplace_Output[i][j] = read_buf;
        //         }
        //     }

        //     r.Laplace_p = Laplace_Output;
        //     result.push_back(r);
        // }

        // if (par().save_type == "csv"){
        //     myfile.open(Laplace_file, std::ofstream::out | std::ofstream::trunc);
        //     for (int i = 0; i < L; i++){
        //         for (int j = 0; j < L; j++){
        //             site = {0, i, j};
        //             peekSite(read_buf, Laplace_p, site);
        //             myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
        //             myfile << ",";
        //             myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
        //             myfile << "\n";
        //         }
        //     }
        //     myfile.close();
        // }

        result.push_back(r);

    }

    if (par().save_type == "h5"){
        saveResult(par().base_filename, "FL", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_FourierLaplaceQmax_hpp_
