#ifndef Hadrons_MScalarSUN_FourierLaplace_hpp_
#define Hadrons_MScalarSUN_FourierLaplace_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         FourierLaplace                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class FourierLaplacePar: Serializable
{
public:
    typedef std::pair<std::string, std::string> OpPair;
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourierLaplacePar,
                                    std::vector<OpPair>, op,
                                    std::string, save_type,
                                    bool, debug,
                                    std::string, base_filename);
};

class FourierLaplaceResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourierLaplaceResult,
                                    std::string, sink,
                                    std::string, source,
                                    std::vector<std::vector<Complex>>, correlator_p,
                                    std::vector<std::vector<Complex>>, correlator_x,
                                    std::vector<std::vector<Complex>>, Laplace_p);
};

template <typename SImpl>
class TFourierLaplace: public Module<FourierLaplacePar>
{
public:
    typedef typename SImpl::Field         Field;
    typedef typename SImpl::ComplexField  ComplexField;
public:
    // constructor
    TFourierLaplace(const std::string name);
    // destructor
    virtual ~TFourierLaplace(void) {};
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

MODULE_REGISTER_TMP(FourierLaplaceSU2, TFourierLaplace<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceSU3, TFourierLaplace<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceSU4, TFourierLaplace<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceSU5, TFourierLaplace<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(FourierLaplaceSU6, TFourierLaplace<ScalarNxNAdjImplR<6>>, MScalarSUN);


/******************************************************************************
 *                 TFourierLaplace implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TFourierLaplace<SImpl>::TFourierLaplace(const std::string name)
: Module<FourierLaplacePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TFourierLaplace<SImpl>::getInput(void)
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
std::vector<std::string> TFourierLaplace<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TFourierLaplace<SImpl>::setup(void)
{
    const unsigned int nd = env().getDim().size();

    envTmpLat(ComplexField, "ftBuf");
    envTmpLat(ComplexField, "correlator_p");
    envTmpLat(ComplexField, "Laplace_p");
    envTmpLat(ComplexField, "Laplace_temp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TFourierLaplace<SImpl>::execute(void)
{
    const unsigned int                           nd      = env().getNd();
    const unsigned int                           nt      = env().getDim().back();
    const unsigned int                           nop     = par().op.size();
    double                                       partVol = 1.;
    const int                                    offset = 1;
    std::set<std::string>                        ops;
    std::vector<FourierLaplaceResult>            result;
    FFT                                          fft(envGetGrid(Field));
    std::vector<int>                             site(nd, 0.);
    Complex                         read_buf;
    std::vector<int>                             qt(nd,0);
    int                                          L = nt;
    FourierLaplaceResult                         r;
    std::string                                  traj = std::to_string(vm().getTrajectory());
    std::ofstream                                myfile;
    std::string                                  Laplace_file;
    std::string                                  Correlator_x_file;
    std::string                                  Correlator_p_file;
    std::vector<std::vector<Complex>>  Output_vector(L);

    envGetTmp(ComplexField, ftBuf);
    envGetTmp(ComplexField, correlator_p);
    envGetTmp(ComplexField, Laplace_p);
    envGetTmp(ComplexField, Laplace_temp);

    // This vector of vectors will be used to save the results into .h5 files
    for (int i = 0; i < L; i++) {
        std::vector<Complex> e(L, 0);
        Output_vector[i] = e;
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

        fft.FFT_all_dim(correlator_p, op1, FFT::forward);
        fft.FFT_all_dim(ftBuf , op2, FFT::forward);

        // Calculate the momentum space value of the two-point function
        correlator_p *= adj(ftBuf);

        if (par().save_type == "h5"){
            // Save two dimensions of this correlator
            for (int i = 0; i < L; i++){
                for (int j = 0; j < L; j++){
                    site = {0, i, j};
                    peekSite(read_buf, correlator_p, site);
                    Output_vector[i][j] = read_buf;
                }
            }

            r.correlator_p = Output_vector;
        }

        if (par().save_type == "csv"){
            myfile.open(Correlator_p_file, std::ofstream::out | std::ofstream::trunc);
            for (int i = 0; i < L; i++){
                for (int j = 0; j < L; j++){
                    site = {0, i, j};
                    peekSite(read_buf, correlator_p, site);
                    myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
                    myfile << ",";
                    myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
                    myfile << "\n";
                }
            }
            myfile.close();
        }

        fft.FFT_all_dim(Laplace_p, correlator_p, FFT::backward);

        if (par().save_type == "h5"){
            // Save two dimensions of the position space correlator
            for (int i = 0; i < L; i++){
                for (int j = 0; j < L; j++){
                    site = {0, i, j};
                    peekSite(read_buf, Laplace_p, site);
                    Output_vector[i][j] = read_buf;
                }
            }
            r.correlator_x = Output_vector;
        }

        if (par().save_type == "csv"){
            myfile.open(Correlator_x_file , std::ofstream::out | std::ofstream::trunc);
            for (int i = 0; i < L; i++){
                for (int j = 0; j < L; j++){
                    site = {0, i, j};
                    peekSite(read_buf, Laplace_p, site);
                    myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
                    myfile << ",";
                    myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
                    myfile << "\n";
                }
            }
            myfile.close();
        }

        LaplaceTransform3D(Laplace_p, Laplace_temp, offset, L, par().debug);

        if (par().save_type == "h5"){
            // Save two dimensions of the position space correlator
            for (int i = 0; i < L; i++){
                for (int j = 0; j < L; j++){
                    site = {0, i, j};
                    peekSite(read_buf, Laplace_p, site);
                    Output_vector[i][j] = read_buf;
                }
            }

            r.Laplace_p = Output_vector;
            result.push_back(r);
        }

        if (par().save_type == "csv"){
            myfile.open(Laplace_file, std::ofstream::out | std::ofstream::trunc);
            for (int i = 0; i < L; i++){
                for (int j = 0; j < L; j++){
                    site = {0, i, j};
                    peekSite(read_buf, Laplace_p, site);
                    myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
                    myfile << ",";
                    myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
                    myfile << "\n";
                }
            }
            myfile.close();
        }

    }

    if (par().save_type == "h5"){
        saveResult(par().base_filename, "FL", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_FourierLaplace_hpp_
