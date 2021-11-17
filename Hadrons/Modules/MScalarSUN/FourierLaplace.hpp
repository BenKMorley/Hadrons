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
                                    std::string, base_filename,
                                    std::string, base_dir);
};

class FourierLaplaceResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourierLaplaceResult,
                                    std::string, sink,
                                    std::string, source,
                                    std::vector<std::vector<std::complex<double>>>, correlator_p,
                                    std::vector<std::vector<std::complex<double>>>, correlator_x,
                                    std::vector<std::vector<std::complex<double>>>, Laplace_p);
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
    std::complex<double>                         read_buf;
    std::vector<int>                             qt(nd,0);
    int                                          L = nt;
    FourierLaplaceResult                         r;
    std::string                                  traj = std::to_string(vm().getTrajectory());
    std::ofstream                                myfile;
    std::string                                  Laplace_file;
    std::string                                  Correlator_x_file;
    std::string                                  Correlator_p_file;
    std::vector<std::vector<std::complex<double>>>  Output_vector(L);


    envGetTmp(ComplexField, ftBuf);
    envGetTmp(ComplexField, correlator_p);
    envGetTmp(ComplexField, Laplace_p);
    envGetTmp(ComplexField, Laplace_temp);

    // This vector of vectors will be used to save the results into .h5 files
    for (int i = 0; i < L; i++) {
        std::vector<std::complex<double>> e(L, 0);
        Output_vector[i] = e;
    }

    for (auto &p : par().op){   
        r.sink    = p.first;
        r.source  = p.second;

        LOG(Message) << "  <" << p.first << " " << p.second << ">" << std::endl;

        auto &op1 = envGet(ComplexField, p.first);
        auto &op2 = envGet(ComplexField, p.second);

        fft.FFT_all_dim(correlator_p, op1, FFT::forward);
        fft.FFT_all_dim(ftBuf , op2, FFT::forward);

        // Calculate the momentum space value of the two-point function
        correlator_p *= adj(ftBuf);

        // Save two dimensions of this correlator
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                site = {0, i, j};
                peekSite(read_buf, correlator_p, site);
                Output_vector[i][j] = read_buf;
            }
        }

        Correlator_p_file = par().base_dir + par().base_filename + "_correlator_p_" + p.first + p.second;
        
        // Start a writer in the appropriate directory
        // ResultWriter writer();
        LOG(Message) << "Saving in directory" << par().base_dir << std::endl;

        // Write out the 2D slice of the lattice
        // write(writer, Correlator_p_file, Output_vector);
        LOG(Message) << "Saving to file" << Correlator_p_file << std::endl;

        r.sink    = p.first;
        r.source  = p.second;
        r.correlator_p = Output_vector;

        result.push_back(r);

        // myfile.open(Correlator_p_file, std::ofstream::out | std::ofstream::trunc);
        // for (int i = 0; i < L; i++){
        //     for (int j = 0; j < L; j++){
        //         site = {0, i, j};
        //         peekSite(read_buf, correlator_p, site);
        //         myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
        //         myfile << ",";
        //         myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
        //         myfile << "\n";
        //     }
        // }
        // myfile.close();

        // qt = {0, 0, 1};
        // peekSite(read_buf, correlator_p, qt);
        // LOG(Message) << "At (0, 0, 1) we have a value of " << read_buf << " for result" << std::endl;
        // the position space correlator temporarily
        fft.FFT_all_dim(Laplace_p, correlator_p, FFT::backward);

        // Save two dimensions of the position space correlator
        // Correlator_x_file = par().base_filename + "_correlator_x_" + p.first + p.second + "." + traj + ".csv";
        // myfile.open(Correlator_x_file, std::ofstream::out | std::ofstream::trunc);
        // for (int i = 0; i < L; i++){
        //     for (int j = 0; j < L; j++){
        //         site = {0, i, j};
        //         peekSite(read_buf, Laplace_p, site);
        //         myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
        //         myfile << ",";
        //         myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
        //         myfile << "\n";
        //     }
        // }
        // myfile.close();

        LOG(Message) << "Hello 1" << std::endl;

        LaplaceTransform3D(Laplace_p, Laplace_temp, offset, L);

        // Save two dimensions of the Laplace momentum space correlator
        Laplace_file = par().base_filename + "_Laplace_p_" + p.first + p.second + "." + traj + ".csv";
        // myfile.open(Laplace_file, std::ofstream::out | std::ofstream::trunc);
        // for (int i = 0; i < L; i++){
        //     for (int j = 0; j < L; j++){
        //         site = {0, i, j};
        //         peekSite(read_buf, Laplace_p, site);
        //         myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.real();
        //         myfile << ",";
        //         myfile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << read_buf.imag();
        //         myfile << "\n";
        //     }
        // }
        // myfile.close();
    }
    saveResult(Correlator_p_file, "FL", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_FourierLaplace_hpp_
