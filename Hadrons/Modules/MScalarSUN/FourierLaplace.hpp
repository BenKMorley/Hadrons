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
                                    std::string, base_filename);
};

class FourierLaplaceResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourierLaplaceResult,
                                    std::string, sink,
                                    std::string, source,
                                    std::vector<std::vector<Complex>>, correlator_p);
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
    std::string                                  Laplace_file = par().base_filename + "_Laplace_p" + "." + traj;
    std::string                                  Correlator_x_file = par().base_filename + "_correlator_x" + "." + traj;
    std::string                                  Correlator_p_file = par().base_filename + "_correlator_p" + "." + traj;

    envGetTmp(ComplexField, ftBuf);
    envGetTmp(ComplexField, correlator_p);
    envGetTmp(ComplexField, Laplace_p);
    envGetTmp(ComplexField, Laplace_temp);

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
        myfile.open(Correlator_p_file, std::ofstream::out | std::ofstream::trunc);
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L - 1; j++){
                site = {0, i, j};
                peekSite(read_buf, correlator_p, site);
                myfile << read_buf;
                myfile << ",";
            }
            site = {0, i, L};
            peekSite(read_buf, correlator_p, site);
            myfile << read_buf;
            myfile << "\n";
        }
        myfile.close();

        // qt = {0, 0, 1};
        // peekSite(read_buf, correlator_p, qt);
        // LOG(Message) << "At (0, 0, 1) we have a value of " << read_buf << " for result" << std::endl;

        // This sets Laplace_p to the position space correlator temporarily
        fft.FFT_all_dim(Laplace_p, correlator_p, FFT::backward);

        // Save two dimensions of the position space correlator
        myfile.open(Correlator_x_file, std::ofstream::out | std::ofstream::trunc);
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L - 1; j++){
                site = {0, i, j};
                peekSite(read_buf, Laplace_p, site);
                myfile << read_buf;
                myfile << ",";
            }
            site = {0, i, L};
            peekSite(read_buf, Laplace_p, site);
            myfile << read_buf;
            myfile << "\n";
        }
        myfile.close();

        LaplaceTransform3D(Laplace_p, Laplace_temp, offset, L);

        // Save two dimensions of the Laplace momentum space correlator
        myfile.open(Laplace_file, std::ofstream::out | std::ofstream::trunc);
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L - 1; j++){
                site = {0, i, j};
                peekSite(read_buf, Laplace_p, site);
                myfile << read_buf;
                myfile << ",";
            }
            site = {0, i, L};
            peekSite(read_buf, Laplace_p, site);
            myfile << read_buf;
            myfile << "\n";
        }
        myfile.close();
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_FourierLaplace_hpp_
