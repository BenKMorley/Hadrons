#ifndef Hadrons_MScalarSUN_ThreeDimCorrelator_hpp_
#define Hadrons_MScalarSUN_ThreeDimCorrelator_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>

#define PI           3.14159265358979323846

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         ThreeDimCorrelator                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MScalarSUN)

class ThreeDimCorrelatorPar: Serializable
{
public:
    typedef std::pair<std::string, std::string> OpPair;
    GRID_SERIALIZABLE_CLASS_MEMBERS(ThreeDimCorrelatorPar,
                                    std::vector<OpPair>, op,
                                    std::string, save_type,
                                    bool, debug,
                                    std::string, base_filename,
                                    double, qmax);
};

class ThreeDimCorrelatorResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ThreeDimCorrelatorResult,
                                    std::string, sink,
                                    std::string, source,
                                    std::vector<std::vector<Complex>>, correlator_p,
                                    std::vector<std::vector<Complex>>, correlator_x,
                                    std::vector<std::vector<Complex>>, Laplace_p,
                                    std::vector<std::vector<std::vector<Complex>>>, P_FULL_3D);
};

template <typename SImpl>
class TThreeDimCorrelator: public Module<ThreeDimCorrelatorPar>
{
public:
    typedef typename SImpl::Field         Field;
    typedef typename SImpl::ComplexField  ComplexField;
public:
    // constructor
    TThreeDimCorrelator(const std::string name);
    // destructor
    virtual ~TThreeDimCorrelator(void) {};
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

MODULE_REGISTER_TMP(ThreeDimCorrelatorSU2, TThreeDimCorrelator<ScalarNxNAdjImplR<2>>, MScalarSUN);
MODULE_REGISTER_TMP(ThreeDimCorrelatorSU3, TThreeDimCorrelator<ScalarNxNAdjImplR<3>>, MScalarSUN);
MODULE_REGISTER_TMP(ThreeDimCorrelatorSU4, TThreeDimCorrelator<ScalarNxNAdjImplR<4>>, MScalarSUN);
MODULE_REGISTER_TMP(ThreeDimCorrelatorSU5, TThreeDimCorrelator<ScalarNxNAdjImplR<5>>, MScalarSUN);
MODULE_REGISTER_TMP(ThreeDimCorrelatorSU6, TThreeDimCorrelator<ScalarNxNAdjImplR<6>>, MScalarSUN);


/******************************************************************************
 *                 TThreeDimCorrelator implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename SImpl>
TThreeDimCorrelator<SImpl>::TThreeDimCorrelator(const std::string name)
: Module<ThreeDimCorrelatorPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename SImpl>
std::vector<std::string> TThreeDimCorrelator<SImpl>::getInput(void)
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
std::vector<std::string> TThreeDimCorrelator<SImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename SImpl>
void TThreeDimCorrelator<SImpl>::setup(void)
{
    const unsigned int nd = env().getDim().size();

    envTmpLat(ComplexField, "ftBuf");
    envTmpLat(ComplexField, "correlator_p");
    envTmpLat(ComplexField, "Laplace_p");
    envTmpLat(ComplexField, "Laplace_temp");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TThreeDimCorrelator<SImpl>::execute(void)
{
    const unsigned int                           nd      = env().getNd();
    const unsigned int                           nt      = env().getDim().back();
    const unsigned int                           nop     = par().op.size();
    double                                       partVol = 1.;
    const int                                    offset = 1;
    std::set<std::string>                        ops;
    std::vector<ThreeDimCorrelatorResult>        result;
    FFT                                          fft(envGetGrid(Field));
    std::vector<int>                             site(nd, 0.);
    Complex                                      read_buf;
    std::vector<int>                             qt(nd,0);
    int                                          L = nt;
    ThreeDimCorrelatorResult                     r;
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

        result.push_back(r);

    }

    if (par().save_type == "h5"){
        saveResult(par().base_filename, "FL", result);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_ThreeDimCorrelator_hpp_
