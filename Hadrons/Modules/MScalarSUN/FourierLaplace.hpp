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
                                    std::string, output);
};

class FourierLaplaceResult: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FourierLaplaceResult,
                                    std::string, sink,
                                    std::string, source,
                                    std::vector<std::vector<Complex>>, data);
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
    envTmpLat(ComplexField, "resultBuf");
    envTmpLat(ComplexField, "resultConnectedBuf");
    envTmpLat(ComplexField, "op2ShiftBuf");
    envTmpLat(ComplexField, "windowField");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename SImpl>
void TFourierLaplace<SImpl>::execute(void)
{

    const unsigned int                           nd      = env().getNd();
    const unsigned int                           nt      = env().getDim().back();
    const unsigned int                           nop     = par().op.size();
    double                                       partVol = 1.;

    std::set<std::string>                        ops;
    std::vector<FourierLaplaceResult>            result;
    FFT                                          fft(envGetGrid(Field));
    TComplex                                     buf1, buf2, wbuf1;
    // std::vector<Complex>                         res(nt, 0.);
    Complex                                      read_buf;
    std::vector<int>                             qt(nd,0);
    Complex                                      bufZeroMode(0., 0.);
    int                                          L = nt;
    double                                       r, window_value;
    Coordinate                                   v_shift({0,0,0});


    envGetTmp(ComplexField, ftBuf);
    envGetTmp(ComplexField, resultBuf);
    envGetTmp(ComplexField, resultConnectedBuf);
    envGetTmp(ComplexField, op2ShiftBuf);
    envGetTmp(ComplexField, windowField);

    for (auto &p : par().op)
    {   
        LOG(Message) << "  <" << p.first << " " << p.second << ">" << std::endl;

        //Remove mean of fields (improves tatistics?)
        bufZeroMode = TensorRemove(sum(op1));
        bufZeroMode = TensorRemove(sum(op2));

        auto &op1 = envGet(ComplexField, p.first);
        auto &op2 = envGet(ComplexField, p.second);

        op1 = op1 - bufZeroMode / static_cast<double>(L * L * L);
        op2 = op2 - bufZeroMode / static_cast<double>(L * L * L);

        qt = {0, 0, 0};
        peekSite(read_buf, op1, qt);
        LOG(Message) << "At (0, 0, 0) we have a value of " << read_buf << " for Operator 1" << std::endl;
        peekSite(read_buf, op2, qt);
        LOG(Message) << "At (0, 0, 0) we have a value of " << read_buf << " for Operator 2" << std::endl;

        qt = {0, 0, 1};
        peekSite(read_buf, op1, qt);
        LOG(Message) << "At (0, 0, 1) we have a value of " << read_buf << " for Operator 1" << std::endl;
        peekSite(read_buf, op2, qt);
        LOG(Message) << "At (0, 0, 1) we have a value of " << read_buf << "for Operator 2" << std::endl;

        fft.FFT_all_dim(resultBuf, op1, FFT::forward);
        fft.FFT_all_dim(ftBuf , op2, FFT::forward);

        // Calculate the momentum space value of the two-point function
        resultBuf *= adj(ftBuf);

        qt = {0, 0, 0};
        peekSite(read_buf, resultBuf, qt);
        LOG(Message) << "At (0, 0, 0) we have a value of " << read_buf << " for result" << std::endl;

        qt = {0, 0, 1};
        peekSite(read_buf, resultBuf, qt);
        LOG(Message) << "At (0, 0, 1) we have a value of " << read_buf << " for result" << std::endl;

        fft.FFT_all_dim(ftBuf, resultBuf, FFT::backward);

        ftBuf *= windowField;

        fft.FFT_all_dim(resultBuf, ftBuf, FFT::forward);

        qt = {0, 0, 0};
        peekSite(bufZeroMode, resultBuf, qt);
        

        LOG(Message) << "Saving result..." << std::endl;

    //     FourierLaplaceResult r;
    //     r.sink    = p.first;
    //     r.source  = p.second;

    //     for (int q1 = 0; q1 < L; q1++)
    //     {    
    //         for (int q2 = 0; q2 < L; q2++)
    //         {
    //             qt[2] = q2; 
    //             peekSite(res[q2],resultBuf, qt);
    //         }
    //     }

    //     r.data    = res;
    //     result.push_back(r);

    }

    // saveResult(par().output, "twopt", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MScalarSUN_FourierLaplace_hpp_
