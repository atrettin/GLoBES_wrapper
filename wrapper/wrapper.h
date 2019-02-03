//#include <stdio.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <globes/globes.h>
#include <boost/python.hpp> 
#include <boost/shared_ptr.hpp>

#include "boost/python/extract.hpp"
#include "boost/python/numeric.hpp"

#include <numpy/ndarrayobject.h>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
extern "C" {
  #include <snu.h>
}
class GLoBESCalculator {

  public:
    GLoBESCalculator(std::string);
    void PrintParameters();
    void SetParameters(double, double, double, double, double, double); // th12, th13, th23, deltacp, dm21, dm31
    void SetParametersArr(boost::python::numeric::array&); 
    double VacuumProbability(int,int,int,double,double);
    // Matter probabilities
    double MatterProbability(int,int,int,double,double);
    double MatterProbabilityPrevBaseline(int,int,int,double);
    double MatterProbPDG(int,int,double,double);
    boost::python::object MatterProbPDGArr(boost::python::numeric::array&, boost::python::numeric::array&, boost::python::numeric::array&, boost::python::numeric::array&);
    double CalcPropDist(double);
    void SetEarthModel(boost::python::list&, boost::python::list&);
    void SetManualDensities(boost::python::list&, boost::python::list&);
    void PrintEarthModel();
    int CalcDensityLayers(double, bool);
    double GetBaselineInExperiment(int);
    void SetDepth(double);
    void SetAtmoHeight(double);
    void SetEarthRadius(double);
    void InitSteriles(int);
    void PrintParametersOrder();
    int GetNumOfOscParams();
    void PrintDensityProfile();
  private: 
    void ApplyParameters();
    int sterileMode;
    int CalcPREMCosines();
    double osc_params[12];
    double detDepth, atmoHeight, innerR, outerR, Rearth;
    double RadPREM;
    int PREM_layers, numInnerLayers, numOuterLayers;
    std::vector<double> inRads, outRads, inDens, outDens;
    std::vector<double> PREM_rad;
    std::vector<double> PREM_dens;
    std::vector<double> PREM_cosines;
};

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
BOOST_PYTHON_MODULE(GLoBES)
{
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  import_array();
  using namespace boost::python;
  class_<GLoBESCalculator>("GLoBESCalculator",  init<const std::string>())
    .def("SetParameters", &GLoBESCalculator::SetParameters)
    .def("SetParametersArr", &GLoBESCalculator::SetParametersArr)
    .def("PrintParameters", &GLoBESCalculator::PrintParameters)
    .def("VacuumProbability", &GLoBESCalculator::VacuumProbability)
    .def("CalcPropDist", &GLoBESCalculator::CalcPropDist)
    .def("SetEarthModel", &GLoBESCalculator::SetEarthModel)
    .def("SetEarthModel", &GLoBESCalculator::SetEarthModel)
    .def("PrintEarthModel", &GLoBESCalculator::PrintEarthModel)
    .def("SetManualDensities", &GLoBESCalculator::SetManualDensities)
    .def("MatterProbability", &GLoBESCalculator::MatterProbability)
    .def("MatterProbPDG",&GLoBESCalculator::MatterProbPDG)
    .def("MatterProbPDGArr",&GLoBESCalculator::MatterProbPDGArr)
    .def("GetBaselineInExperiment", &GLoBESCalculator::GetBaselineInExperiment)
    .def("SetDepth"  ,&GLoBESCalculator::SetDepth)
    .def("SetAtmoHeight", &GLoBESCalculator::SetAtmoHeight) 
    .def("SetEarthRadius", &GLoBESCalculator::SetEarthRadius) 
    .def("InitSteriles", &GLoBESCalculator::InitSteriles)
    .def("PrintParametersOrder", &GLoBESCalculator::PrintParametersOrder)
    .def("GetNumOfOscParams", &GLoBESCalculator::GetNumOfOscParams)
    .def("PrintDensityProfile",&GLoBESCalculator::PrintDensityProfile)
    .def("MatterProbabilityPrevBaseline", &GLoBESCalculator::MatterProbabilityPrevBaseline)
  ;
}
