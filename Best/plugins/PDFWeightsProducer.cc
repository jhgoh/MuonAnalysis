#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/View.h"

#include <LHAPDF/LHAPDF.h>

#include <memory>
#include <vector>
#include <string>

using namespace std;

class PDFWeightsProducer : public edm::EDProducer
{
public:
  PDFWeightsProducer(const edm::ParameterSet& pset);
  void beginJob();
  void produce(edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::InputTag genInfoLabel_;
};

PDFWeightsProducer::PDFWeightsProducer(const edm::ParameterSet& pset)
{
  //genInfoLabel_ = pset.getParameter<edm::InputTag>("genEventInfo");
  genInfoLabel_ = edm::InputTag("generator");

  produces<std::vector<double> >();
}

void PDFWeightsProducer::beginJob()
{
  LHAPDF::initPDFSet(1, "cteq66.LHgrid");
}

void PDFWeightsProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<std::vector<double> > weights(new std::vector<double>);
  if ( event.isRealData() )
  {
    event.put(weights);
    return;
  }

  edm::Handle<GenEventInfoProduct> genInfoHandle;
  event.getByLabel(genInfoLabel_, genInfoHandle);

  const float q = genInfoHandle->pdf()->scalePDF;
  const int id1 = genInfoHandle->pdf()->id.first;
  const int id2 = genInfoHandle->pdf()->id.second;
  const double x1 = genInfoHandle->pdf()->x.first;
  const double x2 = genInfoHandle->pdf()->x.second;

  const double xpdf1 = LHAPDF::xfx(1, x1, q, id1);
  const double xpdf2 = LHAPDF::xfx(1, x2, q, id2);
  const double w0 = xpdf1*xpdf2;

  for ( unsigned int i=1, n=LHAPDF::numberPDF(1); i<=n; ++i )
  {
    LHAPDF::usePDFMember(1, i);
    const double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
    const double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
    weights->push_back(xpdf1_new*xpdf2_new/w0);
  }

  event.put(weights);
}

DEFINE_FWK_MODULE(PDFWeightsProducer);

