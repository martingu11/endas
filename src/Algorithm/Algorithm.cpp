#include <Endas/Algorithm/Algorithm.hpp>


using namespace std;
using namespace endas;


SequentialFilter::SequentialFilter()
{ }

SequentialFilter::~SequentialFilter()
{ }

void SequentialFilter::onResult(OnResultFn fn)
{
    mOnResultFn = fn;
}


SequentialEnsembleFilter::SequentialEnsembleFilter()
{ }

SequentialEnsembleFilter::~SequentialEnsembleFilter()
{ }

void SequentialEnsembleFilter::onResult(OnResultFn fn)
{
    mOnResultFn = fn;
}
