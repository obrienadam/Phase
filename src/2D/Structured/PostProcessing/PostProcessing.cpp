#include "PostProcessing.h"

PostProcessing::PostProcessing(const Input &input, const Solver &solver)
    : _viewer(input, solver), _iter(0),
      _fileWriteFrequency(input.postProcessingInput().get<int>(
          "PostProcessing.fileWriteFrequency")) {}

void PostProcessing::compute(Scalar time, bool force) {
  if (_iter++ % _fileWriteFrequency == 0 || force) {
    _viewer.write(time);
  }
}
