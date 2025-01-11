//
// By: Robert Vogel
// Affiliation: Palmer Lab at UCSD
// Date: 2025-01-09
//
//
// Acknowledgment
//
// Code design and original version completed by Robert Vogel,
// reviewed by Claude Sonnet, the AI assistant from Anthropic
// (Jan 2025), with minor recommendations incorporated.
//
//

#include <array>

class Matrix
{
public:
    Matrix(size_t, size_t);
    double operator()(size_t, size_t) const;
    void operator()(size_t, size_t);

    size_t size() const;
    std::array<size_t,2> dims() const;

private:
    const size_t nrow_;
    const size_t mcol_;
    std::unique_ptr<double[]> data_;
    size_t mat_idx_to_array_(size_t, size_t);
}
