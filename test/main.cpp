#include <array>
#include <grid/core.hpp>
#include <iostream>
#include <unit/double.hpp>

int main()
{
    using namespace Unit;

    Grid::GridVector<Complex<Length>, double, 2> hoge{{-1.0, 1.0, 2}, {-1.0, 1.0, 2}};

    double* fuga = (double*)hoge.data();

    for (int i = 0; i < 8; ++i) {
        fuga[i] = 0.1 * (i + 1);
    }

    for (auto& h : hoge) {
        std::cout << h << std::endl;
    }

    return 0;
}
