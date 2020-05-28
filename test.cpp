#include <array>
#include <iostream>
#include <unit/double.hpp>

int main()
{
    using namespace Unit;

    std::array<Complex<Length>, 4> hoge;

    double* fuga = (double*)hoge.data();

    for (int i = 0; i < 8; ++i) {
        fuga[i] = 0.1 * (i + 1);
    }

    for (auto& h : hoge) {
        std::cout << h << std::endl;
    }

    return 0;
}
