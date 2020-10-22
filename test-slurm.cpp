#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <chrono>
#include "io.hh"

double get_random(double lower_bound, double upper_bound) {
    static std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    static std::default_random_engine re;
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());
    return unif(re);
}

int main(int argc, char** argv) {
//     std::vector<std::string> arguments(argv + 1, argv + argc);
//     std::ofstream file(arguments.empty() ? "testing-slurm.dat" : arguments.at(0) + ".dat");
//     file << "Writing this to a file.\n";
//     file << (arguments.empty() ? "testing-slurm.txt" : arguments.at(0) + ".txt");
    double one = get_random(0, 1000);
    double two = get_random(0, 1000);
    double min1 = std::min(one, two);
    double max1 = std::max(one, two);

    std::vector<double> data;
    for (int i = 0; i < 1000; ++i) data.emplace_back(get_random(min1, max1));

    auto const [min, max] = std::minmax_element(begin(data), end(data));

    int bins = 100;
    std::vector<std::pair<double, int>> histogram(bins, {0., 0});
    for (auto i = 0; i < bins; ++i) histogram[i].first = *min + i * (*max - *min) / bins;
//     for (auto elem : histogram) print(elem);
    for (auto point : data) {
//         hnl::print(*min,*max,point);
        int i = static_cast<int>(std::floor(bins * (point - *min) / (*max - *min)));
        if (i < 0 || i >= 100) hnl::print("going to acces", i, point, *max);
        histogram[i].second++;
    }
    for (auto elem : histogram) std::cout << elem.first << ' ' << elem.second << '\n';
    for (int i = 0; i < histogram.size();++i) std::cout << i << ' ' << histogram[i].first << ' ' << histogram[i].second << '\n';
}
