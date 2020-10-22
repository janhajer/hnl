#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <chrono>
#include "io.hh"
#include "math.hh"

double get_random(double lower, double upper) {
    std::uniform_real_distribution<double> distribution(lower, upper);
    std::default_random_engine random_engine;
    random_engine.seed(std::chrono::system_clock::now().time_since_epoch().count());
    return distribution(random_engine);
}

int main() {
    double one = get_random(0, 10);
    double two = get_random(10000, 100000);
    double lower = std::min(one, two);
    double upper = std::max(one, two);

    std::vector<double> data;
    for (int i = 0; i < 1000; ++i) data.emplace_back(get_random(lower, upper));

    auto const [min, max] = std::minmax_element(begin(data), end(data));

    int bins = 100;
    std::vector<std::pair<double, int>> histogram(bins, {0., 0});
    for (auto i = 0; i < bins; ++i) {
        auto step = hnl::log_value(*min, *max, i, bins);
        hnl::print(step);
        histogram[i].first = step;
    }
    hnl::print("second loop");
    for (auto point : data) {
        auto value = hnl::log_step(*min, *max, point, bins);
        int i = static_cast<int>(std::floor(value));
        hnl::print(value, i);
        if (i == bins) --i;
        if (i < 0 || i >= 100) hnl::print("going to acces", i, point, *max);
        histogram[i].second++;
    }
    hnl::print("the histogram");
    for (auto elem : histogram) std::cout << elem.first << ' ' << elem.second << '\n';
}
