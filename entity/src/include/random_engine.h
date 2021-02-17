//
// Created by leonard on 17.02.21.
//

#ifndef ENTITY_RANDOM_ENGINE_H
#define ENTITY_RANDOM_ENGINE_H

template<typename T>
class RandomEngine {

public:
    RandomEngine(T a, T b, double stddev, double mean) : lower_bound(a), upper_bound(b), stddev(stddev), mean(mean) {};
    virtual ~RandomEngine() = 0;
    virtual T SampleUniformInt() = 0;
    virtual T SampleNormalInt() = 0;

protected:

    T lower_bound;
    T upper_bound;

    double stddev;
    double mean;
};

#endif //ENTITY_RANDOM_ENGINE_H
