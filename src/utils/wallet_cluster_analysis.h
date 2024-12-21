#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <limits>
#include <numeric>
#include <set>

class Point {
public:
    std::vector<double> coordinates;

    Point(const std::vector<double>& coords) : coordinates(coords) {}

    double& operator[](int index) {
        if (index < 0 || index >= coordinates.size()) {
            throw std::out_of_range("Index out of range");
        }
        return coordinates[index];
    }

    const double& operator[](int index) const {
        if (index < 0 || index >= coordinates.size()) {
            throw std::out_of_range("Index out of range");
        }
        return coordinates[index];
    }

    int size() const {
        return coordinates.size();
    }

    bool operator==(const Point& other) const {
        return coordinates == other.coordinates;
    }
};

namespace utils {

    double euclideanDistance(const Point& p1, const Point& p2);
    double manhattanDistance(const Point& p1, const Point& p2);
    Point mean(const std::vector<Point>& points);
    std::vector<double> standardDeviation(const std::vector<Point>& points);
    std::vector<double> variance(const std::vector<Point>& points);
    std::vector<std::vector<double>> covarianceMatrix(const std::vector<Point>& points);
    double determinant(const std::vector<std::vector<double>>& matrix);
    double cofactor(const std::vector<std::vector<double>>& matrix, int row, int col);
    std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matrix);
    double gaussianProbability(const Point& point, const Point& mean, const std::vector<std::vector<double>>& covariance);
    std::vector<std::vector<Point>> kMeans(const std::vector<Point>& points, int k, int maxIterations = 100);
    std::vector<std::vector<Point>> hierarchicalClustering(const std::vector<Point>& points, int k);
    std::vector<Point> regionQuery(const std::vector<Point>& points, int index, double eps);
    void expandCluster(const std::vector<Point>& points, int index, const std::vector<Point>& neighbors, int clusterId,
                       std::vector<std::vector<Point>>& clusters, std::vector<bool>& visited, std::vector<int>& clusterIds,
                       double eps, int minPts);
    std::vector<std::vector<Point>> dbscan(const std::vector<Point>& points, double eps, int minPts);
    std::vector<std::vector<Point>> gmm(const std::vector<Point>& points, int k, int maxIterations = 100);
}

#endif // EXAMPLE_H

