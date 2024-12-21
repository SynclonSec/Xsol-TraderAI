#include "wallet_cluster_analysis.cpp"

// Define a Point class to represent data points
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

// Utility functions
namespace utils {

    // Euclidean distance between two points
    double euclideanDistance(const Point& p1, const Point& p2) {
        if (p1.size() != p2.size()) {
            throw std::invalid_argument("Points must be of the same dimension");
        }
        double sum = 0.0;
        for (int i = 0; i < p1.size(); ++i) {
            sum += std::pow(p1[i] - p2[i], 2);
        }
        return std::sqrt(sum);
    }

    // Manhattan distance between two points
    double manhattanDistance(const Point& p1, const Point& p2) {
        if (p1.size() != p2.size()) {
            throw std::invalid_argument("Points must be of the same dimension");
        }
        double sum = 0.0;
        for (int i = 0; i < p1.size(); ++i) {
            sum += std::abs(p1[i] - p2[i]);
        }
        return sum;
    }

    // Mean of a vector of points
    Point mean(const std::vector<Point>& points) {
        if (points.empty()) {
            throw std::invalid_argument("Points vector must not be empty");
        }
        int dim = points[0].size();
        std::vector<double> meanCoords(dim, 0.0);
        for (const auto& point : points) {
            for (int i = 0; i < dim; ++i) {
                meanCoords[i] += point[i];
            }
        }
        for (int i = 0; i < dim; ++i) {
            meanCoords[i] /= points.size();
        }
        return Point(meanCoords);
    }

    // Standard deviation of a vector of points
    std::vector<double> standardDeviation(const std::vector<Point>& points) {
        if (points.empty()) {
            throw std::invalid_argument("Points vector must not be empty");
        }
        int dim = points[0].size();
        std::vector<double> stdDev(dim, 0.0);
        Point meanPoint = mean(points);
        for (const auto& point : points) {
            for (int i = 0; i < dim; ++i) {
                stdDev[i] += std::pow(point[i] - meanPoint[i], 2);
            }
        }
        for (int i = 0; i < dim; ++i) {
            stdDev[i] = std::sqrt(stdDev[i] / points.size());
        }
        return stdDev;
    }

    // Variance of a vector of points
    std::vector<double> variance(const std::vector<Point>& points) {
        if (points.empty()) {
            throw std::invalid_argument("Points vector must not be empty");
        }
        int dim = points[0].size();
        std::vector<double> var(dim, 0.0);
        Point meanPoint = mean(points);
        for (const auto& point : points) {
            for (int i = 0; i < dim; ++i) {
                var[i] += std::pow(point[i] - meanPoint[i], 2);
            }
        }
        for (int i = 0; i < dim; ++i) {
            var[i] /= points.size();
        }
        return var;
    }

    // Covariance matrix of a vector of points
    std::vector<std::vector<double>> covarianceMatrix(const std::vector<Point>& points) {
        if (points.empty()) {
            throw std::invalid_argument("Points vector must not be empty");
        }
        int dim = points[0].size();
        std::vector<std::vector<double>> covMatrix(dim, std::vector<double>(dim, 0.0));
        Point meanPoint = mean(points);
        for (const auto& point : points) {
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    covMatrix[i][j] += (point[i] - meanPoint[i]) * (point[j] - meanPoint[j]);
                }
            }
        }
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                covMatrix[i][j] /= points.size();
            }
        }
        return covMatrix;
    }

    // Determinant of a matrix
    double determinant(const std::vector<std::vector<double>>& matrix);

    // Cofactor of a matrix
    double cofactor(const std::vector<std::vector<double>>& matrix, int row, int col) {
        int n = matrix.size();
        std::vector<std::vector<double>> minor(n - 1, std::vector<double>(n - 1));
        int minorRow = 0, minorCol = 0;
        for (int i = 0; i < n; ++i) {
            if (i == row) continue;
            minorCol = 0;
            for (int j = 0; j < n; ++j) {
                if (j == col) continue;
                minor[minorRow][minorCol] = matrix[i][j];
                ++minorCol;
            }
            ++minorRow;
        }
        return std::pow(-1, row + col) * determinant(minor);
    }

    double determinant(const std::vector<std::vector<double>>& matrix) {
        int n = matrix.size();
        if (n == 1) {
            return matrix[0][0];
        }
        if (n == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }
        double det = 0.0;
        for (int i = 0; i < n; ++i) {
            det += matrix[0][i] * cofactor(matrix, 0, i);
        }
        return det;
    }

    // Inverse of a matrix
    std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& matrix) {
        int n = matrix.size();
        double det = determinant(matrix);
        if (det == 0) {
            throw std::invalid_argument("Matrix is singular");
        }
        std::vector<std::vector<double>> inv(n, std::vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                inv[j][i] = cofactor(matrix, i, j) / det;
            }
        }
        return inv;
    }

    // Gaussian probability density function
    double gaussianProbability(const Point& point, const Point& mean, const std::vector<std::vector<double>>& covariance) {
        int dim = point.size();
        double det = determinant(covariance);
        if (det == 0) {
            throw std::invalid_argument("Covariance matrix is singular");
        }
        std::vector<std::vector<double>> invCovariance = inverse(covariance);
        double exponent = 0.0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                exponent += (point[i] - mean[i]) * invCovariance[i][j] * (point[j] - mean[j]);
            }
        }
        double coefficient = 1.0 / (std::pow(2 * M_PI, dim / 2.0) * std::sqrt(det));
        return coefficient * std::exp(-0.5 * exponent);
    }

    // K-means clustering algorithm
    std::vector<std::vector<Point>> kMeans(const std::vector<Point>& points, int k, int maxIterations = 100) {
        if (points.empty() || k <= 0) {
            throw std::invalid_argument("Invalid input for k-means clustering");
        }
        std::vector<Point> centroids;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, points.size() - 1);

        // Initialize centroids randomly
        for (int i = 0; i < k; ++i) {
            centroids.push_back(points[dis(gen)]);
        }

        std::vector<std::vector<Point>> clusters(k);
        for (int iter = 0; iter < maxIterations; ++iter) {
            // Assign points to the nearest centroid
            for (const auto& point : points) {
                double minDist = std::numeric_limits<double>::max();
                int clusterIndex = 0;
                for (int j = 0; j < k; ++j) {
                    double dist = euclideanDistance(point, centroids[j]);
                    if (dist < minDist) {
                        minDist = dist;
                        clusterIndex = j;
                    }
                }
                clusters[clusterIndex].push_back(point);
            }

            // Update centroids
            for (int j = 0; j < k; ++j) {
                if (!clusters[j].empty()) {
                    centroids[j] = mean(clusters[j]);
                }
                clusters[j].clear();
            }
        }

        // Final assignment of points to clusters
        for (const auto& point : points) {
            double minDist = std::numeric_limits<double>::max();
            int clusterIndex = 0;
            for (int j = 0; j < k; ++j) {
                double dist = euclideanDistance(point, centroids[j]);
                if (dist < minDist) {
                    minDist = dist;
                    clusterIndex = j;
                }
            }
            clusters[clusterIndex].push_back(point);
        }

        return clusters;
    }

    // Hierarchical clustering algorithm (single linkage)
    std::vector<std::vector<Point>> hierarchicalClustering(const std::vector<Point>& points, int k) {
        if (points.empty() || k <= 0) {
            throw std::invalid_argument("Invalid input for hierarchical clustering");
        }
        std::vector<std::vector<Point>> clusters;
        for (const auto& point : points) {
            clusters.push_back({point});
        }

        while (clusters.size() > k) {
            double minDist = std::numeric_limits<double>::max();
            int clusterIndex1 = 0, clusterIndex2 = 0;
            for (int i = 0; i < clusters.size(); ++i) {
                for (int j = i + 1; j < clusters.size(); ++j) {
                    double dist = std::numeric_limits<double>::max();
                    for (const auto& p1 : clusters[i]) {
                        for (const auto& p2 : clusters[j]) {
                            double tempDist = euclideanDistance(p1, p2);
                            if (tempDist < dist) {
                                dist = tempDist;
                            }
                        }
                    }
                    if (dist < minDist) {
                        minDist = dist;
                        clusterIndex1 = i;
                        clusterIndex2 = j;
                    }
                }
            }

            // Merge the closest clusters
            clusters[clusterIndex1].insert(clusters[clusterIndex1].end(), clusters[clusterIndex2].begin(), clusters[clusterIndex2].end());
            clusters.erase(clusters.begin() + clusterIndex2);
        }

        return clusters;
    }

    // Region query for DBSCAN
    std::vector<Point> regionQuery(const std::vector<Point>& points, int index, double eps) {
        std::vector<Point> neighbors;
        for (int i = 0; i < points.size(); ++i) {
            if (i != index && euclideanDistance(points[index], points[i]) <= eps) {
                neighbors.push_back(points[i]);
            }
        }
        return neighbors;
    }

    // Expand cluster for DBSCAN
    void expandCluster(const std::vector<Point>& points, int index, const std::vector<Point>& neighbors, int clusterId,
                       std::vector<std::vector<Point>>& clusters, std::vector<bool>& visited, std::vector<int>& clusterIds,
                       double eps, int minPts) {
        clusters[clusterId].push_back(points[index]);
        clusterIds[index] = clusterId;

        std::set<int> seeds = {index};
        for (const auto& neighbor : neighbors) {
            int neighborIndex = std::find(points.begin(), points.end(), neighbor) - points.begin();
            if (!visited[neighborIndex]) {
                visited[neighborIndex] = true;
                std::vector<Point> neighborNeighbors = regionQuery(points, neighborIndex, eps);
                if (neighborNeighbors.size() >= minPts) {
                    seeds.insert(neighborIndex);
                }
            }
            if (clusterIds[neighborIndex] == -1) {
                clusterIds[neighborIndex] = clusterId;
                clusters[clusterId].push_back(neighbor);
            }
        }

        for (int seedIndex : seeds) {
            std::vector<Point> seedNeighbors = regionQuery(points, seedIndex, eps);
            if (seedNeighbors.size() >= minPts) {
                expandCluster(points, seedIndex, seedNeighbors, clusterId, clusters, visited, clusterIds, eps, minPts);
            }
        }
    }

    // DBSCAN clustering algorithm
    std::vector<std::vector<Point>> dbscan(const std::vector<Point>& points, double eps, int minPts) {
        if (points.empty() || eps <= 0 || minPts <= 0) {
            throw std::invalid_argument("Invalid input for DBSCAN clustering");
        }
        std::vector<std::vector<Point>> clusters;
        std::vector<bool> visited(points.size(), false);
        std::vector<int> clusterIds(points.size(), -1);
        int clusterId = 0;

        for (int i = 0; i < points.size(); ++i) {
            if (!visited[i]) {
                visited[i] = true;
                std::vector<Point> neighbors = regionQuery(points, i, eps);
                if (neighbors.size() < minPts) {
                    clusterIds[i] = -1; // Mark as noise
                } else {
                    clusters.push_back({});
                    expandCluster(points, i, neighbors, clusterId, clusters, visited, clusterIds, eps, minPts);
                    clusterId++;
                }
            }
        }

        return clusters;
    }

    // Gaussian Mixture Model (GMM) clustering algorithm
    std::vector<std::vector<Point>> gmm(const std::vector<Point>& points, int k, int maxIterations = 100) {
        if (points.empty() || k <= 0) {
            throw std::invalid_argument("Invalid input for GMM clustering");
        }
        std::vector<Point> means;
        std::vector<std::vector<std::vector<double>>> covariances;
        std::vector<double> weights(k, 1.0 / k);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, points.size() - 1);

        // Initialize means randomly
        for (int i = 0; i < k; ++i) {
            means.push_back(points[dis(gen)]);
            covariances.push_back(covarianceMatrix({points[dis(gen)]}));
        }

        std::vector<std::vector<Point>> clusters(k);
        for (int iter = 0; iter < maxIterations; ++iter) {
            // E-step: Compute responsibilities
            std::vector<std::vector<double>> responsibilities(points.size(), std::vector<double>(k, 0.0));
            for (int i = 0; i < points.size(); ++i) {
                double sum = 0.0;
                for (int j = 0; j < k; ++j) {
                    double exponent = gaussianProbability(points[i], means[j], covariances[j]);
                    responsibilities[i][j] = weights[j] * exponent;
                    sum += responsibilities[i][j];
                }
                for (int j = 0; j < k; ++j) {
                    responsibilities[i][j] /= sum;
                }
            }

            // M-step: Update means, covariances, and weights
            for (int j = 0; j < k; ++j) {
                double weightSum = 0.0;
                std::vector<double> newMean(points[0].size(), 0.0);
                std::vector<std::vector<double>> newCovariance(points[0].size(), std::vector<double>(points[0].size(), 0.0));
                for (int i = 0; i < points.size(); ++i) {
                    weightSum += responsibilities[i][j];
                    for (int d = 0; d < points[0].size(); ++d) {
                        newMean[d] += responsibilities[i][j] * points[i][d];
                    }
                }
                for (int d = 0; d < points[0].size(); ++d) {
                    newMean[d] /= weightSum;
                }
                means[j] = Point(newMean);
                weights[j] = weightSum / points.size();

                for (int i = 0; i < points.size(); ++i) {
                    for (int d1 = 0; d1 < points[0].size(); ++d1) {
                        for (int d2 = 0; d2 < points[0].size(); ++d2) {
                            newCovariance[d1][d2] += responsibilities[i][j] *
                                                      (points[i][d1] - means[j][d1]) *
                                                      (points[i][d2] - means[j][d2]);
                        }
                    }
                }
                for (int d1 = 0; d1 < points[0].size(); ++d1) {
                    for (int d2 = 0; d2 < points[0].size(); ++d2) {
                        covariances[j][d1][d2] = newCovariance[d1][d2] / weightSum;
                    }
                }
            }
        }

        // Final assignment of points to clusters
        for (const auto& point : points) {
            double maxProb = -std::numeric_limits<double>::max();
            int clusterIndex = 0;
            for (int j = 0; j < k; ++j) {
                double prob = gaussianProbability(point, means[j], covariances[j]);
                if (prob > maxProb) {
                    maxProb = prob;
                    clusterIndex = j;
                }
            }
            clusters[clusterIndex].push_back(point);
        }

        return clusters;
    }
}

// Example usage
int main() {
    try {
        std::vector<Point> points = {
            Point({1.0, 2.0}),
            Point({2.0, 3.0}),
            Point({5.0, 8.0}),
            Point({8.0, 8.0}),
            Point({1.0, 0.5}),
            Point({9.0, 11.0}),
            Point({8.0, 2.0}),
            Point({10.0, 2.0}),
            Point({9.0, 3.0})
        };

        std::vector<std::vector<Point>> kmeansClusters = utils::kMeans(points, 3);
        std::cout << "K-means Clustering:" << std::endl;
        for (int i = 0; i < kmeansClusters.size(); ++i) {
            std::cout << "Cluster " << i + 1 << ": ";
            for (const auto& point : kmeansClusters[i]) {
                std::cout << "(";
                for (int j = 0; j < point.size(); ++j) {
                    std::cout << point[j];
                    if (j < point.size() - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
        }

        std::vector<std::vector<Point>> hierarchicalClusters = utils::hierarchicalClustering(points, 3);
        std::cout << "Hierarchical Clustering:" << std::endl;
        for (int i = 0; i < hierarchicalClusters.size(); ++i) {
            std::cout << "Cluster " << i + 1 << ": ";
            for (const auto& point : hierarchicalClusters[i]) {
                std::cout << "(";
                for (int j = 0; j < point.size(); ++j) {
                    std::cout << point[j];
                    if (j < point.size() - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
        }

        std::vector<std::vector<Point>> dbscanClusters = utils::dbscan(points, 1.5, 3);
        std::cout << "DBSCAN Clustering:" << std::endl;
        for (int i = 0; i < dbscanClusters.size(); ++i) {
            std::cout << "Cluster " << i + 1 << ": ";
            for (const auto& point : dbscanClusters[i]) {
                std::cout << "(";
                for (int j = 0; j < point.size(); ++j) {
                    std::cout << point[j];
                    if (j < point.size() - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
        }

        std::vector<std::vector<Point>> gmmClusters = utils::gmm(points, 3);
        std::cout << "GMM Clustering:" << std::endl;
        for (int i = 0; i < gmmClusters.size(); ++i) {
            std::cout << "Cluster " << i + 1 << ": ";
            for (const auto& point : gmmClusters[i]) {
                std::cout << "(";
                for (int j = 0; j < point.size(); ++j) {
                    std::cout << point[j];
                    if (j < point.size() - 1) {
                        std::cout << ", ";
                    }
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}

