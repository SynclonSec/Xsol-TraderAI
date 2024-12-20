#include <cstring>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace graph {
	class AdjacencyMatrix {
                public:
                        AdjacencyMatrix(size_t n) : m_n(n), m_edges(n, std::vector<uint64_t>(n,0)) {
				
			}
			size_t getSize() const {
				return m_edges.size();
			}

			void clear() {
				memset(&m_edges, 0, getSize());
			}

			void addEdge(size_t from, size_t to, uint64_t weight) {
				if (from >= getSize() || to >= getSize()) {
					throw std::out_of_range("AdjacencyMatrix: Invalid index");
				}

				m_edges[from][to] = m_edges[to][from] = weight;
			}

		    	size_t getEdgeCount() const {
        			size_t count = 0;
        			for (size_t i = 0; i < getSize(); ++i) {
					// j starts from i+1 because we count each edge only once for undirected graphs

            				for (size_t j = i + 1; j < getSize(); ++j) {
						if (m_edges[i][j] != 0) {
                    					++count;
                				}
            				}
        			}
        			return count;
    			}

		private:
			uint64_t m_n;
			std::vector<std::vector<uint64_t>> m_edges;
	};
};
